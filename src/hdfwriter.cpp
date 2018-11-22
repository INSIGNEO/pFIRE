#include "hdfwriter.hpp"

#include <sstream>

#include "indexing.hpp"
#include "image.hpp"
#include "map.hpp"

const std::vector<std::string> HDFWriter::_components({"x", "y", "z"}); 

HDFWriter::HDFWriter(const std::string& filename, const MPI_Comm& comm, bool truncate_existing)
  : _comm(comm), _filename(filename), _file_h(-1)
{
  hid_t plist_h = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_h, _comm, MPI_INFO_NULL);
  // Open file with parallel properties as set in ctr
  if (truncate_existing)
  {
    _file_h = create_or_truncate(filename, plist_h);
  }
  else
  {
    _file_h = create_or_open(filename, plist_h);
  }
  H5Pclose(plist_h);
}

HDFWriter::~HDFWriter()
{
  H5Fclose(_file_h);
}


void HDFWriter::write_image(const Image& image, const std::string& groupname)
{
  MPI_Comm comm = image.comm();

  int rank;
  MPI_Comm_rank(comm, &rank);

  //No need to set max size as want it to be same as given size.
  std::vector<hsize_t> imshape(image.shape().cbegin(), image.shape().cend());
  imshape.resize(image.ndim());
  std::reverse(imshape.begin(), imshape.end());
  hid_t fspace_h = H5Screate_simple(image.ndim(), imshape.data(), nullptr);

  hid_t dset_h = H5Dcreate(_file_h, groupname.c_str(), H5T_NATIVE_DOUBLE, fspace_h,
                           H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  std::vector<hsize_t> offset = image.mpi_get_offset<hsize_t>();
  offset.resize(image.ndim());
  std::reverse(offset.begin(), offset.end());
  std::vector<hsize_t> chunksize = image.mpi_get_chunksize<hsize_t>();
  chunksize.resize(image.ndim());
  std::reverse(chunksize.begin(), chunksize.end());

  H5Sselect_hyperslab(fspace_h, H5S_SELECT_SET, offset.data(), nullptr, chunksize.data(), nullptr);

  hid_t plist_h = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_h, H5FD_MPIO_COLLECTIVE);

  hid_t dspace_h = H5Screate_simple(image.ndim(), chunksize.data(), nullptr);

  const floating *imgdata = image.get_raw_data_ro();
  H5Dwrite(dset_h, H5T_NATIVE_DOUBLE, dspace_h, fspace_h, plist_h, imgdata);
  image.release_raw_data_ro(imgdata);

  H5Dclose(dset_h);
  H5Sclose(fspace_h);
  H5Sclose(dspace_h);
  H5Pclose(plist_h);

}


void HDFWriter::write_map(const Map& map, const std::string& groupname)
{
  MPI_Comm comm = map.comm();

  int rank;
  MPI_Comm_rank(comm, &rank);

  //No need to set max size as want it to be same as given size.
  std::vector<hsize_t> mapshape(map.shape().cbegin(), map.shape().cend());
  mapshape.resize(map.ndim());
  std::reverse(mapshape.begin(), mapshape.end());

  // Set up for collective writes
  hid_t plist_h = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_h, H5FD_MPIO_COLLECTIVE);

  hid_t mgroup_h = H5Gcreate(_file_h, groupname.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (mgroup_h < 0)
  {
    std::ostringstream errstr;
    errstr << "Failed to create group " << groupname << ".";
    throw std::runtime_error(errstr.str());
  }

  std::pair<integer, integer> own_range = map.get_displacement_ownershiprange();

  for(uinteger idx=0; idx<map.ndim(); idx++)
  {
    // Maps contained in group
    std::ostringstream dsetstr;
    dsetstr << groupname << "/" << _components[idx];
    std::string dsetname = dsetstr.str();
    PetscPrintf(_comm, "%s\n", dsetname.c_str());
    // Create dataspace and dataset to hold full data
    hid_t fspace_h = H5Screate_simple(map.ndim(), mapshape.data(), nullptr);
    hid_t dset_h = H5Dcreate(_file_h, dsetname.c_str(), H5T_NATIVE_DOUBLE, fspace_h,
                             H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    hsize_t lowloc = map.size()*idx;
    hsize_t highloc = map.size()*(idx+1)-1;
    lowloc = own_range.first < (integer)lowloc ? lowloc : own_range.first;
    highloc = own_range.second > (integer)highloc ? highloc : own_range.second;

    std::vector<hsize_t> offset(map.ndim(), 0);
    std::vector<hsize_t> chunksize(map.ndim(), 0);
    if (lowloc < highloc)
    {
      offset = unravel(lowloc, mapshape);
      chunksize = unravel(highloc, mapshape);

      std::transform(chunksize.cbegin(), chunksize.cend(), offset.cbegin(), chunksize.begin(),
                     [](hsize_t x, hsize_t y) -> hsize_t{return 1 + x - y;});
    }

    H5Sselect_hyperslab(fspace_h, H5S_SELECT_SET, offset.data(), nullptr, chunksize.data(), nullptr);

    hid_t dspace_h = H5Screate_simple(map.ndim(), chunksize.data(), nullptr);

    const floating *mapdata = map.get_raw_data_ro();
    const floating *mapstart = &mapdata[lowloc - own_range.first];
    H5Dwrite(dset_h, H5T_NATIVE_DOUBLE, dspace_h, fspace_h, plist_h, mapstart);
    map.release_raw_data_ro(mapdata);

    H5Dclose(dset_h);
    H5Sclose(fspace_h);
    H5Sclose(dspace_h);
  }
  H5Gclose(mgroup_h);
  H5Pclose(plist_h);
}

hid_t HDFWriter::create_or_open(std::string filename, hid_t file_props)
{
  // First supress HDF5 error printing
  herr_t (*old_err)(void*);
  void *old_err_data;
  H5Eget_auto1(&old_err, &old_err_data);
  H5Eset_auto1(nullptr, nullptr);

  // Open if available, will get file_h > 0 if successful
  hid_t file_h = H5Fopen(filename.c_str(), H5F_ACC_RDWR, file_props);
  if (file_h < 0)
  {
    file_h = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, file_props);
  }
  if (file_h < 0)
  {
    std::ostringstream err;
    err << "Failed to open output file " << filename << ".";
    throw std::runtime_error(err.str());
  }
  H5Eset_auto1(old_err, old_err_data);

  return file_h;
}

hid_t HDFWriter::create_or_truncate(std::string filename, hid_t file_props)
{
  // First supress HDF5 error printing
  herr_t (*old_err)(void*);
  void *old_err_data;
  H5Eget_auto1(&old_err, &old_err_data);
  H5Eset_auto1(nullptr, nullptr);

  // Open if available, will get file_h > 0 if successful
  hid_t file_h = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, file_props);
  if (file_h < 0)
  {
    std::ostringstream err;
    err << "Failed to open output file " << filename << ".";
    throw std::runtime_error(err.str());
  }
  H5Eset_auto1(old_err, old_err_data);

  return file_h;
}
