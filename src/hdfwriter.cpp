#include "hdfwriter.hpp"

#include <algorithm>
#include <sstream>

#include "image.hpp"
#include "indexing.hpp"
#include "infix_iterator.hpp"
#include "map.hpp"

const std::string HDFWriter::writer_name = "hdf5";
const std::vector<std::string> HDFWriter::extensions = {".h5"};

HDFWriter::HDFWriter(const std::string& filespec, const MPI_Comm& comm)
  : BaseWriter(filespec, comm), h5_filename(filename), h5_groupname(extra_path), _file_h(-1)
{
  // Will open h5 file or throw
  open_or_create_h5();
}

HDFWriter::~HDFWriter()
{
  H5Fclose(_file_h);
}

void HDFWriter::write_image(const Image& image)
{
  // Sanity check communicators
  MPI_Comm comm = image.comm();
  if (comm != _comm)
  {
    std::ostringstream errss;
    errss << "Communicator mismatch between HDFWriter and provided image";
    throw std::runtime_error(errss.str());
  }

  int rank;
  MPI_Comm_rank(comm, &rank);

  // No need to set max size as want it to be same as given size.
  std::vector<hsize_t> imshape(image.shape().cbegin(), image.shape().cend());
  hid_t fspace_h = H5Screate_simple(image.ndim(), imshape.data(), nullptr);

  hid_t dset_h = H5Dcreate(_file_h, h5_groupname.c_str(), H5T_NATIVE_DOUBLE, fspace_h, H5P_DEFAULT,
      H5P_DEFAULT, H5P_DEFAULT);

  std::vector<hsize_t> offset = image.mpi_get_offset<hsize_t>();
  std::vector<hsize_t> chunksize = image.mpi_get_chunksize<hsize_t>();

  H5Sselect_hyperslab(fspace_h, H5S_SELECT_SET, offset.data(), nullptr, chunksize.data(), nullptr);

  hid_t plist_h = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_h, H5FD_MPIO_COLLECTIVE);

  hid_t dspace_h = H5Screate_simple(image.ndim(), chunksize.data(), nullptr);

  Vec_unique rm_data = image.get_raw_data_row_major();
  const floating* imgdata;
  PetscErrorCode perr = VecGetArrayRead(*rm_data, &imgdata);
  CHKERRABORT(_comm, perr);
  H5Dwrite(dset_h, H5T_NATIVE_DOUBLE, dspace_h, fspace_h, plist_h, imgdata);
  perr = VecRestoreArrayRead(*rm_data, &imgdata);
  CHKERRABORT(_comm, perr);

  H5Dclose(dset_h);
  H5Sclose(fspace_h);
  H5Sclose(dspace_h);
  H5Pclose(plist_h);
}

void HDFWriter::write_map(const Map& map)
{
  // Sanity check communicators
  MPI_Comm comm = map.comm();
  if (comm != _comm)
  {
    std::ostringstream errss;
    errss << "Communicator mismatch between HDFWriter and provided image";
    throw std::runtime_error(errss.str());
  }

  int rank;
  MPI_Comm_rank(comm, &rank);

  // No need to set max size as want it to be same as given size.
  std::vector<hsize_t> mapshape(map.shape().cbegin(), map.shape().cend());
  //  mapshape.resize(map.ndim());
  //  std::reverse(mapshape.begin(), mapshape.end());

  // Set up for collective writes
  hid_t plist_h = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_h, H5FD_MPIO_COLLECTIVE);

  hid_t mgroup_h = H5Gcreate(_file_h, h5_groupname.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (mgroup_h < 0)
  {
    std::ostringstream errstr;
    errstr << "Failed to create group " << h5_groupname << ".";
    throw std::runtime_error(errstr.str());
  }

  for (uinteger idx = 0; idx < map.ndim(); idx++)
  {
    // Maps contained in group
    std::ostringstream dsetstr;
    dsetstr << h5_groupname << "/" << _components[idx];
    std::string dsetname = dsetstr.str();
    PetscPrintf(_comm, "%s\n", dsetname.c_str());
    // Create dataspace and dataset to hold full data
    hid_t fspace_h = H5Screate_simple(map.ndim(), mapshape.data(), nullptr);
    hid_t dset_h = H5Dcreate(_file_h, dsetname.c_str(), H5T_NATIVE_DOUBLE, fspace_h, H5P_DEFAULT,
        H5P_DEFAULT, H5P_DEFAULT);

    auto corners = map.get_dmda_local_extents();
    std::vector<hsize_t> offset(corners.first.cbegin(), corners.first.cend());
    std::vector<hsize_t> chunksize(corners.second.cbegin(), corners.second.cend());
    Vec_unique dimdata = map.get_dim_data_dmda_blocked(idx);

    std::ostringstream ofs;
    std::copy(offset.cbegin(), offset.cend(), infix_ostream_iterator<integer>(ofs, ", "));
    std::ostringstream cks;
    std::copy(chunksize.cbegin(), chunksize.cend(), infix_ostream_iterator<integer>(cks, ", "));
    int rank;
    MPI_Comm_rank(_comm, &rank);

    PetscSynchronizedPrintf(
        _comm, "Rank %i: offset: %s, chunksize: %s\n", rank, ofs.str().c_str(), cks.str().c_str());
    PetscSynchronizedFlush(_comm, PETSC_STDOUT);

    H5Sselect_hyperslab(
        fspace_h, H5S_SELECT_SET, offset.data(), nullptr, chunksize.data(), nullptr);

    hid_t dspace_h = H5Screate_simple(map.ndim(), chunksize.data(), nullptr);

    const floating* mapdata;
    PetscErrorCode perr = VecGetArrayRead(*dimdata, &mapdata);
    CHKERRABORT(_comm, perr);
    H5Dwrite(dset_h, H5T_NATIVE_DOUBLE, dspace_h, fspace_h, plist_h, mapdata);
    perr = VecRestoreArrayRead(*dimdata, &mapdata);

    H5Dclose(dset_h);
    H5Sclose(fspace_h);
    H5Sclose(dspace_h);
  }
  H5Gclose(mgroup_h);
  H5Pclose(plist_h);
}

void HDFWriter::open_or_create_h5()
{
  if (_file_h >= 0)
  {
    std::ostringstream err;
    err << "Attempted to reinitialize already open file, this is a pFIRE bug...";
    throw std::runtime_error(err.str());
  }
  // Open file with parallel properties
  hid_t file_props = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(file_props, _comm, MPI_INFO_NULL);

  // First supress HDF5 error printing
  herr_t (*old_err)(void*);
  void* old_err_data;
  H5Eget_auto1(&old_err, &old_err_data);
  H5Eset_auto1(nullptr, nullptr);

  // Open if available, will get file_h > 0 if successful
  if(BaseWriter::check_truncated(h5_filename))
  {
    _file_h = H5Fopen(h5_filename.c_str(), H5F_ACC_RDWR, file_props);
  }
  else
  {
    _file_h = H5Fcreate(h5_filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, file_props);
    BaseWriter::mark_truncated(h5_filename);
  }
  if (_file_h < 0)
  {
    std::ostringstream err;
    err << "Failed to open output file " << h5_filename << ".";
    throw std::runtime_error(err.str());
  }
  // Reset HDF5 error handling
  H5Eset_auto1(old_err, old_err_data);

  H5Pclose(file_props);
}
