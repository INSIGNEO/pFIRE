//
//   Copyright 2019 University of Sheffield
//
//   Licensed under the Apache License, Version 2.0 (the "License");
//   you may not use this file except in compliance with the License.
//   You may obtain a copy of the License at
//
//       http://www.apache.org/licenses/LICENSE-2.0
//
//   Unless required by applicable law or agreed to in writing, software
//   distributed under the License is distributed on an "AS IS" BASIS,
//   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//   See the License for the specific language governing permissions and
//   limitations under the License.

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

void HDFWriter::write_3d_dataset_parallel(uinteger ndim, const std::vector<hsize_t>& fullshape,
    const std::vector<hsize_t>& chunkshape, const std::vector<hsize_t> &offset,
    const std::string& groupname, Vec& datavec)
{
  PetscPrintf(_comm, "Writing dataset: %s:%s\n", h5_filename.c_str(), groupname.c_str());
  // No need to set max size as want it to be same as given size.
  hid_t fspace_h = H5Screate_simple(ndim, fullshape.data(), nullptr);

  hid_t dset_h = H5Dcreate(_file_h, groupname.c_str(), H5T_NATIVE_DOUBLE, fspace_h, H5P_DEFAULT,
      H5P_DEFAULT, H5P_DEFAULT);

  H5Sselect_hyperslab(fspace_h, H5S_SELECT_SET, offset.data(), nullptr, chunkshape.data(), nullptr);

  hid_t plist_h = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_h, H5FD_MPIO_COLLECTIVE);

  hid_t dspace_h = H5Screate_simple(ndim, chunkshape.data(), nullptr);

  const floating* imgdata;
  PetscErrorCode perr = VecGetArrayRead(datavec, &imgdata);
  CHKERRXX(perr);
  H5Dwrite(dset_h, H5T_NATIVE_DOUBLE, dspace_h, fspace_h, plist_h, imgdata);
  perr = VecRestoreArrayRead(datavec, &imgdata);
  CHKERRXX(perr);

  H5Dclose(dset_h);
  H5Sclose(fspace_h);
  H5Sclose(dspace_h);
  H5Pclose(plist_h);
}

void HDFWriter::write_1d_dataset_rank0(hsize_t nval, const std::string& groupname, const floating* databuf)
{
  PetscPrintf(_comm, "Writing dataset: %s:%s\n", h5_filename.c_str(), groupname.c_str());
  int rank;
  MPI_Comm_rank(_comm, &rank);
  // No need to set max size as want it to be same as given size.
  hid_t fspace_h = H5Screate_simple(1, &nval, nullptr);

  hid_t dset_h = H5Dcreate(_file_h, groupname.c_str(), H5T_NATIVE_DOUBLE, fspace_h, H5P_DEFAULT,
      H5P_DEFAULT, H5P_DEFAULT);

  if(rank == 0)
  {
    H5Dwrite(dset_h, H5T_NATIVE_DOUBLE, H5S_ALL, fspace_h, H5P_DEFAULT, databuf);
  }

  H5Dclose(dset_h);
  H5Sclose(fspace_h);
}

std::string HDFWriter::write_image(const Image& image)
{
  // Sanity check communicators
  MPI_Comm comm = image.comm();
  if (comm != _comm)
  {
    std::ostringstream errss;
    errss << "Communicator mismatch between HDFWriter and provided image";
    throw InternalError(errss.str(), __FILE__, __LINE__);
  }

  std::string groupname = h5_groupname;
  if(groupname.empty())
  {
    groupname = "/registered";
  }

  int rank;
  MPI_Comm_rank(comm, &rank);

  std::vector<hsize_t> imgsizehsizet(image.shape().cbegin(), image.shape().cend());

  write_3d_dataset_parallel(image.ndim(), imgsizehsizet, image.mpi_get_chunksize<hsize_t>(),
      image.mpi_get_offset<hsize_t>(), groupname, *image.get_raw_data_row_major());

  std::string filepath = h5_filename + ":" + h5_groupname;
  return filepath;
}

std::string HDFWriter::write_map(const Map& map)
{
  // Sanity check communicators
  MPI_Comm comm = map.comm();
  if (comm != _comm)
  {
    std::ostringstream errss;
    errss << "Communicator mismatch between HDFWriter and provided image";
    throw InternalError(errss.str(), __FILE__, __LINE__);
  }

  std::string groupname = h5_groupname;
  if(groupname.empty())
  {
    groupname = "/map";
  }

  int rank;
  MPI_Comm_rank(comm, &rank);

  hid_t mgroup_h = H5Gcreate(_file_h, groupname.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (mgroup_h < 0)
  {
    std::ostringstream errstr;
    errstr << "Failed to create group " << h5_groupname << ".";
    throw WriterError(errstr.str());
  }

  for (uinteger idx = 0; idx < map.ndim(); idx++)
  {
    // Get dataset local dimensions
    auto corners = map.get_dmda_local_extents();
    std::vector<hsize_t> offset(corners.first.cbegin(), corners.first.cend());
    std::vector<hsize_t> chunksize(corners.second.cbegin(), corners.second.cend());
    /*
    PetscSynchronizedPrintf(_comm, "Rank %i: ofs: %i %i %i, shp %i %i %i\n", rank, offset[0], offset[1], offset[2],
                            chunksize[0], chunksize[1], chunksize[2]);
    PetscSynchronizedFlush(_comm, PETSC_STDOUT);
    */
    std::ostringstream dsetss;
    dsetss << groupname << "/" << _components[idx];
    std::string dsetname = dsetss.str();

    std::vector<hsize_t> mapsizehsizet(map.shape().cbegin(), map.shape().cend());

    write_3d_dataset_parallel(map.ndim(), mapsizehsizet, chunksize, offset, dsetname,
        *map.get_raw_data_row_major(idx));

    dsetss.clear();
    dsetss.str(std::string());
    dsetss << groupname << "/nodes_" << _components[idx];
    dsetname = dsetss.str();

    write_1d_dataset_rank0(map.shape()[idx], dsetname, map.node_locs()[idx].data());
  }
  H5Gclose(mgroup_h);

  std::string filepath = h5_filename + ":" + h5_groupname;
  return filepath;
}

void HDFWriter::open_or_create_h5()
{
  if (_file_h >= 0)
  {
    std::ostringstream err;
    err << "Attempted to reinitialize already open file, this is a pFIRE bug...";
    throw InternalError(err.str(), __FILE__, __LINE__);
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
    throw WriterError(err.str());
  }
  // Reset HDF5 error handling
  H5Eset_auto1(old_err, old_err_data);
  // Clean up old fileprops
  H5Pclose(file_props);
}
