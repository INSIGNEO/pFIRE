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
#include "mapbase.hpp"

const std::string HDFWriter::writer_name = "hdf5";
const std::vector<std::string> HDFWriter::extensions = {".h5"};

HDFWriter::HDFWriter(const std::string& filespec, const MPI_Comm& comm)
  : BaseWriter(filespec, comm), h5_filename(this->filename()), h5_groupname(this->extra_path()), _file_h(-1)
{
  // Will open h5 file or throw
  open_or_create_h5();
}

HDFWriter::~HDFWriter() { H5Fclose(_file_h); }

void HDFWriter::write_3d_dataset_parallel(integer ndim, coord<hsize_t> fullshape, coord<hsize_t> chunkshape,
                                          coord<hsize_t> offset, const std::string& groupname,
                                          const Vec& datavec, hsize_t ndof, hsize_t idof)
{
  PetscPrintf(this->comm(), "Writing dataset: %s:%s\n", h5_filename.c_str(), groupname.c_str());
  // No need to set max size as want it to be same as given size.
  hid_t fspace_h = H5Screate_simple(ndim, fullshape.data(), nullptr);

  hid_t dset_h = H5Dcreate(_file_h, groupname.c_str(), H5T_NATIVE_DOUBLE, fspace_h, H5P_DEFAULT, H5P_DEFAULT,
                           H5P_DEFAULT);

  H5Sselect_hyperslab(fspace_h, H5S_SELECT_SET, offset.data(), nullptr, chunkshape.data(), nullptr);

  hid_t plist_h = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_h, H5FD_MPIO_COLLECTIVE);

  // Now pretend our 4D array is a 3D array and stride over last dimension to extract relevant DOF for
  // writing
  coord<hsize_t> dof_chunkshape = chunkshape;
  dof_chunkshape[2] *= ndof;
  coord<hsize_t> dof_offset = {0, 0, idof};
  coord<hsize_t> dof_stride = {1, 1, ndof};

  hid_t dspace_h = H5Screate_simple(ndim, chunkshape.data(), nullptr);
  H5Sselect_hyperslab(dspace_h, H5S_SELECT_SET, dof_offset.data(), dof_stride.data(), chunkshape.data(),
                      nullptr);

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
  PetscPrintf(this->comm(), "Writing dataset: %s:%s\n", h5_filename.c_str(), groupname.c_str());
  int rank;
  MPI_Comm_rank(this->comm(), &rank);
  // No need to set max size as want it to be same as given size.
  hid_t fspace_h = H5Screate_simple(1, &nval, nullptr);

  hid_t dset_h = H5Dcreate(_file_h, groupname.c_str(), H5T_NATIVE_DOUBLE, fspace_h, H5P_DEFAULT, H5P_DEFAULT,
                           H5P_DEFAULT);

  if (rank == 0)
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
  if (comm != this->comm())
  {
    std::ostringstream errss;
    errss << "Communicator mismatch between HDFWriter and provided image";
    throw InternalError(errss.str(), __FILE__, __LINE__);
  }

  std::string groupname = h5_groupname;
  if (groupname.empty())
  {
    groupname = "/registered";
  }

  int rank;
  MPI_Comm_rank(comm, &rank);

  write_3d_dataset_parallel(image.ndim(), image.shape().reverse(), image.local_shape().reverse(),
                            image.local_offset().reverse(), groupname, image.data_vector(), 1, 0);

  std::string filepath = h5_filename + ":" + h5_groupname;
  return filepath;
}

std::string HDFWriter::write_map(const MapBase& map)
{
  // Need to handle map layouts differently, so first, check if parallel
  if(map.parallel_layout())
  {
    return this->write_map_parallel(map);
  }
  // Otherwise write as for serial map
  return this->write_map_serial(map);

}

std::string HDFWriter::write_map_serial(const MapBase& map)
{
  int rank;
  MPI_Comm_rank(this->comm(), &rank);

  // Sanity check communicators
  MPI_Comm map_comm = map.comm();
  if (map_comm != MPI_COMM_SELF)
  {
    std::ostringstream errss;
    errss << "Map claims to be serial but map.comm() is not MPI_COMM_SELF";
    throw InternalError(errss.str(), __FILE__, __LINE__);
  }


  std::string groupname = h5_groupname;
  if (groupname.empty())
  {
    groupname = "/map";
  }


  hid_t mgroup_h = H5Gcreate(_file_h, groupname.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (mgroup_h < 0)
  {
    std::ostringstream errstr;
    errstr << "Failed to create group " << h5_groupname << ".";
    throw WriterError(errstr.str());
  }

  for (integer idx = 0; idx < map.ndim(); idx++)
  {
    // Get dataset local dimensions
    coord<hsize_t> offset = {0, 0, 0};
    coord<hsize_t> chunksize = {0, 0, 0};

    // Give rank 0 the whole array to write, easier than subdividing and no slower
    if(rank == 0)
    {
      offset = map.local_offset();
      chunksize = map.local_shape();
    }
    /*
    PetscSynchronizedPrintf(this->comm(), "Rank %i: ofs: %i %i %i, shp %i %i %i\n", rank,
    offset[0], offset[1], offset[2], chunksize[0], chunksize[1], chunksize[2]);
    PetscSynchronizedFlush(this->comm(), PETSC_STDOUT);
    */
    std::ostringstream dsetss;
    dsetss << groupname << "/" << BaseWriter::components()[idx];
    std::string dsetname = dsetss.str();

    write_3d_dataset_parallel(map.ndim(), map.shape(), chunksize, offset, dsetname,
                              map.displacement_vector(), map.ndim(), idx + 1);

    dsetss.clear();
    dsetss.str(std::string());
    dsetss << groupname << "/nodes_" << BaseWriter::components()[idx];
    dsetname = dsetss.str();

    write_1d_dataset_rank0(map.shape()[idx], dsetname, map.node_locs()[idx].data());
  }
  H5Gclose(mgroup_h);

  std::string filepath = h5_filename + ":" + h5_groupname;
  return filepath;
}

std::string HDFWriter::write_map_parallel(const MapBase& map)
{
  // Sanity check communicators
  MPI_Comm comm = map.comm();
  if (comm != this->comm())
  {
    std::ostringstream errss;
    errss << "Communicator mismatch between HDFWriter and provided map";
    throw InternalError(errss.str(), __FILE__, __LINE__);
  }

  std::string groupname = h5_groupname;
  if (groupname.empty())
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

  for (integer idx = 0; idx < map.ndim(); idx++)
  {
    // Get dataset local dimensions
    coord<hsize_t> offset = map.local_offset();
    coord<hsize_t> chunksize = map.local_shape();
    /*
    PetscSynchronizedPrintf(this->comm(), "Rank %i: ofs: %i %i %i, shp %i %i %i\n", rank,
    offset[0], offset[1], offset[2], chunksize[0], chunksize[1], chunksize[2]);
    PetscSynchronizedFlush(this->comm(), PETSC_STDOUT);
    */
    std::ostringstream dsetss;
    dsetss << groupname << "/" << BaseWriter::components()[idx];
    std::string dsetname = dsetss.str();

    write_3d_dataset_parallel(map.ndim(), map.shape(), chunksize, offset, dsetname,
                              map.displacement_vector(), map.ndim(), idx + 1);

    dsetss.clear();
    dsetss.str(std::string());
    dsetss << groupname << "/nodes_" << BaseWriter::components()[idx];
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
  H5Pset_fapl_mpio(file_props, this->comm(), MPI_INFO_NULL);

  // First supress HDF5 error printing
  herr_t (*old_err)(void*);
  void* old_err_data;
  H5Eget_auto1(&old_err, &old_err_data);
  H5Eset_auto1(nullptr, nullptr);

  // Open if available, will get file_h > 0 if successful
  if (BaseWriter::check_truncated(h5_filename))
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
