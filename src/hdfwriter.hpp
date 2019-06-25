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

#ifndef HDFWRITER_HPP
#define HDFWRITER_HPP

#include <string>

#include <hdf5.h>
#include <mpi.h>

#include "basewriter.hpp"
#include "types.hpp"

class HDFWriter: public BaseWriter {
public:
  HDFWriter(const std::string& filename, const MPI_Comm& comm);
  ~HDFWriter();

  std::string write_image(const Image& image);
  std::string write_map(const MapBase& map);

  static const std::string writer_name;
  static const std::vector<std::string> extensions;

// Make sure we always choose the correct datatype to read back to
  static hid_t get_hdf5_petsc_scalar();

protected:
  std::string h5_filename;
  std::string h5_groupname;

private:
  hid_t _file_h;

  void open_or_create_h5();

  void write_1d_dataset_rank0(hsize_t nval, const std::string& groupname, const floating* databuf);

  void write_3d_dataset_parallel(integer ndim, coord<hsize_t> fullshape, coord<hsize_t> chunkshape,
                                 coord<hsize_t> offset, const std::string& groupname, const Vec& datavec,
                                 hsize_t ndof, hsize_t idof);

  std::string write_map_parallel(const MapBase& map);
  std::string write_map_serial(const MapBase& map);
};

#endif // HDFWRITER_HPP
