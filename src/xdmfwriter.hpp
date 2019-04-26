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

#ifndef XDMFWRITER_HPP
#define XDMFWRITER_HPP

#include <string>

#include <hdf5.h>
#include <mpi.h>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include "hdfwriter.hpp"
#include "types.hpp"

namespace pt = boost::property_tree;

class XDMFWriter: public HDFWriter {
public:
  XDMFWriter(std::string filespec, const MPI_Comm& comm);
  ~XDMFWriter();

  std::string write_image(const Image& image);
  std::string write_map(const Map& map);

  static const std::string writer_name;
  static const std::vector<std::string> extensions;

protected:
  int rank;
  std::string xdmf_filename;
  std::string xdmf_groupname;
  pt::ptree xdmf_tree;

  static std::string h5name_from_xdmfname(const std::string &filename);

private:
  void prepopulate_xdmf_ptree();
};

#endif // XDMFWRITER_HPP
