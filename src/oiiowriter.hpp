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

#ifndef OIIOWRITER_HPP
#define OIIOWRITER_HPP

#include <string>

#include "basewriter.hpp"

class OIIOWriter : public BaseWriter {
  public:
  OIIOWriter(std::string filename, const MPI_Comm& comm = PETSC_COMM_WORLD);
  ~OIIOWriter() = default;

  void write_image(const Image& image);
  
  void write_map(const Map& map);
};

#endif // OIIOWRITER_HPP
