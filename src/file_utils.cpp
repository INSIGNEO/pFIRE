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

#include "file_utils.hpp"

#include <boost/filesystem.hpp>

#include "exceptions.hpp"

namespace bf = boost::filesystem;

void throw_if_nonexistent(const std::string& filepath)
{
  if (!bf::exists(filepath))
  {
    throw FileNotFoundError(filepath);
  }
}
