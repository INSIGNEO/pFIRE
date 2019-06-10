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

#ifndef MASK_HPP
#define MASK_HPP

#include "types.hpp"

#include "imagebase.hpp"

class Mask : public ImageBase
{
public:

  static std::shared_ptr<Mask> create_filled_mask(const ImageBase& pattern);

protected:

  Mask(const ImageBase& pattern)
    : ImageBase(pattern) 
  {
  }
};

#endif //MASK_HPP
