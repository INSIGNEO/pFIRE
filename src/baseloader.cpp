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

#include "baseloader.hpp"

#include "exceptions.hpp"

std::unique_ptr<BaseLoader::loader_map> BaseLoader::_loaders =
    std::make_unique<BaseLoader::loader_map>();

bool BaseLoader::register_loader(const std::string &name, loader_creator loader)
{
  if (!BaseLoader::_loaders)
  {
    _loaders = std::make_unique<std::map<std::string, BaseLoader::loader_creator>>();
  }
  if (auto it = _loaders->find(name); it == _loaders->cend())
  {
    (*_loaders)[name] = loader;
    return true;
  }
  std::ostringstream err;
  err << "Loader named " << name << " already registered.";
  throw InternalError(err.str(), __FILE__, __LINE__);
}

BaseLoader_unique BaseLoader::find_loader(const std::string &path, MPI_Comm comm)
{
  // Do it EAFP style...
  for (const auto &it : *_loaders)
  {
    try
    {
      BaseLoader_unique loader = it.second(path, comm);
      return loader;
    }
    catch (InvalidLoaderError&)
    {
      continue;
    }
  }
  std::ostringstream err;
  err << "Failed to find suitable loader for \"" << path << "\".";
  throw InvalidLoaderError(err.str());
}
