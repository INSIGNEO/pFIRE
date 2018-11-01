#include "baseloader.hpp"

std::unique_ptr<BaseLoader::loader_map> BaseLoader::_loaders = std::make_unique<BaseLoader::loader_map>();

bool BaseLoader::register_loader(const std::string &name, loader_creator loader)
{
  if(!BaseLoader::_loaders)
  {
    _loaders = std::make_unique<std::map<std::string, BaseLoader::loader_creator> >(); 
  }
  if (auto it = _loaders->find(name); it == _loaders->cend())
  {
    std::cout << "Registered loader \"" << name << "\".\n";
    (*_loaders)[name] = loader;
    return true;
  }
  std::ostringstream err;
  err << "Loader named " << name << " already registered.";
  throw std::runtime_error(err.str());
}

BaseLoader_unique BaseLoader::find_loader(const std::string &path, MPI_Comm comm)
{
  //Do it EAFP style...
  for(const auto &it: *_loaders)
  {
    try
    {
      BaseLoader_unique loader = it.second(path, comm);
      return loader;
    }
    catch (...)
    {
      continue;
    }
  }
  std::ostringstream err;
  err << "Failed to find suitable loader for \"" << path << "\".";
  throw std::runtime_error(err.str());
}
