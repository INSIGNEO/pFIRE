#ifndef BASELOADER_HPP
#define BASELOADER_HPP

#include <functional>

#include "types.hpp"

class BaseLoader{
public:

  using loader_creator = std::function<BaseLoader_unique(
      const std::string& filename, MPI_Comm comm)>;

  using loader_map = std::map<std::string, loader_creator>;

  BaseLoader(const std::string& path, MPI_Comm comm=PETSC_COMM_WORLD)
  : _comm(comm), _path(path), _shape(intvector(3, 0))
  {};

  virtual ~BaseLoader() = default;

  const intvector &shape() const{ return this->_shape;};
  static const loader_map &loaders(){ return *_loaders;};

  virtual void copy_scaled_chunk(floating ***data, const intvector &size, 
                                 const intvector &offset) const = 0;

  static bool register_loader(const std::string &name, loader_creator loader);

  static BaseLoader_unique find_loader(const std::string &name, MPI_Comm comm=PETSC_COMM_WORLD);

protected:

  MPI_Comm _comm;
  std::string _path;
  intvector _shape;

private:

  static std::unique_ptr<loader_map> _loaders;
  
};

template<class T, class U>
void norm_convert(U *outptr, const T *inptr, uinteger size, double min, double max)
{
  U fscale = (U)(max - min);
  U fofs = (U)min;

  for(size_t idx = 0; idx < size; idx++)
  {
    outptr[idx] = (inptr[idx] - fofs)/fscale;
  }
}

#endif //BASELOADER_HPP
