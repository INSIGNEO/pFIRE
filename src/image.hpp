#ifndef IMAGE_HPP
#define IMAGE_HPP

#include<algorithm>
#include<memory>
#include<string>
#include<vector>

#include<mpi.h>

#include<petscdmda.h>
#include<petscvec.h>

#include "types.hpp"
#include "map.hpp"

#include "fd_routines.hpp"
#include "iterator_routines.hpp"

class Map;

class Image{

public:

  explicit Image(std::string filename, MPI_Comm comm = PETSC_COMM_WORLD);
  explicit Image(intvector shape, MPI_Comm comm = PETSC_COMM_WORLD);
  explicit Image(const Image& image);

  Vec_unique gradient(integer dim);

  // Allow RO access to member variables
  // Note that datavec and dmda remain mutable in this way
  MPI_Comm comm() const {return m_comm;}
  integer ndim() const {return m_ndim;}
  const intvector& shape() const { return m_shape;}
  integer size() const { return m_shape[0] * m_shape[1] * m_shape[2];}
  std::shared_ptr<const Vec> global_vec() const{ return m_globalvec;}
  std::shared_ptr<const Vec> local_vec() const{ return m_localvec;}
  std::shared_ptr<DM> dmda() const{ return m_dmda;}

  friend void create_t_matrix(Image &fixed, Image &moved, Map &map);

protected:

  integer m_ndim;
  MPI_Comm m_comm;
  intvector m_shape;
  Vec_shared m_localvec, m_globalvec;
  std::shared_ptr<DM> m_dmda;

  void initialize_dmda(); 
  void initialize_vectors(); 
  void update_local_from_global();
};

#endif
