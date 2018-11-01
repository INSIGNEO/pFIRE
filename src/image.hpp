#ifndef IMAGE_HPP
#define IMAGE_HPP

#include<mpi.h>

#include "types.hpp"

class Image{

public:

  explicit Image(const intvector &shape, MPI_Comm comm=PETSC_COMM_WORLD);

  Vec_unique gradient(integer dim);

  std::unique_ptr<Image> duplicate() const;
  std::unique_ptr<Image> copy() const;

  // Allow RO access to member variables
  // Note that datavec and dmda remain mutable in this way
  MPI_Comm comm() const {return m_comm;}
  uinteger ndim() const {return m_ndim;}
  const intvector& shape() const { return m_shape;}
  integer size() const { return m_shape[0] * m_shape[1] * m_shape[2];}
  std::shared_ptr<const Vec> global_vec() const{ return m_globalvec;}
  std::shared_ptr<const Vec> local_vec() const{ return m_localvec;}
  std::shared_ptr<DM> dmda() const{ return m_dmda;}

  floating normalize();

//  void set_mask(std::shared_ptr<Mask>);
//  void masked_normalize(const Mask& mask);

  void update_local_from_global();

  static std::unique_ptr<Image> load_file(const std::string &filename,
                                          const Image *existing=nullptr,
                                          MPI_Comm comm=PETSC_COMM_WORLD);

  void save_OIIO(std::string filename);


protected:

  explicit Image(const Image& image);
  Image& operator=(const Image& image);

  MPI_Comm m_comm;
  uinteger m_ndim;
  intvector m_shape;
  Vec_shared m_localvec, m_globalvec;
  DM_shared m_dmda;
//  std::shared_ptr<Mask> mask;

  void initialize_dmda(); 
  void initialize_vectors(); 

  Vec_unique scatter_to_zero(Vec &vec);
};

#endif
