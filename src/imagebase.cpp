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

#include "imagebase.hpp"

#include "petscvec.h"

ImageBase::ImageBase(const intcoord& shape, const MPI_Comm& comm)
    : GridVariable(shape, 1, comm)
{
}


Vec_unique ImageBase::difference_global(const ImageBase& other) const
{
  Vec_unique image_difference = create_unique_vec();
  PetscErrorCode perr = VecDuplicate(this->global_vector(), image_difference.get());
  CHKERRXX(perr);
  perr = VecCopy(this->global_vector(), *image_difference);
  CHKERRXX(perr);
  perr = VecAXPY(*image_difference, -1, other.global_vector());

  return image_difference;
}

/*! Normalize the image to maximum value and return scale factor
 *
 */
floating ImageBase::normalize()
{
  floating norm;
  //PetscErrorCode perr = VecNorm(this->global_vector(), NORM_1, &norm);
  PetscErrorCode perr = VecMax(this->global_vector(), nullptr, &norm);
  CHKERRXX(perr);
  //If VecMax is zero probably an empty image, but still don't want NaNs
  //norm = norm == 0 ? 1 : this->num_total_nodes()/norm;
  norm = norm == 0 ? 1 : 1/norm;
  perr = VecScale(this->global_vector(), norm);
  CHKERRXX(perr);
  this->update_local_vector();
  return norm;
}

floating ImageBase::mutual_information(const ImageBase &other)
{
  integer mi_resolution = 50;
  // Mutual information is given by integrating P(X, Y)log(P(X,Y)/(P(X)P(Y)) over all X and all Y,
  // In this case the probabilities can be calculated from a 2D histogram of X against Y.

  // Could also just find 2D and perform summations after comm, but for now do it this way
  floatvector xhist(mi_resolution, 0);
  floatvector yhist(mi_resolution, 0);
  floatvector2d xyhist(mi_resolution, floatvector(mi_resolution, 0));

  integer img_localsize;
  PetscErrorCode perr = VecGetLocalSize(this->global_vector(), &img_localsize);
  CHKERRXX(perr);

  floating max1, max2;
  perr = VecMax(this->global_vector(), nullptr, &max1);CHKERRXX(perr);
  perr = VecMax(other.global_vector(), nullptr, &max2);CHKERRXX(perr);
  floating max = std::max(max1, max2);

  // Only need RO data from PETSc vecs
  floating const *x_data, *y_data;
  perr = VecGetArrayRead(this->global_vector(), &x_data);
  CHKERRXX(perr);
  perr = VecGetArrayRead(other.global_vector(), &y_data);
  CHKERRXX(perr);

  for (integer idx = 0; idx < img_localsize; idx++)
  {
    integer x = std::lround(x_data[idx]/max * static_cast<double>(mi_resolution-1));
    integer y = std::lround(y_data[idx]/max * static_cast<double>(mi_resolution-1));
    xhist[x] += 1;
    yhist[y] += 1;
    xyhist[x][y] += 1;
  }

  int rank;
  MPI_Comm_rank(this->comm(), &rank);
  if(rank == 0)
  {
    MPI_Reduce(MPI_IN_PLACE, xhist.data(), xhist.size(), MPIU_SCALAR, MPI_SUM, 0, this->comm());
    MPI_Reduce(MPI_IN_PLACE, yhist.data(), yhist.size(), MPIU_SCALAR, MPI_SUM, 0, this->comm());
    for(auto &vec1d: xyhist)
    {
      MPI_Reduce(MPI_IN_PLACE, vec1d.data(), vec1d.size(), MPIU_SCALAR, MPI_SUM, 0, this->comm());
    }
  } 
  else
  {
    MPI_Reduce(xhist.data(), xhist.data(), xhist.size(), MPIU_SCALAR, MPI_SUM, 0, this->comm());
    MPI_Reduce(yhist.data(), yhist.data(), yhist.size(), MPIU_SCALAR, MPI_SUM, 0, this->comm());
    for(auto &vec1d: xyhist)
    {
      MPI_Reduce(vec1d.data(), vec1d.data(), vec1d.size(), MPIU_SCALAR, MPI_SUM, 0, this->comm());
    }
  }

  // Need all probability distributions to sum to 1.0
  integer pix_tot = this->num_total_nodes();
  std::transform(xhist.cbegin(), xhist.cend(), xhist.begin(),
                 [pix_tot](floating x) -> floating {return x/pix_tot;});
  std::transform(xhist.cbegin(), xhist.cend(), xhist.begin(),
                 [pix_tot](floating y) -> floating {return y/pix_tot;});
  for(auto &vec1d: xyhist)
  {
    std::transform(vec1d.cbegin(), vec1d.cend(), vec1d.begin(),
                   [pix_tot](floating y) -> floating {return y/pix_tot;});
  }

  // Use Kahan summation to try to minimize error
  floating mi_total(0);
  floating mi_err(0);
  if(rank == 0)
  {
    for(size_t ix(0); ix < xhist.size(); ix++)
    {
      for(size_t iy(0); iy < yhist.size(); iy++)
      {
        if(xyhist[ix][iy] > 0)
        {
          floating mi_part = xyhist[ix][iy] * std::log(xyhist[ix][iy]/(xhist[ix]*yhist[iy]));
          mi_part -= mi_err;
          floating tmp = mi_total + mi_part;
          mi_err = (tmp - mi_total) - mi_part;
          mi_total = tmp;
        }
      }
    }
  }

  MPI_Bcast(&mi_total, 1, MPIU_SCALAR, 0, this->comm());

  return mi_total;
}
