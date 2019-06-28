#include "image.hpp"

#include "fd_routines.hpp"
#include "indexing.hpp"
#include "infix_iterator.hpp"
#include "mapbase.hpp"
#include "mpi_routines.hpp"

/*! Calculate gradient vectors for T matrix
 *
 */
std::vector<Vec_unique> calculate_tmatrix_gradients(const Image& fixed, const Image& moved, integer ndof)
{
  std::vector<Vec_unique> gradients(ndof + 1);

  for (auto& gradvec : gradients)
  {
    gradvec = create_unique_vec();
    PetscErrorCode perr = VecDuplicate(fixed.global_vector(), gradvec.get());
    CHKERRXX(perr);
  }

  fixed.update_local_vector();
  moved.update_local_vector();

  // Calculate intensity term
  PetscErrorCode perr = VecSet(*gradients[0], 1.);
  CHKERRXX(perr);
  // Calculate 1-0.5(f+m)
  const floating multiplier = -0.5; // Average contributions from both images
  perr = VecAXPBYPCZ(*gradients[0], multiplier, multiplier, 1, fixed.global_vector(), moved.global_vector());
  CHKERRXX(perr);

  Vec_unique local_diff_vec = create_unique_vec();
  perr = DMCreateLocalVector(fixed.dmda(), local_diff_vec.get());
  CHKERRXX(perr);
  perr = DMGlobalToLocalBegin(fixed.dmda(), *gradients[0], INSERT_VALUES, *local_diff_vec);
  CHKERRXX(perr);
  perr = DMGlobalToLocalEnd(fixed.dmda(), *gradients[0], INSERT_VALUES, *local_diff_vec);
  CHKERRXX(perr);

  floating gradsum, gradmin, gradmax;
  VecSum(*gradients[0], &gradsum);
  VecMin(*gradients[0], nullptr, &gradmin);
  VecMax(*gradients[0], nullptr, &gradmax);

#ifdef VERBOSE_DEBUG
  PetscPrintf(fixed.comm(), "Gradients [%i] %f, %f, %f, (sum, min, max)\n", 0, gradsum, gradmin, gradmax);
#endif //VERBOSE_DEBUG

  for (integer idx = 0; idx < fixed.ndim(); idx++)
  {
    gradient_existing(fixed.dmda(), *local_diff_vec, *gradients[idx + 1], idx);
    floating gradsum, gradmin, gradmax;
    VecSum(*gradients[idx + 1], &gradsum);
    VecMin(*gradients[idx + 1], nullptr, &gradmin);
    VecMax(*gradients[idx + 1], nullptr, &gradmax);
#ifdef VERBOSE_DEBUG
    PetscPrintf(fixed.comm(), "Gradients [%i] %f, %f, %f, (sum, min, max)\n", idx+1, gradsum, gradmin, gradmax);
#endif //VERBOSE_DEBUG
  }

  return gradients;
}

std::shared_ptr<Image> Image::warp(const MapBase& map) const
{
  std::shared_ptr<Image> new_image = GridVariable::duplicate(*this);

  map.update_local_vector();

  // No guarantee of where any given source pixel is located, so need to locate rank of each each source and
  // add that plus target pointer to suitable lists
  floatptrvector2d target_ptrs(this->commsize());
  floatcoordvector2d global_source_locs(this->commsize());
  intcoordvector local_target_locs;
  floatcoordvector local_source_locs;

  auto owned_range = this->owned_range();
  auto& owned_lo = owned_range.first;
  auto& owned_hi = owned_range.second;

  intcoord curr_loc;
  auto& xx = curr_loc[0];
  auto& yy = curr_loc[1];
  auto& zz = curr_loc[2];

  floating*** target_data;
  PetscErrorCode perr = DMDAVecGetArray(new_image->dmda(), new_image->global_vector(), &target_data);
  CHKERRXX(perr);
  // Iterate over all owned pixels and fill request structures
  for (zz = owned_lo[2]; zz < owned_hi[2]; zz++)
  {
    for (yy = owned_lo[1]; yy < owned_hi[1]; yy++)
    {
      for (xx = owned_lo[0]; xx < owned_hi[0]; xx++)
      {
        local_target_locs.push_back(curr_loc);
      }
    }
  }

  local_source_locs = map.map_local_coordinates(local_target_locs);

  for (size_t iidx = 0; iidx < local_source_locs.size(); iidx++)
  {
    auto& source_loc = local_source_locs[iidx];
    auto& target_loc = local_target_locs[iidx];
    integer loc_rank = this->get_rank_of_loc(source_loc);
    global_source_locs[loc_rank].push_back(source_loc);
    target_ptrs[loc_rank].push_back(&target_data[target_loc[2]][target_loc[1]][target_loc[0]]);
  }

  // scatter as needed
  if (this->commsize() > 1)
  {
    global_source_locs = p2p_vecscatter(global_source_locs, this->comm());
  }

  floatvector2d interpolated_data;
  for (const auto& rank_sources : global_source_locs)
  {
    interpolated_data.push_back(this->interpolate_pixel_data(rank_sources));
  }

  // scatter back as needed
  if (this->commsize() > 1)
  {
    interpolated_data = p2p_vecscatter(interpolated_data, this->comm());
  }

  for (integer irank = 0; irank < this->commsize(); irank++)
  {
    auto& rank_pixel_data = interpolated_data[irank];
    auto& rank_target_ptrs = target_ptrs[irank];
    for (size_t iidx = 0; iidx < rank_pixel_data.size(); iidx++)
    {
      *rank_target_ptrs[iidx] = rank_pixel_data[iidx];
    }
  }

  DMDAVecRestoreArray(new_image->dmda(), new_image->global_vector(), &target_data);
  CHKERRXX(perr);

  return new_image;
}

floatvector Image::interpolate_pixel_data(const floatcoordvector& source_locations) const
{
  // Create and reserve output vector
  floatvector pixel_data;
  pixel_data.reserve(source_locations.size());

  const intcoordvector* offsets;
  if (this->ndim() == 2)
  {
    offsets = &MapBase::node_offset_list_2d;
  }
  else
  {
    offsets = &MapBase::node_offset_list_3d;
  }

  floating*** image_data;
  PetscErrorCode perr = DMDAVecGetArrayRead(this->dmda(), this->local_vector(), &image_data);
  CHKERRXX(perr);
  // Iterate over all coords to interpolate
  for (const auto& loc : source_locations)
  {
    floating pixel_value = 0;
    intcoord floor_loc = floor(loc);
    for (const auto& offset : *offsets)
    {
      intcoord offset_loc = floor_loc + offset;
      floating coeff = calculate_unscaled_basis_coefficient(offset_loc.cbegin(), offset_loc.cend(),
                                                            loc.cbegin());
      if (coeff > 0)
      {
        pixel_value += coeff * image_data[offset_loc[2]][offset_loc[1]][offset_loc[0]];
      }
    }
    pixel_data.push_back(pixel_value);
  }
  perr = DMDAVecRestoreArrayRead(this->dmda(), this->local_vector(), &image_data);
  CHKERRXX(perr);

  return pixel_data;
}
