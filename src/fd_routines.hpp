#ifndef FD_ROUTINES_HPP
#define FD_ROUTINES_HPP

#include <exception>
#include <memory>

#include <petscdmda.h>
#include <petscvec.h>

#include "petsc_helpers.hpp"
#include "types.hpp"

namespace fd
{
Vec_unique gradient_to_global_unique(const DM &dmda, const Vec &localvec, integer dim);

void gradient_existing(const DM &dmda, const Vec &src, Vec &tgt, integer dim);

} // namespace fd

#endif
