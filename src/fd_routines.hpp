#ifndef FD_ROUTINES_HPP
#define FD_ROUTINES_HPP

#include<exception>
#include<memory>

#include<petscdmda.h>
#include<petscvec.h>

#include "types.hpp"
#include "petsc_helpers.hpp"

namespace fd
{

Vec_unique gradient_to_global_unique(const DM &dmda, const Vec &localvec, integer dim);

}

#endif

