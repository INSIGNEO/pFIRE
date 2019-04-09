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

#ifndef TYPEDEFS_HPP
#define TYPEDEFS_HPP

#include <array>
#include <iostream>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include <petscdm.h>
#include <petscksp.h>
#include <petscmat.h>
#include <petscsys.h>
#include <petscvec.h>

#include "petsc_debug.hpp"

// Forward Defs
//
class Image;
class WorkSpace;
class Map;
class Elastic;
class BaseLoader;
class BaseWriter;

// Useful typedefs
using integer = PetscInt;
using floating = PetscScalar;

using intloc = std::array<integer, 3>;
using floatloc = std::array<floating, 3>;

using intvector = std::vector<integer>;
using intvector2d = std::vector<intvector>;
using intvector3d = std::vector<intvector2d>;

using uinteger = unsigned integer;
using uintvector = std::vector<uinteger>;

using floatvector = std::vector<floating>;
using floatvector2d = std::vector<floatvector>;

using stringset = std::set<std::string>;

using BaseLoader_unique = std::unique_ptr<BaseLoader>;
using BaseWriter_unique = std::unique_ptr<BaseWriter>;

//// Self-destructing PETSc objects
// Can apply identical treatment to many PETSc types

//// Vec
// typedef and helpers for unique_ptr
struct VecDeleter {
  void operator()(Vec* v) const
  {
    debug_deletion(*v);
    VecDestroy(v);
    delete v;
  }
};

using Vec_unique = std::unique_ptr<Vec, VecDeleter>;
inline Vec_unique create_unique_vec()
{
  Vec_unique v = Vec_unique(new Vec);
  *v = nullptr;
  return v;
}

// typedef and helper for shared_ptr
using Vec_shared = std::shared_ptr<Vec>;
inline Vec_shared create_shared_vec()
{
  Vec_shared v = Vec_shared(new Vec, VecDeleter());
  *v = nullptr;
  return v;
}

//// Mat
// typedef and helpers for unique_ptr
struct MatDeleter {
  void operator()(Mat* v) const
  {
    debug_deletion(*v);
    MatDestroy(v);
    delete v;
  }
};

using Mat_unique = std::unique_ptr<Mat, MatDeleter>;
inline Mat_unique create_unique_mat()
{
  Mat_unique v = Mat_unique(new Mat);
  *v = nullptr;
  return v;
}

// typedef and helper for shared_ptr
using Mat_shared = std::shared_ptr<Mat>;
inline Mat_shared create_shared_mat()
{
  Mat_shared v = Mat_shared(new Mat, MatDeleter());
  *v = nullptr;
  return v;
}

//// DM
// typedef and helpers for unique_ptr
struct DMDeleter {
  void operator()(DM* v) const
  {
    DMDestroy(v);
    delete v;
  }
};

using DM_unique = std::unique_ptr<DM, DMDeleter>;
inline DM_unique create_unique_dm()
{
  DM_unique v = DM_unique(new DM);
  *v = nullptr;
  return v;
}

// typedef and helper for shared_ptr
using DM_shared = std::shared_ptr<DM>;
inline DM_shared create_shared_dm()
{
  DM_shared v = DM_shared(new DM, DMDeleter());
  *v = nullptr;
  return v;
}

//// IS
// typedef and helpers for unique_ptr
struct ISDeleter {
  void operator()(IS* v) const
  {
    ISDestroy(v);
    delete v;
  }
};

using IS_unique = std::unique_ptr<IS, ISDeleter>;
inline IS_unique create_unique_is()
{
  IS_unique v = IS_unique(new IS);
  *v = nullptr;
  return v;
}

// typedef and helper for shared_ptr
using IS_shared = std::shared_ptr<IS>;
inline IS_shared create_shared_is()
{
  IS_shared v = IS_shared(new IS, ISDeleter());
  *v = nullptr;
  return v;
}

//// VecScatter
// typedef and helpers for unique_ptr
struct VecScatterDeleter {
  void operator()(VecScatter* v) const
  {
    VecScatterDestroy(v);
    delete v;
  }
};

using VecScatter_unique = std::unique_ptr<VecScatter, VecScatterDeleter>;
inline VecScatter_unique create_unique_vecscatter()
{
  VecScatter_unique v = VecScatter_unique(new VecScatter);
  *v = nullptr;
  return v;
}

// typedef and helper for shared_ptr
using VecScatter_shared = std::shared_ptr<VecScatter>;
inline VecScatter_shared create_shared_vecscatter()
{
  VecScatter_shared v = VecScatter_shared(new VecScatter, VecScatterDeleter());
  *v = nullptr;
  return v;
}

//// KSP
// typedef and helpers for unique_ptr
struct KSPDeleter {
  void operator()(KSP* v) const
  {
    KSPDestroy(v);
    delete v;
  }
};

using KSP_unique = std::unique_ptr<KSP, KSPDeleter>;
inline KSP_unique create_unique_ksp()
{
  KSP_unique v = KSP_unique(new KSP);
  *v = nullptr;
  return v;
}

// typedef and helper for shared_ptr
using KSP_shared = std::shared_ptr<KSP>;
inline KSP_shared create_shared_ksp()
{
  KSP_shared v = KSP_shared(new KSP, KSPDeleter());
  *v = nullptr;
  return v;
}

#endif
