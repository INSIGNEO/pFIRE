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

#ifndef PETSC_DEBUG_HPP
#define PETSC_DEBUG_HPP

#include <cstdint>

#ifdef DEBUG_OBJECTS
#define debug_creation(obj, name) do_debug_creation(obj, name, __FILE__, __LINE__)
#define debug_deletion(obj) do_debug_deletion(obj, __FILE__, __LINE__)
#else
#define debug_creation(obj, name) do {} while(0)
#define debug_deletion(obj) do {} while(0)
#endif

template <typename PetscXX>
void do_debug_creation(PetscXX &uncast_obj, std::string name, std::string filename,
                       uint32_t linenum)
{
  PetscObject obj = (PetscObject) uncast_obj;
  MPI_Comm comm;
  PetscObjectSetName(obj, name.c_str());
  PetscObjectGetComm(obj, &comm);
  PetscPrintf(comm, "Created %s (%s:%i)\n", name.c_str(), filename.c_str(), linenum);
}

template <typename PetscXX>
void do_debug_deletion(PetscXX &uncast_obj, std::string filename, uint32_t linenum)
{
  if(uncast_obj == nullptr)
  {
    return;
  }
  PetscObject obj = (PetscObject) uncast_obj;
  MPI_Comm comm;
  char const *name;
  PetscObjectGetName(obj, &name);
  PetscObjectGetComm(obj, &comm);
  PetscPrintf(comm, "Destroyed %s (%s:%i)\n", name, filename.c_str(), linenum);
}

#endif // PETSC_DEBUG_HPP
