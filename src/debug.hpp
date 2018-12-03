#ifndef DEBUG_HPP
#define DEBUG_HPP

#include<string>

#include<petscmat.h>

#include "types.hpp"

void matrix_dbg_print(const MPI_Comm &comm, const Mat &mat, const std::string &name);

#endif //DEBUG_HPP
