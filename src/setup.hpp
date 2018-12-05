#ifndef SETUP_HPP
#define SETUP_HPP

#include<vector>
#include<string>

std::string get_invocation_name(const std::string &argzero);

void pfire_setup(const std::vector<std::string> &petsc_args);

void check_and_warn_odd_comm();

void pfire_teardown();

#endif //SETUP_HPP
