#include<numeric>
#include<cstdio>
#include<memory>

#include<petscsys.h>
#include<petscdmda.h>

#include "types.hpp"
#include "image.hpp"
#include "indexing.hpp"
#include "map.hpp"

void mainflow();

void mainflow(std::string, std::string);

int main(int argc, char **argv)
{
  PetscErrorCode perr = PetscInitialize(&argc, &argv, NULL, NULL);CHKERRABORT(PETSC_COMM_WORLD, perr);

  mainflow(argv[1], argv[2]);

  PetscFinalize();

  return 0;
}

void mainflow(std::string inputimgpath, std::string warpedpath)
{

  std::unique_ptr<Image> inputimg = Image::create_from_image(inputimgpath);

  floatvector nodespacing = {50,50};

  std::unique_ptr<Map> map = std::make_unique<Map>(*inputimg, nodespacing);

  WorkSpace wksp(*inputimg, *map);

  integer M;
  VecGetSize(*map->m_displacements, &M);
  PetscPrintf(PETSC_COMM_WORLD, "Displacement size: %i\n", M);
  intvector idxn(9, 0);
  std::iota(idxn.begin(), idxn.end(), 0);
  floatvector data = {-50, -50, -50,
                      0, 0, 0,
                      50, 50, 50};
//  floatvector data(9, 10.);
  PetscErrorCode perr;
  perr = VecSetValues(*map->m_displacements, 9, idxn.data(), data.data(), INSERT_VALUES);CHKERRXX(perr);
//  perr = VecSet(*map->m_displacements, 10);CHKERRXX(perr);
  perr = VecAssemblyBegin(*map->m_displacements);CHKERRXX(perr);
  perr = VecAssemblyEnd(*map->m_displacements);CHKERRXX(perr);

  std::unique_ptr<Image> warped = map->warp(*inputimg, wksp);
  warped->save_OIIO(warpedpath);
}
