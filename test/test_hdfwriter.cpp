#define BOOST_TEST_MODULE hdfwriter
#include "test_common.hpp"

#include <numeric>

#include <petscdmda.h>
#include <petscerror.h>

#include <hdf5.h>

#include "fd_routines.hpp"
#include "hdfwriter.hpp"
#include "image.hpp"
#include "indexing.hpp"
#include "types.hpp"

struct im {
  im() : image(imgshape, PETSC_COMM_WORLD) { BOOST_TEST_MESSAGE("setup image"); }

  intvector imgshape = {17, 11, 14};
  Image image;
  std::string filename = "testdata.h5";
  std::string groupname = "/testimage";
  std::string writercompound = filename + ":" + groupname;
};

void write_image_rmaj_indices(Image& image)
{
  PetscErrorCode perr;

  // Access data via dmda and insert values
  // in this case insert row major cell index as value
  auto vtx_lo = image.mpi_get_offset<integer>();
  auto vtx_hi = image.mpi_get_chunksize<integer>();
  std::transform(vtx_lo.cbegin(), vtx_lo.cend(), vtx_hi.begin(), vtx_hi.begin(), std::plus<>());

  integer ymult = image.shape()[2];
  integer xmult = image.shape()[1];

  floating*** ptr;
  perr = DMDAVecGetArray(*image.dmda(), *image.global_vec(), &ptr);
  CHKERRXX(perr);
  for (integer xx = vtx_lo[0]; xx < vtx_hi[0]; xx++)
  {
    for (integer yy = vtx_lo[1]; yy < vtx_hi[1]; yy++)
    {
      for (integer zz = vtx_lo[2]; zz < vtx_hi[2]; zz++)
      {
        ptr[zz][yy][xx] = (xx * xmult + yy) * ymult + zz;
      }
    }
  }
  perr = DMDAVecRestoreArray(*image.dmda(), *image.global_vec(), &ptr);
  CHKERRXX(perr);
}

BOOST_FIXTURE_TEST_SUITE(hdfwriter, im)

BOOST_AUTO_TEST_CASE(check_image_ordering)
{
  // Write row-major indices as image values
  write_image_rmaj_indices(image);

  { // enclose for 'GC' purposes
    HDFWriter wtr(writercompound, PETSC_COMM_WORLD);
    wtr.write_image(image);
  }

  hid_t file_h, dset_h;
  herr_t status;
  std::vector<floating> hdf_data(image.size(), 0.0);
  std::vector<floating> cmp_data(image.size(), 0.0);

  // Create comparison data using iota
  std::iota(cmp_data.begin(), cmp_data.end(), 0);

  file_h = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  dset_h = H5Dopen(file_h, groupname.c_str(), H5P_DEFAULT);

  status = H5Dread(dset_h, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, hdf_data.data());
  if(status < 0)
  {
    throw std::runtime_error("Failed to read dataset");
  }

  H5Dclose(dset_h);
  H5Fclose(file_h);

  BOOST_REQUIRE(std::equal(hdf_data.cbegin(), hdf_data.cend(), cmp_data.cbegin(), cmp_data.cend()));
}

BOOST_AUTO_TEST_SUITE_END()
