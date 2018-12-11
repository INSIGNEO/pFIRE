#ifndef UTILS_HPP
#define UTILS_HPP

#include <numeric>

floating mat_size(integer ndim, integer nvox)
{
  floating nd = std::pow(2, ndim);
  std::cout << "intsize" << sizeof(PetscInt) << "\n";
  // matrix index size seems to be 8 even if PetscInt is 4 byte...
  integer petscintsize = 8;
  return nvox * (nd * sizeof(floating) + (nd + 1) * petscintsize);
}

void explain_memory(const intvector &image_size, const intvector &map_size)
{
  if (image_size.size() < 2 || image_size.size() > 3)
  {
    throw std::runtime_error("Image should be 2D or 3D");
  }

  if (map_size.size() != image_size.size())
  {
    throw std::runtime_error("Map and image dimensions should match");
  }

  integer ndim = image_size.size();
  floating nvox = std::accumulate(image_size.begin(), image_size.end(), 1, std::multiplies<>());
  floating nmapnod = std::accumulate(map_size.begin(), map_size.end(), 1, std::multiplies<>());

  integer floatbytes = sizeof(floating);

  integer nimage = 9;
  integer nwarp_mat = 2;
  integer nmap = 2;
  integer nreg_mat = 2;

  floating warp_mat = mat_size(ndim, nvox * (ndim + 1));
  floating image = 9 * nvox * floatbytes;
  floating reg_mat = 2 * mat_size(ndim, nmapnod * (ndim + 1));
  floating map = 2 * nmapnod * floatbytes;

  floating total_bytes = nwarp_mat * warp_mat + nimage * image + nreg_mat * reg_mat + nmap * map;

  std::cout << "Total data memory needed = " << total_bytes << " for " << nvox << " voxels.\n"
            << "Comprising:\n"
            << "\tImage data and working copies: " << nimage << " x " << image << "\n"
            << "\tInterpolation matrices: " << nwarp_mat << " x " << warp_mat << "\n"
            << "\tMap data and working copy: " << nmap << " x " << map << "\n"
            << "\tSolver matrices: " << nreg_mat << " x " << reg_mat << "\n"
            << std::flush;
}

#endif // UTILS_HPP
