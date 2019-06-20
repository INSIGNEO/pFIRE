#ifndef IMAGE_HPP
#define IMAGE_HPP

#include "types.hpp"

#include "imagebase.hpp"

//! Calculate gradient vectors for T matrix
std::vector<Vec_unique> calculate_tmatrix_gradients(const Image& fixed, const Image& moved, integer ndof);

class Image: public ImageBase {
public:
  Image(const intcoord& shape, const MPI_Comm& comm) : ImageBase(shape, comm) {}

  std::shared_ptr<Image> warp(const MapBase& map) const;

  floatvector interpolate_pixel_data(const floatcoordvector& source_locations) const;

  friend std::vector<Vec_unique> calculate_tmatrix_gradients(const Image& fixed, const Image& moved,
                                                             integer ndof);
};

#endif // IMAGE_HPP
