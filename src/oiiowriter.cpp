#include "oiiowriter.hpp"

#include <OpenImageIO/imageio.h>

#include "image.hpp"

OIIOWriter::OIIOWriter(std::string filename, const MPI_Comm& comm)
  : BaseWriter(std::move(filename), comm)
{
}

void OIIOWriter::write_image(const Image& image, const std::string& groupname)
{
  Vec_unique imgvec = create_unique_vec();
  imgvec = image.scatter_to_zero(*imgvec);

  integer rank;
  MPI_Comm_rank(_comm, &rank);
  if (rank == 0)
  {
    OIIO::ImageOutput* img = OIIO::ImageOutput::create(filename);
    if (img == nullptr)
    {
      throw std::runtime_error("Failed to open image output file");
    }
    OIIO::ImageSpec spec(image.shape()[0], image.shape()[1], 1, OIIO::TypeDesc::UINT16);
    img->open(filename, spec);

    floating* pixdata;
    PetscErrorCode perr = VecGetArray(*imgvec, &pixdata);
    CHKERRABORT(_comm, perr);
    img->write_image(OIIO::TypeDesc::DOUBLE, pixdata);
    img->close();
    perr = VecRestoreArray(*imgvec, &pixdata);
    CHKERRABORT(_comm, perr);

    OIIO::ImageOutput::destroy(img);
  }
  MPI_Barrier(_comm);
}

void OIIOWriter::write_map(const Map& map, const std::string& groupname)
{
  throw std::runtime_error("Cannot save map using OIIO.");
}
