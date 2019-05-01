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

#include "oiiowriter.hpp"

#include <OpenImageIO/imageio.h>

#include "image.hpp"
#include "exceptions.hpp"

const std::string OIIOWriter::writer_name = "OpenImageIO";
const std::vector<std::string> OIIOWriter::extensions = {".png",".jpg",".jpeg",".tiff"};

OIIOWriter::OIIOWriter(std::string filename, const MPI_Comm& comm)
  : BaseWriter(std::move(filename), comm)
{
}

std::string OIIOWriter::write_image(const Image& image)
{
  Vec_unique imgvec = image.scatter_to_zero();

  int rank;
  MPI_Comm_rank(_comm, &rank);
  if (rank == 0)
  {
    // This will be either a smart_ptr or a ImageOutput* depending on oiio version
    auto img = OIIO::ImageOutput::create(filename);
    if (img == nullptr)
    {
      throw WriterError("Failed to open image output file");
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

    //Manual ptr management if needed
#if OIIO_VERSION < 10903
    OIIO::ImageOutput::destroy(img);
#endif //OIIO_VERSION < 10903
  }
  MPI_Barrier(_comm);
  return filename;
}

std::string OIIOWriter::write_map(const Map& map __attribute__((unused)))
{
  throw InvalidWriterError("Cannot save map using OIIO.");
}
