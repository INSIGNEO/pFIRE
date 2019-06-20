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

#include "shirtloader.hpp"

#include <numeric>
#include <sstream>

#include "exceptions.hpp"
#include "file_utils.hpp"

const std::string ShIRTLoader::loader_name = "ShIRT";

BaseLoader_unique ShIRTLoader::Create_Loader(const std::string &path, MPI_Comm comm)
{
  return BaseLoader_unique(new ShIRTLoader(path, comm));
}

ShIRTLoader::ShIRTLoader(const std::string &path, MPI_Comm comm)
    : BaseLoader(path, comm), _file_type(undef)
{
  // Read header and then check total file size, if consistent then continue otherwise bail
  int rank;
  MPI_Comm_rank(comm, &rank);
  intcoord headershape;
  
  if (rank == 0)
  {
    MPI_File fh;
    int mpi_err;
    mpi_err = MPI_File_open(MPI_COMM_SELF, path.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
    if (mpi_err != MPI_SUCCESS)
    {
      headershape[0] = -1;
    }
    else
    {
      try
      {
        headershape = read_and_validate_image_header(fh);
      }
      catch (const InvalidLoaderError &)
      {
        try
        {
          headershape = read_and_validate_mask_header(fh);
        }
        catch (const InvalidLoaderError &)
        {
          headershape[0] = -1;
        }
      }
    }
    MPI_File_close(&fh);
  }

  MPI_Bcast(headershape.data(), headershape.size(), MPI_LONG, 0, comm);
  MPI_Bcast(&_file_type, 1, MPI_LONG, 0, comm);

  if (headershape[0] <= 0)
  {
    throw_if_nonexistent(path);
    throw InvalidLoaderError(path);
  }

  this->set_shape(headershape);
}

intcoord ShIRTLoader::read_and_validate_image_header(const MPI_File &fh)
{
  intcoord shape;

  MPI_Offset fsize;
  int mpi_err = MPI_File_get_size(fh, &fsize);
  if (mpi_err == MPI_SUCCESS)
  {
    std::vector<shirt_header_dtype> headerdata(image_header_length, 0);

    MPI_Status read_status;
    MPI_File_seek(fh, 0, MPI_SEEK_SET);
    MPI_File_read(fh, headerdata.data(), image_header_bytes, MPI_BYTE, &read_status);
    int read_count(0);
    MPI_Get_count(&read_status, MPI_BYTE, &read_count);

    // Currently only handle grayscale images
    if (headerdata[3] != 1)
    {
      throw InvalidLoaderError("Color images not currently supported");
    }

    if (mpi_err == MPI_SUCCESS && read_count == image_header_bytes)
    {
      integer imsize =
          std::accumulate(headerdata.cbegin(), headerdata.cend(), 1, std::multiplies<>());
      if (image_header_bytes + imsize * sizeof(image_data_dtype) == static_cast<uinteger>(fsize))
      {
        std::copy_n(headerdata.cbegin(), shape.size(), shape.begin());
      }
      else
      {
        throw InvalidLoaderError("Not a valid shirt image");
      }
    }
  }
  _file_type = image;
  return shape;
}

intcoord ShIRTLoader::read_and_validate_mask_header(const MPI_File &fh)
{
  intcoord shape;

  MPI_Offset fsize;
  int mpi_err = MPI_File_get_size(fh, &fsize);
  if (mpi_err == MPI_SUCCESS)
  {
    std::vector<shirt_header_dtype> headerdata(mask_header_length, 0);

    MPI_Status read_status;
    MPI_File_seek(fh, 0, MPI_SEEK_SET);
    MPI_File_read(fh, headerdata.data(), mask_header_bytes, MPI_BYTE, &read_status);
    int read_count(0);
    MPI_Get_count(&read_status, MPI_BYTE, &read_count);
    if (mpi_err == MPI_SUCCESS && read_count == mask_header_bytes)
    {
      integer imsize =
          std::accumulate(headerdata.cbegin(), headerdata.cend(), 1, std::multiplies<>());
      if (mask_header_bytes + imsize * sizeof(mask_data_dtype) == static_cast<uinteger>(fsize))
      {
        std::copy_n(headerdata.cbegin(), shape.size(), shape.begin());
      }
      else
      {
        throw InvalidLoaderError("Not a valid shirt image");
      }
    }
  }
  _file_type = mask;
  return shape;
}

void ShIRTLoader::copy_scaled_chunk(
    floating ***data, const intvector &chunksize, const intvector &offset) const
{
  std::vector<int> subsize(chunksize.cbegin(), chunksize.cend());
  std::vector<int> starts(offset.cbegin(), offset.cend());

  switch (_file_type)
  {
    case image:
      copy_chunk_image(data, subsize, starts);
      break;
    case mask:
      copy_chunk_mask(data, subsize, starts);
      break;
    default:
      throw InternalError("ShirtLoader attempted to load from non-shirt file type. This is a bug", __FILE__, __LINE__);
  }
}

void ShIRTLoader::copy_chunk_image(
    floating ***data, const std::vector<int> &subsize, const std::vector<int> &starts) const
{
  int dcount = std::accumulate(subsize.cbegin(), subsize.cend(), 1, std::multiplies<>());
  std::vector<int> size(this->shape().cbegin(), this->shape().cend());

  int irank;
  MPI_Comm_rank(this->comm(), &irank);

  std::vector<image_data_dtype> databuf(dcount, 0);

  MPI_Datatype file_layout;
  MPI_Type_create_subarray(
      size.size(), size.data(), subsize.data(), starts.data(), MPI_ORDER_FORTRAN,
      image_data_mpi_type, &file_layout);
  MPI_Type_commit(&file_layout);


  MPI_File fh;
  int mpi_err;
  mpi_err = MPI_File_open(this->comm(), this->path().c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
  if (mpi_err != MPI_SUCCESS)
  {
    throw_if_nonexistent(this->path());
    throw InvalidLoaderError(this->path()); 
  }

  mpi_err = MPI_File_set_view(
      fh, image_header_bytes, image_data_mpi_type, file_layout, "native", MPI_INFO_NULL);
  if (mpi_err != MPI_SUCCESS)
  {
    int rank;
    MPI_Comm_rank(this->comm(), &rank);
    std::ostringstream errss;
    std::string mpi_err_str;
    int str_size;
    mpi_err_str.resize(MPI_MAX_ERROR_STRING);
    MPI_Error_string(mpi_err, mpi_err_str.data(), &str_size);
    mpi_err_str.resize(str_size);
    errss << "[Rank " << rank << "] Failed to set view: " << mpi_err_str;
    throw InternalError(errss.str(), __FILE__, __LINE__);
  }

  MPI_Status read_status;
  MPI_File_read_all(fh, databuf.data(), dcount, image_data_mpi_type, &read_status);
  int read_count(0);
  MPI_Get_count(&read_status, image_data_mpi_type, &read_count);
  if (mpi_err != MPI_SUCCESS || read_count != dcount)
  {
    throw InternalError("Failed to read data chunk.", __FILE__, __LINE__);
  }

  floating *rankdata = &data[starts[2]][starts[1]][starts[0]];

  for (integer idx = 0; idx < dcount; idx++)
  {
    rankdata[idx] = databuf[idx];
  }

  MPI_File_close(&fh);
  MPI_Type_free(&file_layout);
}

void ShIRTLoader::copy_chunk_mask(
    floating ***data, const std::vector<int> &subsize, const std::vector<int> &starts) const
{
  int dcount = std::accumulate(subsize.cbegin(), subsize.cend(), 1, std::multiplies<>());
  std::vector<int> size(this->shape().cbegin(), this->shape().cend());

  std::vector<mask_data_dtype> databuf(dcount, 0);

  MPI_Datatype file_layout;
  MPI_Type_create_subarray(
      size.size(), size.data(), subsize.data(), starts.data(), MPI_ORDER_FORTRAN,
      mask_data_mpi_type, &file_layout);
  MPI_Type_commit(&file_layout);

  MPI_File fh;
  int mpi_err;
  mpi_err = MPI_File_open(this->comm(), this->path().c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
  if (mpi_err != MPI_SUCCESS)
  {
    throw_if_nonexistent(this->path());
    throw InvalidLoaderError(this->path()); 
  }

  mpi_err = MPI_File_set_view(
      fh, mask_header_bytes, mask_data_mpi_type, file_layout, "native", MPI_INFO_NULL);
  if (mpi_err != MPI_SUCCESS)
  {
    int rank;
    MPI_Comm_rank(this->comm(), &rank);
    std::ostringstream errss;
    std::string mpi_err_str;
    int str_size;
    mpi_err_str.resize(MPI_MAX_ERROR_STRING);
    MPI_Error_string(mpi_err, mpi_err_str.data(), &str_size);
    mpi_err_str.resize(str_size);
    errss << "[Rank " << rank << "] Failed to set view: " << mpi_err_str;
    throw InternalError(errss.str(), __FILE__, __LINE__);
  }

  MPI_Status read_status;
  MPI_File_read_all(fh, databuf.data(), dcount, mask_data_mpi_type, &read_status);
  int read_count(0);
  MPI_Get_count(&read_status, mask_data_mpi_type, &read_count);
  if (mpi_err != MPI_SUCCESS || read_count != dcount)
  {
    throw InternalError("Failed to read data chunk.", __FILE__, __LINE__);
  }

  floating *rankdata = &data[starts[2]][starts[1]][starts[0]];

  for (integer idx = 0; idx < dcount; idx++)
  {
    rankdata[idx] = databuf[idx];
  }

  MPI_File_close(&fh);
  MPI_Type_free(&file_layout);
}
