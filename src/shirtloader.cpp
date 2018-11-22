#include "shirtloader.hpp"

#include <sstream>
#include <numeric>

const std::string ShIRTLoader::loader_name = "ShIRT";

BaseLoader_unique ShIRTLoader::Create_Loader(const std::string &path, MPI_Comm comm)
{
  return BaseLoader_unique(new ShIRTLoader(path, comm));
}


ShIRTLoader::ShIRTLoader(const std::string& path,  MPI_Comm comm)
  : BaseLoader(path, comm), _file_type(undef)
{
  //Read header and then check total file size, if consistent then continue otherwise bail
  int rank;
  MPI_Comm_rank(comm, &rank);

  MPI_File fh;
  int mpi_err;
  mpi_err = MPI_File_open(_comm, path.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
  if (mpi_err != MPI_SUCCESS){
    std::ostringstream errss;
    errss << "Failed to open file " << path << ".";
    throw std::runtime_error(errss.str());
  }
  
  if (rank == 0)
  {
    MPI_File fh;
    int mpi_err;
    mpi_err = MPI_File_open(MPI_COMM_SELF, path.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
    if (mpi_err != MPI_SUCCESS){
      std::ostringstream errss;
      errss << "Failed to open file " << path << ".";
      throw std::runtime_error(errss.str());
    }
    try
    {
      _shape = read_and_validate_image_header(fh);
    }
    catch (const std::runtime_error&)
    {
      try
      {
        _shape = read_and_validate_mask_header(fh);
      }
      catch (const std::runtime_error&)
      {
        _shape[0] = -1;
      }
    }
  MPI_File_close(&fh);
  }

  MPI_Bcast(_shape.data(), _shape.size(), MPI_LONG, 0, comm);
  MPI_Bcast(&_file_type, 1, MPI_LONG, 0, comm);

  if (_shape[0] <= 0)
  {
    std::ostringstream errss;
    errss << path << " is not a valid ShIRT image.";
    throw std::runtime_error(errss.str());
  }
}


intvector ShIRTLoader::read_and_validate_image_header(const MPI_File &fh)
{
  intvector shape(3, 0);

  MPI_Offset fsize;
  int mpi_err = MPI_File_get_size(fh, &fsize);
  if (mpi_err == MPI_SUCCESS)
  {
    std::vector<shirt_header_dtype> headerdata(image_header_length, 0);

    MPI_Status read_status;
    mpi_err = MPI_File_read(fh, headerdata.data(), image_header_bytes, MPI_BYTE, &read_status);
    int read_count(0);
    MPI_Get_count(&read_status, MPI_BYTE, &read_count);

    //Currently only handle grayscale images
    if(headerdata[3] != 1)
    {
      throw std::runtime_error("Color images not currently supported");
    }

    if(mpi_err == MPI_SUCCESS && read_count == image_header_bytes)
    {
      integer imsize = std::accumulate(headerdata.cbegin(), headerdata.cend(), 1,
                                       std::multiplies<>());
      if (image_header_bytes  + imsize * sizeof(image_data_dtype) == fsize)
      {
        std::copy_n(headerdata.cbegin(), shape.size(), shape.begin()); 
      }
      else
      {
        throw std::runtime_error("Not a valid shirt image");
      }
    }
  }
  _file_type = image;
  return shape;
}

intvector ShIRTLoader::read_and_validate_mask_header(const MPI_File &fh)
{
  intvector shape(3, 0);

  MPI_Offset fsize;
  int mpi_err = MPI_File_get_size(fh, &fsize);
  if (mpi_err == MPI_SUCCESS)
  {
    std::vector<shirt_header_dtype> headerdata(mask_header_length, 0);

    MPI_Status read_status;
    mpi_err = MPI_File_read(fh, headerdata.data(), mask_header_bytes, MPI_BYTE, &read_status);
    int read_count(0);
    MPI_Get_count(&read_status, MPI_BYTE, &read_count);
    if(mpi_err == MPI_SUCCESS && read_count == mask_header_bytes)
    {
      integer imsize = std::accumulate(headerdata.cbegin(), headerdata.cend(), 1,
                                       std::multiplies<>());
      if (mask_header_bytes  + imsize * sizeof(mask_data_dtype) == fsize)
      {
        std::copy_n(headerdata.cbegin(), shape.size(), shape.begin()); 
      }
      else
      {
        throw std::runtime_error("Not a valid shirt image");
      }
    }
  }
  _file_type = mask;
  return shape;
}

void ShIRTLoader::copy_scaled_chunk(floating ***data, const intvector &chunksize,
                                  const intvector &offset) const
{
  std::vector<int> subsize(chunksize.cbegin(), chunksize.cend());
  std::vector<int> starts(offset.cbegin(), offset.cend());

  switch(_file_type)
  {
    case image: 
      copy_chunk_image(data, subsize, starts);
      break;
    case mask : 
      copy_chunk_mask(data, subsize, starts);
      break;
    default: 
      throw std::runtime_error("Attempted to load from non-shirt file type. This is a bug");
  }
}

void ShIRTLoader::copy_chunk_image(floating ***data, const std::vector<int> &subsize,
                                   const std::vector<int> &starts) const
{
  int dcount = std::accumulate(subsize.cbegin(), subsize.cend(), 1, std::multiplies<>());
  std::vector<int> size(_shape.cbegin(), _shape.cend());

  std::vector<image_data_dtype> databuf(dcount, 0);

  MPI_Datatype file_layout;
  MPI_Type_create_subarray(size.size(), size.data(), subsize.data(), starts.data(),
                           MPI_ORDER_FORTRAN, image_data_mpi_type, &file_layout);
  MPI_Type_commit(&file_layout);

  MPI_File fh;
  int mpi_err;
  mpi_err = MPI_File_open(_comm, _path.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
  if (mpi_err != MPI_SUCCESS)
  {
    std::ostringstream errss;
    errss << "Failed to open file " << _path << ".";
    throw std::runtime_error(errss.str());
  }

  mpi_err = MPI_File_set_view(fh, image_header_bytes, image_data_mpi_type, file_layout,
                              "native", MPI_INFO_NULL);
  if (mpi_err != MPI_SUCCESS)
  {
    int rank;
    MPI_Comm_rank(_comm, &rank);
    std::ostringstream errss;
    std::string mpi_err_str;
    int str_size;
    mpi_err_str.resize(MPI_MAX_ERROR_STRING);
    MPI_Error_string(mpi_err, mpi_err_str.data(), &str_size);
    mpi_err_str.resize(str_size);
    errss << "[Rank " << rank << "] Failed to set view: " << mpi_err_str;
    throw std::runtime_error(errss.str());
  }

  MPI_Status read_status;
  MPI_File_read_all(fh, databuf.data(), dcount, image_data_mpi_type, &read_status);
  int read_count(0);
  MPI_Get_count(&read_status, image_data_mpi_type, &read_count);
  if(mpi_err != MPI_SUCCESS || read_count != dcount)
  {
    throw std::runtime_error("Failed to read data chunk.");
  }

  floating  *rankdata = &data[starts[2]][starts[1]][starts[0]];

  for(integer idx=0; idx<dcount; idx++)
  {
    rankdata[idx] = databuf[idx];
  }

  MPI_File_close(&fh);
  MPI_Type_free(&file_layout);
}

void ShIRTLoader::copy_chunk_mask(floating ***data, const std::vector<int> &subsize,
                                   const std::vector<int> &starts) const
{
  int dcount = std::accumulate(subsize.cbegin(), subsize.cend(), 1, std::multiplies<>());
  std::vector<int> size(_shape.cbegin(), _shape.cend());

  std::vector<mask_data_dtype> databuf(dcount, 0);

  MPI_Datatype file_layout;
  MPI_Type_create_subarray(size.size(), size.data(), subsize.data(), starts.data(),
                           MPI_ORDER_FORTRAN, mask_data_mpi_type, &file_layout);
  MPI_Type_commit(&file_layout);

  MPI_File fh;
  int mpi_err;
  mpi_err = MPI_File_open(_comm, _path.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
  if (mpi_err != MPI_SUCCESS)
  {
    std::ostringstream errss;
    errss << "Failed to open file " << _path << ".";
    throw std::runtime_error(errss.str());
  }

  mpi_err = MPI_File_set_view(fh, mask_header_bytes, mask_data_mpi_type, file_layout,
                              "native", MPI_INFO_NULL);
  if (mpi_err != MPI_SUCCESS)
  {
    int rank;
    MPI_Comm_rank(_comm, &rank);
    std::ostringstream errss;
    std::string mpi_err_str;
    int str_size;
    mpi_err_str.resize(MPI_MAX_ERROR_STRING);
    MPI_Error_string(mpi_err, mpi_err_str.data(), &str_size);
    mpi_err_str.resize(str_size);
    errss << "[Rank " << rank << "] Failed to set view: " << mpi_err_str;
    throw std::runtime_error(errss.str());
  }

  MPI_Status read_status;
  MPI_File_read_all(fh, databuf.data(), dcount, mask_data_mpi_type, &read_status);
  int read_count(0);
  MPI_Get_count(&read_status, mask_data_mpi_type, &read_count);
  if(mpi_err != MPI_SUCCESS || read_count != dcount)
  {
    throw std::runtime_error("Failed to read data chunk.");
  }

  floating  *rankdata = &data[starts[2]][starts[1]][starts[0]];

  for(integer idx=0; idx<dcount; idx++)
  {
    rankdata[idx] = databuf[idx];
  }

  MPI_File_close(&fh);
  MPI_Type_free(&file_layout);
}
