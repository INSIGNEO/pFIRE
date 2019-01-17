#include "dcmloader.hpp"

#include <dcmtk/config/osconfig.h>
#include <dcmtk/dcmdata/dctagkey.h>
#include <dcmtk/dcmdata/dctk.h>
#include <dcmtk/dcmimgle/dcmimage.h>

#include "exceptions.hpp"
#include "file_utils.hpp"

constexpr uinteger METADATA_MAX_BYTES = 4 * 1024;

constexpr unsigned int TG_IMG(0x0028);
constexpr unsigned int TE_ROWS(0x0010);
constexpr unsigned int TE_COLS(0x0011);
constexpr unsigned int TE_FRAMES(0x0008);

const std::string DCMLoader::loader_name = "DICOM";

BaseLoader_unique DCMLoader::Create_Loader(const std::string &path, MPI_Comm comm)
{
  return BaseLoader_unique(new DCMLoader(path, comm));
}

DCMLoader::DCMLoader(const std::string &path, MPI_Comm comm)
    : BaseLoader(path, comm), _datafile(DcmFileFormat())
{
  dcmtk::log4cplus::Logger rootLogger = dcmtk::log4cplus::Logger::getRoot();
  rootLogger.setLogLevel(OFLogger::OFF_LOG_LEVEL);

  OFCondition status = _datafile.loadFile(
      path.c_str(), EXS_Unknown, EGL_noChange, METADATA_MAX_BYTES, ERM_autoDetect);

  if (status.bad())
  {
    throw_if_nonexistent(path);
    throw InvalidLoaderError(path);
  }

  DcmDataset *dataset = _datafile.getDataset();

  long tmp;
  status = dataset->findAndGetLongInt(DcmTagKey(TG_IMG, TE_ROWS), tmp);
  if (status.bad())
  {
    throw std::runtime_error("Failed to read image shape data");
  }
  this->_shape[0] = tmp;
  status = dataset->findAndGetLongInt(DcmTagKey(TG_IMG, TE_COLS), tmp);
  if (status.bad())
  {
    throw std::runtime_error("Failed to read image shape data");
  }
  this->_shape[1] = tmp;
  status = dataset->findAndGetLongInt(DcmTagKey(TG_IMG, TE_FRAMES), tmp);
  if (status.bad())
  {
    throw std::runtime_error("Failed to read image shape data");
  }
  this->_shape[2] = tmp;
}

void DCMLoader::copy_scaled_chunk(
    floating ***data, const intvector &size, const intvector &offset) const
{
  DcmDataset *dataset = this->_datafile.getDataset();
  DicomImage img(
      dataset, EXS_Unknown, (unsigned long)(0), (unsigned long)offset[2],
      (unsigned long)size[2]); // NB frames in img are reindexed from the first loaded frame

  uinteger bitdepth = img.getDepth();
  double dmin, dmax;
  img.getMinMaxValues(
      dmin, dmax, 1); // mode 1 gets absolute min/max, 0 gets actual from loaded frames

  for (integer slice_idx = offset[2]; slice_idx < offset[2] + size[2]; slice_idx++)
  {
    for (integer col_idx = offset[1]; col_idx < offset[1] + size[1]; col_idx++)
    {
      integer flat_idx = (col_idx + offset[1]) * size[0] + offset[0];
      integer frame_idx = slice_idx - offset[2];
      switch (bitdepth)
      {
        case 8:
        {
          uint8_t *rawdata = (uint8_t *)img.getOutputData(8, slice_idx, 0);
          norm_convert(data[frame_idx][col_idx], rawdata + flat_idx, size[0], dmin, dmax);
          break;
        }
        case 16:
        {
          uint16_t *rawdata = (uint16_t *)img.getOutputData(16, slice_idx, 0);
          norm_convert(data[frame_idx][col_idx], rawdata + flat_idx, size[0], dmin, dmax);
          break;
        }
        default:
        {
          throw std::runtime_error("Unhandled image bit depth.");
          break;
        }
      }
    }
  }
}
