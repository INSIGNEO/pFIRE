#include<numeric>
#include<cstdio>
#include<memory>

#include<OpenImageIO/imageio.h>

#include "types.hpp"

void mainflow();

int main(int argc, char **argv){

  integer fullx = 100, fully = 100;
  integer cropx = 10, cropy = 10;
  integer offsetx = 45, offsety = 45;

  std::string fname(argv[1]);

  OIIO::ImageOutput *out = OIIO::ImageOutput::create(fname);
  OIIO::ImageSpec spec(cropx, cropy, 1, OIIO::TypeDesc::UINT8);
  spec.full_x = 0;
  spec.full_y = 0;
  spec.full_width = fullx;
  spec.full_height = fully;
  spec.x = offsetx;
  spec.y = offsety;

  return 0;
}

