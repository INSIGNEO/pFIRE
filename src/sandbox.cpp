#include<numeric>
#include<cstdio>
#include<memory>

#include<OpenImageIO/imagecache.h>

#include "types.hpp"

void mainflow();

int main(int argc, char **argv){

  std::string fname(argv[1]);

  size_t x = 50;
  size_t y = 50;

  ImageCache_unique cache = create_unique_imagecache();

  OIIO::ImageSpec spec;
  bool ok = cache->get_imagespec (OIIO::ustring(fname), spec);
  if(ok)
  {
    std::cout << fname << " resolution is " << spec.width << "x"  << spec.height << "x"
              << spec.depth
              << " with " << spec.nchannels << " channels"
              << "(";
    for(auto &cn: spec.channelnames){ std::cout << cn << ", ";}
    std::cout << ")\n";
  }

  size_t nc = 1;//spec.nchannels;

  float pixels[nc*x*y];
  ok = cache->get_pixels(OIIO::ustring(fname),
                         0, 0, 0, x, 0, y, 0, 1, 0, 1,
                         OIIO::TypeDesc::FLOAT, pixels,
                         OIIO::AutoStride, OIIO::AutoStride, OIIO::AutoStride);

  if(ok)
  {
    for(size_t c = 0; c < nc; c++)
    {
      for(size_t i = 0; i < x; i++)
      {
        for(size_t j = 0; j < y; j++)
        {
          std::cout << pixels[(y*i + j)*nc + c] << " ";
        }
        std::cout << std::endl;
      }
    }
  }
  else
  {
    std::cout << "Failed to open " << fname << std::endl;
  }

  return 0;
}

