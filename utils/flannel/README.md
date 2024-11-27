# Flannel

It is a set of scripts made by Phil Tooley manipulate image data for ShIRT.

## Installation


```
# Clone the repository:
git clone

# Create a python virtual environment to isolate installation (optional)
python3 -m venv
source venv/bin/activate
 
# Install flannel from its source folder to avoid pip picks its own flannel package
cd flannel/ 
pip3 install .

   
```

Please note that there exist other packages named Flannel tha are also included in 
Pip packages. 

# Python API

```
image_to_shirt()
```

reads input from command line "/path/to/imagefile.ext"
Output: imagefile.image


```
shirt_to_image()
```
Reads input from command line  "/path/to/imgname.image"
Output: imagename.png"


```
image_to_mask()
```   
reads /path/to/imagefile.ext"
Output: imagefile.mask"
 



# Command line interface
Binaries converting image of supported formats to/from Shirt image, or generate mask.
Output files are named after the input file just replacing their extension. 



Convert image of supported extension to Shirt .image format
```
image2shirt  /path/to/imagefile.ext
```

Convert Shirt image to PNG format:
```
shirt2image /path/to/imgname.image
```


Generate mask for input Shirt image. Output file is named after input replacing extension to .mask
```
image2shirtmask  /path/to/imagefile.ext

```