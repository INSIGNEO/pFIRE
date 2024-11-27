from functools import partial
import numpy as np
import h5py
import sys
import skimage.io as skio


def __LINE__():
    return str(sys._getframe(1).f_lineno)




def load_pfire_image(imagepath):
    """Load pFIRE image
    """
    print("\nload_pfire_image filepath={}".format(imagepath) )
    imagepath+=":/registered"  # FIXME data group should be read from xdmf file
    
    filename, group = [x.strip() for x in imagepath.split(':')]
    print("\nload_pfire_image filename={}".format(filename) )
    print("\nload_pfire_image group={}".format(group ))
    #print("load_pfire_image() line="+__LINE__())
            
    if filename.endswith(".xdmf"):    # FIXME read filename from xdmf
        filename += ".h5"
    
        reader = XdmfReader.New()
        domain = read(filename)
        print (domain)
        
    print("\nload_pfire_image filename={}".format(filename) )
    print("load_pfire_image() line="+__LINE__())
    
    with h5py.File(filename, 'r') as fh:
        print(list(fh.keys()))
        imgdata = np.asarray(fh['registered'])

    return imgdata



filename="registered.xdmf.h5"
fh=h5py.File(filename, 'r') 


print("\n open "+filename)
lista = list(fh.keys())

print(*lista)

for element in lista:
    print(element)
        
    
    
print(list(fh.keys()))
imgdata = np.asarray(fh['registered'])
type(imgdata)
print("\n type of imagedata =" + str(type(imgdata)))

print("\n shape of imagedata =" , imgdata.shape)
print("\n registered image\n =" , imgdata)


fh.close()

filename="map.xdmf.h5"
print("\n open "+filename)
fh=h5py.File(filename, 'r') 


for element in list(fh.keys()):
    print(element)
    
data = np.asarray(fh['/map/x'])
print("\n shape of data =" , data.shape)


##########################
imagepath="registered.xdmf"
image_data= load_pfire_image(imagepath)

print ("\n image_data= \n " )
print (image_data)

#datax=map.xdmf.h5:/map/x


 #  <Attribute Name="/map" AttributeType="Vector" Center="Node">
 #               <DataItem ItemType="Function" Dimensions="27 21 1 2" Function="join($0, $1)"/>
 #              <DataItem Name="dx" Format="HDF" Dimensions="27 21 1">map.xdmf.h5:/map/x</DataItem>
 #               <DataItem Name="dy" Format="HDF" Dimensions="27 21 1">map.xdmf.h5:/map/y</DataItem>
 #