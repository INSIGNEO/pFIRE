#!/usr/bin/env python

from setuptools import setup, find_packages

setup(name='flannel',
      version='0.1',
      description='Routines for interacting with the original ShIRT',
      author='Phil Tooley, INSIGNEO',
      author_email='phil.tooley@sheffield.ac.uk',
      url='insigneo.org',
      package_dir={'':'src'},
      packages=find_packages('src'),
      install_requires=['numpy', 'scikit-image'],
      entry_points={
        'console_scripts' : [
          'shirt2image = flannel.helpers.image:shirt_to_image',
          'image2shirt = flannel.helpers.image:image_to_shirt',
          'image2shirtmask = flannel.helpers.image:image_to_mask',
          ],
        'gui_scripts' : []
      }
)
