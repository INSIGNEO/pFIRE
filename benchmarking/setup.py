#!/usr/bin/env python

from setuptools import setup, find_packages

setup(name='pfire_benchmarking',
      version='0.1',
      description='pFIRE-ShIRT cross validation tools',
      author='Phil Tooley, INSIGNEO',
      author_email='phil.tooley@sheffield.ac.uk',
      url='insigneo.org',
      packages=find_packages(''),
      include_package_data=True,
      entry_points={
          'console_scripts':  [
              'pfire-shirt-compare = pfire_benchmarking.validate:main',
              ],
      })
