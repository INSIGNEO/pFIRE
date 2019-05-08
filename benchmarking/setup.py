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
      install_requires=['numpy', 'scipy', 'h5py', 'matplotlib', 'flannel',
                        'configObj'],
      entry_points={
          'console_scripts':  [
              'pfire-shirt-compare = pfire_benchmarking.cross_validate:main',
              'pfire-regression-compare = pfire_benchmarking.regression_validate:main',
              'pfire-integration-test = pfire_benchmarking.__main__:main'
              ],
      })
