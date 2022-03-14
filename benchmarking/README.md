## pFIRE Benchmark suite

The benchmark suite provides a framework to perform regression and comparison tests of pFire/ShiRT over a test cases organised in folders.


# Installation

The suite is written in Python and can be installed in the standard way.
We recoomend to use virtual environments as below:




```bash
	# Create virtual environment
	
	python3 -m venv pyvenv3
	
	
	# Activate
	
	source pyvenv3/bin/activate
	
	# Install dependencies
	
	pip install -r requirements.txt
	
	python setup install .
	
	# add folder to pythonpath
	
	export PYTHONPATH=<path_to>/pFIRE/benchmarking/
	
	
	# When not needed any more deactivate the virtual environment
	
	deactivate


```


# Usage


*  create a testconf file in each test folder
*  pfire_benchmark  <path_to_>/mytest.testconf



