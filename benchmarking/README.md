## pFIRE Benchmark suite

The benchmark suite provides a framework to perform regression and comparison tests of pFire/ShiRT over a test cases organised in folders.


# Installation

The suite is written in Python and can be installed in the standard way.
We recommend to use virtual environments as below:



```bash
	# Create virtual environment
	
	cd pFIRE/benchmarking
	
	python3 -m venv pyvenv3
	
	
	# Activate
	
	source pyvenv3/bin/activate
	
	# Upgrade pip	
	pip install --upgrade pip
		
	# Install dependencies
	
	pip install -r requirements.txt
	
		##python setup.py install .
	pip install .
	
	# add folder to pythonpath
	
	export PYTHONPATH=`pwd`/:$PYTHONPATH
	
	
	../../pFIRE/testdata/integration/regression
	
	# When not needed any more deactivate the virtual environment
	
	deactivate


```


add pfire binary to PATH

export PATH=/<path-to-pfire-binary>/:$PATH

# Usage

```
pfire-integration-test --pfire-executable=</path-to-pfire>/<pfire-executable>  <base directory containing test configuration file> 
```

It is optional to provide the pFIRE executable filename. Apptainer allows to use the container as an executable (a wrapper around a container binary configured un the Runscript section)
 
 
*  create a testconf file in each test folder
*  pfire_benchmark  <path_to_>/mytest.testconf




TODO add testconf syntax and use cases




