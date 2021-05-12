# pFIRE documentation

### In code documentation
This project is mainly writte in C++. In code documentation is generated using 
Doxygen documentation that can be generated running:


```
# In pFIRE source folder run
doxygen

```

Doxygen is a tool to generate in code documentatation: [doxygen](https://www.doxygen.nl)

### User and developer documentation

All the remaining documentation is located in the `/doc` subfoldeer and is 
written using simple text files in Markdown and Restructured text format.
[Sphinx](https://www.sphinx-doc.org/en/master/) with Readthedocs template are used 
to generate documentation both in HTML and PDF format.

Documentation here provided is automatically generated and made available at:

Create Python virtual environment and install Sphinx dependencies:



```

# Create and Activate python environment
python -m venv venv_folder
source venv/bin/activate

# Install required python software

pip install -r requirements.txt

# Generated html documentation
sphinx-build -b html . build

# Generated pdf documentation
sphinx-build -b pdf . build

# When finished deactivate the Python environment
deactivate

```

Documentation is available in `doc/build/html`

### Online Documentation

pFire documentation is available at the following locations:

(Main read the docs documentation)[https://insigneo.github.io/pFIRE/docs.html]
or

(Insigneo pFIRE github documentation)[https://insigneo.github.io/pFIRE/docs.html#]




