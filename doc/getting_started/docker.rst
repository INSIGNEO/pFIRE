===================================
Using pFIRE in Docker
===================================

Docker is a software that creates a package (called image) containing an application, its library dependencies and a minimal execution environment (e.g. a minimal Ubuntu distribution). These images can be distributed and run by any Docker software installed on the user machine. In the following are provided the instructions to execute pFIRE on a user machine with Docker. It as been tested on *nix,  OsX  machines and on Windows10 throught the Windows Linux Subsystem.


There two ways to use pFIRE in a docker:
   * download a ready made image from pfire dockerhub channel
   * build the docker image from the source code	


Download a container image ready made from Dockerhub
--------------------------------------------------------

On the Dockerhub portal are hosted some versions of pFIRE already packaged in a Docker image.
Follow the instructions on the web page to download and use: 

_PETSc test 


Build yourself from dockerfile
---------------------------------
If you have installed Docker with command line support you can generate the docker image from the source docker container file.
Instructions and docker files organised by destination platform are available here:

.. _PETSc: https://www.mcs.anl.gov/petsc/

.. `_pFIRE: docker container source code`: 
	https://github.com/insigneo-pfire/docker-pfire






