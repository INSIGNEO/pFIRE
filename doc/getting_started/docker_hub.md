# pFIRE on Dockerhub #

Docker is a software that creates a package (called image) containing an application, its library dependencies and a  minimal execution environment (e.g. a minimal Ubuntu distribution). These images can be distributed and run by any Docker software installed on the user machine. 
In the following are provided the instructions to execute pFIRE on a user machine with Docker. It as been tested on *nix and OsX machines and Windows10 through Windows Linux Subsystem and (Powerwshell).


### Installation procedure
* install Docker following instructions:
https://www.docker.com/products/docker-desktop

* Download docker image

```bash
docker pull insigneopfire/pfire:v0.4.0

# rename image
docker image tag insigneopfire/pfire:v0.4.0 pfire_image

```

### Usage instruction binding a local data folder



Create a folder test_docker containing the input images to be registered  (fixed.png, moved.png).

Create a textual file  (pfire.cfg) containing the pFIRE execution parameters as below 

```
fixed = fixed.png
moved = moved.png
nodespacing = 10

```

Now you can run pFire using Docker:

```bash

 docker run -d -it --name dockerised_pfire -v "$(pwd)"/test_docker:/test_docker  pfire_image

 docker exec -it  -w  /test_docker dockerised_pfire  pfire  pfire.cfg

```

pFIRE should show its running output and results are available in you local test_docker folder.

```bash
fixed.png  
moved.png  
pfire.cfg 
map.xdmf  
map.xdmf.h5   
registered.xdmf  
registered.xdmf.h5
```

