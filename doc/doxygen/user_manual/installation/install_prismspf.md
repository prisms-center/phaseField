# Installing PRISMS-PF {#install_prismspf}

## Traditional Installation (Non-Docker)
### Downloading PRISMS-PF
PRISMS-PF is available for download at our [GitHub site](https://github.com/prisms-center/phaseField). The recommended method for downloading the PRISMS-PF source code is through git. Using a Linux/Unix terminal, go to the directory where you want PRISMS-PF to be located. To clone the repository, type:
```
$ git clone https://github.com/prisms-center/phaseField.git
```
Git will then download the PRISMS-PF source code into a directory entitled ''phaseField''. A resource for learning to use Git can be found [here](https://git-scm.com/book).

If you prefer not to use Git, a zip file containing the PRISMS-PF source code can be downloaded [here](https://github.com/prisms-center/phaseField/releases).

### Compiling the Core Library
To compile the core PRISMS-PF library, enter the ```phaseField``` directory:
```
$ cd phaseField
```

Next, use CMake to generate a makefile and then compile the library:
```
$ cmake .
$ make -j 8
```
In the last command, the number following ```-j``` is the number of threads used during compilation. The general rule of thumb is to pick a number 1.5x the number of processors available. This process should take a minute or two and will compile both the ''debug'' and ''release'' versions of the core PRISMS-PF library.

Note: This step can be skipped and the core library will be compiled the first time that you compile one of the applications. However, the compilation will be on a single thread and thus will take much longer than if you follow the instructions above and use multiple threads.

## Docker Installation
One option when installing PRISMS-PF is to install it and all of its prerequisites using Docker. Docker is a container platform that can be used to run code identically on any OS X, Linux, or Windows machine. If you aren't familiar with Docker, you can get a quick overview [here](https://opensource.com/resources/what-docker).

1. [Install Docker](https://docs.docker.com/install/)
2. Create a directory on your computer to hold the PRISMS-PF applications folder (e.g. ```$ mkdir ~/DockerWorkspace```) and go to that directory (e.g. ```$ cd ~/DockerWorkspace```)
3. Clone the PRISMS-PF repository (```$ git clone https://github.com/prisms-center/phaseField```)
4. Get the PRISMS-PF Docker image (```$ docker pull prismspf/prismspf:latest```)
5. Now, launch the Docker container, linking the PRISMS-PF applications folder in the container to your local one: (e.g.  ```$ docker run -ti -v ~/DockerWorkspace/phaseField/applications:/home/dealii/phaseField/applications prismspf/prismspf:latest```)

Now you should have a working PRISMS-PF environment. The core library will already be compiled, but you'll need to compile any applications that you want to use. Any changes that you make or output files that you create should be found at ~/DockerWorkspace on your host computer.

To get out of the container, type:
```
$ exit
```

**Note:** If you only want to update (and use) the new Docker image, follow steps 4 and 5 above. If you want to remove the previous image, simply type: 
```
$ docker rmi stvdwtt/prismspf:prismspf
```
or
```
$ docker rmi prismspf/prismspf:2.2
```
