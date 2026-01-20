/** \page installation Installation

Welcome to the installation guide for PRISMS-PF! This is arguably the most difficult part
for most of our users, but we hope to make it as easy as possible. If you run into any
issues, please feel free to \ref contact us.

There are many ways to install PRISMS-PF, and we will try to cover the most common ones
here.

<table>
  <tr>
    <th>Method</th>
    <th>Difficulty</th>
    <th>Notes</th>
  </tr>
  <tr>
    <td>Docker</td>
    <td>Low</td>
    <td>
      Docker is potentially the easiest way to install PRISMS-PF. However, it is not
performant and unsuitable for large simulations.
    </td>
  </tr>
  <tr>
    <td>Source</td>
    <td>Moderate to High</td>
    <td>From source is the most flexible and performant way to install PRISMS-PF. However,
it does require some basic knowledge of build systems like CMake</td>
  </tr>
</table>

Once you have decided on a method, you can follow the instructions below:
<ul class="nav-submenu">
  <li><a href="docker.html">Install with Docker</a></li>
  <li><a href="source.html">Install from Source</a></li>
</ul>

*/

/** \page docker Installing with Docker
As mentioned in the \ref installation page, Docker is potentially the easiest way to
install PRISMS-PF. However, it is not performant and unsuitable for large simulations.
We'll get into those reasons with a brief introduction to Docker.

\subsection docker_introduction Introduction to Docker
Docker is a platform that allows you to run applications in containers. What are
containers? Containers are typically stripped down versions of virtual machines that allow
you to run an application in an isolated and reproducible environment. This custom
environment is defined by a Dockerfile, which is a text file that contains all the
instructions to build the container image. The image can then be run on any machine that
has Docker installed, regardless of the underlying operating system.

Importantly for PRISMS-PF, Docker allows you to run our software without needing to
install any dependencies. All you have to do is download the Docker image! Depending on
your operating system and the difference (if any) between the your hardware and the
emulated hardware, the running of a Docker image may add some overhead, which will slow
down your simulations. Additionally, with the Docker image, you're stuck with the image we
provide. If you want to change the compiler, the configuration of dependencies, or any
other things, you should consider \ref source instead. If you would like to learn more
about Docker, [here](https://opensource.com/resources/what-docker) is a nice resource.

\subsection docker_installation Installation with Docker
To install PRISMS-PF with Docker, you will need to have Docker installed on your machine.
You can find instructions for installing Docker on the official [Docker
website](https://docs.docker.com/get-started/get-docker/).

After you have installed Docker, find a working directory where we can put your files.
This is where we will clone PRISMS-PF so you can interact with the files and applications
there.
```
mkdir ~/prisms_pf_docker
git clone https://github.com/prisms-center/phaseField.git
```

Next, we will pull the docker image with
```
docker pull prismspf/prismspf:latest
```

Finally, we will launch an interactive container with
```
docker run -ti -v
~/prisms_pf_docker/phaseField/applications:/home/dealii/phaseField/applications
prismspf/prismspf:latest
```
This will link your local applications directory (the one in `prisms_workspace`) to the
one in the Docker image. If you plan to modify the core library, you should link one
directory higher to preserve your changes. In other words,
```
docker run -ti -v
~/prisms_pf_docker/phaseField:/home/dealii/phaseField prismspf/prismspf:latest
```
You can then run the applications in the container as you would normally. For example,
to run the `allen_cahn_explicit` application, you can use the following commands:
```
cd allen_cahn_explicit
cmake .
make
mpirun -n 1 ./main
```

*/

/** \page source Installing from Source

\subsection installation_prerequisites Prerequisites
Before you can install PRISMS-PF, you will need to have some prerequisites installed on
your machine. These include:
- Git (to clone the repository)
- A C++ compiler with C++20 support (we recommend GCC or Clang)
- CMake (version 3.13.4 or higher)
- LAPACK (for linear algebra operations)
- MPI (we recommend OpenMPI or MPICH)
- p4est (for adaptive octree meshing)
- VTK (for file I/O) **[optional]**
- Caliper (for profiling) **[optional]**
- deal.II (version 9.6 or higher)

<div class="tabbed">

- <b class="tab-title">Linux</b>
You can install these prerequisites however you like. The
first five are typically available on most OSes through your package manager. For example,
on Ubuntu, you can install them with:
```
sudo apt-get install git lsb-release git subversion wget bc libgmp-dev build-essential
autoconf automake cmake libtool gfortran python3 libboost-all-dev zlib1g-dev openmpi-bin
openmpi-common libopenmpi-dev libblas3 libblas-dev liblapack3 liblapack-dev
libsuitesparse-dev
```
To install p4est and deal.II, we recommend following the instructions on the [deal.II
installation guide](https://dealii.org/current_release/download/).
We suggest building
deal.II and p4est using deal.II's automatic installation script
[candi](https://github.com/dealii/candi).

```
git clone https://github.com/dealii/candi.git
cd candi
```
Modify the `candi.cfg` file, commenting out any PACKAGES other than dealii and p4est.
Then run:
```
./candi.sh
```
or
```
./candi.sh -j <nprocs>
```
to compile with multiple processors (much faster). (Note that we have found that
specifying too many processors can sometimes lead to build failures, so if you run into
issues, try reducing the number of processors specified.)

This will install both p4est and deal.II in a local directory (by
default,`$HOME/dealii-candi/`). Be sure to permanently set the `DEAL_II_DIR` environment
variable to point to the deal.II installation.

To install VTK, we recommend installing from source.
```
git clone https://gitlab.kitware.com/vtk/vtk.git vtk-clone
mkdir vtk-build
mkdir vtk-install
cd vtk-build
cmake ../vtk-clone -DCMAKE_BUILD_TYPE=Release -DVTK_GROUP_ENABLE_QT=OFF
-DVTK_GROUP_ENABLE_MPI=YES -DVTK_USE_MPI=ON -DVTK_BUILD_TESTING=OFF
-DCMAKE_INSTALL_PREFIX=../vtk-install

make -j <numprocs> install
```
Be sure to permanently set the `VTK_DIR` environment variable to point to the VTK
installation (e.g.,`$HOME/vtk-install/lib/cmake/vtk-9.xx`).

- <b class="tab-title">MacOS</b> This is the content of tab 2


- <b class="tab-title">Windows</b> **PRISMS-PF has no plans to support Windows.** If you
are on Windows, we recommend using Docker, another virtual machine, or Windows Subsystem
for Linux (WSL) to run a Linux environment. Please follow the Linux instructions if you
use any of those options.

</div>

\subsection install_prismspf Installing PRISMS-PF

Once you have all the prerequisites installed, you can clone the PRISMS-PF repository
from GitHub, configure the build with CMake, and compile the code:
```
git clone https://github.com/prisms-center/phaseField.git
cd phaseField
mkdir build
cd build
cmake -DPRISMS_PF_WITH_VTK=ON -DCMAKE_INSTALL_PREFIX=<prisms-install-dir> ..
make -j <nprocs>
```
This will build the PRISMS-PF library in both debug and release modes.
Next, set the `PRISMS_PF_DIR` environment variable to point to the PRISMS-PF installation
directory (the directory containing the `lib` and `include` directories).
You can now compile and run any of the applications in the `applications` directory. For
example, to compile and run the `allen_cahn_explicit` application, you can use the
following commands:
```
cd applications/allen_cahn_explicit
cmake .
make
mpirun -n <nprocs> ./main
```
(Note: there is a known bug with compiling application using multiple processors for the
first time. This will be fixed in a future release. After compiling once, the issue should
go away.)


*/
