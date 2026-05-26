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
    <td>Moderate</td>
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
You can then run the applications in the container as you would normally. See the \ref
tutorial for more information.

*/

/** \page source Installing from Source

\subsection installation_prerequisites Prerequisites
Before you can install PRISMS-PF, you will need to have some prerequisites installed on
your machine. These include:
- A C++ compiler with C++20 support
- CMake (3.25+)
- MPI
- p4est (for adaptive octree meshing)
- deal.II (9.6.0+)

You can install these prerequisites however you like. The
first three are typically available on most OS through your package manager. For example,
on Ubuntu, you can install them with:
```
sudo apt-get install git lsb-release git subversion \
wget bc libgmp-dev build-essential autoconf automake \
cmake libtool gfortran python3 libboost-all-dev \
zlib1g-dev openmpi-bin openmpi-common libopenmpi-dev \
libblas3 libblas-dev liblapack3 liblapack-dev \
libsuitesparse-dev
```
To install p4est and deal.II, we recommend following the instructions on the [deal.II
installation guide](https://dealii.org/current_release/download/). Building deal.II and
p4est using deal.II's automatic installation script
[candi](https://github.com/dealii/candi) is by far the best option for most users.

```
git clone https://github.com/dealii/candi.git
cd candi
```
Modify the `candi.cfg` file, commenting out any PACKAGES other than dealii and p4est.
Also, turn on the native optimizations option. You should have something like the
following.

@code
# Option {ON|OFF}: Enable machine-specific optimizations (e.g. -march=native)?
NATIVE_OPTIMIZATIONS=ON
@endcode

@code
# These packages determine the active components of deal.II:
PACKAGES="${PACKAGES} once:p4est"
PACKAGES="${PACKAGES} dealii"
@endcode

@remark If you are building deal.II on a different architecture than you plan to run
PRISMS-PF, you can specify the architecture in the `DEAL_II_CONFOPTS=""` section with
`-march=<arch>`.

@remark If you have trouble with deal.II auto-detecting packages or with the system
installation of Boost try adding the following.
`DEAL_II_CONFOPTS="-DDEAL_II_ALLOW_AUTODETECTION=OFF -DDEAL_II_FORCE_BUNDLED_BOOST=ON"`

Then perform the configuration, compilation, and installation with
```
./candi.sh
```
or
```
./candi.sh -j <nprocs>
```
to compile with multiple processors (much faster).

@note We have found that specifying too many processors can sometimes lead to build
failures, so if you run into issues, try reducing the number of processors specified. If
this doesn't work \ref contact us.

This will install both p4est and deal.II in a local directory (by
default,`$HOME/dealii-candi/`). You can override the installation path with the prefix
command (`./candi.sh -p /path/to/install`). Be sure to permanently set the `DEAL_II_DIR`
environment variable to point to the deal.II installation. The candi script will prompt
you with steps on how to do this.

\subsection install_prismspf Installing PRISMS-PF

Once you have all the prerequisites installed, you can clone the PRISMS-PF repository
from GitHub, configure, build, and install the code.

First, clone the repository and navigate to the main directory.
```
git clone https://github.com/prisms-center/phaseField.git
cd phaseField
```

Next, we'll configure the library. We have a table of configuration options below.
Additionally, if you have CMake 3.28+ you can use our `CMakePresets.json`.


\subsubsection cmake_gen_config General Configuration
<table>
  <tr>
    <th></th>
    <th>Default</th>
    <th>Description</th>
  </tr>
  <tr>
    <td>`PRISMS_PF_AUTODETECTION`</td>
    <td>ON</td>
    <td>
    Try to detect \ref cmake_opt_dependencies that are in path
    </td>
  </tr>
  <tr>
    <td>`PRISMS_PF_UNIT_TESTS`</td>
    <td>OFF</td>
    <td>
    Build the unit tests
    </td>
  </tr>
  <tr>
    <td>`PRISMS_PF_REGRESSION_TESTS`</td>
    <td>OFF</td>
    <td>
    Build the regression tests
    </td>
  </tr>
  <tr>
    <td>`PRISMS_PF_PERFORMANCE_TESTS`</td>
    <td>OFF</td>
    <td>
    Build the performance tests
    </td>
  </tr>
  <tr>
    <td>`PRISMS_PF_EXAMPLES`</td>
    <td>OFF</td>
    <td>
    Build the examples
    </td>
  </tr>
  <tr>
    <td>`PRISMS_PF_DOCS`</td>
    <td>OFF</td>
    <td>
    Build the documentation
    </td>
  </tr>
  <tr>
    <td>`PRISMS_PF_CLANG_TIDY`</td>
    <td>OFF</td>
    <td>
    Run clang-tidy during the build stage
    </td>
  </tr>
</table>

\subsubsection cmake_feature_config Feature Configuration
<table>
  <tr>
    <th></th>
    <th>Default</th>
    <th>Description</th>
  </tr>
  <tr>
    <td>`PRISMS_PF_64BIT_INDICES`</td>
    <td>OFF</td>
    <td>
    Use 64-bit indices for large scale simulations (>2.147 billion DoFs). This relies on
the fact that deal.II is also configured with this option on.

Leaving this on when unnecessary will result in worse performance.
    </td>
  </tr>
  <tr>
    <td>`PRISMS_PF_ADDITIONAL_DEGREES`</td>
    <td>OFF</td>
    <td>
    Compile the core library with element degrees 3+. When turned off, explicit templates
will only be compiled for elements degrees 1 and 2. Turning it on will allow access to
higher element degrees on the application level, but increases compile time.

Note that deal.II is typically configured with element degrees up to 6. If you want higher
element degrees, you must configure deal.II accordingly.
    </td>
  </tr>
  <tr>
    <td>`PRISMS_PF_THREADS`</td>
    <td>ON</td>
    <td>
    Use threading when possible for MPI processes in `Release` mode. When turned off, each
MPI process will only use one thread.
    </td>
  </tr>
  <tr>
    <td>`PRISMS_PF_GPU`</td>
    <td>OFF</td>
    <td>
    <span style="color:red">[Experimental]</span> Use GPU backends instead of CPU. This is
still under active development
    </td>
  </tr>
</table>

\subsubsection cmake_opt_dependencies Optional Dependencies
<table>
  <tr>
    <th></th>
    <th>Default</th>
    <th>Description</th>
  </tr>
  <tr>
    <td>`PRISMS_PF_WITH_VTK`</td>
    <td>OFF</td>
    <td>
    Enable additional VTK file I/O options
    </td>
  </tr>
  <tr>
    <td>`PRISMS_PF_WITH_HDF5`</td>
    <td>OFF</td>
    <td>
    Enable HDF5 file I/O options
    </td>
  </tr>
  <tr>
    <td>`PRISMS_PF_WITH_CALIPER`</td>
    <td>OFF</td>
    <td>
    Enable Caliper as a profiler
    </td>
  </tr>
</table>

First, set the `PRISMS_PF_DIR` environment variable to point to where you want to install
PRISMS-PF. You can add this to your rc file to have permanent access to `$PRISMS_PF_DIR`.
```
export PRISMS_PF_DIR="path/to/install"
```

The most basic configuration is as follows:
```
cmake -DCMAKE_BUILD_TYPE=DebugRelease -B build
cmake --build build -j <nprocs>
cmake --install build --prefix $PRISMS_PF_DIR
```
where `<nprocs>` is how processors to build in parallel with.

This will build the PRISMS-PF library in both `Debug` and `Release` modes (something that
deal.II does that we think is nice). If you'd like to use typical CMake build types you
can also do that.

You can now compile and run any of the applications in the `applications` directory. For
example, to compile and run the explicit `allen_cahn` application, you can use the
following commands:
```
cd applications/allen_cahn/explicit
cmake -DCMAKE_BUILD_TYPE=Release .
make
mpirun -n <nprocs> main
```

Unlike before, there is no `DebugRelease` build type on the application level. As an
application you must choose between `Debug` and `Release`. `Debug` is the default build
type for applications.

@remark Let's talk about build types! `Debug` is slooww, concerningly so sometimes.
However, it's for good reason. We have a bunch of assertions that are being checked! If
you have something wrong with your application, you'll get a meaningful error that will
tell you how to fix it. This is good when developing applications. **Always run in `Debug`
when developing applications**. If you're playing around with model parameters or want to
run non-test systems use `Release`. However, if something goes wrong, you may get an
unintelligible error message or a SEGFAULT.
<br> <br>
Let me repeat myself. **Always run in `Debug` when developing applications**
<br>
\- @landinjm & @fractalsbyx

*/
