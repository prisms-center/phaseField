# Installing CMake and deal.II {#install_prereqs}

Before downloading and installing PRISMS-PF itself, one should install CMake and deal.II. For MacOS or Linux users, it is likely that your distribution already contains CMake. To test if you have it type:
```
$ cmake --version
```

## Installing CMake
CMake can be downloaded [here](https://cmake.org/download), which also has installation instructions. It is also available through package managers such as homebrew, yum, and apt.

Once installed, be sure to add CMake to your path by entering: 
```
$ export PATH=/path/to/cmake/Contents/bin:$PATH
```

filling in the path to the installation of CMake (for example, on Mac OS, the default installation path is /Applications/CMake.app/Contents/bin). For convenience, we recommend adding this line to your bash profile. Your bash profile can be opened via the following terminal command
```
$ vi ~/.profile
```

## Installing deal.II
deal.II can be downloaded [here](https://www.dealii.org/current_release/download/), which also has installation instructions. The deal.II installation process depends on your operating system. The current PRISMS-PF version has been tested with deal.II versions 9.2.0 through 9.5.1. We recommend using the latest version of deal.II, if possible. Other versions are available [here](https://github.com/dealii/dealii/releases).

There are three general approaches to installing deal.II: using a binary package, installing from source, and using Docker. For general use on a personal computer, any of these are appropriate. When available, we recommend installing deal.II using a binary package, this should be the most straightforward approach. If a binary package isn't available, installation using Docker is also fairly straightforward. The third option is to install deal.II and its prerequisites from source, which is more difficult than the other options. When using PRISMS-PF on a HPC cluster we recommend installing from source to make full use of the optimized libraries installed on the cluster.

The process for installing deal.II using these three approaches on Mac OS X, Linux, and Windows is described below.

### Binary Package (Mac OS X)
A binary package for the latest version of deal.II (and other versions) is available here. The binary package includes all the libraries that you need (excluding CMake). To install deal.II open the .dmg file and drag ''deal.II.app'' into your Applications folder. Then, open ''deal.II.app'', which will open a Terminal window and set the necessary environment variables. You may need to modify your security settings so that your operating system will let you open the application (in System Preferences, select Security & Privacy, then under the general tab choose to allow apps downloaded from anywhere; don't click ''Open deal.II.app anyway'', we have found that doing so doesn't actually launch deal.II.app).

In some cases, deal.II.app will not open, possibly freezing at the ''Verifying...'' stage. If this happens, right click on it and select ''Show Package Contents''. Then go to ''/Contents/Resources/bin/'' and launch ''dealii-terminal''. This should launch the Terminal window that sets the environment variables.

### Binary Package (Linux)
Binary packages for deal.II are available for Gentoo, Debian (stretch or newer), Ubuntu (16.10 or newer), and Arch Linux. Please refer to https://www.dealii.org/download.html for details.

### Docker (Mac OS X, Linux, Windows)
Docker is a container platform that can be used to run code identically on any OS X, Linux, or Windows machine. The [PRISMS-PF Docker image](https://hub.docker.com/repository/docker/prismspf/prismspf/) contains a full instantiation of PRISMS-PF and all of its prerequisites. The steps to use it are given in on the [Installing PRISMS-PF page](#install_prismspf). Docker images for deal.II are also posted on [DockerHub](https://hub.docker.com/r/dealii/dealii/). These can be used directly, if so desired.

### Installing from Source (Mac OS X, Linux)
#### Manual Installation
Before building deal.II, you will need to install two libraries, MPI and p4est. Make sure you have both libraries installed before installing deal.II, or else you will have to reinstall deal.II! Many computing clusters already have these libraries installed. MPI libraries (such as MPICH2) are often available through package managers (e.g. yum, apt, pkg_add, port, brew), installing it using your package manager of choice is likely the simplest option. If not, you can find binaries and the source code at http://www.mpich.org/downloads/, as well as instructions for installation. The p4est library may also be available through your package manager. If not, you can download the latest release tarball at http://www.p4est.org/. To install, unzip the archive, move to the p4est root directory and enter the following commands at the command line:
```
$ ./configure --enable-mpi; make; make install
```
Users have found that this approach works better than the setup script ''p4est-setup.sh'' that is on the deal.ii website. More detailed installation instructions can be found in the ''INSTALL'' and ''README'' files in the p4est root directory.

With MPI and p4est installed, you are ready to install deal.II. Open the archive you downloaded from the deal.II website. From the root deal.II directory enter the following commands: 
```
$ mkdir build
$ cd build
$ cmake -DP4EST_DIR=/path/to/installation/dir -DDEAL_II_WITH_P4EST=ON
-DDEAL_II_WITH_MPI=ON -DCMAKE_INSTALL_PREFIX=/path/to/install/dir
$ make install
```

where ''/path/to/install/dir'' is the path to where you want to install deal.II. After installation is complete (which can take take up to an hour), open the file ''summary.log''. The second section of the log file lists the configured features. Ensure that DEAL_II_WITH_MPI and DEAL_II_WITH_P4EST say either ''set up with bundled packages'' or ''set up with external dependences''. If either is listed as ''OFF'', then deal.II was unable to find MPI or p4est. If this is the case, double check that MPI and p4est were installed correctly.

### Installation using candi
An alternative route is to use the [candi (compile and install)](https://github.com/koecher/candi) script to install deal.II and its prerequisites. If it works, it should be much simpler than the traditional route described above. While some users have reported success using this approach, it hasn't worked for everyone.
