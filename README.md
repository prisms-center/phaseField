![](logo_v2.png)

[![GitHub Linux](https://github.com/prisms-center/phaseField/actions/workflows/linux.yml/badge.svg)](https://github.com/prisms-center/phaseField/actions/workflows/linux.yml)
[![Check Links](https://github.com/prisms-center/phaseField/actions/workflows/links.yml/badge.svg)](https://github.com/prisms-center/phaseField/actions/workflows/links.yml)

[![License: LGPL v2.1](https://img.shields.io/badge/License-lgpl-blue.svg)](https://www.gnu.org/licenses/lgpl-2.1)

## Useful links:

[Website](https://prisms-center.github.io/phaseField/) \
[Code](https://github.com/prisms-center/phaseField) \
[Manual](https://prisms-center.github.io/phaseField/doxygen/) \
[Discussions](https://github.com/prisms-center/phaseField/discussions) \
[Training slides](https://github.com/prisms-center/PRISMS-PF_Training_Materials) \
[PFHub phase-field community](https://pages.nist.gov/pfhub/) \
[DOI](https://zenodo.org/badge/latestdoi/22602327)

## What is PRISMS-PF?

PRISMS-PF is a powerful, massively parallel finite element code for conducting phase field and other related simulations of microstructural evolution.  The phase field method is commonly used for predicting the evolution in microstructures under a wide range of conditions and material systems. PRISMS-PF provides a simple interface for solving customizable systems of partial differential equations of the type commonly found in phase field models, and has pre-built application modules, including precipitate evolution, grain growth, and solidification.

With PRISMS-PF, you have access to adaptive meshing and parallelization with near-ideal scaling for over a thousand processors. Moreover, the matrix-free framework from the deal.II library allows much larger than simulations than typical finite element programs – PRISMS-PF has been used for simulations with over one billion degrees of freedom. PRISMS-PF also provides performance competitive with, and often exceeding single-purpose codes.

This code is developed by the PRedictive Integrated Structural Materials Science (PRISMS) Center at University of Michigan which is supported by the U.S. Department of Energy (DOE), Office of Basic Energy Sciences, Division of Materials Sciences and Engineering under Award #DE-SC0008637.

## Quick start guide

For detailed install and usage instructions check out the [Manual](https://prisms-center.github.io/phaseField/doxygen/). A (very) abbreviated version is given below.

### Configure and install

Install the required dependencies:
- CMake [3.25+]
- LAPACK
- MPI
- p4est
- deal.II [9.6.0+]

Configure, compile, and install the main library. We recommend setting the `PRISMS_PF_DIR` variable into your rc file so you have permanent access to it.
```
export PRISMS_PF_DIR="path/to/where/you/want/to/install"
```

```
cmake -B build
cmake --build build -j <nprocs>
cmake --install build --prefix=$PRISMS_PF_DIR
```

Here `<nprocs>` denotes the number of threads you want to use to compile the library.

For access to our executable utilities, you should also add `bin` to path. You can also add this into your rc file to have permanent access to it.
```
export PATH=$PATH:$PRISMS_PF_DIR/bin
```

### Running an application

After building the core library, our example applications can be compiled with similar commands. Make sure you navigate to an application before executing these commands.
```
cmake -DCMAKE_BUILD_TYPE=<type> -B build
cmake --build build -j <nprocs>
```

Here `<type>` is the build type. There are two options `Debug` and `Release`. `Release` contains more optimizations than `Debug`, but often has less meaningful error messages.

Executables can be run in serial or parallel.
```
build/main
```
```
mpirun -n <nprocs> build/main
```

### Visualization

Outputs are most commonly visualized in the following open source applications:
- VisIt (https://visit-dav.github.io/visit-website/)
- Paraview (https://www.paraview.org/)

## Citing PRISMS-PF

Please cite [the following reference](https://www.nature.com/articles/s41524-020-0298-5) when discussing PRISMS-PF in a publication:

S. DeWitt, S. Rudraraju, D. Montiel, W.B. Andrews, and K. Thornton. PRISMS-PF: A general framework for phase-field modeling with a matrix-free finite element method. _npj Computuational Materials_ __6__, 29 (2020).

If you would like to cite a specific release of PRISMS-PF, please use the following format:

PRISMS-PF, v3.0.0 (2026). Available from https://github.com/prisms-center/phaseField. DOI: 10.5281/zenodo.19261983.

For DOI information of other releases, refer to [this site](https://doi.org/10.5281/zenodo.2583307).

## Contributors
Thanks to everyone who has contributed to the project!

<a href="https://github.com/prisms-center/phaseField/graphs/contributors">
  <img src="https://contrib.rocks/image?repo=prisms-center/phaseField" />
</a>

Made with [contrib.rocks](https://contrib.rocks).

## License

Please see the LICENSE file for details.
