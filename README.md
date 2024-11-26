![](logo_v2.png)

[![GitHub Linux](https://github.com/prisms-center/phaseField/actions/workflows/linux.yml/badge.svg)](https://github.com/prisms-center/phaseField/actions/workflows/linux.yml)
[![Clang-Format](https://github.com/prisms-center/phaseField/actions/workflows/clang-format.yml/badge.svg)](https://github.com/prisms-center/phaseField/actions/workflows/clang-format.yml)
[![Clang-Tidy](https://github.com/prisms-center/phaseField/actions/workflows/clang-tidy.yml/badge.svg)](https://github.com/prisms-center/phaseField/actions/workflows/clang-tidy.yml)

[![License: LGPL v2.1](https://img.shields.io/badge/License-lgpl-blue.svg)](https://www.gnu.org/licenses/lgpl-2.1)
[![DOI](https://zenodo.org/badge/22602327.svg)](https://zenodo.org/badge/latestdoi/22602327)

## Useful links:

[PRISMS-PF Website](https://prisms-center.github.io/phaseField/) <br>
[Code repository](https://github.com/prisms-center/phaseField) <br>
[User manual (with installation instructions)](https://prisms-center.github.io/phaseField/doxygen_files/manual.html) <br>
[User registration link](http://goo.gl/forms/GXo7Im8p2Y) <br>
[User forum](https://groups.google.com/forum/#!forum/prisms-pf-users) <br>
[Training slides/exercises](https://goo.gl/BBTkJ8) <br>
[PFHub phase-field community](https://pages.nist.gov/pfhub/)

## What is PRISMS-PF?

PRISMS-PF is a powerful, massively parallel finite element code for conducting phase field and other related simulations of microstructural evolution.  The phase field method is commonly used for predicting the evolution if microstructures under a wide range of conditions and material systems. PRISMS-PF provides a simple interface for solving customizable systems of partial differential equations of the type commonly found in phase field models, and has 24 pre-built application modules, including for precipitate evolution, grain growth, and solidification.

With PRISMS-PF, you have access to adaptive meshing and parallelization with near-ideal scaling for over a thousand processors. Moreover, the matrix-free framework from the deal.II library allows much larger than simulations than typical finite element programs – PRISMS-PF has been used for simulations with over one billion degrees of freedom. PRISMS-PF also provides performance competitive with or exceeding single-purpose codes. For example, even without enabling the mesh adaptivity features in PRISMS-PF, it has been demonstrated to be over 6x faster than an equivalent finite difference code.

This code is developed by the PRedictive Integrated Structural Materials Science (PRISMS) Center
at University of Michigan which is supported by the U.S. Department of Energy (DOE), Office of Basic Energy Sciences, Division of Materials Sciences and Engineering under Award #DE-SC0008637.

## Citing PRISMS-PF

Please cite [the following reference](https://www.nature.com/articles/s41524-020-0298-5) when discussing PRISMS-PF in a publication:

S. DeWitt, S. Rudraraju, D. Montiel, W.B. Andrews, and K. Thornton. PRISMS-PF: A general framework for phase-field modeling with a matrix-free finite element method. _npj Computuational Materials_ __6__, 29 (2020).

If additionally you would like to cite a specific release of PRISMS-PF, please use the following format:

PRISMS-PF, v2.4.0 (2024). Available from https://github.com/prisms-center/phaseField. DOI: 10.5281/zenodo.14026472.

For DOI information for other releases, please refer to [this site](https://doi.org/10.5281/zenodo.14026472).

## Quick start guide:

For detailed instructions on how to download and use PRISMS-PF, please consult the [PRISMS-PF User Manual](https://prisms-center.github.io/phaseField/doxygen_files/manual.html). A (very) abbreviated version of the instructions is given below.

### Install:

Install CMake, p4est, and deal.II (version 9.6 recommended).

Clone the PRISMS-PF GitHub repository and navigate its folder.
```bash
git clone https://github.com/prisms-center/phaseField.git
cd phaseField
```
Configure and compile the main library.
```bash
cmake . && make -j <nprocs>
```
here `<nprocs>` denotes the number of threads you want to use to compile the library.

### Running a pre-built application:

Please refer to the [Running a PRISMS-PF Example App](https://prisms-center.github.io/phaseField/doxygen_files/running_apps.html) for full details including instructions for visualization of the results.

Examples of various phase field models are located under the
applications directory. The easiest way to get started on the code is to
run the example apps in this folder.

The example apps are intended to serve as (1) Demonstration of the
capabilities of this library, (2) Provide a framework for
further development of specialized/advanced applications by
users.

Entering the following commands will run one of the pre-built example applications (the Cahn-Hilliard spinodal decomposition application in this case):
```bash
cd applications/cahnHilliard
cmake .
make -j <nprocs>
```
This will generate two executable files: `main` and `main-debug`. Debug and release are compiler configurations. Debug mode is slower, but contains less optimiziations and more meaningful error messages. This makes it ideal for application/model code development. Release mode has less "safety features" and meaningful error messages, with more optimizations (faster runtime).

Debug execution (serial runs):
```bash
$ ./main-debug
```
Release execution (parallel runs):
```bash
$ mpirun -np <nprocs> ./main
```

### Visualization:

Output of the primal fields fields is in standard vtk
format (parallel:*.pvtu, serial:*.vtu files) which can be visualized with the
following open source applications:

1. VisIt (https://wci.llnl.gov/simulation/computer-codes/visit/downloads)
2. Paraview (http://www.paraview.org/download/)

## Version information:

This version of the code, v2.4, contains moderate changes from v2.3. It was released in November 2024. See [version_changes.md](version_changes.md) for details.

## License:

GNU Lesser General Public License (LGPL). Please see the file
LICENSE for details.

## Further information, questions, issues and bugs:

+ prisms-pf-users@googlegroups.com (user forum)
+ prisms-pf@umich.edu  (developer email list)
