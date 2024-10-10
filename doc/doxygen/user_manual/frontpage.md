# PRISMS-PF User Manual {#manual}

## PRISMS-PF Overview
PRISMS-PF is an open-source finite element simulation framework with a focus on solving equations from phase field models. It is being developed as part of the PRISMS Center at the University of Michigan. PRISMS-PF is built on top of the popular and powerful deal.II finite element library. PRISMS-PF was developed to make high performance phase field simulation capabilities available to the wider scientific community with an easy-to-use and easy-to-learn interface. Benchmark tests show that even without its adaptive meshing capabilities PRISMS-PF is competitive with specialized finite difference codes written in Fortran with MPI parallelization in terms performance,  and in some cases runs several times faster. PRISMS-PF can be used to solve a wide variety of systems of coupled time-dependent partial differential equations (e.g. the diffusion equation) and time-independent partial differential equations (e.g. Poisson's equation). With PRISMS-PF, you have access to adaptive meshing and parallelization that has been demonstrated on over a thousand processors. Moreover, the matrix-free framework from deal.II allows much larger than typical finite element programs -- PRISMS-PF has been used for calculations with over one \emph{billion} degrees of freedom.

This user guide starts with instructions for downloading PRISMS-PF and install the necessary prerequisites. Next are instructions for running the built-in example applications and visualizing the results. A detailed look at the input files follows, and the guide finishes with some instructions on how to create your own PRISMS-PF applications.

If you run into any issues not covered by this guide, or have any questions, please send an email to the PRISMS-PF email list: prisms-pf-users@googlegroups.com. If you'd like to contact the PRISMS-PF developers directly, send an email to: prismsphaseField.dev@umich.edu.

PRISMS-PF was developed as part of the PRedictive Integrated Structural Materials Science (PRISMS) Center at University of Michigan which is supported by the U.S. Department of Energy (DOE), Office of Basic Energy Sciences, Division of Materials Sciences and Engineering under Award \#DE-SC0008637.

## Structure of the PRISMS-PF User Guide
@subpage install

@subpage running_apps

@subpage app_structure

@subpage input_file

@subpage app_files

@subpage tests

@subpage custom_apps
