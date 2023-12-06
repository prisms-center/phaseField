---
:layout: default
---

[![PRISMS-PF Logo](../assets/logo.png)](https://prisms-center.github.io/phaseField/)

# Troubleshooting/FAQ

***
These are some of the most common issues users encounter while using PRISMS-PF 

## CMake Errors 

- **CMake cannot find deal.II**
  
  Typing <span style="font-family:Menlo">"cmake ."</span> gives you the following error
  
  > CMake Error at CMakeLists.txt:8 (FIND_PACKAGE):
  > Could not find a package configuration file provided by "deal.II"
  > (requested version 9.2.0) with any of the following names:
  
  > deal.IIConfig.cmake
  > deal.ii-config.cmake
  
  > Add the installation prefix of "deal.II" to CMAKE_PREFIX_PATH or set
  > "deal.II_DIR" to a directory containing one of the above files.  If
  > "deal.II" provides a separate development package or SDK, be sure it has
  > been installed
  >
  > -- Configuring incomplete, errors occurred! 
  
  Either deal.II has not been installed or the deal.II library has not been found. Make sure deal.II is installed and that the installation directory is included in the path variable. Type <span style="font-family:Menlo">"echo $PATH"</span> to see the directories currently on your path variable. Watch [this video](https://www.youtube.com/watch?v=7bBfigZsveM) for instructions on how to set the correct path.
  
  

## Compilation Errors

- **make:  No rule to make target `release'. Stop.**

  If you find this error when you type <span style="font-family:Menlo">"make release"</span> or <span style="font-family:Menlo">"make debug"</span>, it means a Makefile has not been found in the current directory. Please ensure that the CMake has executed with no errors.

## Runtime Errors

- Failure to run the code  

  Error

  >  HYDU_create_process (utils/launch/launch.c:73): execvp error on file main (No such file or directory)

​	This error can occur whenever the code is run in parallel with the line

​	$ mpirun -n <number of processors>  main

​	To fix it make sure to include "./" before "main" to indicate that the file main is in the current directory, e.g.

​	$ mpirun -n  4 ./main

​	
