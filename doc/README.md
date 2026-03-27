# PRISMS-PF Documentation
The prebuilt documentation can be found here:

https://prisms-center.github.io/phaseField/doxygen/manual.html

## Pull Git Submodules
Make sure Doxygen Awesome has been pulled.
```
git submodule update --init --recursive
```

## Install Doxygen
Doxygen can be installed from https://www.doxygen.nl/

## Build
After Doxygen is installed, simply compile like normal.
```
cmake .
```
```
make
```

## Visualize
Launch a local web server using Python
```
python -m http.server
```

## Licensing
The documentation is covered with the same license as the PRISMS-PF library.
