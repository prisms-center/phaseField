# Tests {#tests}
The 'tests' folder in the root directory contains some tests to ensure that PRISMS-PF is working properly. A battery of automated tests can be found in the 'automated_tests' directory. These can be run by running the python script 'run_automatic_tests.py'. This will run a series of unit and regression tests, printing a summary of the results to the screen and writing to the file 'test_results.txt'. This python file can also be used as an example of automating PRISMS-PF simulations.

The 'unit_tests' directory contains a driver for a series of unit tests. To run them enter:
```
$ cmake .
$ make debug
$ mpirun -n 1 main
```
A summary of the tests will print to the screen.
