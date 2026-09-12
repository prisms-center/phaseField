# How to Contribute to PRISMS-PF
Thank you for taking an interest in contributing to PRISMS-PF! These are some notes that may be useful for contributors. Feel free to expand this document if you think there is anything else that should be known.

## Formatting Rules
PRISMS-PF has rules for text formatting. We format with [clang-format](https://clang.llvm.org/docs/ClangFormat.html). We recommend using an autoformatter with your developing environment. If you are using VSCode, installing the [Clang-Format extension](https://marketplace.visualstudio.com/items?itemName=xaver.clang-format) and turning on **format on save** in settings may be useful. You will also need a clang-format installation from your cli.

To ensure that your commits will pass our requirements, you can run pre-commit inside the repo, which will format the code properly and catch spelling errors. Your changes may need to be staged first.
```bash
# One-time installation:
python3 -m pip install pre_commit
# From any directory in the repo:
python3 -m pre_commit run --all-files
# Repeat until all-passing
```

## Adding an application
You can start to create a new application by copying an existing one. There are a couple of things you should make sure to do though.
First, make sure to rename the project in the new application's `CMakeLists.txt`.
```cmake
# ...
# Name the application
set(APPLICATION_NAME "my_app_name")
# ...
```
The other step is to add the relative path to your new application directory to the list at the top of [`applications/CMakeLists.txt`](/applications/CMakeLists.txt).
```cmake
set(
  # ...
  coupled_allen_cahn_cahn_hilliard
  mechanics/boundary_value_problem
  mechanics/eshelby_inclusion
  my_app_directory
  # ...
)
```

## Style Guidelines
TODO
<!--
Separate stuff into headers and source when appropriate.
Modify CMakeLists.txt in src.
Don't be narrow-minded when adding features. Try to keep things general.
USE IRREDUCIBLE REPRESENTATIONS.
Use c++ idioms.
Catch possible exceptions / do assertions with meaningful output (debug/release).
Template things where appropriate (never assume double, be SIMD-friendly too).
Do inline docs.
For apps, don't use templates if not appropriate (benchmarks for example).
For apps, don't hardcode things that should not be hard-coded.
-->
## Doxygen Documentation
TODO
<!--
Install doxygen.
Make sure submodule is present.
Nav to /docs.
Run cmake build.
Start server with python command.
-->
### Inline Docs
TODO
<!--
Add a brief description to things in headers.
-->
### Pages
TODO
<!--
Add .h files to dirs in doc/doxygen
Add them to DoxygenLayout.xml.in
-->
## Unit Tests
TODO
<!--
Add .cc files of tests.
How to create tests.
cmake instructions.
-->
