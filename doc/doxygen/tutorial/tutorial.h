/** \page tutorial Tutorial

Welcome to the tutorial pages!

Here we will give a recap of our quick start guide, showing how to build and run a
PRISMS-PF application and make basic changes to suit your needs. After reviewing
this page, see the next page on the \ref structure "structure of a PRISMS-PF application".
For developers, please check out our \ref dev_docs "developer's guide" to learn how to
build this documentation locally and contribute to PRISMS-PF.

\section getting_started Getting Started
Before starting, you will need to install PRISMS-PF and its dependencies.
For more information, see our \ref installation page.

Browse our example \ref applications and select one that solves a problem similar
to your problem of interest. For this tutorial, we will use the explicit Allen-Cahn
application.

Each application includes a minimum of three components. The first is a subclass
of \ref PDEOperatorBase, which contains the equations, initial conditions, and boundary
conditions (called `CustomPDE` in the Allen-Cahn app, in `custom_pde.h`). The second
part is a `main.cc` file containing the main function, which sets up the user inputs,
instantiates CustomPDE, creates the \ref Problem, and invokes `Problem.solve()`.
The final part is a parameters file, typically called `parameters.prm`, which allows
specifying different conditions under which to run the simulation. Our example
applications also come with a markdown file containing a description of the equations
being solved and the derivation of the
[weak form](https://en.wikipedia.org/wiki/Weak_formulation).

To run an application, you must first compile the code ("build the application"). To do
this, navigate to the application's directory in your terminal and invoke cmake:
@code{bash}
cmake . -DCMAKE_BUILD_TYPE=Debug
make
@endcode
After compiling the application, you can run in serial (one processor) by simply invoking
the `main` executable and specifying the parameters file using the `-i` option flag:
@code{bash}
./main -i parameters.prm
@endcode
To run with multiple processors using MPI, use `mpirun`:
@code{bash}
mpirun -n 16 ./main -i parameters.prm
@endcode

By default, the output files are saved in a directory called `solutions` in the same
directory as the executable. You can change this and many other settings by editing the
parameters file. See the \ref parameters page for more information on the available
settings. For example, you can change the size of the finite element mesh and the
simulation time step in the parameters file.

To implement your own app, you can start by copying the folder of an application
solving a problem similar to yours (or copy the blank application). Then, edit the `main.cc`
file to set up fields, solve blocks (a.k.a. solvers), and equations. For more information
about the structure of an application, see
\ref structure "Structure of a PRISMS-PF application".

*/

/** \page structure Structure of a PRISMS-PF application

Welcome to the PRISMS-PF application structure tutorial. This page explains the four main
components you need to set up a simulation using PRISMS-PF and how to put them together.
It will cover how to declare \ref fields "fields", how to set up \ref solvers "solvers"
for those fields, where to write your \ref equations "equations", and where to set
PRISMS-PF \ref user_settings "settings". It will then cover how to put these together to
solve a PDE \ref problem "problem".

\subsection fields Fields

The first step to setting up a simulation in PRISMS-PF is declaring all the fields you
wish to evolve. This is done by constructing an array (a C++ `std::vector`) containing
FieldAttributes objects. Each FieldAttributes object has two important attributes:
`name` and `field_type`. The `name` must be a unique alphanumeric identifier for each
field. The `field_type` is the tensor rank of the associated field. It is `Scalar` by
default, but can also be `Vector`. Be sure to take note of the position of each field in
the array, this is what we will refer to as its "field index" (C++ is zero-indexed).

Example:
\code{.cpp}
std::vector<FieldAttributes> fields
= {
    FieldAttributes("c", Scalar), // index 0
    FieldAttributes("u", Vector), // index 1
    FieldAttributes("phi"), // Scalar by default
    FieldAttributes("T")
  };
\endcode

\subsection solvers Solvers

The next step is to declare how you wish to solve each field and in what order groups
of fields are solved. This is done by constructing an array (a C++ `std::vector`)
containing SolveBlock objects. Each SolveBlock object has six important attributes:

\li \c id: A unique integer to identify this solver in equations.
\li \c solve_type: One of the following solve methods: \c Explicit, \c Linear, \c Newton,
or \c Constant. See \ref SolveType for details.
\li \c solve_timing: When to solve each field using the equations.
\par
<ul>
<li> \c Initialized: On increment 0, uses the initial conditions functionality to
explicitly initialize the fields. Solves normally otherwise.
<li> \c Uninitialized: Solves normally with the equations on every increment.
<li> \c PostProcess: Only solves on output increments.
<li> \c NucleationRate: Only solves on nucleation increments and output increments.
</ul>
\li \c field_indices: A set of the field indices of the fields that will be solved by this
solver.
\li \c dependencies_rhs and \c dependencies_lhs: A map that associates each field needed
for the equations on the rhs or lhs side of the solve with the types of dependency (\c
value, \c gradient, \c old, \c trial). This can be quickly populated using the function \c
make_dependency_set.

\subsection equations Equations

To write PDE equations in PRISMS-PF, one must implement methods of a derived class of
PDEOperatorBase. In the provided examples, this class is known as `CustomPDE`. This is
also where you should declare your model constants and settings. There are four virtual
functions that you may override to implement your equations:

\li \c equations_rhs: Equations for the rhs of explicit and implicit solves. See
equations for details.
\li \c equations_lhs: Equations for the lhs of implicit equations. See equations for
details.
\li \c set_initial_condition: Function to set initial conditions by equations for fields
that need them.
\li \c set_dirichlet: Function to set Dirichlet boundary conditions for fields that need
them.

\subsection user_settings PRISMS-PF Settings

PRISMS-PF simulation settings are set up by constructing a UserInputParameters object.
This contains settings for your triangulation, time-stepping, adaptive meshing, simulation
outputs, and more. You can set these settings manually in your code, or by reading in a
parameters file. Detailed documentation of options in UserInputParameters.

Example:
\code{.cpp}
// Reading in settings from a parameter file.
UserInputParameters<dim> user_inputs("parameters.prm");
\endcode

\subsection problem Run a simulation

To run a simulation with PRISMS-PF, create a Problem object from your fields, solve
groups, and, custom PDE operator. You will also need to provide another object called
PFTools, though you can leave this empty unless you are using the /ref nucleation
functions of PRISMS-PF. Lastly, call `problem.run()` to execute the simulation.

\code{.cpp}
std::vector<FieldAttributes> fields({...}); // Fill in
std::vector<SolveBlock> solvers({...}); // Fill in
CustomPDE<dim, degree, double> pde_operator(user_inputs, pf_tools);
UserInputParameters<dim>       user_inputs(parameters_filename);
PhaseFieldTools<dim>           pf_tools;

Problem<dim, degree, double>   problem(fields,
                                       solvers,
                                       user_inputs,
                                       pf_tools,
                                       pde_operator);
problem.run();
\endcode

*/

/** \page dev_docs Developers

There's a lot of institutional knowledge that comes along with being a developer of
PRISMS-PF. Here, we hope to cover all the bits and pieces that you'll need to adhere to
our developer standards.

\subsection Documentation

There are a few components of our documentation:
- Doxygen-style inline code documentation (core library and applications)
- Doxygen-style manual documentation
- GitHub pages website

Both the inline code documentation and manual documentation are bundled together in what
we call the PRISMS-PF User Manual. As a developer (and user) this is what you'll be using
and changing the most. To visualize the User Manual locally you'll need Doxygen and
Graphviz. On Ubuntu these can be installed with
@code
sudo apt-get install doxygen graphviz
@endcode
Additionally, we use Doxygen Awesome as a submodule to make everything look pretty. The
submodules can be downloaded with
@code
git submodule update --init --recursive
@endcode
@warning If you don't have Doxygen Awesome the website will not render correctly.

With all the prerequisites, we just need to build the documentation. Newer versions of
CMake can use the `CMakePresets.json`
@code
cmake --workflow --preset docs
@endcode
or without the preset
@code
cmake -DPRISMS_PF_DOCS=ON -DCMAKE_BUILD_TYPE=Debug -B build
cmake --build build/docs --target doc -j <nprocs>
@endcode

Once the build is complete you can visualize the local documentation by navigating to the
doxygen folder in the build directory and starting a python HTTP server.
@code
cd build/docs/doc/<branch_name>
python3 -m http.server 8000
@endcode
You can open the `http://localhost:8000/` on whatever web browser you use.

@note Some features of the User Manual don't render locally (e.g., the version selector).


 */
