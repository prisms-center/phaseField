/** \page structure Structure of a PRISMS-PF application

Welcome to the PRISMS-PF application structure tutorial. This page explains the four main
components you need to set up a simulation using PRISMS-PF and how to put them together.
It will cover how to declare [fields](#fields), how to set up [solvers](#solvers) for
those fields, where to write your [equations](#equations), and where to set PRISMS-PF
[settings](#user_inputs). It will then cover how to put these together to solve a PDE
[problem](#problem).

\subsection fields Fields
The first step to setting up a simulation in PRISMS-PF is declaring all the fields you
wish to evolve. This is done by constructing an array (a C++ std::vector) containing
FieldAttributes objects. Each FieldAttributes object has two important attributes: `name`,
and `field_type`. The `name` must be a unique alphanumeric identifier for each field. The
`field_type` is the tensor rank of the associated field. It is `Scalar` by default, but
can also be `Vector`. Importantly, be sure to take note of the position of each field in
the array, this is what we will refer to as its "field index" (c++ is zero-indexed).

```c++
std::vector<FieldAttrubutes> fields
= {
    FieldAttributes("c", Scalar), // index 0
    FieldAttributes("u", Vector), // index 1...
    FieldAttributes("phi"), // Scalar by default
    FieldAttributes("T")
  };
```

\subsection solvers Solvers
The next step is to declare how you wish to solve each field and in what order groups
of fields are solved. This is done by constructing an array (a C++ std::vector) containing
SolveBlock objects. Each SolveBlock object has six important attributes:

- `id`: A unique integer to identify this solver in equations.
- `solve_type`: One of the following solve methods: `Explicit`, `Linear`, `Newton`, or
  `Constant`. See \ref SolveType for details.
- `solve_timing`: When to solve each field using the equations.
  - `Initialized`: On increment 0, uses the initial conditions functionality to explicitly
    initialize the fields. Solves normally otherwise.
  - `Uninitialized`: Solves normally with the equations on every increment.
  - `PostProcess`: Only solves on output increments.
  - `NucleationRate`: Only solves on nucleation increments and output increments.
- `field_indices`: A set of the field indices of the fields that will be solved by this
  solver.
- `dependencies_rhs` and `dependencies_lhs`: A map that associates each field needed for
the equations on the rhs or lhs side of the solve with the types of dependency
  (`value`, `gradient`, `old`, `trial`). This can be quickly populated using the
  function `make_dependency_set`.

\subsection equations Equations
To write PDE equations in PRISMS-PF, one must implement methods of a derived class of
`PDEOperatorBase`. In the provided examples, this class is known as `CustomPDE`.
This is also where you should declare your model constants and settings.
There are four virtual functions that you may override to implement your equations:

- `equations_rhs`: Equations for the rhs of explicit and implicit solves. See \ref
  Equations for details.
- `equations_lhs`: Equations for the lhs of implicit equations. See \ref Equations for
  details.
- `set_initial_condition`: Function to set initial conditions by equations for fields that
  need them.
- `set_dirichlet`: Function to set Dirichlet boundary conditions for fields
that need them.

\subsection user_inputs PRISMS-PF Settings
PRISMS-PF simulation settings are set up by constructing a `UserInputParameters` object.
This contains settings for your triangulation, time-stepping, adaptive meshing, simulation
outputs, and more. You can set these settings manually in your code, or by reading in a
parameters file. Detailed documentation of options in `UserInputParameters`.
```c++
// Reading in settings from a parameter file.
UserInputParameters<dim> user_inputs("parameters.prm");
```

\subsection problem Run a simulation

To run a simulation with PRISMS-PF, create a `Problem` object from your fields, solve
groups, and, custom PDE operator. You will also need to provide another object called
PFTools, though you can leave this empty unless you are using the /ref nucleation
functions of PRISMS-PF. Lastly, call problem.run() to execute the simulation.

```c++
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
```
*/
