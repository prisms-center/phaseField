// ===========================================================================
// FUNCTION FOR INITIAL CONDITIONS
// ===========================================================================
#include <cmath>
#include <iostream>
using namespace std;

template <int dim, int degree>
void
customPDE<dim, degree>::setInitialCondition(const dealii::Point<dim> &p,
                                            const unsigned int        index,
                                            double                   &scalar_IC,
                                            dealii::Vector<double>   &vector_IC)
{
  // ---------------------------------------------------------------------
  // ENTER THE INITIAL CONDITIONS HERE FOR SCALAR FIELDS
  // ---------------------------------------------------------------------
  // Enter the function describing conditions for the fields at point "p".
  // Use "if" statements to set the initial condition for each variable
  // according to its variable index

  // The initial condition is a set of overlapping circles/spheres defined
  // by a hyperbolic tangent function. The center of each circle/sphere is
  // given by "center" and its radius is given by "radius".

  /*
   List of scalar variables
   Index: 0; var. nAnodic
   Index: 1; var. muAnodic
   Index: 2; var. nCathodic
   Index: 3; var. muCathodic
   Index: 4; var. psi
   Index: 5; var. mupsi
   Index: 6; var. Phi
   Index: 7; var. irxn
  */

  // xCoordinate of the point considered
  double posx = p[0];

  // yCoordinate of the point considered
  double posy = p[1];

  // zCoordinate of the point considered
  // double posz = p[2];

  const double SYSTEM_WIDTH  = userInputs.domain_size[0];
  const double SYSTEM_HEIGHT = userInputs.domain_size[1];

  double hyperTan = tanh((abs(posx - 1 * SYSTEM_WIDTH) - 1 * cathodeThickness) / deltaV);
  double height   = 0.5 * SYSTEM_HEIGHT;

  double orderParaCathodic = 0.5 * (1 - hyperTan);
  double orderParaAnodic   = 0.5 * (1 + hyperTan);
  double domainParaPsi     = 0.5 * (1 + tanh((posy - height) / deltaV));

  orderParaCathodic = orderParaCathodic * (1 - domainParaPsi);
  orderParaAnodic   = orderParaAnodic * (1 - domainParaPsi);

  if (index == 0)
    { // order parameter for the anodic phase, nAnodic

      scalar_IC = orderParaAnodic;
    }
  else if (index == 2)
    { // order parameter for the cathodic phase, nCathodic

      scalar_IC = orderParaCathodic;
    }
  else if (index == 4)
    { // electrolyte domain parameter, psi

      scalar_IC = domainParaPsi;
    }
  else if (index == 6)
    { // electrolyte potential Phi

      scalar_IC = initialGuessForPhi;
    }
  else
    { // all other variable fields

      scalar_IC = 0.0;
    }
  // ---------------------------------------------------------------------
}

// ===========================================================================
// FUNCTION FOR NON-UNIFORM DIRICHLET BOUNDARY CONDITIONS
// ===========================================================================

template <int dim, int degree>
void
customPDE<dim, degree>::setNonUniformDirichletBCs(const dealii::Point<dim> &p,
                                                  const unsigned int        index,
                                                  const unsigned int        direction,
                                                  const double              time,
                                                  double                   &scalar_BC,
                                                  dealii::Vector<double>   &vector_BC)
{
  // --------------------------------------------------------------------------
  // ENTER THE NON-UNIFORM DIRICHLET BOUNDARY CONDITIONS HERE
  // --------------------------------------------------------------------------
  // Enter the function describing conditions for the fields at point "p".
  // Use "if" statements to set the boundary condition for each variable
  // according to its variable index. This function can be left blank if there
  // are no non-uniform Dirichlet boundary conditions. For BCs that change in
  // time, you can access the current time through the variable "time". The
  // boundary index can be accessed via the variable "direction", which starts
  // at zero and uses the same order as the BC specification in parameters.in
  // (i.e. left = 0, right = 1, bottom = 2, top = 3, front = 4, back = 5).
  /*
if (index == 2){
  if (direction == 1){
      double x=p[0];
      double y=p[1];
      scalar_BC=sin(y/7.0);
  }
}
  */
  // -------------------------------------------------------------------------
}
