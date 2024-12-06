
#ifndef PExtern_HH
#define PExtern_HH

#include "../PField.hh"
#include "../PFunction.hh"
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <vector>

// In future, might have more complicated OutType,
//   so make all have 'void' return and pass everything by reference

extern "C"
{
  // Functions for using a PSimpleBase externally (say Python or Fortran)
  //   written for VarContainer=double*, OutType=double, hence 'dsd' in function names

  void
  PSimpleFunction_dsd_new(char *name, PRISMS::PSimpleBase<double *, double> *&f);

  void
  PSimpleFunction_dsd_delete(PRISMS::PSimpleBase<double *, double> *&f);

  void
  PSimpleFunction_dsd_name(PRISMS::PSimpleBase<double *, double> *f, char *name);

  void
  PSimpleFunction_dsd_calc(PRISMS::PSimpleBase<double *, double> *f,
                           double                                *var,
                           double                                &val);

  void
  PSimpleFunction_dsd_get(PRISMS::PSimpleBase<double *, double> *f, double &val);

  // Functions for using a PSimpleBase externally (say Python or Fortran)
  //   written for VarContainer=double, OutType=double, hence 'dd' in function names

  void
  PSimpleFunction_dd_new(char *name, PRISMS::PSimpleBase<double, double> *&f);

  void
  PSimpleFunction_dd_delete(PRISMS::PSimpleBase<double, double> *&f);

  void
  PSimpleFunction_dd_name(PRISMS::PSimpleBase<double, double> *f, char *name);

  void
  PSimpleFunction_dd_calc(PRISMS::PSimpleBase<double, double> *f,
                          double                               var,
                          double                              &val);

  void
  PSimpleFunction_dd_get(PRISMS::PSimpleBase<double, double> *f, double &val);

  // Functions for using a PFuncBase externally (say Python or Fortran)
  //   written for VarContainer=double*, OutType=double, hence 'dsd' in function names

  void
  PFunction_dsd_new(char *name, PRISMS::PFuncBase<double *, double> *&f);

  void
  PFunction_dsd_delete(PRISMS::PFuncBase<double *, double> *&f);

  void
  PFunction_dsd_name(PRISMS::PFuncBase<double *, double> *f, char *name);

  void
  PFunction_dsd_size(PRISMS::PFuncBase<double *, double> *f, int &size);

  void
  PFunction_dsd_var_name(PRISMS::PFuncBase<double *, double> *f, int i, char *var_name);

  void
  PFunction_dsd_var_description(PRISMS::PFuncBase<double *, double> *f,
                                int                                  i,
                                char                                *var_description);

  // void PFunction_dsd_simplefunc(PRISMS::PFuncBase<double*,double>* f,
  // PSimpleBase<double*, double>* &simplefunc); void
  // PFunction_dsd_grad_simplefunc(PRISMS::PFuncBase<double*,double>* f, int di,
  // PSimpleBase<double*, double>* &simplefunc); void
  // PFunction_dsd_hess_simplefunc(PRISMS::PFuncBase<double*,double>* f, int di, int dj,
  // PSimpleBase<double*, double>* &simplefunc);

  void
  PFunction_dsd_calc(PRISMS::PFuncBase<double *, double> *f, double *var, double &val);

  void
  PFunction_dsd_calc_grad(PRISMS::PFuncBase<double *, double> *f,
                          double                              *var,
                          int                                  di,
                          double                              &val);

  void
  PFunction_dsd_calc_hess(PRISMS::PFuncBase<double *, double> *f,
                          double                              *var,
                          int                                  di,
                          int                                  dj,
                          double                              &val);

  void
  PFunction_dsd_eval(PRISMS::PFuncBase<double *, double> *f, double *var);

  void
  PFunction_dsd_eval_grad(PRISMS::PFuncBase<double *, double> *f, double *var, int di);

  void
  PFunction_dsd_eval_hess(PRISMS::PFuncBase<double *, double> *f,
                          double                              *var,
                          int                                  di,
                          int                                  dj);

  void
  PFunction_dsd_get(PRISMS::PFuncBase<double *, double> *f, double &val);

  void
  PFunction_dsd_get_grad(PRISMS::PFuncBase<double *, double> *f, int di, double &val);

  void
  PFunction_dsd_get_hess(PRISMS::PFuncBase<double *, double> *f,
                         int                                  di,
                         int                                  dj,
                         double                              &val);

  // Functions for using a PBasisSetBase externally (say Python or Fortran)
  //   written for InType=double, OutType=double, hence 'dd' in function names

  void
  PBasisSet_dd_new(char *name, PRISMS::PBasisSetBase<double, double> *&b, int N);

  void
  PBasisSet_dd_delete(PRISMS::PBasisSetBase<double, double> *&b);

  void
  PBasisSet_dd_name(PRISMS::PBasisSetBase<double, double> *b, char *name);

  void
  PBasisSet_dd_description(PRISMS::PBasisSetBase<double, double> *b, char *description);

  void
  PBasisSet_dd_size(PRISMS::PBasisSetBase<double, double> *b, int &size);

  void
  PBasisSet_dd_resize(PRISMS::PBasisSetBase<double, double> *b, int N);

  void
  PBasisSet_dd_max_size(PRISMS::PBasisSetBase<double, double> *b, int &max_size);

  // void PBasisSet_dd_basis_function(PRISMS::PFuncBase<double*,double>* b, int term,
  // PFuncBase<double,double>* &f);

  void
  PBasisSet_dd_calc(PRISMS::PBasisSetBase<double, double> *b,
                    int                                    term,
                    double                                 var,
                    double                                &val);

  void
  PBasisSet_dd_calc_grad(PRISMS::PBasisSetBase<double, double> *b,
                         int                                    term,
                         double                                 var,
                         double                                &val);

  void
  PBasisSet_dd_calc_hess(PRISMS::PBasisSetBase<double, double> *b,
                         int                                    term,
                         double                                 var,
                         double                                &val);

  void
  PBasisSet_dd_eval(PRISMS::PBasisSetBase<double, double> *b, double var);

  void
  PBasisSet_dd_eval_grad(PRISMS::PBasisSetBase<double, double> *b, double var);

  void
  PBasisSet_dd_eval_hess(PRISMS::PBasisSetBase<double, double> *b, double var);

  void
  PBasisSet_dd_get(PRISMS::PBasisSetBase<double, double> *b, int term, double &val);

  void
  PBasisSet_dd_get_grad(PRISMS::PBasisSetBase<double, double> *b, int term, double &val);

  void
  PBasisSet_dd_get_hess(PRISMS::PBasisSetBase<double, double> *b, int term, double &val);

  // void PBasisSet_dd_getall(PRISMS::PBasisSetBase<double,double>* b, const double*
  // &val); void PBasisSet_dd_getall_grad(PRISMS::PBasisSetBase<double,double>* b, const
  // double* &val); void PBasisSet_dd_getall_hess(PRISMS::PBasisSetBase<double,double>* b,
  // const double* &val);

  // Functions for using a PSeriesFunction externally (say Python or Fortran)
  //   written for InType=double, OutType=double, VarContainer=double*,
  //   IndexContainer=int*, hence 'dsis'

  // Functions for using a PSeriesFunction externally (say Python or Fortran)
  //   written for InType=double, OutType=double, VarContainer=double*,
  //   IndexContainer=int*, hence 'dsis'

  void
  PSeriesFunction_dsis_new(PRISMS::PSeriesFunction<double, double, double *, int *> *&f);

  void
  PSeriesFunction_dsis_setnew(
    PRISMS::PSeriesFunction<double, double, double *, int *> *&f,
    PRISMS::PBasisSetBase<double, double>                    **basis_set,
    int                                                        order);

  void
  PSeriesFunction_dsis_delete(
    PRISMS::PSeriesFunction<double, double, double *, int *> *&f);

  void
  PSeriesFunction_dsis_clear(PRISMS::PSeriesFunction<double, double, double *, int *> *f);

  void
  PSeriesFunction_dsis_set(PRISMS::PSeriesFunction<double, double, double *, int *> *f,
                           PRISMS::PBasisSetBase<double, double> **basis_set,
                           int                                     order);

  void
  PSeriesFunction_dsis_order(PRISMS::PSeriesFunction<double, double, double *, int *> *f,
                             int &order);

  void
  PSeriesFunction_dsis_volume(PRISMS::PSeriesFunction<double, double, double *, int *> *f,
                              int &volume);

  void
  PSeriesFunction_dsis_dim(PRISMS::PSeriesFunction<double, double, double *, int *> *f,
                           int                                                       i,
                           int                                                      &dim);

  void
  PSeriesFunction_dsis_get_linear_coeff(
    PRISMS::PSeriesFunction<double, double, double *, int *> *f,
    int                                                       i,
    double                                                   &coeff);

  void
  PSeriesFunction_dsis_get_tensor_coeff(
    PRISMS::PSeriesFunction<double, double, double *, int *> *f,
    int                                                      *term,
    double                                                   &coeff);

  void
  PSeriesFunction_dsis_set_linear_coeff(
    PRISMS::PSeriesFunction<double, double, double *, int *> *f,
    int                                                       i,
    double                                                    coeff);

  void
  PSeriesFunction_dsis_set_tensor_coeff(
    PRISMS::PSeriesFunction<double, double, double *, int *> *f,
    int                                                      *term,
    double                                                    coeff);

  void
  PSeriesFunction_dsis_linear_index(
    PRISMS::PSeriesFunction<double, double, double *, int *> *f,
    int                                                      *term,
    int                                                      &linear_index);

  void
  PSeriesFunction_dsis_tensor_indices(
    PRISMS::PSeriesFunction<double, double, double *, int *> *f,
    int                                                       linear_index,
    int                                                      *term);

  // ----------------------------------------------------------
  // Use these functions if you want to evaluate a single value

  void
  PSeriesFunction_dsis_calc(PRISMS::PSeriesFunction<double, double, double *, int *> *f,
                            double                                                   *var,
                            double &val);

  void
  PSeriesFunction_dsis_calc_grad(
    PRISMS::PSeriesFunction<double, double, double *, int *> *f,
    double                                                   *var,
    int                                                       di,
    double                                                   &val);

  void
  PSeriesFunction_dsis_calc_hess(
    PRISMS::PSeriesFunction<double, double, double *, int *> *f,
    double                                                   *var,
    int                                                       di,
    int                                                       dj,
    double                                                   &val);

  // ----------------------------------------------------------
  // Use these functions to evaluate several values, then use 'get' methods to access
  // results

  void
  PSeriesFunction_dsis_eval(PRISMS::PSeriesFunction<double, double, double *, int *> *f,
                            double *var);

  void
  PSeriesFunction_dsis_eval_grad(
    PRISMS::PSeriesFunction<double, double, double *, int *> *f,
    double                                                   *var);

  void
  PSeriesFunction_dsis_eval_hess(
    PRISMS::PSeriesFunction<double, double, double *, int *> *f,
    double                                                   *var);

  void
  PSeriesFunction_dsis_get(PRISMS::PSeriesFunction<double, double, double *, int *> *f,
                           double                                                   &val);

  void
  PSeriesFunction_dsis_get_grad(
    PRISMS::PSeriesFunction<double, double, double *, int *> *f,
    int                                                       di,
    double                                                   &val);

  void
  PSeriesFunction_dsis_get_hess(
    PRISMS::PSeriesFunction<double, double, double *, int *> *f,
    int                                                       di,
    int                                                       dj,
    double                                                   &val);

  // ----------------------------------------------------------
  // Functions for evaluating basis functions & their derivatives:

  // Use these functions if you want to evaluate a single value

  // use basis index and term index for individual basis function

  void
  PSeriesFunction_dsis_calc_basis(
    PRISMS::PSeriesFunction<double, double, double *, int *> *f,
    int                                                       bindex,
    int                                                       term,
    double                                                   *var,
    double                                                   &val);

  void
  PSeriesFunction_dsis_calc_basis_grad(
    PRISMS::PSeriesFunction<double, double, double *, int *> *f,
    int                                                       bindex,
    int                                                       term,
    double                                                   *var,
    double                                                   &val);

  void
  PSeriesFunction_dsis_calc_basis_hess(
    PRISMS::PSeriesFunction<double, double, double *, int *> *f,
    int                                                       bindex,
    int                                                       term,
    double                                                   *var,
    double                                                   &val);

  // or use tensor indices to evaluate basis function product
  void
  PSeriesFunction_dsis_calc_tensor_basis(
    PRISMS::PSeriesFunction<double, double, double *, int *> *f,
    int                                                      *term,
    double                                                   *var,
    double                                                   &val);

  void
  PSeriesFunction_dsis_calc_tensor_basis_grad(
    PRISMS::PSeriesFunction<double, double, double *, int *> *f,
    int                                                      *term,
    double                                                   *var,
    int                                                       di,
    double                                                   &val);

  void
  PSeriesFunction_dsis_calc_tensor_basis_hess(
    PRISMS::PSeriesFunction<double, double, double *, int *> *f,
    int                                                      *term,
    double                                                   *var,
    int                                                       di,
    int                                                       dj,
    double                                                   &val);

  // ----------------------------------------------------------
  // Use these functions to evaluate all basis functions,
  //   then use following methods to access results.

  void
  PSeriesFunction_dsis_eval_basis_all(
    PRISMS::PSeriesFunction<double, double, double *, int *> *f,
    double                                                   *var);

  void
  PSeriesFunction_dsis_eval_basis(
    PRISMS::PSeriesFunction<double, double, double *, int *> *f,
    double                                                   *var,
    int                                                       i);

  void
  PSeriesFunction_dsis_eval_basis_grad_all(
    PRISMS::PSeriesFunction<double, double, double *, int *> *f,
    double                                                   *var);

  void
  PSeriesFunction_dsis_eval_basis_grad(
    PRISMS::PSeriesFunction<double, double, double *, int *> *f,
    double                                                   *var,
    int                                                       i);

  void
  PSeriesFunction_dsis_eval_basis_hess_all(
    PRISMS::PSeriesFunction<double, double, double *, int *> *f,
    double                                                   *var);

  void
  PSeriesFunction_dsis_eval_basis_hess(
    PRISMS::PSeriesFunction<double, double, double *, int *> *f,
    double                                                   *var,
    int                                                       i);

  // use basis index and term index for individual basis function
  void
  PSeriesFunction_dsis_get_basis(
    PRISMS::PSeriesFunction<double, double, double *, int *> *f,
    int                                                       bindex,
    int                                                       term,
    double                                                   &val);

  void
  PSeriesFunction_dsis_get_basis_grad(
    PRISMS::PSeriesFunction<double, double, double *, int *> *f,
    int                                                       bindex,
    int                                                       term,
    double                                                   &val);

  void
  PSeriesFunction_dsis_get_basis_hess(
    PRISMS::PSeriesFunction<double, double, double *, int *> *f,
    int                                                       bindex,
    int                                                       term,
    double                                                   &val);

  // or use tensor indices to evaluate basis function product
  void
  PSeriesFunction_dsis_get_tensor_basis(
    PRISMS::PSeriesFunction<double, double, double *, int *> *f,
    int                                                      *term,
    double                                                   &val);

  void
  PSeriesFunction_dsis_get_tensor_basis_grad(
    PRISMS::PSeriesFunction<double, double, double *, int *> *f,
    int                                                      *term,
    int                                                       di,
    double                                                   &val);

  void
  PSeriesFunction_dsis_get_tensor_basis_hess(
    PRISMS::PSeriesFunction<double, double, double *, int *> *f,
    int                                                      *term,
    int                                                       di,
    int                                                       dj,
    double                                                   &val);

  // Functions for using constructing a 2D PRISMS::Body externally (say Python or
  // Fortran),
  //   allowing access to PFields
  //   written for Coordinate=double*, DIM=2

  void
  Body2D_new(char *vtkfile, PRISMS::Body<double *, 2> *&b);

  void
  Body2D_delete(PRISMS::Body<double *, 2> *&b);

  // Functions for using a 2D scalar PField externally (say Python or Fortran), as a
  // PFunction.
  //   From a Body pointer, returns a pointer to a PFuncBase
  //   written for Coordinate=double*, OutType=double, DIM=2

  void
  ScalarField2D(char                                 *name,
                PRISMS::Body<double *, 2>            *b,
                PRISMS::PFuncBase<double *, double> *&f);

  // Functions for using constructing a 3D PRISMS::Body externally (say Python or
  // Fortran),
  //   allowing access to PFields
  //   written for Coordinate=double*, DIM=3

  void
  Body3D_new(char *vtkfile, PRISMS::Body<double *, 3> *&b);

  void
  Body3D_delete(PRISMS::Body<double *, 3> *&b);

  // Functions for using a 3D scalar PField externally (say Python or Fortran), as a
  // PFunction.
  //   From a Body pointer, returns a pointer to a PFuncBase
  //   written for Coordinate=double*, OutType=double, DIM=3

  void
  ScalarField3D(char                                 *name,
                PRISMS::Body<double *, 3>            *b,
                PRISMS::PFuncBase<double *, double> *&f);
}

#endif