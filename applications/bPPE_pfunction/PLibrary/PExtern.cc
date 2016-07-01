#ifndef PExtern_CC
#define PExtern_CC

#include<cstring>
#include<iostream>
#include<vector>
#include<cstdlib>

#include "IntegrationTools/PExtern.hh"
#include "IntegrationTools/PFunction.hh"
#include "IntegrationTools/PField.hh"

// In future, might have more complicated OutType, 
//   so make all have 'void' return and pass everything by reference


extern "C"
{
    // Functions for using a PSimpleBase externally (say Python or Fortran)
    //   written for VarContainer=double*, OutType=double, hence 'dsd' in function names
    
    void PSimpleFunction_dsd_new(char* name, PRISMS::PSimpleBase<double*,double>* &f)
    { 
        PRISMS::PLibrary::checkout(std::string(name), f);
    }
    
    void PSimpleFunction_dsd_delete(PRISMS::PSimpleBase<double*,double>* &f)
    { 
        delete f;
        f = NULL;
    }
    
    void PSimpleFunction_dsd_name(PRISMS::PSimpleBase<double*,double> *f, char* name)
    {
        std::strcpy(name, f->name().c_str());
    }
    
    void PSimpleFunction_dsd_calc( PRISMS::PSimpleBase<double*,double> *f, double* var, double &val)
    {
        val = (*f)(var);
    }
    
    void PSimpleFunction_dsd_get( PRISMS::PSimpleBase<double*,double> *f, double  &val)
    {
        val = (*f)();
    }
    
    
    
    
    
    // Functions for using a PFuncBase externally (say Python or Fortran)
    //   written for VarContainer=double*, OutType=double, hence 'dsd' in function names
    
    void PFunction_dsd_new(char* name, PRISMS::PFuncBase<double*,double>* &f)
    { 
        PRISMS::PLibrary::checkout(std::string(name), f);
    }
    
    void PFunction_dsd_delete(PRISMS::PFuncBase<double*,double>* &f)
    { 
        delete f;
        f = NULL;
    }
    
    void PFunction_dsd_name(PRISMS::PFuncBase<double*,double>* f, char* name)
    {
        std::strcpy(name, f->name().c_str());
    }
    
    void PFunction_dsd_size(PRISMS::PFuncBase<double*,double>* f, int &size)
    {
        size = f->size();
    }
    
    void PFunction_dsd_var_name(PRISMS::PFuncBase<double*,double>* f, int i, char* var_name)
    {
        std::strcpy(var_name, f->var_name(i).c_str());
    }
    
    void PFunction_dsd_var_description(PRISMS::PFuncBase<double*,double>* f, int i, char* var_description)
    {
        std::strcpy(var_description, f->var_description(i).c_str());
    }
    
    //void PFunction_dsd_simplefunc(PRISMS::PFuncBase<double*,double> *f, PSimpleBase<double*, double> *simplefunc);
    //void PFunction_dsd_grad_simplefunc(PRISMS::PFuncBase<double*,double> *f, int *di, PSimpleBase<double*, double> *simplefunc);
    //void PFunction_dsd_hess_simplefunc(PRISMS::PFuncBase<double*,double> *f, int *di, int *dj, PSimpleBase<double*, double> *simplefunc);
    
    void PFunction_dsd_calc(PRISMS::PFuncBase<double*,double>* f, double* var, double &val)
    {
        val = (*f)(var);
    }
    
    void PFunction_dsd_calc_grad(PRISMS::PFuncBase<double*,double>* f, double* var, int di, double &val)
    {
        val = (*f).grad(var, di);
    }
    
    void PFunction_dsd_calc_hess(PRISMS::PFuncBase<double*,double>* f, double* var, int di, int dj, double &val)
    {
        val = (*f).hess(var, di, dj);
    }
    
    void PFunction_dsd_eval(PRISMS::PFuncBase<double*,double>* f, double* var)
    {
        (*f)(var);
    }
    
    void PFunction_dsd_eval_grad(PRISMS::PFuncBase<double*,double>* f, double* var, int di)
    {
        (*f).grad(var, di);
    }
    
    void PFunction_dsd_eval_hess(PRISMS::PFuncBase<double*,double>* f, double* var, int di, int dj)
    {
        (*f).hess(var, di, dj);
    }
    
    void PFunction_dsd_get(PRISMS::PFuncBase<double*,double>* f, double &val)
    {
        val = (*f)();
    }
    
    void PFunction_dsd_get_grad(PRISMS::PFuncBase<double*,double>* f, int di, double &val)
    {
        val = (*f).grad(di);
    }
    
    void PFunction_dsd_get_hess(PRISMS::PFuncBase<double*,double>* f, int di, int dj, double &val)
    {
        val = (*f).hess(di, dj);
    }
    
    
    
    // Functions for using constructing a 2D PRISMS::Body externally (say Python or Fortran),
    //   allowing access to PFields
    //   written for Coordinate=double*, OutType=double, DIM=2
    
    void Body2D_new(char* vtkfile, PRISMS::Body<double*,2>* &b)
    {
        b = new PRISMS::Body<double*,2>();
        (*b).read_vtk(std::string(vtkfile));
    };
    
    void Body2D_delete(PRISMS::Body<double*,2>* &b)
    {
        delete b;
        b = NULL;
    };
    
    
    // Functions for using a 2D scalar PField externally (say Python or Fortran), as a PFunction.
    //   From a Body pointer, returns a pointer to a PFuncBase
    //   written for Coordinate=double*, OutType=double, DIM=2
    //   don't delete this! it will be deleted by deleting the Body
    
    void ScalarField2D(char* name, PRISMS::Body<double*,2>* b, PRISMS::PFuncBase<double*,double>* &f)
    {
        f = &((*b).find_scalar_field(std::string(name)));
    };
    
    
    // Functions for using constructing a 3D PRISMS::Body externally (say Python or Fortran),
    //   allowing access to PFields
    //   written for Coordinate=double*, OutType=double, DIM=3
    
    void Body3D_new(char* vtkfile, PRISMS::Body<double*,3>* &b)
    {
        b = new PRISMS::Body<double*,3>();
        (*b).read_vtk(std::string(vtkfile));
    };
    
    void Body3D_delete(PRISMS::Body<double*,3>* &b)
    {
        delete b;
        b = NULL;
    };
    
    
    // Functions for using a 3D scalar PField externally (say Python or Fortran), as a PFunction.
    //   From a Body pointer, returns a pointer to a PFuncBase
    //   written for Coordinate=double*, OutType=double, DIM=2
    //   don't delete this! it will be deleted by deleting the Body
    
    void ScalarField3D(char* name, PRISMS::Body<double*,3>* b, PRISMS::PFuncBase<double*,double>* &f)
    {
        f = &((*b).find_scalar_field(std::string(name)));
    };
    
    
}


#endif
