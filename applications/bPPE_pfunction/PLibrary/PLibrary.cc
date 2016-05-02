// created: 2016-4-27 10:59:25
// version: issue-1
// url: https://github.com/prisms-center/IntegrationTools
// commit: 34f71553d3b1b3cdbbc6a02cfb8be793a16396b5

#ifndef PLIBRARY_CC
#define PLIBRARY_CC

#include<cstring>
#include<stdexcept>
#include "pfunct_faV.hh"
#include "pfunct_McV.hh"
#include "PLibrary.hh"

namespace PRISMS
{

        void PLibrary::checkout( std::string name, PSimpleFunction< dealii::VectorizedArray<double>, double > &simplefunc)
        {
            typedef PSimpleFunction< dealii::VectorizedArray<double>, double > psf;
            if( name == "pfunct_faV_f") { simplefunc = psf( pfunct_faV_f< dealii::VectorizedArray<double> >() ); return;}
            if( name == "pfunct_faV_grad_0") { simplefunc = psf( pfunct_faV_grad_0< dealii::VectorizedArray<double> >() ); return;}
            if( name == "pfunct_faV_hess_0_0") { simplefunc = psf( pfunct_faV_hess_0_0< dealii::VectorizedArray<double> >() ); return;}
            if( name == "pfunct_McV_f") { simplefunc = psf( pfunct_McV_f< dealii::VectorizedArray<double> >() ); return;}
            else throw std::runtime_error( "PSimpleFunction< dealii::VectorizedArray<double>, double > " + name + " was not found in the PLibrary");
        }
        void PLibrary::checkout( std::string name, PSimpleFunction< double*, double > &simplefunc)
        {
            typedef PSimpleFunction< double*, double > psf;
            if( name == "pfunct_faV_f") { simplefunc = psf( pfunct_faV_f< double* >() ); return;}
            if( name == "pfunct_faV_grad_0") { simplefunc = psf( pfunct_faV_grad_0< double* >() ); return;}
            if( name == "pfunct_faV_hess_0_0") { simplefunc = psf( pfunct_faV_hess_0_0< double* >() ); return;}
            if( name == "pfunct_McV_f") { simplefunc = psf( pfunct_McV_f< double* >() ); return;}
            else throw std::runtime_error( "PSimpleFunction< double*, double > " + name + " was not found in the PLibrary");
        }

        void PLibrary::checkout( std::string name, PFunction< dealii::VectorizedArray<double>, double > &func)
        {
            typedef PFunction< dealii::VectorizedArray<double>, double > pf;
            if( name == "pfunct_faV") { func = pf( pfunct_faV< dealii::VectorizedArray<double> >() ); return;}
            if( name == "pfunct_McV") { func = pf( pfunct_McV< dealii::VectorizedArray<double> >() ); return;}
            throw std::runtime_error( "PFunction< dealii::VectorizedArray<double>, double > " + name + " was not found in the PLibrary");
        }
        void PLibrary::checkout( std::string name, PFunction< double*, double > &func)
        {
            typedef PFunction< double*, double > pf;
            if( name == "pfunct_faV") { func = pf( pfunct_faV< double* >() ); return;}
            if( name == "pfunct_McV") { func = pf( pfunct_McV< double* >() ); return;}
            throw std::runtime_error( "PFunction< double*, double > " + name + " was not found in the PLibrary");
        }




        void PLibrary::checkout( std::string name, PSimpleBase< dealii::VectorizedArray<double>, double > *&simplefunc)
        {
            if( name == "pfunct_faV_f") { simplefunc = new pfunct_faV_f< dealii::VectorizedArray<double> >(); return;}
            if( name == "pfunct_faV_grad_0") { simplefunc = new pfunct_faV_grad_0< dealii::VectorizedArray<double> >(); return;}
            if( name == "pfunct_faV_hess_0_0") { simplefunc = new pfunct_faV_hess_0_0< dealii::VectorizedArray<double> >(); return;}
            if( name == "pfunct_McV_f") { simplefunc = new pfunct_McV_f< dealii::VectorizedArray<double> >(); return;}
            throw std::runtime_error( "PSimpleBase< dealii::VectorizedArray<double>, double > " + name + " was not found in the PLibrary");
        }
        void PLibrary::checkout( std::string name, PSimpleBase< double*, double > *&simplefunc)
        {
            if( name == "pfunct_faV_f") { simplefunc = new pfunct_faV_f< double* >(); return;}
            if( name == "pfunct_faV_grad_0") { simplefunc = new pfunct_faV_grad_0< double* >(); return;}
            if( name == "pfunct_faV_hess_0_0") { simplefunc = new pfunct_faV_hess_0_0< double* >(); return;}
            if( name == "pfunct_McV_f") { simplefunc = new pfunct_McV_f< double* >(); return;}
            throw std::runtime_error( "PSimpleBase< double*, double > " + name + " was not found in the PLibrary");
        }

        void PLibrary::checkout( std::string name, PFuncBase< dealii::VectorizedArray<double>, double > *&func)
        {
            if( name == "pfunct_faV") { func = new pfunct_faV< dealii::VectorizedArray<double> >(); return;}
            if( name == "pfunct_McV") { func = new pfunct_McV< dealii::VectorizedArray<double> >(); return;}
            throw std::runtime_error( "PFuncBase< dealii::VectorizedArray<double>, double > " + name + " was not found in the PLibrary");
        }
        void PLibrary::checkout( std::string name, PFuncBase< double*, double > *&func)
        {
            if( name == "pfunct_faV") { func = new pfunct_faV< double* >(); return;}
            if( name == "pfunct_McV") { func = new pfunct_McV< double* >(); return;}
            throw std::runtime_error( "PFuncBase< double*, double > " + name + " was not found in the PLibrary");
        }




}


#endif
