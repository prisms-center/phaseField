// created: 2016-4-20 11:47:45
// version: master
// url: https://github.com/prisms-center/IntegrationTools
// commit: 15edfe1511ccbdb131d055b89133fd543c7c0607

#ifndef PLIBRARY_CC
#define PLIBRARY_CC

#include<cstring>
#include<stdexcept>
#include<vector>
#include "pfunct_faV.hh"
#include "PLibrary.hh"

namespace PRISMS
{

        void PLibrary::checkout( std::string name, PSimpleFunction< std::vector<double>, double > &simplefunc)
        {
            typedef PSimpleFunction< std::vector<double>, double > psf;
            if( name == "pfunct_faV_f") { simplefunc = psf( pfunct_faV_f< std::vector<double> >() ); return;}
            if( name == "pfunct_faV_grad_0") { simplefunc = psf( pfunct_faV_grad_0< std::vector<double> >() ); return;}
            if( name == "pfunct_faV_hess_0_0") { simplefunc = psf( pfunct_faV_hess_0_0< std::vector<double> >() ); return;}
            else throw std::runtime_error( "PSimpleFunction< std::vector<double>, double > " + name + " was not found in the PLibrary");
        }
        void PLibrary::checkout( std::string name, PSimpleFunction< double*, double > &simplefunc)
        {
            typedef PSimpleFunction< double*, double > psf;
            if( name == "pfunct_faV_f") { simplefunc = psf( pfunct_faV_f< double* >() ); return;}
            if( name == "pfunct_faV_grad_0") { simplefunc = psf( pfunct_faV_grad_0< double* >() ); return;}
            if( name == "pfunct_faV_hess_0_0") { simplefunc = psf( pfunct_faV_hess_0_0< double* >() ); return;}
            else throw std::runtime_error( "PSimpleFunction< double*, double > " + name + " was not found in the PLibrary");
        }

        void PLibrary::checkout( std::string name, PFunction< std::vector<double>, double > &func)
        {
            typedef PFunction< std::vector<double>, double > pf;
            if( name == "pfunct_faV") { func = pf( pfunct_faV< std::vector<double> >() ); return;}
            throw std::runtime_error( "PFunction< std::vector<double>, double > " + name + " was not found in the PLibrary");
        }
        void PLibrary::checkout( std::string name, PFunction< double*, double > &func)
        {
            typedef PFunction< double*, double > pf;
            if( name == "pfunct_faV") { func = pf( pfunct_faV< double* >() ); return;}
            throw std::runtime_error( "PFunction< double*, double > " + name + " was not found in the PLibrary");
        }




        void PLibrary::checkout( std::string name, PSimpleBase< std::vector<double>, double > *&simplefunc)
        {
            if( name == "pfunct_faV_f") { simplefunc = new pfunct_faV_f< std::vector<double> >(); return;}
            if( name == "pfunct_faV_grad_0") { simplefunc = new pfunct_faV_grad_0< std::vector<double> >(); return;}
            if( name == "pfunct_faV_hess_0_0") { simplefunc = new pfunct_faV_hess_0_0< std::vector<double> >(); return;}
            throw std::runtime_error( "PSimpleBase< std::vector<double>, double > " + name + " was not found in the PLibrary");
        }
        void PLibrary::checkout( std::string name, PSimpleBase< double*, double > *&simplefunc)
        {
            if( name == "pfunct_faV_f") { simplefunc = new pfunct_faV_f< double* >(); return;}
            if( name == "pfunct_faV_grad_0") { simplefunc = new pfunct_faV_grad_0< double* >(); return;}
            if( name == "pfunct_faV_hess_0_0") { simplefunc = new pfunct_faV_hess_0_0< double* >(); return;}
            throw std::runtime_error( "PSimpleBase< double*, double > " + name + " was not found in the PLibrary");
        }

        void PLibrary::checkout( std::string name, PFuncBase< std::vector<double>, double > *&func)
        {
            if( name == "pfunct_faV") { func = new pfunct_faV< std::vector<double> >(); return;}
            throw std::runtime_error( "PFuncBase< std::vector<double>, double > " + name + " was not found in the PLibrary");
        }
        void PLibrary::checkout( std::string name, PFuncBase< double*, double > *&func)
        {
            if( name == "pfunct_faV") { func = new pfunct_faV< double* >(); return;}
            throw std::runtime_error( "PFuncBase< double*, double > " + name + " was not found in the PLibrary");
        }




}


#endif
