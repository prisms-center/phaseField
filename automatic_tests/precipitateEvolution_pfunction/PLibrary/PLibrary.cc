// created: 2017-1-11 12:28:14
// version: issue-1
// url: https://github.com/prisms-center/IntegrationTools
// commit: 34f71553d3b1b3cdbbc6a02cfb8be793a16396b5

#ifndef PLIBRARY_CC
#define PLIBRARY_CC

#include "PLibrary.hh"

#include "pfunct_McV.hh"
#include "pfunct_Mn1V.hh"
#include "pfunct_Mn2V.hh"
#include "pfunct_Mn3V.hh"
#include "pfunct_faV.hh"
#include "pfunct_fbV.hh"

#include <cstring>
#include <stdexcept>

namespace PRISMS
{

  void
  PLibrary::checkout(std::string name, PSimpleFunction<double *, double> &simplefunc)
  {
    typedef PSimpleFunction<double *, double> psf;
    if (name == "pfunct_faV_f")
      {
        simplefunc = psf(pfunct_faV_f<double *>());
        return;
      }
    if (name == "pfunct_faV_grad_0")
      {
        simplefunc = psf(pfunct_faV_grad_0<double *>());
        return;
      }
    if (name == "pfunct_faV_hess_0_0")
      {
        simplefunc = psf(pfunct_faV_hess_0_0<double *>());
        return;
      }
    if (name == "pfunct_fbV_f")
      {
        simplefunc = psf(pfunct_fbV_f<double *>());
        return;
      }
    if (name == "pfunct_fbV_grad_0")
      {
        simplefunc = psf(pfunct_fbV_grad_0<double *>());
        return;
      }
    if (name == "pfunct_fbV_hess_0_0")
      {
        simplefunc = psf(pfunct_fbV_hess_0_0<double *>());
        return;
      }
    if (name == "pfunct_McV_f")
      {
        simplefunc = psf(pfunct_McV_f<double *>());
        return;
      }
    if (name == "pfunct_Mn1V_f")
      {
        simplefunc = psf(pfunct_Mn1V_f<double *>());
        return;
      }
    if (name == "pfunct_Mn2V_f")
      {
        simplefunc = psf(pfunct_Mn2V_f<double *>());
        return;
      }
    if (name == "pfunct_Mn3V_f")
      {
        simplefunc = psf(pfunct_Mn3V_f<double *>());
        return;
      }
    else
      throw std::runtime_error("PSimpleFunction< double*, double > " + name +
                               " was not found in the PLibrary");
  }

  void
  PLibrary::checkout(std::string name, PFunction<double *, double> &func)
  {
    typedef PFunction<double *, double> pf;
    if (name == "pfunct_faV")
      {
        func = pf(pfunct_faV<double *>());
        return;
      }
    if (name == "pfunct_fbV")
      {
        func = pf(pfunct_fbV<double *>());
        return;
      }
    if (name == "pfunct_McV")
      {
        func = pf(pfunct_McV<double *>());
        return;
      }
    if (name == "pfunct_Mn1V")
      {
        func = pf(pfunct_Mn1V<double *>());
        return;
      }
    if (name == "pfunct_Mn2V")
      {
        func = pf(pfunct_Mn2V<double *>());
        return;
      }
    if (name == "pfunct_Mn3V")
      {
        func = pf(pfunct_Mn3V<double *>());
        return;
      }
    throw std::runtime_error("PFunction< double*, double > " + name +
                             " was not found in the PLibrary");
  }

  void
  PLibrary::checkout(std::string name, PSimpleBase<double *, double> *&simplefunc)
  {
    if (name == "pfunct_faV_f")
      {
        simplefunc = new pfunct_faV_f<double *>();
        return;
      }
    if (name == "pfunct_faV_grad_0")
      {
        simplefunc = new pfunct_faV_grad_0<double *>();
        return;
      }
    if (name == "pfunct_faV_hess_0_0")
      {
        simplefunc = new pfunct_faV_hess_0_0<double *>();
        return;
      }
    if (name == "pfunct_fbV_f")
      {
        simplefunc = new pfunct_fbV_f<double *>();
        return;
      }
    if (name == "pfunct_fbV_grad_0")
      {
        simplefunc = new pfunct_fbV_grad_0<double *>();
        return;
      }
    if (name == "pfunct_fbV_hess_0_0")
      {
        simplefunc = new pfunct_fbV_hess_0_0<double *>();
        return;
      }
    if (name == "pfunct_McV_f")
      {
        simplefunc = new pfunct_McV_f<double *>();
        return;
      }
    if (name == "pfunct_Mn1V_f")
      {
        simplefunc = new pfunct_Mn1V_f<double *>();
        return;
      }
    if (name == "pfunct_Mn2V_f")
      {
        simplefunc = new pfunct_Mn2V_f<double *>();
        return;
      }
    if (name == "pfunct_Mn3V_f")
      {
        simplefunc = new pfunct_Mn3V_f<double *>();
        return;
      }
    throw std::runtime_error("PSimpleBase< double*, double > " + name +
                             " was not found in the PLibrary");
  }

  void
  PLibrary::checkout(std::string name, PFuncBase<double *, double> *&func)
  {
    if (name == "pfunct_faV")
      {
        func = new pfunct_faV<double *>();
        return;
      }
    if (name == "pfunct_fbV")
      {
        func = new pfunct_fbV<double *>();
        return;
      }
    if (name == "pfunct_McV")
      {
        func = new pfunct_McV<double *>();
        return;
      }
    if (name == "pfunct_Mn1V")
      {
        func = new pfunct_Mn1V<double *>();
        return;
      }
    if (name == "pfunct_Mn2V")
      {
        func = new pfunct_Mn2V<double *>();
        return;
      }
    if (name == "pfunct_Mn3V")
      {
        func = new pfunct_Mn3V<double *>();
        return;
      }
    throw std::runtime_error("PFuncBase< double*, double > " + name +
                             " was not found in the PLibrary");
  }

} // namespace PRISMS

#endif
