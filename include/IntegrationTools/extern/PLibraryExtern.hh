#ifndef PLIBRARY_HH
#define PLIBRARY_HH

#include "../PFunction.hh"
#include <cstring>
#include <vector>

namespace PRISMS
{

  /// Library where you can find functions and basis sets
  ///
  namespace PLibrary
  {
    // Use these functions to checkout objects which manage their own memory

    void
    checkout(std::string name, PSimpleFunction<double, double> &simplefunc);
    void
    checkout(std::string name, PSimpleFunction<std::vector<double>, double> &simplefunc);
    void
    checkout(std::string name, PSimpleFunction<double *, double> &simplefunc);

    void
    checkout(std::string name, PFunction<std::vector<double>, double> &func);
    void
    checkout(std::string name, PFunction<double *, double> &func);

    void
    checkout(std::string name, PBasisSet<double, double> &basis_set, int N);

    // Use these functions to checkout new 'Base' objects which the user must delete

    void
    checkout(std::string name, PSimpleBase<double, double> *&simplefunc);
    void
    checkout(std::string name, PSimpleBase<std::vector<double>, double> *&simplefunc);
    void
    checkout(std::string name, PSimpleBase<double *, double> *&simplefunc);

    void
    checkout(std::string name, PFuncBase<std::vector<double>, double> *&func);
    void
    checkout(std::string name, PFuncBase<double *, double> *&func);

    void
    checkout(std::string name, PBasisSetBase<double, double> *&basis_set, int N);
  } // namespace PLibrary

} // namespace PRISMS

#endif
