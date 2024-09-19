// created: 2017-1-11 12:28:14
// version: issue-1
// url: https://github.com/prisms-center/IntegrationTools
// commit: 34f71553d3b1b3cdbbc6a02cfb8be793a16396b5

#ifndef PLIBRARY_HH
#define PLIBRARY_HH

#include "../../../include/IntegrationTools/PFunction.hh"
#include "../../../include/IntegrationTools/PPieceWise.hh"
#include <cstring>

namespace PRISMS
{

  /// Library where you can find functions and basis sets
  ///
  namespace PLibrary
  {
    // Use these functions to checkout objects which manage their own memory

    void
    checkout(std::string name, PSimpleFunction<double *, double> &simplefunc);

    void
    checkout(std::string name, PFunction<double *, double> &func);

    // Use these functions to checkout new 'Base' objects which the user must delete

    void
    checkout(std::string name, PSimpleBase<double *, double> *&simplefunc);

    void
    checkout(std::string name, PFuncBase<double *, double> *&func);

  } // namespace PLibrary

} // namespace PRISMS

#endif
