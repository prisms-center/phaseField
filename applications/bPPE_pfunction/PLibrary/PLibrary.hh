// created: 2016-4-20 11:47:45
// version: master
// url: https://github.com/prisms-center/IntegrationTools
// commit: 15edfe1511ccbdb131d055b89133fd543c7c0607

#ifndef PLIBRARY_HH
#define PLIBRARY_HH

#include<cstring>
#include<vector>
#include "IntegrationTools/PFunction.hh"
#include "IntegrationTools/PPieceWise.hh"

namespace PRISMS
{

    /// Library where you can find functions and basis sets
    ///
    namespace PLibrary
    {
        // Use these functions to checkout objects which manage their own memory

        void checkout( std::string name, PSimpleFunction< std::vector<double>, double > &simplefunc);
        void checkout( std::string name, PSimpleFunction< double*, double > &simplefunc);

        void checkout( std::string name, PFunction< std::vector<double>, double > &func);
        void checkout( std::string name, PFunction< double*, double > &func);




        // Use these functions to checkout new 'Base' objects which the user must delete

        void checkout( std::string name, PSimpleBase< std::vector<double>, double > *&simplefunc);
        void checkout( std::string name, PSimpleBase< double*, double > *&simplefunc);

        void checkout( std::string name, PFuncBase< std::vector<double>, double > *&func);
        void checkout( std::string name, PFuncBase< double*, double > *&func);



    }

}


#endif
