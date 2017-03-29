
#ifndef PSimpleFunction_HH
#define PSimpleFunction_HH

#include<cstring>
#include<iostream>
#include<vector>
#include<cstdlib>

#include "./PSimpleBase.hh"

namespace PRISMS
{
    template< class VarContainer, class OutType>
    class PSimpleFunction 
    {
    private:
        PSimpleBase<VarContainer,OutType> *ptr;

    public:

        std::string name() const
        {
            return (*ptr).name();
        }
        std::string csrc() const
        {
            return (*ptr).csrc();
        }
        std::string sym() const
        {
            return (*ptr).sym();
        }
        std::string latex() const
        {
            return (*ptr).latex();
        }
        

        // ----------------------------------------------------------
        // Use this function if you want to evaluate,
        //   return and store result
        OutType operator()(const VarContainer &var)
        {
            return (*ptr)(var);
        }

        // ----------------------------------------------------------
        // Then use 'get' methods to access results later
        void eval(const VarContainer &var)
        {
            (*ptr)(var);
        }

        OutType operator()() const
        {
            return (*ptr)();
        }

        // PFunction unique members ------------------------------------------

        PSimpleFunction& operator=(const PSimpleFunction &RHS)
        {
            if(ptr != NULL)
                delete ptr;
            ptr = RHS.ptr->clone();
            return *this;
        }

        template<class T> 
        PSimpleFunction& operator=(const T &RHS)
        {
            RHS.is_derived_from_PSimpleBase();

            if(ptr != NULL)
                delete ptr;
            ptr = RHS.clone();
            return *this;
        }

        // If you use this, PSimpleFunction becomes the 'owner' of the function RHS points to
        //    and it will delete it
        PSimpleFunction& set( PSimpleBase<VarContainer,OutType> *RHS)
        {
            if(RHS == NULL)
            {
                std::cout << "Error in PSimpleFunction::set. RHS == NULL." << std::endl;
                exit(1);
            }
            if(ptr != NULL)
                delete ptr;
            ptr = RHS;
            return *this;
        }
        
        PSimpleFunction()
        {
            ptr = NULL;
        }

        PSimpleFunction(const PSimpleFunction &RHS)
        {
            if( RHS.ptr != NULL)
                ptr = RHS.ptr->clone();
            else 
                ptr = NULL;
        }
        
        template<class T> PSimpleFunction(const T &RHS)
        {
            RHS.is_derived_from_PSimpleBase();

            ptr = RHS.clone();
            
        }

        ~PSimpleFunction()
        {
            if(ptr != NULL)
                delete ptr;
        }

    
    };

}


#endif