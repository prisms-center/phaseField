
#ifndef PFunction_HH
#define PFunction_HH

#include<cstring>
#include<iostream>
#include<vector>
#include<cstdlib>

#include "./PSimpleFunction.hh"
#include "./PFuncBase.hh"

namespace PRISMS
{    
    
    /// A class that contains a ptr to a PFuncBase object
    ///   - like a smart ptr class
    ///   - same interface as PFuncBase
    ///   - allows for using PFuncBase objects polymorphically 
    ///     without dereferencing and without worrying about new/delete
    ///   - VarContainer is either a scalar or something whose elements can be accessed with operator[]
    ///    
    ///   example: MyFuncA, MyFuncB, MyFuncC, etc. are defined:
    ///     template< class VarContainer>
    ///     MyFuncX : public PFuncBase<VarContainer, double>
    ///
    ///   // Then you can do things like this:
    ///   
    ///   MyFuncA<std::vector<double> > my_func_a;
    ///   MyFuncB<std::vector<double> > my_func_b;
    ///
    ///   PFuncBase<std::vector<double>, double>* my_func_c_ptr;
    ///   my_func_c_ptr = new MyFuncC<std::vector<double>, double >();
    ///
    ///   PFunction<std::vector<double>, double > f, g, h;
    ///
    ///   f = my_func_a;
    ///   f = my_func_b;
    ///   g.set(my_func_c_ptr->clone());
    ///   h.set(my_func_c_ptr);
    ///   double result = f(3.0) + g(4.0) + h(5.0);
    ///   
    ///   - No deletions are used in this example.  
    ///     PFunction::set makes PFunction the 'owner' of the MyFuncC object and it will delete it.
    ///
    template< class VarContainer, class OutType>
    class PFunction 
    {
    public:

        std::string name() const
        {
            return (*ptr).name();
        }
        int size() const
        {
            return (*ptr).size();
        }
        std::vector<std::string> var_name()
        {
            return (*ptr).var_name();
        }
        std::string var_name(int i)
        {
            return (*ptr).var_name(i);
        }
        std::vector<std::string> var_description()
        {
            return (*ptr).var_description();
        }
        std::string var_description(int i)
        {
            return (*ptr).var_description(i);
        }
        
        PSimpleFunction<VarContainer, OutType> simplefunction() const
        {
            return (*ptr).simplefunction();
        }
        
        PSimpleFunction<VarContainer, OutType> grad_simplefunction(int di) const
        {
            return (*ptr).grad_simplefunction(di);
        }
        
        PSimpleFunction<VarContainer, OutType> hess_simplefunction(int di, int dj) const
        {
            return (*ptr).hess_simplefunction(di, dj);
        }

        // ----------------------------------------------------------
        // Use these functions if you want to evaluate a single value
        OutType operator()(const VarContainer &var)
        {
            return (*ptr)(var);
        }
        OutType grad(const VarContainer &var, int di)
        {
            return (*ptr).grad(var, di);
        }
        OutType hess(const VarContainer &var, int di, int dj)
        {
            return (*ptr).hess(var, di, dj);
        }

        // ----------------------------------------------------------
        // Use these functions to evaluate several values, then use 'get' methods to access results
        void eval(const VarContainer &var)
        {
            return (*ptr).eval(var);
        }
        void eval_grad(const VarContainer &var)
        {
            return (*ptr).eval_grad(var);
        }
        void eval_hess(const VarContainer &var)
        {
            return (*ptr).eval_hess(var);
        }

        OutType operator()() const
        {
            return (*ptr)();
        }
        OutType grad(int di) const
        {
            return (*ptr).grad(di);
        }
        OutType hess(int di, int dj) const
        {
            return (*ptr).hess(di, dj);
        }


        // PFunction unique members ------------------------------------------

        PFunction &operator=(const PFunction &RHS)
        {
            if(ptr != NULL)
                delete ptr;
            ptr = RHS.ptr->clone();
            return *this;
        }

        template<class T> PFunction &operator=(const T &RHS)
        {
            RHS.is_derived_from_PFuncBase();

            if(ptr != NULL)
                delete ptr;
            ptr = RHS.clone();
            return *this;
        }

        // If you use this, PFunction becomes the 'owner' of the function RHS points to
        //    and it will delete it
        PFunction &set( PFuncBase<VarContainer,OutType> *RHS)
        {
            if(RHS == NULL)
            {
                std::cout << "Error in PFunction::set. RHS == NULL." << std::endl;
                exit(1);
            }
            if(ptr != NULL)
                delete ptr;
            ptr = RHS;
            return *this;
        }

        PFunction()
        {
            ptr = NULL;
        }

        PFunction(const PFunction &RHS)
        {
            if( RHS.ptr != NULL)
                ptr = RHS.ptr->clone();
            else 
                ptr = NULL;
        }
        
        template<class T> PFunction(const T &RHS)
        {
            RHS.is_derived_from_PFuncBase();

            ptr = RHS.clone();
            
        }

        ~PFunction()
        {
            if(ptr != NULL)
                delete ptr;
        }

    private:
        PFuncBase<VarContainer,OutType> *ptr;

    };

}


#endif