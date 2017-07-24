
#ifndef SimplePiece_HH
#define SimplePiece_HH

#include<iostream>
#include<vector>

#include "../pfunction/PSimpleFunction.hh"

namespace PRISMS
{   
    /// Class to contain a SimpleFunction and the piece in which it is valid.
    /// 
    ///   This can be evaluated in or out of the piece in which it is declared valid
    ///
    template< class VarContainer, class OutType>
    class SimplePiece : public PSimpleBase<VarContainer, OutType> 
    {
        protected:
        
        mutable PSimpleFunction<VarContainer, OutType> _expr;
        mutable std::vector< PSimpleFunction<VarContainer, bool> > _condition;
        typedef typename std::vector< PSimpleFunction<VarContainer, bool> >::size_type size_type;
        
        public:
        
        SimplePiece( const PSimpleFunction<VarContainer, OutType> &expr, const std::vector< PSimpleFunction<VarContainer, bool> > &condition)
        {
            _expr = expr;
            _condition = condition;
            this->_name = _expr.name();
        }
        
        virtual std::string csrc() const
        {
            std::string str = _expr.csrc();
            for( size_type i=0; i<_condition.size(); i++)
            {
                if( i == 0)
                {
                    str += " if " + _condition[i].csrc();
                }
                else if( i == _condition.size()-1 )
                {
                    str += " and " + _condition[i].csrc();
                }
                else
                {
                    str += ", " + _condition[i].csrc();
                }
            }
            return str;
        }
        
        virtual std::string sym() const
        {
            std::string str = _expr.sym();
            for( size_type i=0; i<_condition.size(); i++)
            {
                if( i == 0)
                {
                    str += " if " + _condition[i].sym();
                }
                else if( i == _condition.size()-1 )
                {
                    str += " and " + _condition[i].sym();
                }
                else
                {
                    str += ", " + _condition[i].sym();
                }
            }
            return str;
        }
        
        virtual std::string latex() const
        {
            std::string str = _expr.latex();
            for( size_type i=0; i<_condition.size(); i++)
            {
                if( i == 0)
                {
                    str += " & \\mbox{ if } " + _condition[i].latex();
                }
                else if( i == _condition.size()-1 )
                {
                    str += " \\mbox{ and } " + _condition[i].latex();
                }
                else
                {
                    str += " \\mbox{, } " + _condition[i].sym();
                }
            }
            return str;
        }
        
        virtual SimplePiece<VarContainer, OutType>* clone() const
        {
            return new SimplePiece<VarContainer, OutType>(*this);
        }
        
        bool in_piece( const VarContainer &var) const
        {
            for( size_type i=0; i<_condition.size(); i++)
            {
                if( !_condition[i](var) )
                    return false;
            }
            return true;
        }
        
        PSimpleFunction<VarContainer, OutType> expr() const
        {
            return _expr;
        }
        
        std::vector< PSimpleFunction<VarContainer, bool> > condition() const
        {
            return _condition;
        }
        
        private:
        
        /// This will return '_expr' evaluated anywhere. Must check in_piece first. We
        ///    don't check it here to avoid double checking when evaluating PPieceWiseSimpleBase
        ///
        virtual OutType eval( const VarContainer &var) const
        { 
            return _expr(var);
        }
    };

}


#endif