
#ifndef Piece_HH
#define Piece_HH

#include "../pfunction/PFunction.hh"
#include "./SimplePiece.hh"
#include <iostream>
#include <vector>

namespace PRISMS
{

  /// Class to contain a Function and the piece in which it is valid.
  ///
  ///   This can be evaluated in or out of the piece in which it is declared valid
  ///
  template <class VarContainer, class OutType>
  class Piece : public PFuncBase<VarContainer, OutType>
  {
  protected:
    PFunction<VarContainer, OutType>                         _expr;
    mutable std::vector<PSimpleFunction<VarContainer, bool>> _condition;

    typedef
      typename std::vector<PSimpleFunction<VarContainer, bool>>::size_type cond_size_type;

  public:
    typedef typename PFuncBase<VarContainer, OutType>::size_type size_type;

    Piece(const PFunction<VarContainer, OutType>                 &expr,
          const std::vector<PSimpleFunction<VarContainer, bool>> &condition)
    {
      _expr                  = expr;
      _condition             = condition;
      this->_name            = _expr.name();
      this->_var_name        = _expr.var_name();
      this->_var_description = _expr.var_description();
    }

    bool
    in_piece(const VarContainer &var) const
    {
      for (cond_size_type i = 0; i < _condition.size(); i++)
        {
          if (!_condition[i](var))
            return false;
        }
      return true;
    }

    PFunction<VarContainer, OutType>
    expr() const
    {
      return _expr;
    }

    std::vector<PSimpleFunction<VarContainer, bool>>
    condition() const
    {
      return _condition;
    }

    SimplePiece<VarContainer, OutType>
    simplepiece() const
    {
      return SimplePiece<VarContainer, OutType>(_expr.simplefunction(), _condition);
    }

    SimplePiece<VarContainer, OutType>
    grad_simplepiece(size_type di) const
    {
      return SimplePiece<VarContainer, OutType>(_expr.grad_simplefunction(di),
                                                _condition);
    }

    SimplePiece<VarContainer, OutType>
    hess_simplepiece(size_type di, size_type dj) const
    {
      return SimplePiece<VarContainer, OutType>(_expr.hess_simplefunction(di, dj),
                                                _condition);
    }

    virtual Piece<VarContainer, OutType> *
    clone() const
    {
      return new Piece<VarContainer, OutType>(*this);
    }

    virtual PSimpleFunction<VarContainer, OutType>
    simplefunction() const
    {
      return PSimpleFunction<VarContainer, OutType>(
        SimplePiece<VarContainer, OutType>(_expr.simplefunction(), _condition));
    }

    virtual PSimpleFunction<VarContainer, OutType>
    grad_simplefunction(size_type di) const
    {
      return PSimpleFunction<VarContainer, OutType>(
        SimplePiece<VarContainer, OutType>(_expr.grad_simplefunction(di), _condition));
    }

    virtual PSimpleFunction<VarContainer, OutType>
    hess_simplefunction(size_type di, size_type dj) const
    {
      return PSimpleFunction<VarContainer, OutType>(
        SimplePiece<VarContainer, OutType>(_expr.hess_simplefunction(di, dj),
                                           _condition));
    }

    // ----------------------------------------------------------
    // Use these functions if you want to evaluate a single value

    /// These will return '_expr' evaluated anywhere. Must check in_piece first. We
    ///    don't check it here to avoid double checking when evaluating PPieceWiseFuncBase
    ///
    virtual OutType
    operator()(const VarContainer &var)
    {
      return _expr(var);
    }

    virtual OutType
    grad(const VarContainer &var, size_type di)
    {
      return _expr.grad(var, di);
    }

    virtual OutType
    hess(const VarContainer &var, size_type di, size_type dj)
    {
      return _expr.hess(var, di, dj);
    }

    // ----------------------------------------------------------
    // Use these functions to evaluate several values, then use 'get' methods to access
    // results
    virtual void
    eval(const VarContainer &var)
    {
      _expr(var);
    }

    virtual void
    eval_grad(const VarContainer &var)
    {
      _expr.eval_grad(var);
    }

    virtual void
    eval_hess(const VarContainer &var)
    {
      _expr.eval_hess(var);
    }

    /// These don't recheck the domain
    virtual OutType
    operator()() const
    {
      return _expr();
    }

    virtual OutType
    grad(size_type di) const
    {
      return _expr.grad(di);
    }

    virtual OutType
    hess(size_type di, size_type dj) const
    {
      return _expr.hess(di, dj);
    }
  };

} // namespace PRISMS

#endif