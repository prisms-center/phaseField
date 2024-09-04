
#ifndef PPieceWiseSimpleBase_HH
#define PPieceWiseSimpleBase_HH

#include "../pfunction/PSimpleBase.hh"
#include "./SimplePiece.hh"
#include <iostream>
#include <stdexcept>
#include <vector>

namespace PRISMS
{
  /// Class to define a PieceWise SimpleFunction
  ///
  ///   Contains a vector of 'SimplePiece'. Throws a domain_error if it
  ///   is evaluated outside of the valid domain of any piece.
  ///
  template <class VarContainer, class OutType>
  class PPieceWiseSimpleBase : public PSimpleBase<VarContainer, OutType>
  {
  public:
    mutable std::vector<SimplePiece<VarContainer, OutType>> _piece;

    PPieceWiseSimpleBase()
    {}

    PPieceWiseSimpleBase(const std::vector<SimplePiece<VarContainer, OutType>> &piece)
    {
      _piece = piece;
    }

    virtual std::string
    csrc() const
    {
      std::string str = "";
      for (int i = 0; i < _piece.size(); i++)
        {
          if (i == 0)
            {
              str += _piece[i].csrc();
            }
          else if (i == _piece.size() - 1)
            {
              str += "; and " + _piece[i].csrc();
            }
          else
            {
              str += "; " + _piece[i].csrc();
            }
        }
      return str;
    }

    virtual std::string
    sym() const
    {
      std::string str = "";
      for (int i = 0; i < _piece.size(); i++)
        {
          if (i == 0)
            {
              str += _piece[i].sym();
            }
          else if (i == _piece.size() - 1)
            {
              str += "; and " + _piece[i].sym();
            }
          else
            {
              str += "; " + _piece[i].sym();
            }
        }
      return str;
    }

    virtual std::string
    latex() const
    {
      std::string str = "";
      for (int i = 0; i < _piece.size(); i++)
        {
          if (i == 0)
            {
              str += "\\left\\{ \\begin{array}{ll} " + _piece[i].latex();
            }
          else if (i == _piece.size() - 1)
            {
              str += " \\\\ " + _piece[i].latex();
            }
          else
            {
              str += " \\\\ " + _piece[i].latex();
            }
        }

      str += " \\end{array} \\right.";
      return str;
    }

    virtual PPieceWiseSimpleBase<VarContainer, OutType> *
    clone() const
    {
      return new PPieceWiseSimpleBase<VarContainer, OutType>(*this);
    }

    bool
    in_piece(const VarContainer &var) const
    {
      for (int i = 0; i < _piece.size(); i++)
        {
          if (_piece[i].in_piece(var))
            return true;
        }
      return false;
    }

    int
    piece(const VarContainer &var)
    {
      for (int i = 0; i < _piece.size(); i++)
        {
          if (_piece[i].in_piece(var))
            return i;
        }

      throw std::domain_error("PPieceWiseSimpleBase: Not in any piece");
    }

  private:
    virtual OutType
    eval(const VarContainer &var) const
    {
      for (int i = 0; i < _piece.size(); i++)
        {
          if (_piece[i].in_piece(var))
            return _piece[i](var);
        }

      throw std::domain_error("PPieceWiseSimpleBase: Not in any piece");
    }
  };
} // namespace PRISMS

#endif