#ifndef Coordinate_HH
#define Coordinate_HH

#include <ostream>

namespace PRISMS
{

  /**
   * \brief A class for a coordinate, templated by dimension. This is a possible option
   * anyplace 'Coordinate' class template is used but it is not the only option. Any class
   * that implements 'Coordinate::operator[]()' should work.
   */
  template <int dim>
  class Coordinate
  {
  public:
    /**
     * \brief Get the number of dimensions.
     */
    [[nodiscard]] int
    size() const
    {
      return dim;
    }

    /**
     * \brief Get the ith coordinate.
     */
    float &
    operator[](int i)
    {
      return coordinate[i];
    }

    /**
     * \brief Get a const reference to the ith coordinate.
     */
    const float &
    operator[](int i) const
    {
      return coordinate[i];
    }

    /**
     * \brief Friend function for output stream.
     */
    template <int D>
    friend std::ostream &
    operator<<(std::ostream &outstream, const Coordinate<D> &coord);

  private:
    float coordinate[dim];
  };

  template <int D>
  std::ostream &
  operator<<(std::ostream &outstream, const Coordinate<D> &coord)
  {
    for (int i = 0; i < coord.size(); i++)
      {
        outstream << coord[i];
        if (i < coord.size() - 1)
          {
            outstream << " ";
          }
      }
    return outstream;
  }
} // namespace PRISMS

#endif