#ifndef Bin_HH
#define Bin_HH

#include <cmath>
#include <field_input/IntegrationTools/datastruc/PNDArray.hh>
#include <stdexcept>

/**
 * \brief A class for binning things.
 *
 * \tparam T typename thats being binned.
 * \tparam Coordinate class.
 */
template <class T, class Coordinate>
class Bin
{
private:
  PNDArray<std::vector<T>> _item;

  std::vector<double> _min;  // min coord
  std::vector<double> _incr; // histogram spacing along each direction (uniform along
                             // each dimension)
  std::vector<int>    _N;    // number of bins along each direction
  std::vector<double> _max;  // max coord
  std::vector<int>    _indices;

public:
  Bin() = default;

  /// Construct a Bin
  ///   'min': minimum value of each coordinate component
  ///   'incr': bin spacing along each direction (this is uniform along each dimension)
  ///   'N': number of bins along each direction
  ///
  /// For example, to bin the range (0->10, 0->20, 10->100), with size 1 bin spacing:
  ///    Bin<T>( {0,0,10}, {1,1,1}, {10, 20, 90})
  ///
  Bin(const std::vector<double> &min, std::vector<double> &incr, std::vector<int> &N)
  {
    _min  = min;
    _incr = incr;
    _N    = N;
    _max  = std::vector<double>(_min.size());
    for (unsigned int i = 0; i < _min.size(); i++)
      _max[i] = _min[i] + _incr[i] * _N[i];
    _indices = _N;
    _item.resize(_N);
  }

  void
  clear()
  {
    (*this) = Bin<T, Coordinate>();
  }

  std::vector<double> &
  min()
  {
    return _min;
  }

  std::vector<double> &
  max()
  {
    return _max;
  }

  /// Add a new item to the bin containing a given coordinate
  void
  add(const T &newitem, const Coordinate &coord)
  {
    indices(coord, _indices);
    _item(_indices).push_back(newitem);
  }

  /// Add a new item to the bin containing a given coordinate,
  ///   if an equivalent item is not already present
  /// return 'true' if added succesfully, 'false' if not
  void
  add_once(const T &newitem, const Coordinate &coord)
  {
    indices(coord, _indices);
    std::vector<T> &singlebin = _item(_indices);
    for (int i = 0; i < singlebin.size(); i++)
      {
        if (singlebin[i] == newitem)
          {
            return;
          }
      }

    singlebin.push_back(newitem);
  }

  /// Add a new item to all bins that fall in cuboid defined by 'min' and 'max'
  /// Coordinates
  template <class PCoord>
  void
  add_range(const T &newitem, const PCoord &min, const PCoord &max)
  {
    int              i;
    std::vector<int> index_min(_item.order());
    std::vector<int> index_max(_item.order());
    indices(min, index_min);
    indices(max, index_max);

    _indices = index_min;

    while (true)
      {
        _item(_indices).push_back(newitem);

        i = 0;
        while (_indices[i] == index_max[i])
          {
            i++;
            if (i == _item.order())
              {
                return;
              }
          }

        _indices[i]++;

        for (i = i - 1; i >= 0; i--)
          {
            _indices[i] = index_min[i];
          }
      }
  }

  std::vector<T> &
  contents(const Coordinate &coord)
  {
    indices(coord, _indices);
    return _item(_indices);
  }

  // maximum size of any bin
  int
  max_size()
  {
    unsigned int max = 0;
    for (int i = 0; i < _item.volume(); i++)
      {
        if (_item(i).size() > max)
          {
            max = _item(i).size();
          }
      }
    return max;
  }

private:
  /// Set 'term' to be the indices into '_item' PNDArray of the bin that contains
  /// 'coord'
  ///   Return 'false' if unsuccesful, 'true' if succesful
  template <class PCoord>
  void
  indices(const PCoord &coord, std::vector<int> &term)
  {
    for (int i = 0; i < _item.order(); i++)
      {
        if ((coord[i] < _min[i]) || (coord[i] > _max[i]))
          throw std::domain_error("Invalid coord, out of bin range");
        term[i] = std::floor((coord[i] - _min[i]) / _incr[i]);
      }
  }

  /// Set 'term' to be the indices into '_item' PNDArray of the bin that contains
  /// 'coord'
  ///   Return 'false' if unsuccesful, 'true' if succesful
  template <class PCoord>
  void
  max_indices(const PCoord &coord, std::vector<int> &term)
  {
    for (int i = 0; i < _item.order(); i++)
      {
        if ((coord[i] < _min[i]) || (coord[i] > _max[i]))
          throw std::domain_error("Invalid coord, out of bin range");

        if (std::floor((coord[i] - _min[i]) / _incr[i]) ==
            (coord[i] - _min[i]) / _incr[i])
          {
            term[i] = std::floor((coord[i] - _min[i]) / _incr[i]);
            term[i]--;
          }
        else
          {
            term[i] = std::floor((coord[i] - _min[i]) / _incr[i]);
          }
      }
  }
};

#endif
