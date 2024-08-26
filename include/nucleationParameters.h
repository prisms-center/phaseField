#ifndef INCLUDE_NUCLEATIONPARAMETERS_H_
#define INCLUDE_NUCLEATIONPARAMETERS_H_

template <int dim>
class nucleationParameters
{
public:
  unsigned int                   var_index;
  std::vector<double>            semiaxes;
  std::vector<double>            ellipsoid_rotation;
  std::vector<double>            freeze_semiaxes;
  double                         no_nucleation_border_thickness;
  double                         hold_time;
  dealii::Tensor<2, dim, double> rotation_matrix;

  nucleationParameters(unsigned int        _var_index,
                       std::vector<double> _semiaxes,
                       std::vector<double> _freeze_semiaxes,
                       std::vector<double> _ellipsoid_rotation,
                       double              _hold_time,
                       double              _no_nucleation_border_thickness)
  {
    var_index                      = _var_index;
    semiaxes                       = _semiaxes;
    freeze_semiaxes                = _freeze_semiaxes;
    ellipsoid_rotation             = _ellipsoid_rotation;
    hold_time                      = _hold_time;
    no_nucleation_border_thickness = _no_nucleation_border_thickness;

    set_rotation_matrix();
  };

  void
  set_rotation_matrix()
  {
    // Rotation conventions:
    // Rx refers to rotations about the x axis, Ry refers to rotations about the
    // y axis, and Rz refers to rotations about the z axis. A positive rotation
    // angle about the z axis corresponds to a clockwise rotation of the
    // particle in a 2D calculation.

    double degrees_to_rad = std::acos(0.0) / 90.0;

    dealii::Tensor<2, dim, double> Rx, Ry, Rz;

    Rx[0][0] = 1.0;
    Rx[1][1] = std::cos(ellipsoid_rotation.at(0) * degrees_to_rad);

    Ry[0][0] = std::cos(ellipsoid_rotation.at(1) * degrees_to_rad);
    Ry[1][1] = 1.0;

    Rz[0][0] = std::cos(ellipsoid_rotation.at(2) * degrees_to_rad);
    Rz[1][0] = std::sin(ellipsoid_rotation.at(2) * degrees_to_rad);
    Rz[0][1] = -std::sin(ellipsoid_rotation.at(2) * degrees_to_rad);
    Rz[1][1] = std::cos(ellipsoid_rotation.at(2) * degrees_to_rad);

    if (dim == 3)
      {
        Rx[1][2] = -std::sin(ellipsoid_rotation.at(0) * degrees_to_rad);
        Rx[2][1] = std::sin(ellipsoid_rotation.at(0) * degrees_to_rad);
        Rx[2][2] = std::cos(ellipsoid_rotation.at(0) * degrees_to_rad);

        Ry[0][2] = std::sin(ellipsoid_rotation.at(1) * degrees_to_rad);
        Ry[2][0] = -std::sin(ellipsoid_rotation.at(1) * degrees_to_rad);
        Ry[2][2] = std::cos(ellipsoid_rotation.at(1) * degrees_to_rad);

        Rz[2][2] = 1.0;
      }

    rotation_matrix = Rx * Ry * Rz; // Note: these are tensor multiplications
  };
};

#endif
