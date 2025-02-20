// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#ifndef element_volume_h
#define element_volume_h

#include <deal.II/base/aligned_vector.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/vectorization.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * \brief Compute the element volume.
 */
template <int dim, int degree, typename number>
class elementVolume
{
public:
  /**
   * \brief Constructor.
   */
  elementVolume() = default;

  /**
   * \brief Destructor.
   */
  ~elementVolume() = default;

  /**
   * \brief Initialize.
   */
  void
  initialize(std::shared_ptr<dealii::MatrixFree<dim, number>> _data);

  /**
   * \brief Compute element volume for the triangulation
   */
  void
  compute_element_volume(const dealii::FESystem<dim> &fe);

  /**
   * \brief Vector that stores element volumes
   */
  dealii::AlignedVector<dealii::VectorizedArray<number>> element_volume;

private:
  /**
   * \brief Matrix-free object.
   */
  std::shared_ptr<dealii::MatrixFree<dim, number>> data;
};

template <int dim, int degree, typename number>
void
elementVolume<dim, degree, number>::initialize(
  std::shared_ptr<dealii::MatrixFree<dim, number>> _data)
{
  data = _data;
}

template <int dim, int degree, typename number>
void
elementVolume<dim, degree, number>::compute_element_volume(
  const dealii::FESystem<dim> &fe)
{
  // Get the number of cell batches. Note this is the same as the cell range in
  // cell_loop()
  const unsigned int n_cells = data->n_cell_batches();

  // Resize vector
  element_volume.resize(n_cells);

  // Set quadrature rule and FEValues to update the JxW values
  dealii::QGaussLobatto<dim> quadrature(degree + 1);
  dealii::FEValues<dim>      fe_values(fe, quadrature, dealii::update_JxW_values);

  // Get the number of quadrature points
  const unsigned int num_quad_points = quadrature.size();

  // Loop over the cells and each lane in the vectorized array
  for (unsigned int cell = 0; cell < n_cells; cell++)
    {
      for (unsigned int lane = 0; lane < data->n_active_entries_per_cell_batch(cell);
           lane++)
        {
          // Get the iterator for the current cell
          auto cell_iterator = data->get_cell_iterator(cell, lane);

          // Reinitialize the cell
          fe_values.reinit(cell_iterator);

          // Initialize volume to 0 for the current cell
          double cell_volume = 0.0;

          // Sum up the JxW values at each quadrature point to compute the element volume
          // in 3D or area in 2D.
          for (unsigned int q_point = 0; q_point < num_quad_points; ++q_point)
            {
              cell_volume += fe_values.JxW(q_point);
            }

          // Store the element volume
          element_volume[cell][lane] = cell_volume;
        }
    }
}

PRISMS_PF_END_NAMESPACE

#endif