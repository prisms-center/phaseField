// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <deal.II/base/exceptions.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/vectorization.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_update_flags.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <prismspf/core/conditional_ostreams.h>
#include <prismspf/core/matrix_free_handler.h>
#include <prismspf/core/multigrid_info.h>
#include <prismspf/core/timer.h>

#include <prismspf/utilities/element_volume.h>

#include <prismspf/config.h>

#include <memory>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree, typename number>
void
ElementVolume<dim, degree, number>::initialize(
  std::shared_ptr<dealii::MatrixFree<dim, number, dealii::VectorizedArray<number>>> _data)
{
  data = std::move(_data);
}

template <unsigned int dim, unsigned int degree, typename number>
void
ElementVolume<dim, degree, number>::compute_element_volume()
{
  Assert(data != nullptr, dealii::ExcNotInitialized());
  // Get the number of cell batches. Note this is the same as the cell range in
  // cell_loop()
  const unsigned int n_cells = data->n_cell_batches();

  // Resize vector
  element_volume.resize(n_cells);

  // Set quadrature rule and FEValues to update the JxW values
  const dealii::QGaussLobatto<dim> quadrature(degree + 1);
  dealii::FEValues<dim>            fe_values(data->get_dof_handler().get_fe(),
                                  quadrature,
                                  dealii::update_JxW_values);

  // Get the number of quadrature points
  const unsigned int num_quad_points = quadrature.size();

  // Loop over the cells and each lane in the vectorized array
  for (unsigned int cell = 0; cell < n_cells; cell++)
    {
      for (unsigned int lane = 0; lane < data->n_active_entries_per_cell_batch(cell);
           lane++)
        {
          // Get the iterator for the cell
          auto cell_iterator = data->get_cell_iterator(cell, lane);

          // Reinitialize the cell
          fe_values.reinit(cell_iterator);

          // Initialize volume to 0 for the cell
          number cell_volume = 0.0;

          // Sum up the JxW values at each quadrature point to compute the element volume
          // in 3D or area in 2D.
          for (unsigned int q_point = 0; q_point < num_quad_points; ++q_point)
            {
              cell_volume += fe_values.JxW(q_point);
            }

          // Store the element volume
          std::cout << cell_volume << std::endl;
          element_volume[cell][lane] = cell_volume;
        }
    }
}

template <unsigned int dim, unsigned int degree, typename number>
const dealii::AlignedVector<dealii::VectorizedArray<number>> &
ElementVolume<dim, degree, number>::get_element_volumes() const
{
  Assert(data != nullptr, dealii::ExcNotInitialized());
  Assert(!element_volume.empty(),
         dealii::ExcMessage("The vector of element volumes is empty. Make sure to call "
                            "`compute_element_volume` beforehand"));
  return element_volume;
};

template <unsigned int dim, unsigned int degree, typename number>
const dealii::VectorizedArray<number> &
ElementVolume<dim, degree, number>::get_element_volume(unsigned cell) const
{
  Assert(data != nullptr, dealii::ExcNotInitialized());
  Assert(!element_volume.empty(),
         dealii::ExcMessage("The vector of element volumes is empty. Make sure to call "
                            "`compute_element_volume` beforehand"));
  Assert(element_volume.size() > cell,
         dealii::ExcIndexRange(cell, 0, element_volume.size()));
  return element_volume[cell];
};

template <unsigned int dim, unsigned int degree, typename number>
ElementVolumeContainer<dim, degree, number>::ElementVolumeContainer(MGInfo<dim> &mg_info)
  : element_volume()
  , multigrid_element_volume(0, 0)
{
  // If we have multigrid, we need to set the min and max levels, and reinitialize the
  // multigrid element volumes.
  if (mg_info.has_multigrid())
    {
      min_level = mg_info.get_mg_min_level();
      max_level = mg_info.get_mg_max_level();
      multigrid_element_volume.resize(min_level, max_level);
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
ElementVolumeContainer<dim, degree, number>::initialize(
  const MatrixFreeContainer<dim, number> &matrix_free_container)
{
  Timer::start_section("reinitialize element volumes");

  ConditionalOStreams::pout_base() << "initializing element volumes...\n" << std::flush;
  element_volume.initialize(matrix_free_container.get_matrix_free());
  if (multigrid_element_volume.n_levels() > 1)
    {
      for (unsigned int level = min_level; level <= max_level; ++level)
        {
          ConditionalOStreams::pout_base()
            << "initializing multgrid element volumes at level " << level << "...\n"
            << std::flush;
          multigrid_element_volume[level].initialize(
            matrix_free_container.get_mg_matrix_free(level));
        }
    }

  Timer::end_section("reinitialize element volumes");
}

template <unsigned int dim, unsigned int degree, typename number>
void
ElementVolumeContainer<dim, degree, number>::compute_element_volume()
{
  Timer::start_section("compute element volumes");

  ConditionalOStreams::pout_base() << "computing element volumes...\n" << std::flush;
  element_volume.compute_element_volume();
  if (multigrid_element_volume.n_levels() > 1)
    {
      for (unsigned int level = min_level; level <= max_level; ++level)
        {
          ConditionalOStreams::pout_base()
            << "computing multgrid element volumes at level " << level << "...\n"
            << std::flush;
          multigrid_element_volume[level].compute_element_volume();
        }
    }

  Timer::end_section("compute element volumes");
}

template <unsigned int dim, unsigned int degree, typename number>
void
ElementVolumeContainer<dim, degree, number>::recompute_element_volume()
{
  this->compute_element_volume();
}

#include "utilities/element_volume.inst"

PRISMS_PF_END_NAMESPACE
