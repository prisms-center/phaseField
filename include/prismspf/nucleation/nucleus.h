// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/point.h>

#include <prismspf/core/types.h>

#include <prismspf/config.h>

#include <mpi.h>
#include <ostream>
#include <string>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * @brief This class contains mutable utilities for phase field problems.
 */
template <unsigned int dim>
struct Nucleus
{
public:
  /**
   * @brief Constructor.
   */
  Nucleus() = default;

  Nucleus(const unsigned int       &_field_index,
          const dealii::Point<dim> &_location,
          const double             &_seed_time,
          const unsigned int       &_seed_increment)
    : field_index(_field_index)
    , location(_location)
    , seed_time(_seed_time)
    , seed_increment(_seed_increment)
  {}

  dealii::Point<dim, dealii::VectorizedArray<double>>
  location_vectorized() const;

  unsigned int       field_index = 0;
  dealii::Point<dim> location;
  double             seed_time      = 0.0;
  unsigned int       seed_increment = 0;

  static MPI_Datatype
  mpi_datatype();
};

template <unsigned int dim>
inline dealii::Point<dim, dealii::VectorizedArray<double>>
Nucleus<dim>::location_vectorized() const
{
  dealii::Point<dim, dealii::VectorizedArray<double>> result;
  for (unsigned int d = 0; d < dim; ++d)
    {
      result[d] = dealii::VectorizedArray<double>(location[d]);
    }
  return result;
}

template <unsigned int dim, typename OStream>
OStream &
operator<<(OStream &ost, const Nucleus<dim> &nucleus)
{
  ost << "Field Index: " << std::to_string(nucleus.field_index)
      << " ; Location: " << nucleus.location
      << " ; Seed Time: " << std::to_string(nucleus.seed_time)
      << " ; Seed Increment: " << std::to_string(nucleus.seed_increment);
  return ost;
}

template <unsigned int dim>
inline MPI_Datatype
Nucleus<dim>::mpi_datatype()
{
  static MPI_Datatype MPI_NUCLEUS = MPI_DATATYPE_NULL;
  if (MPI_NUCLEUS != MPI_DATATYPE_NULL)
    {
      return MPI_NUCLEUS;
    }

  constexpr int block_lengths[4] = {1, static_cast<int>(dim), 1, 1};
  MPI_Datatype  types[4]         = {MPI_UNSIGNED, MPI_DOUBLE, MPI_DOUBLE, MPI_UNSIGNED};
  MPI_Aint      displacements[4];

  Nucleus<dim> dummy;
  MPI_Aint     base;
  MPI_Get_address(&dummy, &base);

  MPI_Get_address(&dummy.field_index, &displacements[0]);
  MPI_Get_address(&dummy.location[0], &displacements[1]);
  MPI_Get_address(&dummy.seed_time, &displacements[2]);
  MPI_Get_address(&dummy.seed_increment, &displacements[3]);

  for (int i = 0; i < 4; ++i)
    {
      displacements[i] -= base;
    }

  MPI_Type_create_struct(4, block_lengths, displacements, types, &MPI_NUCLEUS);
  MPI_Type_commit(&MPI_NUCLEUS);

  return MPI_NUCLEUS;
}

PRISMS_PF_END_NAMESPACE