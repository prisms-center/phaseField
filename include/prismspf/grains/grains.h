// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <prismspf/grains/flood_filler.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * This class converts lists of grains and the vertices inside the grains to a
 * simplified representation (currently spheres, other representations may be
 * added later). Currently, assumptions are made that the elements are
 * rectangular prisms. If not, a valid representation will still be made, but
 * the centroid might be suboptimally placed.
 */
template <unsigned int dim>
class SimplifiedGrainRepresentation
{
public:
  /**
   * Constructor. Creates the simplified representation of a grain from its
   * GrainSet. This sets all of the internal data members and is the only way
   * to set the radius and center of the grain. The members order_parameter_id
   * and old_order_parameter_id are initialized to the same value.
   */
  explicit SimplifiedGrainRepresentation(const GrainSet<dim> &grain_set);

  /**
   * Getter for the grain center/centroid.
   */
  [[nodiscard]] dealii::Point<dim>
  get_center() const;

  /**
   * Getter for the grain radius.
   */
  [[nodiscard]] double
  get_radius() const;

  /**
   * Getter for the grain id.
   */
  [[nodiscard]] unsigned int
  get_grain_id() const;

  /**
   * Setter for the grain id.
   */
  void
  set_grain_id(unsigned int _grain_id);

  /**
   * Getter for the order parameter id.
   */
  [[nodiscard]] unsigned int
  get_order_parameter_id() const;

  /**
   * Setter for the order parameter id.
   */
  void
  set_order_parameter_id(unsigned int _order_parameter_id);

  /**
   * Getter for the old value of the order parameter id (used in the
   * transferGrainIds method of the SimplifiedGrainManipulator class).
   */
  [[nodiscard]] unsigned int
  get_old_order_parameter_id() const;

  /**
   * Setter for the distance from this grain to the nearest grain with the same
   * order parameter.
   */
  void
  set_distance_to_neighbor(double dist);

  /**
   * Getter for the distance from this grain to the nearest grain with the same
   * order parameter
   */
  [[nodiscard]] double
  get_distance_to_neighbor() const;

private:
  /**
   * The center of the circle/sphere that represents the grain.
   */
  dealii::Point<dim> center;

  /**
   * The radius of the circle/sphere that represents the grain.
   */
  double radius;

  /**
   * The id of grain, which is used to distinguish between grains and give them
   * different properties.
   */
  unsigned int grain_id;

  /**
   * The variable index of the order parameter containing this grain.
   */
  unsigned int order_parameter_id;

  /**
   * The variable index of the order parameter containing this grain before
   * reassignment (which is performed outside this class).
   */
  unsigned int old_order_parameter_id;

  /**
   * The distance from this grain to the nearest grain with the same order
   * parameter. This value is used to determine the cutoff distance for
   * tranfering the grain between order parameters.
   */
  double distance_to_neighbor_sharing_op;
};

PRISMS_PF_END_NAMESPACE