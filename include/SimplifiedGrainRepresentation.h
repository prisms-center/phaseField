#ifndef INCLUDE_SIMPLIFIEDGRAINREPRESENTATION_H_
#define INCLUDE_SIMPLIFIEDGRAINREPRESENTATION_H_

#include "FloodFiller.h"

/**
 * This class converts lists of grains and the vertices inside the grains to a
 * simplified representation (currently spheres, other representations may be
 * added later). Currently, assumptions are made that the elements are
 * rectangular prisms. If not, a valid representation will still be made, but
 * the centroid might be suboptimally placed.
 */
template <int dim>
class SimplifiedGrainRepresentation
{
public:
  /**
   * Constructor. Creates the simplified representation of a grain from its
   * GrainSet. This sets all of the internal data members and is the only way
   * to set the radius and center of the grain. The members order_parameter_id
   * and old_order_parameter_id are initialized to the same value.
   */
  SimplifiedGrainRepresentation(const GrainSet<dim> &grain_set);

  /**
   * Getter for the grain center/centroid.
   */
  dealii::Point<dim>
  getCenter() const;

  /**
   * Getter for the grain radius.
   */
  double
  getRadius() const;

  /**
   * Getter for the grain id.
   */
  unsigned int
  getGrainId() const;

  /**
   * Setter for the grain id.
   */
  void
  setGrainId(unsigned int id);

  /**
   * Getter for the order parameter id.
   */
  unsigned int
  getOrderParameterId() const;

  /**
   * Setter for the order parameter id.
   */
  void
  setOrderParameterId(unsigned int id);

  /**
   * Getter for the old value of the order parameter id (used in the
   * transferGrainIds method of the SimplifiedGrainManipulator class).
   */
  unsigned int
  getOldOrderParameterId() const;

  /**
   * Setter for the distance from this grain to the nearest grain with the same
   * order parameter.
   */
  void
  setDistanceToNeighbor(double dist);

  /**
   * Getter for the distance from this grain to the nearest grain with the same
   * order parameter
   */
  double
  getDistanceToNeighbor() const;

protected:
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

/**
 * This class contains methods to make changes to lists of
 * SimplifiedGrainRepresentation objects to aid in order parameter remapping.
 * Currently, this is a stub class containing two unrelated methods and no
 * internal members. If this stays the case, this class may be absorbed into
 * another one.
 */
template <int dim>
class SimplifiedGrainManipulator
{
public:
  /**
   * This method checks for collisions between SimplifiedGrainRepresentation
   * objects with the same order parameter and reassigns them, if needed.
   */
  void
  reassignGrains(std::vector<SimplifiedGrainRepresentation<dim>> &grain_representations,
                 double                                           buffer_distance,
                 std::vector<unsigned int> order_parameter_id_list);

  /**
   * This method checks the centers of two lists of
   * SimplifiedGrainRepresentation objects from different times in the
   * simulation and reassigns the grain ids so that they consistently refer to
   * the same grains.
   */
  void
  transferGrainIds(
    const std::vector<SimplifiedGrainRepresentation<dim>> &old_grain_representations,
    std::vector<SimplifiedGrainRepresentation<dim>> &new_grain_representations) const;

protected:
};

#endif
