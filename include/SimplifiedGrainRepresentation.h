#ifndef INCLUDE_SIMPLIFIEDGRAINREPRESENTATION_H_
#define INCLUDE_SIMPLIFIEDGRAINREPRESENTATION_H_

/**
* This class converts lists of grains and the vertices inside the grains to a simplified representation (currently spheres, other representations may be added later). Currently, assumptions are made that the elements are rectangular prisms. If not, a valid representation will still be made, but the centroid might be suboptimally placed.
*/
template <int dim>
class SimplifiedGrainRepresentation
{
public:
    SimplifiedGrainRepresentation(const GrainSet<dim> & grain_set);


    dealii::Point<dim> getCenter() const;
    double getRadius() const;

    unsigned int getGrainId() const;
    void setGrainId(unsigned int id);

    unsigned int getOrderParameterId() const;
    void setOrderParameterId(unsigned int id);

    unsigned int getOldOrderParameterId() const;
    void setOldOrderParameterId(unsigned int id);

protected:
    dealii::Point<dim> center;
    double radius;
    unsigned int grain_id;
    unsigned int order_parameter_id;
    unsigned int old_order_parameter_id;
};

/**
* This class...
*/
template <int dim>
class SimplifiedGrainManipulator
{
public:
    void reassignGrains(std::vector<SimplifiedGrainRepresentation<dim>> & grain_representations, double buffer_distance, std::vector<unsigned int> order_parameter_id_list);

    void transferGrainIds(const std::vector<SimplifiedGrainRepresentation<dim>> & old_grain_representations, std::vector<SimplifiedGrainRepresentation<dim>> & new_grain_representations) const;


protected:

};

#endif
