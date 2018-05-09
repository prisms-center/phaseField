#ifndef INCLUDE_FLOODFILLER_H_
#define INCLUDE_FLOODFILLER_H_

#include "dealIIheaders.h"

template <int dim>
class GrainSet
{
public:
    void setGrainIndex(unsigned int _grain_index){grain_index = _grain_index;};

    unsigned int getGrainIndex(){return grain_index;};

    void addVertexList(std::vector<dealii::Point<dim>> _vertices){list_of_vertices.push_back(_vertices);};

    std::vector<std::vector<dealii::Point<dim>>> getVertexList(){return list_of_vertices;};

private:
    unsigned int grain_index;
    std::vector<std::vector<dealii::Point<dim>>> list_of_vertices;
};



/**
* This class uses a recursive flood filling algorithm to find connected bodies in a field, given a threshold
*/
template <int dim, int degree>
class FloodFiller
{
public:
    FloodFiller(FESystem<dim> & _fe, QGaussLobatto<dim> _quadrature): quadrature(_quadrature), num_quad_points(_quadrature.size()), dofs_per_cell(_fe.dofs_per_cell){
        fe = & _fe;
    };

    void calcGrainSets(FESystem<dim> & fe, dealii::DoFHandler<dim> &dof_handler, vectorType* solution_field, double threshold, std::vector<GrainSet<dim>> & grain_sets);
protected:

    template <typename T>
    void recursiveFloodFill(T di, T di_end, vectorType* solution_field, double threshold, unsigned int & grain_index, std::vector<GrainSet<dim>> & grain_sets, bool & grain_assigned);

    void communicateGrainSets(std::vector<GrainSet<dim>> & grain_sets);

    void sendUpdate (int procno, std::vector<GrainSet<dim>> & grain_sets) const;

    void receiveUpdate (int procno, std::vector<GrainSet<dim>> & grain_sets) const;

    void broadcastUpdate (int broadcastProc, int thisProc, std::vector<GrainSet<dim>> & grain_sets) const;

    /**
    * Checks to see if grains found on different processors are parts of a larger grain. If so, it merges the grain_sets entries.
    */
    void mergeSplitGrains (std::vector<GrainSet<dim>> & grain_sets) const;

    QGaussLobatto<dim>  quadrature;
    const unsigned int   num_quad_points;
    const unsigned int   dofs_per_cell;

    FESystem<dim> * fe;
};

#endif
