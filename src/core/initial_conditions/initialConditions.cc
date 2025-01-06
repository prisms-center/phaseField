#include <cmath>
#include <core/exceptions.h>
#include <core/initial_conditions/initialConditions.h>
#include <core/matrixFreePDE.h>
#include <field_input/IntegrationTools/PField.hh>
#include <grains/OrderParameterRemapper.h>

template <int dim>
class InitialConditionPField : public Function<dim>
{
public:
  unsigned int   index;
  Vector<double> values;

  using ScalarField = PRISMS::PField<double *, double, dim>;
  ScalarField &inputField;

  InitialConditionPField(const unsigned int _index, ScalarField &_inputField)
    : Function<dim>(1)
    , index(_index)
    , inputField(_inputField)
  {}

  [[nodiscard]] double
  value(const Point<dim>                   &p,
        [[maybe_unused]] const unsigned int component = 0) const override
  {
    double scalar_IC = NAN;

    double coord[dim];
    for (unsigned int i = 0; i < dim; i++)
      {
        coord[i] = p(i);
      }

    scalar_IC = inputField(coord);

    return scalar_IC;
  }
};

// methods to apply initial conditions
template <int dim, int degree>
void
MatrixFreePDE<dim, degree>::applyInitialConditions()
{
  if (userInputs.load_grain_structure)
    {
      // Create the dummy field
      dealii::LinearAlgebra::distributed::Vector<double> grain_index_field;

      // Clear the order parameter fields
      unsigned int op_list_index = 0;
      for (unsigned int var_index = 0; var_index < var_attributes.size(); var_index++)
        {
          if (op_list_index < userInputs.variables_for_remapping.size())
            {
              if (var_index == userInputs.variables_for_remapping.at(op_list_index))
                {
                  *solutionSet[var_index] = 0.0;
                  op_list_index++;
                }
            }
        }

      simplified_grain_representations.clear();

      // Get the index of one of the scalar fields
      unsigned int scalar_field_index = 0;
      for (const auto &[index, variable] : var_attributes)
        {
          if (variable.var_type == SCALAR)
            {
              scalar_field_index = index;
              break;
            }
        }

      matrixFreeObject.initialize_dof_vector(grain_index_field, scalar_field_index);

      // Declare the PField types and containers
      typedef PRISMS::PField<double *, double, dim> ScalarField;
      typedef PRISMS::Body<double *, dim>           Body;
      Body                                          body;

      // Create the filename of the the file to be loaded

      // Load the data from the file using a PField
      std::string filename = userInputs.grain_structure_filename;
      filename += ".vtk";

      // new section added for the choice of unstructured mesh and rectilinear mesh
      if (userInputs.load_vtk_file_type == "UNSTRUCTURED")
        {
          body.read_vtk(filename);
        }
      else if (userInputs.load_vtk_file_type == "RECTILINEAR")
        {
          body.read_RL_vtk(filename);
        }
      else
        {
          pcout << "Error in vtk file type: Use either UNSTRUCTURED OR RECTILINEAR\n";
          abort();
        } // new section ends

      ScalarField &id_field =
        body.find_scalar_field(userInputs.grain_structure_variable_name);

      pcout << "Applying PField initial condition...\n";

      VectorTools::interpolate(*dofHandlersSet[scalar_field_index],
                               InitialConditionPField<dim>(0, id_field),
                               grain_index_field);

      grain_index_field.update_ghost_values();

      // Get the max and min grain ids
      auto         max_id = (unsigned int) grain_index_field.linfty_norm();
      unsigned int min_id = 0;

      // Now locate all of the grains and create simplified representations of
      // them
      QGaussLobatto<dim>       quadrature2(degree + 1);
      FloodFiller<dim, degree> flood_filler(*FESet.at(scalar_field_index), quadrature2);

      pcout << "Locating the grains...\n";
      std::vector<GrainSet<dim>> grain_sets;
      for (unsigned int id = min_id; id < max_id + 1; id++)
        {
          pcout << "Locating grain " << id << "...\n";

          std::vector<GrainSet<dim>> grain_sets_single_id;

          flood_filler.calcGrainSets(*FESet.at(scalar_field_index),
                                     *dofHandlersSet_nonconst.at(scalar_field_index),
                                     &grain_index_field,
                                     (double) id - userInputs.order_parameter_threshold,
                                     (double) id + userInputs.order_parameter_threshold,
                                     min_id,
                                     0,
                                     grain_sets_single_id);

          for (unsigned int g = 0; g < grain_sets_single_id.size(); g++)
            {
              grain_sets_single_id.at(g).setGrainIndex(id);
            }

          grain_sets.insert(grain_sets.end(),
                            grain_sets_single_id.begin(),
                            grain_sets_single_id.end());
        }

      pcout << "Generating simplified representations of the grains...\n";
      for (unsigned int g = 0; g < grain_sets.size(); g++)
        {
          SimplifiedGrainRepresentation<dim> simplified_grain_representation(
            grain_sets.at(g));

          if (dim == 2)
            {
              pcout << "Grain: " << simplified_grain_representation.getGrainId() << " "
                    << simplified_grain_representation.getOrderParameterId()
                    << " Center: " << simplified_grain_representation.getCenter()(0)
                    << " " << simplified_grain_representation.getCenter()(1)
                    << "  Radius: " << simplified_grain_representation.getRadius()
                    << std::endl;
            }
          else
            {
              pcout << "Grain: " << simplified_grain_representation.getGrainId() << " "
                    << simplified_grain_representation.getOrderParameterId()
                    << " Center: " << simplified_grain_representation.getCenter()(0)
                    << " " << simplified_grain_representation.getCenter()(1) << " "
                    << simplified_grain_representation.getCenter()(2)
                    << "  Radius: " << simplified_grain_representation.getRadius()
                    << std::endl;
            }

          simplified_grain_representations.push_back(simplified_grain_representation);
        }

      // Delete grains with very small radii that correspond to a single
      // interfacial element
      for (unsigned int g = 0; g < simplified_grain_representations.size(); g++)
        {
          if (simplified_grain_representations.at(g).getRadius() <
              userInputs.min_radius_for_loading_grains)
            {
              simplified_grain_representations.erase(
                simplified_grain_representations.begin() + g);
              g--;
            }
        }

      pcout << "Reassigning the grains to new order parameters...\n";
      SimplifiedGrainManipulator<dim> simplified_grain_manipulator;
      simplified_grain_manipulator.reassignGrains(simplified_grain_representations,
                                                  userInputs.buffer_between_grains,
                                                  userInputs.variables_for_remapping);

      pcout << "After reassignment: " << std::endl;
      for (unsigned int g = 0; g < simplified_grain_representations.size(); g++)
        {
          if (dim == 2)
            {
              pcout << "Grain: " << simplified_grain_representations.at(g).getGrainId()
                    << " " << simplified_grain_representations.at(g).getOrderParameterId()
                    << " Center: "
                    << simplified_grain_representations.at(g).getCenter()(0) << " "
                    << simplified_grain_representations.at(g).getCenter()(1) << std::endl;
            }
          else
            {
              pcout << "Grain: " << simplified_grain_representations.at(g).getGrainId()
                    << " " << simplified_grain_representations.at(g).getOrderParameterId()
                    << " Center: "
                    << simplified_grain_representations.at(g).getCenter()(0) << " "
                    << simplified_grain_representations.at(g).getCenter()(1) << " "
                    << simplified_grain_representations.at(g).getCenter()(2) << std::endl;
            }
        }

      pcout << "Placing the grains in their new order parameters...\n";
      OrderParameterRemapper<dim> order_parameter_remapper;
      order_parameter_remapper.remap_from_index_field(
        simplified_grain_representations,
        &grain_index_field,
        solutionSet,
        *dofHandlersSet_nonconst.at(scalar_field_index),
        FESet.at(scalar_field_index)->dofs_per_cell);

      // Smooth the order parameters according to Fick's 2nd Law
      // In the time cycle below, we evolve the weak form of Eq.:
      // field_i^(n+1)=field_i^(n)+D*dt*Laplacian(field_i^(n)) Old
      // dt_for_smoothing double dt_for_smoothing =
      // dealii::GridTools::minimal_cell_diameter(triangulation)/1000.0; NEW
      // selection of dt_for_smoothing based on dt_max=(dx^2)/2*dim*D, where
      // D=1, addording to the von Neumann stability condition We chose
      // dt_for_smoothing=0.25*dt_max
      double min_dx = dealii::GridTools::minimal_cell_diameter(triangulation) /
                      (1.0 * degree * std::sqrt(1.0 * dim));
      double dt_for_smoothing = 0.25 * min_dx * min_dx / (1.0 * dim);

      op_list_index = 0;
      for (unsigned int fieldIndex = 0; fieldIndex < fields.size(); fieldIndex++)
        {
          if (op_list_index < userInputs.variables_for_remapping.size())
            {
              if (fieldIndex == userInputs.variables_for_remapping.at(op_list_index))
                {
                  for (unsigned int cycle = 0;
                       cycle < userInputs.num_grain_smoothing_cycles;
                       cycle++)
                    {
                      // Calculates the Laplace RHS and stores the information
                      // in residualSet
                      computeLaplaceRHS(fieldIndex);
                      if (fields[fieldIndex].type == SCALAR)
                        {
                          unsigned int invM_size = invMscalar.locally_owned_size();
                          for (unsigned int dof = 0;
                               dof < solutionSet[fieldIndex]->locally_owned_size();
                               ++dof)
                            {
                              solutionSet[fieldIndex]->local_element(dof) =
                                solutionSet[fieldIndex]->local_element(dof) -
                                invMscalar.local_element(dof % invM_size) *
                                  residualSet[fieldIndex]->local_element(dof) *
                                  dt_for_smoothing;
                            }
                        }
                      else if (fields[fieldIndex].type == VECTOR)
                        {
                          unsigned int invM_size = invMvector.locally_owned_size();
                          for (unsigned int dof = 0;
                               dof < solutionSet[fieldIndex]->locally_owned_size();
                               ++dof)
                            {
                              solutionSet[fieldIndex]->local_element(dof) =
                                solutionSet[fieldIndex]->local_element(dof) -
                                invMvector.local_element(dof % invM_size) *
                                  residualSet[fieldIndex]->local_element(dof) *
                                  dt_for_smoothing;
                            }
                        }

                      solutionSet[fieldIndex]->update_ghost_values();
                    }

                  op_list_index++;
                }
            }
        }
    }

  // Create map of vtk files that are read in
  std::unordered_map<std::string, std::vector<size_t>> file_field_map;
  for (size_t i = 0; i < var_attributes.size(); ++i)
    {
      if (userInputs.load_ICs[i])
        {
          file_field_map[userInputs.load_file_name[i]].push_back(i);
        }
    }
  // Read out unique filenames
  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0 && !file_field_map.empty())
    {
      std::cout << "Unique VTK input files: " << std::endl;
      for (const auto &pair : file_field_map)
        {
          std::cout << pair.first << ", ";
        }
    }
  // Read in each vtk once and apply initial conditions
  using ScalarField = PRISMS::PField<double *, double, dim>;
  using Body        = PRISMS::Body<double *, dim>;
  Body body;

  for (const auto &pair : file_field_map)
    {
      bool                using_parallel_files = false;
      std::string         filename             = pair.first;
      std::vector<size_t> index_list           = pair.second;
      // For parallel file capability
      for (const auto &index : index_list)
        {
          using_parallel_files =
            using_parallel_files || userInputs.load_parallel_file[index];
        }
      if (using_parallel_files)
        {
          std::ostringstream conversion;
          conversion << Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
          filename.append("." + conversion.str() + ".vtk");
        }
      else
        {
          filename.append(".vtk");
        }

      std::cout << "Reading " << filename << "\n";
      // Load the data from the file using a PField
      // new section added for the choice of unstructured mesh and rectilinear mesh
      if (userInputs.load_vtk_file_type == "UNSTRUCTURED")
        {
          body.read_vtk(filename);
        }
      else if (userInputs.load_vtk_file_type == "RECTILINEAR")
        {
          body.read_RL_vtk(filename);
        }
      else
        {
          pcout << "Error in vtk file type: Use either UNSTRUCTURED OR RECTILINEAR\n";
          abort();
        } // new section ends

      for (const auto &index : index_list)
        {
          std::string var_name = userInputs.load_field_name[index];

          // Find the scalar field in the file
          ScalarField &field = body.find_scalar_field(var_name);

          if (var_attributes.at(index).var_type == SCALAR)
            {
              pcout << "Applying PField initial condition for "
                    << userInputs.load_field_name[index] << "...\n";
              VectorTools::interpolate(*dofHandlersSet[index],
                                       InitialConditionPField<dim>(index, field),
                                       *solutionSet[index]);
            }
          else
            {
              AssertThrow(false,
                          ExcMessage("PRISMS-PF Error: We do not support the loading of "
                                     "vector fields from vtks at this moment."));
            }
        }
    }

  unsigned int op_list_index = 0;
  for (const auto &[var_index, variable] : var_attributes)
    {
      bool is_remapped_op = false;
      if (op_list_index < userInputs.variables_for_remapping.size())
        {
          if (var_index == userInputs.variables_for_remapping.at(op_list_index))
            {
              is_remapped_op = true;
              op_list_index++;
            }
        }

      if (!is_remapped_op || !userInputs.load_grain_structure)
        {
          if (userInputs.load_ICs[var_index] == false)
            {
              pcout << "Applying non-PField initial condition...\n";

              if (variable.var_type == SCALAR)
                {
                  VectorTools::interpolate(*dofHandlersSet[var_index],
                                           InitialCondition<dim, degree>(var_index,
                                                                         userInputs,
                                                                         this),
                                           *solutionSet[var_index]);
                }
              else if (variable.var_type == VECTOR)
                {
                  VectorTools::interpolate(*dofHandlersSet[var_index],
                                           InitialConditionVector<dim, degree>(var_index,
                                                                               userInputs,
                                                                               this),
                                           *solutionSet[var_index]);
                }
            }
          pcout << "Application of initial conditions for field number " << var_index
                << " complete \n";
        }
    }
}
