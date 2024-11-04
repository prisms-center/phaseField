#ifndef Body_HH
#define Body_HH

#include "./Mesh.hh"
#include "./PField.hh"
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

namespace PRISMS
{

  /// A class for a Body: a combination of Mesh and Field(s))
  ///
  template <class Coordinate, int DIM>
  class Body
  {
  public:
    Mesh<Coordinate, DIM> mesh;

    std::vector<PField<Coordinate, double, DIM>> scalar_field;

    // std::vector< PField<Coordinate, std::vector<double>, DIM > > vector_field;

    // std::vector< PField<Coordinate, Tensor<double>, DIM > > tensor_field;

    // ----------------------------------------------------------
    // Constructors
    Body() {};

    /// Read from a 2D vtk file
    ///   For now:
    ///      only ASCII files
    ///      only rectilinear grids (though output as UNSTRUCTURED_GRID)
    ///      only (2d) Quad elements
    ///
    void
    read_vtk(const std::string &vtkfile)
    {
      std::cout << "Begin reading unstructured vtk file" << std::endl;

      // read in vtk file here
      std::ifstream infile_mesh(vtkfile.c_str());

      // read mesh info
      mesh.read_vtk(infile_mesh);

      std::ifstream infile(vtkfile.c_str());

      // read point data
      std::istringstream ss;
      std::string        str, name, type, line;
      int                numcomp;
      unsigned long int  Npoints, u, p;

      while (!infile.eof())
        {
          std::getline(infile, line);

          if (line[0] == 'P')
            {
              if (line.size() > 9 && line.substr(0, 10) == "POINT_DATA")
                {
                  // std::cout << line << "\n";
                  ss.clear();
                  ss.str(line);
                  ss >> str >> Npoints;
                }
            }

          if (line[0] == 'S')
            {
              if (line.size() > 6 && line.substr(0, 7) == "SCALARS")
                {
                  ss.clear();
                  ss.str(line);
                  ss >> str >> name >> type >> numcomp;

                  // read LOOKUP_TABLE line
                  std::getline(infile, line);

                  // read data
                  std::cout << "begin reading data" << std::endl;

                  std::vector<double> gid(Npoints);
                  for (unsigned int i = 0; i < Npoints; i++)
                    {
                      infile >> gid[i];
                      // std::cout << data[i] << std::endl;
                    }
                  std::cout << "  done" << std::endl;

                  // construct field
                  std::vector<std::string> var_name(DIM);
                  std::vector<std::string> var_description(DIM);

                  if (DIM == 2)
                    {
                      var_name[0]        = "x";
                      var_description[0] = "x coordinate";
                      var_name[1]        = "y";
                      var_description[1] = "y coordinate";
                    }
                  if (DIM > 2)
                    {
                      var_name[2]        = "z";
                      var_description[2] = "z coordinate";
                    }

                  std::cout << "Construct PField '" << name << "'" << std::endl;
                  scalar_field.push_back(PField<Coordinate, double, DIM>(name,
                                                                         var_name,
                                                                         var_description,
                                                                         mesh,
                                                                         gid,
                                                                         0.0));
                  std::cout << "  done" << std::endl;

                  gid.clear();
                  //
                }
            }
        }

      infile.close();
    }

    void
    read_RL_vtk(const std::string &vtkfile)
    {
      std::cout << "Begin reading vtk file" << std::endl;

      // read in vtk file here
      std::ifstream infile_mesh(vtkfile.c_str());

      // read mesh info
      mesh.read_RL_vtk(infile_mesh);

      std::ifstream infile(vtkfile.c_str());

      // read point data
      std::istringstream ss;
      std::string        str, name, type, line;
      int                numcomp;
      unsigned long int  N_points, Npoints_x, Npoints_y, Npoints_z, Npoints, u, p;

      while (!infile.eof())
        {
          std::getline(infile, line);

          if (line[0] == 'X')
            {
              if (line.size() > 12 && line.substr(0, 13) == "X_COORDINATES")
                {
                  // read header line
                  // std::cout << line << "\n";
                  ss.clear();
                  ss.str(line);
                  ss >> str >> Npoints_x >> type;

                  std::cout << "Read X_COORDINATES: " << Npoints_x << std::endl;

                  // std::cout << "  reserve OK" << std::endl;
                }
            }
          if (line[0] == 'Y')
            {
              if (line.size() > 12 && line.substr(0, 13) == "Y_COORDINATES")
                {
                  // read header line
                  // std::cout << line << "\n";
                  ss.clear();
                  ss.str(line);
                  ss >> str >> Npoints_y >> type;

                  // read points

                  std::cout << "Read Y_COORDINATES: " << Npoints_y << std::endl;

                  // std::cout << "  reserve OK" << std::endl;
                }
            }
          if (line[0] == 'Z')
            {
              if (line.size() > 12 && line.substr(0, 13) == "Z_COORDINATES")
                {
                  // read header line
                  // std::cout << line << "\n";
                  ss.clear();
                  ss.str(line);
                  ss >> str >> Npoints_z >> type;

                  // read points

                  std::cout << "Read Z_COORDINATES: " << Npoints_z << std::endl;

                  // std::cout << "  reserve OK" << std::endl;
                }
            }

          if (line[0] == 'S')
            {
              if (line.size() > 6 && line.substr(0, 7) == "SCALARS")
                {
                  ss.clear();
                  ss.str(line);
                  ss >> str >> name >> type >> numcomp;

                  // read LOOKUP_TABLE line
                  std::getline(infile, line);

                  // read data
                  std::cout << "begin reading data" << std::endl;

                  N_points = (Npoints_x) * (Npoints_y) * (Npoints_z);

                  if (DIM > 2)
                    {
                      Npoints = 8 * (Npoints_x - 1) * (Npoints_y - 1) * (Npoints_z - 1);
                    }
                  if (DIM == 2)
                    {
                      Npoints = 4 * (Npoints_x - 1) * (Npoints_y - 1);
                    }

                  std::vector<float> data(
                    N_points); // NB: data and gid are of different size
                  std::vector<double> gid(Npoints); // gid is the grainID for each node
                                                    // unsigned int gid[Npoints];

                  for (unsigned int i = 0; i < N_points; i++)
                    {
                      infile >> data[i];
                      // std::cout << data[i] << std::endl;
                    }

                  std::cout << "beginning grain_id stencil" << std::endl;

                  u = 0;
                  if (DIM > 2)
                    {
                      for (unsigned int i = 0; i < (Npoints_z - 1); i++)
                        {
                          for (unsigned int j = 0; j < (Npoints_y - 1); j++)
                            {
                              for (unsigned int k = 0; k < (Npoints_x - 1); k++)
                                {
                                  p = k + j * Npoints_x + i * Npoints_x * Npoints_y;

                                  gid[u]     = data[p];
                                  gid[u + 1] = data[p + 1];
                                  gid[u + 2] = data[p + Npoints_x];
                                  gid[u + 3] = data[p + Npoints_x + 1];
                                  gid[u + 4] = data[p + Npoints_x * Npoints_y];
                                  gid[u + 5] = data[p + Npoints_x * Npoints_y + 1];
                                  gid[u + 6] =
                                    data[p + Npoints_x + Npoints_x * Npoints_y];
                                  gid[u + 7] =
                                    data[p + Npoints_x + Npoints_x * Npoints_y + 1];

                                  u += 8;
                                }
                            }
                        }
                    }

                  if (DIM == 2)
                    {
                      for (unsigned int j = 0; j < (Npoints_y - 1); j++)
                        {
                          for (unsigned int k = 0; k < (Npoints_x - 1); k++)
                            {
                              p = k + j * Npoints_x;

                              gid[u]     = data[p];
                              gid[u + 1] = data[p + 1];
                              gid[u + 2] = data[p + Npoints_x];
                              gid[u + 3] = data[p + Npoints_x + 1];

                              u += 4;
                            }
                        }
                    }

                  data.clear();
                  // MPI_Barrier(MPI_COMM_WORLD);

                  // construct field
                  std::vector<std::string> var_name(DIM);
                  std::vector<std::string> var_description(DIM);

                  if (DIM == 2)
                    {
                      var_name[0]        = "x";
                      var_description[0] = "x coordinate";
                      var_name[1]        = "y";
                      var_description[1] = "y coordinate";
                    }
                  if (DIM > 2)
                    {
                      var_name[2]        = "z";
                      var_description[2] = "z coordinate";
                    }

                  std::cout << "Construct PField '" << name << "'" << std::endl;
                  scalar_field.push_back(PField<Coordinate, double, DIM>(name,
                                                                         var_name,
                                                                         var_description,
                                                                         mesh,
                                                                         gid,
                                                                         0.0));
                  std::cout << "  done" << std::endl;

                  gid.clear();
                  //
                }
            }
        }

      infile.close();
    }

    PField<Coordinate, double, DIM> &
    find_scalar_field(std::string name)
    {
      for (unsigned int i = 0; i < scalar_field.size(); i++)
        {
          if (scalar_field[i].name() == name)
            return scalar_field[i];
        }
      throw std::invalid_argument("Could not find scalar_field named '" + name + "'");
    }
  };
} // namespace PRISMS

#endif
