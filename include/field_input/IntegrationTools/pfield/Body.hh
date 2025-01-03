#ifndef Body_HH
#define Body_HH

#include <field_input/IntegrationTools/pfield/Mesh.hh>
#include <field_input/IntegrationTools/pfield/PField.hh>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

namespace PRISMS
{

  /**
   * \brief  A class for a Body-a combination of Mesh and Field(s)-used in file read-in.
   */
  template <class Coordinate, int DIM>
  class Body
  {
  public:
    Mesh<Coordinate, DIM> mesh;

    std::vector<PField<Coordinate, double, DIM>> scalar_field;

    /**
     * \brief Constructor.
     */
    Body() = default;

    /**
     * \brief Read from a 2D/3D vtk file with quad/hex elements
     */
    void
    read_vtk(const std::string &vtkfile)
    {
      // Read vtk file
      std::ifstream infile(vtkfile.c_str());
      if (!infile.is_open())
        {
          throw std::runtime_error("Could not open VTK file: " + vtkfile);
        }

      // Read mesh info
      mesh.read_vtk(infile);

      // Read point data
      std::istringstream ss;
      std::string        str;
      std::string        name;
      std::string        type;
      std::string        line;
      int                numcomp = 0;
      unsigned long int  Npoints = 0;

      while (std::getline(infile, line))
        {
          // Check for Point data
          if (line.rfind("POINT_DATA", 0) == 0)
            {
              ss.clear();
              ss.str(line);
              ss >> str >> Npoints;
            }

          // Check for Scalar data
          if (line.rfind("SCALARS", 0) == 0)
            {
              ss.clear();
              ss.str(line);
              ss >> str >> name >> type >> numcomp;

              // Skip LOOKUP_TABLE line that follows SCALAR line
              std::getline(infile, line);

              // Read data
              std::vector<double> gid(Npoints);
              for (auto &value : gid)
                {
                  infile >> value;
                }

              // Check dim
              switch (DIM)
                {
                  case 1:
                    throw std::runtime_error(
                      "1D dimensional vtk read-in is not currently supported.");
                    break;
                  case 2:
                    break;
                  case 3:
                    break;
                  default:
                    throw std::runtime_error("Invalid dimension for vtk read-in.");
                }

              // Construct field
              scalar_field.emplace_back(name, mesh, gid, 0.0);
            }
        }
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

                  std::cout << "Construct PField '" << name << "'" << std::endl;
                  scalar_field.push_back(
                    PField<Coordinate, double, DIM>(name, mesh, gid, 0.0));
                  std::cout << "  done" << std::endl;

                  gid.clear();
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
            {
              return scalar_field[i];
            }
        }
      throw std::invalid_argument("Could not find scalar_field named '" + name + "'");
    }
  };

  template class Body<double *, 2>;
  template class Body<double *, 3>;

} // namespace PRISMS

#endif
