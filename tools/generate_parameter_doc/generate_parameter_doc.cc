// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <prismspf/user_inputs/user_input_parameters.h>

#include <fstream>

int
main(int argc, char *argv[])
{
  prismspf::UserInputParameters<3> user_inputs;
  dealii::ParameterHandler         parameter_handler;
  std::ofstream                    xml_file("parameter_doc.xml");
  std::ofstream                    json_file("parameter_doc.json");
  std::ofstream                    latex_file("parameter_doc.tex");
  std::ofstream                    prm_file("parameter_doc.prm");
  user_inputs.declare_parameters(parameter_handler, 1);
  parameter_handler.print_parameters(xml_file, dealii::ParameterHandler::XML);
  parameter_handler.print_parameters(json_file, dealii::ParameterHandler::JSON);
  parameter_handler.print_parameters(latex_file, dealii::ParameterHandler::LaTeX);
  parameter_handler.print_parameters(prm_file, dealii::ParameterHandler::PRM);

  return 0;
}
