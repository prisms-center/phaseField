#ifndef INCLUDE_NUCLEATIONPARAMETERS_H_
#define INCLUDE_NUCLEATIONPARAMETERS_H_

template<int dim>
class nucleationParameters
{
public:
    unsigned int var_index;
    std::vector<double> semiaxes;
	std::vector<double> freeze_semiaxes;
	double no_nucleation_border_thickness;
	double hold_time;
	//double min_spacing; // Only enforced for nuclei placed during the same time step
	//double order_parameter_cutoff;
	//unsigned int steps_between_attempts;
};

#endif
