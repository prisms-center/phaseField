#ifndef temporal_discretization_h
#define temporal_discretization_h

#include <core/conditional_ostreams.h>

/**
 * \brief Class that holds temporal discretization parameters.
 */
class temporalDiscretization
{
public:
  /**
   * \brief Constructor.
   */
  temporalDiscretization() = default;

  /**
   * \brief Destructor.
   */
  ~temporalDiscretization() = default;

  /**
   * \brief Print parameters to summary.log
   */
  void
  print_parameter_summary() const;

  // Timestep
  double dt = 0.0;

  // Final time
  double final_time = 0.0;

  // Total number of increments
  uint total_increments = 0;
};

inline void
temporalDiscretization::print_parameter_summary() const
{
  conditionalOStreams::pout_summary()
    << "================================================\n"
    << "\tTemporal Discretization\n"
    << "================================================\n"
    << "Timestep: " << dt << "\n"
    << "Total increments: " << total_increments << "\n"
    << "Final time: " << final_time << "\n\n"
    << std::flush;
}

#endif