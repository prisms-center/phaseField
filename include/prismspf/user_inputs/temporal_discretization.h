#ifndef temporal_discretization_h
#define temporal_discretization_h

#include <prismspf/config.h>
#include <prismspf/core/conditional_ostreams.h>

PRISMS_PF_BEGIN_NAMESPACE

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
  unsigned int total_increments = 0;

  // The current increment
  unsigned int increment = 0;

  // The current time
  double time = 0.0;
};

inline void
temporalDiscretization::print_parameter_summary() const
{
  prisms::conditionalOStreams::pout_summary()
    << "================================================\n"
    << "  Temporal Discretization\n"
    << "================================================\n"
    << "Timestep: " << dt << "\n"
    << "Total increments: " << total_increments << "\n"
    << "Final time: " << final_time << "\n\n"
    << std::flush;
}

PRISMS_PF_END_NAMESPACE

#endif