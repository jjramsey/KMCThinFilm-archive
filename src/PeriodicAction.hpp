#ifndef POST_GLOBAL_TIME_STEP_ACTION
#define POST_GLOBAL_TIME_STEP_ACTION

#include <boost/function.hpp>

#include "SimulationState.hpp"
#include "Lattice.hpp"

/*! \file
  \brief Defines the signature for the function object called periodically during a simulation
 */

namespace KMCThinFilm {
  /*! Signature for the function object called periodically during a simulation. */
  typedef boost::function<void (const SimulationState &, Lattice &)> PeriodicAction;
}

#endif /* POST_GLOBAL_TIME_STEP_ACTION */
