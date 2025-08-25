#ifndef SIMULATION_STATE_HPP
#define SIMULATION_STATE_HPP

// Note: This header file is documented via Doxygen
// <http://www.doxygen.org>. Comments for Doxygen begin with '/*!' or
// '//!', and descriptions of functions, class and member functions
// occur *before* their corresponding class declarations and function
// prototypes.

/*! \file
  \brief Defines the SimulationState class
 */

#include "KMC_Config.hpp"

namespace KMCThinFilm {

  /*! Data structure holding the current state of the simulation.

      This data structure provides an interface to function objects
      with the KMCThinFilm::EventExecutor and
      KMCThinFilm::PeriodicAction signatures for accessing such things
      as the current elapsed time of the simulation, the number of
      local events that have occurred, and so on.

      \see EventExecutor PeriodicAction
   */
  class SimulationState {
    friend class Simulation;
  public:
    SimulationState()
      : elapsed_time_(0), t_stop_(0), maxTime_(0), 
#if KMC_PARALLEL
	t_sector_(0),
#endif
	num_local_events_(0), num_global_steps_(0)
    {}

    /*! The amount of <EM>simulated</EM> time (not the actual wall
      clock time) that has passed since the beginning of the first
      simulation run. */
    double elapsedTime() const {
#if KMC_PARALLEL
      // This means that I'll have to make sure t_sector_ is zero
      // outside a looping over sectors.
      return elapsed_time_ + t_sector_;
#else
      return elapsed_time_;
#endif
    }

    /*! The maximum amount of simulation time (<EM>not</EM> wall-clock
      time) alloted for all the simulation runs. */
    double maxTime() const {return maxTime_;}

    /*! The current value by which the global elasped time is
        incremented in a parallel simulation.

        This returns zero in a serial simulation.
    */
    double globalTimeIncrement() const {return t_stop_;}

    /*! Number of events that have occurred on a particular
      processor. */
    unsigned long long numLocalEvents() const {return num_local_events_;}

    /*! Number of times the global clock has been incremented.

      In a serial simulation, this returns the same value as numLocalEvents().
     */
    unsigned long long numGlobalSteps() const {return num_global_steps_;}
    
  private:
    double elapsed_time_, t_stop_, maxTime_;
#if KMC_PARALLEL
    double t_sector_;
#endif
    unsigned long long num_local_events_, num_global_steps_;
  };

}

#endif /* SIMULATION_STATE_HPP */
