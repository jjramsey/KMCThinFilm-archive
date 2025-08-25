#ifndef EVENT_EXECUTOR_HPP
#define EVENT_EXECUTOR_HPP

#include <vector>

#include <boost/function.hpp>

#include "CellInds.hpp"
#include "CellsToChange.hpp"
#include "Lattice.hpp"
#include "SimulationState.hpp"

/*! \file
  \brief Defines the possible signatures for a function object called when an event executes.
 */

namespace KMCThinFilm {

  /*! Signature for a function object to execute an event where "auto-tracking" is used.

    In auto-tracking, cells affected by the event are determined
    automatically from the offsets used to calculate propensities and
    the positions of cells directly changed by the executed event,
    which are logged automatically when Lattice::setInt() and
    Lattice::setFloat() are called. This mode of tracking may lead to
    redundant propensity calculations. */
  typedef boost::function<void (const CellInds &, const SimulationState &, Lattice &)> EventExecutorAutoTrack;

  /*! Signature for a function object to execute an event where "semi-manual" tracking is used. 

   In semi-manual tracking, the user specifies the relative positions
   of the cells changed by an event. Using this information along with
   the offsets used to calculate propensities can reduce or eliminate
   redundant propensity calculations.*/
  typedef boost::function<void (const CellInds &,
                                const SimulationState &,
                                const Lattice &,
                                std::vector<CellsToChange> &)> EventExecutorSemiManualTrack;
}

#endif /* EVENT_EXECUTOR_HPP */
