#ifndef SIMULATION_HPP
#define SIMULATION_HPP

#include "KMC_Config.hpp"

#if KMC_PARALLEL
#include <mpi.h>
#endif

#include <boost/noncopyable.hpp>
#include <boost/scoped_ptr.hpp>

#include "CellCenteredGroupPropensities.hpp"
#include "EventExecutorGroup.hpp"
#include "PeriodicAction.hpp"
#include "CellNeighOffsets.hpp"
#include "RandNumGen.hpp"
#include "TimeIncrSchemeVars.hpp"
#include "IdsOfSolvers.hpp"

/*! \file
  \brief Defines the Simulation class.
 */

namespace KMCThinFilm {

  class Lattice;

  /*! Class for setting up and running a simulation. 

    Typical usage for this class is something like this:

    \code
    KMC_MAKE_LATTICE_INTVAL_ENUM(My,
                                 GaSiteStatus,
                                 AsSiteStatus);    

    LatticeParams latParams;

    latParams.numIntsPerCell = MyIntVal::SIZE;
    latParams.globalPlanarDims[0] = latParams.globalPlanarDims[1] = myGlobalSize(...);

    #if KMC_PARALLEL
    latParams.ghostExtent[0] = latParams.ghostExtent[1] = myGhostExtent(...);
    #endif

    Simulation sim(latParams);

    sim.setSolver(SolverId::BINARY_TREE);

    RandNumGenSharedPtr rng(new MyRandNumGen(...));    
    sim.setRNG(rng);
    
    TimeIncr::SchemeVars schemeVars;

    // Set values of schemeVars

    sim.setTimeIncrScheme(schemeVars);

    KMC_MAKE_ID_ENUM(MyOverLatticeEvents,
                     DEPOSITION_Ga,
                     DEPOSITION_As2);

    KMC_MAKE_ID_ENUM(MyCellCenteredEvents,
                     GaHOP,
                     As2Dimer_TO_AsCrystal, ... other events ... );


    sim.reserveOverLatticeEvents(MyOverLatticeEvents::SIZE);

    sim.addOverLatticeEvent(MyOverLatticeEvents::DEPOSITION_Ga,
                            myGaDepRate(...), GaDepositionExecute());

    sim.addOverLatticeEvent(MyOverLatticeEvents::DEPOSITION_As2,
                            myAsDepRate(...), AsDepositionExecute());

    sim.addCellCenteredEventGroup(MyCellCenteredEventGroups::HOPS, ...);

    // More possible events defined

    sim.run(approxDepTime);
    sim.removeOverLatticeEvent(MyOverLatticeEvents::DEPOSITION_Ga);
    sim.removeOverLatticeEvent(MyOverLatticeEvents::DEPOSITION_As2);
    sim.run(2.0*approxDepTime);

    \endcode

   */  
  class Simulation : private boost::noncopyable {
  public:

    /*! Constructor to initialize the simulation */
    explicit Simulation(const LatticeParams & paramsForLattice /*!<
                                                                  Parameters
                                                                  for
                                                                  the
                                                                  lattice
                                                                  used
                                                                  internally
                                                                  by
                                                                  the
                                                                  simulation. */);

    //! \cond HIDE_FROM_DOXYGEN
    ~Simulation();
    //! \endcond

    /*! Number of processors over which the lattice in the simulation is distributed.

      \see Lattice::nProcs()
     */
    int nProcs() const;

    /*! When the preprocessor variable KMC_PARALLEL is zero, this
      always returns zero. Otherwise, this is the MPI rank of the
      local process, where the communicator for that process is given
      by Simulation::comm().

      \see Lattice::procID()
    */
    int procID() const;


    /*! The MPI communicator used by the lattice in the simulation for
        the interchange of ghost lattice cells, if the preprocessor
        variable KMC_PARALLEL is non-zero. <STRONG>Not available in
        the serial version of the ARL KMCThinFilm library.</STRONG>

      This is <EM>not</EM> the same as the
      <VAR>latticeCommInitial</VAR> member of the LatticeParams object
      used to construct the lattice.

      \see Lattice::comm()
     */
#if KMC_PARALLEL
    const MPI_Comm & comm() const;
#endif

    /*! The number of processors along in-plane dimension <VAR>dim</VAR>. 

      \see Lattice::procPerDim()
     */
    int procPerDim(int dim) const;

    /*! Coordinates for a processor in an MPI Cartesian topology
      along in-plane dimension <VAR>dim</VAR> for a given processor. 

      \see Lattice::commCoord()
    */
    int commCoord(int dim) const;

    /*! Obtain limiting values of the local coordinates of the cells
        used in the simulation lattice. 

	\see Lattice::getLocalPlanarBBox()
    */
    void getLatticeLocalPlanarBBox(bool wGhost /*!< Indicates whether minimum
						 and maximum values for
						 lattice cell coordinates
						 include the coordinates of
						 ghost lattice cells. */,
				   LatticePlanarBBox & bbox /*!< Bounding box
							      of the local
							      lattice
							      coordinates. */) const;

    /*! Obtain limiting values of the local lattice coordinates for a
      particular sector or sublattice. 

      \see Lattice::getSectorPlanarBBox()
    */
    void getLatticeSectorPlanarBBox(int sectNum /*!< Sector number, ranging from 0 to numSectors() - 1 */,
				    LatticePlanarBBox & bbox /*!<
							       Bounding
							       box of
							       the
							       lattice
							       coordinates
							       within
							       the
							       sector. */) const;

    /*! Obtain limiting values of the global lattice coordinates. 

      \see Lattice::getGlobalPlanarBBox()
     */
    void getLatticeGlobalPlanarBBox(LatticePlanarBBox & bbox /*!< Bounding
							       box of the
							       global lattice
							       coordinates. */) const;

    /*! The amount of <EM>simulated</EM> time (not the actual wall
      clock time) that has passed since the beginning of the first
      simulation run. 

      \see SimulationState::elapsedTime()
    */
    double elapsedTime() const;

    /*! Number of events that have occurred on a particular
      processor.

      \see SimulationState::numLocalEvents()
    */
    unsigned long long numLocalEvents() const;

    /*! Number of times the global clock has been incremented. For
        serial simulations, this yields the same value as
        numLocalEvents(). 

	\see SimulationState::numGlobalSteps()
    */
    unsigned long long numGlobalSteps() const;

    /*! Sets the number of possible cell-centered event groups for the
        simulation to <VAR>numGroups</VAR>, and the number of
        individual possible cell-centered events to
        <VAR>numTotEvents</VAR>.

      Should be called before a call to addCellCenteredEventGroup(), and
      the number of calls to addCellCenteredEventGroup() should be the same
      as <VAR>numGroups</VAR>. */
    void reserveCellCenteredEventGroups(int numGroups, int numTotEvents);

    /*! Adds a possible cell-centered event group to the simulation. */
    void addCellCenteredEventGroup(int eventGroupId /*!< Unique integer ID of event group */,
                                   const CellNeighOffsets & cno /*!<
                                                                  Offsets
                                                                  used
                                                                  to
                                                                  find
                                                                  cells
                                                                  used
                                                                  to
                                                                  determine
                                                                  the
                                                                  event
                                                                  propensities */,
                                   CellCenteredGroupPropensities propensities /*!<
                                                                            Function
                                                                            or
                                                                            function
                                                                            object
                                                                            used to
                                                                            determine
                                                                            the propensities
                                                                            of this group of events.
                                                                            */,
                                   const EventExecutorGroup & eventExecutorGroup /*!< Group of
                                                                                   function
                                                                                   objects,
                                                                                   one of which
                                                                                   runs if
                                                                                   this
                                                                                   event is
                                                                                   executed. */);

    /*! Changes a possible cell-centered event group that has been
        previously added to the simulation.

	\see addCellCenteredEventGroup()
    */
    void changeCellCenteredEventGroup(int eventGroupId /*!< Unique integer ID of event group*/,
                                      const CellNeighOffsets & cno /*!<
                                                                  Offsets
                                                                  used
                                                                  to
                                                                  find
                                                                  cells
                                                                  used
                                                                  to
                                                                  determine
                                                                  the
                                                                  event
                                                                  propensities */,
                                      CellCenteredGroupPropensities propensities /*!<
                                                                            Function
                                                                            or
                                                                            function
                                                                            object
                                                                            used to
                                                                            determine
                                                                            the propensities
                                                                            of this group of events.
                                                                            */,
                                      const EventExecutorGroup & eventExecutorGroup /*!< Group of
                                                                                   function
                                                                                   objects,
                                                                                   one of which
                                                                                   runs if
                                                                                   this
                                                                                   event is
                                                                                   executed. */);

    /*! Removes a previously added possible cell-centered event group from the simulation. 

	\see addCellCenteredEventGroup() changeCellCenteredEventGroup()
    */
    void removeCellCenteredEventGroup(int eventGroupId /*!< Integer ID of event group to be removed. */);

    /*! Sets the number of "over-lattice" possible events for the
        simulation, that is, events that occur at a random site over
        the lattice (e.g. deposition of an atom), to <VAR>num</VAR>.

	Should be called before a call to addOverLatticeEvent(), and
	the number of calls to addOverLatticeEvent() should be the same
	as <VAR>num</VAR>. */
    void reserveOverLatticeEvents(int num);

    /*! Adds a possible "over-lattice" event to the simulation, that
        is, one that may occur at a random site over the lattice
        (e.g. deposition of an atom). Here, the cells affected by the
        event are determined by "auto-tracking."

	If this event occurs, the CellInds argument of
	<VAR>eventExecutor</VAR> will have its first two components, i
	& j, set to a random in-plane position within the current
	sector, and its third component, k, set to one less than the
	height of the lattice.

        \see EventExecutorAutoTrack
     */
    void addOverLatticeEvent(int eventId /*!< Unique integer ID of event. */,
			     double propensityPerUnitArea /*!< Propensity per unit area (e.g. deposition flux) of event */,
			     EventExecutorAutoTrack eventExecutor /*!< Function
                                                                    or
                                                                    function
                                                                    object
                                                                    that
                                                                    runs if
                                                                    this
                                                                    event is
                                                                    executed. */);

    /*! Adds a possible "over-lattice" event to the simulation, that
        is, one that may occur at a random site over the lattice
        (e.g. deposition of an atom). Here, the cells affected by the
        event are determined by "semi-manual" tracking.

	If this event occurs, the CellInds argument of
	<VAR>eventExecutor</VAR> will have its first two components, i
	& j, set to a random in-plane position within the current
	sector, and its third component, k, set to one less than the
	height of the lattice.

        \see EventExecutorSemiManualTrack
     */
    void addOverLatticeEvent(int eventId /*!< Unique integer ID of event. */,
			     double propensityPerUnitArea /*!< Propensity per unit area (e.g. deposition flux) of event */,
			     EventExecutorSemiManualTrack eventExecutor /*!< Function
                                                                          or
                                                                          function
                                                                          object
                                                                          that
                                                                          runs if
                                                                          this
                                                                          event is
                                                                          executed. */,
                             const std::vector<CellNeighOffsets> & cnoVec /*!<
                                                                            Relative positions of cells 
                                                                            directly changed by event. */);

    /*! Changes a possible "over-lattice" event that has been
        previously added to the simulation.  Here, the cells affected
        by the event are determined by "auto-tracking."

	\see addOverLatticeEvent() EventExecutorAutoTrack
    */
    void changeOverLatticeEvent(int eventId /*!< Integer ID of event to be changes. */,
				double propensityPerUnitArea /*!< Propensity per unit area (e.g. deposition flux) of event */,
				EventExecutorAutoTrack eventExecutor /*!< Function
                                                              or
                                                              function
                                                              object
                                                              that
                                                              runs if
                                                              this
                                                              event is
                                                              executed. */);

    /*! Changes a possible "over-lattice" event that has been
        previously added to the simulation. Here, the cells affected
        by the event are determined by "semi-manual" tracking.

	\see addOverLatticeEvent() EventExecutorSemiManualTrack
    */
    void changeOverLatticeEvent(int eventId /*!< Integer ID of event to be changes. */,
				double propensityPerUnitArea /*!< Propensity per unit area (e.g. deposition flux) of event */,
				EventExecutorSemiManualTrack eventExecutor /*!< Function
                                                                          or
                                                                          function
                                                                          object
                                                                          that
                                                                          runs if
                                                                          this
                                                                          event is
                                                                          executed. */,
                                const std::vector<CellNeighOffsets> & cnoVec /*!<
                                                                               Relative positions of cells 
                                                                               directly changed by event. */);

    /*! Removes a previously added possible "over-lattice" from the simulation. 

	\see addOverLatticeEvent() changeOverLatticeEvent()
    */
    void removeOverLatticeEvent(int eventId /*!< Integer ID of event to be removed. */);

    /*! [<STRONG>ADVANCED</STRONG>] If <VAR>doTrack</VAR> is true,
        store indices of the lattice cells changed by the periodic
        actions that occur.

       If a periodic action actually <EM>changes</EM> the lattice
       (i.e. calls Lattice::setInt(), Lattice::setFloat(), or
       Lattice::addPlanes()), then by default, the event list is
       rebuilt, which is an O(N) operation, with N being the number of
       lattice cells. This function changes this behavior. Since the
       indices of changed lattice cells are actually stored, only the
       entries in the event list affected by the changed lattice cells
       need to be changed. Note, though, that if the periodic actions
       change a large number of lattice sites, then the default
       behavior is probably more optimal, or at least less
       memory-hungry.
     */
    void trackCellsChangedByPeriodicActions(bool doTrack);

     /*! Sets the number of time-periodic actions for the
       simulation to <VAR>num</VAR>.

      Should be called before a call to addTimePeriodicAction(), and
      the number of calls to addTimePeriodicAction() should be the
      same as <VAR>num</VAR>. */
    void reserveTimePeriodicActions(int num);

    /*! Adds a time-periodic action to the simulation. 

      Time-periodic actions are actions that occur every P units of
      simulation time, with P being the time period.
     */
    void addTimePeriodicAction(int actionId /*!< Unique integer ID for the action */,
			       PeriodicAction action /*!< Function or
                                                        function
                                                        object
                                                        executed by
                                                        this
                                                        action. */,
			       double periodOfTime /*!< The period of the action */,
			       bool doAtSimEnd /*!< If true, this
                                                  action is done at
                                                  the end of a
                                                  simulation run as
                                                  well. */);
    
    /*! Changes a time-periodic action that had been previously added
        to the simulation.

	\see addTimePeriodicAction()
    */
    void changeTimePeriodicAction(int actionId /*!< Integer ID of the action to be changed */,
				  PeriodicAction action /*!< Function or
                                                        function
                                                        object
                                                        executed by
                                                        this
                                                        action. */,
				  double periodOfTime /*!< The period of the action */,
				  bool doAtSimEnd /*!< If true, this
                                                  action is done at
                                                  the end of a
                                                  simulation run as
                                                  well. */);

    /*! Removes a previously-added time-periodic action from the simulation.

	\see addTimePeriodicAction() changeTimePeriodicAction()
    */
    void removeTimePeriodicAction(int actionId /*!< Integer ID of the action to be removed */);

    /*! Sets the number of step-periodic actions for the
       simulation to <VAR>num</VAR>.

      Should be called before a call to addStepPeriodicAction(), and
      the number of calls to addStepPeriodicAction() should be the
      same as <VAR>num</VAR>. */
    void reserveStepPeriodicActions(int num);

    /*! Adds a step-periodic action to the simulation. 

      Step-periodic actions are actions that occur every P time steps of the
      simulation, with P being the period.
     */
    void addStepPeriodicAction(int actionId /*!< Unique integer ID for the action */,
			       PeriodicAction action /*!< Function or
                                                        function
                                                        object
                                                        executed by
                                                        this
                                                        action. */,
			       int periodOfSteps /*!< The period of the action */,
			       bool doAtSimEnd /*!< If true, this
                                                  action is done at
                                                  the end of a
                                                  simulation run as
                                                  well. */);

    /*! Changes a step-periodic action that had been previously added
        to the simulation.

	\see addStepPeriodicAction()
    */
    void changeStepPeriodicAction(int actionId /*!< Integer ID of the action to be changed */,
				  PeriodicAction action /*!< Function or
                                                        function
                                                        object
                                                        executed by
                                                        this
                                                        action. */,
				  int periodOfSteps /*!< The period of the action */,
				  bool doAtSimEnd /*!< If true, this
                                                  action is done at
                                                  the end of a
                                                  simulation run as
                                                  well. */);

    /*! Removes a previously-added step-periodic action from the simulation.

	\see addStepPeriodicAction() changeStepPeriodicAction()
    */
    void removeStepPeriodicAction(int actionId);

    /*! Sets the type of solver, i.e. the means of storing and choosing events, for the simulation. */
    void setSolver(SolverId::Type sId);

    /*! Sets the parallel time-incrementing scheme.

      If KMC_PARALLEL equals zero, this does nothing. This must be
      called <STRONG>after</STRONG> setSolver() has been called. */
    void setTimeIncrScheme(const TimeIncr::SchemeVars & vars);

    /*! Sets the random-number generator of the simulation.

      This must be called <STRONG>after</STRONG> setSolver() has been called. */
    void setRNG(RandNumGenSharedPtr rng);

    /*! Runs the simulation */
    void run(double runTime /*!< Length of simulation time */);
  private:
    class Impl_;
    boost::scoped_ptr<Impl_> pImpl_;
  };

}

#endif /* SIMULATION_HPP */
