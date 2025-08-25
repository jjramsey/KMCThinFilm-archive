#ifndef SOLVER_HPP
#define SOLVER_HPP

#include "KMC_Config.hpp"

#if KMC_PARALLEL
#include <mpi.h>
#include <vector>
#endif

#include <string>

#include "EventId.hpp"
#include "RandNumGen.hpp"
#include "TimeIncrSchemeVars.hpp"

namespace KMCThinFilm {

  class Lattice;

  class Solver {
  public:

#if KMC_PARALLEL
    Solver(const Lattice * lattice)
      : lattice_(lattice),
	tIncrSchemeName_(TimeIncr::SchemeName::BAD_VALUE)
    {}
#endif

    virtual ~Solver() {};

    void setRNG(RandNumGenSharedPtr rng) {rng_ = rng;}
  
    // This does any preliminary work for adding to the event list, such
    // as clearing any previous contents of the event list.
    virtual void beginBuildingEventList(int numOverLatticeEvents,
                                        int numReservedLatticePlanes) = 0;

    // Used to add entries to the event list being built. Should ONLY be
    // called in between calls to beginBuildingEventList() and
    // endBuildingEventList().
    virtual void addCellCenteredEntryToEventList(const EventId & eId,
						 double propensity,
						 int sectNum) = 0;

    virtual void addOverLatticeEntryToEventList(const EventId & eId,
						double propensity,
						int sectNum) = 0;

    // This does any work needed after entries to event list have all
    // been added.
    virtual void endBuildingEventList() = 0;

    virtual void addOrUpdateCellCenteredEntryToEventList(const EventId & eId,
							 double propensity,
							 int sectNum) = 0;

    virtual void chooseEventIDAndUpdateTime(int sectNum,
					    EventId & chosenEventID,
					    double & time) = 0;

    virtual bool noMoreEvents(int sectNum) const = 0;

#if KMC_PARALLEL
    virtual bool noCellCenteredEvents(int sectNum) const = 0;

    void setTimeIncrScheme(const TimeIncr::SchemeVars & vars);
    TimeIncr::SchemeName::Type getTimeIncrSchemeName() const {return tIncrSchemeName_;}
    bool timeIncrSchemeIsAdaptive() const;
    void updateTStop(const MPI_Comm & comm, double & tStop);
#endif

  protected:
    RandNumGenSharedPtr rng_;
#if KMC_PARALLEL
    // To be used in calculating getLocalMaxAvgPropensityPerInPlaneCell()
    std::vector<int> numSitesPerSector_;
#endif

  private:
    
#if KMC_PARALLEL    
    const Lattice * lattice_;    
    TimeIncr::SchemeName::Type tIncrSchemeName_;

    double nStop_, tStopFixed_, tStopMax_;
    bool timeIncrSchemeIsAdaptive_;

    // Used in calculating time step schemes
    virtual double getLocalMaxAvgPropensityPerPossEvent() const = 0;
    virtual double getLocalMaxSinglePropensity() const = 0;

    typedef void (Solver::*UpdateTStop_)(const MPI_Comm & comm,
					 double & tStop);

    UpdateTStop_ updateTStop_;

    void updateTStopMaxAvgPropensityPerPossEvent_(const MPI_Comm & comm,
						  double & tStop);

    void updateTStopMaxSinglePropensity_(const MPI_Comm & comm,
					 double & tStop);

    void updateTStopFixedValue_(const MPI_Comm & comm,
				double & tStop);

    void getTstopFromLocalPropensities_(double maxLocalPropensity,
					const MPI_Comm & comm,
					double & tStop) const;

#endif

  };

}

#endif /* SOLVER_HPP */
