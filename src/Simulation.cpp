#include "Simulation.hpp"
#include "SimulationState.hpp"
#include "ErrorHandling.hpp"
#include "CallMemberFunction.hpp"
#include "EventId.hpp"
#include "Lattice.hpp"
#include "SolverFactory.hpp"

#include <vector>
#include <deque>
#include <string>
#include <iostream>
#include <set>
#include <map>
#include <cmath>
#include <algorithm>

#include <boost/array.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/variant.hpp>
#include <boost/scoped_ptr.hpp>

#ifdef KMC_AVOID_BOOST_BIMAP
#include "MultiIndexBimap.hpp"
#else
#include <boost/bimap.hpp>
#include <boost/bimap/unordered_set_of.hpp>
#endif

using namespace KMCThinFilm;

struct Simulation::Impl_ {
  Impl_(const LatticeParams & paramsForLattice);

  Lattice lattice_;
  RandNumGenSharedPtr rng_;
  bool tIncrSchemeIsSet_;
  bool rngIsSet_;
  bool runHasBeenExecuted_;

  CellNeighProbe cellNeighProbe_;
  boost::scoped_ptr<Solver> solver_;

#ifdef KMC_AVOID_BOOST_BIMAP
  typedef MultiIndexBimap<int,std::size_t> IdIndexBimap_;
#else
  typedef boost::bimap<boost::bimaps::unordered_set_of<int>,
		       boost::bimaps::unordered_set_of<std::size_t>
		       > IdIndexBimap_;
#endif
  
  IdIndexBimap_ cellCenGroupPropensitiesIdIndexBimap_;
  IdIndexBimap_ overLatticeEventIdIndexBimap_;
  IdIndexBimap_ timePeriodicActionIdIndexBimap_;
  IdIndexBimap_ stepPeriodicActionIdIndexBimap_;

  std::vector<LatticePlanarBBox> sectorPlanarBBox_;
  LatticePlanarBBox localPlanarBBox_;

  SimulationState simState_;

  struct CellCenteredGroupPropensities_ {
    std::vector<CellIndsOffset> cioVec_;
    std::vector<std::size_t> eventVecInds_;
    CellCenteredGroupPropensities propensities_;

    CellCenteredGroupPropensities_(const CellNeighOffsets & cno,
                                   CellCenteredGroupPropensities propensities);
  };

  struct EventExecutorSemiManualTrackInfo_ {
    EventExecutorSemiManualTrack evExec_;
    std::vector<CellsToChange> ctcVec_;
    std::vector<std::vector<CellIndsOffset> > affectedCellOffsets_;

    EventExecutorSemiManualTrackInfo_(EventExecutorSemiManualTrack evExec,
                                      Lattice * lattice, 
                                      const std::vector<CellNeighOffsets> & cnoVec);

  };
  
  typedef boost::variant<EventExecutorAutoTrack,
                         EventExecutorSemiManualTrackInfo_> EventExecutor_;

  class SetEventExecutorAffectedCellOffsets_ :  public boost::static_visitor<> {
  public:

    SetEventExecutorAffectedCellOffsets_(const std::vector<CellIndsOffset> * allReversedOffsets)
      : allReversedOffsets_(allReversedOffsets)
    {}
    
    void operator()(EventExecutorAutoTrack & evExec) { /* Do nothing */ }
    void operator()(EventExecutorSemiManualTrackInfo_ & evExecInfo);

  private:
    const std::vector<CellIndsOffset> * allReversedOffsets_;
  };
                         
  struct ClearEventExecutor_ :  public boost::static_visitor<> {
    
    void operator()(EventExecutorAutoTrack & evExec) {
      evExec.clear();
    }

    void operator()(EventExecutorSemiManualTrackInfo_ & evExecInfo) {
      evExecInfo.evExec_.clear();
      evExecInfo.ctcVec_.clear();
    }

  } clearEventExecutor_;

  class RunEventExecutor_ :  public boost::static_visitor<> {
  public:
    RunEventExecutor_(Impl_ * impl) : impl_(impl) {}
    
    void operator()(EventExecutorAutoTrack & evExec);
    void operator()(EventExecutorSemiManualTrackInfo_ & evExecInfo);

    CellInds ci;
  private:
    Impl_ * impl_;
  } runEventExecutor_;

  struct OverLatticeEvent_ {
    std::vector<double> propensity_;
    EventExecutor_ execute_;
    
    OverLatticeEvent_(double propensityPerUnitArea,
		      const std::vector<LatticePlanarBBox> & sectorPlanarBBox,
		      EventExecutor_ eventExecutor);
  };

  class TimePeriodicAction_ {
  public:
    TimePeriodicAction_(PeriodicAction action,
			double period,
			bool doAtSimEnd)
      : period_(period),
	timeForNextAction_(period),
	doAtSimEnd_(doAtSimEnd),
	action_(action)
    {}

    void doActionIfEndOfPeriod(const SimulationState & simState,
			       Lattice & lattice) {

      if ((simState.elapsed_time_ >= timeForNextAction_) ||
	  (doAtSimEnd_ && (simState.elapsedTime() >= simState.maxTime()))) {

	timeForNextAction_ += period_;

	action_(simState, lattice);
      }

    }

  private:
    double period_, timeForNextAction_;
    bool doAtSimEnd_;
    PeriodicAction action_;
  };

  class StepPeriodicAction_ {
  public:
    StepPeriodicAction_(PeriodicAction action,
			unsigned long long period,
			bool doAtSimEnd)
      : period_(period),
	nStepsForNextAction_(period),
	doAtSimEnd_(doAtSimEnd),
	action_(action)
    {}

    void doActionIfEndOfPeriod(const SimulationState & simState,
			       Lattice & lattice) {

      if ((simState.numGlobalSteps() >= nStepsForNextAction_) ||
	  (doAtSimEnd_ && (simState.elapsedTime() >= simState.maxTime()))) {

	nStepsForNextAction_ += period_;

	action_(simState, lattice);
      }

    }

  private:
    unsigned long long period_, nStepsForNextAction_;
    bool doAtSimEnd_;
    PeriodicAction action_;
  };

  Lattice::TrackType::Type changeTrackingForPeriodicAction_;

  std::vector<CellCenteredGroupPropensities_> cellCenGroupPropensitiesVec_;
  std::vector<EventExecutor_> cellCenEventVec_;
  std::vector<OverLatticeEvent_> overLatticeEventVec_;
  std::vector<TimePeriodicAction_> timePeriodicActionVec_;
  std::vector<StepPeriodicAction_> stepPeriodicActionVec_;

  std::deque<std::size_t> freeCellCenEventVecInds_;

  std::vector<CellIndsOffset> reversedOffsetsVec_;

  std::size_t bimapIdToIndex_(int id,
			      const IdIndexBimap_ & idIndexBimap,
			      const std::string & callingFunc,
			      const std::string & idType) const;

  template <typename T>
  void removeEntryFromVecWithBimap_(int id,
				    IdIndexBimap_ & idIndexBimap,
				    std::vector<T> & myVec,
				    const std::string & callingFunc,
				    const std::string & idType) {

    std::size_t indexToBeOverwritten = bimapIdToIndex_(id,
						       idIndexBimap,
						       callingFunc,
						       idType);

    std::size_t lastIndex = myVec.size()  - 1;

    if (indexToBeOverwritten != lastIndex) {
      myVec.at(indexToBeOverwritten) = myVec.back();
    
      int idToBeRedirected = idIndexBimap.right.at(lastIndex);

      idIndexBimap.left.erase(id);
      idIndexBimap.left.erase(idToBeRedirected);
      idIndexBimap.insert(IdIndexBimap_::value_type(idToBeRedirected, indexToBeOverwritten));
    }
    else {
      idIndexBimap.left.erase(id);
    }

    myVec.pop_back();
  }

  void removeCellCenteredEventGroup_(int eventGroupId, const std::string & callingFunc);

  void mkReversedOffsets_();

  void updateEventAndAddrMapsFromChangedCellInds_(const Lattice::ChangedCellInds & ccInds,
						  const Lattice::OtherCheckedCellInds & ccIndsOther);

  void updateEventAndAddrMapsFromAffectedCellOffsets_(const CellsToChange & ctcVec,
                                                      const std::vector<CellIndsOffset> & affectedCellOffsets,
                                                      const Lattice::OtherCheckedCellInds & ccIndsOther);

  void doPreRunChecks_() const;

  void executeEvent_(const EventId & chosenEventId);

  typedef void (Impl_::*UpdateEventAndAddrMapsAfterPeriodicActions_)();

  UpdateEventAndAddrMapsAfterPeriodicActions_ updateEventAndAddrMapsAfterPeriodicActions_;

  void updateEventAndAddrMapsAfterPeriodicActionsWTrack_();
  void updateEventAndAddrMapsAfterPeriodicActionsNoTrack_();

  void runPeriodicActions_();

  void rebuildEventAndAddrMaps_();

  void run_(double maxTime);

  /*typedef void (Solver::*AddOrUpdateSolverFunc_)(const EventId & eId,
                                           double propensity,
                                           int sectNum);*/

  std::vector<double> tmpPropensitiesVec_; // Temporary "workspace"
                                           // vector for calculating
                                           // propensities.

  template<typename T>
  void doForCellCenteredGroupPropensities_(const CellInds & ci,
                                           int sectNum,
                                           const T & solverFunc) {
  
    for (std::vector<CellCenteredGroupPropensities_>::const_iterator ccGPropItr = cellCenGroupPropensitiesVec_.begin(),
           ccGPropItrEnd = cellCenGroupPropensitiesVec_.end(); ccGPropItr != ccGPropItrEnd; ++ccGPropItr) {
      const std::vector<std::size_t> & eventVecInds = ccGPropItr->eventVecInds_;

      std::size_t eventVecIndsSize = eventVecInds.size();

      // This ensures that all entries in tmpPropensitiesVec_ are initialized to zero.
      tmpPropensitiesVec_.clear();
      tmpPropensitiesVec_.resize(eventVecIndsSize, 0.0);
            
      cellNeighProbe_.attachCellInds(&ci, &(ccGPropItr->cioVec_));
      ccGPropItr->propensities_(cellNeighProbe_, tmpPropensitiesVec_);

      for (std::size_t i = 0; i < eventVecIndsSize; ++i) {
        solverFunc(*solver_, ci, eventVecInds[i], tmpPropensitiesVec_[i], sectNum);
      }
    }

  }

  // Function objects for use with doForCellCenteredGroupPropensities_:

  struct AddOrUpdateCellCenteredEntryToEventList_ {
    void operator()(Solver & solver,
                    const CellInds & ci,
                    std::size_t eventVecInd,
                    double propensity,
                    int sectNum) const {
      solver.addOrUpdateCellCenteredEntryToEventList(EventId(ci, eventVecInd),
                                                     propensity, sectNum);
    }
  } addOrUpdateCellCenteredEntryToEventList_;

  struct AddCellCenteredEntryToEventList_ {
    void operator()(Solver & solver,
                    const CellInds & ci,
                    std::size_t eventVecInd,
                    double propensity,
                    int sectNum) const {

      if (propensity > 0) {
        solver.addCellCenteredEntryToEventList(EventId(ci, eventVecInd),
                                               propensity, sectNum);
      }

    }
  } addCellCenteredEntryToEventList_;

  void updateEventAndAddrMapsFromChangedCell_(int currHeight, const CellInds & ci);

  void updateEventAndAddrMapsFromAffectedCell_(int currHeight,
                                               CellInds & ci /* This
                                                                should *not* 
                                                                be a const reference */);

#if KMC_PARALLEL  
  void updateEventAndAddrMapsAffectedByGhostUpdates_(const std::vector<std::vector<IJK> > & receivedInds);
#endif

};

Simulation::Impl_::CellCenteredGroupPropensities_::CellCenteredGroupPropensities_(const CellNeighOffsets & cno,
                                                                                  CellCenteredGroupPropensities propensities) 
  : cioVec_(cno.numOffsets()),
    propensities_(propensities) {

  for (int i = 0; i < cno.numOffsets(); ++i) {
    cioVec_[i] = cno.getOffset(i);
  }
}

Simulation::Impl_::OverLatticeEvent_::OverLatticeEvent_(double propensityPerUnitArea,
							const std::vector<LatticePlanarBBox> & sectorPlanarBBox,
							EventExecutor_ eventExecutor)
  : execute_(eventExecutor) {
  
  propensity_.reserve(sectorPlanarBBox.size());

  for (std::vector<LatticePlanarBBox>::const_iterator itr = sectorPlanarBBox.begin(),
	 itrEnd = sectorPlanarBBox.end(); itr != itrEnd; ++itr) {
    propensity_.push_back(propensityPerUnitArea*(itr->imaxP1 - itr->imin)*(itr->jmaxP1 - itr->jmin));
  }

}

Simulation::Impl_::EventExecutorSemiManualTrackInfo_::EventExecutorSemiManualTrackInfo_(EventExecutorSemiManualTrack evExec,
                                                                                        Lattice * lattice, 
                                                                                        const std::vector<CellNeighOffsets> & cnoVec) 
  : evExec_(evExec) {

  ctcVec_.reserve(cnoVec.size());

  for (std::vector<CellNeighOffsets>::const_iterator itr = cnoVec.begin(),
         itrEnd = cnoVec.end(); itr != itrEnd; ++itr) {

    int numOffsets = itr->numOffsets();

    ctcVec_.push_back(CellsToChange(lattice, numOffsets));

    for (int i = 0; i < numOffsets; ++i) {
      ctcVec_.back().addOffset(itr->getOffset(i));
    }
    
  }

}

void Simulation::Impl_::SetEventExecutorAffectedCellOffsets_::operator()(EventExecutorSemiManualTrackInfo_ & evExecInfo) {
  
  evExecInfo.affectedCellOffsets_.clear();
  evExecInfo.affectedCellOffsets_.reserve(evExecInfo.ctcVec_.size());
  
  for (std::vector<CellsToChange>::const_iterator ctcItr = evExecInfo.ctcVec_.begin(),
         ctcItrEnd = evExecInfo.ctcVec_.end(); ctcItr != ctcItrEnd; ++ctcItr) {

    const std::vector<CellIndsOffset> & cioVec = ctcItr->getCellIndsOffsetVec();

    // Using an *ordered* set because an unordered set affects the reproducibility of simulations.
    std::set<CellIndsOffset> affectedCellOffsets;

    for (std::vector<CellIndsOffset>::const_iterator rvOffsetItr = allReversedOffsets_->begin(),
           rvOffsetItrEnd = allReversedOffsets_->end();
         rvOffsetItr != rvOffsetItrEnd; ++rvOffsetItr) {
      
      for (std::vector<CellIndsOffset>::const_iterator cioItr = cioVec.begin(),
             cioItrEnd = cioVec.end(); cioItr != cioItrEnd; ++cioItr) {
        affectedCellOffsets.insert(*rvOffsetItr + *cioItr);
      }

    }

    evExecInfo.affectedCellOffsets_.push_back(std::vector<CellIndsOffset>());

    std::vector<CellIndsOffset> & affectedCellOffsetsVec = evExecInfo.affectedCellOffsets_.back();

    affectedCellOffsetsVec.reserve(affectedCellOffsets.size());

    for (std::set<CellIndsOffset>::const_iterator asItr = affectedCellOffsets.begin(),
           asItrEnd = affectedCellOffsets.end(); asItr != asItrEnd; ++asItr) {
      affectedCellOffsetsVec.push_back(*asItr);
    }

  }

}

void Simulation::Impl_::RunEventExecutor_::operator()(EventExecutorAutoTrack & evExec) {  

  impl_->lattice_.trackChanges(Lattice::TrackType::RECORD_CHANGED_CELL_INDS);
  evExec(ci, impl_->simState_, impl_->lattice_);
  impl_->updateEventAndAddrMapsFromChangedCellInds_(impl_->lattice_.getChangedCellInds(),
                                                    impl_->lattice_.getOtherCheckedCellInds());
  impl_->lattice_.trackChanges(Lattice::TrackType::NONE);

}

void Simulation::Impl_::RunEventExecutor_::operator()(EventExecutorSemiManualTrackInfo_ & evExecInfo) {
  impl_->lattice_.trackChanges(Lattice::TrackType::RECORD_ONLY_OTHER_CHANGED_CELL_INDS);
  evExecInfo.evExec_(ci, impl_->simState_, impl_->lattice_, evExecInfo.ctcVec_);

  std::size_t numCenters = evExecInfo.ctcVec_.size();

  for (std::size_t i = 0; i < numCenters; ++i) {
    impl_->updateEventAndAddrMapsFromAffectedCellOffsets_(evExecInfo.ctcVec_[i],
                                                          evExecInfo.affectedCellOffsets_[i],
                                                          impl_->lattice_.getOtherCheckedCellInds());
  }

  impl_->lattice_.trackChanges(Lattice::TrackType::NONE);
}

Simulation::Impl_::Impl_(const LatticeParams & paramsForLattice)
  : lattice_(paramsForLattice),
    tIncrSchemeIsSet_(false),
    rngIsSet_(false),
    runHasBeenExecuted_(false),
    cellNeighProbe_(&lattice_),
    runEventExecutor_(this) {

  sectorPlanarBBox_.resize(lattice_.numSectors());

  lattice_.getLocalPlanarBBox(false, localPlanarBBox_);

  for (int i = 0; i < lattice_.numSectors(); ++i) {
    lattice_.getSectorPlanarBBox(i, sectorPlanarBBox_[i]);
  }

}

std::size_t Simulation::Impl_::bimapIdToIndex_(int id,
					       const IdIndexBimap_ & idIndexBimap,
					       const std::string & callingFunc,
					       const std::string & idType) const {
  std::size_t index = 0;

  try {
    // Not returning within the "try" block in order to mollify the Intel compiler.
    index = idIndexBimap.left.at(id);
  }
  catch (const std::out_of_range & e) {
    exitWithMsg(std::string(e.what()) + "; " + callingFunc + " error: " + idType + " " + 
		boost::lexical_cast<std::string>(id) + " does not exist.");
  }

  return index;
}

void Simulation::Impl_::removeCellCenteredEventGroup_(int eventGroupId, const std::string & callingFunc) {

  std::size_t indexToBeOverwritten = bimapIdToIndex_(eventGroupId,
                                                     cellCenGroupPropensitiesIdIndexBimap_,
                                                     callingFunc,
                                                     "eventGroupId");

  std::size_t lastIndex = cellCenGroupPropensitiesVec_.size()  - 1;

  // Mostly a cut-and-paste from removeEntryFromVecWithBimap_, which
  // is not such a wonderful thing.

  {
    CellCenteredGroupPropensities_ & ccgp = cellCenGroupPropensitiesVec_.at(indexToBeOverwritten);

    for (std::vector<std::size_t>::const_iterator itr = ccgp.eventVecInds_.begin(),
           itrEnd = ccgp.eventVecInds_.end(); itr != itrEnd; ++itr) {
      freeCellCenEventVecInds_.push_back(*itr);
      boost::apply_visitor(clearEventExecutor_, cellCenEventVec_[freeCellCenEventVecInds_.back()]);
    }

    if (indexToBeOverwritten != lastIndex) {

      ccgp = cellCenGroupPropensitiesVec_.back();
    
      int idToBeRedirected = cellCenGroupPropensitiesIdIndexBimap_.right.at(lastIndex);

      cellCenGroupPropensitiesIdIndexBimap_.left.erase(eventGroupId);
      cellCenGroupPropensitiesIdIndexBimap_.left.erase(idToBeRedirected);
      cellCenGroupPropensitiesIdIndexBimap_.insert(IdIndexBimap_::value_type(idToBeRedirected, indexToBeOverwritten));
    }
    else {
      cellCenGroupPropensitiesIdIndexBimap_.left.erase(eventGroupId);
    }
  }
  // Because of the braces surrounding the above, ccgp
  // isn't defined in this scope, so the following pop_back()
  // shouldn't allow ccgp to point to an invalid value.

  cellCenGroupPropensitiesVec_.pop_back();
}

void Simulation::Impl_::mkReversedOffsets_() {

  // Using an *ordered* set because an unordered set affects the reproducibility of simulations.
  std::set<CellIndsOffset> reversedOffsetsSet;

  for (std::vector<CellCenteredGroupPropensities_>::const_iterator ccGPropItr = cellCenGroupPropensitiesVec_.begin(),
         ccGPropItrEnd = cellCenGroupPropensitiesVec_.end(); ccGPropItr != ccGPropItrEnd; ++ccGPropItr) {

    const std::vector<CellIndsOffset> & cioVec = ccGPropItr->cioVec_;
    std::size_t cioVecSize = cioVec.size();

    // i starts at 1 instead of zero in order to avoid the zero offset
    // that's at cioVec[0].
    for (std::size_t i = 1; i < cioVecSize; ++i) {
      reversedOffsetsSet.insert(-(cioVec[i]));
    }

  }

  reversedOffsetsVec_.clear();  
  reversedOffsetsVec_.reserve(reversedOffsetsSet.size());

  for (std::set<CellIndsOffset>::const_iterator itr = reversedOffsetsSet.begin(),
	 itrEnd = reversedOffsetsSet.end(); itr != itrEnd; ++itr) {
    reversedOffsetsVec_.push_back(*itr);
  }

}

void Simulation::Impl_::updateEventAndAddrMapsFromChangedCell_(int currHeight,
                                                               const CellInds & ciOrig) {
  
  CellInds ci(ciOrig);

  lattice_.wrapIndsIfNeeded(ci);

  if ((ci.k >= 0) && (ci.k < currHeight)) {

#if KMC_PARALLEL
    int sectNum;
    bool isGhost = lattice_.addToExportBufferIfNeeded(ci, sectNum);

    if (!isGhost) {
#else
      const int sectNum = 0;
#endif
      doForCellCenteredGroupPropensities_(ci, sectNum,
                                          addOrUpdateCellCenteredEntryToEventList_);
#if KMC_PARALLEL
    }
#endif
  }
}

void Simulation::Impl_::updateEventAndAddrMapsFromAffectedCell_(int currHeight,
                                                                CellInds & ci) {
  lattice_.wrapIndsIfNeeded(ci);

  if ((ci.k >= 0) && (ci.k < currHeight)) {

#if KMC_PARALLEL
    int sectNum = lattice_.sectorOfIndices(ci);

    if (!(sectNum < 0)) {
#else
      const int sectNum = 0;
#endif

      doForCellCenteredGroupPropensities_(ci, sectNum,
                                          addOrUpdateCellCenteredEntryToEventList_);
#if KMC_PARALLEL
    }
#endif
  }
}

void Simulation::Impl_::updateEventAndAddrMapsFromChangedCellInds_(const Lattice::ChangedCellInds & ccInds,
								   const Lattice::OtherCheckedCellInds & ccIndsOther) {

  // Here, I should only rely on ccInds and ccIndsOther having the member
  // functions "size()", "begin()", and "end()".

  CellInds ciPlusOffset;
  int currHeight = lattice_.currHeight();

  for (Lattice::ChangedCellInds::const_iterator ccIndsItr = ccInds.begin(),
	 ccIndsItrEnd = ccInds.end(); ccIndsItr != ccIndsItrEnd;
       ++ccIndsItr) {

    const CellInds & ci = *ccIndsItr;
    updateEventAndAddrMapsFromChangedCell_(currHeight, ci);

    for (std::vector<CellIndsOffset>::const_iterator itr = reversedOffsetsVec_.begin(),
	   itrEnd = reversedOffsetsVec_.end(); itr != itrEnd; ++itr) {

      ciPlusOffset = ci + *itr;
      updateEventAndAddrMapsFromAffectedCell_(currHeight, ciPlusOffset);

    }
  }

  for (Lattice::OtherCheckedCellInds::const_iterator ccIndsItr = ccIndsOther.begin(),
	 ccIndsItrEnd = ccIndsOther.end(); ccIndsItr != ccIndsItrEnd;
       ++ccIndsItr) {

    CellInds ci = *ccIndsItr;
    updateEventAndAddrMapsFromAffectedCell_(currHeight, ci);

  }

}

void Simulation::Impl_::updateEventAndAddrMapsFromAffectedCellOffsets_(const CellsToChange & ctcVec,
                                                                       const std::vector<CellIndsOffset> & affectedCellOffsets,
                                                                       const Lattice::OtherCheckedCellInds & ccIndsOther) {

  int currHeight = lattice_.currHeight();

  const std::vector<CellInds> & changedCellInds = ctcVec.getCellIndsVec();
  
  for (std::vector<CellInds>::const_iterator ciItr = changedCellInds.begin(),
         ciItrEnd = changedCellInds.end(); ciItr != ciItrEnd; ++ciItr) {
    updateEventAndAddrMapsFromChangedCell_(currHeight, *ciItr);
  }

  const CellInds & ciCenter = ctcVec.getCenter();

  for (std::vector<CellIndsOffset>::const_iterator affOffsetItr = affectedCellOffsets.begin(),
         affOffsetItrEnd = affectedCellOffsets.end(); affOffsetItr != affOffsetItrEnd;
       ++affOffsetItr) {

    CellInds ci = ciCenter + *affOffsetItr;
    updateEventAndAddrMapsFromAffectedCell_(currHeight, ci);
  }

  // Here, I should only rely on ccIndsOther having the member
  // functions "size()", "begin()", and "end()".  
  for (Lattice::OtherCheckedCellInds::const_iterator ccIndsItr = ccIndsOther.begin(),
	 ccIndsItrEnd = ccIndsOther.end(); ccIndsItr != ccIndsItrEnd;
       ++ccIndsItr) {

    CellInds ci = *ccIndsItr;
    updateEventAndAddrMapsFromAffectedCell_(currHeight, ci);

  }

}

void Simulation::Impl_::rebuildEventAndAddrMaps_() {
  
  int numSectors = lattice_.numSectors();

#if KMC_PARALLEL
  for (int i = 0; i < numSectors; ++i) {
    lattice_.recvGhosts(i);
  }
#endif

  solver_->beginBuildingEventList(overLatticeEventVec_.size(),
                                  lattice_.planesReserved());

  for (int sectNum = 0; sectNum < numSectors; ++sectNum) {

    int kmaxP1 = lattice_.currHeight();

    CellInds ci;
    //std::vector<double> propensitiesVec;

    for (ci.k = 0; ci.k < kmaxP1; ++(ci.k)) {
      for (ci.i = sectorPlanarBBox_[sectNum].imin; ci.i < sectorPlanarBBox_[sectNum].imaxP1; ++(ci.i)) {
        for (ci.j = sectorPlanarBBox_[sectNum].jmin; ci.j < sectorPlanarBBox_[sectNum].jmaxP1; ++(ci.j)) {

          doForCellCenteredGroupPropensities_(ci,
                                              sectNum,
                                              addCellCenteredEntryToEventList_);

        }
      }
    }

    for (int overLatticeEventIndex = 0; overLatticeEventIndex < overLatticeEventVec_.size();
	 ++overLatticeEventIndex) {

      double p = overLatticeEventVec_[overLatticeEventIndex].propensity_[sectNum];

      if (p > 0) {
	EventId eId(overLatticeEventIndex, sectNum);
	solver_->addOverLatticeEntryToEventList(eId, p, sectNum);
      }
    }
 
  }

  solver_->endBuildingEventList();

}

void Simulation::Impl_::doPreRunChecks_() const {
  bool initError = false;
  std::string initErrStr;

  if (!solver_) {
    initError = true;
    initErrStr += "Must set solver before running simulation.\n";
  }   

  if (!rngIsSet_) {
    initError = true;
    initErrStr += "Must set RNG before running simulation.\n";
  }

#if KMC_PARALLEL
  if (!tIncrSchemeIsSet_) {
    initError = true;
    initErrStr += "Must set time incrementing scheme for the parallel algorithm before running simulation.\n";
  }
#endif

  if (cellCenEventVec_.empty() && overLatticeEventVec_.empty()) {    
    initError = true;
    initErrStr += "Must add events to simulation before running it.\n";
  }

  if (lattice_.currHeight() < 1) {
    initError = true;
    initErrStr += "Lattice must have at least one (possibly empty) plane.\n";
  }

  exitOnCondition(initError, initErrStr);  
}

void Simulation::Impl_::executeEvent_(const EventId & chosenEventID) {
  ++(simState_.num_local_events_);

  if (chosenEventID.isForOverLattice()) {
	    
    int overLatticeEventIndex, sectNum;
    chosenEventID.getEventInfo(overLatticeEventIndex, sectNum);

    OverLatticeEvent_ & chosenEvent = overLatticeEventVec_[overLatticeEventIndex];

    int iExtents = sectorPlanarBBox_[sectNum].imaxP1 - sectorPlanarBBox_[sectNum].imin;
    double iExtentsFrac = rng_->getNumInOpenIntervalFrom0To1();

    int jExtents = sectorPlanarBBox_[sectNum].jmaxP1 - sectorPlanarBBox_[sectNum].jmin;
    double jExtentsFrac = rng_->getNumInOpenIntervalFrom0To1();

    runEventExecutor_.ci.i = sectorPlanarBBox_[sectNum].imin + iExtents*iExtentsFrac;
    runEventExecutor_.ci.j = sectorPlanarBBox_[sectNum].jmin + jExtents*jExtentsFrac;
    runEventExecutor_.ci.k = lattice_.currHeight() - 1;

    if (runEventExecutor_.ci.i == sectorPlanarBBox_[sectNum].imaxP1) { --(runEventExecutor_.ci.i); }
    if (runEventExecutor_.ci.j == sectorPlanarBBox_[sectNum].jmaxP1) { --(runEventExecutor_.ci.j); }

    boost::apply_visitor(runEventExecutor_, chosenEvent.execute_);

  }
  else {

    int cellCenEventIndex;
    chosenEventID.getEventInfo(runEventExecutor_.ci, cellCenEventIndex);

    EventExecutor_ & chosenEvent = cellCenEventVec_[cellCenEventIndex];

    boost::apply_visitor(runEventExecutor_, chosenEvent);

  }
  
}

void Simulation::Impl_::updateEventAndAddrMapsAfterPeriodicActionsWTrack_() {

  const Lattice::ChangedCellInds & ccInds = lattice_.getChangedCellInds();
  const Lattice::OtherCheckedCellInds & ccIndsOther = lattice_.getOtherCheckedCellInds();

#if KMC_PARALLEL
  int numSectors = lattice_.numSectors();

  for (Lattice::ChangedCellInds::const_iterator ccIndsItr = ccInds.begin(),
         ccIndsItrEnd = ccInds.end(); ccIndsItr != ccIndsItrEnd;
       ++ccIndsItr) {
    lattice_.addToExportBufferIfNeeded(*ccIndsItr);
  }

  for (Lattice::OtherCheckedCellInds::const_iterator ccIndsItr = ccIndsOther.begin(),
         ccIndsItrEnd = ccIndsOther.end(); ccIndsItr != ccIndsItrEnd;
       ++ccIndsItr) {
    lattice_.addToExportBufferIfNeeded(*ccIndsItr);    
  }

  // I should only be receiving ghosts from other sites, not sending them.
  lattice_.clearGhostsToSend();

  // Accounting for changes to event lists due to presence of ghosts
  for (int i = 0; i < numSectors; ++i) {
    lattice_.recvGhostsUpdate(i);
    updateEventAndAddrMapsAffectedByGhostUpdates_(lattice_.getReceivedGhostInds());
  }
  
  // Accounting for changes to event lists due to the local effects of
  // the periodic action, now that the ghosts have been taken into
  // account. What's below is similar to what is done by
  // updateEventAndAddrMapsFromChangedCellInds_, except that there are
  // no calls to lattice_.addToExportBufferIfNeeded(), because those
  // calls have already been done.

  CellInds ciPlusOffset;
  int currHeight = lattice_.currHeight();

  for (Lattice::ChangedCellInds::const_iterator ccIndsItr = ccInds.begin(),
         ccIndsItrEnd = ccInds.end(); ccIndsItr != ccIndsItrEnd;
       ++ccIndsItr) {
    
    const CellInds & ci = *ccIndsItr;

    int sectNum = lattice_.sectorOfIndices(ci);
    if (!(sectNum < 0)) {
      // Not adding to export buffer, since that was already done
      doForCellCenteredGroupPropensities_(ci, sectNum,
                                          addOrUpdateCellCenteredEntryToEventList_);
    }

    for (std::vector<CellIndsOffset>::const_iterator oitr = reversedOffsetsVec_.begin(),
           oitrEnd = reversedOffsetsVec_.end(); oitr != oitrEnd; ++oitr) {

      ciPlusOffset = ci + *oitr;
      updateEventAndAddrMapsFromAffectedCell_(currHeight, ciPlusOffset);        
    }

  }

  for (Lattice::OtherCheckedCellInds::const_iterator ccIndsItr = ccIndsOther.begin(),
         ccIndsItrEnd = ccIndsOther.end(); ccIndsItr != ccIndsItrEnd;
       ++ccIndsItr) {
    
    const CellInds & ci = *ccIndsItr;

    int sectNum = lattice_.sectorOfIndices(ci);
    if (!(sectNum < 0)) {
      // Not adding to export buffer, since that was already done
      doForCellCenteredGroupPropensities_(ci, sectNum,
                                          addOrUpdateCellCenteredEntryToEventList_);
    }
    
  }

#else
  updateEventAndAddrMapsFromChangedCellInds_(ccInds, ccIndsOther);
#endif

}

void Simulation::Impl_::updateEventAndAddrMapsAfterPeriodicActionsNoTrack_() {
  rebuildEventAndAddrMaps_();
}

void Simulation::Impl_::runPeriodicActions_() {

  lattice_.trackChanges(changeTrackingForPeriodicAction_);

  for (std::vector<TimePeriodicAction_>::iterator itr = timePeriodicActionVec_.begin(),
	 itrEnd = timePeriodicActionVec_.end(); itr != itrEnd; ++itr) {
    itr->doActionIfEndOfPeriod(simState_, lattice_);
  }

  for (std::vector<StepPeriodicAction_>::iterator itr = stepPeriodicActionVec_.begin(),
	 itrEnd = stepPeriodicActionVec_.end(); itr != itrEnd; ++itr) {    
    itr->doActionIfEndOfPeriod(simState_, lattice_);
  }

  if (lattice_.hasChanged()) {
    KMC_CALL_MEMBER_FUNCTION(*this, updateEventAndAddrMapsAfterPeriodicActions_)();
  }

  lattice_.trackChanges(Lattice::TrackType::NONE);
}

void Simulation::Impl_::run_(double runTime) {

  doPreRunChecks_();

  simState_.maxTime_ += runTime;

  mkReversedOffsets_();

  SetEventExecutorAffectedCellOffsets_ setEventExecutorAffectedCellOffsets_(&reversedOffsetsVec_);

  for (std::vector<EventExecutor_>::iterator itr = cellCenEventVec_.begin(),
         itrEnd = cellCenEventVec_.end(); itr != itrEnd; ++itr) {
    boost::apply_visitor(setEventExecutorAffectedCellOffsets_, *itr);
  }

  for (std::vector<OverLatticeEvent_>::iterator itr = overLatticeEventVec_.begin(),
         itrEnd = overLatticeEventVec_.end(); itr != itrEnd; ++itr) {
    boost::apply_visitor(setEventExecutorAffectedCellOffsets_, itr->execute_);
  }

  // Initializing ciMin_ array used by EventId.  
  EventId::ciMin_[0] = localPlanarBBox_.imin;
  EventId::ciMin_[1] = localPlanarBBox_.jmin;

  // Initializing dimsForFlattening_ array used by EventId.
  EventId::dimsForFlattening_[0] = localPlanarBBox_.imaxP1 - localPlanarBBox_.imin;
  EventId::dimsForFlattening_[1] = localPlanarBBox_.jmaxP1 - localPlanarBBox_.jmin;
  EventId::dimsForFlattening_[2] = cellCenEventVec_.size();

  // Initializing the event maps
  rebuildEventAndAddrMaps_();

  EventId chosenEventID;

#if KMC_PARALLEL
  int numSectors = lattice_.numSectors();

  if ((!runHasBeenExecuted_) && solver_->timeIncrSchemeIsAdaptive()) {

    // Using an integer for a "Boolean" because MPI C interface doesn't support bools.
    int aLocalSectorLacksCellCenteredEvents = 0;
    for (int i = 0; i < numSectors; ++i) {
      if (solver_->noCellCenteredEvents(i)) {
	aLocalSectorLacksCellCenteredEvents = 1;
	break;
      }
    }
    
    int aSectorLacksCellCenteredEvents;
    MPI_Allreduce(&aLocalSectorLacksCellCenteredEvents,
		  &aSectorLacksCellCenteredEvents, 1,
		  MPI_INT, MPI_MAX, lattice_.comm());

    if (aSectorLacksCellCenteredEvents) {

      // This "jumpstarts" the simulation by executing events that
      // introduce cell-centered events into the system, which are then
      // used in adaptive schemes for calculating t_stop_.
      std::vector<double> jumpstartTime(numSectors, 0);

      for (int i = 0; i < numSectors; ++i) {

	lattice_.recvGhostsUpdate(i);
	updateEventAndAddrMapsAffectedByGhostUpdates_(lattice_.getReceivedGhostInds());

	if (!(solver_->noMoreEvents(i))) {
	  solver_->chooseEventIDAndUpdateTime(i, chosenEventID, jumpstartTime[i]);
	  
	  // This is needed for simState_.elapsedTime() to return the
	  // correct value within the loop.
	  simState_.t_sector_ = jumpstartTime[i];

	  executeEvent_(chosenEventID);
	}

	lattice_.sendGhostsUpdate(i);
        updateEventAndAddrMapsAffectedByGhostUpdates_(lattice_.getReceivedLocalInds());
      }

      // The amount of time that has elapsed is the time of the longest event.
      double jumpStartElapsedTimeLocal = *(std::max_element(jumpstartTime.begin(),
							    jumpstartTime.end()));

      double jumpStartElapsedTime;
      MPI_Allreduce(&jumpStartElapsedTimeLocal, &jumpStartElapsedTime,
		    1, MPI_DOUBLE, MPI_MAX, lattice_.comm());

      simState_.elapsed_time_ += jumpStartElapsedTime;
    }
  }

  solver_->updateTStop(lattice_.comm(), simState_.t_stop_);

  exitOnCondition(!(simState_.t_stop_ > 0),
		  "Time step is not greater than zero. Cannot continue with simulation.");
#endif

#if KMC_PARALLEL
  // This compensates for t_sector_ being set to zero at the bottom of
  // the "for" loop below.
  simState_.t_sector_ = 0;
#endif

  // Now doing the actual KMC runs.
  while (simState_.elapsed_time_ < simState_.maxTime_) {

#if KMC_PARALLEL
    
    for (int i = 0; i < numSectors; ++i) {

      lattice_.recvGhostsUpdate(i);
      updateEventAndAddrMapsAffectedByGhostUpdates_(lattice_.getReceivedGhostInds());

      while (true) {

	// Move onto the next sector if there are no events to execute.
        if (solver_->noMoreEvents(i)) {
          break;
        }

	solver_->chooseEventIDAndUpdateTime(i, chosenEventID, simState_.t_sector_);

	if (simState_.t_sector_ <= simState_.t_stop_) {
	  executeEvent_(chosenEventID);
	}
	else {
	  break;
	}

      }

      lattice_.sendGhostsUpdate(i);
      updateEventAndAddrMapsAffectedByGhostUpdates_(lattice_.getReceivedLocalInds());
      
      // Updating t_sector_ *here* rather than at the beginning of the
      // loop so that t_sector_ = 0 outside of the sector loop. (This
      // way, simState_.elapsedTime() will return the correct result
      // outside of a loop over sectors.
      simState_.t_sector_ = 0;
    }

    simState_.elapsed_time_ += simState_.t_stop_;
    ++(simState_.num_global_steps_);

    runPeriodicActions_();

    solver_->updateTStop(lattice_.comm(), simState_.t_stop_);
#else
    
    if (solver_->noMoreEvents(0)) {
      std::cout << "Simulation ran out of events to execute at simulation time = "
		<< simState_.elapsed_time_ << "\n";
      break;
    }

    solver_->chooseEventIDAndUpdateTime(0, chosenEventID, simState_.elapsed_time_);
    executeEvent_(chosenEventID);

    ++(simState_.num_global_steps_);

    runPeriodicActions_();
#endif

  }

  runHasBeenExecuted_ = true;
}

#if KMC_PARALLEL
void Simulation::Impl_::updateEventAndAddrMapsAffectedByGhostUpdates_(const std::vector<std::vector<IJK> > & receivedInds) {

  std::size_t receivedIndsSize = receivedInds.size();

  CellInds ci, ciPlusOffset;
  int currHeight = lattice_.currHeight();

  for (std::size_t bufType = 0; bufType < receivedIndsSize; ++bufType) {
    const std::vector<IJK> & currInds = receivedInds[bufType];

    for (std::vector<IJK>::const_iterator itr = currInds.begin(),
           itrEnd = currInds.end(); itr != itrEnd; ++itr) {

      ci.i = itr->i;
      ci.j = itr->j;
      ci.k = itr->k;      

      int sectNum = lattice_.sectorOfIndices(ci);

      if (!(sectNum < 0)) {
        doForCellCenteredGroupPropensities_(ci, sectNum,
                                            addOrUpdateCellCenteredEntryToEventList_);
      }

      for (std::vector<CellIndsOffset>::const_iterator oitr = reversedOffsetsVec_.begin(),
           oitrEnd = reversedOffsetsVec_.end(); oitr != oitrEnd; ++oitr) {

        ciPlusOffset = ci + *oitr;
        updateEventAndAddrMapsFromAffectedCell_(currHeight, ciPlusOffset);        
      }

    }

  }

}
#endif

Simulation::Simulation(const LatticeParams & paramsForLattice)
  : pImpl_(new Impl_(paramsForLattice)) {

  trackCellsChangedByPeriodicActions(false);
}

Simulation::~Simulation() {}

int Simulation::nProcs() const {return pImpl_->lattice_.nProcs();}
int Simulation::procID() const {return pImpl_->lattice_.procID();}

#if KMC_PARALLEL
const MPI_Comm & Simulation::comm() const {return pImpl_->lattice_.comm();}
#endif

int Simulation::procPerDim(int dim) const {return pImpl_->lattice_.procPerDim(dim);}
int Simulation::commCoord(int dim) const {return pImpl_->lattice_.commCoord(dim);}

void Simulation::getLatticeLocalPlanarBBox(bool wGhost,
				    LatticePlanarBBox & bbox) const {
  pImpl_->lattice_.getLocalPlanarBBox(wGhost, bbox);
}

void Simulation::getLatticeSectorPlanarBBox(int sectNum,
					    LatticePlanarBBox & bbox) const {
  pImpl_->lattice_.getSectorPlanarBBox(sectNum, bbox);
}

void Simulation::getLatticeGlobalPlanarBBox(LatticePlanarBBox & bbox) const {
  pImpl_->lattice_.getGlobalPlanarBBox(bbox);
}

double Simulation::elapsedTime() const {return pImpl_->simState_.elapsed_time_;}
unsigned long long Simulation::numLocalEvents() const {return pImpl_->simState_.num_local_events_;}
unsigned long long Simulation::numGlobalSteps() const {return pImpl_->simState_.num_global_steps_;}

void Simulation::reserveCellCenteredEventGroups(int numGroups, int numTotEvents) {
  pImpl_->cellCenGroupPropensitiesVec_.reserve(numGroups);
  pImpl_->cellCenEventVec_.reserve(numTotEvents);
}

void Simulation::addCellCenteredEventGroup(int eventGroupId,
                                           const CellNeighOffsets & cno,
                                           CellCenteredGroupPropensities propensities,
                                           const EventExecutorGroup & eventExecutorGroup) {

  pImpl_->cellCenGroupPropensitiesIdIndexBimap_.insert(Impl_::IdIndexBimap_::value_type(eventGroupId,
                                                                                        pImpl_->cellCenGroupPropensitiesVec_.size()));

  pImpl_->cellCenGroupPropensitiesVec_.push_back(Impl_::CellCenteredGroupPropensities_(cno,
                                                                                       propensities));

  Impl_::CellCenteredGroupPropensities_ & ccgp = pImpl_->cellCenGroupPropensitiesVec_.back();

  ccgp.eventVecInds_.reserve(eventExecutorGroup.numEventExecutors());

  for (int i = 0; i < eventExecutorGroup.numEventExecutors(); ++i) {

    std::size_t eventVecInd;

    if (pImpl_->freeCellCenEventVecInds_.empty()) {
      eventVecInd = pImpl_->cellCenEventVec_.size();
      pImpl_->cellCenEventVec_.push_back(Impl_::EventExecutor_());
    }
    else {
      eventVecInd = pImpl_->freeCellCenEventVecInds_.front();
      pImpl_->freeCellCenEventVecInds_.pop_front();
    }

    ccgp.eventVecInds_.push_back(eventVecInd);

    switch (eventExecutorGroup.getEventExecutorType(i)) {
    case EventExecutorGroup::EventExecEnum::AUTO:
      pImpl_->cellCenEventVec_[eventVecInd] = eventExecutorGroup.getEventExecutorAutoTrack(i);
      break;
    case EventExecutorGroup::EventExecEnum::SEMIMANUAL:

      pImpl_->cellCenEventVec_[eventVecInd] = Impl_::EventExecutorSemiManualTrackInfo_(eventExecutorGroup.getEventExecutorSemiManualTrack(i),
                                                                                       &(pImpl_->lattice_),
                                                                                       eventExecutorGroup.getEventExecutorOffsetsVec(i));
      break;
    }
    
  }

}

void Simulation::changeCellCenteredEventGroup(int eventGroupId,
                                              const CellNeighOffsets & cno,
                                              CellCenteredGroupPropensities propensities,
                                              const EventExecutorGroup & eventExecutorGroup) {

  pImpl_->removeCellCenteredEventGroup_(eventGroupId, "changeCellCenteredEventGroup");
  addCellCenteredEventGroup(eventGroupId, cno, propensities, eventExecutorGroup);  
}

void Simulation::removeCellCenteredEventGroup(int eventGroupId) {
  pImpl_->removeCellCenteredEventGroup_(eventGroupId, "removeCellCenteredEventGroup");
}

void Simulation::reserveOverLatticeEvents(int num) {
  pImpl_->overLatticeEventVec_.reserve(num);
}

void Simulation::addOverLatticeEvent(int eventId,
				     double propensityPerUnitArea,
				     EventExecutorAutoTrack eventExecutor) {
  pImpl_->overLatticeEventIdIndexBimap_.insert(Impl_::IdIndexBimap_::value_type(eventId,
										pImpl_->overLatticeEventVec_.size()));

  pImpl_->overLatticeEventVec_.push_back(Impl_::OverLatticeEvent_(propensityPerUnitArea,
								  pImpl_->sectorPlanarBBox_,
								  eventExecutor));
}

void Simulation::addOverLatticeEvent(int eventId,
				     double propensityPerUnitArea,
                                     EventExecutorSemiManualTrack eventExecutor,
                                     const std::vector<CellNeighOffsets> & cnoVec) {
  pImpl_->overLatticeEventIdIndexBimap_.insert(Impl_::IdIndexBimap_::value_type(eventId,
										pImpl_->overLatticeEventVec_.size()));

  pImpl_->overLatticeEventVec_.push_back(Impl_::OverLatticeEvent_(propensityPerUnitArea,
								  pImpl_->sectorPlanarBBox_,
								  Impl_::EventExecutorSemiManualTrackInfo_(eventExecutor,
                                                                                                           &(pImpl_->lattice_),
                                                                                                           cnoVec)));
}

void Simulation::changeOverLatticeEvent(int eventId,
					double propensityPerUnitArea,
					EventExecutorAutoTrack eventExecutor) {

  // Entering a bad eventId is a likely error, so I'll try to catch it.
  std::size_t overLatticeEventVecIndex = pImpl_->bimapIdToIndex_(eventId,
								 pImpl_->overLatticeEventIdIndexBimap_,
								 "changeOverLatticeEvent",
								 "eventId");
  
  pImpl_->overLatticeEventVec_.at(overLatticeEventVecIndex) = Impl_::OverLatticeEvent_(propensityPerUnitArea,
										       pImpl_->sectorPlanarBBox_,
										       eventExecutor);
}

void Simulation::changeOverLatticeEvent(int eventId,
					double propensityPerUnitArea,
					EventExecutorSemiManualTrack eventExecutor,
                                        const std::vector<CellNeighOffsets> & cnoVec) {

  // Entering a bad eventId is a likely error, so I'll try to catch it.
  std::size_t overLatticeEventVecIndex = pImpl_->bimapIdToIndex_(eventId,
								 pImpl_->overLatticeEventIdIndexBimap_,
								 "changeOverLatticeEvent",
								 "eventId");
  
  pImpl_->overLatticeEventVec_.at(overLatticeEventVecIndex) = Impl_::OverLatticeEvent_(propensityPerUnitArea,
										       pImpl_->sectorPlanarBBox_,
										       Impl_::EventExecutorSemiManualTrackInfo_(eventExecutor,
                                                                                                                                &(pImpl_->lattice_),
                                                                                                                                cnoVec));
}


void Simulation::removeOverLatticeEvent(int eventId) {

  pImpl_->removeEntryFromVecWithBimap_(eventId,
				       pImpl_->overLatticeEventIdIndexBimap_,
				       pImpl_->overLatticeEventVec_,
				       "removeOverLatticeEvent",
				       "eventId");
}

void Simulation::trackCellsChangedByPeriodicActions(bool doTrack) {
  if (doTrack) {
    pImpl_->changeTrackingForPeriodicAction_ = Lattice::TrackType::RECORD_CHANGED_CELL_INDS;
    pImpl_->updateEventAndAddrMapsAfterPeriodicActions_ = &Impl_::updateEventAndAddrMapsAfterPeriodicActionsWTrack_;
  }
  else {
    pImpl_->changeTrackingForPeriodicAction_ = Lattice::TrackType::CHECK_ONLY_IF_CHANGE_OCCURS;
    pImpl_->updateEventAndAddrMapsAfterPeriodicActions_ = &Impl_::updateEventAndAddrMapsAfterPeriodicActionsNoTrack_;
  }
}

void Simulation::reserveTimePeriodicActions(int num) {
  pImpl_->timePeriodicActionVec_.reserve(num);
}

void Simulation::addTimePeriodicAction(int actionId,
				       PeriodicAction action,
				       double period,
				       bool doAtSimEnd) {

  pImpl_->timePeriodicActionIdIndexBimap_.insert(Impl_::IdIndexBimap_::value_type(actionId,
										  pImpl_->timePeriodicActionVec_.size()));

  pImpl_->timePeriodicActionVec_.push_back(Impl_::TimePeriodicAction_(action, period, doAtSimEnd));
}

void Simulation::changeTimePeriodicAction(int actionId,
					  PeriodicAction action,
					  double period,
					  bool doAtSimEnd) {

  // Entering a bad actionId is a likely error, so I'll try to catch it.
  std::size_t actionVecIndex = pImpl_->bimapIdToIndex_(actionId, 
						       pImpl_->timePeriodicActionIdIndexBimap_,
						       "changePeriodicAction",
						       "actionId");

  pImpl_->timePeriodicActionVec_.at(actionVecIndex) = Impl_::TimePeriodicAction_(action, period, doAtSimEnd);
}

void Simulation::removeTimePeriodicAction(int actionId) {

  pImpl_->removeEntryFromVecWithBimap_(actionId,
				       pImpl_->timePeriodicActionIdIndexBimap_,
				       pImpl_->timePeriodicActionVec_,
				       "removePeriodicAction",
				       "actionId");
  
}

void Simulation::reserveStepPeriodicActions(int num) {
  pImpl_->stepPeriodicActionVec_.reserve(num);
}

void Simulation::addStepPeriodicAction(int actionId,
				       PeriodicAction action,
				       int period,
				       bool doAtSimEnd) {

  pImpl_->stepPeriodicActionIdIndexBimap_.insert(Impl_::IdIndexBimap_::value_type(actionId,
										  pImpl_->stepPeriodicActionVec_.size()));

  pImpl_->stepPeriodicActionVec_.push_back(Impl_::StepPeriodicAction_(action, period, doAtSimEnd));
}

void Simulation::changeStepPeriodicAction(int actionId,
					  PeriodicAction action,
					  int period,
					  bool doAtSimEnd) {

  // Entering a bad actionId is a likely error, so I'll try to catch it.
  std::size_t actionVecIndex = pImpl_->bimapIdToIndex_(actionId, 
						       pImpl_->stepPeriodicActionIdIndexBimap_,
						       "changeStepPeriodicAction",
						       "actionId");

  pImpl_->stepPeriodicActionVec_.at(actionVecIndex) = Impl_::StepPeriodicAction_(action, period, doAtSimEnd);
}

void Simulation::removeStepPeriodicAction(int actionId) {

  pImpl_->removeEntryFromVecWithBimap_(actionId,
				       pImpl_->stepPeriodicActionIdIndexBimap_,
				       pImpl_->stepPeriodicActionVec_,
				       "removeStepPeriodicAction",
				       "actionId");
  
}

void Simulation::setSolver(SolverId::Type sId) {
  pImpl_->solver_.reset(mkSolver(sId, &(pImpl_->lattice_)));

  // Even if the time increment scheme and random number generator
  // have been set before, they'll need to be reset after the solver
  // has been set.
#if KMC_PARALLEL
  pImpl_->tIncrSchemeIsSet_ = false;
#endif
  pImpl_->rngIsSet_ = false;
}

void Simulation::setTimeIncrScheme(const TimeIncr::SchemeVars & vars) {
#if KMC_PARALLEL
  exitOnCondition(!(pImpl_->solver_),
		  "Time increment scheme must be set after the solver has been set");

  pImpl_->solver_->setTimeIncrScheme(vars);
  pImpl_->tIncrSchemeIsSet_ = true;
#else
  std::cerr << "Warning: using setTimeIncrScheme with the serial version of the library does nothing.\n";
#endif
}

void Simulation::setRNG(RandNumGenSharedPtr rng) {
  exitOnCondition(!(pImpl_->solver_),
		  "Random number generator must be set after the solver has been set");

  pImpl_->rng_ = rng;
  pImpl_->solver_->setRNG(rng);
  pImpl_->rngIsSet_ = true;
}

void Simulation::run(double runTime) {
  pImpl_->run_(runTime);
}
