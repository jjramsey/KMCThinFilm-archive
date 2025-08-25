#include "SolverDynamicSchulze.hpp"
#include "Lattice.hpp"
#include "ErrorHandling.hpp"

using namespace KMCThinFilm;

SolverDynamicSchulze::SolverDynamicSchulze(const Lattice * lattice) 
  : 
#if KMC_PARALLEL
  Solver(lattice),
#endif
  propToEventIdList_(lattice->numSectors()),
  propToEventIdListProxy_(lattice->numSectors())
#if KMC_PARALLEL
  , totOverLatticePropensity_(lattice->numSectors()),
  totNumOverLatticeEvents_(lattice->numSectors())
#endif
{}

void SolverDynamicSchulze::beginBuildingEventList(int numOverLatticeEvents,
                                                  int numReservedLatticePlanes) {
  
  for (std::size_t i = 0; i < propToEventIdList_.size(); ++i) {
    propToEventIdList_[i].clear();
    propToEventIdListProxy_[i].clear();
#if KMC_PARALLEL
    totOverLatticePropensity_[i] = 0;
    totNumOverLatticeEvents_[i] = 0;
#endif
  }

  addrMap_.reset(new EventIdMap<EvListItrIndexPair_>(propToEventIdList_.size(),
                                                     numOverLatticeEvents,
                                                     numReservedLatticePlanes,
                                                     EvListItrIndexPair_(propToEventIdList_[0].begin(), -1)));
}

SolverDynamicSchulze::PropToEventIdList_::iterator SolverDynamicSchulze::getPropToEventIdListItr_(double propensity,
                                                                                                  int sectNum) {
  PropToEventIdList_::iterator propToEventIdListItr =
    propToEventIdList_[sectNum].find(propensity);

  if (propToEventIdListItr == propToEventIdList_[sectNum].end()) {
    std::pair<PropToEventIdList_::iterator,bool> itrBoolPair = 
      propToEventIdList_[sectNum].insert(std::make_pair(propensity,
                                                        IndexDequePair_(propToEventIdListProxy_[sectNum].size())));

    propToEventIdListItr = itrBoolPair.first;
    propToEventIdListProxy_[sectNum].push_back(PropToEventIdListProxyEntry_(propToEventIdListItr));
  }

  return propToEventIdListItr;
}

void SolverDynamicSchulze::addCellCenteredEntryToEventList(const EventId & eId,
							   double propensity,
							   int sectNum) {

  PropToEventIdList_::iterator propToEventIdListItr = getPropToEventIdListItr_(propensity,
                                                                               sectNum);

  addItrEntryToEventList_(eId, propToEventIdListItr, sectNum);

}

void SolverDynamicSchulze::addOverLatticeEntryToEventList(const EventId & eId,
							  double propensity,
							  int sectNum) {

  PropToEventIdList_::iterator propToEventIdListItr = getPropToEventIdListItr_(propensity,
                                                                               sectNum);
#if KMC_PARALLEL
  ++(propToEventIdListItr->second.numOverLatticeEvents);
  ++(totNumOverLatticeEvents_[sectNum]);

  totOverLatticePropensity_[sectNum] += propensity;
#endif

  addItrEntryToEventList_(eId, propToEventIdListItr, sectNum);
}

void SolverDynamicSchulze::addItrEntryToEventList_(const EventId & eId,
						   PropToEventIdList_::iterator propToEventIdListItr,
						   int sectNum) {

  std::deque<EventId> & evIdList = propToEventIdListItr->second.eIdDeque;

  std::size_t evListInd = evIdList.size();
  evIdList.push_back(eId);

  std::size_t indexToEIdListProxy = propToEventIdListItr->second.indexToEIdListProxy;
  propToEventIdListProxy_[sectNum][indexToEIdListProxy].recalcPartialSumContrib = true;

  addrMap_->addOrUpdate(eId, EvListItrIndexPair_(propToEventIdListItr, evListInd));
}

void SolverDynamicSchulze::removeFromEventIdList_(std::size_t origInd, std::deque<EventId> & origDeque) {

  if ((origInd + 1) != origDeque.size()) {
    // Replace the entry in origDeque originally associated with origInd
    // with the entry at the rear.
    origDeque[origInd] = origDeque.back();

    // The rear entry is now a duplicate of origDeque[origInd]
    // and should be removed.
    origDeque.pop_back();

    // Update addrMap_ to reflect the replacement.
    addrMap_->getRefToVal(origDeque[origInd]).indexToEId = origInd;
    
  }
  else {
    // If entry in origDeque originally associated with eId *is*
    // the rear entry, then I can just remove the entry.
    origDeque.pop_back();
  }

}

void SolverDynamicSchulze::addOrUpdateCellCenteredEntryToEventList(const EventId & eId,
								   double currPropensity,
								   int sectNum) {

  EvListItrIndexPair_ * evListItrIndexPairPtr = addrMap_->getPtrToVal(eId);

  if ((evListItrIndexPairPtr == NULL) || (evListItrIndexPairPtr->indexToEId < 0)) {

    if (currPropensity > 0) {
      /* If it wasn't in the address map before, but now has a
	 non-zero propensity, then it should be in the event and
	 address maps now. */

      // In this *particular* context, I can use
      // addCellCenteredEntryToEventList without being in between
      // calls to beginBuildingEventList() and endBuildingEventList().
      addCellCenteredEntryToEventList(eId, currPropensity, sectNum);
    }

  }
  else {

    PropToEventIdList_::iterator origItr = evListItrIndexPairPtr->itr;
    std::size_t origInd = evListItrIndexPairPtr->indexToEId;
    
    std::deque<EventId> & origDeque = origItr->second.eIdDeque;

    if (currPropensity > 0) {

      if (origItr->first != currPropensity) {

	/* If eId already has existing entries in the event ID lists
	   and address map, but its propensity has changed, then those
	   entries need to be updated. */

        // Removing eId from its current place in the event ID list.
        // (Note that this *may* involve modifying addrMap_.)
        removeFromEventIdList_(origInd, origDeque);
	std::size_t origIndexToEIdListProxy = origItr->second.indexToEIdListProxy;
        propToEventIdListProxy_[sectNum][origIndexToEIdListProxy].recalcPartialSumContrib = true;

        PropToEventIdList_::iterator newItr =  getPropToEventIdListItr_(currPropensity,
									sectNum);

        std::deque<EventId> & currDeque = newItr->second.eIdDeque;
        evListItrIndexPairPtr->indexToEId = currDeque.size();
        currDeque.push_back(eId);

	std::size_t currIndexToEIdListProxy = newItr->second.indexToEIdListProxy;
        propToEventIdListProxy_[sectNum][currIndexToEIdListProxy].recalcPartialSumContrib = true;

        evListItrIndexPairPtr->itr = newItr;
      }

    }
    else {
      /* If the current propensity is zero, then eId is associated
         with an event that can't happen.

         Therefore, I effectively remove eId from its current place in
         the event ID list and then its associated entry in addrMap_.
      */

      removeFromEventIdList_(origInd, origDeque);
      std::size_t origIndexToEIdListProxy = origItr->second.indexToEIdListProxy;
      propToEventIdListProxy_[sectNum][origIndexToEIdListProxy].recalcPartialSumContrib = true;
      evListItrIndexPairPtr->indexToEId = -1;
    }

    if (origDeque.empty()) {
      // I don't want propensities to be associated with empty event ID lists.
      
      std::size_t origIndexToEIdListProxy = origItr->second.indexToEIdListProxy;
      PropToEventIdListProxy_ & currPropToEventIdListProxy = propToEventIdListProxy_[sectNum];
      
      if ((origIndexToEIdListProxy + 1) != currPropToEventIdListProxy.size()) {
        
        // Replace the entry in currPropToEventIdListProxy originally
        // associated with origIndexToEIdListProxya with the entry at
        // the rear
        currPropToEventIdListProxy[origIndexToEIdListProxy] = currPropToEventIdListProxy.back();

        // Update itr->second.indexToEIdListProxy to reflect the replacement.
        currPropToEventIdListProxy[origIndexToEIdListProxy].itr->second.indexToEIdListProxy = origIndexToEIdListProxy;
      }

      // If the above condition is true, then the last entry is
      // redundant with
      // currPropToEventIdListProxy[origIndexToEIdListProxy] and
      // should be removed. If the above condition is false, then the
      // last entry is the one that was supposed to be removed in the
      // first place.
      currPropToEventIdListProxy.pop_back();

      propToEventIdList_[sectNum].erase(origItr);
    }

  }

}

void SolverDynamicSchulze::chooseEventIDAndUpdateTime(int sectNum,
						      EventId & chosenEventID,
						      double & time) {
  // Choose event to execute (but don't execute it just yet) and
  // update t_sector. Method for choosing the event is a
  // slightly modified version of the one from Schulze, Physical
  // Review E, v. 65, 036704.
    
  PropToEventIdListProxy_ & currPropToEventIdListProxy = propToEventIdListProxy_[sectNum];
  
  double p_s = 0;

  // It's faster to iterate over the PropToEventIdListProxy_ container
  // than the PropToEventIdList_ container.
  for (PropToEventIdListProxy_::iterator itr = currPropToEventIdListProxy.begin(),
         itrEnd = currPropToEventIdListProxy.end(); itr != itrEnd; ++itr) {

    if (itr->recalcPartialSumContrib) {
      PropToEventIdList_::iterator mapItr = itr->itr;
      itr->savedPartialSumContrib = (mapItr->first)*(mapItr->second.eIdDeque.size());
      itr->recalcPartialSumContrib = false;
    }

    p_s += itr->savedPartialSumContrib;
    itr->partialSum = p_s;
  }

  double R = p_s*(rng_->getNumInOpenIntervalFrom0To1());

  PropToEventIdListProxy_::iterator chosenItr = std::lower_bound(currPropToEventIdListProxy.begin(),
                                                                 currPropToEventIdListProxy.end(),
                                                                 R,
                                                                 compPropToEventIdListProxy_);

  double propensityOfChosenEvent = chosenItr->itr->first;
  const std::deque<EventId> & chosenEventIDList = chosenItr->itr->second.eIdDeque;

  // Schulze (Physical Review E, v. 65, 036704) would add a 1 to
  // indForEventList, since he's using 1-based indices instead
  // of the zero-based indices used here.
  std::size_t indForEventList = static_cast<std::size_t>((chosenItr->partialSum - R)/propensityOfChosenEvent);

  // The following condition shouldn't happen since R is
  // supposed to be strictly greater than zero, but might if R
  // is so small that partialSums_[indForPartialSums] - R
  // "equals" partialSums_[indForPartialSums] in floating-point
  // arithmetic.
  //
  if (indForEventList >= static_cast<int>(chosenEventIDList.size())) {
    indForEventList = chosenEventIDList.size() - 1;
  }

  time += -std::log(rng_->getNumInOpenIntervalFrom0To1())/p_s;

  chosenEventID = chosenEventIDList[indForEventList];
}

bool SolverDynamicSchulze::noMoreEvents(int sectNum) const {
  return propToEventIdList_[sectNum].empty();
}

#if KMC_PARALLEL

std::size_t SolverDynamicSchulze::numCellCenteredEvents_(int sectNum) const {

  std::size_t nPossEventsPerSector = 0;
  for (PropToEventIdListProxy_::const_iterator itr = propToEventIdListProxy_[sectNum].begin(),
	 itrEnd = propToEventIdListProxy_[sectNum].end(); itr != itrEnd; ++itr) {
    nPossEventsPerSector += itr->itr->second.eIdDeque.size();
  }

  nPossEventsPerSector -= totNumOverLatticeEvents_[sectNum];

  return nPossEventsPerSector;
}

double SolverDynamicSchulze::totPropensityPerSector_(int sectNum) const {

  const PropToEventIdListProxy_ & currPropToEventIdListProxy = propToEventIdListProxy_[sectNum];

  double p_s = 0;

  // It's faster to iterate over the PropToEventIdListProxy_ container
  // than the PropToEventIdList_ container.
  for (PropToEventIdListProxy_::const_iterator itr = currPropToEventIdListProxy.begin(),
         itrEnd = currPropToEventIdListProxy.end(); itr != itrEnd; ++itr) {
    PropToEventIdList_::iterator mapItr = itr->itr;
    p_s += (itr->recalcPartialSumContrib ? (mapItr->first)*(mapItr->second.eIdDeque.size()) : itr->savedPartialSumContrib);
  }

  return p_s;
}

bool SolverDynamicSchulze::noCellCenteredEvents(int sectNum) const {
  return !(numCellCenteredEvents_(sectNum) > 0);
}

double SolverDynamicSchulze::getLocalMaxAvgPropensityPerPossEvent() const {

  double ps_local_max = 0;

  for (int i = 0; i < propToEventIdListProxy_.size(); ++i) {

    if (propToEventIdListProxy_[i].empty()) {
      continue;
    }

    double p_s = totPropensityPerSector_(i) - totOverLatticePropensity_[i];

    std::size_t nPossEventsPerSector = numCellCenteredEvents_(i);

    if (nPossEventsPerSector > 0) {
      p_s /= nPossEventsPerSector;

      if (p_s > ps_local_max) {
        ps_local_max = p_s;
      }
    }

  }

  return ps_local_max;
}

double SolverDynamicSchulze::getLocalMaxSinglePropensity() const {

  double propensityMaxLocal = 0;

  for (std::size_t i = 0; i < propToEventIdListProxy_.size(); ++i) {
    
    for (PropToEventIdListProxy_::const_iterator itr = propToEventIdListProxy_[i].begin(),
           itrEnd = propToEventIdListProxy_[i].end(); itr != itrEnd; ++itr) {

      PropToEventIdList_::iterator pItr = itr->itr;

      // If all the events in eIdDeque are over-lattice events
      if (pItr->second.eIdDeque.size() ==  pItr->second.numOverLatticeEvents) {
	continue;
      }
      
      double currPropensity = pItr->first;
      
      if (currPropensity > propensityMaxLocal) {
        propensityMaxLocal = currPropensity;
      }

    }
  }

  return propensityMaxLocal;
}

#endif
