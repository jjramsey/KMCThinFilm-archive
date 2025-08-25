#include "SolverBinaryTree.hpp"
#include "Lattice.hpp"

#include <cmath>
#include <algorithm>
#include <boost/array.hpp>

using namespace KMCThinFilm;

SolverBinaryTree::SolverBinaryTree(const Lattice * lattice)
  : 
#if KMC_PARALLEL
  Solver(lattice),
#endif
  events_(lattice->numSectors()),
  treeNodes_(lattice->numSectors()),
  lattice_(lattice)
#if KMC_PARALLEL
  , numOverLatticeEvents_(lattice->numSectors(),0),
  totOverLatticePropensity_(lattice->numSectors(),0)
#endif
{}

void SolverBinaryTree::beginBuildingEventList(int numOverLatticeEvents,
                                              int numReservedLatticePlanes) {

  // WARNING: EventId::dimsForFlattening_ *must* be defined before using this.
  evIdToNodeId_.reset(new EventIdMap<NodeId>(events_.size(),
                                             numOverLatticeEvents,
                                             numReservedLatticePlanes,
                                             -1));

  for (std::size_t i = 0; i < events_.size(); ++i) {
    events_[i].clear();
    treeNodes_[i].clear();

    LatticePlanarBBox sectorBBox;
    lattice_->getSectorPlanarBBox(i, sectorBBox);
    treeNodes_[i].reserve(numReservedLatticePlanes*
                          (sectorBBox.imaxP1 - sectorBBox.imin)*
                          (sectorBBox.jmaxP1 - sectorBBox.jmin)*
                          EventId::dimsForFlattening_[2]);

#if KMC_PARALLEL
    numOverLatticeEvents_[i] = 0;
    totOverLatticePropensity_[i] = 0;
#endif
  }
}

void SolverBinaryTree::appendLeafRaw_(const EventId & eId, double propensity,
				      int sectNum) {

  events_[sectNum].push_back(eId);

  // "Abusing" treeNodes_ as a place for storing propensities for now
  // Later, the vector treeNodes_[sectNum] will be resized and the
  // propensities stored in it moved to the end of the vector.
  treeNodes_[sectNum].push_back(propensity);
}

void SolverBinaryTree::addCellCenteredEntryToEventList(const EventId & eId,
						       double propensity,
						       int sectNum) {
  appendLeafRaw_(eId, propensity, sectNum);
}

void SolverBinaryTree::addOverLatticeEntryToEventList(const EventId & eId,
						      double propensity,
						      int sectNum) {
  appendLeafRaw_(eId, propensity, sectNum);

#if KMC_PARALLEL
  ++(numOverLatticeEvents_[sectNum]);
  totOverLatticePropensity_[sectNum] += propensity;
#endif
}

void SolverBinaryTree::endBuildingEventList() {

  for (int i = 0; i < static_cast<int>(events_.size()); ++i) {
    makeIntNodes_(i);
  }

}

void SolverBinaryTree::addOrUpdateCellCenteredEntryToEventList(const EventId & eId,
							       double currPropensity,
							       int sectNum) {

  NodeId * nodeIdPtr = evIdToNodeId_->getPtrToVal(eId);

  if ((nodeIdPtr == NULL) || (*nodeIdPtr < 0)) {

    /* If it wasn't in the address map before, but now has a
       non-zero propensity, then it should be in the event and
       address maps now. */

    if (currPropensity > 0) {
      appendLeaf_(eId, currPropensity, sectNum);
    }
    
  }
  else {

    NodeId origNodeInd = *nodeIdPtr;
    double & origPropensity = treeNodes_[sectNum][origNodeInd];

    if (currPropensity > 0) {

      if (origPropensity != currPropensity) {
	/* If eId already has existing entries in the event ID lists
           and address map, but its propensity has changed, then those
           entries need to be updated. */

	origPropensity = currPropensity;
	updateAncestorsOfLeafNode_(origNodeInd, sectNum);
      }

    }
    else {

      /* If the current propensity is zero, then eId is associated
         with an event that can't happen.

         Therefore, I remove eId from its current place in the event
         ID list and then its associated entry in the address map.
      */

      std::size_t origLeafInd = origNodeInd - (events_[sectNum].size() - 1);
      EventId & origEventId = events_[sectNum][origLeafInd];

      // "Erasing" the map entry corresponding to the element to be
      // removed by setting it to an invalid (i.e. negative) value.
      *nodeIdPtr = -1;
      
      if ((origLeafInd + 1) != events_[sectNum].size()) {
	// "Removing" the event by replacing it with a copy of the event
	// stored at the last leaf
	origEventId = events_[sectNum].back();
        origPropensity = treeNodes_[sectNum].back();

	evIdToNodeId_->getRefToVal(origEventId) = origNodeInd;
	updateAncestorsOfLeafNode_(origNodeInd, sectNum);
      }

      // Removing the last leaf, which is now redundant with any copy made above.
      
      events_[sectNum].pop_back();
      treeNodes_[sectNum].pop_back();
      
      if (!events_[sectNum].empty()) {

        // At this point, the binary tree is complete but not full,
        // and treeNodes_[sectNum].size() =
        // 2*events_[sectNum].size(). The index of the last non-leaf
        // node, then, is events_[sectNum].size() - 1.
        treeNodes_[sectNum][events_[sectNum].size() - 1] = treeNodes_[sectNum].back();
        treeNodes_[sectNum].pop_back();

        events_[sectNum].push_front(events_[sectNum].back());
        events_[sectNum].pop_back();

        // Now the tree both is full and complete again, and
        // events_[sectNum].size() - 1 is the index of the first leaf.

	evIdToNodeId_->getRefToVal(events_[sectNum].front()) = events_[sectNum].size() - 1;

	updateAncestorsOfLeafNode_(events_[sectNum].size() - 1, sectNum);
      }

    }

  }

}

void SolverBinaryTree::chooseEventIDAndUpdateTime(int sectNum,
                                                  EventId & chosenEventID,
                                                  double & time) {

  double R;
  
  std::size_t chosenChildInd = 0;
  std::size_t numIntNodes = events_[sectNum].size() - 1;

  // The first "if" is factored out from the "while" loop so that R is
  // not calculated unless events_[sectNum].size() > 1.

  if (chosenChildInd < numIntNodes) {      
    R = (treeNodes_[sectNum].front())*(rng_->getNumInOpenIntervalFrom0To1());      
    chooseChildInd_(sectNum, chosenChildInd, R); 
  }

  while (chosenChildInd < numIntNodes) {      
    chooseChildInd_(sectNum, chosenChildInd, R); 
  }

  chosenEventID = events_[sectNum][chosenChildInd - numIntNodes];

  time += -std::log(rng_->getNumInOpenIntervalFrom0To1())/(treeNodes_[sectNum].front());
}

bool SolverBinaryTree::noMoreEvents(int sectNum) const {
  return events_[sectNum].empty();
}

void SolverBinaryTree::appendLeaf_(const EventId & eId,
                                   double propensity, int sectNum) {

  if (!events_[sectNum].empty()) {

    // Note: events_[sectNum].size() - 1 is index of first leaf.
    treeNodes_[sectNum].push_back(treeNodes_[sectNum][events_[sectNum].size() - 1]);

    events_[sectNum].push_back(events_[sectNum].front());
    events_[sectNum].pop_front();

    evIdToNodeId_->getRefToVal(events_[sectNum].back()) = treeNodes_[sectNum].size() - 1;
  }

  evIdToNodeId_->addOrUpdate(eId, treeNodes_[sectNum].size());
  treeNodes_[sectNum].push_back(propensity);
  events_[sectNum].push_back(eId);
   
  // Only need to do the update once, since the two appended leaves
  // are children of the same parent.
  updateAncestorsOfLeafNode_(treeNodes_[sectNum].size() - 1, sectNum);
}

void SolverBinaryTree::makeIntNodes_(int sectNum) {

  std::size_t numPropensities = events_[sectNum].size();
  std::size_t numIntNodes;

  if (numPropensities > 1) {
    
    numIntNodes = numPropensities - 1;

    treeNodes_[sectNum].resize(numIntNodes + numPropensities);  
  
    // Moving the propensity values temporarily stored at the
    // beginning of treeNodes_[sectNum] to the end of it.
    for (std::size_t i = 0; i < numPropensities; ++i) {

      // Have to do the iteration in reverse so that I don't overwrite
      // original propensity values that I haven't used yet.
      std::size_t r = numPropensities - 1 - i;

      treeNodes_[sectNum][r + numIntNodes] = treeNodes_[sectNum][r];
    }

    for (std::size_t i = 0; i < numIntNodes; ++i) {

      // Need to iterate backwards now, since treeNodes_[sectNum][r] needs to have the
      // previous entries of tree properly filled in.
      std::size_t r = numIntNodes - 1 - i;

      std::size_t leftChildInd = 2*r + 1;
      std::size_t rightChildInd = leftChildInd + 1;
      
      treeNodes_[sectNum][r] = treeNodes_[sectNum][leftChildInd] + treeNodes_[sectNum][rightChildInd];
    }
  }
  else {    
    numIntNodes = 0;

    // treeNodes_[sectNum] should already contain the lone propensity
    // value for the one event in the system (if there are any events
    // in the system at all).
  }
  
  for (std::size_t i = 0; i < numPropensities; ++i) {
    evIdToNodeId_->addOrUpdate(events_[sectNum][i], i + numIntNodes);
  }

}

void SolverBinaryTree::updateAncestorsOfLeafNode_(std::size_t nodeInd, int sectNum) {

  std::size_t parentId = nodeInd;

  while (parentId > 0) {
    parentId = (parentId - 1)/2;

    std::size_t leftChildInd = 2*parentId + 1;
    std::size_t rightChildInd = leftChildInd + 1;

    treeNodes_[sectNum][parentId] = treeNodes_[sectNum][leftChildInd] + treeNodes_[sectNum][rightChildInd];
  }

}

#if KMC_PARALLEL

std::size_t SolverBinaryTree::numCellCenteredEvents_(int sectNum) const {
  return (events_[sectNum].size() - numOverLatticeEvents_[sectNum]);
}

bool SolverBinaryTree::noCellCenteredEvents(int sectNum) const {
  return !(numCellCenteredEvents_(sectNum) > 0);
}

double SolverBinaryTree::getLocalMaxAvgPropensityPerPossEvent() const {
  double ps_local_max = 0;
  for (std::size_t i = 0; i < treeNodes_.size(); ++i) {
    std::size_t nPossEventsPerSector = numCellCenteredEvents_(i);
    
    if (nPossEventsPerSector > 0) {

      // If there are no internal nodes for a sector, then there is at
      // most one leaf for that sector, and since nPossEventsPerSector
      // > 0, there is exactly one leaf.
      double p_s = treeNodes_[i].front() - totOverLatticePropensity_[i];

      p_s /= nPossEventsPerSector;

      if (p_s > ps_local_max) {
        ps_local_max = p_s;
      }
    }

  }

  return ps_local_max;
}

double SolverBinaryTree::getLocalMaxSinglePropensity() const {
  double ps_local_max = 0;
  for (std::size_t i = 0; i < events_.size(); ++i) {

    double p_s = 0;

    if (numCellCenteredEvents_(i) > 0) {

      std::size_t nodeIndMin = events_[i].size() - 1;

      for (std::size_t nodeInd = nodeIndMin; 
           nodeInd < treeNodes_[i].size(); ++nodeInd) {
	
	if (events_[i][nodeInd - nodeIndMin].isForOverLattice()) {
	  continue;
	}

	if (treeNodes_[i][nodeInd] > p_s) {
	  p_s = treeNodes_[i][nodeInd];
	}

      }

    }

    if (p_s > ps_local_max) {
      ps_local_max = p_s;
    }
  }

  return ps_local_max;
}

#endif
