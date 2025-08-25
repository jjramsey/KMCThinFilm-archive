#ifndef SOLVER_BINARY_TREE_HPP
#define SOLVER_BINARY_TREE_HPP

#include "Solver.hpp"
#include "EventIdMap.hpp"

#include <vector>
#include <deque>

#include <boost/scoped_ptr.hpp>

namespace KMCThinFilm {

  class Lattice;

  class SolverBinaryTree : public Solver {
  public:
    SolverBinaryTree(const Lattice * lattice);

    virtual void beginBuildingEventList(int numOverLatticeEvents,
                                        int numReservedLatticePlanes);

    virtual void addCellCenteredEntryToEventList(const EventId & eId,
						 double propensity,
						 int sectNum);
    
    virtual void addOverLatticeEntryToEventList(const EventId & eId,
						double propensity,
						int sectNum);
    
    virtual void endBuildingEventList();

    virtual void addOrUpdateCellCenteredEntryToEventList(const EventId & eId,
							 double currPropensity,
							 int sectNum);

    virtual void chooseEventIDAndUpdateTime(int sectNum,
					    EventId & chosenEventID,
					    double & time);

    virtual bool noMoreEvents(int sectNum) const;

#if KMC_PARALLEL
    virtual bool noCellCenteredEvents(int sectNum) const;
#endif

  private:

#if KMC_PARALLEL
    virtual double getLocalMaxAvgPropensityPerPossEvent() const;
    virtual double getLocalMaxSinglePropensity() const;

    std::vector<std::size_t> numOverLatticeEvents_;
    std::vector<double> totOverLatticePropensity_;
#endif

    std::vector<std::deque<EventId> > events_;
    std::vector<std::vector<double> > treeNodes_;
    const Lattice * lattice_;

    // Allows a node ID to be negative in order to indicate an invalid ID.
    typedef std::ptrdiff_t NodeId;

    boost::scoped_ptr<EventIdMap<NodeId> > evIdToNodeId_;

    void makeIntNodes_(int sectNum);

    void updateAncestorsOfLeafNode_(std::size_t nodeInd, int sectNum);

    void appendLeafRaw_(const EventId & eId, double propensity,
			int sectNum);

    void appendLeaf_(const EventId & eId, double propensity,
		     int sectNum);

#if KMC_PARALLEL
    std::size_t numCellCenteredEvents_(int sectNum) const;
#endif

    inline void chooseChildInd_(int sectNum,
                                std::size_t & chosenChildInd,
                                double & R) {

      std::size_t leftChildInd = 2*chosenChildInd + 1;
      
      double leftChildVal = treeNodes_[sectNum][leftChildInd];

      if (R <= leftChildVal) {
        chosenChildInd = leftChildInd;
      }
      else {
        // chosenChildInd is the index of the right child.
        chosenChildInd = leftChildInd + 1;
        R -= leftChildVal;
      }
    }

  };

}


#endif /* SOLVER_BINARY_TREE_HPP */ 
