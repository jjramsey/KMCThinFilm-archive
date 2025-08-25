#ifndef SOLVER_DYNAMIC_SCHULZE_HPP
#define SOLVER_DYNAMIC_SCHULZE_HPP

#include "Solver.hpp"
#include "EventIdMap.hpp"

#include <vector>
#include <deque>
#include <map>
#include <cstddef>

#include <boost/scoped_ptr.hpp>

namespace KMCThinFilm {

  class Lattice;

  class SolverDynamicSchulze : public Solver {
  public:
    SolverDynamicSchulze(const Lattice * lattice);

    virtual void beginBuildingEventList(int numOverLatticeEvents,
                                        int numReservedLatticePlanes);

    virtual void addCellCenteredEntryToEventList(const EventId & eId,
						 double propensity,
						 int sectNum);

    virtual void addOverLatticeEntryToEventList(const EventId & eId,
						double propensity,
						int sectNum);

    virtual void endBuildingEventList() {}

    virtual void addOrUpdateCellCenteredEntryToEventList(const EventId & eId,
							 double propensity,
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

    std::vector<double> totOverLatticePropensity_;
    std::vector<std::size_t> totNumOverLatticeEvents_;
#endif

    struct IndexDequePair_ {
      std::size_t indexToEIdListProxy;
      std::deque<EventId> eIdDeque;

#if KMC_PARALLEL
      std::size_t numOverLatticeEvents;
#endif

      IndexDequePair_(std::size_t ind)
        : indexToEIdListProxy(ind)
#if KMC_PARALLEL
	, numOverLatticeEvents(0)
#endif
      {}
    };

    // Using a std::map rather than a hash map, since if items are
    // inserted during an update, iterators are not invalidated, and I'm
    // storing iterators of a PropToEventIdList_ in addrMap_.
    typedef std::map<double,IndexDequePair_> PropToEventIdList_;
  
    // One instance of PropToEventIdList_ per sector
    std::vector<PropToEventIdList_> propToEventIdList_;

    // This points to an entry in propToEventIdList_ and contains some other stuff.
    struct PropToEventIdListProxyEntry_ {
      PropToEventIdList_::iterator itr;
      
      bool recalcPartialSumContrib;
      double savedPartialSumContrib;
      double partialSum;

      PropToEventIdListProxyEntry_(PropToEventIdList_::iterator myItr)
        : itr(myItr)
      {}
    };

    // This is used (in part) to access propToEventIdList_ indirectly,
    // for use in cases where direct use of propToEventIdList_ would
    // be slower.
    typedef std::deque<PropToEventIdListProxyEntry_> PropToEventIdListProxy_;

    // Used as a comparator for the std::lower_bound function
    struct CompPropToEventIdListProxy_ {
      bool operator()(const PropToEventIdListProxyEntry_ & entry,
                      double R) const {
        return (entry.partialSum < R);
      }
    } compPropToEventIdListProxy_;

    // One instance of PropToEventIdListProxy_ per sector
    std::vector<PropToEventIdListProxy_> propToEventIdListProxy_;

    // Allows a index to be negative in order to indicate that the
    // corresponding event ID is invalid.
    typedef std::ptrdiff_t Index;

    struct EvListItrIndexPair_ {
      PropToEventIdList_::iterator itr;
      Index indexToEId;

      EvListItrIndexPair_(PropToEventIdList_::iterator myItr, Index i)
	: itr(myItr), indexToEId(i)
      {}

    };

    boost::scoped_ptr<EventIdMap<EvListItrIndexPair_> > addrMap_;

    PropToEventIdList_::iterator getPropToEventIdListItr_(double propensity,
                                                          int sectNum);

    void addItrEntryToEventList_(const EventId & eId,
				 PropToEventIdList_::iterator propToEventIdListItr,
				 int sectNum);

    void removeFromEventIdList_(std::size_t origInd, std::deque<EventId> & origDeque);

    double totPropensityPerSector_(int sectNum) const;

#if KMC_PARALLEL
    std::size_t numCellCenteredEvents_(int sectNum) const;
#endif

  };

}

#endif /* SOLVER_DYNAMIC_SCHULZE_HPP */
