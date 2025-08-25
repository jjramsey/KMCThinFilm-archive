#ifndef EVENT_ID_HPP
#define EVENT_ID_HPP

#include "CellInds.hpp"

#include <boost/array.hpp>
#include <boost/lexical_cast.hpp>

namespace KMCThinFilm {

  struct EventId {

    EventId()
      : e1_(0), e2_(0)
    {}

    // This constructor should not be called until dimsForFlattening_ array is defined.
    EventId(const CellInds & ci, int cellCenEventIndex)
      : e1_((ci.i - ciMin_[0]) + dimsForFlattening_[0]*((ci.j - ciMin_[1]) + dimsForFlattening_[1]*cellCenEventIndex)),
	e2_(ci.k)
    {}

    EventId(int overLatticeEventIndex, int sectNum)
      : e1_(overLatticeEventIndex), e2_(-(sectNum + 1))
    {}

    bool isForOverLattice() const {return e2_ < 0;}

    void getEventInfo(int & overLatticeEventIndex, int & sectNum) const {
      overLatticeEventIndex = e1_;
      sectNum = -(e2_ + 1);
    }

    // For hash map
    bool operator==(const EventId & rhs) const {
      return (e1_ == rhs.e1_) && (e2_ == rhs.e2_);
    }

    // For std::map
    bool operator<(const EventId & rhs) const {
      return (e1_ < rhs.e1_) || ((e1_ == rhs.e1_) && (e2_ < rhs.e2_));
    }

    std::string toString() const;

    void getEventInfo(CellInds & ci, int & cellCenEventIndex) const;

    int e1_, e2_;
    static boost::array<int,3> dimsForFlattening_;
    static boost::array<int,2> ciMin_;
  };

  std::size_t hash_value(const EventId & eId);

}

#endif /* EVENT_ID_HPP */
