#include "EventId.hpp"

#include <boost/functional/hash.hpp>

namespace KMCThinFilm {
  
  boost::array<int,3> EventId::dimsForFlattening_ = {{0,0,0}};
  boost::array<int,2> EventId::ciMin_ = {{0,0}};

  std::size_t hash_value(const EventId & eId) {
    std::size_t seed = 0;
    boost::hash_combine(seed, eId.e1_);
    boost::hash_combine(seed, eId.e2_);
    return seed;
  }

}

using namespace KMCThinFilm;

void EventId::getEventInfo(CellInds & ci, int & cellCenEventIndex) const {
  boost::array<int,3> r;

  // Reversing the column-major style flattening. Note deliberate use
  // of the truncation feature of integer division.
  r[0] = e1_/dimsForFlattening_[0];
  r[1] = r[0]/dimsForFlattening_[1];
  r[2] = r[1]/dimsForFlattening_[2];

  ci.i = e1_ - dimsForFlattening_[0]*r[0] + ciMin_[0];
  ci.j = r[0] - dimsForFlattening_[1]*r[1] + ciMin_[1];
  cellCenEventIndex =  r[1] - dimsForFlattening_[2]*r[2];
        
  ci.k = e2_;
}


std::string EventId::toString() const {
  if (isForOverLattice()) {
    int overLatticeEventIndex, sectNum;    
    getEventInfo(overLatticeEventIndex, sectNum);

    return "OverLatticeEvent(Sector=" + boost::lexical_cast<std::string>(sectNum) +
      "; Event index=" + boost::lexical_cast<std::string>(overLatticeEventIndex) + ")";
  }
  else {
    CellInds ci;
    int cellCenEventIndex;
    getEventInfo(ci, cellCenEventIndex);

    return "CellCenteredEvent(Cell indices=[" +
      boost::lexical_cast<std::string>(ci.i) + "," +
      boost::lexical_cast<std::string>(ci.j) + "," +
      boost::lexical_cast<std::string>(ci.k) +
      "]; Event index=" + boost::lexical_cast<std::string>(cellCenEventIndex) + ")";
  }
}
