#ifndef EVENT_ID_MAP_HPP
#define EVENT_ID_MAP_HPP

#include "EventId.hpp"

#include <vector>

namespace KMCThinFilm {
  
  template<typename T>
  class EventIdMap {
  public:

    // WARNING: EventId::dimsForFlattening_ *must* be defined before using this
    EventIdMap(int numSectors, 
               int numOverLatticeEvents,
               int numReservedLatticePlanes,
               const T & defaultVal)
      : defaultVal_(defaultVal) {

      overLatticeEIdMap_.clear();
      overLatticeEIdMap_.reserve(numSectors);

      for (int i = 0; i < numSectors; ++i) {
        overLatticeEIdMap_.push_back(std::vector<T>(numOverLatticeEvents, defaultVal_));
      }

      cellCenteredEIdMap_.clear();
      cellCenteredEIdMap_.reserve(numReservedLatticePlanes);
      cellCenteredEIdMapSize_ = EventId::dimsForFlattening_[0]*EventId::dimsForFlattening_[1]*EventId::dimsForFlattening_[2];
    }

    // This checks if there is a value corresponding to eId. If not it
    // returns NULL. Note that a call to addOrUpdate could invalidate
    // the returned pointer.
    T* getPtrToVal(const EventId & eId) {
      T* outPtr = NULL;

      if (eId.isForOverLattice()) {
        int overLatticeEventIndex, sectNum;
        eId.getEventInfo(overLatticeEventIndex, sectNum);
        
        outPtr = &(overLatticeEIdMap_[sectNum][overLatticeEventIndex]);
      }
      else {
        if ((eId.e2_ < static_cast<int>(cellCenteredEIdMap_.size())) && (eId.e1_ < cellCenteredEIdMapSize_)) {
          outPtr = &(cellCenteredEIdMap_[eId.e2_][eId.e1_]);
        }
      }

      return outPtr;
    }

    // This gets a non-constant reference to a value in the map,
    // presuming that such a reference *exists*.
    T & getRefToVal(const EventId & eId) {

      if (eId.isForOverLattice()) {
        int overLatticeEventIndex, sectNum;
        eId.getEventInfo(overLatticeEventIndex, sectNum);
        
        return overLatticeEIdMap_[sectNum][overLatticeEventIndex];
      }
      else {
        return cellCenteredEIdMap_[eId.e2_][eId.e1_];
      }

    }

    void addOrUpdate(const EventId & eId, const T & val) {
      
      if (eId.isForOverLattice()) {
        int overLatticeEventIndex, sectNum;
        eId.getEventInfo(overLatticeEventIndex, sectNum);
        overLatticeEIdMap_[sectNum][overLatticeEventIndex] = val;
      }
      else {
        
        // This checks if space for eId is already available, and if not,
        // adds it.
        while (eId.e2_ >= static_cast<int>(cellCenteredEIdMap_.size())) {
          cellCenteredEIdMap_.push_back(std::vector<T>(cellCenteredEIdMapSize_, defaultVal_));
        }

        cellCenteredEIdMap_[eId.e2_][eId.e1_] = val;
      }

    }

  private:
    std::vector<std::vector<T> > cellCenteredEIdMap_;
    std::vector<std::vector<T> > overLatticeEIdMap_;
    int cellCenteredEIdMapSize_;
    T defaultVal_;
  };

}

#endif /* EVENT_ID_MAP_HPP */
