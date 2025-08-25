#include "EventExecutorGroup.hpp"
#include "ErrorHandling.hpp"

#include <vector>

#include <boost/lexical_cast.hpp>
#include <boost/variant.hpp>

using namespace KMCThinFilm;

struct EventExecutorGroup::Impl_ {
  Impl_(int numEventsInGroup);

  struct EventExecutorSemiManualTrackInfo_ {
    EventExecutorSemiManualTrack evExec_;
    std::vector<CellNeighOffsets> cnoVec_;

    EventExecutorSemiManualTrackInfo_(EventExecutorSemiManualTrack evExec,
                                      const std::vector<CellNeighOffsets> & cnoVec)
      : evExec_(evExec), cnoVec_(cnoVec)
    {}

  };

  // Note: The implementation of the function getEventExecutorType()
  // assumes that the types in this typedef are in a specific
  // order. It's safe to append a new type, but not to switch types
  // around.
  typedef boost::variant<EventExecutorAutoTrack,
                         EventExecutorSemiManualTrackInfo_> EventExecutorInfo_;

  std::vector<EventExecutorInfo_> eventsInfo_;

  void resetGroup_(int numEventsInGroup);
};

void EventExecutorGroup::Impl_::resetGroup_(int numEventsInGroup) {
  eventsInfo_.clear();
  eventsInfo_.resize(numEventsInGroup);
}


EventExecutorGroup::Impl_::Impl_(int numEventsInGroup) {
  resetGroup_(numEventsInGroup);
}

EventExecutorGroup::EventExecutorGroup(int numEventsInGroup) 
  : pImpl_(new Impl_(numEventsInGroup))
{}

EventExecutorGroup::EventExecutorGroup(const EventExecutorGroup & eeg)
  : pImpl_(new Impl_(*(eeg.pImpl_)))
{}

EventExecutorGroup & EventExecutorGroup::operator=(const EventExecutorGroup & rhs) {
  
  if (this != &rhs) {
    *pImpl_ = *(rhs.pImpl_);
  }
  
  return *this;
}

void EventExecutorGroup::addEventExecutor(int whichEvent,
                                          EventExecutorAutoTrack evExec) {
  
  exitOnCondition((whichEvent < 0) || (whichEvent >= numEventExecutors()),
                  "Event index " + boost::lexical_cast<std::string>(whichEvent) +  " is out of bounds.");

  pImpl_->eventsInfo_[whichEvent] = evExec;
}

void EventExecutorGroup::addEventExecutor(int whichEvent,
                                          EventExecutorSemiManualTrack evExec,
                                          const std::vector<CellNeighOffsets> & cnoVec) {
  
  exitOnCondition((whichEvent < 0) || (whichEvent >= numEventExecutors()),
                  "Event index " + boost::lexical_cast<std::string>(whichEvent) +  " is out of bounds.");

  pImpl_->eventsInfo_[whichEvent] = Impl_::EventExecutorSemiManualTrackInfo_(evExec, cnoVec);
}


void EventExecutorGroup::resetGroup(int numEventsInGroup) {
  pImpl_->resetGroup_(numEventsInGroup);
}

void EventExecutorGroup::clearGroup() {
  pImpl_->eventsInfo_.clear();
}

int EventExecutorGroup::numEventExecutors() const {
  return pImpl_->eventsInfo_.size();
}

EventExecutorGroup::EventExecEnum::Type EventExecutorGroup::getEventExecutorType(int whichEvent) const {

  exitOnCondition((whichEvent < 0) || (whichEvent >= numEventExecutors()),
                  "Event index " + boost::lexical_cast<std::string>(whichEvent) +  " is out of bounds.");
  
  EventExecEnum::Type evExecType;

  // Note: The following "switch" statement assumes that the types
  // within the EventExecutorInfo_ variant typedef are in a specific
  // order.
  switch (pImpl_->eventsInfo_[whichEvent].which()) {
  case 0:
    evExecType = EventExecEnum::AUTO;
    break;
  case 1:
    evExecType = EventExecEnum::SEMIMANUAL;
    break;
  default:
    exitWithMsg("Unknown EventExecutor type!");
  }

  return evExecType;
}

EventExecutorAutoTrack EventExecutorGroup::getEventExecutorAutoTrack(int whichEvent) const {
  return boost::get<EventExecutorAutoTrack>(pImpl_->eventsInfo_[whichEvent]);
}

EventExecutorSemiManualTrack EventExecutorGroup::getEventExecutorSemiManualTrack(int whichEvent) const {
  return boost::get<Impl_::EventExecutorSemiManualTrackInfo_>(pImpl_->eventsInfo_[whichEvent]).evExec_;
}

const std::vector<CellNeighOffsets> & EventExecutorGroup::getEventExecutorOffsetsVec(int whichEvent) const {
  return boost::get<Impl_::EventExecutorSemiManualTrackInfo_>(pImpl_->eventsInfo_[whichEvent]).cnoVec_;
}

EventExecutorGroup::~EventExecutorGroup() {}
