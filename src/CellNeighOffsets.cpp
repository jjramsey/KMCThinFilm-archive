#include "CellNeighOffsets.hpp"
#include "ErrorHandling.hpp"

#include <vector>
#include <boost/lexical_cast.hpp>

using namespace KMCThinFilm;

struct CellNeighOffsets::Impl_ {
  Impl_(int numberOfCells);
  std::vector<CellIndsOffset> offsets_;

  void resetOffsets_(int numberOfCells);
};

void CellNeighOffsets::Impl_::resetOffsets_(int numberOfCells) {
  exitOnCondition(numberOfCells < 1,
		  "There should be at least one offset.");

  offsets_.clear();
  offsets_.resize(numberOfCells);
  offsets_[0] = CellIndsOffset(0,0,0);
}

CellNeighOffsets::Impl_::Impl_(int numberOfCells) {
  resetOffsets_(numberOfCells);
}

CellNeighOffsets::CellNeighOffsets(int numberOfCells) 
  : pImpl_(new Impl_(numberOfCells))
{}

CellNeighOffsets::CellNeighOffsets(const CellNeighOffsets & cno)
  : pImpl_(new Impl_(*(cno.pImpl_)))
{}

CellNeighOffsets & CellNeighOffsets::operator=(const CellNeighOffsets & rhs) {

  if (this != &rhs) {
    *pImpl_ = *(rhs.pImpl_);
  }
  
  return *this;
}

void CellNeighOffsets::addOffset(int whichOffset, 
			     const CellIndsOffset & offset) {

  exitOnCondition(whichOffset == 0,
		  "Offset with index 0 is reserved for the offset with components (0,0,0).");

  pImpl_->offsets_[whichOffset] = offset;
}

void CellNeighOffsets::resetOffsets(int numberOfNewOffsets) {
  pImpl_->resetOffsets_(numberOfNewOffsets);
}

void CellNeighOffsets::clearOffsets() {
  pImpl_->resetOffsets_(1);
}

int CellNeighOffsets::numOffsets() const {
  return pImpl_->offsets_.size();
}

const CellIndsOffset & CellNeighOffsets::getOffset(int whichOffset) const {

  exitOnCondition((whichOffset < 0) || (whichOffset >= numOffsets()),
		  "Offset index " + boost::lexical_cast<std::string>(whichOffset) +  " is out of bounds.");

  return pImpl_->offsets_[whichOffset];
}

CellNeighOffsets::~CellNeighOffsets() {}
