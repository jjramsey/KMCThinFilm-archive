#include "CellsToChange.hpp"
#include "Lattice.hpp"
#include "CellInds.hpp"

using namespace KMCThinFilm;

struct CellsToChange::Impl_ {
  Impl_(Lattice * lattice, int numOffsets);

  Lattice * lattice_;

  std::vector<CellIndsOffset> cioVec_;
  std::vector<CellInds> ciVec_;
};

CellsToChange::Impl_::Impl_(Lattice * lattice, int numOffsets)
  : lattice_(lattice) {
  cioVec_.reserve(numOffsets);
}

CellsToChange::CellsToChange(Lattice * lattice, int numOffsets)
  : pImpl_(new Impl_(lattice, numOffsets))
{}

CellsToChange::CellsToChange(const CellsToChange & ctc)
  : pImpl_(new Impl_(*(ctc.pImpl_)))
{}

CellsToChange & CellsToChange::operator=(const CellsToChange & rhs) {

  if (this != &rhs) {
    *pImpl_ = *(rhs.pImpl_);
  }
  
  return *this;
}

void CellsToChange::setCenter(const CellInds & ci) {
  std::vector<CellInds> & ciVec = pImpl_->ciVec_;

  ciVec.clear();
  ciVec.reserve(pImpl_->cioVec_.size());
  
  for (std::vector<CellIndsOffset>::const_iterator itr = pImpl_->cioVec_.begin(),
         itrEnd = pImpl_->cioVec_.end(); itr != itrEnd; ++itr) {
    ciVec.push_back(ci + *itr);
  }

}

void CellsToChange::addLatticePlanes(int numPlanesToAdd) {
  pImpl_->lattice_->addPlanes(numPlanesToAdd);
}

int CellsToChange::getInt(int whichOffset, int whichInt) const {
  return pImpl_->lattice_->getInt(pImpl_->ciVec_[whichOffset], whichInt);
}

double CellsToChange::getFloat(int whichOffset, int whichFloat) const {
  return pImpl_->lattice_->getFloat(pImpl_->ciVec_[whichOffset], whichFloat);
}

void CellsToChange::setInt(int whichOffset, int whichInt, int val) {
  pImpl_->lattice_->setInt(pImpl_->ciVec_[whichOffset], whichInt, val);
}

void CellsToChange::setFloat(int whichOffset, int whichFloat, double val) {
  pImpl_->lattice_->setFloat(pImpl_->ciVec_[whichOffset], whichFloat, val);
}

void CellsToChange::addOffset(const CellIndsOffset & offset) {
  pImpl_->cioVec_.push_back(offset);
}

void CellsToChange::clear() {
  pImpl_->ciVec_.clear();
  pImpl_->cioVec_.clear();
}

const CellInds & CellsToChange::getCenter() const {
  // This assumes the first entry in cioVec_ is (0,0,0), which, given
  // how the offsets are initialized in the Simulation class, is a
  // valid assumption.
  return pImpl_->ciVec_.front();
}

const CellInds & CellsToChange::getCellInds(int whichOffset) const {
  return pImpl_->ciVec_[whichOffset];
}

const std::vector<CellInds> & CellsToChange::getCellIndsVec() const {
  return pImpl_->ciVec_;
}

const std::vector<CellIndsOffset> & CellsToChange::getCellIndsOffsetVec() const {return pImpl_->cioVec_;}

CellsToChange::~CellsToChange() {}
