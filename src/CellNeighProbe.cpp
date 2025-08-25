#include "CellNeighProbe.hpp"
#include "Lattice.hpp"
#include "ErrorHandling.hpp"

#include <cassert>

using namespace KMCThinFilm;

struct CellNeighProbe::Impl_ {
  Impl_(const Lattice * lattice);

  const Lattice * lattice_;
  const CellInds * ci_;
  const std::vector<CellIndsOffset> * cioVecPtr_;
};

CellNeighProbe::Impl_::Impl_(const Lattice * lattice)
  : lattice_(lattice),
    ci_(NULL),
    cioVecPtr_(NULL)
{}

CellNeighProbe::CellNeighProbe(const Lattice * lattice)
  : pImpl_(new Impl_(lattice))
{}

CellNeighProbe::CellNeighProbe(const CellNeighProbe & ctp)
  : pImpl_(new Impl_(*(ctp.pImpl_)))
{}

CellNeighProbe & CellNeighProbe::operator=(const CellNeighProbe & rhs) {

  if (this != &rhs) {
    *pImpl_ = *(rhs.pImpl_);
  }
  
  return *this;
}

void CellNeighProbe::attachLattice(const Lattice * lattice) {
  pImpl_->lattice_ = lattice;
}

void CellNeighProbe::attachCellInds(const CellInds * ci, 
				    const std::vector<CellIndsOffset> * cioVecPtr) {
  assert(ci != NULL);
  assert(cioVecPtr != NULL);

  pImpl_->ci_ = ci;
  pImpl_->cioVecPtr_ = cioVecPtr;
}

CellToProbe CellNeighProbe::getCellToProbe(int probedCellInd) const {
  return CellToProbe(*(pImpl_->ci_) + (*(pImpl_->cioVecPtr_))[probedCellInd]);
}

double CellNeighProbe::getFloat(const CellToProbe & ctp, int whichFloat) const {
  return pImpl_->lattice_->getFloat(ctp.ci_, whichFloat);
}

int CellNeighProbe::getInt(const CellToProbe & ctp, int whichInt) const {
  return pImpl_->lattice_->getInt(ctp.ci_, whichInt);
}

bool CellNeighProbe::exceedsLatticeHeight(const CellToProbe & ctp) const {
  return ctp.ci_.k >= pImpl_->lattice_->currHeight();
}

bool CellNeighProbe::belowLatticeBottom(const CellToProbe & ctp) const {
  return ctp.ci_.k < 0;
}

CellNeighProbe::~CellNeighProbe() {}
