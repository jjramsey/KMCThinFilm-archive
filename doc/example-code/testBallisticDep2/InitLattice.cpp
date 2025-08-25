#include "InitLattice.hpp"
#include "EventsAndActions.hpp"

using namespace KMCThinFilm;

void SetEmptyCellWithRandColor::operator()(const CellInds & ci,
                                                  const Lattice & lattice,
                                                  std::vector<int> & emptyIntVals,
                                                  std::vector<double> & emptyFloatVals) {
  
  emptyFloatVals[BDFloatVal::COLOR] = rng_->getNumInOpenIntervalFrom0To1();
  
}
