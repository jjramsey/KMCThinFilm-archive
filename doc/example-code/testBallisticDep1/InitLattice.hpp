#ifndef INIT_LATTICE_HPP
#define INIT_LATTICE_HPP

#include <KMCThinFilm/Lattice.hpp>
#include <KMCThinFilm/RandNumGen.hpp>

class SetEmptyCellWithRandColor {
public:
  SetEmptyCellWithRandColor(KMCThinFilm::RandNumGenSharedPtr rng)
    : rng_(rng)
  {}

  void operator()(const KMCThinFilm::CellInds & ci,
                  const KMCThinFilm::Lattice & lattice,
                  std::vector<int> & emptyIntVals,
                  std::vector<double> & emptyFloatVals);
private:
  KMCThinFilm::RandNumGenSharedPtr rng_;
};

#endif /* INIT_LATTICE_HPP */
