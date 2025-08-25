#ifndef EVENTS_AND_ACTIONS_HPP
#define EVENTS_AND_ACTIONS_HPP

#include <KMCThinFilm/CellNeighOffsets.hpp>
#include <KMCThinFilm/CellCenteredGroupPropensities.hpp>
#include <KMCThinFilm/EventExecutor.hpp>
#include <KMCThinFilm/MakeEnum.hpp>
#include <KMCThinFilm/RandNumGen.hpp>

#include <vector>
#include <string>

KMC_MAKE_ID_ENUM(OverLatticeEvents,
                 DEPOSITION);

KMC_MAKE_ID_ENUM(CellCenteredEvents,
                 COLOR_MIXING);

KMC_MAKE_ID_ENUM(PAction,
                 PRINT);

//! [lat enum]
KMC_MAKE_LATTICE_INTVAL_ENUM(BD, IS_OCCUPIED, ACTIVE_ZONE_HEIGHT);

KMC_MAKE_LATTICE_FLOATVAL_ENUM(BD, COLOR);
//! [lat enum]

//! [offset enum]
KMC_MAKE_OFFSET_ENUM(MIX_OFFSET,
                     NORTH, SOUTH, WEST, EAST /* First four neighbors are lateral */, 
                     UP, DOWN);
//! [offset enum]

//! [dep exec]
class DepositionExecute {
public:
  DepositionExecute();

  void operator()(const KMCThinFilm::CellInds & ci,
		  const KMCThinFilm::SimulationState & simState,
		  KMCThinFilm::Lattice & lattice);
private:
  std::vector<KMCThinFilm::CellIndsOffset> neighOffsets_;
};
//! [dep exec]

//! [mix prop]
class ColorMixPropensity {
public:
  ColorMixPropensity(double mixPropPerNeighbor)
    : mixPropPerNeighbor_(mixPropPerNeighbor)
  {} 

  void operator()(const KMCThinFilm::CellNeighProbe & cnp,
                  std::vector<double> & propensityVec) const;
private:
  double mixPropPerNeighbor_;
};
//! [mix prop]

//! [mix exec]
class ColorMixExecute {
public:
  ColorMixExecute(KMCThinFilm::RandNumGenSharedPtr rng_,
                  const KMCThinFilm::CellNeighOffsets * mixCNO,
                  int * numMixes);

  void operator()(const KMCThinFilm::CellInds & ci,
                  const KMCThinFilm::SimulationState & simState,
                  KMCThinFilm::Lattice & lattice);
private:
  KMCThinFilm::RandNumGenSharedPtr rng_;
  const KMCThinFilm::CellNeighOffsets * mixCNO_;
  int * numMixes_;
};
//! [mix exec]

//! [print op]
class PrintPoint3D {
public:
  PrintPoint3D(const std::string & fNameRoot)
    : fNameRoot_(fNameRoot),
      snapShotCntr_(0)
  {}

  void operator()(const KMCThinFilm::SimulationState & simState,
                  KMCThinFilm::Lattice & lattice);

private:
  std::string fNameRoot_;
  int snapShotCntr_;
};
//! [print op]

#endif /* EVENTS_AND_ACTIONS_HPP */
