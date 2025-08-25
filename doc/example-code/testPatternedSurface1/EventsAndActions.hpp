#ifndef EVENTS_AND_ACTIONS_HPP
#define EVENTS_AND_ACTIONS_HPP

/* 
   Comments in this code of the form "//! [...]" are used to assist
   Doxygen in documenting this file.
*/

#include <string>

#include <KMCThinFilm/CellCenteredGroupPropensities.hpp>
#include <KMCThinFilm/EventExecutor.hpp>
#include <KMCThinFilm/MakeEnum.hpp>

#include "PhysicalConstants.hpp"

//! [event and action enums]
KMC_MAKE_ID_ENUM(PSOverLatticeEvents,
		 DEPOSITION);

KMC_MAKE_ID_ENUM(PSCellCenteredEvents,
		 HOP_UP,
		 HOP_DOWN,
		 HOP_LEFT,
		 HOP_RIGHT);

KMC_MAKE_ID_ENUM(PAction,
		 PRINT);
//! [event and action enums]

//! [lattice enum]
KMC_MAKE_LATTICE_INTVAL_ENUM(PS, HEIGHT);
KMC_MAKE_LATTICE_FLOATVAL_ENUM(PS, E_s);
//! [lattice enum]

//! [offset enum]
KMC_MAKE_OFFSET_ENUM(HopOffset,
		     UP, DOWN, LEFT, RIGHT);
//! [offset enum]

//! [dep execute]
void DepositionExecute(const KMCThinFilm::CellInds & ci,
		       const KMCThinFilm::SimulationState & simState,
		       KMCThinFilm::Lattice & lattice);
//! [dep execute]

//! [hop prop]
class HoppingPropensity {
public:
  HoppingPropensity(double E_n, double T) 
    : E_n_(E_n), kBT_(PhysConst::kB*T), k_(kBT_/PhysConst::h)
  {}
  void operator()(const KMCThinFilm::CellNeighProbe & cnp,
                  std::vector<double> & propensityVec) const;
private:
  double E_n_, kBT_, k_;
};
//! [hop prop]

//! [hop exec]
class HoppingExecute {
public:
  void operator()(const KMCThinFilm::CellInds & ci,
		  const KMCThinFilm::SimulationState & simState,
                  const KMCThinFilm::Lattice & lattice,
                  std::vector<KMCThinFilm::CellsToChange> & ctcVec) const;
};
//! [hop exec]

//! [print op]
class PrintASCII {
public:
  PrintASCII(const std::string & fNameRoot)
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
