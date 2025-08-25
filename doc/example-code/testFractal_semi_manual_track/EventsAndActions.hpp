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

//! [event and action enums]
KMC_MAKE_ID_ENUM(FOverLatticeEvents,
		 DEPOSITION);

KMC_MAKE_ID_ENUM(FCellCenteredEvents,
		 HOP_UP,
		 HOP_DOWN,
		 HOP_LEFT,
		 HOP_RIGHT);

KMC_MAKE_ID_ENUM(PAction,
		 PRINT);
//! [event and action enums]

//! [lattice enum]
KMC_MAKE_LATTICE_INTVAL_ENUM(F, HEIGHT);
//! [lattice enum]

/* I have the offsets for hopping include both nearest and
   next-nearest neighbors to make it more obvious what the difference
   is between the FCellCenteredEvents enumeration and the HopOffset
   enumeration. */

//! [offset enum]
KMC_MAKE_OFFSET_ENUM(HopOffset,
		     UP, DOWN, LEFT, RIGHT,
                     RIGHT_UP, RIGHT_DOWN, LEFT_UP, LEFT_DOWN);
//! [offset enum]

//! [dep execute]
void DepositionExecute(const KMCThinFilm::CellInds & ci,
		       const KMCThinFilm::SimulationState & simState,
		       const KMCThinFilm::Lattice & lattice,
                       std::vector<KMCThinFilm::CellsToChange> & ctcVec);
//! [dep execute]

//! [hop prop]
class HoppingPropensity {
public:
  HoppingPropensity(double D) : D_(D) {}
  void operator()(const KMCThinFilm::CellNeighProbe & cnp, 
                  std::vector<double> & propensityVec) const;
private:
  double D_;
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
