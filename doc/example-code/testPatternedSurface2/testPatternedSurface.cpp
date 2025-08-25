/* 
   Any comments in this code of the form "//! [...]" are used to
   assist Doxygen in documenting this file.  */

//! [headers]
#include <string>
#include <fstream>

#include <boost/array.hpp>

#include <KMCThinFilm/Simulation.hpp>
#include <KMCThinFilm/RandNumGenMT19937.hpp>
#include <KMCThinFilm/ErrorHandling.hpp>

#include "EventsAndActions.hpp"
//! [headers]

//! [using decl]
using namespace KMCThinFilm;
//! [using decl]

int main() {

  //! [hardcoded parameters]
  double F = 0.0033, maxCoverage = 0.15, E_n = 0.18, T = 390;
  int numDomains = 16;
  unsigned int seed = 42;
  SolverId::Type sId = SolverId::DYNAMIC_SCHULZE;
  //! [hardcoded parameters]

  //! [reading in E_s]
  DblArray2D E_s;
  
  std::ifstream patternFile("../singleDomain.dat");
  boost::array<int,2> E_s_extents;

  patternFile >> E_s_extents[0] >> E_s_extents[1];
  E_s.resize(E_s_extents);

  int E_s_i, E_s_j;
  double E_s_val;

  while (patternFile >> E_s_i >> E_s_j >> E_s_val) {
    E_s[E_s_i][E_s_j] = E_s_val;
  }
  //! [reading in E_s]
  
  //! [deposition time]
  double approxDepTime = maxCoverage/F;
  //! [deposition time]

  //! [initializing simulation]
  LatticeParams latParams;
  latParams.numIntsPerCell = PSIntVal::SIZE;
  latParams.numFloatsPerCell = PSFloatVal::SIZE;
  latParams.globalPlanarDims[0] = numDomains*(E_s.shape()[0]);
  latParams.globalPlanarDims[1] = numDomains*(E_s.shape()[1]);

  Simulation sim(latParams);
  //! [initializing simulation]

  //! [setting solver and rng]
  sim.setSolver(sId);

  RandNumGenSharedPtr rng(new RandNumGenMT19937(seed));

  sim.setRNG(rng);
  //! [setting solver and rng]

  //! [adding overlattice events]
  sim.reserveOverLatticeEvents(PSOverLatticeEvents::SIZE);
  sim.addOverLatticeEvent(PSOverLatticeEvents::DEPOSITION,
			  F, DepositionExecute);
  //! [adding overlattice events]
  
  //! [making cellneighoffsets]
  CellNeighOffsets hopCNO(HopOffset::SIZE);

  hopCNO.addOffset(HopOffset::UP,    CellIndsOffset(0,+1,0));
  hopCNO.addOffset(HopOffset::DOWN,  CellIndsOffset(0,-1,0));
  hopCNO.addOffset(HopOffset::LEFT,  CellIndsOffset(-1,0,0));
  hopCNO.addOffset(HopOffset::RIGHT, CellIndsOffset(+1,0,0));
  //! [making cellneighoffsets]

  //! [adding cellcentered events]
  sim.reserveCellCenteredEventGroups(1, PSCellCenteredEvents::SIZE);
  
  EventExecutorGroup hopExecs(PSCellCenteredEvents::SIZE);

  std::vector<CellNeighOffsets> tmpExecCNO;
  tmpExecCNO.reserve(1);
  tmpExecCNO.push_back(CellNeighOffsets(2));
  tmpExecCNO.back().addOffset(1, hopCNO.getOffset(HopOffset::LEFT));
  
  hopExecs.addEventExecutor(PSCellCenteredEvents::HOP_LEFT,
                            HoppingExecute(), tmpExecCNO);

  tmpExecCNO.back().resetOffsets(2);
  tmpExecCNO.back().addOffset(1, hopCNO.getOffset(HopOffset::RIGHT));
  hopExecs.addEventExecutor(PSCellCenteredEvents::HOP_RIGHT,
                            HoppingExecute(), tmpExecCNO);

  tmpExecCNO.back().resetOffsets(2);
  tmpExecCNO.back().addOffset(1, hopCNO.getOffset(HopOffset::UP));
  hopExecs.addEventExecutor(PSCellCenteredEvents::HOP_UP,
                            HoppingExecute(), tmpExecCNO);

  tmpExecCNO.back().resetOffsets(2);
  tmpExecCNO.back().addOffset(1, hopCNO.getOffset(HopOffset::DOWN));
  hopExecs.addEventExecutor(PSCellCenteredEvents::HOP_DOWN,
                            HoppingExecute(), tmpExecCNO);

  sim.addCellCenteredEventGroup(1, hopCNO, HoppingPropensity(&E_s, E_n, T), hopExecs);
  
  //! [adding cellcentered events]

  //! [adding printer]
  sim.reserveTimePeriodicActions(PAction::SIZE);
  sim.addTimePeriodicAction(PAction::PRINT,
			    PrintASCII("snapshot"),
			    0.05*approxDepTime, true);
  //! [adding printer]

  //! [running simulation]
  sim.run(approxDepTime);
  sim.removeOverLatticeEvent(PSOverLatticeEvents::DEPOSITION);
  sim.run(0.1*approxDepTime);
  //! [running simulation]

  return 0;
}
