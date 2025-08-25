/* 
   Any comments in this code of the form "//! [...]" are used to
   assist Doxygen in documenting this file.  */

//! [headers]
#include <KMCThinFilm/Simulation.hpp>
#include <KMCThinFilm/RandNumGenMT19937.hpp>

#include "EventsAndActions.hpp"
//! [headers]

//! [using decl]
using namespace KMCThinFilm;
//! [using decl]

int main() {

  //! [hardcoded parameters]
  double F = 1, DoverF = 1e5, maxCoverage = 4;
  int domainSize = 256;
  unsigned int seed = 42;
  SolverId::Type sId = SolverId::DYNAMIC_SCHULZE;
  //! [hardcoded parameters]
  
  //! [deposition time]
  double approxDepTime = maxCoverage/F;
  //! [deposition time]

  //! [initializing simulation]
  LatticeParams latParams;
  latParams.numIntsPerCell = FIntVal::SIZE;
  latParams.globalPlanarDims[0] = latParams.globalPlanarDims[1] = domainSize;

  Simulation sim(latParams);
  //! [initializing simulation]

  //! [setting solver and rng]
  sim.setSolver(sId);

  RandNumGenSharedPtr rng(new RandNumGenMT19937(seed));

  sim.setRNG(rng);
  //! [setting solver and rng]

  //! [adding overlattice events]
  std::vector<CellNeighOffsets> tmpExecCNO;
  tmpExecCNO.reserve(1);

  sim.reserveOverLatticeEvents(FOverLatticeEvents::SIZE);

  tmpExecCNO.push_back(CellNeighOffsets(1));
  sim.reserveOverLatticeEvents(FOverLatticeEvents::SIZE);
  sim.addOverLatticeEvent(FOverLatticeEvents::DEPOSITION,
			  F, DepositionExecute, tmpExecCNO);
  //! [adding overlattice events]
  
  //! [making cellneighoffsets]
  CellNeighOffsets hopCNO(HopOffset::SIZE);

  hopCNO.addOffset(HopOffset::UP,    CellIndsOffset(0,+1,0));
  hopCNO.addOffset(HopOffset::DOWN,  CellIndsOffset(0,-1,0));
  hopCNO.addOffset(HopOffset::LEFT,  CellIndsOffset(-1,0,0));
  hopCNO.addOffset(HopOffset::RIGHT, CellIndsOffset(+1,0,0));

  hopCNO.addOffset(HopOffset::RIGHT_UP,   CellIndsOffset(+1,+1,0));
  hopCNO.addOffset(HopOffset::RIGHT_DOWN, CellIndsOffset(+1,-1,0));
  hopCNO.addOffset(HopOffset::LEFT_UP,    CellIndsOffset(-1,+1,0));
  hopCNO.addOffset(HopOffset::LEFT_DOWN,  CellIndsOffset(-1,-1,0));
  //! [making cellneighoffsets]

  //! [adding cellcentered events]
  sim.reserveCellCenteredEventGroups(1,FCellCenteredEvents::SIZE);

  EventExecutorGroup hopExecs(FCellCenteredEvents::SIZE);
  tmpExecCNO.clear();
  tmpExecCNO.push_back(CellNeighOffsets(2));
  tmpExecCNO.back().addOffset(1, hopCNO.getOffset(HopOffset::LEFT));

  hopExecs.addEventExecutor(FCellCenteredEvents::HOP_LEFT,
                            HoppingExecute(), tmpExecCNO);

  //! [setting execcno for rightward hop]
  tmpExecCNO.back().resetOffsets(2);
  tmpExecCNO.back().addOffset(1, hopCNO.getOffset(HopOffset::RIGHT));
  hopExecs.addEventExecutor(FCellCenteredEvents::HOP_RIGHT,
                            HoppingExecute(), tmpExecCNO);
  //! [setting execcno for rightward hop]

  tmpExecCNO.back().resetOffsets(2);
  tmpExecCNO.back().addOffset(1, hopCNO.getOffset(HopOffset::UP));
  hopExecs.addEventExecutor(FCellCenteredEvents::HOP_UP,
                            HoppingExecute(), tmpExecCNO);

  tmpExecCNO.back().resetOffsets(2);
  tmpExecCNO.back().addOffset(1, hopCNO.getOffset(HopOffset::DOWN));
  hopExecs.addEventExecutor(FCellCenteredEvents::HOP_DOWN,
                            HoppingExecute(), tmpExecCNO);

  sim.addCellCenteredEventGroup(1, hopCNO,
                                HoppingPropensity(DoverF*F),
                                hopExecs);
  //! [adding cellcentered events]

  //! [adding printer]
  sim.reserveTimePeriodicActions(PAction::SIZE);
  sim.addTimePeriodicAction(PAction::PRINT,
			    PrintASCII("snapshot"),
			    0.05*approxDepTime, true);
  //! [adding printer]

  //! [running simulation]
  sim.run(approxDepTime);
  sim.removeOverLatticeEvent(FOverLatticeEvents::DEPOSITION);
  sim.run(0.1*approxDepTime);
  //! [running simulation]

  return 0;
}
