#include <iostream>

#include <KMCThinFilm/Simulation.hpp>
#include <KMCThinFilm/RandNumGenMT19937.hpp>
#include <KMCThinFilm/MakeEnum.hpp>

#include "EventsAndActions.hpp"
#include "InitLattice.hpp"

using namespace KMCThinFilm;

int main() {

  // Parameters used in the simulation
  double F = 1, maxCoverage = 4;
  int domainSize = 100;
  unsigned int seed = 42;
  SolverId::Type sId = SolverId::DYNAMIC_SCHULZE;

  double approxDepTime = maxCoverage/F;

  //! [init sim]
  RandNumGenSharedPtr rng(new RandNumGenMT19937(seed));

  LatticeParams latParams;
  latParams.numIntsPerCell = BDIntVal::SIZE;
  latParams.numFloatsPerCell = BDFloatVal::SIZE;
  latParams.globalPlanarDims[0] = latParams.globalPlanarDims[1] = domainSize;
  latParams.setEmptyCellVals = SetEmptyCellWithRandColor(rng);
  latParams.numPlanesToReserve = 100;

  Simulation sim(latParams);
  //! [init sim]

  sim.setSolver(sId);

  sim.setRNG(rng);

  sim.reserveOverLatticeEvents(OverLatticeEvents::SIZE);
  sim.addOverLatticeEvent(OverLatticeEvents::DEPOSITION,
                          F, DepositionExecute());

  //! [add cellcen events]
  CellNeighOffsets mixCNO(MIX_OFFSET::SIZE);
  mixCNO.addOffset(MIX_OFFSET::NORTH, CellIndsOffset(+1, 0, 0));
  mixCNO.addOffset(MIX_OFFSET::SOUTH, CellIndsOffset(-1, 0, 0));
  mixCNO.addOffset(MIX_OFFSET::WEST,  CellIndsOffset( 0,-1, 0));
  mixCNO.addOffset(MIX_OFFSET::EAST,  CellIndsOffset( 0,+1, 0));
  mixCNO.addOffset(MIX_OFFSET::UP,    CellIndsOffset( 0, 0,+1));
  mixCNO.addOffset(MIX_OFFSET::DOWN,  CellIndsOffset( 0, 0,-1));

  int numMixes;

  EventExecutorGroup mixExec(CellCenteredEvents::SIZE);

  mixExec.addEventExecutor(CellCenteredEvents::COLOR_MIXING,
                           ColorMixExecute(rng, &mixCNO, &numMixes));

  sim.reserveCellCenteredEventGroups(1, CellCenteredEvents::SIZE);
  sim.addCellCenteredEventGroup(1, mixCNO,
                                ColorMixPropensity(10*F),
                                mixExec);
  //! [add cellcen events]
  
  sim.reserveTimePeriodicActions(PAction::SIZE);
  sim.addTimePeriodicAction(PAction::PRINT,
                            PrintPoint3D("snapshot"),
                            0.05*approxDepTime, true);

  sim.run(approxDepTime);

  //! [print nummixes]
  std::cout << "Number of color mixes = " << numMixes << "\n";
  //! [print nummixes]

  return 0;
}
