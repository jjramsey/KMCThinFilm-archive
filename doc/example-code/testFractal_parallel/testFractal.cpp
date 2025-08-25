/* 
   Any comments in this code of the form "//! [...]" are used to
   assist Doxygen in documenting this file.  */

//! [headers]
#include <KMCThinFilm/Simulation.hpp>

#if KMC_PARALLEL
#include <KMCThinFilm/RandNumGenDCMT.hpp>
#else
#include <KMCThinFilm/RandNumGenMT19937.hpp>
#endif

#include "EventsAndActions.hpp"

#include <boost/lexical_cast.hpp>
//! [headers]

//! [using decl]
using namespace KMCThinFilm;
//! [using decl]

int main(int argc, char *argv[]) {

  //! [mpi init]
#if KMC_PARALLEL
  MPI_Init(&argc, &argv);
#endif
  //! [mpi init]

  //! [hardcoded parameters]
  double F = 1, DoverF = 1e5, maxCoverage = 4;
  int domainSize = 256;
  unsigned int seedGlobal = 42;
  SolverId::Type sId = SolverId::DYNAMIC_SCHULZE;

  TimeIncr::SchemeVars schemeVars;
  schemeVars.setSchemeName(TimeIncr::SchemeName::MAX_AVG_PROPENSITY_PER_POSS_EVENT);

  schemeVars.setSchemeParam(TimeIncr::SchemeParam::NSTOP, 
#ifndef BAD_NSTOP
                            1
#else
                            100
#endif
                            );
  //! [hardcoded parameters]
  
  //! [deposition time]
  double approxDepTime = maxCoverage/F;
  //! [deposition time]

  //! [initializing simulation]
  LatticeParams latParams;
  latParams.numIntsPerCell = FIntVal::SIZE;
  latParams.globalPlanarDims[0] = latParams.globalPlanarDims[1] = domainSize;
  latParams.ghostExtent[0] = latParams.ghostExtent[1] = 1;

#ifdef USE_COMPACT_DECOMP
  latParams.parallelDecomp = LatticeParams::COMPACT;
#endif

  Simulation sim(latParams);
  //! [initializing simulation]

  sim.setSolver(sId);

  //! [setting rng]
#if KMC_PARALLEL
  RandNumGenSharedPtr rng(new RandNumGenDCMT(sim.procID(),
                                             seedGlobal,
                                             123*sim.procID() + 456,
                                             RandNumGenDCMT::P521));
#else
  RandNumGenSharedPtr rng(new RandNumGenMT19937(seedGlobal));
#endif

  sim.setRNG(rng);
  //! [setting rng]

  //! [setting parallel time incr scheme]
  sim.setTimeIncrScheme(schemeVars);
  //! [setting parallel time incr scheme]

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

  tmpExecCNO.back().resetOffsets(2);
  tmpExecCNO.back().addOffset(1, hopCNO.getOffset(HopOffset::RIGHT));
  hopExecs.addEventExecutor(FCellCenteredEvents::HOP_RIGHT,
                            HoppingExecute(), tmpExecCNO);

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
			    PrintASCII("outFile_ProcCoords" + 
                                       boost::lexical_cast<std::string>(sim.commCoord(0)) + "_" +
                                       boost::lexical_cast<std::string>(sim.commCoord(1)) + "_snapshot"),
			    0.05*approxDepTime, true);
  //! [adding printer]

  //! [running simulation]
  sim.run(approxDepTime);
  sim.removeOverLatticeEvent(FOverLatticeEvents::DEPOSITION);
  sim.run(0.1*approxDepTime);
  //! [running simulation]

  //! [mpi finalize]
#if KMC_PARALLEL
  MPI_Finalize();
#endif
  //! [mpi finalize]

  return 0;
}
