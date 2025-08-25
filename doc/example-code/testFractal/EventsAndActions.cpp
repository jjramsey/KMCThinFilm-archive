/* 
   Comments in this code of the form "//! [...]" are used to assist
   Doxygen in documenting this file.
*/

#include "EventsAndActions.hpp"

#include <fstream>
#include <boost/lexical_cast.hpp>

#include <KMCThinFilm/ErrorHandling.hpp>

using namespace KMCThinFilm;

//! [dep execute]
void DepositionExecute(const CellInds & ci,
		       const SimulationState & simState,
		       Lattice & lattice) {

  int currVal = lattice.getInt(ci, FIntVal::HEIGHT);
  lattice.setInt(ci, FIntVal::HEIGHT, currVal + 1);

}
//! [dep execute]

//! [hop prop]
void HoppingPropensity::operator()(const CellNeighProbe & cnp, 
                                   std::vector<double> & propensityVec) const {

  int currHeight = cnp.getInt(cnp.getCellToProbe(HopOffset::SELF), FIntVal::HEIGHT);

  if (currHeight > 0) {
    if ((currHeight > cnp.getInt(cnp.getCellToProbe(HopOffset::UP), FIntVal::HEIGHT))
        && (currHeight > cnp.getInt(cnp.getCellToProbe(HopOffset::DOWN), FIntVal::HEIGHT))
        && (currHeight > cnp.getInt(cnp.getCellToProbe(HopOffset::LEFT), FIntVal::HEIGHT))
        && (currHeight > cnp.getInt(cnp.getCellToProbe(HopOffset::RIGHT), FIntVal::HEIGHT))
        && (currHeight > cnp.getInt(cnp.getCellToProbe(HopOffset::RIGHT_UP), FIntVal::HEIGHT))
        && (currHeight > cnp.getInt(cnp.getCellToProbe(HopOffset::RIGHT_DOWN), FIntVal::HEIGHT))
        && (currHeight > cnp.getInt(cnp.getCellToProbe(HopOffset::LEFT_UP), FIntVal::HEIGHT))
        && (currHeight > cnp.getInt(cnp.getCellToProbe(HopOffset::LEFT_DOWN), FIntVal::HEIGHT))) {

      propensityVec[FCellCenteredEvents::HOP_LEFT] = 
        propensityVec[FCellCenteredEvents::HOP_RIGHT] = 
        propensityVec[FCellCenteredEvents::HOP_UP] = 
        propensityVec[FCellCenteredEvents::HOP_DOWN] = D_;
    }
  }

}
//! [hop prop]

//! [hop exec constructor]
HoppingExecute::HoppingExecute(FCellCenteredEvents::Type hopDir) {

  switch (hopDir) {
  case FCellCenteredEvents::HOP_UP:
    jump_i_ = 0;
    jump_j_ = 1;
    break;
  case FCellCenteredEvents::HOP_DOWN:
    jump_i_ = 0;
    jump_j_ = -1;
    break;
  case FCellCenteredEvents::HOP_LEFT:
    jump_i_ = -1;
    jump_j_ = 0;
    break;
  case FCellCenteredEvents::HOP_RIGHT:
    jump_i_ = 1;
    jump_j_ = 0;
    break;
  default:
    exitWithMsg("Bad hop direction!");
  }
}
//! [hop exec constructor]

//! [hop exec op]
void HoppingExecute::operator()(const CellInds & ci,
				const SimulationState & simState,
				Lattice & lattice) const {

  CellInds ciTo(ci.i + jump_i_, ci.j + jump_j_, ci.k);

  int currFrom = lattice.getInt(ci, FIntVal::HEIGHT);
  int currTo = lattice.getInt(ciTo, FIntVal::HEIGHT);
  
  lattice.setInt(ci, FIntVal::HEIGHT, currFrom - 1);
  lattice.setInt(ciTo, FIntVal::HEIGHT, currTo + 1);
}
//! [hop exec op]

//! [print op]
void PrintASCII::operator()(const SimulationState & simState, Lattice & lattice) {

  ++snapShotCntr_;

  std::string fName = fNameRoot_ + boost::lexical_cast<std::string>(snapShotCntr_) + ".dat";
    
  std::ofstream outFile(fName.c_str());

  LatticePlanarBBox localPlanarBBox;
  lattice.getLocalPlanarBBox(false, localPlanarBBox);

  int iminGlobal = localPlanarBBox.imin;
  int jminGlobal = localPlanarBBox.jmin;

  // "P1" here is short for "Plus 1".
  int imaxP1Global = localPlanarBBox.imaxP1;
  int jmaxP1Global = localPlanarBBox.jmaxP1;

  outFile << "# " << iminGlobal << " " << imaxP1Global << " " << jminGlobal << " " << jmaxP1Global
	  << " time:" << simState.elapsedTime() << "\n";

  CellInds ci; ci.k = 0;
  for (ci.i = localPlanarBBox.imin; ci.i < localPlanarBBox.imaxP1; ++(ci.i)) {
    for (ci.j = localPlanarBBox.jmin; ci.j < localPlanarBBox.jmaxP1; ++(ci.j)) {
      outFile << ci.i << " " << ci.j << " "
	      << lattice.getInt(ci, FIntVal::HEIGHT) << "\n";
    }
  }

  outFile.close();
}
//! [print op]
