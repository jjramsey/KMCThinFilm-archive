/* 
   Comments in this code of the form "//! [...]" are used to assist
   Doxygen in documenting this file.
*/

#include "EventsAndActions.hpp"

#include <vector>
#include <fstream>
#include <cmath>

#include <boost/lexical_cast.hpp>

#include <KMCThinFilm/ErrorHandling.hpp>

using namespace KMCThinFilm;

//! [dep execute]
void DepositionExecute(const CellInds & ci,
		       const SimulationState & simState,
		       Lattice & lattice) {

  int currVal = lattice.getInt(ci, PSIntVal::HEIGHT);
  lattice.setInt(ci, PSIntVal::HEIGHT, currVal + 1);

}
//! [dep execute]

//! [hop prop]
void HoppingPropensity::operator()(const CellNeighProbe & cnp,
                                   std::vector<double> & propensityVec) const {

  KMCThinFilm::CellToProbe ctpSelf = cnp.getCellToProbe(HopOffset::SELF);

  int currHeight = cnp.getInt(ctpSelf, PSIntVal::HEIGHT);

  if (currHeight > 0) {
    
    int n = 0;

    if (cnp.getInt(cnp.getCellToProbe(HopOffset::UP), PSIntVal::HEIGHT) >= currHeight) {
      ++n;
    }

    if (cnp.getInt(cnp.getCellToProbe(HopOffset::DOWN), PSIntVal::HEIGHT) >= currHeight) {
      ++n;
    }

    if (cnp.getInt(cnp.getCellToProbe(HopOffset::LEFT), PSIntVal::HEIGHT) >= currHeight) {
      ++n;
    }

    if (cnp.getInt(cnp.getCellToProbe(HopOffset::RIGHT), PSIntVal::HEIGHT) >= currHeight) {
      ++n;
    }

    double E = cnp.getFloat(ctpSelf, PSFloatVal::E_s) + n*E_n_;

    double p = k_*std::exp(-E/kBT_);
    
    for (int i = 0; i < PSCellCenteredEvents::SIZE; ++i) {
      propensityVec[i] = p;
    }

  }
}
//! [hop prop]

//! [hop exec op]
void HoppingExecute::operator()(const CellInds & ci,
				const SimulationState & simState,
                                const Lattice & lattice,
                                std::vector<CellsToChange> & ctcVec) const {

  CellsToChange & ctc = ctcVec[0];
  ctc.setCenter(ci);

  int currFrom = lattice.getInt(ci, PSIntVal::HEIGHT);
  int currTo = ctc.getInt(1, PSIntVal::HEIGHT);
  
  ctc.setInt(0, PSIntVal::HEIGHT, currFrom - 1);
  ctc.setInt(1, PSIntVal::HEIGHT, currTo + 1);
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
	      << lattice.getInt(ci, PSIntVal::HEIGHT) << "\n";
    }
  }

  outFile.close();
}
//! [print op]
