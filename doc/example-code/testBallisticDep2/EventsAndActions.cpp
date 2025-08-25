#include "EventsAndActions.hpp"

#include <fstream>

#include <boost/lexical_cast.hpp>

using namespace KMCThinFilm;

//! [dep exec constructor]
DepositionExecute::DepositionExecute(const LatticePlanarBBox & planarBBox)
: activeZoneHeights_(new IntArray2D) {

  // The resize() member function will initialize the elements of this
  // array to zero.
  activeZoneHeights_->resize(boost::extents[planarBBox.imaxP1 - planarBBox.imin]
                             [planarBBox.jmaxP1 - planarBBox.jmin]);
  
  neighOffsets_.reserve(4);
  neighOffsets_.push_back(CellIndsOffset( 0,-1));
  neighOffsets_.push_back(CellIndsOffset( 0,+1));
  neighOffsets_.push_back(CellIndsOffset(-1, 0));
  neighOffsets_.push_back(CellIndsOffset(+1, 0));

}
//! [dep exec constructor]

//! [dep exec op]
void DepositionExecute::operator()(const KMCThinFilm::CellInds & ci,
                                   const KMCThinFilm::SimulationState & simState,
                                   const KMCThinFilm::Lattice & lattice,
                                   std::vector<KMCThinFilm::CellsToChange> & ctcVec) {
  
  /* Ballistic deposition algorithm for cubic lattices from Meakin and
     Krug, Physical Review A, vol. 46, num. 6, pp. 3390-3399 (1992).*/

  int kDepAtom = (*activeZoneHeights_)[ci.i][ci.j];

  CellInds ciTo(ci.i, ci.j, kDepAtom);

  CellsToChange & ctc = ctcVec[0];

  ctc.setCenter(ciTo);
  
  ctc.addLatticePlanes(ciTo.k - ci.k);
  ctc.setInt(0, BDIntVal::IS_OCCUPIED, 1);

  for (std::vector<CellIndsOffset>::const_iterator offsetItr = neighOffsets_.begin(),
         offsetItrEnd = neighOffsets_.end(); offsetItr != offsetItrEnd; ++offsetItr) {

    CellInds neigh = ciTo + *offsetItr;

    int i = lattice.wrapI(neigh);
    int j = lattice.wrapJ(neigh);

    if ((*activeZoneHeights_)[i][j] < kDepAtom) {
      (*activeZoneHeights_)[i][j] = kDepAtom;
    }
  }

  ++((*activeZoneHeights_)[ci.i][ci.j]);
}
//! [dep exec op]

//! [mix prop op]
void ColorMixPropensity::operator()(const KMCThinFilm::CellNeighProbe & cnp,
                                    std::vector<double> & propensityVec) const {
  
  if (cnp.getInt(cnp.getCellToProbe(MIX_OFFSET::SELF), BDIntVal::IS_OCCUPIED)) {

    int numNeighs = 0;

    // Visiting four lateral neighbors
    for (int whichOffset = 1; whichOffset <= 4; ++whichOffset) {

      if (cnp.getInt(cnp.getCellToProbe(whichOffset), BDIntVal::IS_OCCUPIED)) {
        ++numNeighs;
      }

    }

    CellToProbe downCell = cnp.getCellToProbe(MIX_OFFSET::DOWN);
    if (cnp.belowLatticeBottom(downCell) || cnp.getInt(downCell, BDIntVal::IS_OCCUPIED)) {
      ++numNeighs;
    }

    CellToProbe upCell = cnp.getCellToProbe(MIX_OFFSET::UP);
    if ((!cnp.exceedsLatticeHeight(upCell)) && cnp.getInt(upCell, BDIntVal::IS_OCCUPIED)) {
      ++numNeighs;
    }

    propensityVec[CellCenteredEvents::COLOR_MIXING] = mixPropPerNeighbor_*numNeighs;

  }
}
//! [mix prop op]

//! [mix exec ctor]
ColorMixExecute::ColorMixExecute(RandNumGenSharedPtr rng,
                                 const CellNeighOffsets * mixCNO,
                                 int * numMixes)
  : rng_(rng), 
    mixCNO_(mixCNO),
    numMixes_(numMixes) {

  *numMixes_ = 0;

}
//! [mix exec ctor]

//! [mix exec op]
void ColorMixExecute::operator()(const KMCThinFilm::CellInds & ci,
                                 const KMCThinFilm::SimulationState & simState,
                                 const KMCThinFilm::Lattice & lattice,
                                 std::vector<KMCThinFilm::CellsToChange> & ctcVec) {

  double color = lattice.getFloat(ci, BDFloatVal::COLOR);
  int numColors = 1;

  // Visiting four lateral neighbors
  for (int whichOffset = 1; whichOffset <= 4; ++whichOffset) {

    CellInds ciNeigh = ci + mixCNO_->getOffset(whichOffset);

    if (lattice.getInt(ciNeigh, BDIntVal::IS_OCCUPIED)) {
      color += lattice.getFloat(ciNeigh, BDFloatVal::COLOR);
      ++numColors;      
    }
      
  }

  CellInds ciDown = ci +  mixCNO_->getOffset(MIX_OFFSET::DOWN);

  bool belowLatticeBottom = (ciDown.k < 0);
  bool ciDownIsOccupied = (belowLatticeBottom || lattice.getInt(ciDown, BDIntVal::IS_OCCUPIED));

  if (ciDownIsOccupied) {
    color += (belowLatticeBottom ?
              (rng_->getNumInOpenIntervalFrom0To1()) :
              lattice.getFloat(ciDown, BDFloatVal::COLOR));
    ++numColors;
  }

  

  CellInds ciUp = ci +  mixCNO_->getOffset(MIX_OFFSET::UP);

  if ((ciUp.k < lattice.currHeight()) && lattice.getInt(ciUp, BDIntVal::IS_OCCUPIED)) {

    color += lattice.getFloat(ciUp, BDFloatVal::COLOR);
    ++numColors;

  }

  CellsToChange & ctc = ctcVec[0];

  ctc.setCenter(ci);
  ctc.setFloat(0, BDFloatVal::COLOR, color/numColors);

  ++(*numMixes_);

}
//! [mix exec op]

//! [print op]
void PrintPoint3D::operator()(const SimulationState & simState, Lattice & lattice) {
  ++snapShotCntr_;

  std::string fName = fNameRoot_ + boost::lexical_cast<std::string>(snapShotCntr_) + ".3D";
    
  std::ofstream outFile(fName.c_str());

  LatticePlanarBBox localPlanarBBox;
  lattice.getLocalPlanarBBox(false, localPlanarBBox);

  outFile << "x y z value\n";

  CellInds ci;
  int kMaxP1 = lattice.currHeight();

  for (ci.k = 0; ci.k < kMaxP1; ++(ci.k)) {
    for (ci.i = localPlanarBBox.imin; ci.i < localPlanarBBox.imaxP1; ++(ci.i)) {
      for (ci.j = localPlanarBBox.jmin; ci.j < localPlanarBBox.jmaxP1; ++(ci.j)) {

	if (lattice.getInt(ci, BDIntVal::IS_OCCUPIED) > 0) {
	  outFile << ci.i << " " << ci.j << " "
		  << ci.k << " " << lattice.getFloat(ci, BDFloatVal::COLOR) << "\n";
	}
      }
    }
  }

  outFile.close();
}
//! [print op]
