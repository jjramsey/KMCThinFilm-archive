#include "InitLattice.hpp"
#include "EventsAndActions.hpp"

#include <fstream>

#include <boost/array.hpp>
#include <boost/lexical_cast.hpp>

using namespace KMCThinFilm;

void InitLatticeFromFile::operator()(Lattice & lattice) const {
  
  std::ifstream inpFile(inpFName_.c_str());

  boost::array<int,2> globalPlanarDims;

  inpFile >> globalPlanarDims[0] >> globalPlanarDims[1];

  LatticePlanarBBox globalPlanarBBox;
  lattice.getGlobalPlanarBBox(globalPlanarBBox);

  exitOnCondition((globalPlanarDims[0] != (globalPlanarBBox.imaxP1 - globalPlanarBBox.imin)) ||
                  (globalPlanarDims[1] != (globalPlanarBBox.jmaxP1 - globalPlanarBBox.jmin)),
                  "Mismatch in lattice dimensions and dimensions of strain eng. den. array.");

  lattice.addPlanes(1);

  LatticePlanarBBox localPlanarBBox;
  lattice.getLocalPlanarBBox(false, localPlanarBBox);

  int maxLineNumP1 = (globalPlanarDims[0])*(globalPlanarDims[1]);

  CellInds ci; ci.k = 0;
  for (int lineNum = 0; lineNum < maxLineNumP1; ++lineNum) {
    
    int i, j;
    double E_s;

    inpFile >> i >> j >> E_s;

    if ((i >= localPlanarBBox.imin) && (i < localPlanarBBox.imaxP1) && 
        (j >= localPlanarBBox.jmin) && (j < localPlanarBBox.jmaxP1)) {

      ci.i = i;
      ci.j = j;
        
      lattice.setFloat(ci, PSFloatVal::E_s, E_s);
    }

  }
    
  inpFile.close();

}
