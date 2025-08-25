#include "CellInds.hpp"

namespace KMCThinFilm {

  CellInds operator+(const CellInds & ci, const CellIndsOffset & offset) {
    return CellInds(ci.i + offset.i, ci.j + offset.j, ci.k + offset.k);
  }

  CellIndsOffset operator+(const CellIndsOffset & offset1, const CellIndsOffset & offset2) {
    return CellIndsOffset(offset1.i + offset2.i,
			  offset1.j + offset2.j,
			  offset1.k + offset2.k);
  }

  CellIndsOffset CellIndsOffset::operator-() const {
    return CellIndsOffset(-i, -j, -k);
  }

}
