#ifndef CELL_INDS_HPP
#define CELL_INDS_HPP

#include "IJK.hpp"

// Note: This header file is documented via Doxygen
// <http://www.doxygen.org>. Comments for Doxygen begin with '/*!' or
// '//!', and descriptions of functions, class and member functions
// occur *before* their corresponding class declarations and function
// prototypes.

/*! \file
  \brief Defines the CellInds and CellInds offset classes, and their respective operators.
 */

namespace KMCThinFilm {


  /*! Indices of a lattice cell.

    This class has three <EM>public</EM> data members, i, j, and k,
    which represent the first, second, and third indices. If the
    lattice has primitive lattice vectors ai, aj, and ak, then the
    location of the lattice cell may be taken as i*ai + j*aj + k*ak.

    \see CellIndsOffset KMCThinFilm::operator+()
   */
  class CellInds : public IJK {
  public:
    
    /*! Constructs an empty CellInds. */
    CellInds() {i = j = k = 0;}

    /*! Constructs a CellInds instance from a doublet or triplet of
        integers.*/
    CellInds(int ii /*!< First index of a lattice cell */, 
	     int jj /*!< Second index of lattice cell */, 
	     int kk = 0 /*!< Third index of lattice cell */) {
      i = ii; j = jj; k = kk;
    }

  };

  /*! An offset to indices of some lattice cell.

    Like the CellInds class, it has three <EM>public</EM> data
    members, i, j, and k, which represent the offsets to be added to a
    given set of lattice cell indices. For example, one may write

    \code
    CellInds ci(1,2,3);
    CellInds offset(1,1,-1);

    CellInds ci_w_offset = ci + offset;
    \endcode

    where the indices of ci_w_offset now equal 2, 3, and 2.

    \see CellInds KMCThinFilm::operator+()
   */
  class CellIndsOffset : public IJK {
  public:

    /*! Constructs an empty CellIndsOffset. */
    CellIndsOffset() {i = j = k = 0;}

    /*! Constructs a CellIndsOffset instance from a doublet or triplet of
        integers.*/
    CellIndsOffset(int ii /*!< Offset to the first index of a lattice cell */, 
		   int jj /*!< Offset to the second index of lattice cell */, 
		   int kk = 0 /*!< Offset to the third index of lattice cell */) {
      i = ii; j = jj; k = kk;
    }

    /*! Unary minus operator. Returns a sign-reversed copy of the offset. */
    CellIndsOffset operator-() const;
  };

  /*! Returns a CellInds instance with the indices ci.i + offset.i, ci.j +
    offset.j, and ci.k + offset.k.*/
  CellInds operator+(const CellInds & ci, const CellIndsOffset & offset);

  /*! Returns a CellIndsOffset instance with the indices offset1.i +
    offset2.i, offset1.j + offset2.j, and offset1.k + offset2.k.*/
  CellIndsOffset operator+(const CellIndsOffset & offset1, const CellIndsOffset & offset2);

}

#endif /* CELL_INDS_HPP */
