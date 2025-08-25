#ifndef CELL_NEIGH_OFFSETS_HPP
#define CELL_NEIGH_OFFSETS_HPP

#include "CellInds.hpp"
#include <boost/scoped_ptr.hpp>

/*! \file
  \brief Defines the CellNeighOffsets class.
*/

namespace KMCThinFilm {

  /*! Collection of CellIndsOffset objects that are used to determine
      an event's propensity.

      Each offset in the collection is associated with an integer
      ID. No matter what, this collection always contains a
      CellIndsOffset object with the indices (0,0,0), and this offset
      has the integer ID of zero.
  */
  class CellNeighOffsets {
  public:

    /*! Constructs a CellNeighOffsets object. */
    explicit CellNeighOffsets(int numberOfOffsets /*!< Number of offsets to be added. */);

    //! \cond HIDE_FROM_DOXYGEN
    CellNeighOffsets(const CellNeighOffsets & cno);
    CellNeighOffsets & operator=(const CellNeighOffsets & rhs);
    //! \endcond
    
    /*! Adds an additional offset.

      The argument <VAR>whichOffset</VAR> is not allowed to be zero,
      since that integer ID is already associated with the offset
      (0,0,0).
     */
    void addOffset(int whichOffset /*!< Integer ID of offset. Must
                                      be >= 1 and less than numOffsets() */,
		   const CellIndsOffset & offset /*!< The offset to be added. */);

    /*! Removes all previous offsets and allocates memory for
        <VAR>numberOfNewOffsets</VAR> CellIndsOffset objects. */
    void resetOffsets(int numberOfNewOffsets);

    /*! Removes all previously added offsets.

      Note that no new offsets may be added until resetOffsets() is called. */
    void clearOffsets();

    /*! Returns the offset with the integer ID
        <VAR>whichOffset</VAR>. Note that, here,
        <VAR>whichOffset</VAR> <EM>is</EM> allowed to be zero. */
    const CellIndsOffset & getOffset(int whichOffset /*!< Integer ID of offset. */ ) const;

    /*! Number of offsets stored. */
    int numOffsets() const;

    //! \cond HIDE_FROM_DOXYGEN
    ~CellNeighOffsets();
    //! \endcond

  private:
    class Impl_;
    boost::scoped_ptr<Impl_> pImpl_;
  };

}

#endif /* CELL_NEIGH_OFFSETS_HPP */
