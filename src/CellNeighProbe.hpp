#ifndef CELL_NEIGH_PROBE_HPP
#define CELL_NEIGH_PROBE_HPP

#include <vector>
#include <boost/scoped_ptr.hpp>

#include "CellInds.hpp"

/*!\file
  \brief Defines the CellNeighProbe and CellToProbe classes
 */

namespace KMCThinFilm {

  class Lattice;

  /*! A largely opaque representation of a lattice cell for use with
      the CellNeighProbe class.
   */
  class CellToProbe {
    friend class CellNeighProbe;
  public:

    //! \cond HIDE_FROM_DOXYGEN
    CellToProbe()
      : ci_((CellInds())) /* Extra parens are to avoid C++'s "most
			     vexing parse", but I'm not sure if
			     they're needed. */
    {}
    //! \endcond

    /*! Returns the indices of the lattice cell pointed to by an
        instance of CellToProbe. */
    const CellInds & inds() {return ci_;}
  private:
    explicit CellToProbe(const CellInds & ci)
      : ci_(ci)
    {}

    CellInds ci_;
  };

  /*! A class used by a CellCenteredPropensity function object to
      probe the states of neighboring lattice cells in order to
      determine a cell-centered event's propensity.

      Each neighboring lattice cell is identified by one of the
      integer IDs used to label the offsets inside a CellNeighOffsets
      object, and the indices of this neighboring lattice cell are
      that of the current cell (about which an event is centered) plus
      the offset corresponding to its integer ID.
   */
  class CellNeighProbe {
  public:

    /*! [<STRONG>ADVANCED</STRONG>] Constructs a CellNeighProbe and
        attaches it to a Lattice object. Mainly for use in testing and
        debugging. */
    explicit CellNeighProbe(const Lattice * lattice = NULL);

    //! \cond HIDE_FROM_DOXYGEN
    CellNeighProbe(const CellNeighProbe & ctp);
    CellNeighProbe & operator=(const CellNeighProbe & rhs);
    //! \endcond

    /*! [<STRONG>ADVANCED</STRONG>] Attaches a Lattice object,
        replacing whatever lattice object that was attached
        before. Mainly for use in testing and debugging. */
    void attachLattice(const Lattice * lattice);    

    /*! [<STRONG>ADVANCED</STRONG>] Attaches a vector of
        CellIndsOffset objects (generated from a CellNeighOffsets
        instance). Mainly for use in testing and debugging. */
    void attachCellInds(const CellInds * ci /*!< Pointer to the
					      indices of the cell
					      about which an event is
					      centered. */, 
			const std::vector<CellIndsOffset> * cioVecPtr /*!<
                                                                         Pointer
                                                                         to
                                                                         the
                                                                         vector
                                                                         of
                                                                         CellIndsOffset
                                                                         objects
                                                                         that
                                                                         define
                                                                         the
                                                                         neighboring
                                                                         lattice
                                                                         cells
                                                                         used
                                                                         to
                                                                         calculate
                                                                         the
                                                                         propensity
                                                                         of
                                                                         the
                                                                         event. */);

    /*! Retrieves the lattice cell to be probed from its integer ID. */
    CellToProbe getCellToProbe(int probedCellInd /*!< Integer ID of
                                                    the lattice cell
                                                    (and its
                                                    corresponding
                                                    offset). */) const;

    /*! Retrieves one of the floating-point values (if present) from
        the lattice cell pointed to by <VAR>ctp</VAR>. The parameter
        <VAR>whichFloat</VAR> indicates which of the floating-point
        values to return. */
    double getFloat(const CellToProbe & ctp, int whichFloat) const;

    /*! Retrieves one of the integer values (if present) from the
        lattice cell pointed to by <VAR>ctp</VAR>. The parameter
        <VAR>whichInt</VAR> indicates which of the integer values to
        return. */
    int getInt(const CellToProbe & ctp, int whichInt) const;

    /*! Indicates whether the height coordinate of the (possibly
        non-existent) lattice cell pointed to by <VAR>ctp</VAR>
        exceeds the current height of the lattice. If true, it implies
        that this lattice cell does not actually exist yet in the
        simulation. */
    bool exceedsLatticeHeight(const CellToProbe & ctp) const;

    /*! Indicates whether the height coordinate of the (possibly
        non-existent) lattice cell pointed to by <VAR>ctp</VAR> is
        less than zero. If true, it implies that this lattice cell
        does not and will not actually exist in the simulation. */
    bool belowLatticeBottom(const CellToProbe & ctp) const;
    
    //! \cond HIDE_FROM_DOXYGEN
    ~CellNeighProbe();
    //! \endcond

  private:
    class Impl_;
    boost::scoped_ptr<Impl_> pImpl_;
  };
  
}

#endif /* CELL_NEIGH_PROBE_HPP */
