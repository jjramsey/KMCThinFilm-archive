#ifndef LATTICE_HPP
#define LATTICE_HPP

#include "KMC_Config.hpp"

#include <vector>
#include <deque>
#include <set>

#if KMC_PARALLEL
#include <mpi.h>
#endif

#include <boost/array.hpp>
#include <boost/function.hpp>
#include <boost/noncopyable.hpp>
#include <boost/scoped_ptr.hpp>

#include "CellInds.hpp"

// Note: This header file is documented via Doxygen
// <http://www.doxygen.org>. Comments for Doxygen begin with '/*!' or
// '//!', and descriptions of functions, class and member functions
// occur *before* their corresponding class declarations and function
// prototypes.

/*! \file
  \brief Defines the Lattice class and related parameter objects and typedefs.
 */

namespace KMCThinFilm {

  class Lattice; // Need to forward declare this for the sake of the
		 // following typedefs.

  /*! Signature of a function or function object to be used to set the
      integers and floating-point numbers in an empty lattice cell to
      their proper values.

      Note that when this function is used, the vectors
      <VAR>emptyIntVals</VAR> and <VAR>emptyFloatVals</VAR> have
      already had their sizes set to the number of integer values and
      number of floating-point values per lattice cell,
      respectively. */
  typedef boost::function<void (const CellInds &,
				const Lattice &,
				std::vector<int> & emptyIntVals,
				std::vector<double> & emptyFloatVals)> SetEmptyCellVals;

  /*! Bounding box of in-plane values of lattice cell coordinates. */
  struct LatticePlanarBBox {
    int imin /*! Minimum value of the first lattice cell coordinate. */,
      imaxP1 /*! One more than the maximum value of the first lattice cell coordinate */,
      jmin /*! Minimum value of the second lattice cell coordinate */,
      jmaxP1 /*! One more than the maximum value of the second lattice cell coordinate */;
  };

  /*! Signature of a function or function object to be used to
      initialize a lattice. */
  typedef boost::function<void (Lattice & lattice)> LatticeInitializer;

  /*! Simple function object used to add empty planes to the lattice,
      especially during initialization of the lattice. */
  class AddEmptyPlanes {
  public:    
    /*! Constructor, defines number of lattice planes to add. */
    AddEmptyPlanes(int numPlanesToAdd);

    /*! Overload of operator() that performs the actual adding of planes. */
    void operator()(Lattice & lattice) const;
  private:
    int numPlanesToAdd_;
  };

  /*! Parameter object used in constructing a Lattice object */
  struct LatticeParams {

    /*! Methods of parallel decomposition */
    enum ParallelDecomp {
      COMPACT /*!< The lattice is decomposed such that the perimeter
                 of the part of the lattice owned by each process is a
                 minimum. */,
      ROW /*!< The lattice is decomposed into
             strips of size globalPlanarDims[0]/N by
             globalPlanarDims[1], where N is the number of
             processors. */
    };

    LatticeParams()
      : numIntsPerCell(0), 
	numFloatsPerCell(0),
        numPlanesToReserve(1),
	parallelDecomp(ROW),        
        noAddingPlanesDuringSimulation(false)
#if KMC_PARALLEL
      ,	latticeCommInitial(MPI_COMM_WORLD)
#endif
    {      
      globalPlanarDims[0] = globalPlanarDims[1] = ghostExtent[0] = ghostExtent[1] = 0;
      latInit = AddEmptyPlanes(1);
    }

    boost::array<int,2> globalPlanarDims /*! Array of length two,
					   indicating the number of
					   lattice cells along the
					   in-plane edges of the
					   domain of the lattice. */,

      ghostExtent /*! Array of length two, indicating the size (in
		    lattice cells) of the ghost region along each
		    in-plane edge of the lattice domain in a parallel
		    simulation. <STRONG>Has no effect in a serial
		    simulation.</STRONG> */;

    int numIntsPerCell /*! Size of the integer array at each
			 lattice cell. Defaults to zero. */,
      numFloatsPerCell /*! Size of the floating-point array at each
			 lattice cell. Defaults to zero. */;

    int numPlanesToReserve /*! Number of planes (or at least pointers
                               for those planes) for which to reserve
                               space in memory. Note that this is not
                               the number of planes actually added to
                               a lattice. Actually adding planes is
                               done by the Lattice::addPlanes() member
                               function. The relationship between this
                               parameter and Lattice::addPlanes() is
                               somewhat analogous to that of the
                               std::vector member functions reserve()
                               and push_back(). If this parameter is
                               too small, it is not an error, but it
                               may affect performance somewhat. */;

    ParallelDecomp parallelDecomp /*! Indicates the method of parallel
                                      decomposition in a parallel KMC
                                      simulation. <STRONG>Has no
                                      effect in a serial
                                      simulation.</STRONG>*/;

    LatticeInitializer latInit /*! Function or function object used to
				 initialize the lattice. Defaults to
				 an AddEmptyPlanes object that adds a
				 single lattice plane. */;

    bool noAddingPlanesDuringSimulation /*! When the preprocessor
                                           variable KMC_PARALLEL is
                                           non-zero, this indicates
                                           that no planes will be
                                           added during the
                                           simulation, in order to
                                           avoid unnecessary parallel
                                           communication when this is
                                           the case. Planes, then, may
                                           only be added during
                                           lattice
                                           initialization. <STRONG>Has no
                                           effect in a serial
                                           simulation.</STRONG> */;

#if KMC_PARALLEL
    MPI_Comm latticeCommInitial /*! When the preprocessor variable
				   KMC_PARALLEL is non-zero, this is
				   the MPI communicator used to
				   initialize an instance of the
				   Lattice class, and it defaults to
				   MPI_COMM_WORLD. <STRONG>Not available in
				   the serial version of the
				   ARL KMCThinFilm library.</STRONG> */;
#endif

    SetEmptyCellVals setEmptyCellVals /*! This is empty by default,
					and does not need to be
					defined if the integer and
					floating-point numbers in an
					empty lattice cell are all
					supposed to be zero. If set,
					this function is used to set
					the integers and
					floating-point numbers in an
					empty lattice cell to their
					proper values. Note that when
					this function is used, the
					vectors
					<VAR>emptyIntVals</VAR> and
					<VAR>emptyFloatVals</VAR> have
					already had their sizes set to
					the number of integer values
					and number of floating-point
					values per lattice cell,
					respectively.*/;
  };

  /*! A lattice that may be distributed over several processors.

    This lattice is designed to be simple but general. It is
    represented as a three-dimensional array, which can not only be
    used for square or simple cubic lattices, but also lattices of
    other physical shapes, such as the part of layer <VAR>k</VAR> of a
    simple hexagonal lattice shown below:

    \verbatim
      (0,1,k)-----(1,1,k)----(2,1,k)----(3,1,k)---(4,1,k)
        / \        / \        / \        / \        /
       /   \      /   \      /   \      /   \      /
      /     \    /     \    /     \    /     \    /
     /       \  /       \  /       \  /       \  /
    /_________\/_________\/_________\/_________\/
  (0,0,k)   (1,0,k)   (2,0,k)     (3,0,k)   (4,0,k)
    \endverbatim

    For further flexibility, each site may effectively contain an
    array of integers and/or an array of double-precision floating
    point numbers. The lengths of these arrays is the same for all
    sites on the lattice.

  */
  class Lattice : private boost::noncopyable {
    friend class Simulation;
  public:

    /*! [<STRONG>ADVANCED</STRONG>] Constructor of a Lattice object

      Unless one is debugging or testing an application that uses
      ARL KMCThinFilm, there is little point in calling this directly,
      since when a Simulation object is constructed, it builds a
      Lattice object internally for its own use.
     */
    Lattice(const LatticeParams & paramsForLattice);
    
    /*! Number of processors over which the lattice is distributed. */
    int nProcs() const;

    /*! When the preprocessor variable KMC_PARALLEL is zero, this
      always returns zero. Otherwise, this is the MPI rank of the
      local process, where the communicator for that process is given
      by Lattice::comm(). */
    int procID() const;

    /*! The MPI communicator used by the lattice for interchange of
        ghosts, if the preprocessor variable KMC_PARALLEL is
        non-zero. <STRONG>Not available in the serial version of the
        ARL KMCThinFilm library.</STRONG>

      This is <EM>not</EM> the same as the
      <VAR>latticeCommInitial</VAR> member of the LatticeParams object
      used to construct the lattice. While it contains the same number
      of processors as <VAR>latticeCommInitial</VAR>, the lines of
      code

      \code
      int rank;
      LatticeParams latParams;
      MPI_Comm_rank(latParams.latticeCommInitial, &rank);
      \endcode

      and

      \code
      // Lattice object named "lattice" defined above
      int rank;
      MPI_Comm_rank(lattice.comm(), &rank);
      \endcode

      may not return the same value for <VAR>rank</VAR>.
     */

#if KMC_PARALLEL
    const MPI_Comm & comm() const;
#endif

    /*! If the preprocessor variable KMC_PARALLEL is zero, then this
      always returns 1. Otherwise, this is the number of sectors (or
      sublattices) used in the approximate parallel Kinetic Monte
      Carlo algorithm used by this library.

      For a discussion of what these sectors are, see Y. Shim and
      J. G. Amar, "Semirigorous synchronous sublattice algorithm for
      parallel Kinetic Monte Carlo simulation of thin film growth",
      Physical Review B, vol. 71, 125432 (2005) or S. Plimpton et al.,
      "Crossing the Mesoscale No-Man's Land via Parallel Kinetic Monte
      Carlo", Sandia National Laboratories Technical Report
      SAND2009-6226 (2009).

    */
    int numSectors() const;

    /*! The number of processors along in-plane dimension <VAR>dim</VAR>.

      For example, if the lattice is partitioned as follows,

      \verbatim
      *-----*-----*-----*
      |  3  |  4  |  5  |
      *-----*-----*-----*
      |  0  |  1  |  2  |
      *-----*-----*-----*
      \endverbatim

      then procPerDim(0) returns a value of 3, and procPerDim(2)
      returns a value of 2.

      \see commCoord()
    */
    int procPerDim(int dim) const;

    /*! Coordinates for a processor in an MPI Cartesian topology
        along in-plane dimension <VAR>dim</VAR> for a given processor.

      For example, if the lattice is partitioned as follows,

      \verbatim
      *-----*-----*-----*
      |  3  |  4  |  5  |
      *-----*-----*-----*
      |  0  |  1  |  2  |
      *-----*-----*-----*
      \endverbatim

      where the numbers 0, 1, 2, etc., denote processor ranks, then if
      the MPI rank of a processor is 4, then commCoord(0) returns a
      value of 1, and commCoord(0) returns a value of 2. The
      coordinates for all processor ranks in the above partitioned
      lattice are shown below:

      \verbatim
      *-----*-----*-----*
      |(1,0)|(1,1)|(1,2)|
      *-----*-----*-----*
      |(0,0)|(0,1)|(0,2)|
      *-----*-----*-----*
      \endverbatim

      Note that the coordinates described here are <EM>not</EM> the
      same as the coordinates of a lattice cell.

      \see procPerDim()
     */
    int commCoord(int dim) const;

    /*! Thickness of the ghost region along in-plane dimension <VAR>dim</VAR>.

      Dimension 0 is the dimension along which the first lattice cell
      index varies, and dimension 1 is the dimension along which the
      second lattice cell index varies.

      \see commCoord()
     */
    int ghostExtent(int dim) const;

    /*! Adds zero or more layers to a given lattice, increasing the
        maximum value of the third lattice cell coordinate.

      Note that if <VAR>numPlanesToAdd</VAR> is non-positive, this
      member function does nothing. This behavior may be used to
      conditionally add lattice planes. For example,

      \code
      // Lattice object named "lattice" defined above

      CellInds ci;

      // Misc. calcs ...

      ci.k = someFunction(...);
      kMax = lattice.currHeight() - 1;

      lattice.addPlanes(ci.k - kMax); // Only adds a lattice plane if ci.k exceeds kMax

      lattice.setInt(ci, whichInt, anotherFunction(...));
      \endcode

     */
    void addPlanes(int numPlanesToAdd);

    /*! [<STRONG>ADVANCED</STRONG>] Reserves space in memory for
        adding planes to a lattice<SUP>*</SUP>, but does not actually
        add any additional planes. <STRONG>Almost certainly
        unnecessary if LatticeParams::numPlanesToReserve has
        been set to an appropriate value</STRONG>.

        The relationship between reservePlanes() and addPlanes() is
        somewhat analogous to that of the std::vector member functions
        reserve() and push_back(). The former only allocates space in
        memory, while the latter is responsible for appending actual
        values (i.e. lattice planes). Note that if
        <VAR>numTotalPlanesToReserve</VAR> is too small, it is not an
        error, though it may affect performance somewhat.

        <SUP>*</SUP>Or at least pointers for those planes.

        \see addPlanes()
     */
    void reservePlanes(int numTotalPlanesToReserve);

    /*! Returns the number of planes reserved for the lattice (but not
        necessarily added to the lattice yet).
       
       Note that while this should be at least the number of planes
       that were explicitly reserved, it may be more than that.
     */
    int planesReserved() const;

    /*! [<STRONG>ADVANCED</STRONG>] Update <EM>all</EM> the ghost
        lattice cells from data received from other processors. 

        This member function does not keep track of which lattice
        sites are actually changed during the update, unlike
        recvGhostsUpdate().
    */
    void recvGhosts(int sectNum /*!< Sector number, ranging from 0 to numSectors() - 1 */);

    /*! [<STRONG>ADVANCED</STRONG>] Send data from <EM>all</EM> the
        ghost lattice cells to update the off-process cells to which
        the ghosts correspond.

        This member function does not keep track of which lattice
        sites are actually changed during the update, unlike
        sendGhostsUpdate().
    */
    void sendGhosts(int sectNum /*!< Sector number, ranging from 0 to numSectors() - 1 */);

    /*! Obtain limiting values of the local lattice coordinates.

      "BBox" here is short for "bounding box."

      A typical use for this member function would be something like this:

      \code
      // Lattice object named "lattice" defined above

      bool includeGhosts = myBooleanFunc(...);

      LatticePlanarBBox localPlanarBBox;      
      lattice.getLocalPlanarBBox(includeGhosts, localPlanarBBox);
      kMaxP1 = lattice.currHeight();

      CellInds ci;
      for (ci.k = 0; ci.k < kMaxP1; ++(ci.k)) {
         for (ci.i = localPlanarBBox.imin; ci.i < localPlanarBBox.imaxP1; ++(ci.i)) {
            for (ci.j = localPlanarBBox.jmin; ci.j < localPlanarBBox.jmaxP1; ++(ci.j)) {
      
               // Some calcs ...
            
               lattice.setInt(ci, WHICH_INT, someFunc(...));               
               
               // Other calcs ...

           }
         }
      }
      \endcode
      
      \see getSectorPlanarBBox() getGlobalPlanarBBox()
     */
    void getLocalPlanarBBox(bool wGhost /*!< Indicates whether minimum
                                           and maximum values for
                                           lattice cell coordinates
                                           include the coordinates of
                                           ghost lattice cells. */,
                            LatticePlanarBBox & bbox /*!< Bounding box
                                                        of the local
                                                        lattice
                                                        coordinates. */) const;

    /*! Obtain limiting values of the local lattice coordinates for a
        particular sector or sublattice.

      "BBox" here is short for "bounding box."

      A typical use for this member function would be something like this:

      \code
      // Lattice object named "lattice" defined above

      for (int sectNum = 0; i < lattice.numSectors(); ++i) {

         lattice.recvGhosts(sectNum);
      
	 LatticePlanarBBox sectorPlanarBBox;
         lattice.getSectorPlanarBBox(sectNum, sectorPlanarBBox);
         kMaxP1 = lattice.currHeight();

         CellInds ci;
	 for (ci.k = 0; ci.k < kMaxP1; ++(ci.k)) {
            for (ci.i = sectorPlanarBBox.imin; ci.i < sectorPlanarBBox.imaxP1; ++(ci.i)) {
               for (ci.j = sectorPlanarBBox.jmin; ci.j < sectorPlanarBBox.jmaxP1; ++(ci.j)) {
                  // Some calcs ...
            
                  lattice.setFloat(ci, FLT_VAL1, someFunc(...));

                  // Other calcs ...

              }
            }
         }

         lattice.sendGhosts(sectNum);
      }
      \endcode

      Note that unlike getLocalPlanarBBox(), the minimum and maximum
      lattice coordinate values do not include the coordinates of
      ghost lattice cells.
      
      \see getLocalPlanarBBox() getGlobalPlanarBBox() numSectors()
    */
    void getSectorPlanarBBox(int sectNum /*!< Sector number, ranging from 0 to numSectors() - 1 */,
                             LatticePlanarBBox & bbox /*!< Bounding
                                                        box of the
                                                        lattice
                                                        coordinates
                                                        within the
                                                        sector. */) const;
    
    /*! Obtain limiting values of the global lattice coordinates.

      "BBox" here is short for "bounding box."

      Note that this function in general <EM>cannot</EM> be used to
      determine lattice coordinate values for use with getInt(),
      setInt(), getFloat(), and setFloat(), since those functions
      require local lattice coordinates. To obtain limiting values for
      local lattice coordinates, use getLocalPlanarBBox() or
      getSectorPlanarBBox() instead.

      \see getLocalPlanarBBox() getSectorPlanarBBox()
     */
    void getGlobalPlanarBBox(LatticePlanarBBox & bbox /*!< Bounding
                                                        box of the
                                                        global lattice
                                                        coordinates. */) const;

    /*! The current number of lattice planes.

      This is also one larger than the maximum value of the third lattice coordinate.

      \see getLocalPlanarBBox() getSectorPlanarBBox()
     */
    int currHeight() const;

    /*! Length of the integer array at each lattice cell. */
    int nIntsPerCell() const;

    /*! Length of the double-precision floating-point array at each lattice cell. */
    int nFloatsPerCell() const;

    /*! Returns the value of integer array element <VAR>whichInt</VAR> at lattice cell indices <VAR>ci</VAR>.

      The value of <VAR>whichInt</VAR> ranges from 0 to nIntsPerCell() - 1.

      \see getFloat() setFloat() setInt()
     */
    int getInt(const CellInds & ci, int whichInt) const;

    /*! Returns the value of double-precision array element <VAR>whichFloat</VAR> at lattice cell indices <VAR>ci</VAR>.

      The value of <VAR>whichFloat</VAR> ranges from 0 to nFloatsPerCell() - 1.

      \see getInt() setFloat() setInt()
     */
    double getFloat(const CellInds & ci, int whichFloat) const;

    /*! Sets to <VAR>val</VAR> the value of integer array element <VAR>whichInt</VAR> at lattice cell indices <VAR>ci</VAR>.

      The value of <VAR>whichInt</VAR> ranges from 0 to nIntsPerCell() - 1.

      \see getFloat() getInt() setFloat()
     */
    void setInt(const CellInds & ci, int whichInt, int val);

    /*! Sets to <VAR>val</VAR> the value of double-precision array element <VAR>whichFloat</VAR> at lattice cell indices <VAR>ci</VAR>.

      The value of <VAR>whichFloat</VAR> ranges from 0 to nFloatsPerCell() - 1.

      \see getFloat() getInt() setInt()
     */
    void setFloat(const CellInds & ci, int whichInt, double val);

    /*! Returns a wrapped version of ci.i to account for periodic
        boundary conditions.

      This is likely to only be useful in serial simulations.

      Furthermore, the member functions getInt(),
      getFloat(), setInt(), and setFloat() <EM>already</EM> account
      for periodic boundary conditions, so this function is mostly
      needed for circumstances where one is doing arithmetic of cell
      indices.
     */
    int wrapI(const CellInds & ci) const;

    /*! Returns a wrapped version of ci.j to account for periodic
        boundary conditions.

      This is likely to only be useful in serial simulations, though
      it perhaps may be of use if row-based decomposition is used in
      parallel simulations.

      Furthermore, the member functions getInt(),
      getFloat(), setInt(), and setFloat() <EM>already</EM> account
      for periodic boundary conditions, so this function is mostly
      needed for circumstances where one is doing arithmetic of cell
      indices.
     */
    int wrapJ(const CellInds & ci) const;

    /*! Returns the sector to which a set of cell indices <VAR>ci</VAR> belongs.
      
      If <VAR>ci</VAR> is within a ghost region, then this function
      returns -1.

      In serial, this always returns zero.      
    */
    int sectorOfIndices(const CellInds & ci) const;

    /*! [<STRONG>ADVANCED</STRONG>] Updates ghost lattice cells from
        data received from other processors, but only updates the
        cells that are supposed to have actually changed. Processors
        use addToExportBufferIfNeeded() to mark the cells that
        will be received by a call to this function. */
    void recvGhostsUpdate(int sectNum /*!< Sector number, ranging from 0 to numSectors() - 1 */);

    /*! [<STRONG>ADVANCED</STRONG>] Sends data from the ghost lattice
        cells to update the (usually off-process) cells to which the
        ghosts correspond, but only updates the cells that are
        supposed to have actually changed. Processors use
        addToExportBufferIfNeeded() to mark the cells that will be
        received by a call to this function. */
    void sendGhostsUpdate(int sectNum /*!< Sector number, ranging from 0 to numSectors() - 1 */);

    /*! [<STRONG>ADVANCED</STRONG>] Returns the lists of the indices
        of cells that are changed by the latest call to
        recvGhostsUpdate().*/
    const std::vector<std::vector<IJK> > & getReceivedGhostInds() const;

    /*! [<STRONG>ADVANCED</STRONG>] Returns the lists of the indices
        of cells that are changed by the latest call to
        sendGhostsUpdate(). */
    const std::vector<std::vector<IJK> > & getReceivedLocalInds() const;

    /*! [<STRONG>ADVANCED</STRONG>] Clears any ghosts that would be
        sent by a call to sendGhostsUpdate() */
    void clearGhostsToSend();

    /*! [<STRONG>ADVANCED</STRONG>] Adds <VAR>ci</VAR> to the lists of
        CellInds objects to be sent to other processors, provided that
        <VAR>ci</VAR> is either in a ghost region or in the part of
        the boundary of the lattice that is sent to other
        processors. Returns true if <VAR>ci</VAR> is in a ghost
        region.

        This should typically be used between a pair of calls to
        \link KMCThinFilm::Lattice::recvGhostsUpdate recvGhostsUpdate(i)\endlink 
        and \link KMCThinFilm::Lattice::sendGhostsUpdate sendGhostsUpdate(i)\endlink,
        where the value of <VAR>i</VAR> is the same for both calls. If
        ghosts are never exported, then clearGhostsToSend() should be
        used in place of sendGhostsUpdate(), and calls to
        addToExportBufferIfNeeded() may be called before any calls to
        recvGhostsUpdate().
    */
    bool addToExportBufferIfNeeded(const CellInds & ci);

    //! \cond HIDE_FROM_DOXYGEN
    ~Lattice(); // Destructor should be public so that "smart"
		// pointers to the Lattice class can delete it.

#if KMC_PARALLEL
    bool addToExportBufferIfNeeded(const CellInds & ci, int & sectNum);

    // This is only needed for debugging.
    int exportVal(const CellInds & ci) const;
#endif
    //! \endcond

  private:

    // Used in Simulation class

    //! \cond HIDE_FROM_DOXYGEN    
    struct TrackType {
      enum Type {
	NONE,
	CHECK_ONLY_IF_CHANGE_OCCURS,
	RECORD_CHANGED_CELL_INDS,
        RECORD_ONLY_OTHER_CHANGED_CELL_INDS,
      };
    };
    //! \endcond

    // Using a std::set turns out to be faster than a
    // boost::unordered_set for the purposes of this library (probably
    // because there are usually so few elements in ChangedCellInds).
    typedef std::set<CellInds> ChangedCellInds;

    typedef std::deque<CellInds> OtherCheckedCellInds;

    void trackChanges(TrackType::Type trackType);
    bool hasChanged() const;

    const ChangedCellInds & getChangedCellInds() const;
    const OtherCheckedCellInds & getOtherCheckedCellInds() const;

    void wrapIndsIfNeeded(CellInds & ci) const;

    class Impl_;
    boost::scoped_ptr<Impl_> pImpl_;
  };

}

#endif /* LATTICE_HPP */
