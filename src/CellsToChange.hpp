#ifndef CELLS_TO_CHANGE_HPP
#define CELLS_TO_CHANGE_HPP

#include <boost/scoped_ptr.hpp>
#include <vector>

namespace KMCThinFilm {

  class Lattice;
  class CellInds;
  class CellIndsOffset;

  /*! Class used to change the lattice when using a function object
      satisfying the EventExecutorSemiManualTrack signature to execute
      an event.

      Instances of this class contain offsets that define the relative
      positions of cells that would be directly changed by the
      event. These relative positions can be translated to actual
      positions by setting the <EM>center</EM> of a CellsToChange
      instance, that is, the absolute coordinates of the cell
      associated with the offset \f$(0,0,0)\f$ stored within the
      instance. For example, if the offsets are \f$(0,0,0)\f$,
      \f$(0,\pm1,0)\f$, \f$(\pm1,0,0)\f$, and \f$(0,0,\pm1)\f$, and
      the center is \f$(a,b,c)\f$, then the absolute coordinates of
      the cells changed by the event would be \f$(a,b,c)\f$,
      \f$(a,b\pm1,c)\f$, \f$(a\pm1,b,c)\f$, and
      \f$(a,b,c\pm1)\f$. Once the absolute coordinates of the cells
      are defined, the values in the lattice associated with these
      cells can be changed via setInt() and setFloat().

      The offsets are specified through the <VAR>cnoVec</VAR> argument
      in the following member functions:

      - EventExecutorGroup::addEventExecutor(int whichEvent, EventExecutorSemiManualTrack evExec, const std::vector<CellNeighOffsets> & cnoVec)
      - Simulation::addOverLatticeEvent(int, double, EventExecutorSemiManualTrack evExec, const std::vector<CellNeighOffsets> & cnoVec)

      The offsets specified in the first element of <VAR>cnoVec</VAR>
      correspond to the offsets stored in the first element of the
      argument of type std::vector<CellsToChange>& in a
      EventExecutorSemiManualTrack function object, the offsets
      specified in the second element of <VAR>cnoVec</VAR> correspond
      to the offsets stored in the second element of the argument of
      type std::vector<CellsToChange>& in that same function object,
      etc.

   */
  class CellsToChange {
    friend class Simulation;
  public:

    /*! Sets to <VAR>ci</VAR> the absolute coordinates of the cell
        associated with the offset \f$(0,0,0)\f$ stored within a
        CellsToChange object. */
    void setCenter(const CellInds & ci);

    /*! Returns the absolute coordinates that become associated with
        the offset ID <VAR>whichOffset</VAR> once setCenter() has
        been called. 

        Obviously, setCenter() <STRONG>MUST</STRONG> be called before
        this function.
    */
    const CellInds & getCellInds(int whichOffset) const;

    /*! Returns the value of integer array element <VAR>whichInt</VAR>
        at the lattice cell associated with offset ID
        <VAR>whichOffset</VAR>.

        Given a lattice object instance <TT>lattice</TT> and a
        CellsToChange instance <TT>ctc</TT>, this function is
        essentially equivalent to
        <TT>lattice.getInt(ctc.getCellInds(whichOffset),
        whichInt)</TT>

        setCenter() <STRONG>MUST</STRONG> be called before this
        function.

      \see Lattice::getInt()
     */
    int getInt(int whichOffset, int whichInt) const;

     /*! Returns the value of double-precision array element
        <VAR>whichFloat</VAR> at the lattice cell associated with
        offset ID <VAR>whichOffset</VAR>.

        Given a lattice object instance <TT>lattice</TT> and a
        CellsToChange instance <TT>ctc</TT>, this function is
        essentially equivalent to
        <TT>lattice.getFloat(ctc.getCellInds(whichOffset),
        whichFloat)</TT>.

        setCenter() <STRONG>MUST</STRONG> be called before this
        function.

      \see Lattice::getFloat()
     */
    double getFloat(int whichOffset, int whichFloat) const;

    /*! Sets the value of integer array element <VAR>whichInt</VAR> at
        lattice cell associated with offset ID <VAR>whichOffset</VAR>
        to <VAR>val</VAR>.

        setCenter() <STRONG>MUST</STRONG> be called before this
        function.

      \see Lattice::setInt()
    */
    void setInt(int whichOffset, int whichInt, int val);

    /*! Sets the value of double-precision array element
        <VAR>whichFloat</VAR> at lattice cell associated with offset
        ID <VAR>whichOffset</VAR> to <VAR>val</VAR>.

        setCenter() <STRONG>MUST</STRONG> be called before this
        function.

      \see Lattice::setFloat()
    */
    void setFloat(int whichOffset, int whichFloat, double val);

    /*! Adds zero or more layers to a given lattice, increasing the
        maximum value of the third lattice cell coordinate.

        Note that if <VAR>numPlanesToAdd</VAR> is non-positive, this
        member function does nothing.

        \see Lattice::addPlanes()
     */
    void addLatticePlanes(int numPlanesToAdd);

    //! \cond HIDE_FROM_DOXYGEN
    CellsToChange(const CellsToChange & ctc);
    CellsToChange & operator=(const CellsToChange & rhs);
    ~CellsToChange();
    //! \endcond

  private:
    
    // For use by Simulation class
    explicit CellsToChange(Lattice * lattice, int numOffsets);
    void addOffset(const CellIndsOffset & offset);
    const CellInds & getCenter() const;
    const std::vector<CellIndsOffset> & getCellIndsOffsetVec() const;
    const std::vector<CellInds> & getCellIndsVec() const;
    void clear();

    class Impl_;
    boost::scoped_ptr<Impl_> pImpl_;
  };

}

#endif /* CELLS_TO_CHANGE_HPP */
