#include "Lattice.hpp"
#include "ErrorHandling.hpp"
#include "CallMemberFunction.hpp"

#include "wrapInd.hpp"

#include <deque>
#include <cassert>
#include <stdexcept>
#include <algorithm>

#include <boost/multi_array.hpp>
#include <boost/lexical_cast.hpp>

using namespace KMCThinFilm;

struct Lattice::Impl_ {

  Impl_(const LatticeParams & paramsForLattice,
	Lattice * const self);

  Lattice * const self_; // Back-pointer to actual Lattice class
  int nIntsPerCell_, nFloatsPerCell_;
  boost::array<int,2> globalPlanarDims_, ghostExtent_;
  SetEmptyCellVals setEmptyCellVals_;
  LatticeParams::ParallelDecomp parallelDecomp_;

  boost::array<int,2> commCoords_, procsPerDim_,
    globalOffset_, globalOffsetMinusGhostExtent_, extentWGhost_;

  typedef boost::multi_array<int, 3> IntArray_;
  typedef boost::multi_array<double, 3> FloatArray_;

  struct IntFloatArrayPair_ {
    IntArray_ intArray_;
    FloatArray_ floatArray_;

    IntFloatArrayPair_(const boost::array<int,2> & extentWGhost,
		       const boost::array<int,2> & globalOffsetMinusGhostExtent, 
		       int nIntsPerCell, int nFloatsPerCell);
  };
  
  std::vector<IntFloatArrayPair_> lattice_;

  typedef void (Impl_::*SetInt_)(const CellInds & ci, int whichInt, int val);
  typedef void (Impl_::*SetFloat_)(const CellInds & ci, int whichFloat, double val);

  typedef void (Impl_::*AppendPlane_)();

  SetInt_ setInt_;
  SetFloat_ setFloat_;

  AppendPlane_ appendPlaneOnly_; // This points to either
				 // appendPlaneNoExplicitEmpty_ or
				 // appendPlaneWithExplicitEmpty_.

  AppendPlane_ appendPlane_; // This points to appendPlaneOnly_,
			     // appendPlaneAndRecordThatChangeOccurred_,
			     // or
			     // appendPlaneAndRecordChangedCellInds_.

  bool latticeModified_;
  ChangedCellInds changedCellInds_;  
  OtherCheckedCellInds otherCheckedCellInds_;

  int nProcs_, procID_;
  
  /* These functions are here to minimize the dependence on the types
     of ChangedCellInds and OtherCheckedCellInds. If the type of
     ChangedCellInds or OtherCheckedCellInds changes (e.g. to a deque
     or a set or a vector), then only these functions and a typedef or
     two in Lattice.hpp may need to change. */
  void addToChangedCellInds_(const CellInds & ci) {changedCellInds_.insert(ci);} // for set
  void addToOtherCheckedCellInds_(const CellInds & ci) {otherCheckedCellInds_.push_back(ci);}
  
  void setIntOnly_(const CellInds & ci, int whichInt, int val);
  void setIntAndRecordThatChangeOccurred_(const CellInds & ci, int whichInt, int val);
  void setIntAndRecordChangedCellInds_(const CellInds & ci, int whichInt, int val);

  void setFloatOnly_(const CellInds & ci, int whichFloat, double val);
  void setFloatAndRecordThatChangeOccurred_(const CellInds & ci, int whichFloat, double val);
  void setFloatAndRecordChangedCellInds_(const CellInds & ci, int whichFloat, double val);

  void addPlanes_(int numPlanesToAdd);

  void appendPlaneNoExplicitEmpty_();
  void appendPlaneWithExplicitEmpty_();
  void appendPlaneFake_();

  void appendPlaneAndRecordThatChangeOccurred_();
  void appendPlaneAndRecordChangedCellInds_();

#if !KMC_PARALLEL
  void wrapBothInds_(int & i, int & j) const {
    i = wrapInd(i, globalPlanarDims_[0]);
    j = wrapInd(j, globalPlanarDims_[1]);
  }
#endif

#if KMC_PARALLEL
  MPI_Datatype ijkDType_;

  typedef int (Impl_::*SetNewLatticeHeight_)(int currLocalHeight);
  SetNewLatticeHeight_ setNewLatticeHeight_;

  int setNewLatticeHeightActual_(int currLocalHeight);
  int setNewLatticeHeightFake_(int currLocalHeight);

  typedef void (Impl_::*WrapInds_)(int & i, int & j) const;

  WrapInds_ wrapIndsIfNeeded_;
  
  void noWrap_(int & i, int & j) const {}

  void wrapJOnly_(int & i, int & j) const {
    j = wrapInd(j, globalPlanarDims_[1]);
  }

  // For 1-D decomposition
  struct ProcIdCart1_ {
    enum Type {
      Back, Fwd, SIZE
    };
  };

  struct ExportVal1_ {
    enum Type {
      NONE, Back, Fwd, SIZE
    };
  };
  
  // For 2-D decomposition
  struct ProcIdCart2_ {
    enum Type {
      BackRow, FwdRow, FwdCol, BackCol,
      BackRowBackCol, FwdRowBackCol, BackRowFwdCol, FwdRowFwdCol,
      SIZE
    };
  };

  struct ExportVal2_ {
    enum Type {
      NONE,
      BackRow1,
      FwdRow1,
      BackRow2,
      FwdRow2,
      BackRowOverlap,
      FwdRowOverlap,
      BackCol1,
      FwdCol1,
      BackCol2,
      FwdCol2,
      BackColOverlap,
      FwdColOverlap,
      BackRowBackCol,
      FwdRowFwdCol,
      BackRowFwdCol,
      FwdRowBackCol,
      SIZE
    };
  };

  // For row decomposition

  struct BufferTypeRow_ {
    enum Type {ROW, SIZE};
  };

  int extentInt_, extentFloat_;
  
  struct SendRecvGhostRowInfo_ {
    int sendRank, recvRank, sendPlanarCoords, recvPlanarCoords;
  };
  
  std::vector<SendRecvGhostRowInfo_> sendRecvGhostRowInfo_;

  void setUpMPICartForRowPartitions_(MPI_Comm latticeCommInitial, std::vector<int> & pIdNeighs);
  void setUpSectoringForRowPartitions_(const std::vector<int> & pIdNeighs);

  void recvGhostsRow_(int sectNum, int newLatticeHeight);
  void sendGhostsRow_(int sectNum, int newLatticeHeight);

  void recvGhostsUpdateRow_(int sectNum);
  void sendGhostsUpdateRow_(int sectNum);

  // For compact decomposition

  struct BufferTypeCompact_ {
    enum Type {ROW, COL, CORNER, SIZE};
  };

  MPI_Datatype ghostCornerInt_, ghostCornerFloat_;
  boost::array<MPI_Datatype,2> ghostHalfRowInt_, ghostHalfColInt_,
    ghostHalfRowFloat_, ghostHalfColFloat_;

  struct SendRecvGhostCompactInfo_ {
    int sendCornerRank, recvCornerRank,
      sendHalfRowRank, recvHalfRowRank,
      sendHalfColRank, recvHalfColRank;

    boost::array<int,2> sendCornerPlanarCoords, recvCornerPlanarCoords,
      sendHalfRowPlanarCoords, recvHalfRowPlanarCoords,
      sendHalfColPlanarCoords, recvHalfColPlanarCoords;

    int whichHalfRow, whichHalfCol;
  };

  std::vector<SendRecvGhostCompactInfo_> sendRecvGhostCompactInfo_;
  
  struct ProcsPerDim_ {
    int nprocs_x, nprocs_y;
    ProcsPerDim_(int x, int y) : nprocs_x(x), nprocs_y(y) {}
  };

  double calcHalfPerimeter_(const ProcsPerDim_ & ppd, double dim_x, double dim_y) const;
  void getPossProcsPerDim_(int nProcs, std::deque<ProcsPerDim_> & dequeOfFactors) const;
  void findProcsPerDimCompact_(int nProcs, double dim_x, double dim_y, boost::array<int,2> & ppd) const;

  void setUpMPICartForCompactPartitions_(MPI_Comm latticeCommInitial, std::vector<int> & pIdNeighs);
  void setUpSectoringForCompactPartitions_(const std::vector<int> & pIdNeighs);

  void recvGhostsCompact_(int sectNum, int newLatticeHeight);
  void sendGhostsCompact_(int sectNum, int newLatticeHeight);

  void recvGhostsUpdateCompact_(int sectNum);
  void sendGhostsUpdateCompact_(int sectNum);

  // For sector decomposition in general

  std::vector<int> iminS_, imaxP1S_, jminS_, jmaxP1S_;

  int nSectors_;

  MPI_Comm latticeCommCart_;

  boost::array<int,2> localDims_;

  void getLocalPlanarBBox_(bool wGhost,
			   int & imin, int & imaxP1, 
			   int & jmin, int & jmaxP1) const;

  typedef void (Impl_::*SendRecvGhosts_)(int sectNum, int newLatticeHeight);
  typedef void (Impl_::*SendRecvGhostsUpdate_)(int sectNum);

  SendRecvGhosts_ sendGhosts_, recvGhosts_;
  SendRecvGhostsUpdate_ sendGhostsUpdate_, recvGhostsUpdate_;

  // Buffers use IJK rather than CellInds because only the former is a POD
  std::vector<std::vector<IJK> > ghostSendIndsBuffer_, ghostRecvIndsBuffer_, localRecvIndsBuffer_;
  std::vector<std::vector<int> > ghostSendIntBuffer_, ghostRecvIntBuffer_, localRecvIntBuffer_;
  std::vector<std::vector<double> > ghostSendFloatBuffer_, ghostRecvFloatBuffer_, localRecvFloatBuffer_;

  std::vector<std::vector<std::vector<IJK> > > localSendIndsBuffer_;
  std::vector<std::vector<std::vector<int> > > localSendIntBuffer_;
  std::vector<std::vector<std::vector<double> > > localSendFloatBuffer_;

  std::vector<std::vector<LatticePlanarBBox> > ghostRecvBufferBounds_, localRecvBufferBounds_;

  struct CellParInfo_ {
    // Using "signed char" rather than "int" to save some space.
    signed char sectNum, exportVal;
  };

  boost::multi_array<CellParInfo_,2> cellParInfoArray_;

  struct LocalBufferLocs_ {
    std::size_t nBuffers;
    enum { nBuffersMax = 3 };
    boost::array<int,nBuffersMax> sectNum, bufType;
  };

  std::vector<signed char> exportValToGhostBufType_;
  std::vector<LocalBufferLocs_> exportValToLocalBufferLocs_;

  void setLatticeValsFromRecvBuffers_(std::vector<std::vector<IJK> > & recvIndsBuffer,
                                      const std::vector<std::vector<int> > & recvIntBuffer,
                                      const std::vector<std::vector<double> > & recvFloatBuffer,
                                      const std::vector<LatticePlanarBBox> & recvBufferBounds);

  void setSendIntFloatBuffers_(const std::vector<std::vector<IJK> > & sendIndsBuffer,
                               std::vector<std::vector<int> > & sendIntBuffer,
                               std::vector<std::vector<double> > & sendFloatBuffer);

  typedef void (Impl_::*ReExportIfNeeded_)();
  ReExportIfNeeded_ reExportIfNeeded_;

  void reExportIfNeededActual_();
  void reExportIfNeededFake_() {}

#else
  // Used for getReceivedGhostInds() and getReceivedLocalInds() in serial mode.
  std::vector<std::vector<IJK> > dummyNullVector_;
#endif

};

Lattice::Impl_::Impl_(const LatticeParams & paramsForLattice,
		      Lattice * const self)
  : self_(self),
    nIntsPerCell_(paramsForLattice.numIntsPerCell),
    nFloatsPerCell_(paramsForLattice.numFloatsPerCell),
    globalPlanarDims_(paramsForLattice.globalPlanarDims),
    ghostExtent_(paramsForLattice.ghostExtent),
    setEmptyCellVals_(paramsForLattice.setEmptyCellVals),
    parallelDecomp_(paramsForLattice.parallelDecomp)
 {
  
  exitOnCondition((nIntsPerCell_ < 1) && (nFloatsPerCell_ < 1),
                  "There must be at least one integer or floating point number defined per lattice site.");

  exitOnCondition((globalPlanarDims_[0] <= 0) || (globalPlanarDims_[1] <= 0),
		  "One of the in-plane dimensions of the lattice is less than or equal to zero (or was not set).");

  if (setEmptyCellVals_) {
    appendPlaneOnly_ = &Impl_::appendPlaneWithExplicitEmpty_;
  }
  else {
    appendPlaneOnly_ = &Impl_::appendPlaneNoExplicitEmpty_;
  }
  // Note: After paramsForLattice.latInit() has been run,
  // appendPlaneOnly_ may be changed to appendPlaneFake_.

  assert(!(paramsForLattice.numPlanesToReserve < 0));
  // Need to call this before any instance of lattice_.capacity() is
  // called.
  lattice_.reserve(paramsForLattice.numPlanesToReserve);

#if KMC_PARALLEL
  setNewLatticeHeight_ = &Impl_::setNewLatticeHeightActual_;
  // Note: After paramsForLattice.latInit() has been run,
  // setNewLatticeHeight_ may be changed to setNewLatticeHeightFake_.

  MPI_Comm_size(paramsForLattice.latticeCommInitial, &nProcs_);

  std::vector<int> pIdNeighs;

  switch (paramsForLattice.parallelDecomp) {
  case LatticeParams::COMPACT:
    setUpMPICartForCompactPartitions_(paramsForLattice.latticeCommInitial, pIdNeighs);
    break;
  case LatticeParams::ROW:
    ghostExtent_[1] = 0; // Need to define this before calculating
			 // extentWGhost_[1].
    setUpMPICartForRowPartitions_(paramsForLattice.latticeCommInitial, pIdNeighs);
    break;
  }
#else
  nProcs_ = 1;
  procID_ = 0;
  procsPerDim_[0] = procsPerDim_[1] = 1;
  commCoords_[0] = commCoords_[1] = 0;

  ghostExtent_[0] = ghostExtent_[1] = 0;
#endif
  
  // Need to define localDims_ before setting up sectoring
  for (std::size_t i = 0; i < 2; ++i) {
#if KMC_PARALLEL
    double localDimUsual = globalPlanarDims_[i]/procsPerDim_[i];

    localDims_[i] = (commCoords_[i] == (procsPerDim_[i] - 1) ? 
                     localDimUsual + (globalPlanarDims_[i] % procsPerDim_[i]) :
                     localDimUsual);

    globalOffset_[i] = (commCoords_[i])*localDimUsual;
    extentWGhost_[i] = localDims_[i] + 2*ghostExtent_[i];
    globalOffsetMinusGhostExtent_[i] = globalOffset_[i] - ghostExtent_[i];
#else
    globalOffset_[i] = globalOffsetMinusGhostExtent_[i] = 0;
    extentWGhost_[i] = globalPlanarDims_[i];
#endif
  }

#if KMC_PARALLEL

  // cellParInfoArray_ must be sized and initialized before sectoring
  // is set up. The extra padding is so that sectorOfIndices can be
  // used on some indices that are outside a ghost region but still
  // may be probed by functions in the Simulation class().
  cellParInfoArray_.resize(boost::extents[localDims_[0] + 4*ghostExtent_[0]][localDims_[1] + 4*ghostExtent_[1]]);
  boost::array<int,2> globalOffsetMinusTwiceGhostExtent;
  for (std::size_t i = 0; i < 2; ++i) {
    globalOffsetMinusTwiceGhostExtent[i] = globalOffset_[i] - 2*ghostExtent_[i];
  }
  cellParInfoArray_.reindex(globalOffsetMinusTwiceGhostExtent);  
  CellParInfo_ dummyCellParInfo = {-2, -2};
  std::fill(cellParInfoArray_.data(),
            cellParInfoArray_.data() + cellParInfoArray_.num_elements(),
            dummyCellParInfo);

  switch (paramsForLattice.parallelDecomp) {
  case LatticeParams::COMPACT:
    setUpSectoringForCompactPartitions_(pIdNeighs);
    break;
  case LatticeParams::ROW:
    setUpSectoringForRowPartitions_(pIdNeighs);
    break;
  }

  int blockcount = 3;
  MPI_Aint displ = 0;
  MPI_Datatype type = MPI_INT;
  MPI_Datatype ijkDTypeTmp;

  MPI_Type_create_struct(1, &blockcount, &displ, &type,
                         &ijkDTypeTmp);
  MPI_Type_commit(&ijkDTypeTmp);
    
  // Ensuring that ijkDType_ has the correct size when used to
  // communicate an array of CellInds.
  MPI_Type_create_resized(ijkDTypeTmp, 0, sizeof(CellInds),
                          &ijkDType_);
  MPI_Type_commit(&ijkDType_);
  MPI_Type_free(&ijkDTypeTmp);
#endif
}

Lattice::Impl_::IntFloatArrayPair_::IntFloatArrayPair_(const boost::array<int,2> & extentWGhost,
						       const boost::array<int,2> & globalOffsetMinusGhostExtent, 
						       int nIntsPerCell, int nFloatsPerCell) 
  : intArray_(boost::extents[0][0][0], boost::c_storage_order()),
    floatArray_(boost::extents[0][0][0], boost::c_storage_order()) {

  boost::array<int,3> localStart;
  localStart[0] = globalOffsetMinusGhostExtent[0]; 
  localStart[1] = globalOffsetMinusGhostExtent[1];
  localStart[2] = 0;

  // According to the documentation for Boost MultiArray, the resize()
  // member function of intArray_ and floatArray_ should initialize
  // any new elements of these arrays with their respective default
  // constructors int() and double(), so these arrays should already
  // be initialized to zero.

  if (nIntsPerCell > 0) {
    intArray_.resize(boost::extents[extentWGhost[0]][extentWGhost[1]][nIntsPerCell]);
    intArray_.reindex(localStart);
  }

  if (nFloatsPerCell > 0) {
    floatArray_.resize(boost::extents[extentWGhost[0]][extentWGhost[1]][nFloatsPerCell]);
    floatArray_.reindex(localStart);
  }
}

void Lattice::Impl_::setIntOnly_(const CellInds & ci, int whichInt, int val) {
  int iWrapped = ci.i;
  int jWrapped = ci.j;

#if KMC_PARALLEL
  KMC_CALL_MEMBER_FUNCTION(*this, wrapIndsIfNeeded_)(iWrapped, jWrapped);
#else
  wrapBothInds_(iWrapped, jWrapped);
#endif

  lattice_[ci.k].intArray_[iWrapped][jWrapped][whichInt] = val;
}

void Lattice::Impl_::setFloatOnly_(const CellInds & ci, int whichFloat, double val) {
  int iWrapped = ci.i;
  int jWrapped = ci.j;

#if KMC_PARALLEL
  KMC_CALL_MEMBER_FUNCTION(*this, wrapIndsIfNeeded_)(iWrapped, jWrapped);
#else
  wrapBothInds_(iWrapped, jWrapped);
#endif

  lattice_[ci.k].floatArray_[iWrapped][jWrapped][whichFloat] = val;
}

void Lattice::Impl_::setIntAndRecordThatChangeOccurred_(const CellInds & ci, int whichInt, int val) {
  latticeModified_ = true;
  setIntOnly_(ci, whichInt, val);
}

void Lattice::Impl_::setFloatAndRecordThatChangeOccurred_(const CellInds & ci, int whichFloat, double val) {
  latticeModified_ = true;
  setFloatOnly_(ci, whichFloat, val);
}

void Lattice::Impl_::setIntAndRecordChangedCellInds_(const CellInds & ci, int whichInt, int val) {
  setIntAndRecordThatChangeOccurred_(ci, whichInt, val);
  addToChangedCellInds_(ci);
}

void Lattice::Impl_::setFloatAndRecordChangedCellInds_(const CellInds & ci, int whichFloat, double val) {
  setFloatAndRecordThatChangeOccurred_(ci, whichFloat, val);
  addToChangedCellInds_(ci);
}

void Lattice::Impl_::addPlanes_(int numPlanesToAdd) {
  for (int i = 0; i < numPlanesToAdd; ++i) {
    KMC_CALL_MEMBER_FUNCTION(*this, appendPlane_)();
  }
}

void Lattice::Impl_::appendPlaneFake_() {
  abortWithMsg("Additional planes are not to be added during the simulation when noAddingPlanesDuringSimulation = true.");
}

void Lattice::Impl_::appendPlaneNoExplicitEmpty_() {
  lattice_.push_back(IntFloatArrayPair_(extentWGhost_, globalOffsetMinusGhostExtent_,
					nIntsPerCell_, nFloatsPerCell_));
}

void Lattice::Impl_::appendPlaneWithExplicitEmpty_() {

  // Note that the third argument is the height of the lattice
  // *before* adding a new lattice plane, i.e., the maximum value that
  // ci.k should have.
  CellInds ci(0,0, lattice_.size());

  appendPlaneNoExplicitEmpty_();
  
  int imin, imaxP1, jmin, jmaxP1;
#if KMC_PARALLEL
  // imin, imaxP1, jmin, jmaxP1 include ghosts.
  getLocalPlanarBBox_(true, imin, imaxP1, jmin, jmaxP1);
#else
  imin = jmin = 0;
  imaxP1 = globalPlanarDims_[0];
  jmaxP1 = globalPlanarDims_[1];
#endif

  std::vector<int> emptyIntVals(nIntsPerCell_);
  std::vector<double> emptyFloatVals(nFloatsPerCell_);

  for (ci.i = imin; ci.i < imaxP1; ++(ci.i)) {
    for (ci.j = jmin; ci.j < jmaxP1; ++(ci.j)) {

      setEmptyCellVals_(ci, *self_, emptyIntVals, emptyFloatVals);
      assert(nIntsPerCell_ == static_cast<int>(emptyIntVals.size()));
      assert(nFloatsPerCell_ == static_cast<int>(emptyFloatVals.size()));

      for (int whichInt = 0; whichInt < nIntsPerCell_; ++whichInt) {
	lattice_.back().intArray_[ci.i][ci.j][whichInt] = emptyIntVals[whichInt];
      }

      for (int whichFloat = 0; whichFloat < nFloatsPerCell_; ++whichFloat) {
	lattice_.back().floatArray_[ci.i][ci.j][whichFloat] = emptyFloatVals[whichFloat];
      }

    }
  }

}

void Lattice::Impl_::appendPlaneAndRecordThatChangeOccurred_() {
  latticeModified_ = true;
  KMC_CALL_MEMBER_FUNCTION(*this, appendPlaneOnly_)();
}

void Lattice::Impl_::appendPlaneAndRecordChangedCellInds_() {

  // Note that the third argument is the height of the lattice
  // *before* adding a new lattice plane, i.e., the maximum value that
  // ci.k should have.
  CellInds ci(0,0, lattice_.size());

  appendPlaneAndRecordThatChangeOccurred_();

  int imin, imaxP1, jmin, jmaxP1;
#if KMC_PARALLEL
  // imin, imaxP1, jmin, jmaxP1 do NOT include ghosts.
  getLocalPlanarBBox_(false, imin, imaxP1, jmin, jmaxP1);
#else
  imin = jmin = 0;
  imaxP1 = globalPlanarDims_[0];
  jmaxP1 = globalPlanarDims_[1];
#endif

  for (ci.i = imin; ci.i < imaxP1; ++(ci.i)) {
    for (ci.j = jmin; ci.j < jmaxP1; ++(ci.j)) {
      addToOtherCheckedCellInds_(ci);
    }
  }

}

#if KMC_PARALLEL
double Lattice::Impl_::calcHalfPerimeter_(const ProcsPerDim_ & ppd, double dim_x, double dim_y) const {
  return dim_x/ppd.nprocs_x + dim_y/ppd.nprocs_y;
}
#endif

#if KMC_PARALLEL
void Lattice::Impl_::getPossProcsPerDim_(int nProcs, std::deque<ProcsPerDim_> & dequeOfFactors) const {
  for (int i = 1; i <= nProcs; ++i) {
    if ((nProcs % i) == 0) {
      dequeOfFactors.push_back(ProcsPerDim_(i, nProcs/i));
    }
  }
}
#endif

#if KMC_PARALLEL
void Lattice::Impl_::findProcsPerDimCompact_(int nProcs, double dim_x, double dim_y, boost::array<int,2> & ppd) const {

  std::deque<ProcsPerDim_> dequeOfFactors;
  getPossProcsPerDim_(nProcs, dequeOfFactors);

  int dequeOfFactorsSize = dequeOfFactors.size();
  
  double minHalfPerim = calcHalfPerimeter_(dequeOfFactors[0], dim_x, dim_y);
  int procPerDimIndexAtMin = 0;

  // Starting the iteration at 1, *not* zero  
  for (int i = 1; i < dequeOfFactorsSize; ++i) {
    double halfPerim = calcHalfPerimeter_(dequeOfFactors[i], dim_x, dim_y);

    if (halfPerim < minHalfPerim) {
      minHalfPerim = halfPerim;
      procPerDimIndexAtMin = i;
    }
  }

  ppd[0] = dequeOfFactors[procPerDimIndexAtMin].nprocs_x;
  ppd[1] = dequeOfFactors[procPerDimIndexAtMin].nprocs_y;
  
}
#endif

#if KMC_PARALLEL
int Lattice::Impl_::setNewLatticeHeightActual_(int currLocalHeight) {
  
  int latticeHeightMax;

  // Make sure that the portions of the lattice on all processors have
  // the same number of layers.
  MPI_Allreduce(&currLocalHeight, &latticeHeightMax, 1,
                MPI_INT, MPI_MAX, latticeCommCart_);

  addPlanes_(latticeHeightMax - currLocalHeight); /* If
						    latticeHeightMax
						    minus
						    latticeHeight
						    equals 0, no new
						    planes are
						    added. */

  return latticeHeightMax;
}

int Lattice::Impl_::setNewLatticeHeightFake_(int currLocalHeight) {
  // This is pretty much a no-op
  return currLocalHeight;
}

#endif

#if KMC_PARALLEL
void Lattice::Impl_::setUpMPICartForCompactPartitions_(MPI_Comm latticeCommInitial, std::vector<int> & pIdNeighs) {
  findProcsPerDimCompact_(nProcs_, globalPlanarDims_[0], globalPlanarDims_[1], procsPerDim_);
  
  std::vector<int> isPeriodic(2,1);
  MPI_Cart_create(latticeCommInitial, 2, &procsPerDim_[0],
		  &isPeriodic[0], 1, &latticeCommCart_);

  MPI_Comm_rank(latticeCommCart_, &procID_);
  MPI_Cart_coords(latticeCommCart_, procID_, 2, &commCoords_[0]);

  pIdNeighs.resize(ProcIdCart2_::SIZE);

  MPI_Cart_shift(latticeCommCart_, 0, 1, &pIdNeighs[ProcIdCart2_::BackRow], &pIdNeighs[ProcIdCart2_::FwdRow]);
  MPI_Cart_shift(latticeCommCart_, 1, 1, &pIdNeighs[ProcIdCart2_::BackCol], &pIdNeighs[ProcIdCart2_::FwdCol]);

  boost::array<int,2> backRowBackCol(commCoords_), fwdRowBackCol(commCoords_),
    backRowFwdCol(commCoords_), fwdRowFwdCol(commCoords_);

  --(backRowBackCol[0]); --(backRowBackCol[1]);
  ++(fwdRowBackCol[0]); --(fwdRowBackCol[1]);
  --(backRowFwdCol[0]); ++(backRowFwdCol[1]);
  ++(fwdRowFwdCol[0]); ++(fwdRowFwdCol[1]);
  
  MPI_Cart_rank(latticeCommCart_, &backRowBackCol[0], &pIdNeighs[ProcIdCart2_::BackRowBackCol]);
  MPI_Cart_rank(latticeCommCart_, &fwdRowBackCol[0], &pIdNeighs[ProcIdCart2_::FwdRowBackCol]);
  MPI_Cart_rank(latticeCommCart_, &backRowFwdCol[0], &pIdNeighs[ProcIdCart2_::BackRowFwdCol]);
  MPI_Cart_rank(latticeCommCart_, &fwdRowFwdCol[0], &pIdNeighs[ProcIdCart2_::FwdRowFwdCol]);  
}
#endif

#if KMC_PARALLEL
void Lattice::Impl_::setUpSectoringForCompactPartitions_(const std::vector<int> & pIdNeighs) {
  nSectors_ = 4;

  wrapIndsIfNeeded_ = &Impl_::noWrap_;

  sendGhosts_ = &Impl_::sendGhostsCompact_;
  recvGhosts_ = &Impl_::recvGhostsCompact_;

  sendGhostsUpdate_ = &Impl_::sendGhostsUpdateCompact_;
  recvGhostsUpdate_ = &Impl_::recvGhostsUpdateCompact_;

  reExportIfNeeded_ = &Impl_::reExportIfNeededActual_;

  iminS_.resize(nSectors_);
  imaxP1S_.resize(nSectors_);
  jminS_.resize(nSectors_);
  jmaxP1S_.resize(nSectors_);

  ghostRecvBufferBounds_.resize(nSectors_);
  localRecvBufferBounds_.resize(nSectors_);

  for (int i = 0; i < nSectors_; ++i) {
    ghostRecvBufferBounds_[i].resize(BufferTypeCompact_::SIZE);
    localRecvBufferBounds_[i].resize(BufferTypeCompact_::SIZE);
  }
      
  boost::array<int,2> firstHalfLocalDims, secondHalfLocalDims, maxHalfLocalDims;

  for (std::size_t i = 0; i < 2; ++i) {
    firstHalfLocalDims[i] = localDims_[i]/2;
    secondHalfLocalDims[i] = localDims_[i] - firstHalfLocalDims[i];

    maxHalfLocalDims[i] = ((firstHalfLocalDims[i] > secondHalfLocalDims[i])
                           ? firstHalfLocalDims[i]
                           :secondHalfLocalDims[i]);
  }

  ghostSendIndsBuffer_.resize(BufferTypeCompact_::SIZE);
  ghostSendIntBuffer_.resize(BufferTypeCompact_::SIZE);
  ghostSendFloatBuffer_.resize(BufferTypeCompact_::SIZE);

  std::size_t reservedRow = maxHalfLocalDims[1]*lattice_.capacity();
  std::size_t reservedCol = maxHalfLocalDims[0]*lattice_.capacity();
  std::size_t reservedCorner = ghostExtent_[0]*ghostExtent_[1]*lattice_.capacity();

  ghostSendIndsBuffer_[BufferTypeCompact_::ROW].reserve(reservedRow);  
  ghostSendIndsBuffer_[BufferTypeCompact_::COL].reserve(reservedCol);
  ghostSendIndsBuffer_[BufferTypeCompact_::CORNER].reserve(reservedCorner);

  ghostSendIntBuffer_[BufferTypeCompact_::ROW].reserve(nIntsPerCell_*reservedRow);
  ghostSendIntBuffer_[BufferTypeCompact_::COL].reserve(nIntsPerCell_*reservedCol);
  ghostSendIntBuffer_[BufferTypeCompact_::CORNER].reserve(nIntsPerCell_*reservedCorner);

  ghostSendFloatBuffer_[BufferTypeCompact_::ROW].reserve(nFloatsPerCell_*reservedRow);
  ghostSendFloatBuffer_[BufferTypeCompact_::COL].reserve(nFloatsPerCell_*reservedCol);
  ghostSendFloatBuffer_[BufferTypeCompact_::CORNER].reserve(nFloatsPerCell_*reservedCorner);

  localSendIndsBuffer_.resize(nSectors_);
  localSendIntBuffer_.resize(nSectors_);
  localSendFloatBuffer_.resize(nSectors_);

  for (int i = 0; i < nSectors_; ++i) {

    localSendIndsBuffer_[i].resize(BufferTypeCompact_::SIZE);
    localSendIndsBuffer_[i][BufferTypeCompact_::ROW].reserve(reservedRow);
    localSendIndsBuffer_[i][BufferTypeCompact_::COL].reserve(reservedCol);
    localSendIndsBuffer_[i][BufferTypeCompact_::CORNER].reserve(reservedCorner);

    localSendIntBuffer_[i].resize(BufferTypeCompact_::SIZE);
    localSendIntBuffer_[i][BufferTypeCompact_::ROW].reserve(nIntsPerCell_*reservedRow);
    localSendIntBuffer_[i][BufferTypeCompact_::COL].reserve(nIntsPerCell_*reservedCol);
    localSendIntBuffer_[i][BufferTypeCompact_::CORNER].reserve(nIntsPerCell_*reservedCorner);

    localSendFloatBuffer_[i].resize(BufferTypeCompact_::SIZE);
    localSendFloatBuffer_[i][BufferTypeCompact_::ROW].reserve(nFloatsPerCell_*reservedRow);
    localSendFloatBuffer_[i][BufferTypeCompact_::COL].reserve(nFloatsPerCell_*reservedCol);
    localSendFloatBuffer_[i][BufferTypeCompact_::CORNER].reserve(nFloatsPerCell_*reservedCorner);

  }

  ghostRecvIndsBuffer_.resize(BufferTypeCompact_::SIZE);
  localRecvIndsBuffer_.resize(BufferTypeCompact_::SIZE);
  ghostRecvIntBuffer_.resize(BufferTypeCompact_::SIZE);
  localRecvIntBuffer_.resize(BufferTypeCompact_::SIZE);
  ghostRecvFloatBuffer_.resize(BufferTypeCompact_::SIZE);
  localRecvFloatBuffer_.resize(BufferTypeCompact_::SIZE);

  sendRecvGhostCompactInfo_.resize(nSectors_);

  boost::array<LatticePlanarBBox, ExportVal2_::SIZE> ghostExportBBox, localExportBBox;

  // ghostExportBBox[ExportVal2_::NONE] is just dummy values; I just
  // need to make sure that imin = imaxP1 so that the following "for"
  // loop works.
  ghostExportBBox[ExportVal2_::NONE].imin = ghostExportBBox[ExportVal2_::NONE].imaxP1 = 0;
  
  ghostExportBBox[ExportVal2_::BackRow1].imin = globalOffset_[0] - ghostExtent_[0];
  ghostExportBBox[ExportVal2_::BackRow1].jmin = globalOffset_[1];
  ghostExportBBox[ExportVal2_::BackRow1].imaxP1 = ghostExportBBox[ExportVal2_::BackRow1].imin + ghostExtent_[0];
  // Note that the third addend is -ghostExtent_[1] not
  // +ghostExtent_[1], to take into account the space taken by
  // ghostExportBBox[ExportVal2_::BackRowOverlap].
  ghostExportBBox[ExportVal2_::BackRow1].jmaxP1 = ghostExportBBox[ExportVal2_::BackRow1].jmin + firstHalfLocalDims[1] - ghostExtent_[1];

  ghostExportBBox[ExportVal2_::BackRow2].imin = globalOffset_[0] - ghostExtent_[0];
  // Note that the third addend is +ghostExtent_[1] not
  // -ghostExtent_[1], to take into account the space taken by
  // ghostExportBBox[ExportVal2_::BackRowOverlap].
  ghostExportBBox[ExportVal2_::BackRow2].jmin = globalOffset_[1] + firstHalfLocalDims[1] + ghostExtent_[1];
  ghostExportBBox[ExportVal2_::BackRow2].imaxP1 = ghostExportBBox[ExportVal2_::BackRow2].imin + ghostExtent_[0];
  ghostExportBBox[ExportVal2_::BackRow2].jmaxP1 = ghostExportBBox[ExportVal2_::BackRow2].jmin + secondHalfLocalDims[1] - ghostExtent_[1];

  ghostExportBBox[ExportVal2_::BackRowOverlap].imin = globalOffset_[0] - ghostExtent_[0];
  ghostExportBBox[ExportVal2_::BackRowOverlap].jmin = globalOffset_[1] + firstHalfLocalDims[1] - ghostExtent_[1];
  ghostExportBBox[ExportVal2_::BackRowOverlap].imaxP1 = ghostExportBBox[ExportVal2_::BackRowOverlap].imin + ghostExtent_[0];
  ghostExportBBox[ExportVal2_::BackRowOverlap].jmaxP1 = ghostExportBBox[ExportVal2_::BackRowOverlap].jmin + 2*ghostExtent_[1];

  ghostExportBBox[ExportVal2_::FwdRow1] = ghostExportBBox[ExportVal2_::BackRow1];
  ghostExportBBox[ExportVal2_::FwdRow1].imin += localDims_[0] + ghostExtent_[0];
  ghostExportBBox[ExportVal2_::FwdRow1].imaxP1 += localDims_[0] + ghostExtent_[0];
  
  ghostExportBBox[ExportVal2_::FwdRow2] = ghostExportBBox[ExportVal2_::BackRow2];
  ghostExportBBox[ExportVal2_::FwdRow2].imin += localDims_[0] + ghostExtent_[0];
  ghostExportBBox[ExportVal2_::FwdRow2].imaxP1 += localDims_[0] + ghostExtent_[0];

  ghostExportBBox[ExportVal2_::FwdRowOverlap] = ghostExportBBox[ExportVal2_::BackRowOverlap];
  ghostExportBBox[ExportVal2_::FwdRowOverlap].imin += localDims_[0] + ghostExtent_[0];
  ghostExportBBox[ExportVal2_::FwdRowOverlap].imaxP1 += localDims_[0] + ghostExtent_[0];

  ghostExportBBox[ExportVal2_::BackCol1].imin = globalOffset_[0];
  ghostExportBBox[ExportVal2_::BackCol1].jmin = globalOffset_[1] - ghostExtent_[1];
  // Note that the third addend is -ghostExtent_[0] not
  // +ghostExtent_[0], to take into account the space taken by
  // ghostExportBBox[ExportVal2_::BackColOverlap].
  ghostExportBBox[ExportVal2_::BackCol1].imaxP1 = ghostExportBBox[ExportVal2_::BackCol1].imin + firstHalfLocalDims[0] - ghostExtent_[0];
  ghostExportBBox[ExportVal2_::BackCol1].jmaxP1 = ghostExportBBox[ExportVal2_::BackCol1].jmin + ghostExtent_[1];

  // Note that the third addend is +ghostExtent_[0] not
  // -ghostExtent_[0], to take into account the space taken by
  // ghostExportBBox[ExportVal2_::BackRowOverlap].
  ghostExportBBox[ExportVal2_::BackCol2].imin = globalOffset_[0] + firstHalfLocalDims[0] + ghostExtent_[0];
  ghostExportBBox[ExportVal2_::BackCol2].jmin = globalOffset_[1] - ghostExtent_[0];
  ghostExportBBox[ExportVal2_::BackCol2].imaxP1 = ghostExportBBox[ExportVal2_::BackCol2].imin + secondHalfLocalDims[0] - ghostExtent_[0];
  ghostExportBBox[ExportVal2_::BackCol2].jmaxP1 = ghostExportBBox[ExportVal2_::BackCol2].jmin + ghostExtent_[1];

  ghostExportBBox[ExportVal2_::BackColOverlap].imin = globalOffset_[0] + firstHalfLocalDims[0] - ghostExtent_[0];
  ghostExportBBox[ExportVal2_::BackColOverlap].jmin = globalOffset_[1] - ghostExtent_[0];
  ghostExportBBox[ExportVal2_::BackColOverlap].imaxP1 = ghostExportBBox[ExportVal2_::BackColOverlap].imin + 2*ghostExtent_[0];
  ghostExportBBox[ExportVal2_::BackColOverlap].jmaxP1 = ghostExportBBox[ExportVal2_::BackColOverlap].jmin + ghostExtent_[1];

  ghostExportBBox[ExportVal2_::FwdCol1] = ghostExportBBox[ExportVal2_::BackCol1];
  ghostExportBBox[ExportVal2_::FwdCol1].jmin += localDims_[1] + ghostExtent_[1];
  ghostExportBBox[ExportVal2_::FwdCol1].jmaxP1 += localDims_[1] + ghostExtent_[1];

  ghostExportBBox[ExportVal2_::FwdCol2] = ghostExportBBox[ExportVal2_::BackCol2];
  ghostExportBBox[ExportVal2_::FwdCol2].jmin += localDims_[1] + ghostExtent_[1];
  ghostExportBBox[ExportVal2_::FwdCol2].jmaxP1 += localDims_[1] + ghostExtent_[1];

  ghostExportBBox[ExportVal2_::FwdColOverlap] = ghostExportBBox[ExportVal2_::BackColOverlap];
  ghostExportBBox[ExportVal2_::FwdColOverlap].jmin += localDims_[1] + ghostExtent_[1];
  ghostExportBBox[ExportVal2_::FwdColOverlap].jmaxP1 += localDims_[1] + ghostExtent_[1];

  ghostExportBBox[ExportVal2_::BackRowBackCol].imin = globalOffset_[0] - ghostExtent_[0];
  ghostExportBBox[ExportVal2_::BackRowBackCol].jmin = globalOffset_[1] - ghostExtent_[1];
  ghostExportBBox[ExportVal2_::BackRowBackCol].imaxP1 = ghostExportBBox[ExportVal2_::BackRowBackCol].imin + ghostExtent_[0];
  ghostExportBBox[ExportVal2_::BackRowBackCol].jmaxP1 = ghostExportBBox[ExportVal2_::BackRowBackCol].jmin + ghostExtent_[1];

  ghostExportBBox[ExportVal2_::FwdRowFwdCol].imin = globalOffset_[0] + localDims_[0];
  ghostExportBBox[ExportVal2_::FwdRowFwdCol].jmin = globalOffset_[1] + localDims_[1];
  ghostExportBBox[ExportVal2_::FwdRowFwdCol].imaxP1 = ghostExportBBox[ExportVal2_::FwdRowFwdCol].imin + ghostExtent_[0];
  ghostExportBBox[ExportVal2_::FwdRowFwdCol].jmaxP1 = ghostExportBBox[ExportVal2_::FwdRowFwdCol].jmin + ghostExtent_[1];

  ghostExportBBox[ExportVal2_::BackRowFwdCol].imin = ghostExportBBox[ExportVal2_::BackRowBackCol].imin;
  ghostExportBBox[ExportVal2_::BackRowFwdCol].jmin = ghostExportBBox[ExportVal2_::FwdRowFwdCol].jmin;
  ghostExportBBox[ExportVal2_::BackRowFwdCol].imaxP1 = ghostExportBBox[ExportVal2_::BackRowFwdCol].imin + ghostExtent_[0];
  ghostExportBBox[ExportVal2_::BackRowFwdCol].jmaxP1 = ghostExportBBox[ExportVal2_::BackRowFwdCol].jmin + ghostExtent_[1];

  ghostExportBBox[ExportVal2_::FwdRowBackCol].imin = ghostExportBBox[ExportVal2_::FwdRowFwdCol].imin;
  ghostExportBBox[ExportVal2_::FwdRowBackCol].jmin = ghostExportBBox[ExportVal2_::BackRowBackCol].jmin;
  ghostExportBBox[ExportVal2_::FwdRowBackCol].imaxP1 = ghostExportBBox[ExportVal2_::FwdRowBackCol].imin + ghostExtent_[0];
  ghostExportBBox[ExportVal2_::FwdRowBackCol].jmaxP1 = ghostExportBBox[ExportVal2_::FwdRowBackCol].jmin + ghostExtent_[1];

  localExportBBox[ExportVal2_::NONE].imin = globalOffset_[0] + ghostExtent_[0];
  localExportBBox[ExportVal2_::NONE].jmin = globalOffset_[1] + ghostExtent_[1];
  localExportBBox[ExportVal2_::NONE].imaxP1 = globalOffset_[0] + localDims_[0] - ghostExtent_[0];
  localExportBBox[ExportVal2_::NONE].jmaxP1 = globalOffset_[1] + localDims_[1] - ghostExtent_[1];

  localExportBBox[ExportVal2_::BackRow1].imin = globalOffset_[0];
  localExportBBox[ExportVal2_::BackRow1].jmin = globalOffset_[1] + ghostExtent_[1];
  localExportBBox[ExportVal2_::BackRow1].imaxP1 = localExportBBox[ExportVal2_::BackRow1].imin + ghostExtent_[0];
  localExportBBox[ExportVal2_::BackRow1].jmaxP1 = localExportBBox[ExportVal2_::BackRow1].jmin + firstHalfLocalDims[1] - 2*ghostExtent_[1];

  localExportBBox[ExportVal2_::BackRow2].imin = globalOffset_[0];
  localExportBBox[ExportVal2_::BackRow2].jmin = globalOffset_[1] + firstHalfLocalDims[1] + ghostExtent_[1];
  localExportBBox[ExportVal2_::BackRow2].imaxP1 = localExportBBox[ExportVal2_::BackRow2].imin + ghostExtent_[0];
  localExportBBox[ExportVal2_::BackRow2].jmaxP1 = localExportBBox[ExportVal2_::BackRow2].jmin + secondHalfLocalDims[1] - 2*ghostExtent_[1];
  
  localExportBBox[ExportVal2_::BackRowOverlap] = ghostExportBBox[ExportVal2_::BackRowOverlap];
  localExportBBox[ExportVal2_::BackRowOverlap].imin += ghostExtent_[0];
  localExportBBox[ExportVal2_::BackRowOverlap].imaxP1 += ghostExtent_[0];

  localExportBBox[ExportVal2_::FwdRow1] = localExportBBox[ExportVal2_::BackRow1];
  localExportBBox[ExportVal2_::FwdRow1].imin += localDims_[0] - ghostExtent_[0];
  localExportBBox[ExportVal2_::FwdRow1].imaxP1 += localDims_[0] - ghostExtent_[0];
  
  localExportBBox[ExportVal2_::FwdRow2] = localExportBBox[ExportVal2_::BackRow2];
  localExportBBox[ExportVal2_::FwdRow2].imin += localDims_[0] - ghostExtent_[0];
  localExportBBox[ExportVal2_::FwdRow2].imaxP1 += localDims_[0] - ghostExtent_[0];

  localExportBBox[ExportVal2_::FwdRowOverlap] = localExportBBox[ExportVal2_::BackRowOverlap];
  localExportBBox[ExportVal2_::FwdRowOverlap].imin += localDims_[0] - ghostExtent_[0];
  localExportBBox[ExportVal2_::FwdRowOverlap].imaxP1 += localDims_[0] - ghostExtent_[0];

  localExportBBox[ExportVal2_::BackCol1].imin = globalOffset_[0] + ghostExtent_[0];
  localExportBBox[ExportVal2_::BackCol1].jmin = globalOffset_[1];
  localExportBBox[ExportVal2_::BackCol1].imaxP1 = localExportBBox[ExportVal2_::BackCol1].imin + firstHalfLocalDims[0] - 2*ghostExtent_[0];
  localExportBBox[ExportVal2_::BackCol1].jmaxP1 = localExportBBox[ExportVal2_::BackCol1].jmin + ghostExtent_[1];

  localExportBBox[ExportVal2_::BackCol2].imin = globalOffset_[0] + firstHalfLocalDims[0] + ghostExtent_[0];
  localExportBBox[ExportVal2_::BackCol2].jmin = globalOffset_[1];
  localExportBBox[ExportVal2_::BackCol2].imaxP1 = localExportBBox[ExportVal2_::BackCol2].imin + secondHalfLocalDims[0] - 2*ghostExtent_[0];
  localExportBBox[ExportVal2_::BackCol2].jmaxP1 = localExportBBox[ExportVal2_::BackCol2].jmin + ghostExtent_[1];

  localExportBBox[ExportVal2_::BackColOverlap].imin = globalOffset_[0] + firstHalfLocalDims[0] - ghostExtent_[0];
  localExportBBox[ExportVal2_::BackColOverlap].jmin = globalOffset_[1];
  localExportBBox[ExportVal2_::BackColOverlap].imaxP1 = localExportBBox[ExportVal2_::BackColOverlap].imin + 2*ghostExtent_[0];
  localExportBBox[ExportVal2_::BackColOverlap].jmaxP1 = localExportBBox[ExportVal2_::BackColOverlap].jmin + ghostExtent_[1];

  localExportBBox[ExportVal2_::FwdCol1] = localExportBBox[ExportVal2_::BackCol1];
  localExportBBox[ExportVal2_::FwdCol1].jmin += localDims_[1] - ghostExtent_[1];
  localExportBBox[ExportVal2_::FwdCol1].jmaxP1 += localDims_[1] - ghostExtent_[1];

  localExportBBox[ExportVal2_::FwdCol2] = localExportBBox[ExportVal2_::BackCol2];
  localExportBBox[ExportVal2_::FwdCol2].jmin += localDims_[1] - ghostExtent_[1];
  localExportBBox[ExportVal2_::FwdCol2].jmaxP1 += localDims_[1] - ghostExtent_[1];

  localExportBBox[ExportVal2_::FwdColOverlap] = localExportBBox[ExportVal2_::BackColOverlap];
  localExportBBox[ExportVal2_::FwdColOverlap].jmin += localDims_[1] - ghostExtent_[1];
  localExportBBox[ExportVal2_::FwdColOverlap].jmaxP1 += localDims_[1] - ghostExtent_[1];

  localExportBBox[ExportVal2_::BackRowBackCol].imin = globalOffset_[0];
  localExportBBox[ExportVal2_::BackRowBackCol].jmin = globalOffset_[1];
  localExportBBox[ExportVal2_::BackRowBackCol].imaxP1 = localExportBBox[ExportVal2_::BackRowBackCol].imin + ghostExtent_[0];
  localExportBBox[ExportVal2_::BackRowBackCol].jmaxP1 = localExportBBox[ExportVal2_::BackRowBackCol].jmin + ghostExtent_[1];

  localExportBBox[ExportVal2_::FwdRowFwdCol].imin = globalOffset_[0] + localDims_[0] - ghostExtent_[0];
  localExportBBox[ExportVal2_::FwdRowFwdCol].jmin = globalOffset_[1] + localDims_[1] - ghostExtent_[1];
  localExportBBox[ExportVal2_::FwdRowFwdCol].imaxP1 = localExportBBox[ExportVal2_::FwdRowFwdCol].imin + ghostExtent_[0];
  localExportBBox[ExportVal2_::FwdRowFwdCol].jmaxP1 = localExportBBox[ExportVal2_::FwdRowFwdCol].jmin + ghostExtent_[1];

  localExportBBox[ExportVal2_::BackRowFwdCol].imin = localExportBBox[ExportVal2_::BackRowBackCol].imin;
  localExportBBox[ExportVal2_::BackRowFwdCol].jmin = localExportBBox[ExportVal2_::FwdRowFwdCol].jmin;
  localExportBBox[ExportVal2_::BackRowFwdCol].imaxP1 = localExportBBox[ExportVal2_::BackRowFwdCol].imin + ghostExtent_[0];
  localExportBBox[ExportVal2_::BackRowFwdCol].jmaxP1 = localExportBBox[ExportVal2_::BackRowFwdCol].jmin + ghostExtent_[1];

  localExportBBox[ExportVal2_::FwdRowBackCol].imin = localExportBBox[ExportVal2_::FwdRowFwdCol].imin;
  localExportBBox[ExportVal2_::FwdRowBackCol].jmin = localExportBBox[ExportVal2_::BackRowBackCol].jmin;
  localExportBBox[ExportVal2_::FwdRowBackCol].imaxP1 = localExportBBox[ExportVal2_::FwdRowBackCol].imin + ghostExtent_[0];
  localExportBBox[ExportVal2_::FwdRowBackCol].jmaxP1 = localExportBBox[ExportVal2_::FwdRowBackCol].jmin + ghostExtent_[1];
  
  for (int exportVal = 0; exportVal < ExportVal2_::SIZE; ++exportVal) {
    
    for (int i = ghostExportBBox[exportVal].imin; i < ghostExportBBox[exportVal].imaxP1; ++i) {
      for (int j = ghostExportBBox[exportVal].jmin; j < ghostExportBBox[exportVal].jmaxP1; ++j) {
        cellParInfoArray_[i][j].sectNum = -1;
        cellParInfoArray_[i][j].exportVal = exportVal;
      }
    }

    for (int i = localExportBBox[exportVal].imin; i < localExportBBox[exportVal].imaxP1; ++i) {
      for (int j = localExportBBox[exportVal].jmin; j < localExportBBox[exportVal].jmaxP1; ++j) {
        cellParInfoArray_[i][j].exportVal = exportVal;
      }
    }

  }

  // Sector 0
  // ========

  // Corner
  sendRecvGhostCompactInfo_[0].sendCornerRank = pIdNeighs[ProcIdCart2_::FwdRowFwdCol];
  sendRecvGhostCompactInfo_[0].recvCornerRank = pIdNeighs[ProcIdCart2_::BackRowBackCol];

  sendRecvGhostCompactInfo_[0].sendCornerPlanarCoords[0] = localDims_[0] - ghostExtent_[0] + globalOffset_[0];
  sendRecvGhostCompactInfo_[0].recvCornerPlanarCoords[0] =                -ghostExtent_[0] + globalOffset_[0];

  sendRecvGhostCompactInfo_[0].sendCornerPlanarCoords[1] = localDims_[1] - ghostExtent_[1] + globalOffset_[1];
  sendRecvGhostCompactInfo_[0].recvCornerPlanarCoords[1] =                -ghostExtent_[1] + globalOffset_[1];

  localRecvBufferBounds_[0][BufferTypeCompact_::CORNER].imin = sendRecvGhostCompactInfo_[0].sendCornerPlanarCoords[0];
  localRecvBufferBounds_[0][BufferTypeCompact_::CORNER].jmin = sendRecvGhostCompactInfo_[0].sendCornerPlanarCoords[1];

  localRecvBufferBounds_[0][BufferTypeCompact_::CORNER].imaxP1 = localRecvBufferBounds_[0][BufferTypeCompact_::CORNER].imin + ghostExtent_[0];
  localRecvBufferBounds_[0][BufferTypeCompact_::CORNER].jmaxP1 = localRecvBufferBounds_[0][BufferTypeCompact_::CORNER].jmin + ghostExtent_[1];

  ghostRecvBufferBounds_[0][BufferTypeCompact_::CORNER].imin = sendRecvGhostCompactInfo_[0].recvCornerPlanarCoords[0];
  ghostRecvBufferBounds_[0][BufferTypeCompact_::CORNER].jmin = sendRecvGhostCompactInfo_[0].recvCornerPlanarCoords[1];

  ghostRecvBufferBounds_[0][BufferTypeCompact_::CORNER].imaxP1 = ghostRecvBufferBounds_[0][BufferTypeCompact_::CORNER].imin + ghostExtent_[0];
  ghostRecvBufferBounds_[0][BufferTypeCompact_::CORNER].jmaxP1 = ghostRecvBufferBounds_[0][BufferTypeCompact_::CORNER].jmin + ghostExtent_[1];  

  // Partial row
  sendRecvGhostCompactInfo_[0].whichHalfRow = 0;

  sendRecvGhostCompactInfo_[0].sendHalfRowRank = pIdNeighs[ProcIdCart2_::FwdRow];
  sendRecvGhostCompactInfo_[0].recvHalfRowRank = pIdNeighs[ProcIdCart2_::BackRow];

  sendRecvGhostCompactInfo_[0].sendHalfRowPlanarCoords[0] = localDims_[0] - ghostExtent_[0] + globalOffset_[0];
  sendRecvGhostCompactInfo_[0].recvHalfRowPlanarCoords[0] =                -ghostExtent_[0] + globalOffset_[0];

  sendRecvGhostCompactInfo_[0].sendHalfRowPlanarCoords[1] = sendRecvGhostCompactInfo_[0].recvHalfRowPlanarCoords[1] =  globalOffset_[1];

  localRecvBufferBounds_[0][BufferTypeCompact_::ROW].imin = sendRecvGhostCompactInfo_[0].sendHalfRowPlanarCoords[0];
  localRecvBufferBounds_[0][BufferTypeCompact_::ROW].jmin = sendRecvGhostCompactInfo_[0].sendHalfRowPlanarCoords[1];

  localRecvBufferBounds_[0][BufferTypeCompact_::ROW].imaxP1 = localRecvBufferBounds_[0][BufferTypeCompact_::ROW].imin + ghostExtent_[0];
  localRecvBufferBounds_[0][BufferTypeCompact_::ROW].jmaxP1 = localRecvBufferBounds_[0][BufferTypeCompact_::ROW].jmin + firstHalfLocalDims[1] + ghostExtent_[1];

  ghostRecvBufferBounds_[0][BufferTypeCompact_::ROW].imin = sendRecvGhostCompactInfo_[0].recvHalfRowPlanarCoords[0];
  ghostRecvBufferBounds_[0][BufferTypeCompact_::ROW].jmin = sendRecvGhostCompactInfo_[0].recvHalfRowPlanarCoords[1];

  ghostRecvBufferBounds_[0][BufferTypeCompact_::ROW].imaxP1 = ghostRecvBufferBounds_[0][BufferTypeCompact_::ROW].imin + ghostExtent_[0];
  ghostRecvBufferBounds_[0][BufferTypeCompact_::ROW].jmaxP1 = ghostRecvBufferBounds_[0][BufferTypeCompact_::ROW].jmin + firstHalfLocalDims[1] + ghostExtent_[1];
  
  // Partial column
  sendRecvGhostCompactInfo_[0].whichHalfCol = 0;

  sendRecvGhostCompactInfo_[0].sendHalfColRank = pIdNeighs[ProcIdCart2_::FwdCol];
  sendRecvGhostCompactInfo_[0].recvHalfColRank = pIdNeighs[ProcIdCart2_::BackCol];
  
  sendRecvGhostCompactInfo_[0].sendHalfColPlanarCoords[0] = sendRecvGhostCompactInfo_[0].recvHalfColPlanarCoords[0] = globalOffset_[0];

  sendRecvGhostCompactInfo_[0].sendHalfColPlanarCoords[1] = localDims_[1] - ghostExtent_[1] + globalOffset_[1];
  sendRecvGhostCompactInfo_[0].recvHalfColPlanarCoords[1] =                -ghostExtent_[1] + globalOffset_[1];

  localRecvBufferBounds_[0][BufferTypeCompact_::COL].imin = sendRecvGhostCompactInfo_[0].sendHalfColPlanarCoords[0];
  localRecvBufferBounds_[0][BufferTypeCompact_::COL].jmin = sendRecvGhostCompactInfo_[0].sendHalfColPlanarCoords[1];

  localRecvBufferBounds_[0][BufferTypeCompact_::COL].imaxP1 = localRecvBufferBounds_[0][BufferTypeCompact_::COL].imin + firstHalfLocalDims[0] + ghostExtent_[0];
  localRecvBufferBounds_[0][BufferTypeCompact_::COL].jmaxP1 = localRecvBufferBounds_[0][BufferTypeCompact_::COL].jmin + ghostExtent_[1];

  ghostRecvBufferBounds_[0][BufferTypeCompact_::COL].imin = sendRecvGhostCompactInfo_[0].recvHalfColPlanarCoords[0];
  ghostRecvBufferBounds_[0][BufferTypeCompact_::COL].jmin = sendRecvGhostCompactInfo_[0].recvHalfColPlanarCoords[1];

  ghostRecvBufferBounds_[0][BufferTypeCompact_::COL].imaxP1 = ghostRecvBufferBounds_[0][BufferTypeCompact_::COL].imin + firstHalfLocalDims[0] + ghostExtent_[0];
  ghostRecvBufferBounds_[0][BufferTypeCompact_::COL].jmaxP1 = ghostRecvBufferBounds_[0][BufferTypeCompact_::COL].jmin + ghostExtent_[1];

  // Sector boundaries
  iminS_[0] = globalOffset_[0];
  jminS_[0] = globalOffset_[1];

  imaxP1S_[0] = firstHalfLocalDims[0] + globalOffset_[0];
  jmaxP1S_[0] = firstHalfLocalDims[1] + globalOffset_[1];

  // Sector 1
  // ========

  // Corner
  sendRecvGhostCompactInfo_[1].sendCornerRank = pIdNeighs[ProcIdCart2_::BackRowFwdCol];
  sendRecvGhostCompactInfo_[1].recvCornerRank = pIdNeighs[ProcIdCart2_::FwdRowBackCol];

  sendRecvGhostCompactInfo_[1].sendCornerPlanarCoords[0] = globalOffset_[0];
  sendRecvGhostCompactInfo_[1].recvCornerPlanarCoords[0] = localDims_[0] + globalOffset_[0];

  sendRecvGhostCompactInfo_[1].sendCornerPlanarCoords[1] = localDims_[1] - ghostExtent_[1] + globalOffset_[1];
  sendRecvGhostCompactInfo_[1].recvCornerPlanarCoords[1] =                -ghostExtent_[1] + globalOffset_[1];

  localRecvBufferBounds_[1][BufferTypeCompact_::CORNER].imin = sendRecvGhostCompactInfo_[1].sendCornerPlanarCoords[0];
  localRecvBufferBounds_[1][BufferTypeCompact_::CORNER].jmin = sendRecvGhostCompactInfo_[1].sendCornerPlanarCoords[1];

  localRecvBufferBounds_[1][BufferTypeCompact_::CORNER].imaxP1 = localRecvBufferBounds_[1][BufferTypeCompact_::CORNER].imin + ghostExtent_[0];
  localRecvBufferBounds_[1][BufferTypeCompact_::CORNER].jmaxP1 = localRecvBufferBounds_[1][BufferTypeCompact_::CORNER].jmin + ghostExtent_[1];

  ghostRecvBufferBounds_[1][BufferTypeCompact_::CORNER].imin = sendRecvGhostCompactInfo_[1].recvCornerPlanarCoords[0];
  ghostRecvBufferBounds_[1][BufferTypeCompact_::CORNER].jmin = sendRecvGhostCompactInfo_[1].recvCornerPlanarCoords[1];

  ghostRecvBufferBounds_[1][BufferTypeCompact_::CORNER].imaxP1 = ghostRecvBufferBounds_[1][BufferTypeCompact_::CORNER].imin + ghostExtent_[0];
  ghostRecvBufferBounds_[1][BufferTypeCompact_::CORNER].jmaxP1 = ghostRecvBufferBounds_[1][BufferTypeCompact_::CORNER].jmin + ghostExtent_[1];  

  // Partial row
  sendRecvGhostCompactInfo_[1].whichHalfRow = 0;

  sendRecvGhostCompactInfo_[1].sendHalfRowRank = pIdNeighs[ProcIdCart2_::BackRow];
  sendRecvGhostCompactInfo_[1].recvHalfRowRank = pIdNeighs[ProcIdCart2_::FwdRow];

  sendRecvGhostCompactInfo_[1].sendHalfRowPlanarCoords[0] = globalOffset_[0];
  sendRecvGhostCompactInfo_[1].recvHalfRowPlanarCoords[0] = localDims_[0] + globalOffset_[0];

  sendRecvGhostCompactInfo_[1].sendHalfRowPlanarCoords[1] = sendRecvGhostCompactInfo_[1].recvHalfRowPlanarCoords[1] = globalOffset_[1];

  localRecvBufferBounds_[1][BufferTypeCompact_::ROW].imin = sendRecvGhostCompactInfo_[1].sendHalfRowPlanarCoords[0];
  localRecvBufferBounds_[1][BufferTypeCompact_::ROW].jmin = sendRecvGhostCompactInfo_[1].sendHalfRowPlanarCoords[1];

  localRecvBufferBounds_[1][BufferTypeCompact_::ROW].imaxP1 = localRecvBufferBounds_[1][BufferTypeCompact_::ROW].imin + ghostExtent_[0];
  localRecvBufferBounds_[1][BufferTypeCompact_::ROW].jmaxP1 = localRecvBufferBounds_[1][BufferTypeCompact_::ROW].jmin + firstHalfLocalDims[1] + ghostExtent_[1];

  ghostRecvBufferBounds_[1][BufferTypeCompact_::ROW].imin = sendRecvGhostCompactInfo_[1].recvHalfRowPlanarCoords[0];
  ghostRecvBufferBounds_[1][BufferTypeCompact_::ROW].jmin = sendRecvGhostCompactInfo_[1].recvHalfRowPlanarCoords[1];

  ghostRecvBufferBounds_[1][BufferTypeCompact_::ROW].imaxP1 = ghostRecvBufferBounds_[1][BufferTypeCompact_::ROW].imin + ghostExtent_[0];
  ghostRecvBufferBounds_[1][BufferTypeCompact_::ROW].jmaxP1 = ghostRecvBufferBounds_[1][BufferTypeCompact_::ROW].jmin + firstHalfLocalDims[1] + ghostExtent_[1];

  // Partial column
  sendRecvGhostCompactInfo_[1].whichHalfCol = 1;

  sendRecvGhostCompactInfo_[1].sendHalfColRank = pIdNeighs[ProcIdCart2_::FwdCol];
  sendRecvGhostCompactInfo_[1].recvHalfColRank = pIdNeighs[ProcIdCart2_::BackCol];

  sendRecvGhostCompactInfo_[1].sendHalfColPlanarCoords[0] = sendRecvGhostCompactInfo_[1].recvHalfColPlanarCoords[0] = firstHalfLocalDims[0] - ghostExtent_[0] + globalOffset_[0];

  sendRecvGhostCompactInfo_[1].sendHalfColPlanarCoords[1] = localDims_[1] - ghostExtent_[1] + globalOffset_[1];
  sendRecvGhostCompactInfo_[1].recvHalfColPlanarCoords[1] =                -ghostExtent_[1] + globalOffset_[1];

  localRecvBufferBounds_[1][BufferTypeCompact_::COL].imin = sendRecvGhostCompactInfo_[1].sendHalfColPlanarCoords[0];
  localRecvBufferBounds_[1][BufferTypeCompact_::COL].jmin = sendRecvGhostCompactInfo_[1].sendHalfColPlanarCoords[1];

  localRecvBufferBounds_[1][BufferTypeCompact_::COL].imaxP1 = localRecvBufferBounds_[1][BufferTypeCompact_::COL].imin + secondHalfLocalDims[0] + ghostExtent_[0];
  localRecvBufferBounds_[1][BufferTypeCompact_::COL].jmaxP1 = localRecvBufferBounds_[1][BufferTypeCompact_::COL].jmin + ghostExtent_[1];

  ghostRecvBufferBounds_[1][BufferTypeCompact_::COL].imin = sendRecvGhostCompactInfo_[1].recvHalfColPlanarCoords[0];
  ghostRecvBufferBounds_[1][BufferTypeCompact_::COL].jmin = sendRecvGhostCompactInfo_[1].recvHalfColPlanarCoords[1];

  ghostRecvBufferBounds_[1][BufferTypeCompact_::COL].imaxP1 = ghostRecvBufferBounds_[1][BufferTypeCompact_::COL].imin + secondHalfLocalDims[0] + ghostExtent_[0];
  ghostRecvBufferBounds_[1][BufferTypeCompact_::COL].jmaxP1 = ghostRecvBufferBounds_[1][BufferTypeCompact_::COL].jmin + ghostExtent_[1];

  // Sector boundaries
  iminS_[1] = firstHalfLocalDims[0] + globalOffset_[0];
  jminS_[1] = globalOffset_[1];

  imaxP1S_[1] = localDims_[0] + globalOffset_[0];
  jmaxP1S_[1] = firstHalfLocalDims[1] + globalOffset_[1];

  // Sector 2
  // ========

  // Corner
  sendRecvGhostCompactInfo_[2].sendCornerRank = pIdNeighs[ProcIdCart2_::BackRowBackCol];
  sendRecvGhostCompactInfo_[2].recvCornerRank = pIdNeighs[ProcIdCart2_::FwdRowFwdCol];

  sendRecvGhostCompactInfo_[2].sendCornerPlanarCoords[0] = globalOffset_[0];
  sendRecvGhostCompactInfo_[2].recvCornerPlanarCoords[0] = localDims_[0] + globalOffset_[0];

  sendRecvGhostCompactInfo_[2].sendCornerPlanarCoords[1] = globalOffset_[1];
  sendRecvGhostCompactInfo_[2].recvCornerPlanarCoords[1] = localDims_[1] + globalOffset_[1];

  localRecvBufferBounds_[2][BufferTypeCompact_::CORNER].imin = sendRecvGhostCompactInfo_[2].sendCornerPlanarCoords[0];
  localRecvBufferBounds_[2][BufferTypeCompact_::CORNER].jmin = sendRecvGhostCompactInfo_[2].sendCornerPlanarCoords[1];

  localRecvBufferBounds_[2][BufferTypeCompact_::CORNER].imaxP1 = localRecvBufferBounds_[2][BufferTypeCompact_::CORNER].imin + ghostExtent_[0];
  localRecvBufferBounds_[2][BufferTypeCompact_::CORNER].jmaxP1 = localRecvBufferBounds_[2][BufferTypeCompact_::CORNER].jmin + ghostExtent_[1];

  ghostRecvBufferBounds_[2][BufferTypeCompact_::CORNER].imin = sendRecvGhostCompactInfo_[2].recvCornerPlanarCoords[0];
  ghostRecvBufferBounds_[2][BufferTypeCompact_::CORNER].jmin = sendRecvGhostCompactInfo_[2].recvCornerPlanarCoords[1];

  ghostRecvBufferBounds_[2][BufferTypeCompact_::CORNER].imaxP1 = ghostRecvBufferBounds_[2][BufferTypeCompact_::CORNER].imin + ghostExtent_[0];
  ghostRecvBufferBounds_[2][BufferTypeCompact_::CORNER].jmaxP1 = ghostRecvBufferBounds_[2][BufferTypeCompact_::CORNER].jmin + ghostExtent_[1];

  // Partial row
  sendRecvGhostCompactInfo_[2].whichHalfRow = 1;

  sendRecvGhostCompactInfo_[2].sendHalfRowRank = pIdNeighs[ProcIdCart2_::BackRow];
  sendRecvGhostCompactInfo_[2].recvHalfRowRank = pIdNeighs[ProcIdCart2_::FwdRow];

  sendRecvGhostCompactInfo_[2].sendHalfRowPlanarCoords[0] = globalOffset_[0];
  sendRecvGhostCompactInfo_[2].recvHalfRowPlanarCoords[0] = localDims_[0] + globalOffset_[0];

  sendRecvGhostCompactInfo_[2].sendHalfRowPlanarCoords[1] = sendRecvGhostCompactInfo_[2].recvHalfRowPlanarCoords[1] = firstHalfLocalDims[1] - ghostExtent_[1] + globalOffset_[1];

  localRecvBufferBounds_[2][BufferTypeCompact_::ROW].imin = sendRecvGhostCompactInfo_[2].sendHalfRowPlanarCoords[0];
  localRecvBufferBounds_[2][BufferTypeCompact_::ROW].jmin = sendRecvGhostCompactInfo_[2].sendHalfRowPlanarCoords[1];

  localRecvBufferBounds_[2][BufferTypeCompact_::ROW].imaxP1 = localRecvBufferBounds_[2][BufferTypeCompact_::ROW].imin + ghostExtent_[0];
  localRecvBufferBounds_[2][BufferTypeCompact_::ROW].jmaxP1 = localRecvBufferBounds_[2][BufferTypeCompact_::ROW].jmin + secondHalfLocalDims[1] + ghostExtent_[1];

  ghostRecvBufferBounds_[2][BufferTypeCompact_::ROW].imin = sendRecvGhostCompactInfo_[2].recvHalfRowPlanarCoords[0];
  ghostRecvBufferBounds_[2][BufferTypeCompact_::ROW].jmin = sendRecvGhostCompactInfo_[2].recvHalfRowPlanarCoords[1];

  ghostRecvBufferBounds_[2][BufferTypeCompact_::ROW].imaxP1 = ghostRecvBufferBounds_[2][BufferTypeCompact_::ROW].imin + ghostExtent_[0];
  ghostRecvBufferBounds_[2][BufferTypeCompact_::ROW].jmaxP1 = ghostRecvBufferBounds_[2][BufferTypeCompact_::ROW].jmin + secondHalfLocalDims[1] + ghostExtent_[1];

  // Partial column
  sendRecvGhostCompactInfo_[2].whichHalfCol = 1;

  sendRecvGhostCompactInfo_[2].sendHalfColRank = pIdNeighs[ProcIdCart2_::BackCol];
  sendRecvGhostCompactInfo_[2].recvHalfColRank = pIdNeighs[ProcIdCart2_::FwdCol];

  sendRecvGhostCompactInfo_[2].sendHalfColPlanarCoords[0] = sendRecvGhostCompactInfo_[2].recvHalfColPlanarCoords[0] = firstHalfLocalDims[0] - ghostExtent_[0] + globalOffset_[0];

  sendRecvGhostCompactInfo_[2].sendHalfColPlanarCoords[1] = globalOffset_[1];
  sendRecvGhostCompactInfo_[2].recvHalfColPlanarCoords[1] = localDims_[1] + globalOffset_[1];

  localRecvBufferBounds_[2][BufferTypeCompact_::COL].imin = sendRecvGhostCompactInfo_[2].sendHalfColPlanarCoords[0];
  localRecvBufferBounds_[2][BufferTypeCompact_::COL].jmin = sendRecvGhostCompactInfo_[2].sendHalfColPlanarCoords[1];

  localRecvBufferBounds_[2][BufferTypeCompact_::COL].imaxP1 = localRecvBufferBounds_[2][BufferTypeCompact_::COL].imin + secondHalfLocalDims[0] + ghostExtent_[0];
  localRecvBufferBounds_[2][BufferTypeCompact_::COL].jmaxP1 = localRecvBufferBounds_[2][BufferTypeCompact_::COL].jmin + ghostExtent_[1];

  ghostRecvBufferBounds_[2][BufferTypeCompact_::COL].imin = sendRecvGhostCompactInfo_[2].recvHalfColPlanarCoords[0];
  ghostRecvBufferBounds_[2][BufferTypeCompact_::COL].jmin = sendRecvGhostCompactInfo_[2].recvHalfColPlanarCoords[1];

  ghostRecvBufferBounds_[2][BufferTypeCompact_::COL].imaxP1 = ghostRecvBufferBounds_[2][BufferTypeCompact_::COL].imin + secondHalfLocalDims[0] + ghostExtent_[0];
  ghostRecvBufferBounds_[2][BufferTypeCompact_::COL].jmaxP1 = ghostRecvBufferBounds_[2][BufferTypeCompact_::COL].jmin + ghostExtent_[1];

  // Sector boundaries
  iminS_[2] = firstHalfLocalDims[0] + globalOffset_[0];
  jminS_[2] = firstHalfLocalDims[1] + globalOffset_[1];

  imaxP1S_[2] = localDims_[0] + globalOffset_[0];
  jmaxP1S_[2] = localDims_[1] + globalOffset_[1];

  // Sector 3
  // ========

  // Corner
  sendRecvGhostCompactInfo_[3].sendCornerRank = pIdNeighs[ProcIdCart2_::FwdRowBackCol];
  sendRecvGhostCompactInfo_[3].recvCornerRank = pIdNeighs[ProcIdCart2_::BackRowFwdCol];

  sendRecvGhostCompactInfo_[3].sendCornerPlanarCoords[0] = localDims_[0] - ghostExtent_[0] + globalOffset_[0];
  sendRecvGhostCompactInfo_[3].recvCornerPlanarCoords[0] =                -ghostExtent_[0] + globalOffset_[0];

  sendRecvGhostCompactInfo_[3].sendCornerPlanarCoords[1] = globalOffset_[1];
  sendRecvGhostCompactInfo_[3].recvCornerPlanarCoords[1] = localDims_[1] + globalOffset_[1];

  localRecvBufferBounds_[3][BufferTypeCompact_::CORNER].imin = sendRecvGhostCompactInfo_[3].sendCornerPlanarCoords[0];
  localRecvBufferBounds_[3][BufferTypeCompact_::CORNER].jmin = sendRecvGhostCompactInfo_[3].sendCornerPlanarCoords[1];

  localRecvBufferBounds_[3][BufferTypeCompact_::CORNER].imaxP1 = localRecvBufferBounds_[3][BufferTypeCompact_::CORNER].imin + ghostExtent_[0];
  localRecvBufferBounds_[3][BufferTypeCompact_::CORNER].jmaxP1 = localRecvBufferBounds_[3][BufferTypeCompact_::CORNER].jmin + ghostExtent_[1];

  ghostRecvBufferBounds_[3][BufferTypeCompact_::CORNER].imin = sendRecvGhostCompactInfo_[3].recvCornerPlanarCoords[0];
  ghostRecvBufferBounds_[3][BufferTypeCompact_::CORNER].jmin = sendRecvGhostCompactInfo_[3].recvCornerPlanarCoords[1];

  ghostRecvBufferBounds_[3][BufferTypeCompact_::CORNER].imaxP1 = ghostRecvBufferBounds_[3][BufferTypeCompact_::CORNER].imin + ghostExtent_[0];
  ghostRecvBufferBounds_[3][BufferTypeCompact_::CORNER].jmaxP1 = ghostRecvBufferBounds_[3][BufferTypeCompact_::CORNER].jmin + ghostExtent_[1];  

  // Partial row
  sendRecvGhostCompactInfo_[3].whichHalfRow = 1;

  sendRecvGhostCompactInfo_[3].sendHalfRowRank = pIdNeighs[ProcIdCart2_::FwdRow];
  sendRecvGhostCompactInfo_[3].recvHalfRowRank = pIdNeighs[ProcIdCart2_::BackRow];

  sendRecvGhostCompactInfo_[3].sendHalfRowPlanarCoords[0] = localDims_[0] - ghostExtent_[0] + globalOffset_[0];
  sendRecvGhostCompactInfo_[3].recvHalfRowPlanarCoords[0] =                -ghostExtent_[0] + globalOffset_[0];

  sendRecvGhostCompactInfo_[3].sendHalfRowPlanarCoords[1] = sendRecvGhostCompactInfo_[3].recvHalfRowPlanarCoords[1] = firstHalfLocalDims[1] - ghostExtent_[1] + globalOffset_[1];

  localRecvBufferBounds_[3][BufferTypeCompact_::ROW].imin = sendRecvGhostCompactInfo_[3].sendHalfRowPlanarCoords[0];
  localRecvBufferBounds_[3][BufferTypeCompact_::ROW].jmin = sendRecvGhostCompactInfo_[3].sendHalfRowPlanarCoords[1];

  localRecvBufferBounds_[3][BufferTypeCompact_::ROW].imaxP1 = localRecvBufferBounds_[3][BufferTypeCompact_::ROW].imin + ghostExtent_[0];
  localRecvBufferBounds_[3][BufferTypeCompact_::ROW].jmaxP1 = localRecvBufferBounds_[3][BufferTypeCompact_::ROW].jmin + secondHalfLocalDims[1] + ghostExtent_[1];

  ghostRecvBufferBounds_[3][BufferTypeCompact_::ROW].imin = sendRecvGhostCompactInfo_[3].recvHalfRowPlanarCoords[0];
  ghostRecvBufferBounds_[3][BufferTypeCompact_::ROW].jmin = sendRecvGhostCompactInfo_[3].recvHalfRowPlanarCoords[1];

  ghostRecvBufferBounds_[3][BufferTypeCompact_::ROW].imaxP1 = ghostRecvBufferBounds_[3][BufferTypeCompact_::ROW].imin + ghostExtent_[0];
  ghostRecvBufferBounds_[3][BufferTypeCompact_::ROW].jmaxP1 = ghostRecvBufferBounds_[3][BufferTypeCompact_::ROW].jmin + secondHalfLocalDims[1] + ghostExtent_[1];

  // Partial column
  sendRecvGhostCompactInfo_[3].whichHalfCol = 0;

  sendRecvGhostCompactInfo_[3].sendHalfColRank = pIdNeighs[ProcIdCart2_::BackCol];
  sendRecvGhostCompactInfo_[3].recvHalfColRank = pIdNeighs[ProcIdCart2_::FwdCol];

  sendRecvGhostCompactInfo_[3].sendHalfColPlanarCoords[0] = sendRecvGhostCompactInfo_[3].recvHalfColPlanarCoords[0] = globalOffset_[0];

  sendRecvGhostCompactInfo_[3].sendHalfColPlanarCoords[1] = globalOffset_[1];
  sendRecvGhostCompactInfo_[3].recvHalfColPlanarCoords[1] = localDims_[1] + globalOffset_[1];

  localRecvBufferBounds_[3][BufferTypeCompact_::COL].imin = sendRecvGhostCompactInfo_[3].sendHalfColPlanarCoords[0];
  localRecvBufferBounds_[3][BufferTypeCompact_::COL].jmin = sendRecvGhostCompactInfo_[3].sendHalfColPlanarCoords[1];

  localRecvBufferBounds_[3][BufferTypeCompact_::COL].imaxP1 = localRecvBufferBounds_[3][BufferTypeCompact_::COL].imin + firstHalfLocalDims[0] + ghostExtent_[0];
  localRecvBufferBounds_[3][BufferTypeCompact_::COL].jmaxP1 = localRecvBufferBounds_[3][BufferTypeCompact_::COL].jmin + ghostExtent_[1];

  ghostRecvBufferBounds_[3][BufferTypeCompact_::COL].imin = sendRecvGhostCompactInfo_[3].recvHalfColPlanarCoords[0];
  ghostRecvBufferBounds_[3][BufferTypeCompact_::COL].jmin = sendRecvGhostCompactInfo_[3].recvHalfColPlanarCoords[1];

  ghostRecvBufferBounds_[3][BufferTypeCompact_::COL].imaxP1 = ghostRecvBufferBounds_[3][BufferTypeCompact_::COL].imin + firstHalfLocalDims[0] + ghostExtent_[0];
  ghostRecvBufferBounds_[3][BufferTypeCompact_::COL].jmaxP1 = ghostRecvBufferBounds_[3][BufferTypeCompact_::COL].jmin + ghostExtent_[1];

  // Sector boundaries
  iminS_[3] = globalOffset_[0];
  jminS_[3] = firstHalfLocalDims[1] + globalOffset_[1];

  imaxP1S_[3] = firstHalfLocalDims[0] + globalOffset_[0];
  jmaxP1S_[3] = localDims_[1] + globalOffset_[1];


  for (int sectNum = 0; sectNum < nSectors_; ++sectNum) {
    for (int i = iminS_[sectNum]; i < imaxP1S_[sectNum]; ++i) {
      for (int j = jminS_[sectNum]; j < jmaxP1S_[sectNum]; ++j) {
        cellParInfoArray_[i][j].sectNum = sectNum;
      }
    }
  }

  exportValToGhostBufType_.resize(ExportVal2_::SIZE);
  exportValToLocalBufferLocs_.resize(ExportVal2_::SIZE);

  exportValToGhostBufType_[ExportVal2_::NONE] = -42; // Garbage value, should never get used.

  exportValToGhostBufType_[ExportVal2_::BackRow1] = BufferTypeCompact_::ROW;
  exportValToGhostBufType_[ExportVal2_::FwdRow1] = BufferTypeCompact_::ROW;
  exportValToGhostBufType_[ExportVal2_::BackRow2] = BufferTypeCompact_::ROW;
  exportValToGhostBufType_[ExportVal2_::FwdRow2] = BufferTypeCompact_::ROW;

  exportValToGhostBufType_[ExportVal2_::BackRowOverlap] = BufferTypeCompact_::ROW;
  exportValToGhostBufType_[ExportVal2_::FwdRowOverlap] = BufferTypeCompact_::ROW;

  exportValToGhostBufType_[ExportVal2_::BackCol1] = BufferTypeCompact_::COL;
  exportValToGhostBufType_[ExportVal2_::BackCol2] = BufferTypeCompact_::COL;
  exportValToGhostBufType_[ExportVal2_::FwdCol1] = BufferTypeCompact_::COL;
  exportValToGhostBufType_[ExportVal2_::FwdCol2] = BufferTypeCompact_::COL;

  exportValToGhostBufType_[ExportVal2_::BackColOverlap] = BufferTypeCompact_::COL;
  exportValToGhostBufType_[ExportVal2_::FwdColOverlap] = BufferTypeCompact_::COL;

  exportValToGhostBufType_[ExportVal2_::BackRowBackCol] = BufferTypeCompact_::CORNER;
  exportValToGhostBufType_[ExportVal2_::FwdRowFwdCol] = BufferTypeCompact_::CORNER;
  exportValToGhostBufType_[ExportVal2_::BackRowFwdCol] = BufferTypeCompact_::CORNER;
  exportValToGhostBufType_[ExportVal2_::FwdRowBackCol] = BufferTypeCompact_::CORNER;

  exportValToLocalBufferLocs_[ExportVal2_::NONE].nBuffers = 0;

  exportValToLocalBufferLocs_[ExportVal2_::BackRow1].nBuffers = 1;
  exportValToLocalBufferLocs_[ExportVal2_::BackRow1].sectNum[0] = 1; /* Note
                                                                        that this sector 
                                                                        is *not* the sector 
                                                                        that BackRow1 is located in, 
                                                                        but rather the sector responsible 
                                                                        for exporting BackRow1 off-proc.*/
  exportValToLocalBufferLocs_[ExportVal2_::BackRow1].bufType[0] = BufferTypeCompact_::ROW;

  exportValToLocalBufferLocs_[ExportVal2_::FwdRow1].nBuffers = 1;
  exportValToLocalBufferLocs_[ExportVal2_::FwdRow1].sectNum[0] = 0;
  exportValToLocalBufferLocs_[ExportVal2_::FwdRow1].bufType[0] = BufferTypeCompact_::ROW;

  exportValToLocalBufferLocs_[ExportVal2_::BackRow2].nBuffers = 1;
  exportValToLocalBufferLocs_[ExportVal2_::BackRow2].sectNum[0] = 2;
  exportValToLocalBufferLocs_[ExportVal2_::BackRow2].bufType[0] = BufferTypeCompact_::ROW;

  exportValToLocalBufferLocs_[ExportVal2_::FwdRow2].nBuffers = 1;
  exportValToLocalBufferLocs_[ExportVal2_::FwdRow2].sectNum[0] = 3;
  exportValToLocalBufferLocs_[ExportVal2_::FwdRow2].bufType[0] = BufferTypeCompact_::ROW;

  exportValToLocalBufferLocs_[ExportVal2_::BackRowOverlap].nBuffers = 2;
  exportValToLocalBufferLocs_[ExportVal2_::BackRowOverlap].sectNum[0] = 1;
  exportValToLocalBufferLocs_[ExportVal2_::BackRowOverlap].bufType[0] = BufferTypeCompact_::ROW;
  exportValToLocalBufferLocs_[ExportVal2_::BackRowOverlap].sectNum[1] = 2;
  exportValToLocalBufferLocs_[ExportVal2_::BackRowOverlap].bufType[1] = BufferTypeCompact_::ROW;

  exportValToLocalBufferLocs_[ExportVal2_::FwdRowOverlap].nBuffers = 2;
  exportValToLocalBufferLocs_[ExportVal2_::FwdRowOverlap].sectNum[0] = 0;
  exportValToLocalBufferLocs_[ExportVal2_::FwdRowOverlap].bufType[0] = BufferTypeCompact_::ROW;
  exportValToLocalBufferLocs_[ExportVal2_::FwdRowOverlap].sectNum[1] = 3;
  exportValToLocalBufferLocs_[ExportVal2_::FwdRowOverlap].bufType[1] = BufferTypeCompact_::ROW;

  exportValToLocalBufferLocs_[ExportVal2_::BackCol1].nBuffers = 1;
  exportValToLocalBufferLocs_[ExportVal2_::BackCol1].sectNum[0] = 3;
  exportValToLocalBufferLocs_[ExportVal2_::BackCol1].bufType[0] = BufferTypeCompact_::COL;

  exportValToLocalBufferLocs_[ExportVal2_::FwdCol1].nBuffers = 1;
  exportValToLocalBufferLocs_[ExportVal2_::FwdCol1].sectNum[0] = 0;
  exportValToLocalBufferLocs_[ExportVal2_::FwdCol1].bufType[0] = BufferTypeCompact_::COL;

  exportValToLocalBufferLocs_[ExportVal2_::BackCol2].nBuffers = 1;
  exportValToLocalBufferLocs_[ExportVal2_::BackCol2].sectNum[0] = 2;
  exportValToLocalBufferLocs_[ExportVal2_::BackCol2].bufType[0] = BufferTypeCompact_::COL;

  exportValToLocalBufferLocs_[ExportVal2_::FwdCol2].nBuffers = 1;
  exportValToLocalBufferLocs_[ExportVal2_::FwdCol2].sectNum[0] = 1;
  exportValToLocalBufferLocs_[ExportVal2_::FwdCol2].bufType[0] = BufferTypeCompact_::COL;

  exportValToLocalBufferLocs_[ExportVal2_::BackColOverlap].nBuffers = 2;
  exportValToLocalBufferLocs_[ExportVal2_::BackColOverlap].sectNum[0] = 3;
  exportValToLocalBufferLocs_[ExportVal2_::BackColOverlap].bufType[0] = BufferTypeCompact_::COL;
  exportValToLocalBufferLocs_[ExportVal2_::BackColOverlap].sectNum[1] = 2;
  exportValToLocalBufferLocs_[ExportVal2_::BackColOverlap].bufType[1] = BufferTypeCompact_::COL;

  exportValToLocalBufferLocs_[ExportVal2_::FwdColOverlap].nBuffers = 2;
  exportValToLocalBufferLocs_[ExportVal2_::FwdColOverlap].sectNum[0] = 0;
  exportValToLocalBufferLocs_[ExportVal2_::FwdColOverlap].bufType[0] = BufferTypeCompact_::COL;
  exportValToLocalBufferLocs_[ExportVal2_::FwdColOverlap].sectNum[1] = 1;
  exportValToLocalBufferLocs_[ExportVal2_::FwdColOverlap].bufType[1] = BufferTypeCompact_::COL;

  exportValToLocalBufferLocs_[ExportVal2_::BackRowBackCol].nBuffers = 3;
  exportValToLocalBufferLocs_[ExportVal2_::BackRowBackCol].sectNum[0] = 1;
  exportValToLocalBufferLocs_[ExportVal2_::BackRowBackCol].sectNum[1] = 3;
  exportValToLocalBufferLocs_[ExportVal2_::BackRowBackCol].sectNum[2] = 2;
  exportValToLocalBufferLocs_[ExportVal2_::BackRowBackCol].bufType[0] = BufferTypeCompact_::ROW;
  exportValToLocalBufferLocs_[ExportVal2_::BackRowBackCol].bufType[1] = BufferTypeCompact_::COL;
  exportValToLocalBufferLocs_[ExportVal2_::BackRowBackCol].bufType[2] = BufferTypeCompact_::CORNER;

  exportValToLocalBufferLocs_[ExportVal2_::FwdRowFwdCol].nBuffers = 3;
  exportValToLocalBufferLocs_[ExportVal2_::FwdRowFwdCol].sectNum[0] = 3;
  exportValToLocalBufferLocs_[ExportVal2_::FwdRowFwdCol].sectNum[1] = 1;
  exportValToLocalBufferLocs_[ExportVal2_::FwdRowFwdCol].sectNum[2] = 0;
  exportValToLocalBufferLocs_[ExportVal2_::FwdRowFwdCol].bufType[0] = BufferTypeCompact_::ROW;
  exportValToLocalBufferLocs_[ExportVal2_::FwdRowFwdCol].bufType[1] = BufferTypeCompact_::COL;
  exportValToLocalBufferLocs_[ExportVal2_::FwdRowFwdCol].bufType[2] = BufferTypeCompact_::CORNER;

  exportValToLocalBufferLocs_[ExportVal2_::BackRowFwdCol].nBuffers = 3;
  exportValToLocalBufferLocs_[ExportVal2_::BackRowFwdCol].sectNum[0] = 2;
  exportValToLocalBufferLocs_[ExportVal2_::BackRowFwdCol].sectNum[1] = 0;
  exportValToLocalBufferLocs_[ExportVal2_::BackRowFwdCol].sectNum[2] = 1;
  exportValToLocalBufferLocs_[ExportVal2_::BackRowFwdCol].bufType[0] = BufferTypeCompact_::ROW;
  exportValToLocalBufferLocs_[ExportVal2_::BackRowFwdCol].bufType[1] = BufferTypeCompact_::COL;
  exportValToLocalBufferLocs_[ExportVal2_::BackRowFwdCol].bufType[2] = BufferTypeCompact_::CORNER;

  exportValToLocalBufferLocs_[ExportVal2_::FwdRowBackCol].nBuffers = 3;
  exportValToLocalBufferLocs_[ExportVal2_::FwdRowBackCol].sectNum[0] = 0;
  exportValToLocalBufferLocs_[ExportVal2_::FwdRowBackCol].sectNum[1] = 2;
  exportValToLocalBufferLocs_[ExportVal2_::FwdRowBackCol].sectNum[2] = 3;
  exportValToLocalBufferLocs_[ExportVal2_::FwdRowBackCol].bufType[0] = BufferTypeCompact_::ROW;
  exportValToLocalBufferLocs_[ExportVal2_::FwdRowBackCol].bufType[1] = BufferTypeCompact_::COL;
  exportValToLocalBufferLocs_[ExportVal2_::FwdRowBackCol].bufType[2] = BufferTypeCompact_::CORNER;

  // MPI datatypes
  // =============

  MPI_Datatype contigInt, contigFloat;

  MPI_Type_contiguous(nIntsPerCell_, MPI_INT, &contigInt);
  MPI_Type_contiguous(nFloatsPerCell_, MPI_DOUBLE, &contigFloat);
  MPI_Type_commit(&contigInt);
  MPI_Type_commit(&contigFloat);

  MPI_Type_vector(ghostExtent_[0], ghostExtent_[1], extentWGhost_[1], contigInt, &ghostCornerInt_);
  MPI_Type_vector(ghostExtent_[0], ghostExtent_[1], extentWGhost_[1], contigFloat, &ghostCornerFloat_);
  MPI_Type_commit(&ghostCornerInt_);
  MPI_Type_commit(&ghostCornerFloat_);

  MPI_Type_vector(ghostExtent_[0], firstHalfLocalDims[1] + ghostExtent_[1], extentWGhost_[1], contigInt, &ghostHalfRowInt_[0]);
  MPI_Type_vector(ghostExtent_[0], firstHalfLocalDims[1] + ghostExtent_[1], extentWGhost_[1], contigFloat, &ghostHalfRowFloat_[0]);
  MPI_Type_vector(ghostExtent_[0], secondHalfLocalDims[1] + ghostExtent_[1], extentWGhost_[1], contigInt, &ghostHalfRowInt_[1]);
  MPI_Type_vector(ghostExtent_[0], secondHalfLocalDims[1] + ghostExtent_[1], extentWGhost_[1], contigFloat, &ghostHalfRowFloat_[1]);

  MPI_Type_vector(firstHalfLocalDims[0] + ghostExtent_[0], ghostExtent_[1], extentWGhost_[1], contigInt, &ghostHalfColInt_[0]);
  MPI_Type_vector(firstHalfLocalDims[0] + ghostExtent_[0], ghostExtent_[1], extentWGhost_[1], contigFloat, &ghostHalfColFloat_[0]);
  MPI_Type_vector(secondHalfLocalDims[0] + ghostExtent_[0], ghostExtent_[1], extentWGhost_[1], contigInt, &ghostHalfColInt_[1]);
  MPI_Type_vector(secondHalfLocalDims[0] + ghostExtent_[0], ghostExtent_[1], extentWGhost_[1], contigFloat, &ghostHalfColFloat_[1]);

  for (int i = 0; i < 2; ++i) {
    MPI_Type_commit(&ghostHalfRowInt_[i]);
    MPI_Type_commit(&ghostHalfRowFloat_[i]);
    MPI_Type_commit(&ghostHalfColInt_[i]);
    MPI_Type_commit(&ghostHalfColFloat_[i]);
  }

}
#endif

#if KMC_PARALLEL
void Lattice::Impl_::setUpMPICartForRowPartitions_(MPI_Comm latticeCommInitial, std::vector<int> & pIdNeighs) {
  procsPerDim_[0] = nProcs_;
  procsPerDim_[1] = 1;

  int isPeriodic = 1;

  MPI_Cart_create(latticeCommInitial, 1, &procsPerDim_[0],
		  &isPeriodic, 1, &latticeCommCart_);

  MPI_Comm_rank(latticeCommCart_, &procID_);
  MPI_Cart_coords(latticeCommCart_, procID_, 1, &commCoords_[0]);
  commCoords_[1] = 0;
  
  pIdNeighs.resize(ProcIdCart1_::SIZE);

  MPI_Cart_shift(latticeCommCart_, 0, 1, &pIdNeighs[ProcIdCart1_::Back], &pIdNeighs[ProcIdCart1_::Fwd]);
}
#endif

#if KMC_PARALLEL
void Lattice::Impl_::setUpSectoringForRowPartitions_(const std::vector<int> & pIdNeighs) {
  nSectors_ = 2;

  wrapIndsIfNeeded_ = &Impl_::wrapJOnly_;

  sendGhosts_ = &Impl_::sendGhostsRow_;
  recvGhosts_ = &Impl_::recvGhostsRow_;

  sendGhostsUpdate_ = &Impl_::sendGhostsUpdateRow_;
  recvGhostsUpdate_ = &Impl_::recvGhostsUpdateRow_;

  reExportIfNeeded_ = &Impl_::reExportIfNeededFake_;

  iminS_.resize(nSectors_);
  imaxP1S_.resize(nSectors_);
  jminS_.resize(nSectors_);
  jmaxP1S_.resize(nSectors_);

  ghostRecvBufferBounds_.resize(nSectors_);
  localRecvBufferBounds_.resize(nSectors_);

  for (int i = 0; i < nSectors_; ++i) {
    ghostRecvBufferBounds_[i].resize(BufferTypeRow_::SIZE);
    localRecvBufferBounds_[i].resize(BufferTypeRow_::SIZE);
  }

  int firstHalfLocalDims = localDims_[0]/2;

  std::size_t reservedRow = localDims_[1]*lattice_.capacity();

  ghostSendIndsBuffer_.resize(BufferTypeRow_::SIZE);
  ghostSendIndsBuffer_[BufferTypeRow_::ROW].reserve(reservedRow);

  ghostSendIntBuffer_.resize(BufferTypeRow_::SIZE);
  ghostSendIntBuffer_[BufferTypeRow_::ROW].reserve(nIntsPerCell_*reservedRow);

  ghostSendFloatBuffer_.resize(BufferTypeRow_::SIZE);
  ghostSendFloatBuffer_[BufferTypeRow_::ROW].reserve(nFloatsPerCell_*reservedRow);

  localSendIndsBuffer_.resize(nSectors_);
  localSendIntBuffer_.resize(nSectors_);
  localSendFloatBuffer_.resize(nSectors_);

  for (int i = 0; i < nSectors_; ++i) {

    localSendIndsBuffer_[i].resize(BufferTypeRow_::SIZE);
    localSendIndsBuffer_[i][BufferTypeRow_::ROW].reserve(reservedRow);

    localSendIntBuffer_[i].resize(BufferTypeRow_::SIZE);
    localSendIntBuffer_[i][BufferTypeRow_::ROW].reserve(nIntsPerCell_*reservedRow);

    localSendFloatBuffer_[i].resize(BufferTypeRow_::SIZE);
    localSendFloatBuffer_[i][BufferTypeRow_::ROW].reserve(nFloatsPerCell_*reservedRow);

  }

  ghostRecvIndsBuffer_.resize(BufferTypeRow_::SIZE);
  localRecvIndsBuffer_.resize(BufferTypeRow_::SIZE);
  ghostRecvIntBuffer_.resize(BufferTypeRow_::SIZE);
  localRecvIntBuffer_.resize(BufferTypeRow_::SIZE);
  ghostRecvFloatBuffer_.resize(BufferTypeRow_::SIZE);
  localRecvFloatBuffer_.resize(BufferTypeRow_::SIZE);

  boost::array<LatticePlanarBBox, ExportVal1_::SIZE> ghostExportBBox, localExportBBox;

  // ghostExportBBox[ExportVal1_::NONE] is just dummy values; I just
  // need to make sure that imin = imaxP1 so that the following "for"
  // loop works.
  ghostExportBBox[ExportVal1_::NONE].imin = ghostExportBBox[ExportVal1_::NONE].imaxP1 = 0;

  ghostExportBBox[ExportVal1_::Back].imin = globalOffset_[0] - ghostExtent_[0];
  ghostExportBBox[ExportVal1_::Back].jmin = globalOffset_[1];
  ghostExportBBox[ExportVal1_::Back].imaxP1 = ghostExportBBox[ExportVal1_::Back].imin + ghostExtent_[0];
  ghostExportBBox[ExportVal1_::Back].jmaxP1 = ghostExportBBox[ExportVal1_::Back].jmin + localDims_[1];

  ghostExportBBox[ExportVal1_::Fwd].imin = globalOffset_[0] + localDims_[0];
  ghostExportBBox[ExportVal1_::Fwd].jmin = globalOffset_[1];
  ghostExportBBox[ExportVal1_::Fwd].imaxP1 = ghostExportBBox[ExportVal1_::Fwd].imin + ghostExtent_[0];
  ghostExportBBox[ExportVal1_::Fwd].jmaxP1 = ghostExportBBox[ExportVal1_::Fwd].jmin + localDims_[1];

  localExportBBox[ExportVal1_::NONE].imin = globalOffset_[0] + ghostExtent_[0];
  localExportBBox[ExportVal1_::NONE].jmin = globalOffset_[1];
  localExportBBox[ExportVal1_::NONE].imaxP1 = globalOffset_[0] + localDims_[0] - ghostExtent_[0];
  localExportBBox[ExportVal1_::NONE].jmaxP1 = globalOffset_[1] + localDims_[1];

  localExportBBox[ExportVal1_::Back].imin = globalOffset_[0];
  localExportBBox[ExportVal1_::Back].jmin = globalOffset_[1];
  localExportBBox[ExportVal1_::Back].imaxP1 = localExportBBox[ExportVal1_::Back].imin + ghostExtent_[0];
  localExportBBox[ExportVal1_::Back].jmaxP1 = localExportBBox[ExportVal1_::Back].jmin + localDims_[1];

  localExportBBox[ExportVal1_::Fwd].imin = globalOffset_[0] + localDims_[0] - ghostExtent_[0];
  localExportBBox[ExportVal1_::Fwd].jmin = globalOffset_[1];
  localExportBBox[ExportVal1_::Fwd].imaxP1 = localExportBBox[ExportVal1_::Fwd].imin + ghostExtent_[0];
  localExportBBox[ExportVal1_::Fwd].jmaxP1 = localExportBBox[ExportVal1_::Fwd].jmin + localDims_[1];

  for (int exportVal = 0; exportVal < ExportVal1_::SIZE; ++exportVal) {
    
    for (int i = ghostExportBBox[exportVal].imin; i < ghostExportBBox[exportVal].imaxP1; ++i) {
      for (int j = ghostExportBBox[exportVal].jmin; j < ghostExportBBox[exportVal].jmaxP1; ++j) {
        cellParInfoArray_[i][j].sectNum = -1;
        cellParInfoArray_[i][j].exportVal = exportVal;
      }
    }

    for (int i = localExportBBox[exportVal].imin; i < localExportBBox[exportVal].imaxP1; ++i) {
      for (int j = localExportBBox[exportVal].jmin; j < localExportBBox[exportVal].jmaxP1; ++j) {
        cellParInfoArray_[i][j].exportVal = exportVal;
      }
    }

  }

  sendRecvGhostRowInfo_.resize(nSectors_);

  // Sector 0
  // ========

  sendRecvGhostRowInfo_[0].sendRank = pIdNeighs[ProcIdCart1_::Fwd];
  sendRecvGhostRowInfo_[0].recvRank = pIdNeighs[ProcIdCart1_::Back];

  sendRecvGhostRowInfo_[0].sendPlanarCoords = localDims_[0] - ghostExtent_[0] + globalOffset_[0];
  sendRecvGhostRowInfo_[0].recvPlanarCoords = -ghostExtent_[0] + globalOffset_[0];

  localRecvBufferBounds_[0][BufferTypeRow_::ROW].imin = sendRecvGhostRowInfo_[0].sendPlanarCoords;
  localRecvBufferBounds_[0][BufferTypeRow_::ROW].jmin = 0;

  localRecvBufferBounds_[0][BufferTypeRow_::ROW].imaxP1 = localRecvBufferBounds_[0][BufferTypeRow_::ROW].imin + ghostExtent_[0];
  localRecvBufferBounds_[0][BufferTypeRow_::ROW].jmaxP1 = localRecvBufferBounds_[0][BufferTypeRow_::ROW].jmin + localDims_[1];
  
  ghostRecvBufferBounds_[0][BufferTypeRow_::ROW].imin = sendRecvGhostRowInfo_[0].recvPlanarCoords;
  localRecvBufferBounds_[0][BufferTypeRow_::ROW].jmin = 0;

  ghostRecvBufferBounds_[0][BufferTypeRow_::ROW].imaxP1 = ghostRecvBufferBounds_[0][BufferTypeRow_::ROW].imin + ghostExtent_[0];
  ghostRecvBufferBounds_[0][BufferTypeRow_::ROW].jmaxP1 = ghostRecvBufferBounds_[0][BufferTypeRow_::ROW].jmin + localDims_[1];

  iminS_[0] = globalOffset_[0];
  jminS_[0] = globalOffset_[1];

  imaxP1S_[0] = firstHalfLocalDims + globalOffset_[0];
  jmaxP1S_[0] = localDims_[1] + globalOffset_[1];

  // Sector 1
  // ========

  sendRecvGhostRowInfo_[1].sendRank = pIdNeighs[ProcIdCart1_::Back];
  sendRecvGhostRowInfo_[1].recvRank = pIdNeighs[ProcIdCart1_::Fwd];

  sendRecvGhostRowInfo_[1].sendPlanarCoords = globalOffset_[0];
  sendRecvGhostRowInfo_[1].recvPlanarCoords = localDims_[0] + globalOffset_[0];

  localRecvBufferBounds_[1][BufferTypeRow_::ROW].imin = sendRecvGhostRowInfo_[1].sendPlanarCoords;
  localRecvBufferBounds_[1][BufferTypeRow_::ROW].jmin = 0;

  localRecvBufferBounds_[1][BufferTypeRow_::ROW].imaxP1 = localRecvBufferBounds_[1][BufferTypeRow_::ROW].imin + ghostExtent_[0];
  localRecvBufferBounds_[1][BufferTypeRow_::ROW].jmaxP1 = localRecvBufferBounds_[1][BufferTypeRow_::ROW].jmin + localDims_[1];
  
  ghostRecvBufferBounds_[1][BufferTypeRow_::ROW].imin = sendRecvGhostRowInfo_[1].recvPlanarCoords;
  localRecvBufferBounds_[1][BufferTypeRow_::ROW].jmin = 0;

  ghostRecvBufferBounds_[1][BufferTypeRow_::ROW].imaxP1 = ghostRecvBufferBounds_[1][BufferTypeRow_::ROW].imin + ghostExtent_[0];
  ghostRecvBufferBounds_[1][BufferTypeRow_::ROW].jmaxP1 = ghostRecvBufferBounds_[1][BufferTypeRow_::ROW].jmin + localDims_[1];

  iminS_[1] = firstHalfLocalDims + globalOffset_[0];
  jminS_[1] = globalOffset_[1];
  imaxP1S_[1] = localDims_[0] + globalOffset_[0];
  jmaxP1S_[1] = localDims_[1] + globalOffset_[1];

  extentInt_ = ghostExtent_[0]*localDims_[1]*nIntsPerCell_;
  extentFloat_ = ghostExtent_[0]*localDims_[1]*nFloatsPerCell_;

  for (int sectNum = 0; sectNum < nSectors_; ++sectNum) {
    for (int i = iminS_[sectNum]; i < imaxP1S_[sectNum]; ++i) {
      for (int j = jminS_[sectNum]; j < jmaxP1S_[sectNum]; ++j) {
        cellParInfoArray_[i][j].sectNum = sectNum;
      }
    }
  }
  
  exportValToGhostBufType_.resize(ExportVal1_::SIZE);
  exportValToLocalBufferLocs_.resize(ExportVal1_::SIZE);

  exportValToGhostBufType_[ExportVal1_::NONE] = -42; // Garbage value, should never get used.  
  exportValToGhostBufType_[ExportVal1_::Back] = BufferTypeRow_::ROW;  
  exportValToGhostBufType_[ExportVal1_::Fwd] = BufferTypeRow_::ROW;
  
  exportValToLocalBufferLocs_[ExportVal1_::NONE].nBuffers = 0;

  exportValToLocalBufferLocs_[ExportVal1_::Back].nBuffers = 1;
  exportValToLocalBufferLocs_[ExportVal1_::Back].sectNum[0] = 1;
  exportValToLocalBufferLocs_[ExportVal1_::Back].bufType[0] = BufferTypeRow_::ROW;
  
  exportValToLocalBufferLocs_[ExportVal1_::Fwd].nBuffers = 1;
  exportValToLocalBufferLocs_[ExportVal1_::Fwd].sectNum[0] = 0;
  exportValToLocalBufferLocs_[ExportVal1_::Fwd].bufType[0] = BufferTypeRow_::ROW;

}
#endif

#if KMC_PARALLEL
void Lattice::Impl_::recvGhostsCompact_(int sectNum, int newLatticeHeight) {

  const SendRecvGhostCompactInfo_ & srInfo = sendRecvGhostCompactInfo_[sectNum];

  for (int i = 0; i < newLatticeHeight; ++i) {

    const int tagCornerInt = 50, tagHalfRowInt = 51, tagHalfColInt = 52;
    const int tagCornerFloat = 60, tagHalfRowFloat = 61, tagHalfColFloat = 62;
  
    //MPI_Barrier(latticeCommCart_);

    boost::array<MPI_Status, 12> status;
    boost::array<MPI_Request, 12> request;

    int requestInd = 0;

    if (nIntsPerCell_ > 0) {      

      MPI_Isend(&(lattice_[i].intArray_[srInfo.sendCornerPlanarCoords[0]][srInfo.sendCornerPlanarCoords[1]][0]), 1,
              ghostCornerInt_, srInfo.sendCornerRank, tagCornerInt, latticeCommCart_, &request[requestInd++]);
      MPI_Irecv(&(lattice_[i].intArray_[srInfo.recvCornerPlanarCoords[0]][srInfo.recvCornerPlanarCoords[1]][0]), 1,
              ghostCornerInt_, srInfo.recvCornerRank, tagCornerInt, latticeCommCart_, &request[requestInd++]);

      MPI_Isend(&(lattice_[i].intArray_[srInfo.sendHalfRowPlanarCoords[0]][srInfo.sendHalfRowPlanarCoords[1]][0]), 1,
              ghostHalfRowInt_[srInfo.whichHalfRow], srInfo.sendHalfRowRank, tagHalfRowInt, latticeCommCart_, &request[requestInd++]);
      MPI_Irecv(&(lattice_[i].intArray_[srInfo.recvHalfRowPlanarCoords[0]][srInfo.recvHalfRowPlanarCoords[1]][0]), 1,
              ghostHalfRowInt_[srInfo.whichHalfRow], srInfo.recvHalfRowRank, tagHalfRowInt, latticeCommCart_, &request[requestInd++]);

      MPI_Isend(&(lattice_[i].intArray_[srInfo.sendHalfColPlanarCoords[0]][srInfo.sendHalfColPlanarCoords[1]][0]), 1,
              ghostHalfColInt_[srInfo.whichHalfCol], srInfo.sendHalfColRank, tagHalfColInt, latticeCommCart_, &request[requestInd++]);
      MPI_Irecv(&(lattice_[i].intArray_[srInfo.recvHalfColPlanarCoords[0]][srInfo.recvHalfColPlanarCoords[1]][0]), 1,
              ghostHalfColInt_[srInfo.whichHalfCol], srInfo.recvHalfColRank, tagHalfColInt, latticeCommCart_, &request[requestInd++]);
    }

    if (nFloatsPerCell_ > 0) {
      
      MPI_Isend(&(lattice_[i].floatArray_[srInfo.sendCornerPlanarCoords[0]][srInfo.sendCornerPlanarCoords[1]][0]), 1,
              ghostCornerFloat_, srInfo.sendCornerRank, tagCornerFloat, latticeCommCart_, &request[requestInd++]);
      MPI_Irecv(&(lattice_[i].floatArray_[srInfo.recvCornerPlanarCoords[0]][srInfo.recvCornerPlanarCoords[1]][0]), 1,
              ghostCornerFloat_, srInfo.recvCornerRank, tagCornerFloat, latticeCommCart_, &request[requestInd++]);

      MPI_Isend(&(lattice_[i].floatArray_[srInfo.sendHalfRowPlanarCoords[0]][srInfo.sendHalfRowPlanarCoords[1]][0]), 1,
              ghostHalfRowFloat_[srInfo.whichHalfRow], srInfo.sendHalfRowRank, tagHalfRowFloat, latticeCommCart_, &request[requestInd++]);
      MPI_Irecv(&(lattice_[i].floatArray_[srInfo.recvHalfRowPlanarCoords[0]][srInfo.recvHalfRowPlanarCoords[1]][0]), 1,
              ghostHalfRowFloat_[srInfo.whichHalfRow], srInfo.recvHalfRowRank, tagHalfRowFloat, latticeCommCart_, &request[requestInd++]);

      MPI_Isend(&(lattice_[i].floatArray_[srInfo.sendHalfColPlanarCoords[0]][srInfo.sendHalfColPlanarCoords[1]][0]), 1,
              ghostHalfColFloat_[srInfo.whichHalfCol], srInfo.sendHalfColRank, tagHalfColFloat, latticeCommCart_, &request[requestInd++]);
      MPI_Irecv(&(lattice_[i].floatArray_[srInfo.recvHalfColPlanarCoords[0]][srInfo.recvHalfColPlanarCoords[1]][0]), 1,
              ghostHalfColFloat_[srInfo.whichHalfCol], srInfo.recvHalfColRank, tagHalfColFloat, latticeCommCart_, &request[requestInd++]);
    }

    MPI_Waitall(requestInd, &request[0], &status[0]);
  }
}
#endif

#if KMC_PARALLEL
void Lattice::Impl_::recvGhostsUpdateCompact_(int sectNum) {

  const SendRecvGhostCompactInfo_ & srInfo = sendRecvGhostCompactInfo_[sectNum];

  const int tagCornerInt = 50, tagHalfRowInt = 51, tagHalfColInt = 52;
  const int tagCornerFloat = 60, tagHalfRowFloat = 61, tagHalfColFloat = 62;
  const int tagCornerInds = 70, tagHalfRowInds = 71, tagHalfColInds = 72; 
  
  //MPI_Barrier(latticeCommCart_);

  boost::array<MPI_Status, 18> status;
  boost::array<MPI_Request, 18> request;

  int requestInd = 0;

#ifndef MPI_ISEND_NOT_CONST_CORRECT /* If the MPI implementation doesn't allow 
                                      MPI_Isend to have a constant first argument
                                      (even though it should), don't use the "const"
                                      modifier. */
  const
#endif
  std::vector<std::vector<IJK> > & localSendIndsBuffer = localSendIndsBuffer_[sectNum];

  boost::array<int, BufferTypeCompact_::SIZE> localSendIndsBufferSize, ghostRecvIndsBufferSize;

  for (std::size_t i = 0; i < BufferTypeCompact_::SIZE; ++i) {
    localSendIndsBufferSize[i] = localSendIndsBuffer[i].size();
  }
  
  MPI_Isend(&localSendIndsBufferSize[BufferTypeCompact_::CORNER], 1, MPI_INT, 
            srInfo.sendCornerRank, tagCornerInds, latticeCommCart_, &request[requestInd++]);
  MPI_Irecv(&ghostRecvIndsBufferSize[BufferTypeCompact_::CORNER], 1, MPI_INT, 
            srInfo.recvCornerRank, tagCornerInds, latticeCommCart_, &request[requestInd++]);
  
  MPI_Isend(&localSendIndsBufferSize[BufferTypeCompact_::ROW], 1, MPI_INT, 
            srInfo.sendHalfRowRank, tagHalfRowInds, latticeCommCart_, &request[requestInd++]);
  MPI_Irecv(&ghostRecvIndsBufferSize[BufferTypeCompact_::ROW], 1, MPI_INT, 
            srInfo.recvHalfRowRank, tagHalfRowInds, latticeCommCart_, &request[requestInd++]);
  
  MPI_Isend(&localSendIndsBufferSize[BufferTypeCompact_::COL], 1, MPI_INT, 
            srInfo.sendHalfColRank, tagHalfColInds, latticeCommCart_, &request[requestInd++]);
  MPI_Irecv(&ghostRecvIndsBufferSize[BufferTypeCompact_::COL], 1, MPI_INT, 
            srInfo.recvHalfColRank, tagHalfColInds, latticeCommCart_, &request[requestInd++]);

  MPI_Waitall(requestInd, &request[0], &status[0]);

  for (std::size_t i = 0; i < BufferTypeCompact_::SIZE; ++i) {
    ghostRecvIndsBuffer_[i].resize(ghostRecvIndsBufferSize[i]);
  }

  requestInd = 0; // Resetting requestInd
  
  MPI_Isend(&(localSendIndsBuffer[BufferTypeCompact_::CORNER][0]), localSendIndsBufferSize[BufferTypeCompact_::CORNER],
            ijkDType_, srInfo.sendCornerRank, tagCornerInds, latticeCommCart_, &request[requestInd++]);
  MPI_Irecv(&(ghostRecvIndsBuffer_[BufferTypeCompact_::CORNER][0]), ghostRecvIndsBufferSize[BufferTypeCompact_::CORNER],
            ijkDType_, srInfo.recvCornerRank, tagCornerInds, latticeCommCart_, &request[requestInd++]);
  
  MPI_Isend(&(localSendIndsBuffer[BufferTypeCompact_::ROW][0]), localSendIndsBufferSize[BufferTypeCompact_::ROW],
            ijkDType_, srInfo.sendHalfRowRank, tagHalfRowInds, latticeCommCart_, &request[requestInd++]);
  MPI_Irecv(&(ghostRecvIndsBuffer_[BufferTypeCompact_::ROW][0]), ghostRecvIndsBufferSize[BufferTypeCompact_::ROW],
            ijkDType_, srInfo.recvHalfRowRank, tagHalfRowInds, latticeCommCart_, &request[requestInd++]);
  
  MPI_Isend(&(localSendIndsBuffer[BufferTypeCompact_::COL][0]), localSendIndsBufferSize[BufferTypeCompact_::COL], 
            ijkDType_, srInfo.sendHalfColRank, tagHalfColInds, latticeCommCart_, &request[requestInd++]);
  MPI_Irecv(&(ghostRecvIndsBuffer_[BufferTypeCompact_::COL][0]), ghostRecvIndsBufferSize[BufferTypeCompact_::COL], 
            ijkDType_, srInfo.recvHalfColRank, tagHalfColInds, latticeCommCart_, &request[requestInd++]);

  if (nIntsPerCell_ > 0) {

    for (std::size_t i = 0; i < BufferTypeCompact_::SIZE; ++i) {
      ghostRecvIntBuffer_[i].resize(nIntsPerCell_*ghostRecvIndsBufferSize[i]);
    }
    
#ifndef MPI_ISEND_NOT_CONST_CORRECT
    const
#endif
    std::vector<std::vector<int> > & localSendIntBuffer = localSendIntBuffer_[sectNum];
  
    MPI_Isend(&(localSendIntBuffer[BufferTypeCompact_::CORNER][0]), nIntsPerCell_*localSendIndsBufferSize[BufferTypeCompact_::CORNER],
              MPI_INT, srInfo.sendCornerRank, tagCornerInt, latticeCommCart_, &request[requestInd++]);
    MPI_Irecv(&(ghostRecvIntBuffer_[BufferTypeCompact_::CORNER][0]), nIntsPerCell_*ghostRecvIndsBufferSize[BufferTypeCompact_::CORNER],
              MPI_INT, srInfo.recvCornerRank, tagCornerInt, latticeCommCart_, &request[requestInd++]);
  
    MPI_Isend(&(localSendIntBuffer[BufferTypeCompact_::ROW][0]), nIntsPerCell_*localSendIndsBufferSize[BufferTypeCompact_::ROW],
              MPI_INT, srInfo.sendHalfRowRank, tagHalfRowInt, latticeCommCart_, &request[requestInd++]);
    MPI_Irecv(&(ghostRecvIntBuffer_[BufferTypeCompact_::ROW][0]), nIntsPerCell_*ghostRecvIndsBufferSize[BufferTypeCompact_::ROW],
              MPI_INT, srInfo.recvHalfRowRank, tagHalfRowInt, latticeCommCart_, &request[requestInd++]);
  
    MPI_Isend(&(localSendIntBuffer[BufferTypeCompact_::COL][0]), nIntsPerCell_*localSendIndsBufferSize[BufferTypeCompact_::COL], 
              MPI_INT, srInfo.sendHalfColRank, tagHalfColInt, latticeCommCart_, &request[requestInd++]);
    MPI_Irecv(&(ghostRecvIntBuffer_[BufferTypeCompact_::COL][0]), nIntsPerCell_*ghostRecvIndsBufferSize[BufferTypeCompact_::COL], 
              MPI_INT, srInfo.recvHalfColRank, tagHalfColInt, latticeCommCart_, &request[requestInd++]);
  }

  if (nFloatsPerCell_ > 0) {

    for (std::size_t i = 0; i < BufferTypeCompact_::SIZE; ++i) {
      ghostRecvFloatBuffer_[i].resize(nFloatsPerCell_*ghostRecvIndsBufferSize[i]);
    }    

#ifndef MPI_ISEND_NOT_CONST_CORRECT
    const
#endif
    std::vector<std::vector<double> > & localSendFloatBuffer = localSendFloatBuffer_[sectNum];
  
    MPI_Isend(&(localSendFloatBuffer[BufferTypeCompact_::CORNER][0]), nFloatsPerCell_*localSendIndsBufferSize[BufferTypeCompact_::CORNER],
              MPI_DOUBLE, srInfo.sendCornerRank, tagCornerFloat, latticeCommCart_, &request[requestInd++]);
    MPI_Irecv(&(ghostRecvFloatBuffer_[BufferTypeCompact_::CORNER][0]), nFloatsPerCell_*ghostRecvIndsBufferSize[BufferTypeCompact_::CORNER],
              MPI_DOUBLE, srInfo.recvCornerRank, tagCornerFloat, latticeCommCart_, &request[requestInd++]);
  
    MPI_Isend(&(localSendFloatBuffer[BufferTypeCompact_::ROW][0]), nFloatsPerCell_*localSendIndsBufferSize[BufferTypeCompact_::ROW],
              MPI_DOUBLE, srInfo.sendHalfRowRank, tagHalfRowFloat, latticeCommCart_, &request[requestInd++]);
    MPI_Irecv(&(ghostRecvFloatBuffer_[BufferTypeCompact_::ROW][0]), nFloatsPerCell_*ghostRecvIndsBufferSize[BufferTypeCompact_::ROW],
              MPI_DOUBLE, srInfo.recvHalfRowRank, tagHalfRowFloat, latticeCommCart_, &request[requestInd++]);
  
    MPI_Isend(&(localSendFloatBuffer[BufferTypeCompact_::COL][0]), nFloatsPerCell_*localSendIndsBufferSize[BufferTypeCompact_::COL], 
              MPI_DOUBLE, srInfo.sendHalfColRank, tagHalfColFloat, latticeCommCart_, &request[requestInd++]);
    MPI_Irecv(&(ghostRecvFloatBuffer_[BufferTypeCompact_::COL][0]), nFloatsPerCell_*ghostRecvIndsBufferSize[BufferTypeCompact_::COL], 
              MPI_DOUBLE, srInfo.recvHalfColRank, tagHalfColFloat, latticeCommCart_, &request[requestInd++]);
  }

  MPI_Waitall(requestInd, &request[0], &status[0]);

  for (int i = 0; i < BufferTypeCompact_::SIZE; ++i) {
    localSendIndsBuffer_[sectNum][i].clear();
    localSendIntBuffer_[sectNum][i].clear();
    localSendFloatBuffer_[sectNum][i].clear();
  }

}
#endif

#if KMC_PARALLEL
void Lattice::Impl_::sendGhostsCompact_(int sectNum, int newLatticeHeight) {

  const SendRecvGhostCompactInfo_ & srInfo = sendRecvGhostCompactInfo_[sectNum];

  for (int i = 0; i < newLatticeHeight; ++i) {

    const int tagCornerInt = 50, tagHalfRowInt = 51, tagHalfColInt = 52;
    const int tagCornerFloat = 60, tagHalfRowFloat = 61, tagHalfColFloat = 62;
  
    //MPI_Barrier(latticeCommCart_);

    boost::array<MPI_Status, 12> status;
    boost::array<MPI_Request, 12> request;

    int requestInd = 0;

    if (nIntsPerCell_ > 0) {      

      MPI_Irecv(&(lattice_[i].intArray_[srInfo.sendCornerPlanarCoords[0]][srInfo.sendCornerPlanarCoords[1]][0]), 1,
              ghostCornerInt_, srInfo.sendCornerRank, tagCornerInt, latticeCommCart_, &request[requestInd++]);
      MPI_Isend(&(lattice_[i].intArray_[srInfo.recvCornerPlanarCoords[0]][srInfo.recvCornerPlanarCoords[1]][0]), 1,
              ghostCornerInt_, srInfo.recvCornerRank, tagCornerInt, latticeCommCart_, &request[requestInd++]);

      MPI_Irecv(&(lattice_[i].intArray_[srInfo.sendHalfRowPlanarCoords[0]][srInfo.sendHalfRowPlanarCoords[1]][0]), 1,
              ghostHalfRowInt_[srInfo.whichHalfRow], srInfo.sendHalfRowRank, tagHalfRowInt, latticeCommCart_, &request[requestInd++]);
      MPI_Isend(&(lattice_[i].intArray_[srInfo.recvHalfRowPlanarCoords[0]][srInfo.recvHalfRowPlanarCoords[1]][0]), 1,
              ghostHalfRowInt_[srInfo.whichHalfRow], srInfo.recvHalfRowRank, tagHalfRowInt, latticeCommCart_, &request[requestInd++]);

      MPI_Irecv(&(lattice_[i].intArray_[srInfo.sendHalfColPlanarCoords[0]][srInfo.sendHalfColPlanarCoords[1]][0]), 1,
              ghostHalfColInt_[srInfo.whichHalfCol], srInfo.sendHalfColRank, tagHalfColInt, latticeCommCart_, &request[requestInd++]);
      MPI_Isend(&(lattice_[i].intArray_[srInfo.recvHalfColPlanarCoords[0]][srInfo.recvHalfColPlanarCoords[1]][0]), 1,
              ghostHalfColInt_[srInfo.whichHalfCol], srInfo.recvHalfColRank, tagHalfColInt, latticeCommCart_, &request[requestInd++]);
    }

    if (nFloatsPerCell_ > 0) {
      
      MPI_Irecv(&(lattice_[i].floatArray_[srInfo.sendCornerPlanarCoords[0]][srInfo.sendCornerPlanarCoords[1]][0]), 1,
              ghostCornerFloat_, srInfo.sendCornerRank, tagCornerFloat, latticeCommCart_, &request[requestInd++]);
      MPI_Isend(&(lattice_[i].floatArray_[srInfo.recvCornerPlanarCoords[0]][srInfo.recvCornerPlanarCoords[1]][0]), 1,
              ghostCornerFloat_, srInfo.recvCornerRank, tagCornerFloat, latticeCommCart_, &request[requestInd++]);

      MPI_Irecv(&(lattice_[i].floatArray_[srInfo.sendHalfRowPlanarCoords[0]][srInfo.sendHalfRowPlanarCoords[1]][0]), 1,
              ghostHalfRowFloat_[srInfo.whichHalfRow], srInfo.sendHalfRowRank, tagHalfRowFloat, latticeCommCart_, &request[requestInd++]);
      MPI_Isend(&(lattice_[i].floatArray_[srInfo.recvHalfRowPlanarCoords[0]][srInfo.recvHalfRowPlanarCoords[1]][0]), 1,
              ghostHalfRowFloat_[srInfo.whichHalfRow], srInfo.recvHalfRowRank, tagHalfRowFloat, latticeCommCart_, &request[requestInd++]);

      MPI_Irecv(&(lattice_[i].floatArray_[srInfo.sendHalfColPlanarCoords[0]][srInfo.sendHalfColPlanarCoords[1]][0]), 1,
              ghostHalfColFloat_[srInfo.whichHalfCol], srInfo.sendHalfColRank, tagHalfColFloat, latticeCommCart_, &request[requestInd++]);
      MPI_Isend(&(lattice_[i].floatArray_[srInfo.recvHalfColPlanarCoords[0]][srInfo.recvHalfColPlanarCoords[1]][0]), 1,
              ghostHalfColFloat_[srInfo.whichHalfCol], srInfo.recvHalfColRank, tagHalfColFloat, latticeCommCart_, &request[requestInd++]);
    }

    MPI_Waitall(requestInd, &request[0], &status[0]);
  }

}
#endif



#if KMC_PARALLEL
void Lattice::Impl_::sendGhostsUpdateCompact_(int sectNum) {

  const SendRecvGhostCompactInfo_ & srInfo = sendRecvGhostCompactInfo_[sectNum];

  const int tagCornerInt = 50, tagHalfRowInt = 51, tagHalfColInt = 52;
  const int tagCornerFloat = 60, tagHalfRowFloat = 61, tagHalfColFloat = 62;
  const int tagCornerInds = 70, tagHalfRowInds = 71, tagHalfColInds = 72; 
  
  //MPI_Barrier(latticeCommCart_);

  boost::array<MPI_Status, 18> status;
  boost::array<MPI_Request, 18> request;

  int requestInd = 0;

  boost::array<int, BufferTypeCompact_::SIZE> localRecvIndsBufferSize, ghostSendIndsBufferSize;

  for (std::size_t i = 0; i < BufferTypeCompact_::SIZE; ++i) {
    ghostSendIndsBufferSize[i] = ghostSendIndsBuffer_[i].size();
  }
  
  MPI_Irecv(&localRecvIndsBufferSize[BufferTypeCompact_::CORNER], 1, MPI_INT, 
            srInfo.sendCornerRank, tagCornerInds, latticeCommCart_, &request[requestInd++]);
  MPI_Isend(&ghostSendIndsBufferSize[BufferTypeCompact_::CORNER], 1, MPI_INT, 
            srInfo.recvCornerRank, tagCornerInds, latticeCommCart_, &request[requestInd++]);
  
  MPI_Irecv(&localRecvIndsBufferSize[BufferTypeCompact_::ROW], 1, MPI_INT, 
            srInfo.sendHalfRowRank, tagHalfRowInds, latticeCommCart_, &request[requestInd++]);
  MPI_Isend(&ghostSendIndsBufferSize[BufferTypeCompact_::ROW], 1, MPI_INT, 
            srInfo.recvHalfRowRank, tagHalfRowInds, latticeCommCart_, &request[requestInd++]);
  
  MPI_Irecv(&localRecvIndsBufferSize[BufferTypeCompact_::COL], 1, MPI_INT, 
            srInfo.sendHalfColRank, tagHalfColInds, latticeCommCart_, &request[requestInd++]);
  MPI_Isend(&ghostSendIndsBufferSize[BufferTypeCompact_::COL], 1, MPI_INT, 
            srInfo.recvHalfColRank, tagHalfColInds, latticeCommCart_, &request[requestInd++]);

  MPI_Waitall(requestInd, &request[0], &status[0]);

  for (std::size_t i = 0; i < BufferTypeCompact_::SIZE; ++i) {
    localRecvIndsBuffer_[i].resize(localRecvIndsBufferSize[i]);
  }

  requestInd = 0; // Resetting requestInd
  
  MPI_Irecv(&(localRecvIndsBuffer_[BufferTypeCompact_::CORNER][0]), localRecvIndsBufferSize[BufferTypeCompact_::CORNER],
            ijkDType_, srInfo.sendCornerRank, tagCornerInds, latticeCommCart_, &request[requestInd++]);
  MPI_Isend(&(ghostSendIndsBuffer_[BufferTypeCompact_::CORNER][0]), ghostSendIndsBufferSize[BufferTypeCompact_::CORNER],
            ijkDType_, srInfo.recvCornerRank, tagCornerInds, latticeCommCart_, &request[requestInd++]);
  
  MPI_Irecv(&(localRecvIndsBuffer_[BufferTypeCompact_::ROW][0]), localRecvIndsBufferSize[BufferTypeCompact_::ROW],
            ijkDType_, srInfo.sendHalfRowRank, tagHalfRowInds, latticeCommCart_, &request[requestInd++]);
  MPI_Isend(&(ghostSendIndsBuffer_[BufferTypeCompact_::ROW][0]), ghostSendIndsBufferSize[BufferTypeCompact_::ROW],
            ijkDType_, srInfo.recvHalfRowRank, tagHalfRowInds, latticeCommCart_, &request[requestInd++]);
  
  MPI_Irecv(&(localRecvIndsBuffer_[BufferTypeCompact_::COL][0]), localRecvIndsBufferSize[BufferTypeCompact_::COL], 
            ijkDType_, srInfo.sendHalfColRank, tagHalfColInds, latticeCommCart_, &request[requestInd++]);
  MPI_Isend(&(ghostSendIndsBuffer_[BufferTypeCompact_::COL][0]), ghostSendIndsBufferSize[BufferTypeCompact_::COL], 
            ijkDType_, srInfo.recvHalfColRank, tagHalfColInds, latticeCommCart_, &request[requestInd++]);

  if (nIntsPerCell_ > 0) {

    for (std::size_t i = 0; i < BufferTypeCompact_::SIZE; ++i) {
      localRecvIntBuffer_[i].resize(nIntsPerCell_*localRecvIndsBufferSize[i]);
    }
  
    MPI_Irecv(&(localRecvIntBuffer_[BufferTypeCompact_::CORNER][0]), nIntsPerCell_*localRecvIndsBufferSize[BufferTypeCompact_::CORNER],
              MPI_INT, srInfo.sendCornerRank, tagCornerInt, latticeCommCart_, &request[requestInd++]);
    MPI_Isend(&(ghostSendIntBuffer_[BufferTypeCompact_::CORNER][0]), nIntsPerCell_*ghostSendIndsBufferSize[BufferTypeCompact_::CORNER],
              MPI_INT, srInfo.recvCornerRank, tagCornerInt, latticeCommCart_, &request[requestInd++]);
  
    MPI_Irecv(&(localRecvIntBuffer_[BufferTypeCompact_::ROW][0]), nIntsPerCell_*localRecvIndsBufferSize[BufferTypeCompact_::ROW],
              MPI_INT, srInfo.sendHalfRowRank, tagHalfRowInt, latticeCommCart_, &request[requestInd++]);
    MPI_Isend(&(ghostSendIntBuffer_[BufferTypeCompact_::ROW][0]), nIntsPerCell_*ghostSendIndsBufferSize[BufferTypeCompact_::ROW],
              MPI_INT, srInfo.recvHalfRowRank, tagHalfRowInt, latticeCommCart_, &request[requestInd++]);
  
    MPI_Irecv(&(localRecvIntBuffer_[BufferTypeCompact_::COL][0]), nIntsPerCell_*localRecvIndsBufferSize[BufferTypeCompact_::COL], 
              MPI_INT, srInfo.sendHalfColRank, tagHalfColInt, latticeCommCart_, &request[requestInd++]);
    MPI_Isend(&(ghostSendIntBuffer_[BufferTypeCompact_::COL][0]), nIntsPerCell_*ghostSendIndsBufferSize[BufferTypeCompact_::COL], 
              MPI_INT, srInfo.recvHalfColRank, tagHalfColInt, latticeCommCart_, &request[requestInd++]);
  }

  if (nFloatsPerCell_ > 0) {

    for (std::size_t i = 0; i < BufferTypeCompact_::SIZE; ++i) {
      localRecvFloatBuffer_[i].resize(nFloatsPerCell_*localRecvIndsBufferSize[i]);
    }
  
    MPI_Irecv(&(localRecvFloatBuffer_[BufferTypeCompact_::CORNER][0]), nFloatsPerCell_*localRecvIndsBufferSize[BufferTypeCompact_::CORNER],
              MPI_DOUBLE, srInfo.sendCornerRank, tagCornerFloat, latticeCommCart_, &request[requestInd++]);
    MPI_Isend(&(ghostSendFloatBuffer_[BufferTypeCompact_::CORNER][0]), nFloatsPerCell_*ghostSendIndsBufferSize[BufferTypeCompact_::CORNER],
              MPI_DOUBLE, srInfo.recvCornerRank, tagCornerFloat, latticeCommCart_, &request[requestInd++]);
  
    MPI_Irecv(&(localRecvFloatBuffer_[BufferTypeCompact_::ROW][0]), nFloatsPerCell_*localRecvIndsBufferSize[BufferTypeCompact_::ROW],
              MPI_DOUBLE, srInfo.sendHalfRowRank, tagHalfRowFloat, latticeCommCart_, &request[requestInd++]);
    MPI_Isend(&(ghostSendFloatBuffer_[BufferTypeCompact_::ROW][0]), nFloatsPerCell_*ghostSendIndsBufferSize[BufferTypeCompact_::ROW],
              MPI_DOUBLE, srInfo.recvHalfRowRank, tagHalfRowFloat, latticeCommCart_, &request[requestInd++]);
  
    MPI_Irecv(&(localRecvFloatBuffer_[BufferTypeCompact_::COL][0]), nFloatsPerCell_*localRecvIndsBufferSize[BufferTypeCompact_::COL], 
              MPI_DOUBLE, srInfo.sendHalfColRank, tagHalfColFloat, latticeCommCart_, &request[requestInd++]);
    MPI_Isend(&(ghostSendFloatBuffer_[BufferTypeCompact_::COL][0]), nFloatsPerCell_*ghostSendIndsBufferSize[BufferTypeCompact_::COL], 
              MPI_DOUBLE, srInfo.recvHalfColRank, tagHalfColFloat, latticeCommCart_, &request[requestInd++]);
  }

  MPI_Waitall(requestInd, &request[0], &status[0]);

  for (int i = 0; i < BufferTypeCompact_::SIZE; ++i) {
    ghostSendIndsBuffer_[i].clear();
    ghostSendIntBuffer_[i].clear();
    ghostSendFloatBuffer_[i].clear();
  }

}
#endif

#if KMC_PARALLEL
void Lattice::Impl_::recvGhostsRow_(int sectNum, int newLatticeHeight) {

  const SendRecvGhostRowInfo_ & srInfo = sendRecvGhostRowInfo_[sectNum];

  for (int i = 0; i < newLatticeHeight; ++i) {
    
    const int tagInt = 50, tagFloat = 51;

    //MPI_Barrier(latticeCommCart_);

    boost::array<MPI_Status, 4> status;
    boost::array<MPI_Request, 4> request;

    int requestInd = 0;

    if (nIntsPerCell_ > 0) {
      MPI_Isend(&(lattice_[i].intArray_[srInfo.sendPlanarCoords][0][0]), extentInt_,
		MPI_INT, srInfo.sendRank, tagInt, latticeCommCart_, &request[requestInd++]);
      MPI_Irecv(&(lattice_[i].intArray_[srInfo.recvPlanarCoords][0][0]), extentInt_,
		MPI_INT, srInfo.recvRank, tagInt, latticeCommCart_, &request[requestInd++]);
    }

    if (nFloatsPerCell_ > 0) {
      MPI_Isend(&(lattice_[i].floatArray_[srInfo.sendPlanarCoords][0][0]), extentFloat_,
		MPI_DOUBLE, srInfo.sendRank, tagFloat, latticeCommCart_, &request[requestInd++]);
      MPI_Irecv(&(lattice_[i].floatArray_[srInfo.recvPlanarCoords][0][0]), extentFloat_,
		MPI_DOUBLE, srInfo.recvRank, tagFloat, latticeCommCart_, &request[requestInd++]);            
    }

    MPI_Waitall(requestInd, &request[0], &status[0]);

  }

}
#endif

#if KMC_PARALLEL
void Lattice::Impl_::recvGhostsUpdateRow_(int sectNum) {
  const SendRecvGhostRowInfo_ & srInfo = sendRecvGhostRowInfo_[sectNum];

  const int tagInt = 50, tagFloat = 51, tagInds = 52;

  //MPI_Barrier(latticeCommCart_);

  boost::array<MPI_Status, 6> status;
  boost::array<MPI_Request, 6> request;

  int requestInd = 0;

  std::vector<IJK> & localSendIndsBuffer = localSendIndsBuffer_[sectNum][BufferTypeRow_::ROW];
  std::vector<IJK> & ghostRecvIndsBuffer = ghostRecvIndsBuffer_[BufferTypeRow_::ROW];

  int localSendIndsBufferSize = localSendIndsBuffer.size(), ghostRecvIndsBufferSize;

  MPI_Isend(&localSendIndsBufferSize, 1, MPI_INT, srInfo.sendRank, tagInds,
            latticeCommCart_, &request[requestInd++]);
  MPI_Irecv(&ghostRecvIndsBufferSize, 1, MPI_INT, srInfo.recvRank, tagInds,
            latticeCommCart_, &request[requestInd++]);

  MPI_Waitall(requestInd, &request[0], &status[0]);

  ghostRecvIndsBuffer.resize(ghostRecvIndsBufferSize);

  requestInd = 0; // Resetting requestInd

  MPI_Isend(&localSendIndsBuffer[0], localSendIndsBufferSize, ijkDType_, 
            srInfo.sendRank, tagInds, latticeCommCart_, &request[requestInd++]);
  MPI_Irecv(&ghostRecvIndsBuffer[0], ghostRecvIndsBufferSize, ijkDType_,
            srInfo.recvRank, tagInds, latticeCommCart_, &request[requestInd++]);

  if (nIntsPerCell_ > 0) {

    ghostRecvIntBuffer_[BufferTypeRow_::ROW].resize(nIntsPerCell_*ghostRecvIndsBufferSize);

    MPI_Isend(&(localSendIntBuffer_[sectNum][BufferTypeRow_::ROW][0]),
              nIntsPerCell_*localSendIndsBufferSize, MPI_INT, 
              srInfo.sendRank, tagInt, latticeCommCart_, &request[requestInd++]);
    MPI_Irecv(&ghostRecvIntBuffer_[BufferTypeRow_::ROW][0],
              nIntsPerCell_*ghostRecvIndsBufferSize, MPI_INT,
              srInfo.recvRank, tagInt, latticeCommCart_, &request[requestInd++]);
  }

  if (nFloatsPerCell_ > 0) {

    ghostRecvFloatBuffer_[BufferTypeRow_::ROW].resize(nFloatsPerCell_*ghostRecvIndsBufferSize);

    MPI_Isend(&(localSendFloatBuffer_[sectNum][BufferTypeRow_::ROW][0]),
              nFloatsPerCell_*localSendIndsBufferSize, MPI_DOUBLE, 
              srInfo.sendRank, tagFloat, latticeCommCart_, &request[requestInd++]);
    MPI_Irecv(&ghostRecvFloatBuffer_[BufferTypeRow_::ROW][0],
              nFloatsPerCell_*ghostRecvIndsBufferSize, MPI_DOUBLE,
              srInfo.recvRank, tagFloat, latticeCommCart_, &request[requestInd++]);
  }
  
  MPI_Waitall(requestInd, &request[0], &status[0]);

  localSendIndsBuffer.clear();
}
#endif

#if KMC_PARALLEL
void Lattice::Impl_::sendGhostsRow_(int sectNum, int newLatticeHeight) {

  const SendRecvGhostRowInfo_ & srInfo = sendRecvGhostRowInfo_[sectNum];

  for (int i = 0; i < newLatticeHeight; ++i) {
    
    const int tagInt = 50, tagFloat = 51;

    //MPI_Barrier(latticeCommCart_);

    boost::array<MPI_Status, 4> status;
    boost::array<MPI_Request, 4> request;

    int requestInd = 0;

    if (nIntsPerCell_ > 0) {
      MPI_Irecv(&(lattice_[i].intArray_[srInfo.sendPlanarCoords][0][0]), extentInt_,
		MPI_INT, srInfo.sendRank, tagInt, latticeCommCart_, &request[requestInd++]);
      MPI_Isend(&(lattice_[i].intArray_[srInfo.recvPlanarCoords][0][0]), extentInt_,
		MPI_INT, srInfo.recvRank, tagInt, latticeCommCart_, &request[requestInd++]);
    }

    if (nFloatsPerCell_ > 0) {
      MPI_Irecv(&(lattice_[i].floatArray_[srInfo.sendPlanarCoords][0][0]), extentFloat_,
		MPI_DOUBLE, srInfo.sendRank, tagFloat, latticeCommCart_, &request[requestInd++]);
      MPI_Isend(&(lattice_[i].floatArray_[srInfo.recvPlanarCoords][0][0]), extentFloat_,
		MPI_DOUBLE, srInfo.recvRank, tagFloat, latticeCommCart_, &request[requestInd++]);            
    }

    MPI_Waitall(requestInd, &request[0], &status[0]);

  }

}
#endif

#if KMC_PARALLEL
void Lattice::Impl_::sendGhostsUpdateRow_(int sectNum) {
  const SendRecvGhostRowInfo_ & srInfo = sendRecvGhostRowInfo_[sectNum];

  const int tagInt = 50, tagFloat = 51, tagInds = 52;

  //MPI_Barrier(latticeCommCart_);

  boost::array<MPI_Status, 6> status;
  boost::array<MPI_Request, 6> request;

  int requestInd = 0;

  std::vector<IJK> & localRecvIndsBuffer = localRecvIndsBuffer_[BufferTypeRow_::ROW];
  std::vector<IJK> & ghostSendIndsBuffer = ghostSendIndsBuffer_[BufferTypeRow_::ROW];

  int localRecvIndsBufferSize, ghostSendIndsBufferSize = ghostSendIndsBuffer.size();

  MPI_Irecv(&localRecvIndsBufferSize, 1, MPI_INT, srInfo.sendRank, tagInds,
            latticeCommCart_, &request[requestInd++]);
  MPI_Isend(&ghostSendIndsBufferSize, 1, MPI_INT, srInfo.recvRank, tagInds,
            latticeCommCart_, &request[requestInd++]);

  MPI_Waitall(requestInd, &request[0], &status[0]);

  localRecvIndsBuffer.resize(localRecvIndsBufferSize);

  requestInd = 0; // Resetting requestInd

  MPI_Irecv(&localRecvIndsBuffer[0], localRecvIndsBufferSize, ijkDType_, 
            srInfo.sendRank, tagInds, latticeCommCart_, &request[requestInd++]);
  MPI_Isend(&ghostSendIndsBuffer[0], ghostSendIndsBufferSize, ijkDType_,
            srInfo.recvRank, tagInds, latticeCommCart_, &request[requestInd++]);

  if (nIntsPerCell_ > 0) {

    localRecvIntBuffer_[BufferTypeRow_::ROW].resize(nIntsPerCell_*localRecvIndsBufferSize);

    MPI_Irecv(&(localRecvIntBuffer_[BufferTypeRow_::ROW][0]),
              nIntsPerCell_*localRecvIndsBufferSize, MPI_INT, 
              srInfo.sendRank, tagInt, latticeCommCart_, &request[requestInd++]);
    MPI_Isend(&ghostSendIntBuffer_[BufferTypeRow_::ROW][0],
              nIntsPerCell_*ghostSendIndsBufferSize, MPI_INT,
              srInfo.recvRank, tagInt, latticeCommCart_, &request[requestInd++]);
  }

  if (nFloatsPerCell_ > 0) {

    localRecvFloatBuffer_[BufferTypeRow_::ROW].resize(nFloatsPerCell_*localRecvIndsBufferSize);

    MPI_Irecv(&(localRecvFloatBuffer_[BufferTypeRow_::ROW][0]),
              nFloatsPerCell_*localRecvIndsBufferSize, MPI_DOUBLE, 
              srInfo.sendRank, tagFloat, latticeCommCart_, &request[requestInd++]);
    MPI_Isend(&ghostSendFloatBuffer_[BufferTypeRow_::ROW][0],
              nFloatsPerCell_*ghostSendIndsBufferSize, MPI_DOUBLE,
              srInfo.recvRank, tagFloat, latticeCommCart_, &request[requestInd++]);
  }
  
  MPI_Waitall(requestInd, &request[0], &status[0]);

  ghostSendIndsBuffer.clear();
}
#endif

#if KMC_PARALLEL
void Lattice::Impl_::getLocalPlanarBBox_(bool wGhost,
					 int & imin, int & imaxP1, 
					 int & jmin, int & jmaxP1) const {
  imin = globalOffset_[0];
  jmin = globalOffset_[1];

  imaxP1 = localDims_[0] + globalOffset_[0];
  jmaxP1 = localDims_[1] + globalOffset_[1];

  if (wGhost) {
    imin -= ghostExtent_[0];
    jmin -= ghostExtent_[1];

    imaxP1 += ghostExtent_[0];
    jmaxP1 += ghostExtent_[1];    
  }

}
#endif

#if KMC_PARALLEL
void Lattice::Impl_::setSendIntFloatBuffers_(const std::vector<std::vector<IJK> > & sendIndsBuffer,
                                             std::vector<std::vector<int> > & sendIntBuffer,
                                             std::vector<std::vector<double> > & sendFloatBuffer) {

  std::size_t sendIndsBufferSize = sendIndsBuffer.size();

  for (std::size_t bufType = 0; bufType < sendIndsBufferSize; ++bufType) {

    const std::vector<IJK> & currInds = sendIndsBuffer[bufType];
    std::vector<int> & currInt = sendIntBuffer[bufType];
    std::vector<double> & currFloat = sendFloatBuffer[bufType];

    currInt.clear();
    currFloat.clear();

    for (std::vector<IJK>::const_iterator itr = currInds.begin(),
           itrEnd = currInds.end(); itr != itrEnd; ++itr) {

      const IJK & ci = *itr;

      for (std::size_t whichInt = 0; whichInt < nIntsPerCell_; ++whichInt) {
        currInt.push_back(lattice_[ci.k].intArray_[ci.i][ci.j][whichInt]);
      }

      for (std::size_t whichFloat = 0; whichFloat < nFloatsPerCell_; ++whichFloat) {
        currFloat.push_back(lattice_[ci.k].floatArray_[ci.i][ci.j][whichFloat]);
      }

    }
    
  }
  
}
#endif

#if KMC_PARALLEL
void Lattice::Impl_::setLatticeValsFromRecvBuffers_(std::vector<std::vector<IJK> > & recvIndsBuffer,
                                                    const std::vector<std::vector<int> > & recvIntBuffer,
                                                    const std::vector<std::vector<double> > & recvFloatBuffer,
                                                    const std::vector<LatticePlanarBBox> & recvBufferBounds) {

  std::size_t recvIndsBufferSize = recvIndsBuffer.size();

  for (std::size_t bufType = 0; bufType < recvIndsBufferSize; ++bufType) {

    std::vector<IJK> & currInds = recvIndsBuffer[bufType];
    const std::vector<int> & currInt = recvIntBuffer[bufType];
    const std::vector<double> & currFloat = recvFloatBuffer[bufType];
    const LatticePlanarBBox & currBounds = recvBufferBounds[bufType];
    
    std::size_t currIndsSize = currInds.size();

    for (std::size_t i = 0; i < currIndsSize; ++i) {

      IJK & ci = currInds[i];

      if (ci.i < currBounds.imin) {
        ci.i += globalPlanarDims_[0];
      }
      else if (ci.i >= currBounds.imaxP1) {
        ci.i -= globalPlanarDims_[0];
      }

      if (ci.j < currBounds.jmin) {
        ci.j += globalPlanarDims_[1];
      }
      else if (ci.j >= currBounds.jmaxP1) {
        ci.j -= globalPlanarDims_[1];
      }

      std::size_t intOffset = nIntsPerCell_*i;
      std::size_t floatOffset = nFloatsPerCell_*i;

      for (std::size_t whichInt = 0; whichInt < nIntsPerCell_; ++whichInt) {
        lattice_[ci.k].intArray_[ci.i][ci.j][whichInt] = currInt[whichInt + intOffset];
      }

      for (std::size_t whichFloat = 0; whichFloat < nFloatsPerCell_; ++whichFloat) {
        lattice_[ci.k].floatArray_[ci.i][ci.j][whichFloat] = currFloat[whichFloat + floatOffset];
      }

    }

  }
}
#endif

#if KMC_PARALLEL
void Lattice::Impl_::reExportIfNeededActual_() {
  int recvIndsBufferSize = localRecvIndsBuffer_.size();

  for (int bufType = 0; bufType < recvIndsBufferSize; ++bufType) {
    std::vector<IJK> & currInds = localRecvIndsBuffer_[bufType];
    
    for (std::vector<IJK>::const_iterator itr = currInds.begin(),
           itrEnd = currInds.end(); itr != itrEnd; ++itr) {

      const IJK & ci = *itr;
      const Impl_::CellParInfo_ & cellParInfo_ = cellParInfoArray_[ci.i][ci.j];
      const Impl_::LocalBufferLocs_ & bufLocs = exportValToLocalBufferLocs_[cellParInfo_.exportVal];

      for (std::size_t i = 0; i < bufLocs.nBuffers; ++i) {
        if (bufLocs.bufType[i] != bufType) {
          localSendIndsBuffer_[bufLocs.sectNum[i]][bufLocs.bufType[i]].push_back(ci);
        }
      }

    }
  }
}
#endif

Lattice::Lattice(const LatticeParams & paramsForLattice)
  : pImpl_(new Impl_(paramsForLattice, this)) {
  trackChanges(TrackType::NONE);
  paramsForLattice.latInit(*this);

#if KMC_PARALLEL
  if (paramsForLattice.noAddingPlanesDuringSimulation) {
    // These should only be changed *after* initialization.
    pImpl_->setNewLatticeHeight_ = &Impl_::setNewLatticeHeightFake_;
    pImpl_->appendPlaneOnly_ = &Impl_::appendPlaneFake_;
  }
#endif
}

Lattice::~Lattice() {}

int Lattice::nProcs() const {return pImpl_->nProcs_;}
int Lattice::procID() const {return pImpl_->procID_;}

#if KMC_PARALLEL
const MPI_Comm & Lattice::comm() const {return pImpl_->latticeCommCart_;}
#endif

int Lattice::procPerDim(int dim) const {assert((dim >=0) && (dim < 2)); return pImpl_->procsPerDim_[dim];}
int Lattice::commCoord(int dim) const {assert((dim >=0) && (dim < 2)); return pImpl_->commCoords_[dim];}

int Lattice::ghostExtent(int dim) const {assert((dim >=0) && (dim < 2)); return pImpl_->ghostExtent_[dim];}

int Lattice::currHeight() const {return pImpl_->lattice_.size();}

int Lattice::nIntsPerCell() const {return pImpl_->nIntsPerCell_;}
int Lattice::nFloatsPerCell() const {return pImpl_->nFloatsPerCell_;}

int Lattice::numSectors() const {
#if KMC_PARALLEL
  return pImpl_->nSectors_;
#else
  return 1;
#endif
}

int Lattice::sectorOfIndices(const CellInds & ci) const {
#if KMC_PARALLEL
  return pImpl_->cellParInfoArray_[ci.i][ci.j].sectNum;
#else
  return 0;
#endif
}

#if KMC_PARALLEL
int Lattice::exportVal(const CellInds & ci) const {
  return pImpl_->cellParInfoArray_[ci.i][ci.j].exportVal;
}
#endif

int Lattice::getInt(const CellInds & ci, int whichInt) const {
  int iWrapped = ci.i;
  int jWrapped = ci.j;

#if KMC_PARALLEL
  KMC_CALL_MEMBER_FUNCTION(*pImpl_, pImpl_->wrapIndsIfNeeded_)(iWrapped, jWrapped);
#else
  pImpl_->wrapBothInds_(iWrapped, jWrapped);
#endif

  return pImpl_->lattice_[ci.k].intArray_[iWrapped][jWrapped][whichInt];
}

double Lattice::getFloat(const CellInds & ci, int whichFloat) const {
  int iWrapped = ci.i;
  int jWrapped = ci.j;

#if KMC_PARALLEL
  KMC_CALL_MEMBER_FUNCTION(*pImpl_, pImpl_->wrapIndsIfNeeded_)(iWrapped, jWrapped);
#else
  pImpl_->wrapBothInds_(iWrapped, jWrapped);
#endif

  return pImpl_->lattice_[ci.k].floatArray_[iWrapped][jWrapped][whichFloat];
}

void Lattice::setInt(const CellInds & ci, int whichInt, int val) {
  KMC_CALL_MEMBER_FUNCTION(*pImpl_, pImpl_->setInt_)(ci, whichInt, val);
}

void Lattice::setFloat(const CellInds & ci, int whichFloat, double val) {
  KMC_CALL_MEMBER_FUNCTION(*pImpl_, pImpl_->setFloat_)(ci, whichFloat, val);
}

void Lattice::addPlanes(int numPlanesToAdd) {
  pImpl_->addPlanes_(numPlanesToAdd);
}

void Lattice::reservePlanes(int numTotalPlanesToReserve) {
  assert(!(numTotalPlanesToReserve < 0));
  pImpl_->lattice_.reserve(numTotalPlanesToReserve);
}

int Lattice::planesReserved() const {
  return pImpl_->lattice_.capacity();
}

void Lattice::recvGhosts(int sectNum) {
#if KMC_PARALLEL
  int newLatticeHeight = KMC_CALL_MEMBER_FUNCTION(*pImpl_, pImpl_->setNewLatticeHeight_)(currHeight());
  KMC_CALL_MEMBER_FUNCTION(*pImpl_, pImpl_->recvGhosts_)(sectNum, newLatticeHeight);
#endif
}

void Lattice::recvGhostsUpdate(int sectNum) {
#if KMC_PARALLEL
  KMC_CALL_MEMBER_FUNCTION(*pImpl_, pImpl_->setNewLatticeHeight_)(currHeight());

  pImpl_->setSendIntFloatBuffers_(pImpl_->localSendIndsBuffer_[sectNum],
                                  pImpl_->localSendIntBuffer_[sectNum],
                                  pImpl_->localSendFloatBuffer_[sectNum]);

  KMC_CALL_MEMBER_FUNCTION(*pImpl_, pImpl_->recvGhostsUpdate_)(sectNum);

  pImpl_->setLatticeValsFromRecvBuffers_(pImpl_->ghostRecvIndsBuffer_,
                                         pImpl_->ghostRecvIntBuffer_,
                                         pImpl_->ghostRecvFloatBuffer_,
                                         pImpl_->ghostRecvBufferBounds_[sectNum]);
#endif
}

void Lattice::sendGhosts(int sectNum) {
#if KMC_PARALLEL
  int newLatticeHeight = KMC_CALL_MEMBER_FUNCTION(*pImpl_, pImpl_->setNewLatticeHeight_)(currHeight());
  KMC_CALL_MEMBER_FUNCTION(*pImpl_, pImpl_->sendGhosts_)(sectNum, newLatticeHeight);
#endif
}

void Lattice::sendGhostsUpdate(int sectNum) {
#if KMC_PARALLEL
  KMC_CALL_MEMBER_FUNCTION(*pImpl_, pImpl_->setNewLatticeHeight_)(currHeight());

  pImpl_->setSendIntFloatBuffers_(pImpl_->ghostSendIndsBuffer_,
                                  pImpl_->ghostSendIntBuffer_,
                                  pImpl_->ghostSendFloatBuffer_);

  KMC_CALL_MEMBER_FUNCTION(*pImpl_, pImpl_->sendGhostsUpdate_)(sectNum);

  pImpl_->setLatticeValsFromRecvBuffers_(pImpl_->localRecvIndsBuffer_,
                                         pImpl_->localRecvIntBuffer_,
                                         pImpl_->localRecvFloatBuffer_,
                                         pImpl_->localRecvBufferBounds_[sectNum]);

  KMC_CALL_MEMBER_FUNCTION(*pImpl_, pImpl_->reExportIfNeeded_)();
#endif
}

void Lattice::getLocalPlanarBBox(bool wGhost,
                                 LatticePlanarBBox & bbox) const {
#if KMC_PARALLEL
  pImpl_->getLocalPlanarBBox_(wGhost, bbox.imin, bbox.imaxP1, bbox.jmin, bbox.jmaxP1);
#else
  getGlobalPlanarBBox(bbox);
#endif
}

void Lattice::getSectorPlanarBBox(int sectNum,
                                  LatticePlanarBBox & bbox) const {
#if KMC_PARALLEL
  assert(sectNum >= 0);
  assert(sectNum < numSectors());

  bbox.imin = pImpl_->iminS_[sectNum];
  bbox.imaxP1 = pImpl_->imaxP1S_[sectNum];

  bbox.jmin = pImpl_->jminS_[sectNum];
  bbox.jmaxP1 = pImpl_->jmaxP1S_[sectNum];
#else
  getGlobalPlanarBBox(bbox);
#endif
}

void Lattice::getGlobalPlanarBBox(LatticePlanarBBox & bbox) const {

  bbox.imin = bbox.jmin = 0;
  
  bbox.imaxP1 = pImpl_->globalPlanarDims_[0];
  bbox.jmaxP1 = pImpl_->globalPlanarDims_[1];
}

void Lattice::trackChanges(TrackType::Type trackType) {

  switch (trackType) {
  case TrackType::CHECK_ONLY_IF_CHANGE_OCCURS:
    pImpl_->setInt_ = &Impl_::setIntAndRecordThatChangeOccurred_;
    pImpl_->setFloat_ = &Impl_::setFloatAndRecordThatChangeOccurred_;

    pImpl_->appendPlane_ = &Impl_::appendPlaneAndRecordThatChangeOccurred_;
    break;

  case TrackType::RECORD_CHANGED_CELL_INDS:
    pImpl_->setInt_ = &Impl_::setIntAndRecordChangedCellInds_;
    pImpl_->setFloat_ = &Impl_::setFloatAndRecordChangedCellInds_;

    pImpl_->appendPlane_ = &Impl_::appendPlaneAndRecordChangedCellInds_;
    break;

  case TrackType::RECORD_ONLY_OTHER_CHANGED_CELL_INDS:
    pImpl_->setInt_ = &Impl_::setIntOnly_;
    pImpl_->setFloat_ = &Impl_::setFloatOnly_;
    
    pImpl_->appendPlane_ = &Impl_::appendPlaneAndRecordChangedCellInds_;
    break;

  case TrackType::NONE:
    pImpl_->setInt_ = &Impl_::setIntOnly_;
    pImpl_->setFloat_ = &Impl_::setFloatOnly_;

    // Note that appendPlaneOnly_ is itself a function pointer.
    pImpl_->appendPlane_ = pImpl_->appendPlaneOnly_;
    break;

  default:
    // Should not be able to get here
    abortWithMsg("Bad track type value");
  }
  
  // Disregarding any calls to setInt or setFloat up to this point
  pImpl_->latticeModified_ = false;
  pImpl_->changedCellInds_.clear();
  pImpl_->otherCheckedCellInds_.clear();
}

bool Lattice::hasChanged() const {return pImpl_->latticeModified_;}

const Lattice::ChangedCellInds & Lattice::getChangedCellInds() const {
  return pImpl_->changedCellInds_;
}

const Lattice::OtherCheckedCellInds & Lattice::getOtherCheckedCellInds() const {
  return pImpl_->otherCheckedCellInds_;
}

void Lattice::wrapIndsIfNeeded(CellInds & ci) const {
#if KMC_PARALLEL
  KMC_CALL_MEMBER_FUNCTION(*pImpl_, pImpl_->wrapIndsIfNeeded_)(ci.i, ci.j);
#else
   pImpl_->wrapBothInds_(ci.i, ci.j);
#endif
}

int Lattice::wrapI(const CellInds & ci) const {
  return wrapInd(ci.i, pImpl_->globalPlanarDims_[0]);
}

int Lattice::wrapJ(const CellInds & ci) const {
  return wrapInd(ci.j, pImpl_->globalPlanarDims_[1]);
}

const std::vector<std::vector<IJK> > & Lattice::getReceivedGhostInds() const {
#if KMC_PARALLEL
  return pImpl_->ghostRecvIndsBuffer_;
#else
  return pImpl_->dummyNullVector_;
#endif
}

const std::vector<std::vector<IJK> > & Lattice::getReceivedLocalInds() const {
#if KMC_PARALLEL
  return pImpl_->localRecvIndsBuffer_;
#else
  return pImpl_->dummyNullVector_;
#endif
}

void Lattice::clearGhostsToSend() {
#if KMC_PARALLEL
  for (std::size_t i = 0; i < pImpl_->ghostSendIndsBuffer_.size(); ++i) {
    pImpl_->ghostSendIndsBuffer_[i].clear();
  }
#endif  
}

#if KMC_PARALLEL
bool Lattice::addToExportBufferIfNeeded(const CellInds & ci, int & sectNum) {
  const Impl_::CellParInfo_ & cellParInfo_ = pImpl_->cellParInfoArray_[ci.i][ci.j];
  sectNum = cellParInfo_.sectNum;

  assert(sectNum != -2);

  bool isGhost = (sectNum < 0);

  if (isGhost) {
    // Note: If isGhost == true, then cellParInfo_.exportVal > 0.
    pImpl_->ghostSendIndsBuffer_[pImpl_->exportValToGhostBufType_[cellParInfo_.exportVal]].push_back(ci);
  }
  else {
    const Impl_::LocalBufferLocs_ & bufLocs = pImpl_->exportValToLocalBufferLocs_[cellParInfo_.exportVal];

    for (std::size_t i = 0; i < bufLocs.nBuffers; ++i) {
      pImpl_->localSendIndsBuffer_[bufLocs.sectNum[i]][bufLocs.bufType[i]].push_back(ci);
    }
  }

  return isGhost;
}
#endif

bool Lattice::addToExportBufferIfNeeded(const CellInds & ci) {
  bool isGhost = false;

#if KMC_PARALLEL
  int sectNum;
  isGhost = addToExportBufferIfNeeded(ci, sectNum);
#endif

  return isGhost;
}

AddEmptyPlanes::AddEmptyPlanes(int numPlanesToAdd)
  : numPlanesToAdd_(numPlanesToAdd)
{}

void AddEmptyPlanes::operator()(Lattice & lattice) const {
  lattice.addPlanes(numPlanesToAdd_);
}
