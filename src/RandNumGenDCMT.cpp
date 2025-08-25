#include <climits>
#include <cfloat>

#include "RandNumGenDCMT.hpp"

using namespace KMCThinFilm;

RandNumGenDCMT::RandNumGenDCMT(int rank, 
			       uint32_t seed_master,
			       uint32_t seed_perProc,
			       Period period)
  : seed_(seed_master) {

  mts_ = get_mt_parameter_id_st(32, period, rank, seed_);
  sgenrand_mt(seed_perProc, mts_);  
}

double RandNumGenDCMT::getNumInOpenIntervalFrom0To1() {

  // If the denominator is UINT_MAX, the first argument of
  // get_mt_parameter_id_st *must* be 32.
  double myRand = static_cast<double>(genrand_mt(mts_))/UINT_MAX;

  // myRand is supposed to be in the *open* interval between 0 and 1;
  if (myRand == 0) {
    myRand += DBL_EPSILON;
  }
  else if (myRand == 1) {
    myRand -= DBL_EPSILON;
  }

  return myRand;
}

RandNumGenDCMT::~RandNumGenDCMT() {
  free_mt_struct(mts_);
}
