#include <cfloat>

#include "RandNumGenRngStreams.hpp"
#include "ErrorHandling.hpp"

using namespace KMCThinFilm;

RandNumGenRngStreams::RandNumGenRngStreams(std::vector<unsigned long> & seed,
					   int rank,
					   long advExp, long advConst) {

  exitOnCondition(seed.size() < 6, "The \"seed\" argument must be of length 6");

  RngStream_SetPackageSeed(&seed[0]);

  rng_ = RngStream_CreateStream("KMCThinFilm");

  for (int i = 0; i < rank; ++i) {
    RngStream_AdvanceState(rng_, advExp, advConst);
  }

}

double RandNumGenRngStreams::getNumInOpenIntervalFrom0To1() {
  double myRand = RngStream_RandU01(rng_);

  // myRand is supposed to be in the *open* interval between 0 and 1;
  if (myRand == 0) {
    myRand += DBL_EPSILON;
  }
  else if (myRand == 1) {
    myRand -= DBL_EPSILON;
  }

  return myRand;
}

RandNumGenRngStreams::~RandNumGenRngStreams() {
  RngStream_DeleteStream(&rng_);
}
