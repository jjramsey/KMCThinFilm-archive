#include "RandNumGenMT19937.hpp"

#include <cfloat>

using namespace KMCThinFilm;

double RandNumGenMT19937::getNumInOpenIntervalFrom0To1() {
  double myRand = uni01_(rng_);

  // myRand is supposed to be in the *open* interval between 0 and 1;
  if (myRand == 0) {
    myRand += DBL_EPSILON;
  }

  return myRand;
}
