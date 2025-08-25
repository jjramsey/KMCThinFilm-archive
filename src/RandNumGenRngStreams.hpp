#ifndef RAND_NUM_GEN_RNG_STREAMS_HPP
#define RAND_NUM_GEN_RNG_STREAMS_HPP

#include <vector>
#include <string>

extern "C" {
#include <RngStream.h>
}

/*!\file
  \brief Defines the RandNumGenRngStreams class
 */

#include "RandNumGen.hpp"

namespace KMCThinFilm {

  /*! Wrapper class for the parallel pseudorandom number generator
      library RngStreams.

      RngStreams is a parallel pseudorandom number generator by Pierre
      L'Ecuyer and Richard Simard. It may be obtained from this web
      page: http://statmath.wu.ac.at/software/RngStreams/
  */
  class RandNumGenRngStreams : public RandNumGen {
  public:

    /*! Constructor, sets the initial state for each of the parallel
        streams of the pseudorandom number generator.

	The vector <VAR>seed</VAR> that seeds the pseudorandom number
	generator must have a length of 6, and it should be the same
	for every process.

	Each parallel pseudorandom number stream is offset by \f$rank
	\times (2^{advExp} + advConst)\f$, if <VAR>advExp</VAR> is
	nonnegative, or \f$rank \times (-2^{-advExp} + advConst)\f$,
	otherwise. The argument <VAR>rank</VAR> should be the MPI rank
	of the processor.
     */
    RandNumGenRngStreams(std::vector<unsigned long> & seed,
			 int rank,
			 long advExp = 127,
			 long advConst = 0);

    /*! Returns a (pseudo)random number between 0 and 1, not including either 0 or 1. */
    virtual double getNumInOpenIntervalFrom0To1();

    virtual ~RandNumGenRngStreams();
  private:
    RngStream rng_;
  };

}

#endif /* RAND_NUM_GEN_RNG_STREAMS_HPP */
