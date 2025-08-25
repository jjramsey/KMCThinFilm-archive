#ifndef RAND_NUM_GEN_MT19937_HPP
#define RAND_NUM_GEN_MT19937_HPP

/*!\file
  \brief Defines the RandNumGenMT19937 class
*/

#include "RandNumGen.hpp"

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>

namespace KMCThinFilm {
  
  /*! Wrapper class for the <EM>serial</EM> Mersenne Twister
      pseudorandom number generator with a period of \f$2^{19937} - 1\f$.

      This is not recommended for parallel simulations.
   */
  class RandNumGenMT19937 : public RandNumGen {
  public:
    /*! Constructor, sets the initial state for the pseudorandom
        number generator. */
    RandNumGenMT19937(unsigned int seed)
      : rng_(seed)
    {}

    /*! Returns a (pseudo)random number between 0 and 1, not including
        either 0 or 1. */
    virtual double getNumInOpenIntervalFrom0To1();

  private:
    boost::random::mt19937 rng_;
    boost::random::uniform_01<double> uni01_;
  };

}

#endif /* RAND_NUM_GEN_MT19937_HPP */
