#ifndef RAND_NUM_GEN_DCMT_HPP
#define RAND_NUM_GEN_DCMT_HPP

extern "C" {
#include <dc.h>
}

#include "RandNumGen.hpp"

/*!\file
  \brief Defines the RandNumGenDCMT class
 */

namespace KMCThinFilm {

  /*! Wrapper class for the parallel pseudorandom number generator
      library DCMT: "Dynamic Creator of Mersenne Twisters."

      DCMT is a parallel pseudorandom number generator by Makoto
      Matsumoto and Takuji Nishimura, and it is based on the serial
      Mersenne Twister pseudorandom number generator. The library may
      be obtained from this web page:
      http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/DC/dc.html
   */
  class RandNumGenDCMT : public RandNumGen {
  public:

    /*! This enumeration forces the <VAR>period</VAR> argument of the
      RandNumGenDCMT constructor to be one of a certain number of
      discrete values.

      When the value of this enumeration is P<EM>n</EM>, where
      <EM>n</EM> is a number, then the period of the pseudorandom
      number stream is \f$2^n - 1\f$. The value of
      <EM>n</EM> is such that the period is a Mersenne prime.
    */
    enum Period {
      P521 = 521,
      P607 = 607,
      P1279 = 1279, 
      P2203 = 2203, 
      P2281 = 2281, 
      P3217 = 3217, 
      P4253 = 4253, 
      P4423 = 4423, 
      P9689 = 9689, 
      P9941 = 9941, 
      P11213 = 11213, 
      P19937 = 19937, 
      P21701 = 21701, 
      P23209 = 23209, 
      P44497 = 44497
    };

    /*! Constructor, sets the initial state for each of the parallel
        streams of the pseudorandom number generator.
     */
    RandNumGenDCMT(int rank /*!< MPI rank of the processor */,
		   uint32_t seed /*!< Global seed for the pseudorandom
                                    number generator */, 
		   uint32_t seed_perProc /*!< Per-processor seed for
                                            the pseudorandom number
                                            generator, should be
                                            different for each process */,
		   Period period = P521 /*!< Enumeration constant
                                           indicating the exponent for
                                           the period of the
                                           pseudorandom number
                                           stream */);

    /*! Returns a (pseudo)random number between 0 and 1, not including either 0 or 1. */
    virtual double getNumInOpenIntervalFrom0To1();

    virtual ~RandNumGenDCMT();
  private:
    uint32_t seed_;
    mt_struct * mts_;
  };

}

#endif /* RAND_NUM_GEN_DCMT_HPP */
