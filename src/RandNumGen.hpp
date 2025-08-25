#ifndef RAND_NUM_GEN_HPP
#define RAND_NUM_GEN_HPP

#include <boost/shared_ptr.hpp>

// Note: This header file is documented via Doxygen
// <http://www.doxygen.org>. Comments for Doxygen begin with '/*!' or
// '//!', and descriptions of functions, class and member functions
// occur *before* their corresponding class declarations and function
// prototypes.

/*! \file
  \brief Defines abstract class RandNumGen
 */

namespace KMCThinFilm {
  
  /*! Abstract class to implement wrapper classes for random number
      generators */
  class RandNumGen {
  public:
    /*! Returns a (pseudo)random number between 0 and 1,
        <STRONG>NOT</STRONG> including either 0 or 1. 

	It is especially important for implementations of this class
	to exclude 0 from this function's return value. The amount of
	time by which a sector's clock is incremented is proportional
	to log(r), where r is the return value of the function, and of
	course, log(0) is undefined.

	It is also ill-adviced to allow 1 to be a return value of this
	function, either, since log(1) = 0, which would imply that an
	event in a KMC simulation took no time at all.
    */
    virtual double getNumInOpenIntervalFrom0To1() = 0;

    //! \cond HIDE_FROM_DOXYGEN
    virtual ~RandNumGen() {}
    //! \endcond

  };

  /*! "Smart" pointer used to point to an implementation of the RandNumGen class. */
  typedef boost::shared_ptr<RandNumGen> RandNumGenSharedPtr;

}

#endif /* RAND_NUM_GEN_HPP */
