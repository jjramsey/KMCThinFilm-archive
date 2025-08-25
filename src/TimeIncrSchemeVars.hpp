#ifndef TIME_INCR_SCHEMES_HPP
#define TIME_INCR_SCHEMES_HPP

#include <string>
#include <boost/scoped_ptr.hpp>

/*! \file
  \brief Defines the TimeIncr::SchemeVars class and related enumerations
 */

namespace KMCThinFilm {

  /*! This namespace wraps a class and enumerations relating to
      specifying parallel time-stepping schemes.*/
  namespace TimeIncr {

    /*! Namespace to enclose the SchemeName::Type enumeration */
    namespace SchemeName {

      /*! Types of parallel time-incrementing schemes */
      enum Type {
	MAX_AVG_PROPENSITY_PER_POSS_EVENT /*!< Scheme where the time
					    step is inversely
					    proportional to the
					    maximum of the average
					    propensities per possible
					    event from all sectors,
					    where the average
					    propensity per possible
					    event is the sum of all
					    propensities in a sector
					    divided by the number of
					    possible events in that
					    sector. This is similar to
					    the adaptive time step
					    scheme used in SPPARKS
					    (http://spparks.sandia.gov/). */,
	MAX_SINGLE_PROPENSITY /*!< Scheme where the time step is
                                 inversely proportional to largest
                                 event propensity available in the
                                 simulation. This is similar to the
                                 adaptive time step scheme suggested
                                 by Shim and Amar in <EM>Physical
                                 Review B</EM>, vol. 71, 125432
                                 (2005).*/,
	FIXED_VALUE /*!< Scheme where the time step is a fixed value. */
	//! \cond HIDE_FROM_DOXYGEN
	, BAD_VALUE
	//! \endcond
      };
    }

    /*! Namespace to enclose the SchemeParam::Type enumeration */
    namespace SchemeParam {

      /*! Parameters used in the parallel time-incrementing schemes */
      enum Type {
	NSTOP /*!< This is a premultiplier used in adaptive time
                 schemes where the time step takes the form
                 \f$t_{stop} = NSTOP/P\f$, where <VAR>P</VAR> is some
                 function of the available event propensities. */,
	TSTOP_MAX /*!< Maximum time step for the adaptive schemes. The
                      time step is set to this value if the step
                      determined by the adaptive scheme would
		      otherwise exceed it. */,
	TSTOP /*!< Time step used for scheme SchemeName::FIXED_VALUE. */
      };
    }

    /*! Type of object used to store parameters for the parallel time stepping schemes.

      Note that in a typical use of a class instance of this, the
      set*() functions will tend to be used rather than the get*()
      functions. However, the get*() functions are exposed for the
      sake of testing and debugging.
     */
    class SchemeVars {
    public:

      //! \cond HIDE_FROM_DOXYGEN
      SchemeVars();

      explicit SchemeVars(const SchemeVars & params);
      SchemeVars & operator=(const SchemeVars & rhs);
      //! \endcond

      /*! Sets the name of the time-stepping scheme to be used. */
      void setSchemeName(SchemeName::Type name);

      /*! Sets a parameter of the time-stepping scheme to be used. */
      void setSchemeParam(SchemeParam::Type paramName /*!< Name of parameter */,
			  double paramVal /*!< Value of parameter */);

      /*! Returns the name of the time-stepping scheme set by setSchemeName(). */
      SchemeName::Type getSchemeName() const;

      /*! Returns the value of the parameter, if it has been set.

	It is not an error if the parameter was not set, but if
	it wasn't, then the return value is likely garbage.
       */
      double getSchemeParamIfAvailable(SchemeParam::Type paramName /*!< Name of parameter */,
				       bool & isAvailable /*!< If this
                                                             is false,
                                                             then the
                                                             parameter
                                                             was not
                                                             set.*/) const;

      /*! Returns the value of the parameter or causes the KMC
          application to terminate with an error message. */
      double getSchemeParamOrDie(SchemeParam::Type paramName /*!< Name of parameter */,
				 const std::string & msgIfDie /*!<
                                                                 Error
                                                                 message
                                                                 printed
                                                                 if
                                                                 the
                                                                 parameter
                                                                 was not set. */) const;

      /*! Returns the value of the parameter if it has been set, and
          returns a default value otherwise. */
      double getSchemeParamOrReturnDefaultVal(SchemeParam::Type paramName /*!< Name of parameter */,
					      double defaultVal /*!
                                                                    The
                                                                    value
                                                                    to
                                                                    be
                                                                    returned
                                                                    if
                                                                    the
                                                                    parameter
                                                                    was
                                                                    not
                                                                    set.*/) const;

      //! \cond HIDE_FROM_DOXYGEN
      ~SchemeVars();
      //! \endcond

    private:
      class Impl_;
      boost::scoped_ptr<Impl_> pImpl_;
    };
    
  }

}

#endif /* TIME_INCR_SCHEMES_HPP */
