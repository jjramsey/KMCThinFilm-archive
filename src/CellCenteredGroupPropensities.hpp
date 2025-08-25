#ifndef PROPENSITY_GROUP_CALC_HPP
#define PROPENSITY_GROUP_CALC_HPP

#include <vector>
#include <boost/function.hpp>

#include "CellNeighProbe.hpp"

/*! \file
    \brief Defines the signature of the function object used to determine the propensities of a group of related events
*/

namespace KMCThinFilm {
  /*! Signature of the function object used to determine the propensities of a group of related events. */
  typedef boost::function<void (const CellNeighProbe &, std::vector<double> &)> CellCenteredGroupPropensities;
}

#endif /* PROPENSITY_GROUP_CALC_HPP */
