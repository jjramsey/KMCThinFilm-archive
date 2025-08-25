#ifndef IDS_OF_SOLVERS_HPP
#define IDS_OF_SOLVERS_HPP

// This file is named IdsOfSolvers.hpp rather than SolverIds.hpp
// because the Doxygen documentation excludes files with the pattern
// "Solver*.hpp".

/*!\file
  \brief Defines the enumeration SolverId::Type
 */

namespace KMCThinFilm {

  /*! Namespace to enclose the SolverId::Type enumeration */
  namespace SolverId {
    
    /*! Types of solvers used in a KMC simulation */
    enum Type {
      DYNAMIC_SCHULZE /*!< A solver implementing an algorithm similar
                         to that described by Schulze in <EM>Physical
                         Review E</EM>, vol. 65, 036704 (2002), but
                         allowing for the number of event propensities
                         to be determined at runtime. This should
                         allow the choosing of an event to scale with
                         the number of unique propensities, rather
                         than the number of events. */,

      BINARY_TREE /*!< A solver where the possible events are stored
                     as leaves in a binary tree, allowing the choosing
                     of an event to scale as \f$O(log_2 N)\f$, where
                     <VAR>N</VAR> is the number of possible
                     events. May be faster than DYNAMIC_SCHULZE in
                     practice. */
    };

  }

}

#endif /* IDS_OF_SOLVERS_HPP */
