#ifndef ERROR_HANDLING_FUNCS_HPP
#define ERROR_HANDLING_FUNCS_HPP

#include <string>

// Note: This header file is documented via Doxygen
// <http://www.doxygen.org>. Comments for Doxygen begin with '/*!' or
// '//!', and descriptions of functions, class and member functions
// occur *before* their corresponding class declarations and function
// prototypes.

/*! \file
  \brief Defines functions that exit the KMC application in case of an error.

  These are free functions that don't belong to any particular class.
*/

namespace KMCThinFilm {

  /*! Exits the program and prints the message string <VAR>msg</VAR>
    on rank 0 of MPI_COMM_WORLD. */
  void exitWithMsg(const std::string & msg);

  /*! Aborts the program and prints the message string <VAR>msg</VAR>
    on all processors. */
  void abortWithMsg(const std::string & msg);
  
  /*! Exits the program and prints the message string <VAR>msg</VAR>
      on rank 0 of MPI_COMM_WORLD if <VAR>condition</VAR> is true.

    This function allows for limited error handling, but only for the
    cases where <VAR>condition</VAR> is tested on all processors and
    <EM>will evaluate to the same value on all processors</EM>. If
    these conditions will not always be satisfied, then another method
    of error handling, e.g. abortOnCondition(), should be
    used instead.

    \see abortOnCondition()
   */
  void exitOnCondition(bool condition, const std::string & msg);

  /*! Aborts the program and prints the message string <VAR>msg</VAR>
      on all processors where <VAR>condition</VAR> is true.

    This function terminates less gracefully than exitOnCondition(),
    but may be used for cases where <VAR>condition</VAR> may not
    necessarily evaluate to the same value on all processors.

    \see exitOnCondition()
   */
  void abortOnCondition(bool condition, const std::string & msg);

}

#endif /* ERROR_HANDLING_FUNCS_HPP */
