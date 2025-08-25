#ifndef SOLVER_FACTORY_HPP
#define SOLVER_FACTORY_HPP

#include "IdsOfSolvers.hpp"
#include "Solver.hpp"

namespace KMCThinFilm {

  class Lattice;

  // This returns a raw pointer to a Solver, which can be captured in
  // the constructor of an appropriate "smart" pointer.
  Solver * mkSolver(SolverId::Type sId, const Lattice * lattice);
  
}

#endif /* SOLVER_FACTORY_HPP */
