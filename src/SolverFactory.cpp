#include "SolverFactory.hpp"

#include "Lattice.hpp"
#include "ErrorHandling.hpp"

#include "SolverDynamicSchulze.hpp"
#include "SolverBinaryTree.hpp"

namespace KMCThinFilm {

  // Right now, just using a "dumb" switch-based factory, since I
  // doubt that I'll have enough new solver types to justify something
  // fancier, like what the "Gang of Four" would do in a design
  // pattern.
  
  Solver * mkSolver(SolverId::Type sId, const Lattice * lattice) {
    
    Solver * solver = NULL;

    switch (sId) {
    case SolverId::DYNAMIC_SCHULZE:
      solver = new SolverDynamicSchulze(lattice);
      break;
    case SolverId::BINARY_TREE:
      solver = new SolverBinaryTree(lattice);
      break;
    default:
      exitWithMsg("Bad SolverId value");
    }

    return solver;
  }
  
}
