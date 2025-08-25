#include "KMC_Config.hpp"
#include "ErrorHandling.hpp"

#if KMC_PARALLEL
#include <mpi.h>
#endif

#include <cstdlib>
#include <iostream>

namespace KMCThinFilm {

  void exitWithMsg(const std::string & msg) {

#if KMC_PARALLEL
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) {
#endif
      std::cerr << msg << "\n";

#if KMC_PARALLEL
    }
    
    MPI_Finalize();
#endif

    std::exit(EXIT_FAILURE);
  }
  
  void exitOnCondition(bool condition, const std::string & msg) {
    
    if (condition) {
      exitWithMsg(msg);
    }
  }

  void abortWithMsg(const std::string & msg) {
    std::cerr << msg << "\n";

#if KMC_PARALLEL
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
#else
    std::exit(EXIT_FAILURE);
#endif
  }

  void abortOnCondition(bool condition, const std::string & msg) {

    if (condition) {
      abortWithMsg(msg);
    }
  }

}

