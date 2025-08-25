#include "wrapInd.hpp"

namespace KMCThinFilm {
  int wrapInd(int i, int dim) {
    int iout = i;

    while (iout < 0) { 
      iout += dim;
    }

    /* Note that I need the greater than _or equals to_ to take into
       account that i is supposed to be within 0 and (dim - 1). */
    while (iout >= dim) {
      iout -= dim;
    }
  
    return iout;
  }
}
