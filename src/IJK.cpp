#include "IJK.hpp"

#include <boost/functional/hash.hpp>

namespace KMCThinFilm {

  std::size_t hash_value(const IJK & ijk) {
    std::size_t seed = 0;
    boost::hash_combine(seed, ijk.i);
    boost::hash_combine(seed, ijk.j);
    boost::hash_combine(seed, ijk.k);
    return seed;
  }

}
