#ifndef IJK_HPP
#define IJK_HPP

#include <cstddef>

// Note: This header file is documented via Doxygen
// <http://www.doxygen.org>. Comments for Doxygen begin with '/*!' or
// '//!', and descriptions of functions, class and member functions
// occur *before* their corresponding class declarations and function
// prototypes.

namespace KMCThinFilm {

  /*! POD struct inherited by CellInds and CellIndsOffset. */
  struct IJK {
    // Note: In order to remain a POD, and thus readily transmissible
    // via MPI, this struct must *never* have user-declared
    // constructors, private or protected non-static data members,
    // base classes, or virtual functions.
    int i, j, k;

    /*! Defines a less-than operation via lexicographic comparison for
        use with ordered maps and sets.

      Given two IJK instances ijk1 and ijk2, ijk1 will be less than
      ijk2 if one of the following conditions is met:

      - ijk1.i < ijk2.i
      - ijk1.i equals ijk2.i, and ijk1.j < ijk2.j
      - ijk1.i equals ijk2.i, ijk1.j equals ijk2.j, and ijk1.k < ijk2.k

      This operator allows an IJK instance to be used as an index in a
      std::map (or element in a std::set). */
    bool operator<(const IJK & rhs) const {
      return (i < rhs.i)
        || ((i == rhs.i) && (j < rhs.j))
        || ((i == rhs.i) && (j == rhs.j) && (k < rhs.k));
    }

    /*! Defines an equality operation for use with unordered maps and sets.

      Two IJK instances ijk1 and ijk2 are equal if their corresponding
      components are equal, i.e. ijk1.i == ijk2.i, ijk1.j == ijk2.j,
      and ijk1.k == ijk2.k.

      This operator, along with the hash_value() function, allows an
      IJK instance to be used as an index in a boost::unordered_map
      (or element in a boost::unordered_set). */
    bool operator==(const IJK & rhs) const {
      return (i == rhs.i) && (j == rhs.j) && (k == rhs.k);
    }

  };

  /*! Hashing function to allow instances of SiteInds to be used as an
      index in a boost::unordered_map. */
  std::size_t hash_value(const IJK & ijk);

}

#endif /* IJK_HPP */
