#ifndef INIT_LATTICE_HPP
#define INIT_LATTICE_HPP

#include <string>

#include <KMCThinFilm/Lattice.hpp>
#include <KMCThinFilm/ErrorHandling.hpp>

class InitLatticeFromFile {
public:
  InitLatticeFromFile(const std::string & inpFName)
    : inpFName_(inpFName)
  {}

  void operator()(KMCThinFilm::Lattice & lattice) const;
private:
  std::string inpFName_;
};

#endif /* INIT_LATTICE_HPP */
