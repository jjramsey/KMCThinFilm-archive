#include "Solver.hpp"
#include "ErrorHandling.hpp"
#include "CallMemberFunction.hpp"

#include <limits>

using namespace KMCThinFilm;

#if KMC_PARALLEL

void Solver::getTstopFromLocalPropensities_(double propensityMaxLocal,
					    const MPI_Comm & comm,
					    double & tStop) const {
  double propensityMax;
  MPI_Allreduce(&propensityMaxLocal, &propensityMax, 1, MPI_DOUBLE, MPI_MAX, comm);

  if (propensityMax > 0) {
    double tStopPoss = nStop_/propensityMax;
    tStop = ((tStopPoss > tStopMax_) ? tStopMax_ : tStopPoss);
  }
}

void Solver::updateTStopMaxAvgPropensityPerPossEvent_(const MPI_Comm & comm,
						      double & tStop) {
  getTstopFromLocalPropensities_(getLocalMaxAvgPropensityPerPossEvent(), comm, tStop);
}

void Solver::updateTStopMaxSinglePropensity_(const MPI_Comm & comm,
					     double & tStop) {
  getTstopFromLocalPropensities_(getLocalMaxSinglePropensity(), comm, tStop);
}

void Solver::updateTStopFixedValue_(const MPI_Comm & comm,
				    double & tStop) {
  tStop = tStopFixed_;
}

bool Solver::timeIncrSchemeIsAdaptive() const {
  return timeIncrSchemeIsAdaptive_;
}

void Solver::updateTStop(const MPI_Comm & comm, double & tStop) {
  KMC_CALL_MEMBER_FUNCTION(*this, updateTStop_)(comm, tStop);
}

void Solver::setTimeIncrScheme(const TimeIncr::SchemeVars & vars) {
  tIncrSchemeName_ = vars.getSchemeName();

  tStopMax_ = vars.getSchemeParamOrReturnDefaultVal(TimeIncr::SchemeParam::TSTOP_MAX,
						   std::numeric_limits<double>::max());

  switch (tIncrSchemeName_) {
  case TimeIncr::SchemeName::MAX_AVG_PROPENSITY_PER_POSS_EVENT:
    {
      nStop_ = vars.getSchemeParamOrDie(TimeIncr::SchemeParam::NSTOP,
					"Parameter NSTOP not set for MAX_AVG_PROPENSITY_PER_POSS_EVENT time increment scheme");

      updateTStop_ = &Solver::updateTStopMaxAvgPropensityPerPossEvent_;
      timeIncrSchemeIsAdaptive_ = true;
    }
    break;
  case TimeIncr::SchemeName::MAX_SINGLE_PROPENSITY:
    {
      nStop_ = vars.getSchemeParamOrDie(TimeIncr::SchemeParam::NSTOP,
					"Parameter NSTOP not set for MAX_SINGLE_PROPENSITY time increment scheme");

      updateTStop_ = &Solver::updateTStopMaxSinglePropensity_;
      timeIncrSchemeIsAdaptive_ = true;
    }
    break;
  case TimeIncr::SchemeName::FIXED_VALUE:
    {
      tStopFixed_ = vars.getSchemeParamOrDie(TimeIncr::SchemeParam::TSTOP,
					     "Parameter TSTOP not set for FIXED_VALUE time increment scheme");

      updateTStop_ = &Solver::updateTStopFixedValue_;
      timeIncrSchemeIsAdaptive_ = false;
    }
    break;
  default:
    exitWithMsg("Bad parallel time increment scheme name");
  }
}
#endif
