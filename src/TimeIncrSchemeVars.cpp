#include "TimeIncrSchemeVars.hpp"
#include "ErrorHandling.hpp"

#include <map>

using namespace KMCThinFilm::TimeIncr;

struct SchemeVars::Impl_ {

  // Making sure the name_ is initialized to a known value
  Impl_()
    : name_(SchemeName::BAD_VALUE)
  {}

  SchemeName::Type name_;

  typedef std::map<SchemeParam::Type, double> ParamMap_;
  ParamMap_ paramMap_;

  bool get_(SchemeParam::Type paramName, double & paramVal) const {
    ParamMap_::const_iterator itr = paramMap_.find(paramName);

    if (itr != paramMap_.end()) {
      paramVal = itr->second;
      return true;
    }
    else {
      return false;
    }

  }

};

SchemeVars::SchemeVars()
  : pImpl_(new Impl_)
{}

SchemeVars::SchemeVars(const SchemeVars & params)
  : pImpl_(new Impl_(*(params.pImpl_)))
{}

SchemeVars & SchemeVars::operator=(const SchemeVars & rhs) {

  if (this != &rhs) {
    *pImpl_ = *(rhs.pImpl_);
  }
  
  return *this;
}

SchemeVars::~SchemeVars() {}

void SchemeVars::setSchemeName(SchemeName::Type name) {
  pImpl_->name_ = name;
}

void SchemeVars::setSchemeParam(SchemeParam::Type paramName, double paramVal) {
  pImpl_->paramMap_[paramName] = paramVal;
}

SchemeName::Type SchemeVars::getSchemeName() const {
  return pImpl_->name_;
}

double SchemeVars::getSchemeParamIfAvailable(SchemeParam::Type paramName,
					     bool & isAvailable) const {
  double retVal = -1.0; // Initializing to some default value. Doesn't
			// matter what the value is.
  isAvailable = pImpl_->get_(paramName, retVal);
  return retVal;
}

double SchemeVars::getSchemeParamOrDie(SchemeParam::Type paramName,
				       const std::string & msgIfDie) const {  
  double retVal;
  KMCThinFilm::exitOnCondition(!(pImpl_->get_(paramName, retVal)), msgIfDie);
  return retVal;
}

double SchemeVars::getSchemeParamOrReturnDefaultVal(SchemeParam::Type paramName,
						    double defaultVal) const {
  double retVal;
  if (pImpl_->get_(paramName, retVal)) {
    return retVal;
  }
  else {
    return defaultVal;
  }  
}
