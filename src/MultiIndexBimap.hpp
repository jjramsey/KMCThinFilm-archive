#ifndef MULTI_INDEX_BIMAP_HPP
#define MULTI_INDEX_BIMAP_HPP

#include <stdexcept>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/hashed_index.hpp>

namespace KMCThinFilm {

  /* This is an implementation of a bimap using Boost.MultiIndex
     rather than Boost.Bimap, which produces (apparently benign)
     warnings in Intel compilers and also causes problems with the TAU
     profiler for some reason.

     Used only if KMC_AVOID_BOOST_BIMAP is defined.
  */

  template<typename LeftType,typename RightType>
  class MultiIndexBimap
  {
    struct LeftMemberSpace_ {
      LeftMemberSpace_(MultiIndexBimap * bimap)
	: bimap_(bimap)
      {}

      const RightType & at(const LeftType & lt) const {return bimap_->left_at_(lt);}

      void erase(const LeftType & lt) {bimap_->left_erase_(lt);}

      MultiIndexBimap * bimap_;
    };

    struct RightMemberSpace_ {
      //friend class MultiIndexBimap;

      RightMemberSpace_(MultiIndexBimap * bimap)
	: bimap_(bimap)
      {}

      const LeftType & at(const RightType & rt) const {return bimap_->right_at_(rt);}

      MultiIndexBimap * bimap_;
    };

  public:
    MultiIndexBimap()
      : left(this), right(this)
    {}

    struct value_type {
      value_type(const LeftType & lt, const RightType & rt)
	: left(lt), right(rt)
      {}
    
      LeftType left;
      RightType right;
    };

    void insert(const value_type & vt) {
      bimap_.insert(vt);
    }

    LeftMemberSpace_ left;
    RightMemberSpace_ right;

  private:

    // Tags for multi-index
    struct Left_{};
    struct Right_{};

    typedef boost::multi_index_container<
      value_type,
      boost::multi_index::indexed_by<
	boost::multi_index::hashed_unique<
	  boost::multi_index::tag<Left_>,
	  boost::multi_index::member<value_type,LeftType,&value_type::left> >,
	    boost::multi_index::hashed_unique<
	      boost::multi_index::tag<Right_>,
	      boost::multi_index::member<value_type,RightType,&value_type::right> >
	>
      > MultiIndexBimap_;

    MultiIndexBimap_ bimap_;

    typedef typename MultiIndexBimap_::template index<Left_>::type::iterator MultiIndexBimap_left_iterator;
    typedef typename MultiIndexBimap_::template index<Right_>::type::iterator MultiIndexBimap_right_iterator;

    const RightType & left_at_(const LeftType & lt) const {
      MultiIndexBimap_left_iterator itr = bimap_.template get<Left_>().find(lt);

      if (itr == bimap_.template get<Left_>().end()) {
	throw std::out_of_range("Left-hand key out of range");
      }

      return itr->right;
    }

    const LeftType & right_at_(const RightType & rt) const {
      MultiIndexBimap_right_iterator itr = bimap_.template get<Right_>().find(rt);

      if (itr == bimap_.template get<Right_>().end()) {
	throw std::out_of_range("Right-hand key out of range");
      }

      return itr->left;
    }

    void left_erase_(const LeftType & lt) {
      bimap_.template get<Left_>().erase(lt);
    }

  public:
    typedef typename MultiIndexBimap_::const_iterator const_iterator;

    const_iterator begin() const {return bimap_.begin();}
    const_iterator end() const {return bimap_.end();}
  };

}

#endif /* MULTI_INDEX_BIMAP_HPP */
