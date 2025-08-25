#ifndef EVENT_EXECUTOR_GROUP_HPP
#define EVENT_EXECUTOR_GROUP_HPP

#include <vector>

#include <boost/scoped_ptr.hpp>

#include "CellInds.hpp"
#include "CellNeighOffsets.hpp"
#include "EventExecutor.hpp"

/*! \file
  \brief Defines the class EventExecutorGroup.
 */

namespace KMCThinFilm {

  /*! Collection of function objects satisfying the
      EventExecutorAutoTrack and EventExecutorSemiManualTrack
      signatures, used to specify a group of possible cell-centered
      events.
   */
  class EventExecutorGroup {
    friend class Simulation;
  public:

    /*! Constructs an EventExecutorGroup object. */
    explicit EventExecutorGroup(int numEventsInGroup /*!< Number of
                                                        function
                                                        objects in
                                                        group that
                                                        execute
                                                        events */);

    //! \cond HIDE_FROM_DOXYGEN
    EventExecutorGroup(const EventExecutorGroup & eeg);
    EventExecutorGroup & operator=(const EventExecutorGroup & rhs);
    //! \endcond

    /*! Adds an EventExecutorAutoTrack object to the group      
     */
    void addEventExecutor(int whichEvent /*!< Integer ID of event in
                                           group. Must be >= 0 and
                                           less than
                                           numEventExecutors()*/,
                          EventExecutorAutoTrack evExec /*!< Function object to execute the event */);

    /*! Adds an EventExecutorSemiManualTrack object and its associated
        offsets to the group.

        The argument <VAR>cnoVec</VAR> indicates the relative
        positions of cells that would be directly changed by executing
        <VAR>evExec</VAR>. If, for example, the event changed cells
        \f$(i_1, j_1, k_1)\f$ and \f$(i_1 + a, j_1 + b, k_1 + c)\f$,
        and also \f$(i_2, j_2, k_2)\f$, \f$(i_2 + d_1, j_2 + e_1, k_2
        + f_1)\f$, and \f$(i_2 + d_2, j_2 + e_2, k_2 + f_2)\f$, where
        the relative positions of \f$(i_1, j_1, k_1)\f$ and \f$(i_2,
        j_2, k_2)\f$ may vary with respect to each other, then the
        first entry of <VAR>cnoVec</VAR> would be the offsets
        associated with \f$(i_1, j_1, k_1)\f$, that is, \f$(0,0,0)\f$
        and \f$(a,b,c)\f$, and the second entry of <VAR>cnoVec</VAR>
        would be the offsets associated with \f$(i_2, j_2, k_2)\f$,
        that is, \f$(0,0,0)\f$, \f$(d_1, e_1, f_1)\f$ and \f$(d_2,
        e_2, f_2)\f$.

    */
    void addEventExecutor(int whichEvent /*!< Integer ID of event in
                                           group. Must be >= 0 and
                                           less than
                                           numEventExecutors()*/,
                          EventExecutorSemiManualTrack evExec /*!< Function object to execute the event */,
                          const std::vector<CellNeighOffsets> & cnoVec  /*!<
                                                                          Relative positions of cells 
                                                                          directly changed by event. */);

    /*! Removes all previous function objects from the group and
        allocates memory for <VAR>numEventsInGroup</VAR> function
        objects. */
    void resetGroup(int numEventsInGroup);

    /*! Removes all previous function objects from the group. */
    void clearGroup();

    /*! Returns the number of function objects stored in the group. */
    int numEventExecutors() const;

    //! \cond HIDE_FROM_DOXYGEN
    ~EventExecutorGroup();
    //! \endcond

  private:
    
    struct EventExecEnum {
      enum Type {AUTO, SEMIMANUAL};
    };
    
    // To be used by Simulation class
    EventExecEnum::Type getEventExecutorType(int whichEvent) const;
    
    EventExecutorAutoTrack getEventExecutorAutoTrack(int whichEvent) const;
    EventExecutorSemiManualTrack getEventExecutorSemiManualTrack(int whichEvent) const;
    const std::vector<CellNeighOffsets> & getEventExecutorOffsetsVec(int whichEvent) const;

    class Impl_;
    boost::scoped_ptr<Impl_> pImpl_;
  };

}

#endif /* EVENT_EXECUTOR_GROUP_HPP */
