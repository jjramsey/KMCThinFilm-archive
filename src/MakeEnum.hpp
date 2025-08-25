#ifndef MAKE_ENUM_HPP
#define MAKE_ENUM_HPP

/*! \file
  \brief Defines convenience macros to set up enumeration constants for a KMC simulation.
 */

/*! Convenience macro for defining enumeration constants for the
  integers defined in each lattice cell.

  For example, the following

  \code
  KMC_MAKE_LATTICE_INTVAL_ENUM(GaInN_, Ga, In, N);
  \endcode

  creates an enumeration with values GaInN_IntVal::Ga,
  GaInN_IntVal::In, and GaInN_IntVal::N, and also defines
  GaInN_IntVal::SIZE, which indicates the number of enumeration
  constants defined.

  Note that this macro only works with preprocessors that do C99-style
  variadic macros.  */
#define KMC_MAKE_LATTICE_INTVAL_ENUM(EnumName, ...)	\
  struct EnumName##IntVal {				\
    enum Type {__VA_ARGS__, SIZE};			\
  }

/*! Convenience macro for defining enumeration constants for the
  floating-point parameters defined in each lattice cell.

  For example, the following

  \code
  KMC_MAKE_LATTICE_FLOATVAL_ENUM(SiCoords, SiX, SiY, SiZ);
  \endcode

  creates an enumeration with values SiCoordsFloatVal::SiX,
  SiCoordsFloatVal::SiY, and SiCoordsFloatVal::SiZ and also defines
  SiCoordsFloatVal::SIZE, which indicates the number of enumeration
  constants defined.

  Note that this macro only works with preprocessors that do C99-style
  variadic macros. */
#define KMC_MAKE_LATTICE_FLOATVAL_ENUM(EnumName, ...)	\
  struct EnumName##FloatVal {			\
    enum Type {__VA_ARGS__, SIZE};			\
  }

/*! Convenience macro for defining enumeration constants for the
    offsets used in a cell propensity calculation.

    For example, the following

    \code
    KMC_MAKE_OFFSET_ENUM(HopOffset,
                         UP, DOWN, LEFT, RIGHT);
    \endcode

    creates an enumeration that not only contains the values
    HopOffset::UP, HopOffset::DOWN, HopOffset::LEFT, and
    HopOffset::RIGHT, but also the values HopOffset::SELF -- <EM>which
    always has an integer value of zero</EM> -- and HopOffset::SIZE,
    which indicates the number of enumeration constants defined.

    Note that this macro only works with preprocessors that do
    C99-style variadic macros.

    \see CellNeighOffsets CellNeighProbe
*/
#define KMC_MAKE_OFFSET_ENUM(EnumName, ...)	\
  struct EnumName {				\
    enum Type {SELF, __VA_ARGS__, SIZE};	\
  }


/*! Convenience macro for defining enumeration constants that identify
    possible events.

    For example, the following

    \code
    KMC_MAKE_ID_ENUM(MyCellCenteredEvents,
                     HOP_UP,
		     HOP_DOWN,
		     HOP_LEFT,
		     HOP_RIGHT);
    \endcode

    creates an enumeration with values MyCellCenteredEvents::HOP_UP,
    MyCellCenteredEvents::HOP_DOWN, MyCellCenteredEvents::HOP_LEFT,
    and MyCellCenteredEvents::HOP_RIGHT and also defines
    MyCellCenteredEvents::SIZE, which indicates the number of enumeration
    constants defined.

    Note that this macro only works with preprocessors that do
    C99-style variadic macros.

    \see Simulation::addCellCenteredEvent()
 */
#define KMC_MAKE_ID_ENUM(EnumName, ...)	\
  struct EnumName {					\
    enum Type {__VA_ARGS__, SIZE};			\
  }

#endif /* MAKE_ENUM_HPP */
