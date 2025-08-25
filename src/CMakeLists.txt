####################
# Cache variables
####################

set(BOOST_ROOT "/usr/local" CACHE PATH "Path to root of Boost installation")

set(KMC_BUILD_SERIAL TRUE CACHE BOOL "Indicates whether to attempt to build the serial version of the library")
set(KMC_BUILD_PARALLEL TRUE CACHE BOOL "Indicates whether to attempt to build the parallel version of the library")

set(KMC_USE_DCMT TRUE CACHE BOOL "Indicates whether to attempt to use the Dynamic Creator Mersenne Twister as a random number generator")
set(DCMT_ROOT "/usr/local/dcmt" CACHE PATH "Path to root of installation of Dynamic Creator Mersenne Twister")

set(KMC_USE_RNGSTREAMS TRUE CACHE BOOL "Indicates whether to attempt to use the RngStreams library as a random number generator")
set(RNGSTREAMS_ROOT "/usr/local/rngStreams" CACHE PATH "Path to root of installation of RngStreams")

set(BUILD_SHARED_LIBS TRUE CACHE BOOL "Indicates whether the library should be static or shared")

# Whether to use rpath
set(KMC_THIN_FILM_USE_RPATH TRUE CACHE BOOL "Indicates if paths to linked libraries are in executables")

# For debugging
set(KMC_AVOID_BOOST_BIMAP FALSE CACHE BOOL "Indicates whether to implement bimaps with Boost.Bimap or an ad hoc implementation based on Boost.MultiIndex")
mark_as_advanced(KMC_AVOID_BOOST_BIMAP)

##################################################
# Rpath variables, so that I don't have to set 
# LD_LIBRARY_PATH before running
##################################################

if (KMC_THIN_FILM_USE_RPATH)
  set(CMAKE_SKIP_BUILD_RPATH FALSE)
  set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
endif()

####################
# Library
####################

set(KMC_CPP_FILES
  CellInds.cpp
  CellsToChange.cpp
  CellNeighProbe.cpp
  ErrorHandling.cpp
  EventExecutorGroup.cpp
  EventId.cpp
  IJK.cpp
  Lattice.cpp
  CellNeighOffsets.cpp
  RandNumGenMT19937.cpp
  Simulation.cpp
  Solver.cpp
  SolverBinaryTree.cpp
  SolverDynamicSchulze.cpp
  SolverFactory.cpp
  TimeIncrSchemeVars.cpp
  wrapInd.cpp)

set(KMC_HPP_FILES
  CellInds.hpp
  CellsToChange.hpp
  CellNeighProbe.hpp
  ErrorHandling.hpp
  EventExecutor.hpp
  EventExecutorGroup.hpp
  IdsOfSolvers.hpp
  IJK.hpp
  Lattice.hpp
  MakeEnum.hpp
  PeriodicAction.hpp
  CellNeighOffsets.hpp
  CellCenteredGroupPropensities.hpp
  RandNumGen.hpp
  RandNumGenMT19937.hpp
  Simulation.hpp
  SimulationState.hpp
  TimeIncrSchemeVars.hpp)

if (KMC_USE_DCMT)
  find_library(DCMT_LIBRARY "dcmt" "${DCMT_ROOT}/lib")

  if (DCMT_LIBRARY)
    set(KMC_CPP_FILES ${KMC_CPP_FILES} RandNumGenDCMT.cpp)
    set(KMC_HPP_FILES ${KMC_HPP_FILES} RandNumGenDCMT.hpp)

    include_directories("${DCMT_ROOT}/include")
  else()
    message("No DCMT library found under directory ${DCMT_ROOT}/lib. Will not use DCMT.")
  endif()
endif()

if (KMC_USE_RNGSTREAMS)
  find_library(RNGSTREAMS_LIBRARY "rngstreams" "${RNGSTREAMS_ROOT}/lib")

  if (RNGSTREAMS_LIBRARY)
    set(KMC_CPP_FILES ${KMC_CPP_FILES} RandNumGenRngStreams.cpp)
    set(KMC_HPP_FILES ${KMC_HPP_FILES} RandNumGenRngStreams.hpp)

    include_directories("${RNGSTREAMS_ROOT}/include")
  else()
    message("No RngStreams library found under directory ${RNGSTREAMS_ROOT}/lib. Will not use RngStreams.")
  endif()
endif()

if (KMC_AVOID_BOOST_BIMAP)
  add_definitions(-DKMC_AVOID_BOOST_BIMAP)
  message("Avoiding use of Boost.Bimap.")
endif()

# If find_package tries to use the CMake that comes with some installations of Boost,
# some things break horribly.
set(Boost_NO_BOOST_CMAKE TRUE)
find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

if (KMC_BUILD_SERIAL)

  add_library(KMCThinFilmSerial ${KMC_CPP_FILES})

  set(KMC_PARALLEL 0)
  configure_file("${PROJECT_SOURCE_DIR}/src/KMC_Config.hpp.in"
    "${PROJECT_BINARY_DIR}/serial/KMC_Config.hpp")

  set_property(TARGET KMCThinFilmSerial
    APPEND
    PROPERTY INCLUDE_DIRECTORIES "${PROJECT_BINARY_DIR}/serial")

  set_target_properties(KMCThinFilmSerial PROPERTIES
    VERSION "${KMC_THIN_FILM_MAJOR_VERSION}.${KMC_THIN_FILM_MINOR_VERSION}.${KMC_THIN_FILM_SUB_VERSION}"
    SOVERSION ${KMC_THIN_FILM_MAJOR_VERSION})

  if (KMC_USE_DCMT)
    target_link_libraries(KMCThinFilmSerial ${DCMT_LIBRARY})
  endif()

  if (KMC_USE_RNGSTREAMS)
    target_link_libraries(KMCThinFilmSerial ${RNGSTREAMS_LIBRARY})
  endif()

  install(TARGETS KMCThinFilmSerial DESTINATION lib)
  install(FILES "${PROJECT_BINARY_DIR}/serial/KMC_Config.hpp"
    DESTINATION include/KMCThinFilm/serial)

endif()

if (KMC_BUILD_PARALLEL)
  find_package(MPI)
  if (MPI_CXX_FOUND)

    add_library(KMCThinFilmParallel ${KMC_CPP_FILES})

    set(KMC_PARALLEL 1)
    configure_file("${PROJECT_SOURCE_DIR}/src/KMC_Config.hpp.in"
      "${PROJECT_BINARY_DIR}/parallel/KMC_Config.hpp")

    set_property(TARGET KMCThinFilmParallel
      APPEND
      PROPERTY INCLUDE_DIRECTORIES "${PROJECT_BINARY_DIR}/parallel")

    set_target_properties(KMCThinFilmParallel PROPERTIES
      VERSION "${KMC_THIN_FILM_MAJOR_VERSION}.${KMC_THIN_FILM_MINOR_VERSION}.${KMC_THIN_FILM_SUB_VERSION}"
      SOVERSION ${KMC_THIN_FILM_MAJOR_VERSION})

    include_directories(${MPI_CXX_INCLUDE_PATH})
    target_link_libraries(KMCThinFilmParallel ${MPI_CXX_LIBRARIES})

    set_property(TARGET KMCThinFilmParallel
      APPEND
      PROPERTY COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS}")

    set_property(TARGET KMCThinFilmParallel
      APPEND
      PROPERTY LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")

    if (KMC_USE_DCMT)
      target_link_libraries(KMCThinFilmParallel ${DCMT_LIBRARY})
    endif()

    if (KMC_USE_RNGSTREAMS)
      target_link_libraries(KMCThinFilmParallel ${RNGSTREAMS_LIBRARY})
    endif()

    install(TARGETS KMCThinFilmParallel DESTINATION lib)
    install(FILES "${PROJECT_BINARY_DIR}/parallel/KMC_Config.hpp"
      DESTINATION include/KMCThinFilm/parallel)

    # Should these "includes" be put outside of the "if (MPI_CXX_FOUND)" statement?
    include(CMakePushCheckState)
    include(CheckCXXSourceCompiles)

    cmake_push_check_state()
    set(CMAKE_REQUIRED_FLAGS ${MPI_CXX_COMPILE_FLAGS})
    set(CMAKE_REQUIRED_INCLUDES ${MPI_CXX_INCLUDE_PATH})
    set(CMAKE_REQUIRED_LIBRARIES ${MPI_CXX_LIBRARIES})

    check_cxx_source_compiles(
      "#include <mpi.h>
int main() {
   const void * buf;
   MPI_Request request;
   MPI_Isend(buf, 0, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &request);
   return 0;
}" MPI_ISEND_IS_CONST_CORRECT)

   if (MPI_ISEND_IS_CONST_CORRECT)
     message("MPI_Isend appears to be const correct.")
   else()
     message("MPI_Isend is not const correct. Working around the issue ...")
     add_definitions(-DMPI_ISEND_NOT_CONST_CORRECT)
   endif()

   cmake_pop_check_state()

  else()
    message("MPI not found, so will not build parallel version of library.")
  endif()
else()
  message("Will not attempt to build parallel version of library.")
endif()

install(FILES
  ${KMC_HPP_FILES}
  DESTINATION include/KMCThinFilm)
