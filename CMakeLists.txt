cmake_minimum_required( VERSION 2.8 )

project(EM-MIX-CNP)

if(NOT CMAKE_BUILD_TYPE)
  message("Type of build not specified. Use 'Release' by default.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
endif(NOT CMAKE_BUILD_TYPE)

FIND_PACKAGE(Boost REQUIRED system random thread date_time COMPONENTS)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")

set( LIBLEMON_ROOT "$ENV{HOME}/lemon" CACHE PATH "Additional search directory for lemon library" )

if( NOT CPLEX_INC_DIR OR NOT CPLEX_LIB_DIR OR NOT CONCERT_INC_DIR OR NOT CONCERT_LIB_DIR )
  set( CPLEX_AUTODETECTION_PREFIXES "$ENV{HOME}/ILOG;/opt/ibm/ILOG;$ENV{HOME}/Applications/IBM/ILOG;/research/compbio/software/ILOG;/research/compbio/software/ILOG" CACHE STRING "CPLEX automatic lookup directories" )
  foreach( CplexInstallPrefix ${CPLEX_AUTODETECTION_PREFIXES} )
    file( GLOB_RECURSE CplexLibraryPaths ${CplexInstallPrefix}/*/cplex/lib/*/libcplex.a )
    if( CplexLibraryPaths )
      list( LENGTH CplexLibraryPaths NumPath )
      if( NumPath GREATER 1 )
        message( WARNING "Found multiple CPLEX installations in ${CplexInstallPrefix}, automatically selecting the first one, please set {CPLEX|CONCERT}_{INC|LIB}_DIR manually if this is not the right selection" )
      else()
        message( "Found CPLEX's installation in ${CPlexInstallPrefix}" )
      endif()
      list( GET CplexLibraryPaths 0 CplexLibraryPath )

      set( MyRegex "${CplexInstallPrefix}/(.*)/cplex/?/lib/(.*)/libcplex.a" )
      string( REGEX MATCH "${MyRegex}" ValidatedCplexLibPath ${CplexLibraryPath} )
      if( NOT ValidatedCplexLibPath )
        message( WARNING "Unable to recognize ${CplexLibraryPath} as a valid CPLEX path" )
      else()
        string( REGEX REPLACE "${MyRegex}" "\\1" CplexHome ${ValidatedCplexLibPath} )
        string( REGEX REPLACE "${MyRegex}" "\\2" CplexArch ${ValidatedCplexLibPath} )

        set( CPLEX_INC_DIR   "${CplexInstallPrefix}/${CplexHome}/cplex/include"            CACHE PATH "CPLEX's include dir" )
        set( CPLEX_LIB_DIR   "${CplexInstallPrefix}/${CplexHome}/cplex/lib/${CplexArch}"   CACHE PATH "CPLEX's library dir" )
        set( CONCERT_INC_DIR "${CplexInstallPrefix}/${CplexHome}/concert/include"          CACHE PATH "Concert's include dir" )
        set( CONCERT_LIB_DIR "${CplexInstallPrefix}/${CplexHome}/concert/lib/${CplexArch}" CACHE PATH "Concert's library dir" )

        break()
      endif()
    else()
      message( "Unable to find CPLEX and Concert installation in ${CplexInstallPrefix}" )
    endif()
  endforeach()

  if( NOT CPLEX_INC_DIR OR NOT CPLEX_LIB_DIR OR NOT CONCERT_INC_DIR OR NOT CONCERT_LIB_DIR )
    message( "Unable to autodeted CPLEX's and Concert's installation directories, please set {CPLEX|CONCERT}_{INC|LIB}_DIR manually" )
  else()
    message( "Autodetected CPLEX's include dir   : ${CPLEX_INC_DIR}" )
    message( "Autodetected CPLEX's library dir   : ${CPLEX_LIB_DIR}" )
    message( "Autodetected Concert's include dir : ${CONCERT_INC_DIR}" )
    message( "Autodetected Concert's library dir : ${CONCERT_LIB_DIR}" )
  endif()
endif()


add_definitions( -DIL_STD )

#COMMENT THIS LINE and UNCOMMENT THE NEXT IF YOU ARE USING A VERSION OF CPLEX < 12.7 ON MAC OSX
#SET(CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LIBRARY "libstdc++")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
#  set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libstdc++" )

include_directories( ${CPLEX_INC_DIR} ${CONCERT_INC_DIR} ${Boost_INCLUDE_DIRS} "${LIBLEMON_ROOT}/include" "src" )
link_directories( ${CPLEX_LIB_DIR} ${CONCERT_LIB_DIR} ${Boost_LIBRARIES} "${LIBLEMON_ROOT}/lib" )

set( CommonLibs
  ilocplex
  cplex
  concert
  pthread
  emon
)

add_executable ( mixcnp
  src/mixcnp.cpp
  src/worker.h
  src/worker.cpp
  src/basecarchitect.h
  src/basecarchitect.cpp
  src/carchitect.h
  src/carchitect.cpp
  src/marchitect.h
  src/marchitect.cpp	
  src/basic_types.h
  src/basic_types.cpp
  src/copynumbertree.cpp
  src/copynumbertree.h
  src/inputinstance.h
  src/inputinstance.cpp
  src/fmcsolution.h
  src/fmcsolution.cpp
  src/manager.h
  src/manager.cpp
  src/refiner.h
  src/refiner.cpp
)
target_link_libraries(mixcnp ${CommonLibs} ${Boost_LIBRARIES})

add_executable ( visualize
  src/visualize.cpp
  src/basic_types.h
  src/basic_types.cpp
  src/copynumbertree.cpp
  src/copynumbertree.h
  src/fmcsolution.h
  src/fmcsolution.cpp
  src/inputinstance.h
  src/inputinstance.cpp
)
target_link_libraries(visualize ${CommonLibs} ${Boost_LIBRARIES})

add_executable ( compare
  src/compare.cpp
  src/comparison.h
  src/comparison.cpp
  src/basic_types.h
  src/basic_types.cpp
  src/copynumbertree.cpp
  src/copynumbertree.h
  src/fmcsolution.h
  src/fmcsolution.cpp
  src/inputinstance.h
  src/inputinstance.cpp
  src/tripletarchitect.cpp
  src/tripletarchitect.h
)
target_link_libraries(compare ${CommonLibs} ${Boost_LIBRARIES})

add_executable ( triplet
  src/triplet.cpp
  src/basic_types.h
  src/basic_types.cpp
  src/copynumbertree.cpp
  src/copynumbertree.h
  src/tripletarchitect.cpp
  src/tripletarchitect.h
)
target_link_libraries(triplet ${CommonLibs} ${Boost_LIBRARIES})

add_executable( check_carchitect EXCLUDE_FROM_ALL 
  test/check_carchitect.cpp 
  test/basic_check.h
  test/basic_check.cpp
  src/basecarchitect.h
  src/basecarchitect.cpp
  src/carchitect.h
  src/carchitect.cpp
  src/basic_types.h
  src/basic_types.cpp
  src/copynumbertree.cpp
  src/copynumbertree.h
  src/inputinstance.h
  src/inputinstance.cpp
)
target_link_libraries( check_carchitect ${CommonLibs} ${Boost_LIBRARIES} )

add_executable( check_marchitect EXCLUDE_FROM_ALL
  test/check_marchitect.cpp
  test/basic_check.h
  test/basic_check.cpp
  src/marchitect.h
  src/marchitect.cpp
  src/basic_types.h
  src/basic_types.cpp
  src/copynumbertree.cpp
  src/copynumbertree.h
  src/inputinstance.h
  src/inputinstance.cpp
)
target_link_libraries( check_marchitect ${CommonLibs} ${Boost_LIBRARIES} )

add_executable( check_worker EXCLUDE_FROM_ALL
  test/check_worker.cpp
  test/basic_check.h
  test/basic_check.cpp
  src/worker.h
  src/worker.cpp
  src/basecarchitect.h
  src/basecarchitect.cpp
  src/carchitect.h
  src/carchitect.cpp
  src/marchitect.h
  src/marchitect.cpp
  src/basic_types.h
  src/basic_types.cpp
  src/copynumbertree.cpp
  src/copynumbertree.h
  src/inputinstance.h
  src/inputinstance.cpp
)
target_link_libraries( check_worker ${CommonLibs} ${Boost_LIBRARIES} )

add_executable( check_manager EXCLUDE_FROM_ALL
  test/check_manager.cpp
  test/basic_check.h
  test/basic_check.cpp
  src/manager.h
  src/manager.cpp
  src/worker.h
  src/worker.cpp
  src/basecarchitect.h
  src/basecarchitect.cpp
  src/carchitect.h
  src/carchitect.cpp
  src/marchitect.h
  src/marchitect.cpp
  src/basic_types.h
  src/basic_types.cpp
  src/copynumbertree.cpp
  src/copynumbertree.h
  src/inputinstance.h
  src/inputinstance.cpp
  src/fmcsolution.h
  src/fmcsolution.cpp
  src/refiner.h
  src/refiner.cpp
)
target_link_libraries( check_manager ${CommonLibs} ${Boost_LIBRARIES} )


enable_testing()
add_test( NAME CheckCArchitect COMMAND check_carchitect )
add_test( NAME CheckMArchitect COMMAND check_marchitect )
add_test( NAME CheckWorker COMMAND check_worker )
add_test( NAME CheckManager COMMAND check_manager )
add_custom_target( check COMMAND ${CMAKE_CTEST_COMMAND} DEPENDS 
  check_carchitect
  check_marchitect
  check_worker
  check_manager
)
