#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.13 FATAL_ERROR)
project(minimal_example)

#----------------------------------------------------------------------------

find_package(ROOT REQUIRED)
find_package(GAPS REQUIRED)
find_package(Boost 1.74 REQUIRED program_options timer log log_setup)

add_definitions(-DUSE_BOOST_PROGRAM_OPTIONS)
add_definitions(-DUSE_BOOST_TIMER)
add_definitions(-DUSE_BOOST_LOG)
add_definitions(-DBOOST_LOG_DYN_LINK)

include(${ROOT_USE_FILE})
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

#----------------------------------------------------------------------------
# Add the executable, and link it to the ROOT/GAPS Libs
#

set(LIB_DIRS ${ROOT_LIBRARIES}
                            ${Geant4_LIBRARIES}
                            ${GAPS_LIBRARIES}/libGAPSCommon.so
                            ${GAPS_LIBRARIES}/libGAPSAnalysis.so
                            ${GAPS_LIBRARIES}/libGAPSReconstruction.so
                            Boost::program_options
                            Boost::timer
                            Boost::log
                            Boost::log_setup)

add_executable(example minimum_example.cpp)
target_include_directories(example PRIVATE ${GAPS_INCLUDE_DIRS}/gaps)
target_link_libraries(example ${LIB_DIRS})

add_executable(Occupancy Occupancy.cpp)
target_include_directories(Occupancy PRIVATE ${GAPS_INCLUDE_DIRS}/gaps)
target_link_libraries(Occupancy ${LIB_DIRS})
