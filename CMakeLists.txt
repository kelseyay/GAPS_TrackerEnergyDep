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

add_executable(exam minimum_example.cpp)
target_include_directories(exam PRIVATE ${GAPS_INCLUDE_DIRS}/gaps)
target_link_libraries(exam ${LIB_DIRS})

add_executable(Occupancy Occupancy.cpp)
target_include_directories(Occupancy PRIVATE ${GAPS_INCLUDE_DIRS}/gaps)
target_link_libraries(Occupancy ${LIB_DIRS})

add_executable(EdepHist FullTrackerEdepHist.cpp)
target_include_directories(EdepHist PRIVATE ${GAPS_INCLUDE_DIRS}/gaps)
target_link_libraries(EdepHist ${LIB_DIRS})

add_executable(TellMe TellMe.cpp)
target_include_directories(TellMe PRIVATE ${GAPS_INCLUDE_DIRS}/gaps)
target_link_libraries(TellMe ${LIB_DIRS})

add_executable(RecTellMe RecTellMe.cpp)
target_include_directories(RecTellMe PRIVATE ${GAPS_INCLUDE_DIRS}/gaps)
target_link_libraries(RecTellMe ${LIB_DIRS})

add_executable(FitHist LGFitHist.cpp)
target_include_directories(FitHist PRIVATE ${GAPS_INCLUDE_DIRS}/gaps)
target_link_libraries(FitHist ${LIB_DIRS})

add_executable(NoCutsEdepHist NoCutsFullTrackerEdepHist.cpp)
target_include_directories(NoCutsEdepHist PRIVATE ${GAPS_INCLUDE_DIRS}/gaps)
target_link_libraries(NoCutsEdepHist ${LIB_DIRS})

add_executable(Open2Droot Open2Droot.cpp)
target_include_directories(Open2Droot PRIVATE ${GAPS_INCLUDE_DIRS}/gaps)
target_link_libraries(Open2Droot ${LIB_DIRS})

add_executable(ArbOpen2Droot Open2DArbroot.cpp)
target_include_directories(ArbOpen2Droot PRIVATE ${GAPS_INCLUDE_DIRS}/gaps)
target_link_libraries(ArbOpen2Droot ${LIB_DIRS})

add_executable(RankPulse PulsedRanking.cpp)
target_include_directories(RankPulse PRIVATE ${GAPS_INCLUDE_DIRS}/gaps)
target_link_libraries(RankPulse ${LIB_DIRS})

add_executable(RankPulse2 PulsedRanking2.cpp)
target_include_directories(RankPulse2 PRIVATE ${GAPS_INCLUDE_DIRS}/gaps)
target_link_libraries(RankPulse2 ${LIB_DIRS})

add_executable(OneHist TrackerOneHist.cpp)
target_include_directories(OneHist PRIVATE ${GAPS_INCLUDE_DIRS}/gaps)
target_link_libraries(OneHist ${LIB_DIRS})

add_executable(MCOneDet MCOneDetHist.cpp)
target_include_directories(MCOneDet PRIVATE ${GAPS_INCLUDE_DIRS}/gaps)
target_link_libraries(MCOneDet ${LIB_DIRS})

add_executable(OneDet OneDetHist.cpp)
target_include_directories(OneDet PRIVATE ${GAPS_INCLUDE_DIRS}/gaps)
target_link_libraries(OneDet ${LIB_DIRS})

add_executable(Zedep ChargeMCEdep.cpp)
target_include_directories(Zedep PRIVATE ${GAPS_INCLUDE_DIRS}/gaps)
target_link_libraries(Zedep ${LIB_DIRS})

add_executable(Zedep2 ChargeMCEdep2.cpp)
target_include_directories(Zedep2 PRIVATE ${GAPS_INCLUDE_DIRS}/gaps)
target_link_libraries(Zedep2 ${LIB_DIRS})

add_executable(Zedep3 ChargeMCEdep3.cpp)
target_include_directories(Zedep3 PRIVATE ${GAPS_INCLUDE_DIRS}/gaps)
target_link_libraries(Zedep3 ${LIB_DIRS})

add_executable(DZedep ChargeDataEdep.cpp)
target_include_directories(DZedep PRIVATE ${GAPS_INCLUDE_DIRS}/gaps)
target_link_libraries(DZedep ${LIB_DIRS})

add_executable(DZedep2 ChargeDataEdep2.cpp)
target_include_directories(DZedep2 PRIVATE ${GAPS_INCLUDE_DIRS}/gaps)
target_link_libraries(DZedep2 ${LIB_DIRS})

add_executable(DZedep3 ChargeDataEdep3.cpp)
target_include_directories(DZedep3 PRIVATE ${GAPS_INCLUDE_DIRS}/gaps)
target_link_libraries(DZedep3 ${LIB_DIRS})

add_executable(Ledep MCEdepLayers.cpp)
target_include_directories(Ledep PRIVATE ${GAPS_INCLUDE_DIRS}/gaps)
target_link_libraries(Ledep ${LIB_DIRS})

add_executable(MCCutEff CutsMCEdep.cpp)
target_include_directories(MCCutEff PRIVATE ${GAPS_INCLUDE_DIRS}/gaps)
target_link_libraries(MCCutEff ${LIB_DIRS})

add_executable(RecCutEff CutsRecEdep.cpp)
target_include_directories(RecCutEff PRIVATE ${GAPS_INCLUDE_DIRS}/gaps)
target_link_libraries(RecCutEff ${LIB_DIRS})

add_executable(OccuOpen OpenOccu.cpp)
target_include_directories(OccuOpen PRIVATE ${GAPS_INCLUDE_DIRS}/gaps)
target_link_libraries(OccuOpen ${LIB_DIRS})

add_executable(TrmvB TrMvsBeta.cpp)
target_include_directories(TrmvB PRIVATE ${GAPS_INCLUDE_DIRS}/gaps)
target_link_libraries(TrmvB ${LIB_DIRS})

add_executable(MCEdepvB MCEdepvsBeta.cpp)
target_include_directories(MCEdepvB PRIVATE ${GAPS_INCLUDE_DIRS}/gaps)
target_link_libraries(MCEdepvB ${LIB_DIRS})

add_executable(MCZvB MCZvsBeta.cpp)
target_include_directories(MCZvB PRIVATE ${GAPS_INCLUDE_DIRS}/gaps)
target_link_libraries(MCZvB ${LIB_DIRS})

add_executable(RecEdepvB RecEdepvsBeta.cpp)
target_include_directories(RecEdepvB PRIVATE ${GAPS_INCLUDE_DIRS}/gaps)
target_link_libraries(RecEdepvB ${LIB_DIRS})

add_executable(PartsRecEdepvB PartsRecEdepvsBeta.cpp)
target_include_directories(PartsRecEdepvB PRIVATE ${GAPS_INCLUDE_DIRS}/gaps)
target_link_libraries(PartsRecEdepvB ${LIB_DIRS})

add_executable(MCPartsEdepvB MCPartsEdepvsBeta.cpp)
target_include_directories(MCPartsEdepvB PRIVATE ${GAPS_INCLUDE_DIRS}/gaps)
target_link_libraries(MCPartsEdepvB ${LIB_DIRS})

add_executable(SimpleHisto SimpleHisto.cpp)
target_include_directories(SimpleHisto PRIVATE ${GAPS_INCLUDE_DIRS}/gaps)
target_link_libraries(SimpleHisto ${LIB_DIRS})

add_executable(OccNoCuts OccuNoCuts.cpp)
target_include_directories(OccNoCuts PRIVATE ${GAPS_INCLUDE_DIRS}/gaps)
target_link_libraries(OccNoCuts ${LIB_DIRS})
