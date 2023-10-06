cmake_minimum_required(VERSION 3.12)
project(warwick-legend VERSION 0.1.0)

# Compiler/build settings
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_BUILD_TYPE RelWithDebInfo)

# CTest if we need it
include(CTest)

# Dependencies
find_package(Geant4 10.6 REQUIRED gdml ui_all vis_all)

find_package(ROOT 6.06 CONFIG REQUIRED )


# Build
add_executable(warwick-legend
  warwick-legend.cc
  src/WLGDActionInitialization.cc
  src/WLGDBiasChangeCrossSection.cc
  src/WLGDBiasMultiParticleChangeCrossSection.cc
  src/WLGDCrystalHit.cc
  src/WLGDCrystalSD.cc
  src/WLGDDetectorConstruction.cc
  src/WLGDEventAction.cc
  src/WLGDPrimaryGeneratorAction.cc
  src/WLGDPetersGammaCascadeReader.cc
  src/WLGDRunAction.cc
  src/WLGDStackingAction.cc
  src/WLGDSteppingAction.cc
  src/WLGDTrackingAction.cc
  src/WLGDTrajectory.cc)
target_include_directories(warwick-legend PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(warwick-legend PRIVATE ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})

# Copy macro needed to run in interactive mode to build directory.
# By default, the macro is assumed to be in the working directory
# where warwick-legend is run from.
configure_file(vis.mac vis.mac COPYONLY)

# Test
if(BUILD_TESTING)
  add_subdirectory(test)
endif()
