 #    Copyright (c) 2010-2017, Delft University of Technology
 #    All rigths reserved
 #
 #    This file is part of the Tudat. Redistribution and use in source and
 #    binary forms, with or without modification, are permitted exclusively
 #    under the terms of the Modified BSD license. You should have received
 #    a copy of the license with this file. If not, please or visit:
 #    http://tudat.tudelft.nl/LICENSE.
 #

# Specify minimum CMake version required.
cmake_minimum_required(VERSION 2.6)

# Specify project name.
project(SatellitePropagatorExamples)

# Load UserSettings.txt
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR})
  message(STATUS "Building ${PROJECT_NAME} standalone.")
  include("${CMAKE_CURRENT_SOURCE_DIR}/UserSettings.txt" OPTIONAL)
else()
  message(STATUS "Building ${PROJECT_NAME} from within ${CMAKE_PROJECT_NAME}.")
  include("${CMAKE_CURRENT_SOURCE_DIR}/UserSettings.txt" OPTIONAL)
  include("${CMAKE_SOURCE_DIR}/UserSettings.txt" OPTIONAL)
  STRING(REGEX REPLACE ${CMAKE_SOURCE_DIR} "" RELATIVE_PROJECT_PATH ${CMAKE_CURRENT_SOURCE_DIR})
  set(RELATIVE_PROJECT_PATH "${RELATIVE_PROJECT_PATH}" CACHE STRING "Relative path wrt to project for function")
  message(STATUS "Relative path (wrt to project): ${RELATIVE_PROJECT_PATH}")
endif()

# Set CMake build-type. If it not supplied by the user (either directly as an argument of through
# the "UserSettings.txt" file, the default built type is "Release".
if((NOT CMAKE_BUILD_TYPE) OR (CMAKE_BUILD_TYPE STREQUAL "Release"))
    message(STATUS "WARNING: building release version!")
    set(CMAKE_BUILD_TYPE Release)
elseif(CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(STATUS "WARNING: building debug version!")
endif()

# Set module path to local CMake scripts.
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMakeModules)

# Set compiler based on preferences (e.g. USE_CLANG) and system.
include(${CMAKE_MODULE_PATH}/compiler.cmake)

# Define the directory with the source code.
set(SRCROOT "${CMAKE_CURRENT_SOURCE_DIR}")

# Define the code root directory.
set(CODEROOT "${CMAKE_CURRENT_SOURCE_DIR}/..")

# Set testing options based on platform.
enable_testing()

# Set lib and bin directories where static libraries and unit tests are built.
if(NOT LIBROOT)
  set(LIBROOT "${CODEROOT}/lib")
endif()
if(NOT BINROOT)
  set(BINROOT "${CODEROOT}/bin")
endif()

# Set the global macros for setting up targets.
macro(setup_executable_target target_name CUSTOM_OUTPUT_PATH)
  set_property(TARGET ${target_name} PROPERTY RUNTIME_OUTPUT_DIRECTORY "${BINROOT}/applications")
  install(TARGETS ${target_name} RUNTIME DESTINATION "${BINROOT}/applications")
endmacro(setup_executable_target)

macro(setup_library_target target_name CUSTOM_OUTPUT_PATH)
  set_property(TARGET ${target_name} PROPERTY LIBRARY_OUTPUT_DIRECTORY "${LIBROOT}")
  set_property(TARGET ${target_name} PROPERTY ARCHIVE_OUTPUT_DIRECTORY "${LIBROOT}")
endmacro(setup_library_target)

macro(setup_unit_test_target target_name CUSTOM_OUTPUT_PATH)
  set_property(TARGET ${target_name} PROPERTY RUNTIME_OUTPUT_DIRECTORY "${BINROOT}/unit_tests")
  get_property(CUSTOM_TEST_PROGRAM_NAME TARGET ${target_name} PROPERTY OUTPUT_NAME)
  add_test("${target_name}" "${BINROOT}/unit_tests/${target_name}")
endmacro(setup_unit_test_target)

# Define the install targets to create a distribution.
if(NOT TUDAT_BUNDLE_DISTRIBUTION_PATH)
    set(TUDAT_BUNDLE_DISTRIBUTION_PATH "${CODEROOT}")
endif(NOT TUDAT_BUNDLE_DISTRIBUTION_PATH)

if(NOT SATELLITE_PROPAGATOR_EXAMPLES_DISTRIBUTION_PATH)
    set(SATELLITE_PROPAGATOR_EXAMPLES_DISTRIBUTION_PATH 
        "${TUDAT_BUNDLE_DISTRIBUTION_PATH}/tudatApplications/satellitePropagatorExamples")
endif(NOT SATELLITE_PROPAGATOR_EXAMPLES_DISTRIBUTION_PATH)

# Install files.
install(DIRECTORY "${SRCROOT}/"
        DESTINATION "${SATELLITE_PROPAGATOR_EXAMPLES_DISTRIBUTION_PATH}/SatellitePropagatorExamples"
        PATTERN ".DS_STORE" EXCLUDE
        PATTERN "CMakeLists.txt.user" EXCLUDE
        PATTERN ".svn" EXCLUDE
        PATTERN ".git" EXCLUDE
        PATTERN ".bzr" EXCLUDE
)

# Include the top-level directories.
include_directories(AFTER
  "${CODEROOT}"
)

# Find Eigen3 library on local system.
find_package(Eigen3 REQUIRED)

# Include Eigen3 directories.
# Set CMake flag to suppress Eigen warnings (platform-dependent solution).
if(NOT APPLE)
  include_directories(SYSTEM AFTER "${EIGEN3_INCLUDE_DIR}")
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${EIGEN3_INCLUDE_DIR}\"")
endif()

# Configure Boost libraries.
if(NOT Boost_USE_STATIC_LIBS)
  set(Boost_USE_STATIC_LIBS ON)
endif()
if(NOT Boost_USE_MULTITHREADED)
  set(Boost_USE_MULTITHREADED ON)
endif()
if(NOT Boost_USE_STATIC_RUNTIME)
  set(Boost_USE_STATIC_RUNTIME ON)
endif()

# Find Boost libraries on local system.
find_package(Boost 1.45.0 
             COMPONENTS thread date_time system unit_test_framework filesystem regex REQUIRED)

# Include Boost directories.
# Set CMake flag to suppress Boost warnings (platform-dependent solution).
if(NOT APPLE)
  include_directories(SYSTEM AFTER "${Boost_INCLUDE_DIRS}")
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${Boost_INCLUDE_DIRS}\"")
endif()

# Find Tudat library on local system.
find_package(Tudat 2.0 REQUIRED)

# Include Tudat directories.
# Set CMake flag to suppress Tudat warnings (platform-dependent solution).
if(NOT APPLE)
  include_directories(SYSTEM AFTER "${TUDAT_INCLUDE_DIR}")
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${TUDAT_INCLUDE_DIR}\"")
endif()

  # Find CSPICE library on local system.
  find_package(Spice)

  # Include CSpice directories.
  if(NOT APPLE)
    include_directories(SYSTEM AFTER "${SPICE_INCLUDE_DIR}")
  else( )
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${SPICE_INCLUDE_DIR}\"")
  endif( )

option(USE_NRLMSISE00 "build Tudat with NRLMSISE-00 enabled" ON)
if(NOT USE_NRLMSISE00)
  message(STATUS "NRLMSISE-00 disabled!")
  add_definitions(-DUSE_NRLMSISE00=0)
else()
  message(STATUS "NRLMSISE-00 enabled!")
  add_definitions(-DUSE_NRLMSISE00=1)
  # Find USE_NRLMSISE00 library on local system.
  find_package(NRLMSISE00)

  # Include NRLMSISE00 directories.
  if(NOT APPLE)
    include_directories(SYSTEM AFTER "${NRLMSISE00_INCLUDE_DIR}")
  else( )
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${NRLMSISE00_INCLUDE_DIR}\"")
  endif( )
endif( )

list(APPEND TUDAT_APPLICATION_EXTERNAL_LIBRARIES "")
list(APPEND TUDAT_APPLICATION_EXTERNAL_INTERFACE_LIBRARIES "")
list(APPEND TUDAT_APPLICATION_ITRS_LIBRARIES "")


if(USE_SOFA)
 list(APPEND TUDAT_APPLICATION_EXTERNAL_LIBRARIES sofa)
 list(APPEND TUDAT_APPLICATION_EXTERNAL_INTERFACE_LIBRARIES tudat_sofa_interface )
 list(APPEND TUDAT_APPLICATION_ITRS_LIBRARIES tudat_earth_orientation )
endif()

if(USE_CSPICE)
  list(APPEND TUDAT_APPLICATION_EXTERNAL_LIBRARIES cspice)
  list(APPEND TUDAT_APPLICATION_EXTERNAL_INTERFACE_LIBRARIES tudat_spice_interface )
endif()

if(USE_NRLMSISE00)
  list(APPEND TUDAT_APPLICATION_EXTERNAL_LIBRARIES nrlmsise00)
endif()


list(APPEND TUDAT_APPLICATION_PROPAGATION_LIBRARIES tudat_simulation_setup tudat_ground_stations tudat_propagators
    tudat_aerodynamics tudat_system_models tudat_geometric_shapes tudat_relativity tudat_gravitation tudat_mission_segments
    tudat_electro_magnetism tudat_propulsion tudat_ephemerides ${TUDAT_APPLICATION_ITRS_LIBRARIES} tudat_numerical_integrators tudat_reference_frames
     tudat_statistics tudat_propagators ${TUDAT_APPLICATION_EXTERNAL_INTERFACE_LIBRARIES} tudat_basic_astrodynamics tudat_root_finders tudat_basic_mathematics
     tudat_input_output tudat_basics ${TUDAT_APPLICATION_EXTERNAL_LIBRARIES})

list(APPEND TUDAT_APPLICATION_ESTIMATION_LIBRARIES tudat_simulation_setup tudat_observation_models tudat_ground_stations tudat_acceleration_partials
    tudat_observation_partials tudat_estimatable_parameters tudat_orbit_determination  tudat_propagators
    tudat_aerodynamics tudat_system_models tudat_geometric_shapes tudat_relativity tudat_gravitation tudat_mission_segments
    tudat_electro_magnetism tudat_propulsion tudat_ephemerides ${TUDAT_APPLICATION_ITRS_LIBRARIES} tudat_numerical_integrators tudat_reference_frames
    tudat_statistics tudat_propagators ${TUDAT_APPLICATION_EXTERNAL_INTERFACE_LIBRARIES} tudat_basic_astrodynamics tudat_root_finders tudat_basic_mathematics
    tudat_input_output tudat_basics ${TUDAT_APPLICATION_EXTERNAL_LIBRARIES})



## INTEGRATION: SLIDE RESULTS
add_executable(po_application_LunarOrbiterPropagationIntegrationSettings "${SRCROOT}/NumericalIntegration/Generation/lunarOrbiterPropagatorIntegratorSettings.cpp")
setup_executable_target(po_application_LunarOrbiterPropagationIntegrationSettings "${SRCROOT}")
target_link_libraries(po_application_LunarOrbiterPropagationIntegrationSettings ${TUDAT_APPLICATION_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES} )

add_executable(po_application_KeplerOrbitErrorTrend "${SRCROOT}/NumericalIntegration/Generation/keplerOrbitTruncationAndRoundingErrorTrend.cpp")
setup_executable_target(po_application_KeplerOrbitErrorTrend "${SRCROOT}")
target_link_libraries(po_application_KeplerOrbitErrorTrend ${TUDAT_APPLICATION_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES} )


## EQUATIONS OF MOTION: SLIDE RESULTS
add_executable(po_application_PerturbedSatellitePropagationElementTypes "${SRCROOT}/EquationsOfMotion/Generation/perturbedSatellitePropagationElementTypes.cpp")
setup_executable_target(po_application_PerturbedSatellitePropagationElementTypes "${SRCROOT}")
target_link_libraries(po_application_PerturbedSatellitePropagationElementTypes ${TUDAT_APPLICATION_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES} )

add_executable(po_application_PerturbedSatellitePropagationElementTypesSingular "${SRCROOT}/EquationsOfMotion/Generation/perturbedSatellitePropagationElementTypesSingular.cpp")
setup_executable_target(po_application_PerturbedSatellitePropagationElementTypesSingular "${SRCROOT}")
target_link_libraries(po_application_PerturbedSatellitePropagationElementTypesSingular ${TUDAT_APPLICATION_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES} )


## ENVIRONMENT MODELS: SLIDE RESULTS
add_executable(po_application_EphemerisInfluence "${SRCROOT}/EnvironmentModels/Generation/planetaryEphemerisInfluence.cpp")
setup_executable_target(po_application_EphemerisInfluence "${SRCROOT}")
target_link_libraries(po_application_EphemerisInfluence ${TUDAT_APPLICATION_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES} )

add_executable(po_application_RotationModelInfluence "${SRCROOT}/EnvironmentModels/Generation/rotationModelInfluence.cpp")
setup_executable_target(po_application_RotationModelInfluence "${SRCROOT}")
target_link_libraries(po_application_RotationModelInfluence ${TUDAT_APPLICATION_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES} )


## ACCELERATION MODELS: SLIDE RESULTS
add_executable(po_application_SphericalHarmonicCase "${SRCROOT}/AccelerationModels/Generation/sphericalHarmonicInfluence.cpp")
setup_executable_target(po_application_SphericalHarmonicCase "${SRCROOT}")
target_link_libraries(po_application_SphericalHarmonicCase ${TUDAT_APPLICATION_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES} )

add_executable(po_application_SphericalHarmonicCaseSolarSystem "${SRCROOT}/AccelerationModels/Generation/sphericalHarmonicInfluenceSolarSystem.cpp")
setup_executable_target(po_application_SphericalHarmonicCaseSolarSystem "${SRCROOT}")
target_link_libraries(po_application_SphericalHarmonicCaseSolarSystem ${TUDAT_APPLICATION_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES} )

add_executable(po_application_PropagationOrigin "${SRCROOT}/AccelerationModels/Generation/propagationOriginInfluence.cpp")
setup_executable_target(po_application_PropagationOrigin "${SRCROOT}")
target_link_libraries(po_application_PropagationOrigin ${TUDAT_APPLICATION_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES} )

add_executable(po_application_MutualGravitationalInteraction "${SRCROOT}/AccelerationModels/Generation/mutualGravitationalInteraction.cpp")
setup_executable_target(po_application_MutualGravitationalInteraction "${SRCROOT}")
target_link_libraries(po_application_MutualGravitationalInteraction ${TUDAT_APPLICATION_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES} )

#add_executable(po_application_MutualExtendedGravitationalInteraction "${SRCROOT}/AccelerationModels/Generation/mutualExtendedGravitationalInteractionPhobos.cpp")
#setup_executable_target(po_application_MutualExtendedGravitationalInteraction "${SRCROOT}")
#target_link_libraries(po_application_MutualExtendedGravitationalInteraction ${TUDAT_APPLICATION_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES} )

add_executable(po_application_PhobosRadiationPressure "${SRCROOT}/AccelerationModels/Generation/phobosOrbitPropagationWithRadiationPressure.cpp")
setup_executable_target(po_application_PhobosRadiationPressure "${SRCROOT}")
target_link_libraries(po_application_PhobosRadiationPressure ${TUDAT_APPLICATION_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES} )

add_executable(po_application_SphericalHarmonicCaseReEnrty "${SRCROOT}/AccelerationModels/Generation/apolloCapsuleEntrySphericalHarmonicInfluence.cpp")
setup_executable_target(po_application_SphericalHarmonicCaseReEnrty "${SRCROOT}")
target_link_libraries(po_application_SphericalHarmonicCaseReEnrty ${TUDAT_APPLICATION_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES} )


## UNCERTAINTY MODELLING: SLIDE RESULTS
add_executable(po_application_InitialStateAbsorptionEntry "${SRCROOT}/UncertaintyModelling/Generation/initialStateAbsorptionOfJ2_Entry.cpp")
setup_executable_target(po_application_InitialStateAbsorptionEntry "${SRCROOT}")
target_link_libraries(po_application_InitialStateAbsorptionEntry ${TUDAT_APPLICATION_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES} )

add_executable(po_application_InitialStateAbsorptionSolarSystem "${SRCROOT}/UncertaintyModelling/Generation/initialStateAbsorptionOfJ2_SolarSystem.cpp")
setup_executable_target(po_application_InitialStateAbsorptionSolarSystem "${SRCROOT}")
target_link_libraries(po_application_InitialStateAbsorptionSolarSystem ${TUDAT_APPLICATION_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES} )

add_executable(po_application_AccelerationDirectionInfluence "${SRCROOT}/UncertaintyModelling/Generation/accelerationDirectionInfluence.cpp")
setup_executable_target(po_application_AccelerationDirectionInfluence "${SRCROOT}")
target_link_libraries(po_application_AccelerationDirectionInfluence ${TUDAT_APPLICATION_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES} )




