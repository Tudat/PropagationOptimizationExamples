/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <Tudat/SimulationSetup/tudatSimulationHeader.h>
#include <Tudat/Astrodynamics/Ephemerides/approximatePlanetPositions.h>
#include <Tudat/External/SpiceInterface/spiceEphemeris.h>

#include "propagationAndOptimization/applicationOutput.h"


//! Execute propagation of orbit of LunarOrbiter around the Earth.
int main()
{
    std::string outputDirectory = tudat_applications::getOutputPath( "EnvironmentModels/" );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            USING STATEMENTS              //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    using namespace tudat;
    using namespace tudat::simulation_setup;
    using namespace tudat::propagators;
    using namespace tudat::numerical_integrators;
    using namespace tudat::orbital_element_conversions;
    using namespace tudat::basic_mathematics;
    using namespace tudat::unit_conversions;
    using namespace tudat::ephemerides;
    using namespace tudat::spice_interface;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////     CREATE ENVIRONMENT AND VEHICLE       //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Load Spice kernels.
    spice_interface::loadStandardSpiceKernels(
    { input_output::getSpiceKernelPath( ) + "de430.bsp" } );

    std::shared_ptr< Ephemeris > approximateEphemeris = std::make_shared< ApproximatePlanetPositions>(
                ApproximatePlanetPositionsBase::BodiesWithEphemerisData::mars );
    std::shared_ptr< Ephemeris > spiceEphemeris = std::make_shared< SpiceEphemeris >(
                "Mars_Barycenter", "Sun", false, false, false );
    //    std::shared_ptr< Ephemeris > spiceEphemerisSun = std::make_shared< SpiceEphemeris >(
    //                    "Sun", "SSB", false, false, false );

    double sunGravitationalParameter = getBodyGravitationalParameter( "Sun" );
    double startTime = -400.0 * physical_constants::JULIAN_YEAR;
    double endTime = 400.0 * physical_constants::JULIAN_YEAR;
    double timeStep = 30.0 * physical_constants::JULIAN_DAY;

    std::map< double, Eigen::Vector6d > spiceStates;
    std::map< double, Eigen::Vector6d > approximateStates;

    double currentTime = startTime;
    while( currentTime < endTime )
    {
        spiceStates[ currentTime ] =
                convertCartesianToKeplerianElements( spiceEphemeris->getCartesianState( currentTime ), sunGravitationalParameter );
        approximateStates[ currentTime ] =
                convertCartesianToKeplerianElements( approximateEphemeris->getCartesianState( currentTime ), sunGravitationalParameter );

//        std::cout<<spiceStates[ currentTime ].transpose( )<<std::endl;
//        std::cout<<approximateStates[ currentTime ].transpose( )<<std::endl;

//        std::cout<<( approximateStates[ currentTime ] - spiceStates[ currentTime ] ).transpose( ) <<std::endl;

        currentTime += timeStep;
    }

    input_output::writeDataMapToTextFile( spiceStates,
                                          "marsSpiceStates.dat",
                                          outputDirectory );
    input_output::writeDataMapToTextFile( approximateStates,
                                          "marsApproximateStates.dat",
                                          outputDirectory );

    return EXIT_SUCCESS;

}
