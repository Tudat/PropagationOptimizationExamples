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

#include "propagationAndOptimization/applicationOutput.h"


//! Execute propagation of orbit of LunarOrbiter around the Earth.
int main()
{
    std::string outputDirectory = tudat_applications::getOutputPath( "AccelerationModels/" );

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
    spice_interface::loadStandardSpiceKernels( );

    //for( unsigned int ephemerisCase = 0; ephemerisCase < 3; ephemerisCase++ )

    // Create body objects.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Earth" );
    bodiesToCreate.push_back( "Moon" );
    bodiesToCreate.push_back( "Sun" );
    bodiesToCreate.push_back( "Jupiter" );
    bodiesToCreate.push_back( "Mars" );
    bodiesToCreate.push_back( "Venus" );

    std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( bodiesToCreate, -3600.0, 14.0 * tudat::physical_constants::JULIAN_DAY + 3600.0 );

    // Create Earth object
    NamedBodyMap bodyMap = createBodies( bodySettings );

    // Create spacecraft object.
    bodyMap[ "LunarOrbiter" ] = boost::make_shared< simulation_setup::Body >( );

    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );
    for( unsigned int originCase = 0; originCase < 4; originCase++ )
    {

        for( unsigned int accelerationCase = 0; accelerationCase < 3; accelerationCase++ )
        {
            std::cout<<"Test case: "<<accelerationCase<<" "<<originCase<<std::endl;
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            ///////////////////////            CREATE ACCELERATIONS          //////////////////////////////////////////////////////
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            // Define propagator settings variables.
            SelectedAccelerationMap accelerationMap;
            std::vector< std::string > bodiesToPropagate;
            std::vector< std::string > centralBodies;

            bodiesToPropagate.push_back( "LunarOrbiter" );
            if( originCase == 0 )
            {
                centralBodies.push_back( "Moon" );
            }
            else if( originCase == 1 )
            {
                centralBodies.push_back( "Earth" );
            }
            else if( originCase == 2 )
            {
                centralBodies.push_back( "Sun" );
            }
            else if( originCase == 3 )
            {
                centralBodies.push_back( "SSB" );
            }


            // Define propagation settings.
            std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfLunarOrbiter;

            if( accelerationCase == 0 )
            {
                accelerationsOfLunarOrbiter[ "Moon" ].push_back( boost::make_shared< AccelerationSettings >(
                                                                     basic_astrodynamics::central_gravity ) );
            }
            else if( accelerationCase == 1 )
            {
                accelerationsOfLunarOrbiter[ "Earth" ].push_back( boost::make_shared< AccelerationSettings >(
                                                                      basic_astrodynamics::central_gravity ) );
                accelerationsOfLunarOrbiter[ "Moon" ].push_back( boost::make_shared< AccelerationSettings >(
                                                                     basic_astrodynamics::central_gravity ) );
                accelerationsOfLunarOrbiter[ "Sun" ].push_back( boost::make_shared< AccelerationSettings >(
                                                                    basic_astrodynamics::central_gravity ) );

            }
            else if( accelerationCase == 2 )
            {
                accelerationsOfLunarOrbiter[ "Earth" ].push_back( boost::make_shared< AccelerationSettings >(
                                                                      basic_astrodynamics::central_gravity ) );
                accelerationsOfLunarOrbiter[ "Moon" ].push_back( boost::make_shared< AccelerationSettings >(
                                                                     basic_astrodynamics::central_gravity ) );
                accelerationsOfLunarOrbiter[ "Sun" ].push_back( boost::make_shared< AccelerationSettings >(
                                                                    basic_astrodynamics::central_gravity ) );
                accelerationsOfLunarOrbiter[ "Jupiter" ].push_back( boost::make_shared< AccelerationSettings >(
                                                                        basic_astrodynamics::central_gravity ) );
                accelerationsOfLunarOrbiter[ "Mars" ].push_back( boost::make_shared< AccelerationSettings >(
                                                                     basic_astrodynamics::central_gravity ) );
                accelerationsOfLunarOrbiter[ "Venus" ].push_back( boost::make_shared< AccelerationSettings >(
                                                                      basic_astrodynamics::central_gravity ) );
            }
            accelerationMap[  "LunarOrbiter" ] = accelerationsOfLunarOrbiter;


            // Create acceleration models and propagation settings.
            basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                        bodyMap, accelerationMap, bodiesToPropagate, centralBodies );

            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            ///////////////////////             CREATE PROPAGATION SETTINGS            ////////////////////////////////////////////
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


            // Set initial conditions for the LunarOrbiter satellite that will be propagated in this simulation.
            // The initial conditions are given in Keplerian elements and later on converted to Cartesian
            // elements.

            // Set Keplerian elements for LunarOrbiter.
            Eigen::Vector6d lunarOrbiterInitialStateInKeplerianElements;
            lunarOrbiterInitialStateInKeplerianElements( semiMajorAxisIndex ) = 2000.0E3;
            lunarOrbiterInitialStateInKeplerianElements( eccentricityIndex ) = 0.05;
            lunarOrbiterInitialStateInKeplerianElements( inclinationIndex ) = convertDegreesToRadians( 85.3 );
            lunarOrbiterInitialStateInKeplerianElements( argumentOfPeriapsisIndex )
                    = convertDegreesToRadians( 235.7 );
            lunarOrbiterInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex )
                    = convertDegreesToRadians( 23.4 );
            lunarOrbiterInitialStateInKeplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 139.87 );

            // Convert LunarOrbiter state from Keplerian elements to Cartesian elements.
            double moonGravitationalParameter = bodyMap.at( "Moon" )->getGravityFieldModel( )->getGravitationalParameter( );
            Eigen::VectorXd systemInitialState = convertKeplerianToCartesianElements(
                        lunarOrbiterInitialStateInKeplerianElements,
                        moonGravitationalParameter ) + spice_interface::getBodyCartesianStateAtEpoch(
                        "Moon", centralBodies.at( 0 ), "ECLIPJ2000", "None", 0.0 );

            // Set simulation end epoch.
            const double simulationStartEpoch = 0.0;
            const double simulationEndEpoch = 14.0 * tudat::physical_constants::JULIAN_DAY;

            boost::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
                    boost::make_shared< TranslationalStatePropagatorSettings< double > >
                    ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState, simulationEndEpoch, cowell );

            boost::shared_ptr< IntegratorSettings< > > integratorSettings = boost::make_shared< IntegratorSettings< > >
                    ( rungeKutta4, simulationStartEpoch, 10.0 );

            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            // Create simulation object and propagate dynamics.
            SingleArcDynamicsSimulator< > dynamicsSimulator(
                        bodyMap, integratorSettings, propagatorSettings );
            std::map< double, Eigen::VectorXd > integrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
            std::map< double, Eigen::VectorXd > integrationResultToPrint;

            for( std::map< double, Eigen::VectorXd >::const_iterator stateIterator = integrationResult.begin( );
                 stateIterator != integrationResult.end( ); stateIterator++ )
            {
                if( originCase == 0 )
                {
                    integrationResultToPrint[ stateIterator->first ] = stateIterator->second;
                }
                else if( originCase == 1 )
                {
                    integrationResultToPrint[ stateIterator->first ] =
                            stateIterator->second - spice_interface::getBodyCartesianStateAtEpoch(
                                "Moon", "Earth", "ECLIPJ2000", "None", stateIterator->first );
                }
                else if( originCase == 2 )
                {
                    integrationResultToPrint[ stateIterator->first ] =
                            stateIterator->second - spice_interface::getBodyCartesianStateAtEpoch(
                                "Moon", "Sun", "ECLIPJ2000", "None", stateIterator->first );
                }
                else if( originCase == 3 )
                {
                    integrationResultToPrint[ stateIterator->first ] =
                            stateIterator->second - spice_interface::getBodyCartesianStateAtEpoch(
                                "Moon", "SSB", "ECLIPJ2000", "None", stateIterator->first );
                }
            }

            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            ///////////////////////        PROVIDE OUTPUT TO CONSOLE AND FILES           //////////////////////////////////////////
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


            // Write satellite propagation history to file.
            input_output::writeDataMapToTextFile( integrationResultToPrint,
                                                  "stateMoonOrbiterOriginCases_" +
                                                  boost::lexical_cast< std::string >( originCase ) + "_" +
                                                  boost::lexical_cast< std::string >( accelerationCase ) +
                                                  ".dat",
                                                  outputDirectory, "",
                                                  std::numeric_limits< double >::digits10,
                                                  std::numeric_limits< double >::digits10,
                                                  "," );

            // Write satellite propagation history to file.
            input_output::writeDataMapToTextFile( integrationResult,
                                                  "stateMoonOrbiterOriginCases_full_" +
                                                  boost::lexical_cast< std::string >( originCase ) + "_" +
                                                  boost::lexical_cast< std::string >( accelerationCase ) +
                                                  ".dat",
                                                  outputDirectory, "",
                                                  std::numeric_limits< double >::digits10,
                                                  std::numeric_limits< double >::digits10,
                                                  "," );

            // Final statement.
            // The exit code EXIT_SUCCESS indicates that the program was successfully executed.

        }
    }
    return EXIT_SUCCESS;

}
