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
    spice_interface::loadStandardSpiceKernels( );

    for( unsigned int ephemerisCase = 0; ephemerisCase < 3; ephemerisCase++ )
    {

        // Create body objects.
        std::vector< std::string > bodiesToCreate;
        bodiesToCreate.push_back( "Earth" );
        bodiesToCreate.push_back( "Moon" );
        bodiesToCreate.push_back( "Sun" );
        bodiesToCreate.push_back( "Jupiter" );
        bodiesToCreate.push_back( "Mars" );
        bodiesToCreate.push_back( "Venus" );

        std::map< std::string, std::shared_ptr< BodySettings > > bodySettings =
                getDefaultBodySettings( bodiesToCreate, -3600.0, 14.0 * tudat::physical_constants::JULIAN_DAY + 3600.0 );


        for( std::map< std::string, std::shared_ptr< BodySettings > >::iterator settingsIterator = bodySettings.begin( );
             settingsIterator != bodySettings.end( ); settingsIterator++ )
        {
            settingsIterator->second->ephemerisSettings->resetFrameOrientation( "ECLIPJ2000" );
            settingsIterator->second->rotationModelSettings->resetOriginalFrame( "ECLIPJ2000" );
            if( ephemerisCase > 0 )
            {
                if( settingsIterator->first != "Sun" )
                {
                    std::shared_ptr< EphemerisSettings > newEphemerisSettings =
                            settingsIterator->second->ephemerisSettings;
                    if( ephemerisCase == 1 )
                    {
                        ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData ephemerisBodyId;
                        if( settingsIterator->first == "Jupiter" )
                        {
                            ephemerisBodyId = ApproximatePlanetPositionsBase::jupiter;
                            newEphemerisSettings = std::make_shared < ApproximatePlanetPositionSettings >(
                                        ephemerisBodyId, false );
                        }
                        else if( settingsIterator->first == "Mars" )
                        {
                            ephemerisBodyId = ApproximatePlanetPositionsBase::jupiter;
                            newEphemerisSettings = std::make_shared < ApproximatePlanetPositionSettings >(
                                        ephemerisBodyId, false );
                        }
                        else if( settingsIterator->first == "Venus" )
                        {
                            ephemerisBodyId = ApproximatePlanetPositionsBase::jupiter;
                            newEphemerisSettings = std::make_shared < ApproximatePlanetPositionSettings >(
                                        ephemerisBodyId, false );
                        }
                    }
                    else if( ephemerisCase == 2 )
                    {
                        if( settingsIterator->first != "Sun" )
                        {
                            double gravitationalParameter;
                            Eigen::Vector6d initialCartesianState;
                            std::string referenceBody;
                            std::string ephemerisBody = settingsIterator->first;
                            if( ephemerisBody == "Mars" || ephemerisBody == "Jupiter" )
                            {
                                ephemerisBody += " Barycenter";
                            }

                            if( settingsIterator->first == "Moon" )
                            {
                                referenceBody = "Earth";
                            }
                            else
                            {
                                referenceBody = "Sun";
                            }

                            gravitationalParameter = getBodyGravitationalParameter( referenceBody ) +
                                    getBodyGravitationalParameter( settingsIterator->first );
                            initialCartesianState = getBodyCartesianStateAtEpoch(
                                        ephemerisBody, referenceBody, "ECLIPJ2000", "NONE", 0.0 );
                            newEphemerisSettings = std::make_shared< KeplerEphemerisSettings >(
                                        convertCartesianToKeplerianElements( initialCartesianState, gravitationalParameter ),
                                        0.0, gravitationalParameter, referenceBody, "ECLIPJ2000" );
                        }
//                        else if( ephemerisCase == 3 )
//                        {
//                            std::string centralBody;
//                            std::string ephemerisBody = settingsIterator->first;

//                            if( settingsIterator->first == "Sun" )
//                            {
//                                centralBody = "SSB";
//                            }
//                            else if( settingsIterator->first == "Moon" )
//                            {
//                                centralBody = "Earth";
//                            }
//                            else
//                            {
//                                centralBody = "Sun";
//                            }

//                            if( ephemerisBody == "Mars" || ephemerisBody == "Jupiter" )
//                            {
//                                ephemerisBody += " Barycenter";
//                            }

//                            newEphemerisSettings = std::make_shared< EphemerisSettings >(
//                                        convertCartesianToKeplerianElements( initialCartesianState, gravitationalParameter ),
//                                        0.0, gravitationalParameter, referenceBody, "ECLIPJ2000" );
//                        }
                    }

                    settingsIterator->second->ephemerisSettings = newEphemerisSettings;
                }
            }
        }

        // Create Earth object
        NamedBodyMap bodyMap = createBodies( bodySettings );

        // Create spacecraft object.
        bodyMap[ "LunarOrbiter" ] = std::make_shared< simulation_setup::Body >( );

        // Finalize body creation.
        setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );

        for( unsigned int accelerationCase = 0; accelerationCase < 3; accelerationCase++ )
        {
            std::cout<<"Test case: "<<accelerationCase<<" "<<ephemerisCase<<std::endl;
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            ///////////////////////            CREATE ACCELERATIONS          //////////////////////////////////////////////////////
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            // Define propagator settings variables.
            SelectedAccelerationMap accelerationMap;
            std::vector< std::string > bodiesToPropagate;
            std::vector< std::string > centralBodies;

            bodiesToPropagate.push_back( "LunarOrbiter" );
            centralBodies.push_back( "Moon" );

            // Define propagation settings.
            std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfLunarOrbiter;

            if( accelerationCase == 0 )
            {
                accelerationsOfLunarOrbiter[ "Moon" ].push_back( std::make_shared< AccelerationSettings >(
                                                                     basic_astrodynamics::central_gravity ) );
            }
            else if( accelerationCase == 1 )
            {
                accelerationsOfLunarOrbiter[ "Earth" ].push_back( std::make_shared< AccelerationSettings >(
                                                                      basic_astrodynamics::central_gravity ) );
                accelerationsOfLunarOrbiter[ "Moon" ].push_back( std::make_shared< AccelerationSettings >(
                                                                     basic_astrodynamics::central_gravity ) );
                accelerationsOfLunarOrbiter[ "Sun" ].push_back( std::make_shared< AccelerationSettings >(
                                                                    basic_astrodynamics::central_gravity ) );

            }
            else if( accelerationCase == 2 )
            {
                accelerationsOfLunarOrbiter[ "Earth" ].push_back( std::make_shared< AccelerationSettings >(
                                                                      basic_astrodynamics::central_gravity ) );
                accelerationsOfLunarOrbiter[ "Moon" ].push_back( std::make_shared< AccelerationSettings >(
                                                                     basic_astrodynamics::central_gravity ) );
                accelerationsOfLunarOrbiter[ "Sun" ].push_back( std::make_shared< AccelerationSettings >(
                                                                    basic_astrodynamics::central_gravity ) );
                accelerationsOfLunarOrbiter[ "Jupiter" ].push_back( std::make_shared< AccelerationSettings >(
                                                                        basic_astrodynamics::central_gravity ) );
                accelerationsOfLunarOrbiter[ "Mars" ].push_back( std::make_shared< AccelerationSettings >(
                                                                     basic_astrodynamics::central_gravity ) );
                accelerationsOfLunarOrbiter[ "Venus" ].push_back( std::make_shared< AccelerationSettings >(
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
                        moonGravitationalParameter );

            // Set simulation end epoch.
            const double simulationStartEpoch = 0.0;
            const double simulationEndEpoch = 14.0 * tudat::physical_constants::JULIAN_DAY;

            std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
                    std::make_shared< TranslationalStatePropagatorSettings< double > >
                    ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState, simulationEndEpoch, encke );

            std::shared_ptr< IntegratorSettings< > > integratorSettings = std::make_shared< IntegratorSettings< > >
                    ( rungeKutta4, simulationStartEpoch, 10.0 );

            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            // Create simulation object and propagate dynamics.
            SingleArcDynamicsSimulator< > dynamicsSimulator(
                        bodyMap, integratorSettings, propagatorSettings );
            std::map< double, Eigen::VectorXd > integrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );

            std::map< double, double > orbitalEnergy;

            double currentSpeed;
            for( std::map< double, Eigen::VectorXd >::const_iterator stateIterator = integrationResult.begin( );
                 stateIterator != integrationResult.end( ); stateIterator++ )
            {
                currentSpeed = stateIterator->second.segment( 3, 3 ).norm( );
                orbitalEnergy[ stateIterator->first ] =
                        currentSpeed * currentSpeed / 2.0 - moonGravitationalParameter / stateIterator->second.segment( 0, 3 ).norm( );
            }

            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            ///////////////////////        PROVIDE OUTPUT TO CONSOLE AND FILES           //////////////////////////////////////////
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            // Write satellite propagation history to file.
            input_output::writeDataMapToTextFile( orbitalEnergy,
                                                  "orbitalEnergyMoonOrbiter_EphemerisInfluence_" +
                                                  boost::lexical_cast< std::string >( ephemerisCase ) + "_" +
                                                  boost::lexical_cast< std::string >( accelerationCase ) +
                                                  ".dat",
                                                  outputDirectory,
                                                  "",
                                                  std::numeric_limits< double >::digits10,
                                                  std::numeric_limits< double >::digits10,
                                                  "," );

            // Write satellite propagation history to file.
            input_output::writeDataMapToTextFile( integrationResult,
                                                  "stateMoonOrbiter_EphemerisInfluence_" +
                                                  boost::lexical_cast< std::string >( ephemerisCase ) + "_" +
                                                  boost::lexical_cast< std::string >( accelerationCase ) +
                                                  ".dat",
                                                  outputDirectory,
                                                  "",
                                                  std::numeric_limits< double >::digits10,
                                                  std::numeric_limits< double >::digits10,
                                                  "," );

            // Final statement.
            // The exit code EXIT_SUCCESS indicates that the program was successfully executed.
        }
    }
    return EXIT_SUCCESS;

}
