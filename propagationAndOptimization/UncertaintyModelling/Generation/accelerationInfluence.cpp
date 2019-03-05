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
    std::string outputDirectory = tudat_applications::getOutputPath( "UncertaintyModelling/" );

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            USING STATEMENTS              //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    using namespace tudat;
    using namespace tudat::simulation_setup;
    using namespace tudat::propagators;
    using namespace tudat::numerical_integrators;
    using namespace tudat::orbital_element_conversions;
    using namespace tudat::basic_mathematics;
    using namespace tudat::gravitation;
    using namespace tudat::numerical_integrators;
    using namespace tudat::basic_astrodynamics;


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////     CREATE ENVIRONMENT AND VEHICLE       //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Load Spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Set simulation time settings.
    const double simulationStartEpoch = 0.0;
    const double simulationEndEpoch = 7.0 * tudat::physical_constants::JULIAN_DAY;

    // Define body settings for simulation.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Sun" );
    bodiesToCreate.push_back( "Earth" );
    bodiesToCreate.push_back( "Moon" );
    bodiesToCreate.push_back( "Mars" );
    bodiesToCreate.push_back( "Venus" );
    bodiesToCreate.push_back( "Jupiter" );
    bodiesToCreate.push_back( "Saturn" );

    // Create body objects.
    std::map< std::string, std::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( bodiesToCreate, simulationStartEpoch - 300.0, simulationEndEpoch + 300.0 );
    for( unsigned int i = 0; i < bodiesToCreate.size( ); i++ )
    {
        bodySettings[ bodiesToCreate.at( i ) ]->ephemerisSettings->resetFrameOrientation( "J2000" );
        bodySettings[ bodiesToCreate.at( i ) ]->rotationModelSettings->resetOriginalFrame( "J2000" );
    }
    NamedBodyMap bodyMap = createBodies( bodySettings );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE VEHICLE            /////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create spacecraft object.
    bodyMap[ "Asterix" ] = std::make_shared< simulation_setup::Body >( );
    bodyMap[ "Asterix" ]->setConstantBodyMass( 400.0 );

    // Create aerodynamic coefficient interface settings.
    double referenceArea = 4.0;
    double aerodynamicCoefficient = 1.2;
    std::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings =
            std::make_shared< ConstantAerodynamicCoefficientSettings >(
                referenceArea, aerodynamicCoefficient * Eigen::Vector3d::UnitX( ), 1, 1 );

    // Create and set aerodynamic coefficients object
    bodyMap[ "Asterix" ]->setAerodynamicCoefficientInterface(
                createAerodynamicCoefficientInterface( aerodynamicCoefficientSettings, "Asterix" ) );

    // Create radiation pressure settings
    double referenceAreaRadiation = 4.0;
    double radiationPressureCoefficient = 1.2;
    std::vector< std::string > occultingBodies;
    occultingBodies.push_back( "Earth" );
    std::shared_ptr< RadiationPressureInterfaceSettings > asterixRadiationPressureSettings =
            std::make_shared< CannonBallRadiationPressureInterfaceSettings >(
                "Sun", referenceAreaRadiation, radiationPressureCoefficient, occultingBodies );

    // Create and set radiation pressure settings
    bodyMap[ "Asterix" ]->setRadiationPressureInterface(
                "Sun", createRadiationPressureInterface(
                    asterixRadiationPressureSettings, "Asterix", bodyMap ) );

    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "J2000" );

    for(  unsigned int integratorCase = 0; integratorCase < 2; integratorCase++ )
    {
        for( unsigned int accelerationCase = 0; accelerationCase < 8; accelerationCase++ )
        {
            std::cout<<"Acc case "<<accelerationCase<<std::endl;
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            ///////////////////////            CREATE ACCELERATIONS          //////////////////////////////////////////////////////
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            // Define propagator settings variables.
            SelectedAccelerationMap accelerationMap;
            std::vector< std::string > bodiesToPropagate;
            std::vector< std::string > centralBodies;

            // Define propagation settings.
            std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfAsterix;

            if( accelerationCase == 0 )
            {
                accelerationsOfAsterix[ "Earth" ].push_back( std::make_shared< AccelerationSettings >(
                                                                 central_gravity ) );
            }
            else
            {
                if( accelerationCase <= 6 )
                {
                    accelerationsOfAsterix[ "Earth" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >( 8, 8 ) );
                }
            }

            if( accelerationCase > 1 )
            {
                accelerationsOfAsterix[ "Sun" ].push_back( std::make_shared< AccelerationSettings >(
                                                               central_gravity ) );
                if( accelerationCase <= 5 )
                {
                    accelerationsOfAsterix[ "Moon" ].push_back( std::make_shared< AccelerationSettings >(
                                                                    central_gravity ) );
                }
            }

            if( accelerationCase > 2 )
            {
                accelerationsOfAsterix[ "Earth" ].push_back( std::make_shared< AccelerationSettings >(
                                                                 aerodynamic ) );
            }

            if( accelerationCase > 3 )
            {
                accelerationsOfAsterix[ "Sun" ].push_back( std::make_shared< AccelerationSettings >(
                                                               cannon_ball_radiation_pressure ) );
            }

            if( accelerationCase > 4 )
            {
                accelerationsOfAsterix[ "Mars" ].push_back( std::make_shared< AccelerationSettings >(
                                                                central_gravity ) );
                accelerationsOfAsterix[ "Venus" ].push_back( std::make_shared< AccelerationSettings >(
                                                                 central_gravity ) );
                accelerationsOfAsterix[ "Saturn" ].push_back( std::make_shared< AccelerationSettings >(
                                                                  central_gravity ) );
                accelerationsOfAsterix[ "Jupiter" ].push_back( std::make_shared< AccelerationSettings >(
                                                                 central_gravity ) );
            }

            if( accelerationCase > 5 )
            {
                accelerationsOfAsterix[ "Moon" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >( 2, 2 ) );
            }

            if( accelerationCase > 6 )
            {
                accelerationsOfAsterix[ "Earth" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >( 32, 32 ) );
            }



            accelerationMap[ "Asterix" ] = accelerationsOfAsterix;
            bodiesToPropagate.push_back( "Asterix" );
            centralBodies.push_back( "Earth" );

            AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                        bodyMap, accelerationMap, bodiesToPropagate, centralBodies );

            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            ///////////////////////             CREATE PROPAGATION SETTINGS            ////////////////////////////////////////////
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            // Set Keplerian elements for Asterix.
            Eigen::Vector6d asterixInitialStateInKeplerianElements;
            asterixInitialStateInKeplerianElements( semiMajorAxisIndex ) = 6800.0E3;
            asterixInitialStateInKeplerianElements( eccentricityIndex ) = 0.01;
            asterixInitialStateInKeplerianElements( inclinationIndex ) = unit_conversions::convertDegreesToRadians( 85.3 );
            asterixInitialStateInKeplerianElements( argumentOfPeriapsisIndex ) = unit_conversions::convertDegreesToRadians( 235.7 );
            asterixInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex ) = unit_conversions::convertDegreesToRadians( 23.4 );
            asterixInitialStateInKeplerianElements( trueAnomalyIndex ) = unit_conversions::convertDegreesToRadians( 139.87 );

            double earthGravitationalParameter = bodyMap.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );
            const Eigen::Vector6d asterixInitialState = convertKeplerianToCartesianElements(
                        asterixInitialStateInKeplerianElements, earthGravitationalParameter );



            TranslationalPropagatorType propagator = cowell;


            std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariables;
            //            dependentVariables.push_back(
            //                        std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
            //                            third_body_central_gravity, "Asterix", "Sun" ) );
            //            dependentVariables.push_back(
            //                        std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
            //                            third_body_central_gravity, "Asterix", "Mars" ) );
            //            dependentVariables.push_back(
            //                        std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
            //                            third_body_central_gravity, "Asterix", "Venus" ) );
            //            dependentVariables.push_back(
            //                        std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
            //                            third_body_central_gravity, "Asterix", "Jupiter" ) );
            //            dependentVariables.push_back(
            //                        std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
            //                            third_body_central_gravity, "Asterix", "Saturn" ) );
            //            dependentVariables.push_back(
            //                        std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
            //                            cannon_ball_radiation_pressure, "Asterix", "Sun" ) );
            //            dependentVariables.push_back(
            //                        std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
            //                            aerodynamic, "Asterix", "Earth" ) );
            //            dependentVariables.push_back(
            //                        std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
            //                            cannon_ball_radiation_pressure, "Asterix", "Sun" ) );

            std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
                    std::make_shared< TranslationalStatePropagatorSettings< double > >
                    ( centralBodies, accelerationModelMap, bodiesToPropagate, asterixInitialState, simulationEndEpoch, propagator,
                      std::make_shared< DependentVariableSaveSettings >( dependentVariables ) );

            const double fixedStepSize = ( integratorCase == 0 ) ? 10.0 : 1.0;
            std::shared_ptr< IntegratorSettings< > > integratorSettings =
                    std::make_shared< IntegratorSettings< > >
                    ( rungeKutta4, 0.0, fixedStepSize );

            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


            // Create simulation object and propagate dynamics.
            SingleArcDynamicsSimulator< > dynamicsSimulator( bodyMap, integratorSettings, propagatorSettings );
            std::map< double, Eigen::VectorXd > integrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );

            input_output::writeDataMapToTextFile(
                        integrationResult, "accelerationModelInfluence" +
                        std::to_string( accelerationCase ) + "_" +
                        std::to_string( integratorCase ) +
                        ".dat", outputDirectory );

        }
    }
    // Final statement.
    // The exit code EXIT_SUCCESS indicates that the program was successfully executed.
    return EXIT_SUCCESS;


}
