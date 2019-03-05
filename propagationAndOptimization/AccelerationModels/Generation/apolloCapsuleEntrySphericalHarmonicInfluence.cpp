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

#include <Tudat/Astrodynamics/Aerodynamics/UnitTests/testApolloCapsuleCoefficients.h>

#include "propagationAndOptimization/applicationOutput.h"


//! Execute propagation of orbits of Apollo during entry.
int main( )
{
    std::string outputPath = tudat_applications::getOutputPath( "AccelerationModels/" );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            USING STATEMENTS              //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    using namespace tudat::ephemerides;
    using namespace tudat::interpolators;
    using namespace tudat::numerical_integrators;
    using namespace tudat::spice_interface;
    using namespace tudat::simulation_setup;
    using namespace tudat::basic_astrodynamics;
    using namespace tudat::orbital_element_conversions;
    using namespace tudat::propagators;
    using namespace tudat::aerodynamics;
    using namespace tudat::basic_mathematics;
    using namespace tudat::input_output;
    using namespace tudat;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE ENVIRONMENT            //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Load Spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Set simulation start epoch.
    const double simulationStartEpoch = 0.0;

    // Set simulation end epoch.
    const double simulationEndEpoch = 3600.0;

    // Set numerical integration fixed step size.
    const double fixedStepSize = 0.1;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE ENVIRONMENT            //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define simulation body settings.
    std::map< std::string, std::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( { "Earth" }, simulationStartEpoch - 10.0 * fixedStepSize,
                                    simulationEndEpoch + 10.0 * fixedStepSize );
    bodySettings[ "Earth" ]->ephemerisSettings = std::make_shared< simulation_setup::ConstantEphemerisSettings >(
                Eigen::Vector6d::Zero( ), "SSB", "J2000" );
    bodySettings[ "Earth" ]->rotationModelSettings->resetOriginalFrame( "J2000" );

    // Create Earth object
    simulation_setup::NamedBodyMap bodyMap = simulation_setup::createBodies( bodySettings );

    // Create vehicle objects.
    bodyMap[ "Apollo" ] = std::make_shared< simulation_setup::Body >( );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE VEHICLE            /////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    // Create vehicle aerodynamic coefficients
    bodyMap[ "Apollo" ]->setAerodynamicCoefficientInterface(
                unit_tests::getApolloCoefficientInterface( ) );
    bodyMap[ "Apollo" ]->setConstantBodyMass( 5.0E3 );

    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "J2000" );

    double angleStep = mathematical_constants::PI / 10.0;

    for( unsigned int simulationCase = 0; simulationCase < 1; simulationCase++ )
    {
        int latitudeCase = 0;
        for( double initialLatitude = -mathematical_constants::PI/2.0 + angleStep / 2.0;
             initialLatitude <= mathematical_constants::PI/2.0;
             initialLatitude += angleStep )
        {
            int longitudeCase = 0;
            for( double initialLongitude = 0.0; initialLongitude < 2.0 * mathematical_constants::PI;
                 initialLongitude += angleStep )
            {

                int headingCase = 0;
                for( double initialHeadingAngle = 0.0; initialHeadingAngle <= mathematical_constants::PI / 2.0 + 0.001;
                     initialHeadingAngle += mathematical_constants::PI / 4.0 )
                {
                    std::cout<<std::setprecision( 16 )<<simulationCase<<" "<<initialLatitude<<" "<<initialLongitude<<" "<<initialHeadingAngle<<std::endl;
                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    ///////////////////////             CREATE ACCELERATIONS            ///////////////////////////////////////////////////
                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                    // Define propagator settings variables.
                    SelectedAccelerationMap accelerationMap;
                    std::vector< std::string > bodiesToPropagate;
                    std::vector< std::string > centralBodies;

                    // Define acceleration model settings.
                    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfApollo;
                    if( simulationCase == 0 )
                    {
                        accelerationsOfApollo[ "Earth" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >( 32, 32 ) );
                    }
                    else if( simulationCase == 1 )
                    {
                        accelerationsOfApollo[ "Earth" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >( 8, 8 ) );
                    }
                    else if( simulationCase == 2 )
                    {
                        accelerationsOfApollo[ "Earth" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >( 8, 8 ) );
                    }
                    else if( simulationCase == 3 )
                    {
                        accelerationsOfApollo[ "Earth" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >( 4, 4 ) );
                    }
                    else if( simulationCase == 4 )
                    {
                        accelerationsOfApollo[ "Earth" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >( 4, 0 ) );
                    }
                    else if( simulationCase == 5 )
                    {
                        accelerationsOfApollo[ "Earth" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >( 3, 0 ) );
                    }
                    else if( simulationCase == 6 )
                    {
                        accelerationsOfApollo[ "Earth" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >( 2, 0 ) );
                    }
                    else if( simulationCase == 7 )
                    {
                        accelerationsOfApollo[ "Earth" ].push_back( std::make_shared< AccelerationSettings >( central_gravity ) );
                    }
                    accelerationsOfApollo[ "Earth" ].push_back( std::make_shared< AccelerationSettings >( aerodynamic ) );
                    accelerationMap[  "Apollo" ] = accelerationsOfApollo;

                    bodiesToPropagate.push_back( "Apollo" );
                    centralBodies.push_back( "Earth" );

                    // Create acceleration models
                    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                                bodyMap, accelerationMap, bodiesToPropagate, centralBodies );

                    // Define constant 30 degree angle of attack
                    double constantAngleOfAttack = 0.0 * mathematical_constants::PI / 180.0;
                    bodyMap.at( "Apollo" )->getFlightConditions( )->getAerodynamicAngleCalculator( )->setOrientationAngleFunctions(
                                boost::lambda::constant( constantAngleOfAttack ) );

                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    ///////////////////////             CREATE PROPAGATION SETTINGS            ////////////////////////////////////////////
                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                    // Set spherical elements for Apollo.
                    Eigen::Vector6d apolloSphericalEntryState;
                    apolloSphericalEntryState( SphericalOrbitalStateElementIndices::radiusIndex ) =
                            spice_interface::getAverageRadius( "Earth" ) + 120.0E3;
                    apolloSphericalEntryState( SphericalOrbitalStateElementIndices::latitudeIndex ) = initialLatitude;
                    apolloSphericalEntryState( SphericalOrbitalStateElementIndices::longitudeIndex ) = longitudeCase;
                    apolloSphericalEntryState( SphericalOrbitalStateElementIndices::speedIndex ) = 7.65E3;
                    apolloSphericalEntryState( SphericalOrbitalStateElementIndices::flightPathIndex ) =
                            -1.25 * mathematical_constants::PI / 180.0;
                    apolloSphericalEntryState( SphericalOrbitalStateElementIndices::headingAngleIndex ) = initialHeadingAngle;

                    // Convert apollo state from spherical elements to Cartesian elements.
                    Eigen::Vector6d systemInitialState = convertSphericalOrbitalToCartesianState(
                                apolloSphericalEntryState );

                    std::shared_ptr< ephemerides::RotationalEphemeris > earthRotationalEphemeris =
                            bodyMap.at( "Earth" )->getRotationalEphemeris( );
                    systemInitialState = transformStateToGlobalFrame( systemInitialState, simulationStartEpoch, earthRotationalEphemeris );

                    // Define list of dependent variables to save.
                    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariablesList;
                    dependentVariablesList.push_back(
                                std::make_shared< SingleDependentVariableSaveSettings >( mach_number_dependent_variable, "Apollo" ) );
                    dependentVariablesList.push_back(
                                std::make_shared< SingleDependentVariableSaveSettings >(
                                    altitude_dependent_variable, "Apollo", "Earth" ) );
                    dependentVariablesList.push_back(
                                std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                                    aerodynamic, "Apollo", "Earth", 1 ) );
                    dependentVariablesList.push_back(
                                std::make_shared< SingleDependentVariableSaveSettings >(
                                    aerodynamic_force_coefficients_dependent_variable, "Apollo" ) );

                    // Create object with list of dependent variables
                    std::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave =
                            std::make_shared< DependentVariableSaveSettings >( dependentVariablesList );

                    // Define termination conditions
                    std::shared_ptr< SingleDependentVariableSaveSettings > terminationDependentVariable =
                            std::make_shared< SingleDependentVariableSaveSettings >(
                                altitude_dependent_variable, "Apollo", "Earth" );

                    std::vector< std::shared_ptr< PropagationTerminationSettings > > terminationSettingsList;
                    terminationSettingsList.push_back(
                                std::make_shared< PropagationDependentVariableTerminationSettings >(
                                    terminationDependentVariable, 25.0E3, true ) );
                    terminationSettingsList.push_back(
                                std::make_shared< PropagationTimeTerminationSettings >(
                                    simulationEndEpoch  ) );
                    std::shared_ptr< PropagationTerminationSettings > totalTerminationSettings =
                            std::make_shared< PropagationHybridTerminationSettings >(
                                terminationSettingsList, true );


                    // Create propagation settings.
                    std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
                            std::make_shared< TranslationalStatePropagatorSettings< double > >
                            ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState,
                              totalTerminationSettings, cowell, dependentVariablesToSave );
                    std::shared_ptr< IntegratorSettings< > > integratorSettings =
                            std::make_shared< RungeKuttaVariableStepSizeSettings< > >
                            ( rungeKuttaVariableStepSize, 0.0, fixedStepSize,
                              RungeKuttaCoefficients::rungeKuttaFehlberg78, 1.0E-6, 3600.0, 1.0E-14, 1.0E-14 );
                    //std::make_shared< IntegratorSettings< > >
                    //( rungeKutta4, simulationStartEpoch, fixedStepSize );


                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


                    // Create simulation object and propagate dynamics.
                    SingleArcDynamicsSimulator< > dynamicsSimulator(
                                bodyMap, integratorSettings, propagatorSettings );

                    std::map< double, Eigen::VectorXd > stateHistoryInertialFrame =
                            dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
                    std::map< double, Eigen::VectorXd > stateHistoryEarthFixedFrame;

                    for( std::map< double, Eigen::VectorXd >::const_iterator stateIterator = stateHistoryInertialFrame.begin( );
                         stateIterator != stateHistoryInertialFrame.end( ); stateIterator++ )
                    {
                        stateHistoryEarthFixedFrame[ stateIterator->first ] =
                                transformStateToTargetFrame(
                                    Eigen::Vector6d( stateIterator->second ), stateIterator->first, earthRotationalEphemeris );
                    }

                    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::VectorXd > > inertialInterpolator =
                            std::make_shared< interpolators::LagrangeInterpolator< double, Eigen::VectorXd > >(
                                stateHistoryInertialFrame, 8 );
                    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::VectorXd > > bodyFixedInterpolator =
                            std::make_shared< interpolators::LagrangeInterpolator< double, Eigen::VectorXd > >(
                                stateHistoryEarthFixedFrame, 8 );


                    std::map< double, Eigen::VectorXd > interpolatedStateHistoryInertialFrame;
                    std::map< double, Eigen::VectorXd > interpolatedStateHistoryEarthFixedFrame;

                    double currentTime = 0;
                    double endTime = stateHistoryInertialFrame.rbegin( )->first;
                    while( currentTime < endTime )
                    {
                        interpolatedStateHistoryInertialFrame[ currentTime ] = inertialInterpolator->interpolate( currentTime );
                        interpolatedStateHistoryEarthFixedFrame[ currentTime ] = bodyFixedInterpolator->interpolate( currentTime );
                        currentTime += 1.0;
                    }

                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    ///////////////////////        PROVIDE OUTPUT TO FILE                        //////////////////////////////////////////
                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


                    // Write Apollo propagation history to file.
                    writeDataMapToTextFile( interpolatedStateHistoryInertialFrame,
                                            "stateReEntrySphericalHarmonicCases_" +
                                            boost::lexical_cast< std::string >( latitudeCase ) + "_" +
                                            boost::lexical_cast< std::string >( longitudeCase ) + "_" +
                                            boost::lexical_cast< std::string >( headingCase ) + "_" +
                                            boost::lexical_cast< std::string >( simulationCase ) +".dat",
                                            outputPath,
                                            "",
                                            std::numeric_limits< double >::digits10,
                                            std::numeric_limits< double >::digits10,
                                            "," );
                    writeDataMapToTextFile( interpolatedStateHistoryEarthFixedFrame,
                                            "bodyFixedStateReEntrySphericalHarmonicCases_" +
                                            boost::lexical_cast< std::string >( latitudeCase ) + "_" +
                                            boost::lexical_cast< std::string >( longitudeCase ) + "_" +
                                            boost::lexical_cast< std::string >( headingCase ) + "_" +
                                            boost::lexical_cast< std::string >( simulationCase ) +".dat",
                                            outputPath,
                                            "",
                                            std::numeric_limits< double >::digits10,
                                            std::numeric_limits< double >::digits10,
                                            "," );
                    writeDataMapToTextFile( dynamicsSimulator.getDependentVariableHistory( ),
                                            "dependentVariablesReEntrySphericalHarmonicCases_" +
                                            boost::lexical_cast< std::string >( latitudeCase ) + "_" +
                                            boost::lexical_cast< std::string >( longitudeCase ) + "_" +
                                            boost::lexical_cast< std::string >( headingCase ) + "_" +
                                            boost::lexical_cast< std::string >( simulationCase ) +".dat",
                                            outputPath,
                                            "",
                                            std::numeric_limits< double >::digits10,
                                            std::numeric_limits< double >::digits10,
                                            "," );
                    headingCase++;
                }
                longitudeCase++;
            }
            latitudeCase++;
        }
    }
    // Final statement.
    // The exit code EXIT_SUCCESS indicates that the program was successfully executed.
    return EXIT_SUCCESS;

}
