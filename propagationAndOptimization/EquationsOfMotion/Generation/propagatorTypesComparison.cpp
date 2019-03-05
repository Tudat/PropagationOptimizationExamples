/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References:
 *      Mooij, E. "Orbit-State Model Selection for Solar-Sailing Mission Optimization."
 *          AIAA/AAS Astrodynamics Specialist Conference. 2012.
 */

#include <Tudat/SimulationSetup/tudatSimulationHeader.h>

#include "Tudat/Mathematics/Statistics/basicStatistics.h"
#include "Tudat/InputOutput/basicInputOutput.h"
#include "Tudat/Basics/utilities.h"

#include "propagationAndOptimization/applicationOutput.h"

//! Execute propagation of orbit of Satellite around the Earth.
int main( )
{
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            USING STATEMENTS              //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    using namespace tudat;
    using namespace tudat::aerodynamics;
    using namespace tudat::basic_astrodynamics;
    using namespace tudat::basic_mathematics;
    using namespace tudat::ephemerides;
    using namespace tudat::estimatable_parameters;
    using namespace tudat::gravitation;
    using namespace tudat::input_output;
    using namespace tudat::numerical_integrators;
    using namespace tudat::orbital_element_conversions;
    using namespace tudat::propagators;
    using namespace tudat::simulation_setup;
    using namespace tudat::unit_conversions;

    using namespace tudat_applications;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            DEFINE TEST CASES             //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Predefine variables
    bool keplerOrbit = true; // toggle use of accelerations (true == no accelerations)
    std::string keplerSuffix = ( keplerOrbit == false ) ? "" : "_Kepler";
    double simulationDuration = 10.0 * physical_constants::JULIAN_DAY;
    double integrationRelativeTolerance = 1.0e-10;
    double integrationAbsoluteTolerance = 1.0e-10;
    double integrationReferenceTolerance = 1.0e-18;
    double integrationConstantTimeStepSize = 10.0;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////     CREATE ENVIRONMENT AND VEHICLE       //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Load Spice kernels
    spice_interface::loadStandardSpiceKernels( );

    // Set simulation time settings
    const double simulationStartEpoch = 7.0 * physical_constants::JULIAN_YEAR + 30.0 * 6.0 * physical_constants::JULIAN_DAY;
    const double simulationEndEpoch = simulationDuration + simulationStartEpoch;

    // Define body settings for simulation
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Earth" );
    bodiesToCreate.push_back( "Sun" );
    bodiesToCreate.push_back( "Moon" );

    // Create body objects
    std::map< std::string, std::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( bodiesToCreate, simulationStartEpoch - 1.0e3, simulationEndEpoch + 1.0e3 );

    for ( unsigned int i = 0; i < bodiesToCreate.size( ); i++ )
    {
        bodySettings[ bodiesToCreate.at( i ) ]->ephemerisSettings->resetFrameOrientation( "J2000" );
        bodySettings[ bodiesToCreate.at( i ) ]->rotationModelSettings->resetOriginalFrame( "J2000" );
    }

    bodySettings[ "Earth" ]->gravityFieldSettings = std::make_shared< FromFileSphericalHarmonicsGravityFieldSettings >( ggm02s );
    bodySettings[ "Earth" ]->atmosphereSettings = std::make_shared< ExponentialAtmosphereSettings >( aerodynamics::earth );
    NamedBodyMap bodyMap = createBodies( bodySettings );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE VEHICLE            /////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create spacecraft object
    bodyMap[ "Satellite" ] = std::make_shared< Body >( );
    const double satelliteMass = 1000.0;
    bodyMap[ "Satellite" ]->setConstantBodyMass( satelliteMass );

    // Set constant aerodynamic drag coefficient
    const double referenceAreaAerodynamic = 37.5;
    const Eigen::Vector3d aerodynamicCoefficients = 2.2 * Eigen::Vector3d::UnitX( ); // only drag coefficient
    std::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings =
            std::make_shared< ConstantAerodynamicCoefficientSettings >( referenceAreaAerodynamic, aerodynamicCoefficients, true, true );

    bodyMap[ "Satellite" ]->setAerodynamicCoefficientInterface(
                createAerodynamicCoefficientInterface( aerodynamicCoefficientSettings, "Satellite" ) );

    // Create radiation pressure settings
    double referenceAreaRadiation = referenceAreaAerodynamic;
    double radiationPressureCoefficient = 1.25;
    std::vector< std::string > occultingBodies;
    occultingBodies.push_back( "Earth" );
    std::shared_ptr< RadiationPressureInterfaceSettings > SatelliteRadiationPressureSettings =
            std::make_shared< CannonBallRadiationPressureInterfaceSettings >(
                "Sun", referenceAreaRadiation, radiationPressureCoefficient, occultingBodies );

    // Create and set radiation pressure settings
    bodyMap[ "Satellite" ]->setRadiationPressureInterface( "Sun", createRadiationPressureInterface(
                                                               SatelliteRadiationPressureSettings, "Satellite", bodyMap ) );

    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "J2000" );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE ACCELERATIONS          //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define propagator settings variables
    SelectedAccelerationMap accelerationMap;
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;

    // Switch between Kepler orbit or perturbed environment
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfSatellite;
    if ( keplerOrbit )
    {
        // Only central gravity
        accelerationsOfSatellite[ "Earth" ].push_back( std::make_shared< AccelerationSettings >( central_gravity ) );
    }
    else
    {
        // Define spherical harmonics, third bodies, solar radiation and aerodynamic forces
        accelerationsOfSatellite[ "Earth" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >( 8, 8 ) );
        for ( unsigned int i = 0; i < bodiesToCreate.size( ); i++ )
        {
            if ( bodiesToCreate.at( i ) != "Earth" )
            {
                accelerationsOfSatellite[ bodiesToCreate.at( i ) ].push_back( std::make_shared< AccelerationSettings >( central_gravity ) );
            }
        }
        accelerationsOfSatellite[ "Sun" ].push_back( std::make_shared< AccelerationSettings >( cannon_ball_radiation_pressure ) );
        accelerationsOfSatellite[ "Earth" ].push_back( std::make_shared< AccelerationSettings >( aerodynamic ) );
    }

    // Add acceleration information
    accelerationMap[ "Satellite" ] = accelerationsOfSatellite;
    bodiesToPropagate.push_back( "Satellite" );
    centralBodies.push_back( "Earth" );

    AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap, accelerationMap, bodiesToPropagate, centralBodies );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE INITIAL CONDITIONS              ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    std::vector< double > eccentricityList = { 0.1, 0.3, 0.6, 0.9, 0.99 };
    for ( unsigned int eccentricityType = 0; eccentricityType < 4; eccentricityType++ )
    {
        std::cout << "Eccentricity: " << eccentricityList.at( eccentricityType ) << std::endl;

        std::string fileSuffix = keplerSuffix + "_e" +
                std::to_string( eccentricityList.at( eccentricityType ) );

        // Set Keplerian initial conditions
        Eigen::Vector6d satelliteInitialStateInKeplerianElements;
        satelliteInitialStateInKeplerianElements( semiMajorAxisIndex ) = 7500.0E3;
        satelliteInitialStateInKeplerianElements( eccentricityIndex ) = eccentricityList.at( 0 );
        satelliteInitialStateInKeplerianElements( inclinationIndex ) = convertDegreesToRadians( 93.0 );
        satelliteInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex ) = convertDegreesToRadians( 158.7 );
        satelliteInitialStateInKeplerianElements( argumentOfPeriapsisIndex ) = convertDegreesToRadians( 23.4 );
        satelliteInitialStateInKeplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 0.0 );

        if( eccentricityType != 0 )
        {
            double newEccentricity = eccentricityList.at( eccentricityType );
            satelliteInitialStateInKeplerianElements( semiMajorAxisIndex ) *=
                    ( 1.0 - satelliteInitialStateInKeplerianElements( eccentricityIndex ) ) /
                    ( 1.0 - newEccentricity );
            satelliteInitialStateInKeplerianElements( eccentricityIndex ) = newEccentricity;
        }
        // Convert to Cartesian elements
        double mainGravitationalParameter = bodyMap.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );
        const Eigen::Vector6d satelliteInitialState = convertKeplerianToCartesianElements(
                    satelliteInitialStateInKeplerianElements, mainGravitationalParameter );

        std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::VectorXd > > benchmarkStateInterpolator;
        {
            // Propagator settings
            std::shared_ptr< TranslationalStatePropagatorSettings< long double > > propagatorSettings =
                    std::make_shared< TranslationalStatePropagatorSettings< long double > >(
                        centralBodies, accelerationModelMap, bodiesToPropagate, satelliteInitialState.cast< long double >( ),
                        simulationEndEpoch, cowell );
            std::shared_ptr< IntegratorSettings< Time > > integratorSettings =
                    std::make_shared< RungeKuttaVariableStepSizeSettingsScalarTolerances< Time > >(
                        simulationStartEpoch, 5.0, RungeKuttaCoefficients::rungeKuttaFehlberg78, 1.0e-9, 5.0,
                        integrationReferenceTolerance, integrationReferenceTolerance );

            std::cout << std::endl << "Benchmark: " << std::endl;

            // Simulate orbit and output computation time
            SingleArcDynamicsSimulator< long double, Time > dynamicsSimulator(
                        bodyMap, integratorSettings, propagatorSettings, true, false, false, false );

            std::map< Time, Eigen::Matrix< long double, Eigen::Dynamic, 1 > > cartesianIntegrationResult =
                    dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
            std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > > cartesianDoubleIntegrationResult;

            utilities::castMatrixMap( cartesianIntegrationResult, cartesianDoubleIntegrationResult );

            benchmarkStateInterpolator = interpolators::createOneDimensionalInterpolator(
                         cartesianDoubleIntegrationResult, std::make_shared< interpolators::LagrangeInterpolatorSettings >( 8 ) );
            // Store number of function evaluations for variable step-size
            std::cout << "Total Number of Function Evaluations: " <<
                      dynamicsSimulator.getCumulativeNumberOfFunctionEvaluations( ).rbegin( )->second  << std::endl;

            // Store propagation time for constant step-size
            std::cout << "Total Propagation Time: " <<
                         dynamicsSimulator.getCumulativeComputationTimeHistory( ).rbegin( )->second << std::endl;

            writeDataMapToTextFile( cartesianDoubleIntegrationResult, "cartesian_benchmark" + fileSuffix + ".dat",
                                    getOutputPath( "EquationsOfMotion/" ) );


        }

        {
            // Propagator settings
            std::shared_ptr< TranslationalStatePropagatorSettings< long double > > propagatorSettings =
                    std::make_shared< TranslationalStatePropagatorSettings< long double > >(
                        centralBodies, accelerationModelMap, bodiesToPropagate, satelliteInitialState.cast< long double >( ),
                        simulationEndEpoch, unified_state_model_exponential_map );
            std::shared_ptr< IntegratorSettings< Time > > integratorSettings =
                    std::make_shared< RungeKuttaVariableStepSizeSettingsScalarTolerances< Time > >(
                        simulationStartEpoch, 5.0, RungeKuttaCoefficients::rungeKuttaFehlberg78, 1.0e-9, 5.0,
                         integrationReferenceTolerance, integrationReferenceTolerance );

            std::cout << std::endl << "Benchmark: " << std::endl;

            // Simulate orbit and output computation time
            SingleArcDynamicsSimulator< long double, Time > dynamicsSimulator(
                        bodyMap, integratorSettings, propagatorSettings, true, false, false, false );

            std::map< Time, Eigen::Matrix< long double, Eigen::Dynamic, 1 > > cartesianIntegrationResult =
                    dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
            std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > > cartesianDoubleIntegrationResult;

            utilities::castMatrixMap( cartesianIntegrationResult, cartesianDoubleIntegrationResult );

            // Store number of function evaluations for variable step-size
            std::cout << "Total Number of Function Evaluations: " <<
                      dynamicsSimulator.getCumulativeNumberOfFunctionEvaluations( ).rbegin( )->second  << std::endl;

            // Store propagation time for constant step-size
            std::cout << "Total Propagation Time: " <<
                         dynamicsSimulator.getCumulativeComputationTimeHistory( ).rbegin( )->second << std::endl;

            std::map< double, Eigen::VectorXd > cartesianBenchmarkInterpolatedResult;

            for( auto stateIterator : cartesianDoubleIntegrationResult )
            {
                cartesianBenchmarkInterpolatedResult[ stateIterator.first ] =
                        benchmarkStateInterpolator->interpolate( stateIterator.first );
            }


            writeDataMapToTextFile( cartesianDoubleIntegrationResult, "cartesian_benchmark_2_" + fileSuffix + ".dat",
                                    getOutputPath( "EquationsOfMotion/" ) );

            writeDataMapToTextFile( cartesianBenchmarkInterpolatedResult, "cartesian_benchmark_interpolated_2_" + fileSuffix + ".dat",
                                    getOutputPath( "EquationsOfMotion/" ) );


        }

        {
            // Propagator settings
            std::shared_ptr< TranslationalStatePropagatorSettings< long double > > propagatorSettings =
                    std::make_shared< TranslationalStatePropagatorSettings< long double > >(
                        centralBodies, accelerationModelMap, bodiesToPropagate, satelliteInitialState.cast< long double >( ),
                        simulationEndEpoch, cowell );
            std::shared_ptr< IntegratorSettings< Time > > integratorSettings =
                    std::make_shared< RungeKuttaVariableStepSizeSettingsScalarTolerances< Time > >(
                        simulationStartEpoch, 5.0, RungeKuttaCoefficients::rungeKuttaFehlberg78, 1.0e-9, 5.0,
                        100.0 * integrationReferenceTolerance, 100.0 * integrationReferenceTolerance );

            std::cout << std::endl << "Benchmark: " << std::endl;

            // Simulate orbit and output computation time
            SingleArcDynamicsSimulator< long double, Time > dynamicsSimulator(
                        bodyMap, integratorSettings, propagatorSettings, true, false, false, false );

            std::map< Time, Eigen::Matrix< long double, Eigen::Dynamic, 1 > > cartesianIntegrationResult =
                    dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
            std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > > cartesianDoubleIntegrationResult;

            utilities::castMatrixMap( cartesianIntegrationResult, cartesianDoubleIntegrationResult );

            // Store number of function evaluations for variable step-size
            std::cout << "Total Number of Function Evaluations: " <<
                      dynamicsSimulator.getCumulativeNumberOfFunctionEvaluations( ).rbegin( )->second  << std::endl;

            // Store propagation time for constant step-size
            std::cout << "Total Propagation Time: " <<
                         dynamicsSimulator.getCumulativeComputationTimeHistory( ).rbegin( )->second << std::endl;

            std::map< double, Eigen::VectorXd > cartesianBenchmarkInterpolatedResult;

            for( auto stateIterator : cartesianDoubleIntegrationResult )
            {
                cartesianBenchmarkInterpolatedResult[ stateIterator.first ] =
                        benchmarkStateInterpolator->interpolate( stateIterator.first );
            }


            writeDataMapToTextFile( cartesianDoubleIntegrationResult, "cartesian_benchmark_3_" + fileSuffix + ".dat",
                                    getOutputPath( "EquationsOfMotion/" ) );

            writeDataMapToTextFile( cartesianBenchmarkInterpolatedResult, "cartesian_benchmark_interpolated_3_" + fileSuffix + ".dat",
                                    getOutputPath( "EquationsOfMotion/" ) );


        }

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////             LOOP OVER PROPAGATORS                  ////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // Vector of function evaluations and times
        std::vector< unsigned int > numberOfFunctionEvaluations;
        std::vector< double > totalPropagationTime;

        // Loop over propagators
        std::vector< string > nameAdditionPropagator = { "_cowell", "_encke", "_kepl", "_equi", "_usm7", "_usm6", "_usmem" };
        std::vector< string > nameAdditionIntegrator = { "_var", "_const", "_var2", "_const2" };
        for ( unsigned int propagatorType = 0; propagatorType < nameAdditionPropagator.size( ); propagatorType++ )
        {
            // Progress
            std::cout << std::endl << "Propagator: " << propagatorType << std::endl;

            // Loop over integrators
            for ( unsigned int integratorType = 0; integratorType < nameAdditionIntegrator.size( ); integratorType++ )
            {
                // Progress
                std::cout << "Integrator: " << integratorType << std::endl;

                ///////////////////////     CREATE SIMULATION SETTINGS          ////////////////////////////////////////////

                // Propagator settings
                std::shared_ptr< TranslationalStatePropagatorSettings< > > propagatorSettings =
                        std::make_shared< TranslationalStatePropagatorSettings< > >(
                            centralBodies, accelerationModelMap, bodiesToPropagate, satelliteInitialState,
                            simulationEndEpoch, static_cast< TranslationalPropagatorType >( propagatorType ) );

                // Integrator settings
                std::shared_ptr< IntegratorSettings< > > integratorSettings;
                // Integrator dependent on loop
                if ( integratorType == 0 || integratorType == 2 )
                {
                    double toleranceMultiplier = ( integratorType >= 2 ) ? 100.0 : 1.0;
                    integratorSettings = std::make_shared< RungeKuttaVariableStepSizeSettingsScalarTolerances< > >(
                                simulationStartEpoch, 100.0, RungeKuttaCoefficients::rungeKuttaFehlberg78, 1.0e-5, 1.0e5,
                                toleranceMultiplier * integrationRelativeTolerance,
                                toleranceMultiplier * integrationAbsoluteTolerance );
                }
                else if ( integratorType == 1 || integratorType == 3 )
                {
                    double stepMultiplier = ( integratorType >= 2 ) ? 5.0 : 1.0;
                    integratorSettings = std::make_shared< IntegratorSettings< > >(
                                rungeKutta4, simulationStartEpoch, stepMultiplier * integrationConstantTimeStepSize );
                }

                ///////////////////////     PROPAGATE ORBIT                     ////////////////////////////////////////////

                SingleArcDynamicsSimulator< > dynamicsSimulator( bodyMap, integratorSettings, propagatorSettings,
                                                                 true, false, false, false );

                // Retrieve results
                std::map< double, Eigen::VectorXd > cartesianIntegrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
                std::map< double, Eigen::VectorXd > cartesianBenchmarkInterpolatedResult;

                for( auto stateIterator : cartesianIntegrationResult )
                {
                    cartesianBenchmarkInterpolatedResult[ stateIterator.first ] =
                            benchmarkStateInterpolator->interpolate( stateIterator.first );
                }

                if ( integratorType == 0 || integratorType == 2 )
                {
                    // Store number of function evaluations for variable step-size
                    numberOfFunctionEvaluations.push_back( dynamicsSimulator.getCumulativeNumberOfFunctionEvaluations( ).rbegin( )->second );
                    std::cout << "Total Number of Function Evaluations: " << numberOfFunctionEvaluations.back( ) << std::endl;
                }
                else
                {
                    // Store propagation time for constant step-size
                    totalPropagationTime.push_back( dynamicsSimulator.getCumulativeComputationTimeHistory( ).rbegin( )->second );
                    std::cout << "Total Propagation Time: " << totalPropagationTime.back( ) << std::endl;
                }


                ///////////////////////     PROVIDE OUTPUT TO FILES             ////////////////////////////////////////////

                // Write perturbed satellite propagation history to file
                writeDataMapToTextFile( cartesianIntegrationResult, "cartesian" + nameAdditionPropagator[ propagatorType ] +
                                        nameAdditionIntegrator[ integratorType ] + fileSuffix + ".dat", getOutputPath( "EquationsOfMotion/" ) );
                writeDataMapToTextFile( cartesianBenchmarkInterpolatedResult, "cartesian_benchmark" + nameAdditionPropagator[ propagatorType ] +
                                        nameAdditionIntegrator[ integratorType ] + fileSuffix + ".dat", getOutputPath( "EquationsOfMotion/" ) );
            }
        }

        // Write function evaluations and times to file
        writeMatrixToFile( utilities::convertStlVectorToEigenVector( numberOfFunctionEvaluations ),
                           "functionEvaluations" + fileSuffix + ".dat", 16, getOutputPath( "EquationsOfMotion/" ) );
        writeMatrixToFile( utilities::convertStlVectorToEigenVector( totalPropagationTime ),
                           "propagationTime"  + fileSuffix + ".dat", 16, getOutputPath( "EquationsOfMotion/" ) );

    }

    // Final statement
    // The exit code EXIT_SUCCESS indicates that the program was successfully executed
    return EXIT_SUCCESS;
}
