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


//! Execute propagation of orbit of Asterix around the Earth.
template< typename StateScalarType = double >
void runSimulations( )
{
    std::string outputDirectory = tudat_applications::getOutputPath( "NumericalIntegration/" );

    bool useLongDoubles = false;
    double toleranceFactor = 1.0;
    std::string fileSuffix = "";
    if( !( sizeof( StateScalarType ) == 8 ) )
    {
        useLongDoubles = true;
        fileSuffix = "_long";
        toleranceFactor = 0.01;
    }
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

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////     CREATE ENVIRONMENT AND VEHICLE       //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Load Spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Create body objects.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Earth" );
    bodiesToCreate.push_back( "Moon" );
    bodiesToCreate.push_back( "Sun" );

    std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( bodiesToCreate );


    // Create Earth object
    NamedBodyMap bodyMap = createBodies( bodySettings );

    // Create spacecraft object.
    bodyMap[ "Asterix" ] = boost::make_shared< simulation_setup::Body >( );
    bodyMap[ "Asterix" ]->setConstantBodyMass( 400.0 );

    // Create aerodynamic coefficient interface settings.
    double referenceArea = 4.0;
    double aerodynamicCoefficient = 1.2;
    boost::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings =
            boost::make_shared< ConstantAerodynamicCoefficientSettings >(
                referenceArea, aerodynamicCoefficient * Eigen::Vector3d::UnitX( ), 1, 1 );

    // Create and set aerodynamic coefficients object
    bodyMap[ "Asterix" ]->setAerodynamicCoefficientInterface(
                createAerodynamicCoefficientInterface( aerodynamicCoefficientSettings, "Asterix" ) );

    // Create radiation pressure settings
    double referenceAreaRadiation = 4.0;
    double radiationPressureCoefficient = 1.2;
    std::vector< std::string > occultingBodies;
    occultingBodies.push_back( "Earth" );
    boost::shared_ptr< RadiationPressureInterfaceSettings > asterixRadiationPressureSettings =
            boost::make_shared< CannonBallRadiationPressureInterfaceSettings >(
                "Sun", referenceAreaRadiation, radiationPressureCoefficient, occultingBodies );

    // Create and set radiation pressure settings
    bodyMap[ "Asterix" ]->setRadiationPressureInterface(
                "Sun", createRadiationPressureInterface(
                    asterixRadiationPressureSettings, "Asterix", bodyMap ) );

    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "Earth", "ECLIPJ2000" );

    for( unsigned int accelerationCase = 0; accelerationCase < 1; accelerationCase++ )
    {
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////            CREATE ACCELERATIONS          //////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // Define propagator settings variables.
        SelectedAccelerationMap accelerationMap;
        std::vector< std::string > bodiesToPropagate;
        std::vector< std::string > centralBodies;

        bodiesToPropagate.push_back( "Asterix" );
        centralBodies.push_back( "Earth" );

        // Define propagation settings.
        std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfAsterix;

        if( accelerationCase == 0 )
        {
            accelerationsOfAsterix[ "Earth" ].push_back( boost::make_shared< AccelerationSettings >(
                                                             basic_astrodynamics::central_gravity ) );
        }
        else if( accelerationCase == 1 )
        {
            accelerationsOfAsterix[ "Earth" ].push_back( boost::make_shared< SphericalHarmonicAccelerationSettings >(
                                                             8, 8 ) );
            accelerationsOfAsterix[ "Moon" ].push_back( boost::make_shared< AccelerationSettings >(
                                                            basic_astrodynamics::central_gravity ) );
            accelerationsOfAsterix[ "Sun" ].push_back( boost::make_shared< AccelerationSettings >(
                                                           basic_astrodynamics::central_gravity ) );
            accelerationsOfAsterix[ "Sun" ].push_back( boost::make_shared< AccelerationSettings >(
                                                           basic_astrodynamics::cannon_ball_radiation_pressure ) );

        }
        accelerationMap[  "Asterix" ] = accelerationsOfAsterix;


        // Create acceleration models and propagation settings.
        basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                    bodyMap, accelerationMap, bodiesToPropagate, centralBodies );



        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////             CREATE PROPAGATION SETTINGS            ////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        std::vector< double > eccentricities;
        if( accelerationCase == 0 )
        {
            eccentricities = { 0.01, 0.05, 0.1, 0.25, 0.5, 0.9, 0.95, 0.99 };
        }
        else
        {
            eccentricities = { 0.01, 0.05, 0.1 };
        }

        for( unsigned int i = 0; i < eccentricities.size( ); i++ )
        {
            for( unsigned int j = 0; j < 4; j++ )
            {
                std::map< double, double > functionEvaluationCounter;
                for( unsigned int k = 0; k < 14; k++ )
                {                                       
                    int lmax = 1;
                    if( accelerationCase != 0 )
                    {
                        lmax = 2;
                    }
                    for( unsigned int l = 0; l < lmax; l++ )
                    {
                        std::cout<<i<<" "<<j<<" "<<k<<" "<<l<<std::endl;

                        // Set initial conditions for the Asterix satellite that will be propagated in this simulation.
                        // The initial conditions are given in Keplerian elements and later on converted to Cartesian
                        // elements.

                        // Set Keplerian elements for Asterix.
                        Eigen::Vector6d asterixInitialStateInKeplerianElements;
                        asterixInitialStateInKeplerianElements( semiMajorAxisIndex ) = 7500.0E3;
                        asterixInitialStateInKeplerianElements( eccentricityIndex ) = eccentricities.at( i );
                        asterixInitialStateInKeplerianElements( inclinationIndex ) = convertDegreesToRadians( 85.3 );
                        asterixInitialStateInKeplerianElements( argumentOfPeriapsisIndex )
                                = convertDegreesToRadians( 235.7 );
                        asterixInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex )
                                = convertDegreesToRadians( 23.4 );
                        asterixInitialStateInKeplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 139.87 );

                        // Convert Asterix state from Keplerian elements to Cartesian elements.
                        double earthGravitationalParameter = bodyMap.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );
                        Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > systemInitialState = convertKeplerianToCartesianElements(
                                    asterixInitialStateInKeplerianElements,
                                    earthGravitationalParameter ).template cast< StateScalarType >( );


                        // Set simulation end epoch.
                        const double simulationStartEpoch = 0.0;
                        const double simulationEndEpoch = 14.0 * tudat::physical_constants::JULIAN_DAY;


                        TranslationalPropagatorType propagatorType = cowell;
                        if( l > 0 )
                        {
                            propagatorType = encke;
                        }
                        boost::shared_ptr< TranslationalStatePropagatorSettings< StateScalarType > > propagatorSettings =
                                boost::make_shared< TranslationalStatePropagatorSettings< StateScalarType > >
                                ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState,
                                  boost::make_shared< PropagationTimeTerminationSettings >( simulationEndEpoch, false ), propagatorType );

                        boost::shared_ptr< IntegratorSettings< > > integratorSettings;
                        if( k == 0 )
                        {
                            double timeStep = std::pow( 2.0,  static_cast< double >( 1.0 + j ) );
                            integratorSettings = boost::make_shared< IntegratorSettings< > >
                                    ( rungeKutta4, simulationStartEpoch, timeStep );
                        }
                        else if( k == 1 )
                        {
                            double tolerance = std::pow( 10.0, static_cast< double >( -15.0 + 2.0 * j ) );
                            integratorSettings = boost::make_shared< RungeKuttaVariableStepSizeSettings< > >(
                                        rungeKuttaVariableStepSize, simulationStartEpoch, 10.0,
                                        numerical_integrators::RungeKuttaCoefficients::rungeKuttaFehlberg45,
                                        std::numeric_limits< double >::epsilon( ), std::numeric_limits< double >::infinity( ),
                                        toleranceFactor * tolerance, toleranceFactor * 1.0E-10 );
                        }
                        else if( k == 2 )
                        {
                            double tolerance = std::pow( 10.0, static_cast< double >( -15.0 + 2.0 * j ) );
                            integratorSettings = boost::make_shared< RungeKuttaVariableStepSizeSettings< > >(
                                        rungeKuttaVariableStepSize, simulationStartEpoch, 10.0,
                                        numerical_integrators::RungeKuttaCoefficients::rungeKuttaFehlberg56,
                                        std::numeric_limits< double >::epsilon( ), std::numeric_limits< double >::infinity( ),
                                        toleranceFactor * tolerance, toleranceFactor * 1.0E-10 );
                        }
                        else if( k == 3 )
                        {
                            double tolerance = std::pow( 10.0, static_cast< double >( -15.0 + 2.0 * j ) );
                            integratorSettings = boost::make_shared< RungeKuttaVariableStepSizeSettings< > >(
                                        rungeKuttaVariableStepSize, simulationStartEpoch, 10.0,
                                        numerical_integrators::RungeKuttaCoefficients::rungeKuttaFehlberg78,
                                        std::numeric_limits< double >::epsilon( ), std::numeric_limits< double >::infinity( ),
                                        toleranceFactor * tolerance, toleranceFactor * 1.0E-10 );
                        }
                        else if( k == 4 )
                        {
                            double tolerance = std::pow( 10.0, static_cast< double >( -15.0 + 2.0 * j ) );
                            integratorSettings = boost::make_shared< tudat::numerical_integrators::BulirschStoerIntegratorSettings< > >(
                                        simulationStartEpoch, 10.0,
                                        numerical_integrators::bulirsch_stoer_sequence, 4,
                                        std::numeric_limits< double >::epsilon( ), std::numeric_limits< double >::infinity( ),
                                        toleranceFactor * tolerance, toleranceFactor * 1.0E-10 );
                        }
                        else if( k == 5 )
                        {
                            double tolerance = std::pow( 10.0, static_cast< double >( -15.0 + 2.0 * j ) );
                            integratorSettings = boost::make_shared< tudat::numerical_integrators::BulirschStoerIntegratorSettings< > >(
                                        simulationStartEpoch, 10.0,
                                        numerical_integrators::bulirsch_stoer_sequence, 6,
                                        std::numeric_limits< double >::epsilon( ), std::numeric_limits< double >::infinity( ),
                                        toleranceFactor * tolerance, toleranceFactor * 1.0E-10 );
                        }
                        else if( k == 6 )
                        {
                            double tolerance = std::pow( 10.0, static_cast< double >( -15.0 + 2.0 * j ) );
                            integratorSettings = boost::make_shared< tudat::numerical_integrators::BulirschStoerIntegratorSettings< > >(
                                        simulationStartEpoch, 10.0,
                                        numerical_integrators::bulirsch_stoer_sequence, 8,
                                        std::numeric_limits< double >::epsilon( ), std::numeric_limits< double >::infinity( ),
                                        toleranceFactor * tolerance, toleranceFactor * 1.0E-10 );
                        }
                        else if( k == 7 )
                        {
                            double tolerance = std::pow( 10.0, static_cast< double >( -15.0 + 2.0 * j ) );
                            integratorSettings = boost::make_shared< tudat::numerical_integrators::BulirschStoerIntegratorSettings< > >(
                                        simulationStartEpoch, 10.0,
                                        numerical_integrators::bulirsch_stoer_sequence, 10,
                                        std::numeric_limits< double >::epsilon( ), std::numeric_limits< double >::infinity( ),
                                        toleranceFactor * tolerance, toleranceFactor * 1.0E-10 );
                        }
                        else if( k == 8 )
                        {
                            double tolerance = std::pow( 10.0, static_cast< double >( -15.0 + 2.0 * j ) );
                            integratorSettings = boost::make_shared< RungeKuttaVariableStepSizeSettings< > >(
                                        rungeKuttaVariableStepSize, simulationStartEpoch, 10.0,
                                        numerical_integrators::RungeKuttaCoefficients::rungeKutta87DormandPrince,
                                        std::numeric_limits< double >::epsilon( ), std::numeric_limits< double >::infinity( ),
                                        100.0 * toleranceFactor * tolerance, toleranceFactor * 1.0E-10 );
                        }
                        else if( k == 9 )
                        {
                            double tolerance = std::pow( 10.0, static_cast< double >( -15.0 + 2.0 * j ) );
                            integratorSettings = boost::make_shared< AdamsBashforthMoultonSettings< > >(
                                        simulationStartEpoch, 10.0,
                                        1.0, std::numeric_limits< double >::infinity( ),
                                        toleranceFactor * tolerance, 1E6 * toleranceFactor * tolerance );
                        }
                        else if( k == 10 )
                        {
                            double tolerance = std::pow( 10.0, static_cast< double >( -15.0 + 2.0 * j ) );
                            integratorSettings = boost::make_shared< AdamsBashforthMoultonSettings< > >(
                                        simulationStartEpoch, 10.0,
                                        1.0, std::numeric_limits< double >::infinity( ),
                                        toleranceFactor * tolerance, 1E6 * toleranceFactor * tolerance, 6, 6 );
                        }
                        else if( k == 11 )
                        {
                            double tolerance = std::pow( 10.0, static_cast< double >( -15.0 + 2.0 * j ) );
                            integratorSettings = boost::make_shared< AdamsBashforthMoultonSettings< > >(
                                        simulationStartEpoch, 10.0,
                                        1.0, std::numeric_limits< double >::infinity( ),
                                        toleranceFactor * tolerance, 1E6 * toleranceFactor * tolerance, 8, 8 );
                        }
                        else if( k == 12 )
                        {
                            double tolerance = std::pow( 10.0, static_cast< double >( -15.0 + 2.0 * j ) );
                            integratorSettings = boost::make_shared< AdamsBashforthMoultonSettings< > >(
                                        simulationStartEpoch, 10.0,
                                        1.0, std::numeric_limits< double >::infinity( ),
                                        toleranceFactor * tolerance, 1E6 * toleranceFactor * tolerance, 10, 10 );
                        }
                        else if( k == 13 )
                        {
                            double timeStep = std::pow( 2.0,  static_cast< double >( 1.0 + j ) );
                            integratorSettings = boost::make_shared< AdamsBashforthMoultonSettings< > >(
                                        simulationStartEpoch, timeStep,
                                        timeStep, timeStep, 1.0, 1.0 );
                        }


                        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                        ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
                        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                        // Create simulation object and propagate dynamics.
                        SingleArcDynamicsSimulator< StateScalarType > dynamicsSimulator(
                                    bodyMap, integratorSettings, propagatorSettings );
                        functionEvaluationCounter[ k ] = dynamicsSimulator.getDynamicsStateDerivative( )->getNumberOfFunctionEvaluations( );
                        std::cout<<"Function evaluations: "<<functionEvaluationCounter[ k ]<<std::endl;
                        std::map< double, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > integrationResult =
                                dynamicsSimulator.getEquationsOfMotionNumericalSolution( );

                        std::map< double, Eigen::VectorXd > integrationError;
                        if( accelerationCase == 0 )
                        {

                            // Compare propagated orbit against numerical result.
                            Eigen::VectorXd analyticalSolution;
                            for( typename std::map< double, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >::const_iterator
                                 resultIterator = integrationResult.begin( );
                                 resultIterator != integrationResult.end( ); resultIterator++ )
                            {
                                analyticalSolution = convertKeplerianToCartesianElements( propagateKeplerOrbit(
                                                                                              asterixInitialStateInKeplerianElements, resultIterator->first, earthGravitationalParameter ),
                                                                                          earthGravitationalParameter );
                                integrationError[ resultIterator->first ] = resultIterator->second.template cast< double >( ) -
                                        analyticalSolution;
                            }
                        }

                        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                        ///////////////////////        PROPAGATE BACKWARDS IN TIME                   //////////////////////////////////////////
                        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                        if( k < 4 )
                        {
                            double propagationEndTime = integrationResult.rbegin( )->first;
                            Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > propagationEndState = integrationResult.rbegin( )->second;

                            systemInitialState = propagationEndState;

                            integratorSettings = boost::make_shared< IntegratorSettings< > >
                                    ( rungeKutta4, propagationEndTime, -std::pow( 10.0, j ) );

                            propagatorSettings =
                                    boost::make_shared< TranslationalStatePropagatorSettings< StateScalarType > >
                                    ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState,
                                      boost::make_shared< PropagationTimeTerminationSettings >( simulationStartEpoch, true ), cowell );

                            // Create simulation object and propagate dynamics.
                            SingleArcDynamicsSimulator< StateScalarType > dynamicsSimulator2(
                                        bodyMap, integratorSettings, propagatorSettings );
                            std::map< double, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > integrationResult2 =
                                    dynamicsSimulator2.getEquationsOfMotionNumericalSolution( );
                            std::map< double, Eigen::VectorXd > integrationError2;

                            // Compare propagated orbit against numerical result.
                            Eigen::VectorXd analyticalSolution;
                            for( typename std::map< double, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >::const_iterator
                                 resultIterator = integrationResult2.begin( );
                                 resultIterator != integrationResult2.end( ); resultIterator++ )
                            {
                                analyticalSolution = convertKeplerianToCartesianElements(
                                            propagateKeplerOrbit(
                                                asterixInitialStateInKeplerianElements, resultIterator->first, earthGravitationalParameter ),
                                            earthGravitationalParameter );
                                integrationError2[ resultIterator->first ] = resultIterator->second.template cast< double >( ) -
                                            analyticalSolution;
                            }

                            // Write satellite propagation history to file.
                            input_output::writeDataMapToTextFile( integrationResult2,
                                                                  "numericalKeplerOrbitBack_eccSett_" + boost::lexical_cast< std::string >( i ) +
                                                                  "_intType"  + boost::lexical_cast< std::string >( k ) +
                                                                  "_intSett"  + boost::lexical_cast< std::string >( j ) +
                                                                  fileSuffix +
                                                                  ".dat",
                                                                  outputDirectory,
                                                                  "",
                                                                  std::numeric_limits< double >::digits10,
                                                                  std::numeric_limits< double >::digits10,
                                                                  "," );

                            // Write satellite propagation history to file.
                            input_output::writeDataMapToTextFile( integrationError2,
                                                                  "numericalKeplerOrbitErrorBack_e_" + boost::lexical_cast< std::string >( i ) +
                                                                  "_intType"  + boost::lexical_cast< std::string >( k ) +
                                                                  "_intSett"  + boost::lexical_cast< std::string >( j ) +
                                                                  fileSuffix +
                                                                  ".dat",
                                                                  outputDirectory,
                                                                  "",
                                                                  std::numeric_limits< double >::digits10,
                                                                  std::numeric_limits< double >::digits10,
                                                                  "," );
                        }


                        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                        ///////////////////////        PROVIDE OUTPUT TO CONSOLE AND FILES           //////////////////////////////////////////
                        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                        // Write satellite propagation history to file.
                        input_output::writeDataMapToTextFile( integrationResult,
                                                              "numericalKeplerOrbit_eccSett_" + boost::lexical_cast< std::string >( i ) +
                                                              "_intType"  + boost::lexical_cast< std::string >( k ) +
                                                              "_intSett"  + boost::lexical_cast< std::string >( j ) +
                                                              "_propSett"  + boost::lexical_cast< std::string >( l ) +
                                                              "_accSett"  + boost::lexical_cast< std::string >( accelerationCase ) +
                                                              fileSuffix +
                                                              ".dat",
                                                              outputDirectory,
                                                              "",
                                                              std::numeric_limits< double >::digits10,
                                                              std::numeric_limits< double >::digits10,
                                                              "," );

                        if( accelerationCase == 0 )
                        {
                            // Write satellite propagation history to file.
                            input_output::writeDataMapToTextFile( integrationError,
                                                                  "numericalKeplerOrbitError_e_" + boost::lexical_cast< std::string >( i ) +
                                                                  "_intType"  + boost::lexical_cast< std::string >( k ) +
                                                                  "_intSett"  + boost::lexical_cast< std::string >( j ) +
                                                                  fileSuffix +
                                                                  ".dat",
                                                                  outputDirectory,
                                                                  "",
                                                                  std::numeric_limits< double >::digits10,
                                                                  std::numeric_limits< double >::digits10,
                                                                  "," );
                        }
                    }
                }
                input_output::writeDataMapToTextFile( functionEvaluationCounter ,
                                                      "functionEvaluations_e_" + boost::lexical_cast< std::string >( i ) +
                                                      "_intSett"  + boost::lexical_cast< std::string >( j ) +
                                                      fileSuffix +
                                                      ".dat",
                                                      outputDirectory,
                                                      "",
                                                      std::numeric_limits< double >::digits10,
                                                      std::numeric_limits< double >::digits10,
                                                      "," );
            }
        }
    }

}

int main()
{
    //runSimulations< double >( );
    runSimulations< double >( );
    runSimulations< long double >( );


    // Final statement.
    // The exit code EXIT_SUCCESS indicates that the program was successfully executed.
    return EXIT_SUCCESS;
}

