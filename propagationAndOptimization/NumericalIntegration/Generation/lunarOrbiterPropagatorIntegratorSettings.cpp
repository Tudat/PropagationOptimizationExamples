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


boost::shared_ptr< IntegratorSettings< > > getIntegratorSettings(
        const int j, const int k, const double simulationStartEpoch, const double toleranceFactor, const double timeStepMultiplier = 1.0 )
{
    boost::shared_ptr< IntegratorSettings< > > integratorSettings;
    if( k == 0 )
    {
        double timeStep = timeStepMultiplier * std::pow( 2.0,  static_cast< double >( 1.0 + j ) );
        integratorSettings = boost::make_shared< IntegratorSettings< > >
                ( rungeKutta4, simulationStartEpoch, timeStep );
    }
    else if( k == 1 )
    {
        double tolerance = std::pow( 10.0, static_cast< double >( -15.0 + 2.0 * j ) );
        integratorSettings = boost::make_shared< RungeKuttaVariableStepSizeSettings< > >(
                    rungeKuttaVariableStepSize, simulationStartEpoch, timeStepMultiplier * 10.0,
                    numerical_integrators::RungeKuttaCoefficients::rungeKuttaFehlberg45,
                    std::numeric_limits< double >::epsilon( ), std::numeric_limits< double >::infinity( ),
                    toleranceFactor * tolerance, toleranceFactor * 1.0E-10 );
    }
    else if( k == 2 )
    {
        double tolerance = std::pow( 10.0, static_cast< double >( -15.0 + 2.0 * j ) );
        integratorSettings = boost::make_shared< RungeKuttaVariableStepSizeSettings< > >(
                    rungeKuttaVariableStepSize, simulationStartEpoch, timeStepMultiplier * 10.0,
                    numerical_integrators::RungeKuttaCoefficients::rungeKuttaFehlberg56,
                    std::numeric_limits< double >::epsilon( ), std::numeric_limits< double >::infinity( ),
                    toleranceFactor * tolerance, toleranceFactor * 1.0E-10 );
    }
    else if( k == 3 )
    {
        double tolerance = std::pow( 10.0, static_cast< double >( -15.0 + 2.0 * j ) );
        integratorSettings = boost::make_shared< RungeKuttaVariableStepSizeSettings< > >(
                    rungeKuttaVariableStepSize, simulationStartEpoch, timeStepMultiplier * 10.0,
                    numerical_integrators::RungeKuttaCoefficients::rungeKuttaFehlberg78,
                    std::numeric_limits< double >::epsilon( ), std::numeric_limits< double >::infinity( ),
                    toleranceFactor * tolerance, toleranceFactor * 1.0E-10 );
    }
    else if( k == 4 )
    {
        double tolerance = std::pow( 10.0, static_cast< double >( -15.0 + 2.0 * j ) );
        integratorSettings = boost::make_shared< tudat::numerical_integrators::BulirschStoerIntegratorSettings< > >(
                    simulationStartEpoch, timeStepMultiplier * 10.0,
                    numerical_integrators::bulirsch_stoer_sequence, 4,
                    std::numeric_limits< double >::epsilon( ), std::numeric_limits< double >::infinity( ),
                    toleranceFactor * tolerance, toleranceFactor * 1.0E-10 );
    }
    else if( k == 5 )
    {
        double tolerance = std::pow( 10.0, static_cast< double >( -15.0 + 2.0 * j ) );
        integratorSettings = boost::make_shared< tudat::numerical_integrators::BulirschStoerIntegratorSettings< > >(
                    simulationStartEpoch, timeStepMultiplier * 10.0,
                    numerical_integrators::bulirsch_stoer_sequence, 6,
                    std::numeric_limits< double >::epsilon( ), std::numeric_limits< double >::infinity( ),
                    toleranceFactor * tolerance, toleranceFactor * 1.0E-10 );
    }
    else if( k == 6 )
    {
        double tolerance = std::pow( 10.0, static_cast< double >( -15.0 + 2.0 * j ) );
        integratorSettings = boost::make_shared< tudat::numerical_integrators::BulirschStoerIntegratorSettings< > >(
                    simulationStartEpoch, timeStepMultiplier * 10.0,
                    numerical_integrators::bulirsch_stoer_sequence, 8,
                    std::numeric_limits< double >::epsilon( ), std::numeric_limits< double >::infinity( ),
                    toleranceFactor * tolerance, toleranceFactor * 1.0E-10 );
    }
    else if( k == 7 )
    {
        double tolerance = std::pow( 10.0, static_cast< double >( -15.0 + 2.0 * j ) );
        integratorSettings = boost::make_shared< tudat::numerical_integrators::BulirschStoerIntegratorSettings< > >(
                    simulationStartEpoch, timeStepMultiplier * 10.0,
                    numerical_integrators::bulirsch_stoer_sequence, 10,
                    std::numeric_limits< double >::epsilon( ), std::numeric_limits< double >::infinity( ),
                    toleranceFactor * tolerance, toleranceFactor * 1.0E-10 );
    }
    else if( k == 8 )
    {
        double tolerance = std::pow( 10.0, static_cast< double >( -13.0 + 2.0 * j ) );
        integratorSettings = boost::make_shared< RungeKuttaVariableStepSizeSettings< > >(
                    rungeKuttaVariableStepSize, simulationStartEpoch, timeStepMultiplier * 10.0,
                    numerical_integrators::RungeKuttaCoefficients::rungeKutta87DormandPrince,
                    std::numeric_limits< double >::epsilon( ), std::numeric_limits< double >::infinity( ),
                    100.0 * toleranceFactor * tolerance, toleranceFactor * 1.0E-10 );
    }
    else if( k == 9 )
    {
        double tolerance = std::pow( 10.0, static_cast< double >( -15.0 + 2.0 * j ) );
        integratorSettings = boost::make_shared< AdamsBashforthMoultonSettings< > >(
                    simulationStartEpoch, timeStepMultiplier * 10.0,
                    1.0, std::numeric_limits< double >::infinity( ),
                    toleranceFactor * tolerance, 1E6 * toleranceFactor * tolerance );
    }
    else if( k == 10 )
    {
        double tolerance = std::pow( 10.0, static_cast< double >( -15.0 + 2.0 * j ) );
        integratorSettings = boost::make_shared< AdamsBashforthMoultonSettings< > >(
                    simulationStartEpoch,  timeStepMultiplier * 10.0,
                    1.0, std::numeric_limits< double >::infinity( ),
                    toleranceFactor * tolerance, 1E6 * toleranceFactor * tolerance, 6, 6 );
    }
    else if( k == 11 )
    {
        double tolerance = std::pow( 10.0, static_cast< double >( -15.0 + 2.0 * j ) );
        integratorSettings = boost::make_shared< AdamsBashforthMoultonSettings< > >(
                    simulationStartEpoch,  timeStepMultiplier * 10.0,
                    1.0, std::numeric_limits< double >::infinity( ),
                    toleranceFactor * tolerance, 1E6 * toleranceFactor * tolerance, 8, 8 );
    }
    else if( k == 12 )
    {
        double tolerance = std::pow( 10.0, static_cast< double >( -15.0 + 2.0 * j ) );
        integratorSettings = boost::make_shared< AdamsBashforthMoultonSettings< > >(
                    simulationStartEpoch,  timeStepMultiplier * 10.0,
                    1.0, std::numeric_limits< double >::infinity( ),
                    toleranceFactor * tolerance, 1E6 * toleranceFactor * tolerance, 10, 10 );
    }
    else if( k == 13 )
    {
        double timeStep =  timeStepMultiplier * std::pow( 2.0,  static_cast< double >( 1.0 + j ) );
        integratorSettings = boost::make_shared< AdamsBashforthMoultonSettings< > >(
                    simulationStartEpoch, timeStep,
                    std::fabs( timeStep ), std::fabs( timeStep ), 1.0, 1.0 );
    }
    return integratorSettings;
}

//! Execute propagation of orbit of spacecraft around the Earth.
/*!
 *
 *  This function performs a numerical integration of equations of motion for various settings, with the goal of providing
 *  output for determining the errors incurred when performing a numerical integration. The test case is a spacecraft in LEO
 *  The code runs multiple for-loop, which modify the settings of the simulation. The following iteration variables are used
 *  in the for loops:
 *
 *  - accelerationCase: Values 0 and 1. If equal to 0: inly Earth point-mass is used (has analytical solution: Kepler orbit).
 *    If equal to 1, Earth SH up to D/O 5/5, Sun and Moon third-body perturbations, and Sun radiation pressure are used.
 *
 *  - i: Iterates over the vector 'eccentricities', different values of i will use different values of the initial eccentricity
 *
 *  - j: Defines the time-step/tolerances of the numerical integrator. Low value of j is for strict tolerance/small step size.
 *    For increasing j, larger step sizes/tolerances are used.
 *
 *  - k: Defines the type of numerical integrator used, as defined below:
 *        0: RK4
 *        1: RKF4(5)
 *        2: RKF5(6)
 *        3: RKF7(8)
 *        4: Bulirsch-Stoer (sequence of size 4)
 *        5: Bulirsch-Stoer (sequence of size 6)
 *        6: Bulirsch-Stoer (sequence of size 8)
 *        7: Bulirsch-Stoer (sequence of size 10)
 *        8: DOPRI8(7)
 *        9: ABM, variable order, variable step-size
 *        10: ABM, 6th order, variable step-size
 *        11: ABM, 8th order, variable step-size
 *        12: ABM, 10th order, variable step-size
 *        13: ABM, variable order, fixed step-size
 *
 *  - l: Propagtor that is used:
 *        0: Cowell
 *        1: Gauss - Kepler
 *        2: Gauss - MEE
 *        3: Encke
 */

template< typename StateScalarType = double >
void runSimulations( )
{
    unsigned int numberOfIntegrators = 14;
    unsigned int numberOfPropagators = 4;
    unsigned int numberOfTolerances = 6;
    bool performForwardsBackwardsIntegration = true;

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
    ///////////////////////     CREATE ENVIRONMENT AND VEHICLE       //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Load Spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Create body objects.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Earth" );
    bodiesToCreate.push_back( "Moon" );
    bodiesToCreate.push_back( "Sun" );

    double simulationStartEpoch = 0.0;
    double simulationEndEpoch = 7.0 * tudat::physical_constants::JULIAN_DAY;

    std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( bodiesToCreate, -1.0E8, 1.0E8 );


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
                                                             5, 5 ) );
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
            eccentricities = { 0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95 };
        }
        else
        {
            eccentricities = { 0.01, 0.1, 0.5, 0.9 };
        }

        for( unsigned int i = 0; i < eccentricities.size( ); i++ )
        {
            std::vector< std::vector< std::map< double, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > > > benchmarkResult;
            benchmarkResult.resize( numberOfIntegrators );

            for( unsigned int j = 0; j < numberOfTolerances; j++ )
            {
                std::vector< std::map< double, double > > functionEvaluationCounter;
                std::vector< std::map< double, Eigen::Vector2d > > forwardBackwardError;
                std::vector< std::map< double, Eigen::Vector7d > > propagatedEndStates;

                functionEvaluationCounter.resize( numberOfPropagators);
                forwardBackwardError.resize( numberOfPropagators);
                propagatedEndStates.resize( numberOfPropagators);
                for( unsigned int k = 0; k < numberOfIntegrators; k++ )
                {
                    unsigned int lmax = 1;
                    if( accelerationCase != 0 )
                    {
                        lmax = numberOfPropagators;
                        benchmarkResult[ k ].resize( numberOfPropagators);
                    }
                    else
                    {
                        lmax = numberOfPropagators - 1;
                        benchmarkResult[ k ].resize( numberOfPropagators);
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
                        if( accelerationCase == 1 )
                        {
                            simulationEndEpoch = 7.0 * tudat::physical_constants::JULIAN_DAY;
                        }
                        else
                        {
                            simulationEndEpoch = 7.0 * tudat::physical_constants::JULIAN_DAY;
                        }

                        TranslationalPropagatorType propagatorType = cowell;

                        if( l == 1 )
                        {
                            propagatorType = gauss_keplerian;
                        }
                        else if( l == 2 )
                        {
                            propagatorType = gauss_modified_equinoctial;
                        }
                        else if( l == 3 )
                        {
                            propagatorType = encke;
                        }

                        boost::shared_ptr< TranslationalStatePropagatorSettings< StateScalarType > > propagatorSettings =
                                boost::make_shared< TranslationalStatePropagatorSettings< StateScalarType > >
                                ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState,
                                  boost::make_shared< PropagationTimeTerminationSettings >( simulationEndEpoch, true ), propagatorType );

                        boost::shared_ptr< IntegratorSettings< > > integratorSettings = getIntegratorSettings(
                                    j, k, simulationStartEpoch, toleranceFactor, 1.0 );

                        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                        ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
                        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                        // Create simulation object and propagate dynamics.
                        SingleArcDynamicsSimulator< StateScalarType > dynamicsSimulator(
                                    bodyMap, integratorSettings, propagatorSettings );

                        functionEvaluationCounter[ l ][ k ] =
                                dynamicsSimulator.getDynamicsStateDerivative( )->getNumberOfFunctionEvaluations( );

                        std::map< double, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > integrationResult =
                                dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
                        Eigen::Vector7d vectorToSave;
                        vectorToSave( 0 ) = integrationResult.rbegin( )->first ;
                        vectorToSave.segment( 1, 6 ) = integrationResult.rbegin( )->second.template cast< double >( );
                        propagatedEndStates[ l ][ k ] = vectorToSave;

//                        // Write satellite propagation history to file.
//                        input_output::writeDataMapToTextFile( integrationResult,
//                                                              "perturbedOrbit_e_" + boost::lexical_cast< std::string >( i ) +
//                                                              "_intType"  + boost::lexical_cast< std::string >( k ) +
//                                                              "_intSett"  + boost::lexical_cast< std::string >( j ) +
//                                                              "_propSett"  + boost::lexical_cast< std::string >( l ) +
//                                                              fileSuffix +
//                                                              ".dat",
//                                                              outputDirectory,
//                                                              "",
//                                                              std::numeric_limits< double >::digits10,
//                                                              std::numeric_limits< double >::digits10,
//                                                              "," );

                        std::map< double, double > integrationError;
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
                                integrationError[ resultIterator->first ] = ( resultIterator->second.template cast< double >( ) -
                                                                              analyticalSolution ).segment( 0, 3 ).norm( );
                            }


                            // Write satellite propagation history to file.
                            input_output::writeDataMapToTextFile( integrationError,
                                                                  "numericalKeplerOrbitError_e_" + boost::lexical_cast< std::string >( i ) +
                                                                  "_intType"  + boost::lexical_cast< std::string >( k ) +
                                                                  "_intSett"  + boost::lexical_cast< std::string >( j ) +
                                                                  "_propSett"  + boost::lexical_cast< std::string >( l ) +
                                                                  fileSuffix +
                                                                  ".dat",
                                                                  outputDirectory,
                                                                  "",
                                                                  std::numeric_limits< double >::digits10,
                                                                  std::numeric_limits< double >::digits10,
                                                                  "," );
                        }
                        else if ( accelerationCase == 1 && j == 0 )
                        {
                            benchmarkResult[ k ][ l ] = integrationResult;

                        }
                        else if( accelerationCase == 1 && j > 0 )
                        {

                            boost::shared_ptr< interpolators::OneDimensionalInterpolator<
                                    double, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > > stateInterpolator =
                                    boost::make_shared< interpolators::LagrangeInterpolator<
                                    double, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > >( integrationResult, 8 );

                            boost::shared_ptr< interpolators::OneDimensionalInterpolator<
                                    double, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > > stateInterpolator2 =
                                    boost::make_shared< interpolators::LagrangeInterpolator<
                                    double, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > >( integrationResult, 6 );

                            std::map< double, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > interpolatedResult;
                            std::map< double, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > interpolatedResult2;


                            for( typename std::map< double, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >::iterator it =
                                 benchmarkResult[ k ][ l ].begin( ); it != benchmarkResult[ k ][ l ].end( ); it++ )
                            {
                                interpolatedResult[ it->first ] = stateInterpolator->interpolate(
                                            it->first ) - it->second;
                                interpolatedResult2[ it->first ] = stateInterpolator2->interpolate(
                                            it->first ) - it->second;
                            }

                            // Write satellite propagation history to file.
                            input_output::writeDataMapToTextFile( interpolatedResult,
                                                                  "interpolatedPerturbedOrbit_e_" + boost::lexical_cast< std::string >( i ) +
                                                                  "_intType"  + boost::lexical_cast< std::string >( k ) +
                                                                  "_intSett"  + boost::lexical_cast< std::string >( j ) +
                                                                  "_propSett"  + boost::lexical_cast< std::string >( l ) +
                                                                  fileSuffix +
                                                                  ".dat",
                                                                  outputDirectory,
                                                                  "",
                                                                  std::numeric_limits< double >::digits10,
                                                                  std::numeric_limits< double >::digits10,
                                                                  "," );

                            // Write satellite propagation history to file.
                            input_output::writeDataMapToTextFile( interpolatedResult2,
                                                                  "interpolatedPerturbedOrbitB_e_" + boost::lexical_cast< std::string >( i ) +
                                                                  "_intType"  + boost::lexical_cast< std::string >( k ) +
                                                                  "_intSett"  + boost::lexical_cast< std::string >( j ) +
                                                                  "_propSett"  + boost::lexical_cast< std::string >( l ) +
                                                                  fileSuffix +
                                                                  ".dat",
                                                                  outputDirectory,
                                                                  "",
                                                                  std::numeric_limits< double >::digits10,
                                                                  std::numeric_limits< double >::digits10,
                                                                  "," );
                        }

                        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                        ///////////////////////        PROPAGATE BACKWARDS IN TIME                   //////////////////////////////////////////
                        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                        if( performForwardsBackwardsIntegration )
                        {
                            double propagationEndTime = integrationResult.rbegin( )->first;
                            Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > propagationEndState = integrationResult.rbegin( )->second;
                            Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > newSystemInitialState = propagationEndState;

                            integratorSettings = getIntegratorSettings(
                                        j, k, propagationEndTime, toleranceFactor, -1.0 );

                            propagatorSettings =
                                    boost::make_shared< TranslationalStatePropagatorSettings< StateScalarType > >
                                    ( centralBodies, accelerationModelMap, bodiesToPropagate, newSystemInitialState,
                                      boost::make_shared< PropagationTimeTerminationSettings >( simulationStartEpoch, true ), propagatorType );

                            // Create simulation object and propagate dynamics.
                            SingleArcDynamicsSimulator< StateScalarType > dynamicsSimulator2(
                                        bodyMap, integratorSettings, propagatorSettings );
                            std::map< double, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > integrationResult2 =
                                    dynamicsSimulator2.getEquationsOfMotionNumericalSolution( );
                            forwardBackwardError[ l ][ k ] =
                                    ( Eigen::Vector2d( ) << integrationResult2.begin( )->first,
                                      ( integrationResult2.begin( )->second - integrationResult.begin( )->second ).segment( 0, 3 ).
                                      template cast< double >( ).norm( ) ).finished( );

//                            // Write satellite propagation history to file.
//                            input_output::writeDataMapToTextFile( integrationResult2,
//                                                                  "perturbedOrbitBackward_e_" + boost::lexical_cast< std::string >( i ) +
//                                                                  "_intType"  + boost::lexical_cast< std::string >( k ) +
//                                                                  "_intSett"  + boost::lexical_cast< std::string >( j ) +
//                                                                  "_propSett"  + boost::lexical_cast< std::string >( l ) +
//                                                                  fileSuffix +
//                                                                  ".dat",
//                                                                  outputDirectory,
//                                                                  "",
//                                                                  std::numeric_limits< double >::digits10,
//                                                                  std::numeric_limits< double >::digits10,
//                                                                  "," );


                            std::cout<<"Forward/backward "<<
                                       integrationResult.begin( )->first<<" "<<
                                       integrationResult2.begin( )->first<<" "<<
                                       propagationEndTime<<" "<<
                                       ( integrationResult2.begin( )->second - integrationResult.begin( )->second ).transpose( )<<" "<<
                                       ( integrationResult2.rbegin( )->second - integrationResult.rbegin( )->second ).transpose( )<<std::endl;

                            std::map< double, double > integrationError2;

                            // Compare propagated orbit against numerical result.

                            if( accelerationCase == 0 )
                            {
                                Eigen::VectorXd analyticalSolution;
                                for( typename std::map< double, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >::const_iterator
                                     resultIterator = integrationResult2.begin( );
                                     resultIterator != integrationResult2.end( ); resultIterator++ )
                                {
                                    analyticalSolution = convertKeplerianToCartesianElements(
                                                propagateKeplerOrbit(
                                                    asterixInitialStateInKeplerianElements, resultIterator->first, earthGravitationalParameter ),
                                                earthGravitationalParameter );
                                    integrationError2[ resultIterator->first ] = ( resultIterator->second.template cast< double >( ) -
                                                                                   analyticalSolution ).segment( 0, 3 ).norm( );
                                }

                                // Write forward/backward error to file.
                                input_output::writeDataMapToTextFile( integrationError2,
                                                                      "numericalKeplerOrbitErrorBack_e_" + boost::lexical_cast< std::string >( i ) +
                                                                      "_intType"  + boost::lexical_cast< std::string >( k ) +
                                                                      "_intSett"  + boost::lexical_cast< std::string >( j ) +
                                                                      "_propSett"  + boost::lexical_cast< std::string >( l ) +
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


                for( int l = 0; l < numberOfPropagators; l++ )
                {
                    input_output::writeDataMapToTextFile( functionEvaluationCounter[ l ] ,
                                                          "functionEvaluations_e_" + boost::lexical_cast< std::string >( i ) +
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
                    if( performForwardsBackwardsIntegration )
                    {
                        input_output::writeDataMapToTextFile( forwardBackwardError[ l ] ,
                                                              "forwardBackwardError_e_" + boost::lexical_cast< std::string >( i ) +
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
                        input_output::writeDataMapToTextFile( propagatedEndStates[ l ] ,
                                                              "propagatedEndState_e_" + boost::lexical_cast< std::string >( i ) +
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
                    }
                }
            }
        }
    }

}

int main()
{
    //runSimulations< double >( );
    runSimulations< long double >( );

    // Final statement.
    // The exit code EXIT_SUCCESS indicates that the program was successfully executed.
    return EXIT_SUCCESS;
}

