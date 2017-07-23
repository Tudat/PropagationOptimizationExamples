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
int main()
{
    std::string outputDirectory = tudat_applications::getOutputPath( "IntergratorAndPropagatorInfluence/" );

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
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "pck00009.tpc" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "de-403-masses.tpc" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "de421.bsp" );


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
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );

    for( unsigned int accelerationCase = 0; accelerationCase < 2; accelerationCase++ )
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
                for( unsigned int k = 0; k < 4; k++ )
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
                        Eigen::VectorXd systemInitialState = convertKeplerianToCartesianElements(
                                    asterixInitialStateInKeplerianElements,
                                    earthGravitationalParameter );


                        // Set simulation end epoch.
                        const double simulationStartEpoch = 0.0;
                        const double simulationEndEpoch = 5.0 * tudat::physical_constants::JULIAN_DAY;


                        TranslationalPropagatorType propagatorType = cowell;
                        if( l > 0 )
                        {
                            propagatorType = encke;
                        }
                        boost::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
                                boost::make_shared< TranslationalStatePropagatorSettings< double > >
                                ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState, simulationEndEpoch, propagatorType );

                        boost::shared_ptr< IntegratorSettings< > > integratorSettings;
                        if( k == 0 )
                        {
                            integratorSettings = boost::make_shared< IntegratorSettings< > >
                                    ( rungeKutta4, simulationStartEpoch, std::pow( 10.0,  static_cast< double >( 1.0 + j ) ) );
                        }
                        else if( k == 1 )
                        {
                            double tolerance = std::pow( 10.0, static_cast< double >( -14.0 + j ) );
                            integratorSettings = boost::make_shared< RungeKuttaVariableStepSizeSettings< > >(
                                        rungeKuttaVariableStepSize, simulationStartEpoch, 10.0,
                                        numerical_integrators::RungeKuttaCoefficients::rungeKuttaFehlberg45,
                                        std::numeric_limits< double >::epsilon( ), std::numeric_limits< double >::infinity( ),
                                        tolerance, tolerance );
                        }
                        else if( k == 2 )
                        {
                            double tolerance = std::pow( 10.0, static_cast< double >( -14.0 + j ) );
                            integratorSettings = boost::make_shared< RungeKuttaVariableStepSizeSettings< > >(
                                        rungeKuttaVariableStepSize, simulationStartEpoch, 10.0,
                                        numerical_integrators::RungeKuttaCoefficients::rungeKuttaFehlberg56,
                                        std::numeric_limits< double >::epsilon( ), std::numeric_limits< double >::infinity( ),
                                        tolerance, tolerance );
                        }
                        else if( k == 3 )
                        {
                            double tolerance = std::pow( 10.0, static_cast< double >( -14.0 + j ) );
                            integratorSettings = boost::make_shared< RungeKuttaVariableStepSizeSettings< > >(
                                        rungeKuttaVariableStepSize, simulationStartEpoch, 10.0,
                                        numerical_integrators::RungeKuttaCoefficients::rungeKuttaFehlberg78,
                                        std::numeric_limits< double >::epsilon( ), std::numeric_limits< double >::infinity( ),
                                        tolerance, tolerance );
                        }

                        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                        ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
                        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                        // Create simulation object and propagate dynamics.
                        SingleArcDynamicsSimulator< > dynamicsSimulator(
                                    bodyMap, integratorSettings, propagatorSettings );
                        std::map< double, Eigen::VectorXd > integrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );

                        std::map< double, Eigen::VectorXd > integrationError;
                        if( accelerationCase == 0 )
                        {

                            // Compare propagated orbit against numerical result.
                            Eigen::VectorXd analyticalSolution;
                            for( std::map< double, Eigen::VectorXd >::const_iterator resultIterator = integrationResult.begin( );
                                 resultIterator != integrationResult.end( ); resultIterator++ )
                            {
                                analyticalSolution = convertKeplerianToCartesianElements( propagateKeplerOrbit(
                                                                                              asterixInitialStateInKeplerianElements, resultIterator->first, earthGravitationalParameter ),
                                                                                          earthGravitationalParameter );
                                integrationError[ resultIterator->first ] = resultIterator->second - analyticalSolution;
                            }
                        }

                        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                        ///////////////////////        PROPAGATE BACKWARDS IN TIME                   //////////////////////////////////////////
                        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                        if( k == 0 && accelerationCase == 0 )
                        {
                            double propagationEndTime = integrationResult.rbegin( )->first;
                            Eigen::VectorXd propagationEndState = integrationResult.rbegin( )->second;

                            systemInitialState = propagationEndState;

                            integratorSettings = boost::make_shared< IntegratorSettings< > >
                                    ( rungeKutta4, propagationEndTime, -std::pow( 10.0, j ) );

                            propagatorSettings =
                                    boost::make_shared< TranslationalStatePropagatorSettings< double > >
                                    ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState, simulationStartEpoch, cowell );

                            // Create simulation object and propagate dynamics.
                            SingleArcDynamicsSimulator< > dynamicsSimulator2(
                                        bodyMap, integratorSettings, propagatorSettings );
                            std::map< double, Eigen::VectorXd > integrationResult2 = dynamicsSimulator2.getEquationsOfMotionNumericalSolution( );
                            std::map< double, Eigen::VectorXd > integrationError2;

                            // Compare propagated orbit against numerical result.
                            Eigen::VectorXd analyticalSolution;
                            for( std::map< double, Eigen::VectorXd >::const_iterator resultIterator = integrationResult2.begin( );
                                 resultIterator != integrationResult2.end( ); resultIterator++ )
                            {
                                analyticalSolution = convertKeplerianToCartesianElements(
                                            propagateKeplerOrbit(
                                                asterixInitialStateInKeplerianElements, resultIterator->first, earthGravitationalParameter ),
                                            earthGravitationalParameter );
                                integrationError2[ resultIterator->first ] = resultIterator->second - analyticalSolution;
                            }

                            // Write satellite propagation history to file.
                            input_output::writeDataMapToTextFile( integrationResult2,
                                                                  "numericalKeplerOrbitBack_eccSett_" + boost::lexical_cast< std::string >( i ) +
                                                                  "_intType"  + boost::lexical_cast< std::string >( k ) +
                                                                  "_intSett"  + boost::lexical_cast< std::string >( j ) +
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


    // Final statement.
    // The exit code EXIT_SUCCESS indicates that the program was successfully executed.
    return EXIT_SUCCESS;
}
