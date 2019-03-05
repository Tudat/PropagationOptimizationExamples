/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <ctime>

#include <Tudat/SimulationSetup/tudatSimulationHeader.h>

#include "propagationAndOptimization/applicationOutput.h"


//! Execute propagation of orbit of spacecraft around the Earth, using an RK4 integrator with a range of time step
/*!
 *  Execute propagation of orbit of spacecraft around the Earth, using only a point mass Earth acceleration model, using an RK4
 *  integrator with a range of time step. Since the orbit has an analytical solution (Kepler orbit), we can determine the
 *  integration error for each numerical integration. The output of this executable is used to visualize the impact of
 *  truncation and rounding errors.
 *
 *  The simulation is run with double precision state/time values, as well as long double state and split (int + long double)
 *  time representation.
 */
template< typename StateScalarType, typename TimeType >
void runIntegrationErrorSimulation( )
{
    std::string outputDirectory = tudat_applications::getOutputPath( "NumericalIntegration/" );

    std::string fileSuffix = "";
    if( !( sizeof( StateScalarType ) == 8 ) )
    {
        fileSuffix = "_long";
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

    std::map< std::string, std::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( bodiesToCreate );


    // Create Earth object
    NamedBodyMap bodyMap = createBodies( bodySettings );

    // Create spacecraft object.
    bodyMap[ "Asterix" ] = std::make_shared< simulation_setup::Body >( );

    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "Earth", "ECLIPJ2000" );

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
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfAsterix;

    accelerationsOfAsterix[ "Earth" ].push_back( std::make_shared< AccelerationSettings >(
                                                     basic_astrodynamics::central_gravity ) );
    accelerationMap[  "Asterix" ] = accelerationsOfAsterix;


    // Create acceleration models and propagation settings.
    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap, accelerationMap, bodiesToPropagate, centralBodies );



    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE PROPAGATION SETTINGS            ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    std::map< double, Eigen::VectorXd > integrationError;
    std::map< double, double > runTimes;
    for( double stepExponent = -2.0; stepExponent <= 3.0; stepExponent += 0.05 )
    {
        clock_t begin = clock();

        // Set initial conditions for the Asterix satellite that will be propagated in this simulation.
        // The initial conditions are given in Keplerian elements and later on converted to Cartesian
        // elements.

        // Set Keplerian elements for Asterix.
        Eigen::Matrix< StateScalarType, 6, 1 > asterixInitialStateInKeplerianElements;
        asterixInitialStateInKeplerianElements( semiMajorAxisIndex ) = 7500.0E3;
        asterixInitialStateInKeplerianElements( eccentricityIndex ) = 0.1;
        asterixInitialStateInKeplerianElements( inclinationIndex ) = convertDegreesToRadians( 85.3 );
        asterixInitialStateInKeplerianElements( argumentOfPeriapsisIndex )
                = convertDegreesToRadians( 235.7 );
        asterixInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex )
                = convertDegreesToRadians( 23.4 );
        asterixInitialStateInKeplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 139.87 );

        // Convert Asterix state from Keplerian elements to Cartesian elements.
        double earthGravitationalParameter = bodyMap.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );
        Eigen::Matrix< StateScalarType, 6, 1 > systemInitialState = convertKeplerianToCartesianElements< StateScalarType >(
                    asterixInitialStateInKeplerianElements,
                    earthGravitationalParameter );


        // Set simulation end epoch.
        const double simulationStartEpoch = 0.0;
        const double simulationEndEpoch = 3.0 * 3600.0;

        std::cout<<"Exponent of time step: "<<stepExponent<<", time step: "<<std::pow( 10, stepExponent )<<", ";

        TranslationalPropagatorType propagatorType = cowell;
        std::shared_ptr< TranslationalStatePropagatorSettings< StateScalarType > > propagatorSettings =
                std::make_shared< TranslationalStatePropagatorSettings< StateScalarType > >
                ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState.template cast< StateScalarType >( ),
                  simulationEndEpoch, propagatorType );

        std::shared_ptr< IntegratorSettings< TimeType > > integratorSettings =
                std::make_shared< IntegratorSettings< TimeType > >
                    ( rungeKutta4, simulationStartEpoch, std::pow( 10, stepExponent ) );

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


        // Create simulation object and propagate dynamics.
        SingleArcDynamicsSimulator< StateScalarType, TimeType > dynamicsSimulator(
                    bodyMap, integratorSettings, propagatorSettings );
        Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > propagationEndState =
                dynamicsSimulator.getEquationsOfMotionNumericalSolution( ).rbegin( )->second;
        double finalPropagationTime =
                dynamicsSimulator.getEquationsOfMotionNumericalSolution( ).rbegin( )->first;
        clock_t end = clock();

        double elapsedSeconds = double(end - begin) / CLOCKS_PER_SEC;


       Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > analyticalSolution =
               convertKeplerianToCartesianElements< StateScalarType >(
                    propagateKeplerOrbit< StateScalarType >(
                        asterixInitialStateInKeplerianElements, finalPropagationTime, earthGravitationalParameter ),
                    earthGravitationalParameter );

        integrationError[ stepExponent ] = ( propagationEndState - analyticalSolution ).template cast< double >( );
        runTimes[ stepExponent ] = elapsedSeconds;
        std::cout<<"final error: "<<( propagationEndState - analyticalSolution ).transpose( )<<std::endl;
        std::cout<<"Run time: "<<elapsedSeconds<<std::endl<<std::endl;


    }

    input_output::writeDataMapToTextFile( integrationError,
                                          "integrationErrorBehaviour" + fileSuffix + ".dat",
                                          outputDirectory,
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    input_output::writeDataMapToTextFile( runTimes,
                                          "integrationRunTime" + fileSuffix + ".dat",
                                          outputDirectory,
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );
}

int main()
{
    //runIntegrationErrorSimulation< double, double >( );
    runIntegrationErrorSimulation< long double, tudat::Time >( );

    // Final statement.
    // The exit code EXIT_SUCCESS indicates that the program was successfully executed.
    return EXIT_SUCCESS;
}
