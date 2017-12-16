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


    // Create body objects.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Earth" );
    bodiesToCreate.push_back( "Moon" );
    bodiesToCreate.push_back( "Sun" );

    // Set simulation end epoch.
    const double simulationStartEpoch = 0.0;
    const double simulationEndEpoch = 14.0 * tudat::physical_constants::JULIAN_DAY;

    std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( bodiesToCreate, simulationStartEpoch - 3600.0, simulationEndEpoch + 3600.0 );

    // Create Earth object
    NamedBodyMap bodyMap = createBodies( bodySettings );

    // Create spacecraft object.
    bodyMap[ "LunarOrbiter" ] = boost::make_shared< simulation_setup::Body >( );

    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "Moon", "ECLIPJ2000" );

    for( unsigned int eccentricityCase = 0; eccentricityCase < 2; eccentricityCase++ )
    {
        for( unsigned int propagationCase = 0; propagationCase < 10; propagationCase++ )
        {
            std::cout<<"Test case: "<<propagationCase<<" "<<std::endl;


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
            std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfLunarOrbiter;

            accelerationsOfLunarOrbiter[ "Moon" ].push_back( boost::make_shared< AccelerationSettings >(
                                                                 basic_astrodynamics::central_gravity ) );

            double empiricalAccelerationNorm = 1.0E-8;
            if( propagationCase == 1 )
            {
                accelerationsOfLunarOrbiter[ "Moon" ].push_back( boost::make_shared< EmpiricalAccelerationSettings >(
                                                                     empiricalAccelerationNorm * Eigen::Vector3d::UnitX( ) ) );
            }
            else if( propagationCase == 2 )
            {
                accelerationsOfLunarOrbiter[ "Moon" ].push_back( boost::make_shared< EmpiricalAccelerationSettings >(
                                                                     empiricalAccelerationNorm * Eigen::Vector3d::UnitY( ) ) );
            }
            else if( propagationCase == 3 )
            {
                accelerationsOfLunarOrbiter[ "Moon" ].push_back( boost::make_shared< EmpiricalAccelerationSettings >(
                                                                     empiricalAccelerationNorm * Eigen::Vector3d::UnitZ( ) ) );
            }
            else if( propagationCase == 4 )
            {
                accelerationsOfLunarOrbiter[ "Moon" ].push_back( boost::make_shared< EmpiricalAccelerationSettings >(
                                                                     Eigen::Vector3d::Zero( ),
                                                                     empiricalAccelerationNorm * Eigen::Vector3d::UnitX( ) ) );
            }
            else if( propagationCase == 5 )
            {
                accelerationsOfLunarOrbiter[ "Moon" ].push_back( boost::make_shared< EmpiricalAccelerationSettings >(
                                                                     Eigen::Vector3d::Zero( ),
                                                                     empiricalAccelerationNorm * Eigen::Vector3d::UnitY( ) ) );
            }
            else if( propagationCase == 6 )
            {
                accelerationsOfLunarOrbiter[ "Moon" ].push_back( boost::make_shared< EmpiricalAccelerationSettings >(
                                                                     Eigen::Vector3d::Zero( ),
                                                                     empiricalAccelerationNorm * Eigen::Vector3d::UnitZ( ) ) );
            }
            else if( propagationCase == 7 )
            {
                accelerationsOfLunarOrbiter[ "Moon" ].push_back( boost::make_shared< EmpiricalAccelerationSettings >(
                                                                     Eigen::Vector3d::Zero( ),
                                                                     Eigen::Vector3d::Zero( ),
                                                                     empiricalAccelerationNorm * Eigen::Vector3d::UnitX( ) ) );
            }
            else if( propagationCase == 8 )
            {
                accelerationsOfLunarOrbiter[ "Moon" ].push_back( boost::make_shared< EmpiricalAccelerationSettings >(
                                                                     Eigen::Vector3d::Zero( ),
                                                                     Eigen::Vector3d::Zero( ),
                                                                     empiricalAccelerationNorm * Eigen::Vector3d::UnitY( ) ) );
            }
            else if( propagationCase == 9 )
            {
                accelerationsOfLunarOrbiter[ "Moon" ].push_back( boost::make_shared< EmpiricalAccelerationSettings >(
                                                                     Eigen::Vector3d::Zero( ),
                                                                     Eigen::Vector3d::Zero( ),
                                                                     empiricalAccelerationNorm * Eigen::Vector3d::UnitZ( ) ) );
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
            if( eccentricityCase == 0 )
            {
                lunarOrbiterInitialStateInKeplerianElements( eccentricityIndex ) = 0.01;
            }
            else if( eccentricityCase == 1 )
            {
                lunarOrbiterInitialStateInKeplerianElements( eccentricityIndex ) = 0.05;
            }
            lunarOrbiterInitialStateInKeplerianElements( inclinationIndex ) = convertDegreesToRadians( 89.3 );
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

            // Define list of dependent variables to save.
            std::vector< boost::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariablesList;
            dependentVariablesList.push_back(
                        boost::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                            basic_astrodynamics::empirical_acceleration, "LunarOrbiter", "Moon", 0 ) );

            // Create object with list of dependent variables
            boost::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave =
                    boost::make_shared< DependentVariableSaveSettings >( dependentVariablesList );

            boost::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
                    boost::make_shared< TranslationalStatePropagatorSettings< double > >
                    ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState, simulationEndEpoch, encke );

            boost::shared_ptr< IntegratorSettings< > > integratorSettings = boost::make_shared< IntegratorSettings< > >
                    ( rungeKutta4, simulationStartEpoch, 10.0 );

            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            // Create simulation object and propagate dynamics.
            SingleArcDynamicsSimulator< > dynamicsSimulator(
                        bodyMap, integratorSettings, propagatorSettings );
            std::map< double, Eigen::VectorXd > integrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );

            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            ///////////////////////        PROVIDE OUTPUT TO CONSOLE AND FILES           //////////////////////////////////////////
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



            // Write satellite propagation history to file.
            input_output::writeDataMapToTextFile( integrationResult,
                                                  "accelerationDirectionInfluence_" +
                                                  boost::lexical_cast< std::string >( propagationCase ) + "_" +
                                                  boost::lexical_cast< std::string >( eccentricityCase ) + ".dat",
                                                  outputDirectory );
            // Final statement.
            // The exit code EXIT_SUCCESS indicates that the program was successfully executed.

        }
    }

    return EXIT_SUCCESS;

}
