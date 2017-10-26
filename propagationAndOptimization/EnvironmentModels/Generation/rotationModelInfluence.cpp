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

    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "moon_pa_de421_1900-2050.bpc" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "moon_080317.tf" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "moon_assoc_pa.tf" );


    for( unsigned int rotationModelCase = 0; rotationModelCase < 3; rotationModelCase++ )
    {
        // Create body objects.
        std::vector< std::string > bodiesToCreate;
        bodiesToCreate.push_back( "Earth" );
        bodiesToCreate.push_back( "Moon" );
        bodiesToCreate.push_back( "Sun" );

        std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings =
                getDefaultBodySettings( bodiesToCreate, -3600.0, 14.0 * tudat::physical_constants::JULIAN_DAY + 3600.0 );

        std::string moonFrame;
        if( rotationModelCase == 0 )
        {
            moonFrame = "IAU_Moon";
        }
        else if( rotationModelCase == 1 )
        {
            moonFrame = "MOON_PA";
        }
        else if( rotationModelCase == 2 )
        {
            moonFrame = "MOON_ME";
        }
        bodySettings[ "Moon" ]->rotationModelSettings = boost::make_shared< RotationModelSettings >(
                    spice_rotation_model, "ECLIPJ2000", moonFrame );
        boost::dynamic_pointer_cast< SphericalHarmonicsGravityFieldSettings >(
                    bodySettings[ "Moon" ]->gravityFieldSettings )->resetAssociatedReferenceFrame( moonFrame );

        std::cout<<"IAU: "<<std::endl<<spice_interface::computeRotationQuaternionBetweenFrames( "ECLIPJ2000", "IAU_MOON", 0.0 ).toRotationMatrix( )<<std::endl;
        std::cout<<"ME : "<<std::endl<<spice_interface::computeRotationQuaternionBetweenFrames( "ECLIPJ2000", "MOON_ME", 0.0 ).toRotationMatrix( )<<std::endl;
        std::cout<<"PA : "<<std::endl<<spice_interface::computeRotationQuaternionBetweenFrames( "ECLIPJ2000", "MOON_PA", 0.0 ).toRotationMatrix( )<<std::endl;

        // Create Earth object
        NamedBodyMap bodyMap = createBodies( bodySettings );

        // Create spacecraft object.
        bodyMap[ "LunarOrbiter" ] = boost::make_shared< simulation_setup::Body >( );

        // Finalize body creation.
        setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );

        for( unsigned int accelerationCase = 0; accelerationCase < 4; accelerationCase++ )
        {
            std::cout<<"Test case: "<<rotationModelCase<<" "<<accelerationCase<<std::endl;
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

            int maximumDegree, maximumOrder;
            if( accelerationCase == 0 )
            {
                maximumDegree = 16;
                maximumOrder = 16;
            }
            else if( accelerationCase == 1 )
            {
                maximumDegree = 4;
                maximumOrder = 4;
            }
            else if( accelerationCase == 2 )
            {
                maximumDegree = 2;
                maximumOrder = 2;
            }
            else if( accelerationCase == 3 )
            {
                maximumDegree = 2;
                maximumOrder = 0;
            }

            accelerationsOfLunarOrbiter[ "Moon" ].push_back( boost::make_shared< SphericalHarmonicAccelerationSettings >(
                                                                 maximumDegree, maximumOrder ) );
            accelerationsOfLunarOrbiter[ "Earth" ].push_back( boost::make_shared< AccelerationSettings >(
                                                                  basic_astrodynamics::central_gravity ) );
            accelerationsOfLunarOrbiter[ "Sun" ].push_back( boost::make_shared< AccelerationSettings >(
                                                                basic_astrodynamics::central_gravity ) );
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
            const double simulationEndEpoch = 4.0 * tudat::physical_constants::JULIAN_DAY;

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

            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            ///////////////////////        PROVIDE OUTPUT TO CONSOLE AND FILES           //////////////////////////////////////////
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



            // Write satellite propagation history to file.
            input_output::writeDataMapToTextFile( integrationResult,
                                                  "stateMoonOrbiterRotationCases_full_" +
                                                  boost::lexical_cast< std::string >( rotationModelCase ) + "_" +
                                                  boost::lexical_cast< std::string >( accelerationCase ) +
                                                  ".dat",
                                                  outputDirectory );
        }
        // Final statement.
        // The exit code EXIT_SUCCESS indicates that the program was successfully executed.

    }

    return EXIT_SUCCESS;

}
