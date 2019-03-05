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
    std::string outputDirectory = tudat_applications::getOutputPath( "AccelerationModels/" );

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

    double simulationStartEpoch = 0.0;
    double simulationEndEpoch = 100.0 * physical_constants::JULIAN_YEAR;

    // Create body objects.
    std::vector< std::string > bodiesToCreate;
    //bodiesToCreate.push_back( "Earth" );
    bodiesToCreate.push_back( "Sun" );
    bodiesToCreate.push_back( "Mercury" );
    bodiesToCreate.push_back( "Venus" );
    bodiesToCreate.push_back( "Earth" );
    bodiesToCreate.push_back( "Moon" );
    bodiesToCreate.push_back( "Mars" );
    bodiesToCreate.push_back( "Jupiter" );
//    bodiesToCreate.push_back( "Saturn" );
//    bodiesToCreate.push_back( "Moon" );

    std::cout<<"T1"<<std::endl;
    std::map< std::string, std::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( bodiesToCreate );
    std::cout<<"T1"<<std::endl;

    bodySettings[ "Sun" ]->gravityFieldSettings = std::make_shared< SphericalHarmonicsGravityFieldSettings >(
                getBodyGravitationalParameter( "Sun" ), 695.7E6,
                ( Eigen::Matrix3d( ) << 1.0, 0.0, 0.0,
                  0.0, 0.0, 0.0,
                  2.0E-7 / calculateLegendreGeodesyNormalizationFactor( 2, 0), 0.0, 0.0 ).finished( ),
                Eigen::Matrix3d::Zero( ), "IAU_Sun" );
    std::cout<<"T1"<<std::endl;


    // Create Earth object
    NamedBodyMap bodyMap = createBodies( bodySettings );
    std::cout<<"T1"<<std::endl;

    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );

    for( unsigned int simulationCase = 0; simulationCase < 2; simulationCase++ )
    {
        std::cout<<"Test case: "<<simulationCase<<std::endl;
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////            CREATE ACCELERATIONS          //////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // Define propagator settings variables.
        SelectedAccelerationMap accelerationMap;
        std::vector< std::string > bodiesToPropagate;
        std::vector< std::string > centralBodies;


        for( unsigned int i = 0; i < bodiesToCreate.size( ); i++ )
        {
            bodiesToPropagate.push_back( bodiesToCreate.at( i ) );
            if( bodiesToPropagate.at( i ) != "Moon" )
            {
                centralBodies.push_back( "SSB" );
            }
            else
            {
                centralBodies.push_back( "Earth" );
            }

            for( unsigned int j = 0; j < bodiesToCreate.size( ); j++ )
            {
                std::cout<<i<<" "<<j<<std::endl;
                if( i != j )
                {
                    if( bodiesToCreate.at( j )  != "Sun" )
                    {
                        accelerationMap[ bodiesToCreate.at( i ) ][ bodiesToCreate.at( j ) ].push_back(
                                    std::make_shared< AccelerationSettings >( basic_astrodynamics::central_gravity ) );
                    }
                    else if( simulationCase == 0 )
                    {
                        accelerationMap[ bodiesToCreate.at( i ) ][ bodiesToCreate.at( j ) ].push_back(
                                std::make_shared< AccelerationSettings >( basic_astrodynamics::central_gravity ) );
                    }
                    else if( simulationCase == 1 )
                    {
                        accelerationMap[ bodiesToCreate.at( i ) ][ bodiesToCreate.at( j ) ].push_back(
                                std::make_shared< SphericalHarmonicAccelerationSettings >( 2, 0 ) );
                    }
                }
                std::cout<<i<<" "<<j<<std::endl;
            }
        }

        // Create acceleration models and propagation settings.
        basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                    bodyMap, accelerationMap, bodiesToPropagate, centralBodies );

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////             CREATE PROPAGATION SETTINGS            ////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


        Eigen::VectorXd systemInitialState = getInitialStatesOfBodies(
                   bodiesToPropagate, centralBodies, bodyMap, simulationStartEpoch );

        std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
                std::make_shared< TranslationalStatePropagatorSettings< double > >
                ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState, simulationEndEpoch, cowell );

        std::shared_ptr< IntegratorSettings< > > integratorSettings =
                std::make_shared< RungeKuttaVariableStepSizeSettings< > >(
                    rungeKuttaVariableStepSize, simulationStartEpoch, 3600.0,
                    numerical_integrators::RungeKuttaCoefficients::rungeKuttaFehlberg78,
                    std::numeric_limits< double >::epsilon( ), std::numeric_limits< double >::infinity( ),
                    1.0E-13, 1.0E-13 );

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // Create simulation object and propagate dynamics.
        SingleArcDynamicsSimulator< > dynamicsSimulator(
                    bodyMap, integratorSettings, propagatorSettings );
        std::map< double, Eigen::VectorXd > integrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );

        std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::VectorXd > > stateInterpolator =
                std::make_shared< interpolators::LagrangeInterpolator< double, Eigen::VectorXd > >(
                    integrationResult, 8 );

        std::map< double, Eigen::VectorXd > integrationResultToSave;
        integrationResultToSave[ integrationResult.begin( )->first ] = integrationResult.begin( )->second;
        double currentTime = simulationStartEpoch + physical_constants::JULIAN_DAY * 2.0;
        int numberOfDataPoints = 5.0E4;
        double timeStep = ( simulationEndEpoch - simulationStartEpoch ) / numberOfDataPoints;
        while( currentTime < simulationEndEpoch - physical_constants::JULIAN_DAY * 2.0 )
        {
            integrationResultToSave[ currentTime ] = stateInterpolator->interpolate( currentTime );
            currentTime += timeStep;
        }

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////        PROVIDE OUTPUT TO CONSOLE AND FILES           //////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



        // Write satellite propagation history to file.
        input_output::writeDataMapToTextFile( integrationResultToSave,
                                              "stateSolarSystemSphericalHarmonicCases_full_" +
                                              boost::lexical_cast< std::string >( simulationCase ) +
                                              ".dat",
                                              outputDirectory, "",
                                              std::numeric_limits< double >::digits10,
                                              std::numeric_limits< double >::digits10,
                                              "," );
    }
    // Final statement.
    // The exit code EXIT_SUCCESS indicates that the program was successfully executed.



    return EXIT_SUCCESS;

}
