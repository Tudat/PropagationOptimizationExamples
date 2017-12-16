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
#include <Tudat/Astrodynamics/OrbitDetermination/determinePostFitParameterInfluence.h>
#include <Tudat/Astrodynamics/Ephemerides/tabulatedEphemeris.h>

#include "propagationAndOptimization/applicationOutput.h"


//! Execute propagation of orbit of Asterix around the Earth.
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
    using namespace tudat::basic_astrodynamics;
    using namespace tudat::unit_conversions;
    using namespace tudat::estimatable_parameters;
    using namespace tudat::ephemerides;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////     CREATE ENVIRONMENT AND VEHICLE       //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    double simulationStartEpoch = 0.0;
    double simulationEndEpoch = physical_constants::JULIAN_YEAR;

    // Load Spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "moon_pa_de421_1900-2050.bpc" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "moon_080317.tf" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "moon_assoc_pa.tf" );

    // Create body objects.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Earth" );
    bodiesToCreate.push_back( "Moon" );
    bodiesToCreate.push_back( "Sun" );
    bodiesToCreate.push_back( "Jupiter" );
    bodiesToCreate.push_back( "Mars" );
    bodiesToCreate.push_back( "Venus" );

    std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( bodiesToCreate, simulationStartEpoch - 86400.0, simulationEndEpoch + 86400.0  );

    bodySettings[ "Moon" ]->rotationModelSettings = boost::make_shared< RotationModelSettings >(
                spice_rotation_model, "ECLIPJ2000", "MOON_PA" );
    boost::dynamic_pointer_cast< SphericalHarmonicsGravityFieldSettings >(
                bodySettings[ "Moon" ]->gravityFieldSettings )->resetAssociatedReferenceFrame( "MOON_PA" );


    // Create Earth object
    NamedBodyMap bodyMap = createBodies( bodySettings );

    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "Earth", "ECLIPJ2000" );

    for( unsigned int accelerationCase = 0; accelerationCase < 2; accelerationCase++ )
    {
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////            CREATE ACCELERATIONS          //////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // Define propagator settings variables.
        SelectedAccelerationMap accelerationMap;
        std::vector< std::string > bodiesToPropagate;
        std::vector< std::string > centralBodies;

        bodiesToPropagate.push_back( "Moon" );
        centralBodies.push_back( "Earth" );

        // Define propagation settings.
        std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfMoon;

        if( accelerationCase == 0 )
        {
            accelerationsOfMoon[ "Earth" ].push_back( boost::make_shared< MutualSphericalHarmonicAccelerationSettings >(
                                                          4, 4, 0, 0 ) );
        }
        else
        {
            accelerationsOfMoon[ "Earth" ].push_back( boost::make_shared< MutualSphericalHarmonicAccelerationSettings >(
                                                          4, 4, 2, 2 ) );

        }
        accelerationsOfMoon[ "Sun" ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );
        accelerationsOfMoon[ "Jupiter" ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );
        accelerationsOfMoon[ "Mars" ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );
        accelerationsOfMoon[ "Venus" ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );

        accelerationMap[ "Moon" ] = accelerationsOfMoon;


        // Create acceleration models and propagation settings.
        basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                    bodyMap, accelerationMap, bodiesToPropagate, centralBodies );

        Eigen::VectorXd systemInitialState = spice_interface::getBodyCartesianStateAtEpoch(
                    "Moon", "Earth", "ECLIPJ2000", "None", simulationStartEpoch );

        boost::shared_ptr< PropagatorSettings< double > > propagatorSettings =
                boost::make_shared< TranslationalStatePropagatorSettings< double > >
                ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState, simulationEndEpoch,
                  cowell );

        boost::shared_ptr< IntegratorSettings< > > integratorSettings =
                boost::make_shared< IntegratorSettings< > >
                ( rungeKutta4, simulationStartEpoch, 1200.0 );

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // Create simulation object and propagate dynamics.
        SingleArcDynamicsSimulator< > dynamicsSimulator(
                    bodyMap, integratorSettings, propagatorSettings );
        std::map< double, Eigen::VectorXd > integrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );

        // Write satellite propagation history to file.
        input_output::writeDataMapToTextFile( integrationResult,
                                              "moonMutualGravitationalAttractionInfluence_" +
                                              boost::lexical_cast< std::string >( accelerationCase ) +
                                              ".dat",
                                              outputDirectory);

        if( accelerationCase == 1 )
        {
            double moonC20 =
                    boost::dynamic_pointer_cast< gravitation::SphericalHarmonicsGravityField >(
                        bodyMap.at( "Moon" )->getGravityFieldModel( ) )->getCosineCoefficients( )( 2, 0 );
            double moonC21 =
                    boost::dynamic_pointer_cast< gravitation::SphericalHarmonicsGravityField >(
                        bodyMap.at( "Moon" )->getGravityFieldModel( ) )->getCosineCoefficients( )( 2, 1 );
            double moonC22 =
                    boost::dynamic_pointer_cast< gravitation::SphericalHarmonicsGravityField >(
                        bodyMap.at( "Moon" )->getGravityFieldModel( ) )->getCosineCoefficients( )( 2, 2 );

            boost::shared_ptr< EstimatableParameterSettings > perturbedParameterSettings;
            perturbedParameterSettings = (
                        boost::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                            2, 0, 2, 2, "Moon", spherical_harmonics_cosine_coefficient_block ) );

            std::pair< boost::shared_ptr< PodOutput< double > >, Eigen::VectorXd > estimationOutput =
                    determinePostfitParameterInfluence< double, double >(
                        bodyMap, integratorSettings, propagatorSettings, perturbedParameterSettings,
                        6.0 * 3600.0, boost::assign::list_of( -moonC20 )( -moonC21 )( -moonC22 ),
                        boost::assign::list_of( 0 )( 1 )( 2 ), 3 );

            std::cout<<"prefit rms: "<<linear_algebra::getVectorEntryRootMeanSquare(
                           estimationOutput.first->residualHistory_.at( 0 ) )<<std::endl;
            std::cout<<"postfit rms: "<<linear_algebra::getVectorEntryRootMeanSquare(
                           estimationOutput.first->residuals_ )<<std::endl;
            input_output::writeMatrixToFile(
                        estimationOutput.first->residualHistory_.at( 0 ), "moonMutualJ2PreFitResiduals.dat", 16,
                        outputDirectory );
            input_output::writeMatrixToFile(
                        estimationOutput.first->residuals_, "moonMutualJ2PostFitResiduals.dat" , 16,
                        outputDirectory  );
        }
    }



    // Final statement.
    // The exit code EXIT_SUCCESS indicates that the program was successfully executed.
    return EXIT_SUCCESS;
}
