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
#include <Tudat/Astrodynamics/Gravitation/triAxialEllipsoidGravity.h>

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
    double simulationEndEpoch = physical_constants::JULIAN_YEAR / 365.0;

    // Load Spice kernels.
    spice_interface::loadStandardSpiceKernels( { input_output::getSpiceKernelPath( ) + "de430_mar097_small.bsp" } );

    // Create body objects.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Earth" );
    bodiesToCreate.push_back( "Sun" );
    bodiesToCreate.push_back( "Mars" );
    bodiesToCreate.push_back( "Venus" );
    bodiesToCreate.push_back( "Phobos" );


    std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( bodiesToCreate );

    Eigen::MatrixXd phobosCosineCoefficients = Eigen::MatrixXd::Zero( 3, 3 );

    phobosCosineCoefficients( 0, 0 ) = 1.0;
    phobosCosineCoefficients( 2, 0 ) = -0.072 / basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 0 );
    phobosCosineCoefficients( 2, 2 ) = -0.048 /basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 2 );


    bodySettings[ "Phobos" ]->gravityFieldSettings = boost::make_shared< SphericalHarmonicsGravityFieldSettings >(
                7.087546066894452E05, gravitation::calculateTriAxialEllipsoidReferenceRadius(
                    13.00E3, 11.39E3, 9.07E3 ), phobosCosineCoefficients,  Eigen::MatrixXd::Zero( 3, 3 ), "IAU_Phobos" );

//    bodySettings[ "Moon" ]->rotationModelSettings = boost::make_shared< RotationModelSettings >(
//                spice_rotation_model, "ECLIPJ2000", "MOON_PA" );
//    boost::dynamic_pointer_cast< SphericalHarmonicsGravityFieldSettings >(
//                bodySettings[ "Moon" ]->gravityFieldSettings )->resetAssociatedReferenceFrame( "MOON_PA" );


    // Create Earth object
    NamedBodyMap bodyMap = createBodies( bodySettings );

    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "Mars", "ECLIPJ2000" );

    for( unsigned int accelerationCase = 0; accelerationCase < 1; accelerationCase++ )
    {
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////            CREATE ACCELERATIONS          //////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // Define propagator settings variables.
        SelectedAccelerationMap accelerationMap;
        std::vector< std::string > bodiesToPropagate;
        std::vector< std::string > centralBodies;

        bodiesToPropagate.push_back( "Phobos" );
        centralBodies.push_back( "Mars" );

        // Define propagation settings.
        std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfPhobos;

        if( accelerationCase == 0 )
        {
            accelerationsOfPhobos[ "Mars" ].push_back( boost::make_shared< MutualExtendedBodySphericalHarmonicAccelerationSettings >(
                                                          2, 2, 2, 2 ) );
        }
        else
        {
            accelerationsOfPhobos[ "Mars" ].push_back( boost::make_shared< MutualSphericalHarmonicAccelerationSettings >(
                                                          2, 2, 2, 2 ) );

        }
        accelerationsOfPhobos[ "Sun" ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );
        accelerationsOfPhobos[ "Earth" ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );
        accelerationsOfPhobos[ "Venus" ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );

        accelerationMap[ "Phobos" ] = accelerationsOfPhobos;


        // Create acceleration models and propagation settings.
        basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                    bodyMap, accelerationMap, bodiesToPropagate, centralBodies );

        Eigen::VectorXd systemInitialState = spice_interface::getBodyCartesianStateAtEpoch(
                    "Phobos", "Mars", "ECLIPJ2000", "None", simulationStartEpoch );

        std::cout<<systemInitialState.transpose( )<<std::endl;

        std::vector< boost::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariables;
        dependentVariables.push_back(
                    boost::make_shared< MutualExtendedSphericalHarmonicAccelerationTermsDependentVariableSaveSettings >(
                        "Phobos", "Mars", 2, 2, 2, 2 ) );


        boost::shared_ptr< PropagatorSettings< double > > propagatorSettings =
                boost::make_shared< TranslationalStatePropagatorSettings< double > >
                ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState, simulationEndEpoch,
                  cowell, boost::make_shared< DependentVariableSaveSettings >( dependentVariables ) );

        boost::shared_ptr< IntegratorSettings< > > integratorSettings =
                        boost::make_shared< RungeKuttaVariableStepSizeSettings< > >
                        ( rungeKuttaVariableStepSize, simulationStartEpoch, 60.0,
                          RungeKuttaCoefficients::rungeKuttaFehlberg78, 60.0, 60.0, 1.0, 1.0 );

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // Create simulation object and propagate dynamics.
        SingleArcDynamicsSimulator< > dynamicsSimulator(
                    bodyMap, integratorSettings, propagatorSettings, true, false, false );
        std::map< double, Eigen::VectorXd > integrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );

        // Write satellite propagation history to file.
        input_output::writeDataMapToTextFile( integrationResult,
                                              "phobosMutualExtendedGravitationalAttractionInfluence_varstep_3_" +
                                              boost::lexical_cast< std::string >( accelerationCase ) +
                                              ".dat",
                                              outputDirectory);

        input_output::writeDataMapToTextFile( dynamicsSimulator.getDependentVariableHistory( ),
                                              "phobosMutualExtendedGravitationalAccelerationTerms_varstep_3_" +
                                              boost::lexical_cast< std::string >( accelerationCase ) +
                                              ".dat",
                                              outputDirectory);

        if( accelerationCase == 0 )
        {


            using namespace observation_models;
            using namespace estimatable_parameters;

            // Check input consistency
            boost::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > translationalPropagatorSettings =
                    boost::dynamic_pointer_cast< propagators::TranslationalStatePropagatorSettings< double > >( propagatorSettings ) ;


            // Getlist of bodies for which the dynamics is to be fit
            std::vector< std::string > observedBodies = translationalPropagatorSettings->bodiesToIntegrate_;

            // Create list of ideal observation settings and initial states to estimate
            std::vector< LinkEnds > linkEndsList;
            ObservationSettingsMap observationSettingsMap;
            std::vector< boost::shared_ptr< EstimatableParameterSettings > > initialStateParameterNames;
            for( unsigned int i = 0; i < observedBodies.size( ); i++ )
            {
                // Add current body to list of observed bodies
                LinkEnds observationLinkEnds;
                observationLinkEnds[ observed_body ] = std::make_pair( observedBodies.at( i ), "" );
                linkEndsList.push_back( observationLinkEnds );
                observationSettingsMap.insert(
                            std::make_pair( observationLinkEnds, boost::make_shared< ObservationSettings >(
                                                position_observable ) ) );

                // Add current body to list of estimated bodies
                initialStateParameterNames.push_back(
                            boost::make_shared< InitialTranslationalStateEstimatableParameterSettings< double > >(
                                observedBodies.at( i ), translationalPropagatorSettings->getInitialStates( ).segment( i * 6, 6 ),
                                translationalPropagatorSettings->centralBodies_.at( i ) ) );
            }

            // Create initial state estimation objects
            boost::shared_ptr< EstimatableParameterSet< double > > initialStateParametersToEstimate =
                    createParametersToEstimate< double >( initialStateParameterNames, bodyMap );

            // Get range over which observations are to be simulated.
            std::pair< double, double > dataTimeInterval =
                    getTabulatedEphemerisSafeInterval( bodyMap.at( observedBodies.at( 0 ) )->getEphemeris( ) );
            double startTime = dataTimeInterval.first;
            double endTime = dataTimeInterval.second;

            // Define list of times at which observations are to be simulated

            double simulatedObservationInterval = 4.0 * 3600.0;

            std::vector< double > baseTimeList;
            double currentTime = startTime;
            while( currentTime < endTime )
            {
                baseTimeList.push_back( currentTime );
                currentTime += simulatedObservationInterval;
            }
            std::map< ObservableType, std::map< LinkEnds, std::pair< std::vector< double >, LinkEndType > > > measurementSimulationInput;
            for( unsigned int i = 0; i < linkEndsList.size( ); i++ )
            {
                measurementSimulationInput[ position_observable ][ linkEndsList.at( i ) ] =
                        std::make_pair( baseTimeList, observed_body );
            }


            // Simulate ideal observations
            typedef Eigen::Matrix< double, Eigen::Dynamic, 1 > ObservationVectorType;
            typedef std::map< LinkEnds, std::pair< ObservationVectorType, std::pair< std::vector< double >, LinkEndType > > > SingleObservablePodInputType;
            typedef std::map< ObservableType, SingleObservablePodInputType > PodInputDataType;
            PodInputDataType observationsAndTimes = simulateObservations< double, double >(
                        measurementSimulationInput, createObservationSimulators(
                            observationSettingsMap, bodyMap ) );


            accelerationsOfMoon.clear( );
            accelerationsOfMoon[ "Earth" ].push_back( boost::make_shared< MutualSphericalHarmonicAccelerationSettings >(
                                                          4, 4, 4, 4 ) );
            accelerationsOfMoon[ "Sun" ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );
            accelerationsOfMoon[ "Jupiter" ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );
            accelerationsOfMoon[ "Mars" ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );
            accelerationsOfMoon[ "Venus" ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );

            accelerationMap.clear( );
            accelerationMap[ "Moon" ] = accelerationsOfMoon;


            // Create acceleration models and propagation settings.
            basic_astrodynamics::AccelerationMap reducedAccelerationModelMap = createAccelerationModelsMap(
                        bodyMap, accelerationMap, bodiesToPropagate, centralBodies );

            Eigen::VectorXd systemInitialState = spice_interface::getBodyCartesianStateAtEpoch(
                        "Moon", "Earth", "ECLIPJ2000", "None", simulationStartEpoch );

            boost::shared_ptr< PropagatorSettings< double > > reducedPropagatorSettings =
                    boost::make_shared< TranslationalStatePropagatorSettings< double > >
                    ( centralBodies, reducedAccelerationModelMap, bodiesToPropagate, systemInitialState, simulationEndEpoch,
                      cowell );

            // Create orbit determination object.
            OrbitDeterminationManager< double, double > orbitDeterminationManager =
                    OrbitDeterminationManager< double, double >(
                        bodyMap, initialStateParametersToEstimate, observationSettingsMap,
                        integratorSettings, reducedPropagatorSettings );

            // Retrieve nominal (e.g. pre-fit) body states
            Eigen::VectorXd nominalBodyStates = initialStateParametersToEstimate->template getFullParameterValues< double >( );

            // Define estimation input
            boost::shared_ptr< PodInput< double, double > > podInput =
                    boost::make_shared< PodInput< double, double > >(
                        observationsAndTimes, initialStateParametersToEstimate->getParameterSetSize( ) );
            podInput->defineEstimationSettings( true, true, false, true, true );


            // Fit nominal dynamics to pertrubed dynamical model
            boost::shared_ptr< PodOutput< double > > podOutput = orbitDeterminationManager.estimateParameters(
                        podInput, boost::make_shared< EstimationConvergenceChecker >( 3 ) );

            std::cout<<initialStateParametersToEstimate->template getFullParameterValues< double >( ) - nominalBodyStates<<std::endl;

            input_output::writeMatrixToFile(
                        podOutput->residualHistory_.at( 0 ), "moonMutualExtendedPreFitResiduals_varstep_3.dat", 16,
                        outputDirectory );
            input_output::writeMatrixToFile(
                        podOutput->residuals_, "moonMutualExtendedPostFitResiduals_varstep_3.dat" , 16,
                        outputDirectory  );
        }
    }



    // Final statement.
    // The exit code EXIT_SUCCESS indicates that the program was successfully executed.
    return EXIT_SUCCESS;
}
