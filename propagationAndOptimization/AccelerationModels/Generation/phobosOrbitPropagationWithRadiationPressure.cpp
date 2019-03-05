/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <Tudat/SimulationSetup/tudatSimulationHeader.h>
#include <Tudat/External/SpiceInterface/spiceInterface.h>

#include "propagationAndOptimization/applicationOutput.h"

void propagatePhobosOrbit(
        const int testCase )
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
    using namespace tudat::gravitation;
    using namespace tudat::numerical_integrators;
    using namespace tudat::spice_interface;


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////     CREATE ENVIRONMENT AND VEHICLE       //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    // Load Spice kernels.
    spice_interface::loadStandardSpiceKernels( );


    // Set simulation time settings.
    const double simulationStartEpoch = 0.0;
    const double simulationEndEpoch = 10.0 * tudat::physical_constants::JULIAN_YEAR;

    // Define body settings for simulation.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Sun" );
    bodiesToCreate.push_back( "Mars" );

    // Create body objects.
    std::map< std::string, std::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( bodiesToCreate, simulationStartEpoch - 86400.0, simulationEndEpoch + 86400.0 );
    for( unsigned int i = 0; i < bodiesToCreate.size( ); i++ )
    {
        bodySettings[ bodiesToCreate.at( i ) ]->ephemerisSettings->resetFrameOrientation( "J2000" );
        bodySettings[ bodiesToCreate.at( i ) ]->rotationModelSettings->resetOriginalFrame( "J2000" );
    }
    NamedBodyMap bodyMap = createBodies( bodySettings );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE VEHICLE            /////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create spacecraft object.
    bodyMap[ "Phobos" ] = std::make_shared< simulation_setup::Body >( );
    bodyMap[ "Phobos" ]->setConstantBodyMass( 1.0E16 );

    // Create radiation pressure settings
    double referenceAreaRadiation = mathematical_constants::PI * 10.0E3 * 10.0E3;
    double radiationPressureCoefficient = 1.2;
    std::vector< std::string > occultingBodies;
    occultingBodies.push_back( "Mars" );
    std::shared_ptr< RadiationPressureInterfaceSettings > phobosRadiationPressureSettings =
            std::make_shared< CannonBallRadiationPressureInterfaceSettings >(
                "Sun", referenceAreaRadiation, radiationPressureCoefficient, occultingBodies );

    // Create and set radiation pressure settings
    bodyMap[ "Phobos" ]->setRadiationPressureInterface(
                "Sun", createRadiationPressureInterface(
                    phobosRadiationPressureSettings, "Phobos", bodyMap ) );


    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "J2000" );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE ACCELERATIONS          //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define propagator settings variables.
    SelectedAccelerationMap accelerationMap;
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;

    // Define propagation settings.
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfPhobos;
    accelerationsOfPhobos[ "Sun" ].push_back( std::make_shared< AccelerationSettings >(
                                                  basic_astrodynamics::central_gravity ) );
    accelerationsOfPhobos[ "Mars" ].push_back( std::make_shared< AccelerationSettings >(
                                                   basic_astrodynamics::central_gravity ) );

    if( testCase == 1 )
    {
        accelerationsOfPhobos[ "Sun" ].push_back( std::make_shared< AccelerationSettings >(
                                                      basic_astrodynamics::cannon_ball_radiation_pressure ) );
    }

    accelerationMap[ "Phobos" ] = accelerationsOfPhobos;
    bodiesToPropagate.push_back( "Phobos" );
    centralBodies.push_back( "Mars" );

    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap, accelerationMap, bodiesToPropagate, centralBodies );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE PROPAGATION SETTINGS            ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Set Keplerian elements for Phobos.0
    Eigen::Vector6d phobosInitialStateInKeplerianElements;
    phobosInitialStateInKeplerianElements( semiMajorAxisIndex ) = 9378.0E3;
    phobosInitialStateInKeplerianElements( eccentricityIndex ) = 0.0151;
    phobosInitialStateInKeplerianElements( inclinationIndex ) = unit_conversions::convertDegreesToRadians( 1.08 );
    phobosInitialStateInKeplerianElements( argumentOfPeriapsisIndex )
            = unit_conversions::convertDegreesToRadians( 235.7 );
    phobosInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex )
            = unit_conversions::convertDegreesToRadians( 23.4 );
    phobosInitialStateInKeplerianElements( trueAnomalyIndex ) = unit_conversions::convertDegreesToRadians( 139.87 );


    double marsGravitationalParameter = bodyMap.at( "Mars" )->getGravityFieldModel( )->getGravitationalParameter( );
    const Eigen::Vector6d phobosInitialState = convertKeplerianToCartesianElements(
                phobosInitialStateInKeplerianElements, marsGravitationalParameter );


    std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< double > >
            ( centralBodies, accelerationModelMap, bodiesToPropagate, phobosInitialState, simulationEndEpoch, encke );

    const double fixedStepSize = 600.0;
    std::shared_ptr< IntegratorSettings< > > integratorSettings =
            std::make_shared< RungeKuttaVariableStepSizeSettings< > >
            ( rungeKuttaVariableStepSize, 0.0, fixedStepSize,
              RungeKuttaCoefficients::CoefficientSets::rungeKuttaFehlberg78, fixedStepSize, fixedStepSize, 1.0, 1.0 );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    // Create simulation object and propagate dynamics.
    SingleArcDynamicsSimulator< > dynamicsSimulator(
                bodyMap, integratorSettings, propagatorSettings, true, false, false );
    std::map< double, Eigen::VectorXd > integrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
    std::map< double, Eigen::VectorXd > keplerianResults;
    std::map< double, Eigen::VectorXd > meeResults;
    std::map< double, Eigen::VectorXd > keplerOrbitResults;

    for( std::map< double, Eigen::VectorXd >::const_iterator resultIterator = integrationResult.begin( ); resultIterator !=
         integrationResult.end( ); resultIterator++ )
    {
        keplerianResults[ resultIterator->first ] = convertCartesianToKeplerianElements(
                    Eigen::Vector6d( resultIterator->second ), marsGravitationalParameter );
        meeResults[ resultIterator->first ]  = convertCartesianToModifiedEquinoctialElements(
                    Eigen::Vector6d( resultIterator->second ), marsGravitationalParameter, false );
        keplerOrbitResults[ resultIterator->first ]  = convertKeplerianToCartesianElements(
                    propagateKeplerOrbit(
                        phobosInitialStateInKeplerianElements, resultIterator->first, marsGravitationalParameter ),
                    marsGravitationalParameter );
    }

    // Write perturbed satellite propagation history to file.
    input_output::writeDataMapToTextFile( integrationResult,
                                          "phobosPropagationHistorySrp" + boost::lexical_cast< std::string >( testCase ) + ".dat",
                                          outputDirectory,
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    // Write perturbed satellite propagation history to file.
    input_output::writeDataMapToTextFile( keplerianResults,
                                          "phobosPropagationKeplerianHistorySrp" + boost::lexical_cast< std::string >( testCase ) + ".dat",
                                          outputDirectory,
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    // Write perturbed satellite propagation history to file.
    input_output::writeDataMapToTextFile( meeResults,
                                          "phobosPropagationMeeHistorySrp" + boost::lexical_cast< std::string >( testCase ) + ".dat",
                                          outputDirectory,
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    // Write perturbed satellite propagation history to file.
    input_output::writeDataMapToTextFile( keplerOrbitResults,
                                          "unperturbedPhobosKeplerianHistorySrp" + boost::lexical_cast< std::string >( testCase ) + ".dat",
                                          outputDirectory,
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );
}

int main( )
{
    propagatePhobosOrbit( 0 );
    propagatePhobosOrbit( 1 );
}
