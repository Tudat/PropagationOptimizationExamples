/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_MAIN


#include <limits>

#include <Tudat/SimulationSetup/tudatEstimationHeader.h>
#include <Tudat/Astrodynamics/Aerodynamics/UnitTests/testApolloCapsuleCoefficients.h>
#include <Tudat/Astrodynamics/OrbitDetermination/determinePostFitParameterInfluence.h>

#include "propagationAndOptimization/applicationOutput.h"

int main( )
{
    std::string outputPath = tudat_applications::getOutputPath( "UncertaintyModelling/" );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            USING STATEMENTS              //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    using namespace tudat;
    using namespace ephemerides;
    using namespace interpolators;
    using namespace numerical_integrators;
    using namespace spice_interface;
    using namespace simulation_setup;
    using namespace basic_astrodynamics;
    using namespace orbital_element_conversions;
    using namespace propagators;
    using namespace aerodynamics;
    using namespace basic_mathematics;
    using namespace input_output;
    using namespace estimatable_parameters;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE ENVIRONMENT            //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Load Spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Set simulation start epoch.
    const double simulationStartEpoch = 0.0;

    // Set simulation end epoch.
    const double simulationEndEpoch = 3100.0;

    // Set numerical integration fixed step size.
    const double fixedStepSize = 1.0;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE ENVIRONMENT            //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define simulation body settings.
    std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( { "Earth" }, simulationStartEpoch - 10.0 * fixedStepSize,
                                    simulationEndEpoch + 10.0 * fixedStepSize );
    bodySettings[ "Earth" ]->ephemerisSettings = boost::make_shared< simulation_setup::ConstantEphemerisSettings >(
                Eigen::Vector6d::Zero( ), "SSB", "J2000" );
    bodySettings[ "Earth" ]->rotationModelSettings->resetOriginalFrame( "J2000" );

    // Create Earth object
    simulation_setup::NamedBodyMap bodyMap = simulation_setup::createBodies( bodySettings );

    double earthC20 =
            boost::dynamic_pointer_cast< gravitation::SphericalHarmonicsGravityField >(
                bodyMap.at( "Earth" )->getGravityFieldModel( ) )->getCosineCoefficients( )( 2, 0 );
    // Create vehicle objects.
    bodyMap[ "Apollo" ] = boost::make_shared< simulation_setup::Body >( );
    bodyMap[ "Apollo" ]->setEphemeris( boost::make_shared< TabulatedCartesianEphemeris< > >(
                                           boost::shared_ptr< interpolators::OneDimensionalInterpolator
                                           < double, Eigen::Vector6d > >( ), "Earth", "J2000" ) );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE VEHICLE            /////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    // Create vehicle aerodynamic coefficients
    bodyMap[ "Apollo" ]->setAerodynamicCoefficientInterface(
                unit_tests::getApolloCoefficientInterface( ) );
    bodyMap[ "Apollo" ]->setConstantBodyMass( 5.0E3 );

    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "J2000" );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE ACCELERATIONS            ///////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define propagator settings variables.
    SelectedAccelerationMap accelerationMap;
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;

    // Define acceleration model settings.
    std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfApollo;
    accelerationsOfApollo[ "Earth" ].push_back( boost::make_shared< SphericalHarmonicAccelerationSettings >( 2, 0 ) );
    accelerationsOfApollo[ "Earth" ].push_back( boost::make_shared< AccelerationSettings >( aerodynamic ) );
    accelerationMap[  "Apollo" ] = accelerationsOfApollo;

    bodiesToPropagate.push_back( "Apollo" );
    centralBodies.push_back( "Earth" );

    // Create acceleration models
    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap, accelerationMap, bodiesToPropagate, centralBodies );

    // Define constant 30 degree angle of attack
    double constantAngleOfAttack = 0.0 * mathematical_constants::PI / 180.0;
    bodyMap.at( "Apollo" )->getFlightConditions( )->getAerodynamicAngleCalculator( )->setOrientationAngleFunctions(
                boost::lambda::constant( constantAngleOfAttack ) );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE PROPAGATION SETTINGS            ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Set spherical elements for Apollo.
    Eigen::Vector6d apolloSphericalEntryState;
    apolloSphericalEntryState( SphericalOrbitalStateElementIndices::radiusIndex ) =
            spice_interface::getAverageRadius( "Earth" ) + 120.0E3;
    apolloSphericalEntryState( SphericalOrbitalStateElementIndices::latitudeIndex ) =
            25.0 * mathematical_constants::PI / 180.0;
    apolloSphericalEntryState( SphericalOrbitalStateElementIndices::longitudeIndex )  =
            25.0 * mathematical_constants::PI / 180.0;
    apolloSphericalEntryState( SphericalOrbitalStateElementIndices::speedIndex ) = 7.65E3;
    apolloSphericalEntryState( SphericalOrbitalStateElementIndices::flightPathIndex ) =
            -1.25 * mathematical_constants::PI / 180.0;
    apolloSphericalEntryState( SphericalOrbitalStateElementIndices::headingAngleIndex )  =
            25.0 * mathematical_constants::PI / 180.0;

    // Convert apollo state from spherical elements to Cartesian elements.
    Eigen::Vector6d systemInitialState = convertSphericalOrbitalToCartesianState(
                apolloSphericalEntryState );

    boost::shared_ptr< ephemerides::RotationalEphemeris > earthRotationalEphemeris =
            bodyMap.at( "Earth" )->getRotationalEphemeris( );
    systemInitialState = transformStateToGlobalFrame( systemInitialState, simulationStartEpoch, earthRotationalEphemeris );

    // Define list of dependent variables to save.
    std::vector< boost::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariablesList;
    dependentVariablesList.push_back(
                boost::make_shared< SingleDependentVariableSaveSettings >( mach_number_dependent_variable, "Apollo" ) );
    dependentVariablesList.push_back(
                boost::make_shared< SingleDependentVariableSaveSettings >(
                    altitude_dependent_variable, "Apollo", "Earth" ) );
    dependentVariablesList.push_back(
                boost::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                    aerodynamic, "Apollo", "Earth", 1 ) );
    dependentVariablesList.push_back(
                boost::make_shared< SingleDependentVariableSaveSettings >(
                    aerodynamic_force_coefficients_dependent_variable, "Apollo" ) );

    // Create object with list of dependent variables
    boost::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave =
            boost::make_shared< DependentVariableSaveSettings >( dependentVariablesList );

    // Define termination conditions
    boost::shared_ptr< SingleDependentVariableSaveSettings > terminationDependentVariable =
            boost::make_shared< SingleDependentVariableSaveSettings >(
                altitude_dependent_variable, "Apollo", "Earth" );
    boost::shared_ptr< PropagationTerminationSettings > terminationSettings =
            boost::make_shared< PropagationDependentVariableTerminationSettings >(
                terminationDependentVariable, 25.0E3, true );

    // Create propagation settings.
    boost::shared_ptr< PropagatorSettings< double > > propagatorSettings =
            boost::make_shared< TranslationalStatePropagatorSettings< double > >
            ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState,
              terminationSettings, cowell, dependentVariablesToSave );
    boost::shared_ptr< IntegratorSettings< > > integratorSettings =
            boost::make_shared< IntegratorSettings< > >
            ( rungeKutta4, simulationStartEpoch, fixedStepSize );

    boost::shared_ptr< EstimatableParameterSettings > perturbedParameterSettings;
    perturbedParameterSettings = (
                boost::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                    2, 0, 2, 0, "Earth", spherical_harmonics_cosine_coefficient_block ) );

    std::pair< boost::shared_ptr< PodOutput< double > >, Eigen::VectorXd > estimationOutput =
            determinePostfitParameterInfluence(
                bodyMap, integratorSettings, propagatorSettings, perturbedParameterSettings,
                1.0, boost::assign::list_of( -earthC20 ), boost::assign::list_of( 0 ), 4 );

    input_output::writeMatrixToFile(
                estimationOutput.first->residualHistory_.at( 0 ), "entryJ2SensitivityPreFitResiduals.dat"  );
    input_output::writeMatrixToFile(
                estimationOutput.first->residuals_, "entryJ2SensitivityPostFitResiduals.dat"  );


    SingleArcDynamicsSimulator< > simulator( bodyMap, integratorSettings, propagatorSettings );

    input_output::writeDataMapToTextFile( simulator.getEquationsOfMotionNumericalSolution( ),
                                          "entryJ2SensitivityNominalTrajectory.dat", outputPath );

    boost::function< double( ) > positionPerturbationFunction =
            statistics::createBoostContinuousRandomVariableGeneratorFunction(
                statistics::normal_boost_distribution, boost::assign::list_of( 0 )( 100.0 ), 0.0 );

    std::cout<<"Parameter difference "<<estimationOutput.second.transpose( )<<std::endl;

    for( unsigned int i = 0; i < 1000; i++ )
    {
        Eigen::VectorXd statePerturbation = Eigen::VectorXd::Zero( 6 );
        for( unsigned int j = 0; j < 3; j++ )
        {
            statePerturbation( j ) = positionPerturbationFunction( );
        }

        std::cout<<"======================== Perturbation "<<i<<"  "<<statePerturbation.transpose( )<<std::endl;


        Eigen::VectorXd newInitialState = systemInitialState + statePerturbation;

        propagatorSettings->resetInitialStates( newInitialState );
        SingleArcDynamicsSimulator< > perturbedSimulator( bodyMap, integratorSettings, propagatorSettings );
        input_output::writeDataMapToTextFile(
                    perturbedSimulator.getEquationsOfMotionNumericalSolution( ),
                    "entryJ2SensitivityPerturbedTrajectory" + boost::lexical_cast< std::string >( i ) + ".dat", outputPath );
        input_output::writeMatrixToFile( statePerturbation, "entryJ2SensitivityPerturbation" +
                                         boost::lexical_cast< std::string >( i ) + ".dat", 16, outputPath );

    }


}
