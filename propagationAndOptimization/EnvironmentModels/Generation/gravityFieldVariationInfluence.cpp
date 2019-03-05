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
using namespace tudat::ephemerides;
using namespace tudat::spice_interface;

//! Function to get tidal deformation model for earth
std::vector< std::shared_ptr< GravityFieldVariationSettings > > getEarthGravityFieldVariationSettings( )
{
    std::vector< std::shared_ptr< GravityFieldVariationSettings > > gravityFieldVariations;

    std::vector< std::string > deformingBodies;
    deformingBodies.push_back( "Moon" );

    std::vector< std::vector< std::complex< double > > > loveNumbers;

    std::vector< std::complex< double > > degreeTwoLoveNumbers_;
    degreeTwoLoveNumbers_.push_back( std::complex< double >( 0.29525, -0.00087 ) );
    degreeTwoLoveNumbers_.push_back( std::complex< double >( 0.29525, -0.00087 ) );
    degreeTwoLoveNumbers_.push_back( std::complex< double >( 0.29525, -0.00087 ) );
    std::vector< std::complex< double > > degreeThreeLoveNumbers_;
    degreeThreeLoveNumbers_.push_back( std::complex< double >( 0.093, 0.0 ) );
    degreeThreeLoveNumbers_.push_back( std::complex< double >( 0.093, 0.0 ) );
    degreeThreeLoveNumbers_.push_back( std::complex< double >( 0.093, 0.0 ) );
    degreeThreeLoveNumbers_.push_back( std::complex< double >( 0.093, 0.0 ) );
    loveNumbers.push_back( degreeTwoLoveNumbers_ );
    loveNumbers.push_back( degreeThreeLoveNumbers_ );


    std::shared_ptr< GravityFieldVariationSettings > moonGravityFieldVariation =
            std::make_shared< BasicSolidBodyGravityFieldVariationSettings >(
                deformingBodies, loveNumbers, 6378137.0 );
    gravityFieldVariations.push_back( moonGravityFieldVariation );

    deformingBodies[ 0 ] = "Sun";
    std::shared_ptr< GravityFieldVariationSettings > sunSingleGravityFieldVariation =
            std::make_shared< BasicSolidBodyGravityFieldVariationSettings >(
                deformingBodies, loveNumbers, 6378137.0 );
    gravityFieldVariations.push_back( sunSingleGravityFieldVariation );

    return gravityFieldVariations;
}

//! Execute propagation of orbit of earthOrbiter around the earth.
int main()
{
    std::string outputDirectory = tudat_applications::getOutputPath( "EnvironmentModels/" );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////     CREATE ENVIRONMENT AND VEHICLE       //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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

    std::map< std::string, std::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( bodiesToCreate, -3600.0, 14.0 * tudat::physical_constants::JULIAN_DAY + 3600.0 );
    bodySettings[ "Earth" ]->gravityFieldVariationSettings = getEarthGravityFieldVariationSettings( );

    // Create earth object
    NamedBodyMap bodyMap = createBodies( bodySettings );

    // Create spacecraft object.
    // Create spacecraft object.
    bodyMap[ "EarthOrbiter" ] = std::make_shared< simulation_setup::Body >( );
    bodyMap[ "EarthOrbiter" ]->setConstantBodyMass( 400.0 );

    // Create aerodynamic coefficient interface settings.
    double referenceArea = 4.0;
    double aerodynamicCoefficient = 1.2;
    std::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings =
            std::make_shared< ConstantAerodynamicCoefficientSettings >(
                referenceArea, aerodynamicCoefficient * Eigen::Vector3d::UnitX( ), 1, 1 );

    // Create and set aerodynamic coefficients object
    bodyMap[ "EarthOrbiter" ]->setAerodynamicCoefficientInterface(
                createAerodynamicCoefficientInterface( aerodynamicCoefficientSettings, "EarthOrbiter" ) );

    // Create radiation pressure settings
    double referenceAreaRadiation = 4.0;
    double radiationPressureCoefficient = 1.2;
    std::vector< std::string > occultingBodies;
    occultingBodies.push_back( "Earth" );
    std::shared_ptr< RadiationPressureInterfaceSettings > asterixRadiationPressureSettings =
            std::make_shared< CannonBallRadiationPressureInterfaceSettings >(
                "Sun", referenceAreaRadiation, radiationPressureCoefficient, occultingBodies );

    // Create and set radiation pressure settings
    bodyMap[ "EarthOrbiter" ]->setRadiationPressureInterface(
                "Sun", createRadiationPressureInterface(
                    asterixRadiationPressureSettings, "EarthOrbiter", bodyMap ) );

    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE ACCELERATIONS          //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define propagator settings variables.
    SelectedAccelerationMap accelerationMap;
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;

    bodiesToPropagate.push_back( "EarthOrbiter" );

    centralBodies.push_back( "Earth" );


    // Define propagation settings.
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfEarthOrbiter;

    accelerationsOfEarthOrbiter[ "Earth" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >(
                                                          16, 16 ) );
    accelerationsOfEarthOrbiter[ "Moon" ].push_back( std::make_shared< AccelerationSettings >(
                                                         basic_astrodynamics::central_gravity ) );
    accelerationsOfEarthOrbiter[ "Sun" ].push_back( std::make_shared< AccelerationSettings >(
                                                        basic_astrodynamics::central_gravity ) );
    accelerationMap[  "EarthOrbiter" ] = accelerationsOfEarthOrbiter;


    // Create acceleration models and propagation settings.
    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap, accelerationMap, bodiesToPropagate, centralBodies );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE PROPAGATION SETTINGS            ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Set Keplerian elements for earthOrbiter.
    Eigen::Vector6d earthOrbiterInitialStateInKeplerianElements;
    earthOrbiterInitialStateInKeplerianElements( semiMajorAxisIndex ) = 6850.0E3;
    earthOrbiterInitialStateInKeplerianElements( eccentricityIndex ) = 0.001;
    earthOrbiterInitialStateInKeplerianElements( inclinationIndex ) = convertDegreesToRadians( 85.3 );
    earthOrbiterInitialStateInKeplerianElements( argumentOfPeriapsisIndex )
            = convertDegreesToRadians( 235.7 );
    earthOrbiterInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex )
            = convertDegreesToRadians( 23.4 );
    earthOrbiterInitialStateInKeplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 139.87 );

    // Convert earthOrbiter state from Keplerian elements to Cartesian elements.
    double earthGravitationalParameter = bodyMap.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );
    Eigen::VectorXd systemInitialState = convertKeplerianToCartesianElements(
                earthOrbiterInitialStateInKeplerianElements,
                earthGravitationalParameter ) + spice_interface::getBodyCartesianStateAtEpoch(
                "Earth", centralBodies.at( 0 ), "ECLIPJ2000", "None", 0.0 );

    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariablesList;
    dependentVariablesList.push_back( std::make_shared< SingleDependentVariableSaveSettings >(
                                          total_gravity_field_variation_acceleration,  "EarthOrbiter", "Earth" ) );
    dependentVariablesList.push_back(
                std::make_shared< SingleVariationSphericalHarmonicAccelerationSaveSettings >(
                    "EarthOrbiter", "Earth", gravitation::basic_solid_body, "Sun" ) );
    dependentVariablesList.push_back(
                std::make_shared< SingleVariationSphericalHarmonicAccelerationSaveSettings >(
                    "EarthOrbiter", "Earth", gravitation::basic_solid_body, "Moon" ) );
    dependentVariablesList.push_back( std::make_shared< SingleVariationSingleTermSphericalHarmonicAccelerationSaveSettings >(
                                          "EarthOrbiter", "Earth", 3, 3, gravitation::basic_solid_body, "Sun" ) );
    dependentVariablesList.push_back( std::make_shared< SingleVariationSingleTermSphericalHarmonicAccelerationSaveSettings >(
                                          "EarthOrbiter", "Earth", 3, 3, gravitation::basic_solid_body, "Moon" ) );

    // Create object with list of dependent variables
    std::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave =
            std::make_shared< DependentVariableSaveSettings >( dependentVariablesList, false );

    // Set simulation end epoch.
    const double simulationStartEpoch = 0.0;
    const double simulationEndEpoch = 4.0 * tudat::physical_constants::JULIAN_DAY;

    std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< double > >
            ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState, simulationEndEpoch, cowell,
              dependentVariablesToSave );

    std::shared_ptr< IntegratorSettings< > > integratorSettings = std::make_shared< IntegratorSettings< > >
            ( rungeKutta4, simulationStartEpoch, 10.0 );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create simulation object and propagate dynamics.
    SingleArcDynamicsSimulator< > dynamicsSimulator(
                bodyMap, integratorSettings, propagatorSettings );
    std::map< double, Eigen::VectorXd > dependentVariableResult = dynamicsSimulator.getDependentVariableHistory( );
    std::map< double, Eigen::VectorXd > sphericalHarmonicVariations;
    std::shared_ptr< gravitation::TimeDependentSphericalHarmonicsGravityField > earthGravityField =
            std::dynamic_pointer_cast< gravitation::TimeDependentSphericalHarmonicsGravityField >(
                bodyMap.at( "Earth" )->getGravityFieldModel( ) );
    for( auto it : dependentVariableResult )
    {
        dynamicsSimulator.getDynamicsStateDerivative( )->computeStateDerivative(
                    it.first, it.second );
        Eigen::MatrixXd cosineVariations = earthGravityField->getTotalCosineCoefficientCorrection( 3, 3 );
        Eigen::MatrixXd sineVariations = earthGravityField->getTotalSineCoefficientCorrection( 3, 3 );
        Eigen::VectorXd currentOutput = Eigen::VectorXd( 12 );
        currentOutput.segment( 0, 3 ) = cosineVariations.block( 2, 0, 1, 3 ).transpose( );
        currentOutput.segment( 3, 4 ) = cosineVariations.block( 3, 0, 1, 4 ).transpose( );
        currentOutput.segment( 7, 2 ) = sineVariations.block( 2, 1, 1, 2 ).transpose( );
        currentOutput.segment( 9, 3 ) = sineVariations.block( 3, 1, 1, 3 ).transpose( );
        sphericalHarmonicVariations[ it.first ] = currentOutput;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////        PROVIDE OUTPUT TO CONSOLE AND FILES           //////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



    // Write satellite propagation history to file.
    input_output::writeDataMapToTextFile( dependentVariableResult,
                                          "stateEarthOrbiterGravityFieldVariations.dat",
                                          outputDirectory );

    input_output::writeDataMapToTextFile( sphericalHarmonicVariations,
                                          "stateEarthOrbiterGravityFieldVariationsCoefficients.dat",
                                          outputDirectory );


    // The exit code EXIT_SUCCESS indicates that the program was successfully executed.

    return EXIT_SUCCESS;
}
