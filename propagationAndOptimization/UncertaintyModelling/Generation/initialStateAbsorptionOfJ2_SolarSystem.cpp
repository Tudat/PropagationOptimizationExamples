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
#include <Tudat/Astrodynamics/OrbitDetermination/determinePostFitParameterInfluence.h>

#include "propagationAndOptimization/applicationOutput.h"

int main( )
{
    std::string outputPath = tudat_applications::getOutputPath( "UncertaintyModelling/" );

    using namespace tudat::simulation_setup;
    using namespace tudat::estimatable_parameters;
    using namespace tudat::propagators;
    using namespace tudat::numerical_integrators;
    using namespace tudat::orbital_element_conversions;
    using namespace tudat::basic_mathematics;
    using namespace tudat::unit_conversions;
    using namespace tudat::ephemerides;
    using namespace tudat::spice_interface;
    using namespace tudat;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////     CREATE ENVIRONMENT AND VEHICLE       //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Load Spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    double simulationStartEpoch = 0.0;
    double simulationEndEpoch = 25.0 * physical_constants::JULIAN_YEAR;

    // Create body objects.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Sun" );
    bodiesToCreate.push_back( "Mercury" );
    bodiesToCreate.push_back( "Venus" );
    bodiesToCreate.push_back( "Earth" );
    bodiesToCreate.push_back( "Mars" );
    bodiesToCreate.push_back( "Jupiter" );
    bodiesToCreate.push_back( "Moon" );


    for( unsigned int i = 0; i < bodiesToCreate.size( ) - 1; i++ )
    {
        std::string targetBody = bodiesToCreate.at( i );

        std::cout<<"************************************* "<<targetBody<<" ************************************"<<std::endl;
        std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings =
                getDefaultBodySettings( bodiesToCreate );

        double sunNormalizedJ2 = 2.0E-7 / calculateLegendreGeodesyNormalizationFactor( 2, 0 );
        bodySettings[ "Sun" ]->gravityFieldSettings = boost::make_shared< SphericalHarmonicsGravityFieldSettings >(
                    getBodyGravitationalParameter( "Sun" ), 695.7E6,
                    ( Eigen::Matrix3d( ) << 1.0, 0.0, 0.0,
                      0.0, 0.0, 0.0,
                      sunNormalizedJ2 , 0.0, 0.0 ).finished( ),
                    Eigen::Matrix3d::Zero( ), "IAU_Sun" );
        bodySettings[ targetBody ]->ephemerisSettings = boost::make_shared< InterpolatedSpiceEphemerisSettings >(
                    simulationStartEpoch, simulationEndEpoch, 300.0, "SSB", "ECLIPJ2000" );

        NamedBodyMap bodyMap = createBodies( bodySettings );

        // Finalize body creation.
        setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );

        { // Define propagator settings variables.

            SelectedAccelerationMap accelerationMap;
            std::vector< std::string > bodiesToPropagate;
            bodiesToPropagate.push_back( targetBody );
            std::vector< std::string > centralBodies;
            if( targetBody != "Moon" )
            {
                centralBodies.push_back( "SSB" );
            }
            else
            {
                centralBodies.push_back( "Earth" );
            }

            for( unsigned int j = 0; j < bodiesToCreate.size( ); j++ )
            {
                if( bodiesToCreate.at( j ) != targetBody )
                {
                    if( bodiesToCreate.at( j )  != "Sun" )
                    {
                        accelerationMap[ targetBody ][ bodiesToCreate.at( j ) ].push_back(
                                    boost::make_shared< AccelerationSettings >( basic_astrodynamics::central_gravity ) );
                    }
                    else
                    {
                        accelerationMap[ targetBody ][ bodiesToCreate.at( j ) ].push_back(
                                    boost::make_shared< SphericalHarmonicAccelerationSettings >( 2, 0 ) );
                    }
                }
            }
            // Create acceleration models and propagation settings.
            basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                        bodyMap, accelerationMap, bodiesToPropagate, centralBodies );


            Eigen::VectorXd systemInitialState = getInitialStatesOfBodies(
                        bodiesToPropagate, centralBodies, bodyMap, simulationStartEpoch );

            boost::shared_ptr< PropagatorSettings< double > > propagatorSettings =
                    boost::make_shared< TranslationalStatePropagatorSettings< double > >
                    ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState, simulationEndEpoch, cowell );

            boost::shared_ptr< IntegratorSettings< > > integratorSettings =
                    boost::make_shared< RungeKuttaVariableStepSizeSettings< > >
                    ( rungeKuttaVariableStepSize, ( simulationStartEpoch ), 12.0 * 3600.0,
                      RungeKuttaCoefficients::CoefficientSets::rungeKuttaFehlberg78,
                      12.0 * 3600.0, 12.0 * 3600.0, 1.0, 1.0 );

            boost::shared_ptr< EstimatableParameterSettings > perturbedParameterSettings;
            perturbedParameterSettings = (
                        boost::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                            2, 0, 2, 0, "Sun", spherical_harmonics_cosine_coefficient_block ) );

            std::pair< boost::shared_ptr< PodOutput< double > >, Eigen::VectorXd > estimationOutput =
                    determinePostfitParameterInfluence(
                        bodyMap, integratorSettings, propagatorSettings, perturbedParameterSettings,
                        6.0 * 3600.0, boost::assign::list_of( -sunNormalizedJ2 ), boost::assign::list_of( 0 ), 4 );

            input_output::writeMatrixToFile(
                        estimationOutput.first->residualHistory_.at( 0 ), "preFitResidualsSunJ2Body" +
                        boost::lexical_cast< std::string >( i ) + ".dat", 16, outputPath );
            input_output::writeMatrixToFile(
                        estimationOutput.first->residuals_, "postFitResidualsSunJ2Body" +
                        boost::lexical_cast< std::string >( i ) + ".dat", 16, outputPath );

            std::cout<<"Parameter difference "<<estimationOutput.second.transpose( )<<std::endl;
        }
    }
}
