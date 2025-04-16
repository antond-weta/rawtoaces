// Copyright Contributors to the rawtoaces project.
// SPDX-License-Identifier: Apache-2.0
// https://github.com/AcademySoftwareFoundation/rawtoaces

#include <OpenImageIO/unittest.h>

#include "../src/core/spectral_solver.hpp"
#include "../src/core/matrix_solver.hpp"
#include "../src/core/constants.hpp"
#include "../src/util/data_adapters/spectral_loader.hpp"

#include "test_data.hpp"

using namespace std;
using namespace rta;

#if defined( WIN32 ) || defined( WIN64 )
const std::string separator = ";";

int setenv( const char *name, const char *value, int overwrite )
{
    return _putenv_s( name, value );
}

#else
const std::string separator = ":";
#endif

const std::string data_path = "../../../data/" + separator + "../../data/";

const int verbosity = 2;

void check_dimensions( const rta::core::SpectralData &data, size_t c )
{
    std::cerr << "check_dimensions" << std::endl;
    const auto &d = data.data;
    OIIO_CHECK_EQUAL( d.size(), 1 );

    const auto &m = d.at( "main" );
    OIIO_CHECK_EQUAL( m.size(), c );

    for ( size_t i = 0; i < c; i++ )
    {
        const auto &s = m[i];
        OIIO_CHECK_EQUAL( s.shape.first, 380.0 );
        OIIO_CHECK_EQUAL( s.shape.last, 780.0 );
        OIIO_CHECK_EQUAL( s.shape.step, 5.0 );
        OIIO_CHECK_EQUAL( s.values.size(), 81 );
    }
}

// Check if the camera contains the Nikon D200 sensitivity data.
void check_camera( const rta::core::SpectralData &camera )
{
    std::cerr << "check_camera" << std::endl;
    check_dimensions( camera, 3 );

    for ( size_t i = 0; i < 3; i++ )
    {
        const auto &s = camera.data.at( "main" )[i];
        for ( size_t j = 0; j < 81; j++ )
            OIIO_CHECK_EQUAL_THRESH(
                s.values[j], rta::test::nikon_d200_sensitivity[j][i], 1e-30 );
    }
}

// Check if the illuminant contains the ISO7589 studio tungsten light spectrum.
void check_illuminant( const rta::core::SpectralData &illuminant )
{
    std::cerr << "check_illuminant" << std::endl;
    check_dimensions( illuminant, 1 );

    const auto &s = illuminant.data.at( "main" )[0];

    for ( size_t j = 0; j < 81; j++ )
        OIIO_CHECK_EQUAL_THRESH(
            s.values[j], rta::test::illuminant_ISO7589_SPD[j], 1e-30 );
}

// Check if the cmf object contains the CIE-1931 XYZ colour matching functions.
// The data file contains 360:830:1 samples, so the code also checks
// if resampling to 380:780:5 works correctly.
void check_cmf( const rta::core::SpectralData &cmf )
{
    std::cerr << "check_cmf" << std::endl;
    check_dimensions( cmf, 3 );

    for ( size_t i = 0; i < 3; i++ )
    {
        const auto &s = cmf.data.at( "main" )[i];
        for ( size_t j = 0; j < 81; j++ )
            OIIO_CHECK_EQUAL_THRESH(
                s.values[j], rta::test::colour_matching_XYZ[j][i], 1e-30 );
    }
}

// Check if the illuminant contains the ISO7589 studio tungsten light spectrum
// scaled to the Nikon D200 curves.
void check_scale_illuminant( const rta::core::SpectralData &illuminant )
{
    std::cerr << "check_scale_illuminant" << std::endl;
    check_dimensions( illuminant, 1 );

    double sum = 0;
    for ( size_t i = 0; i < 81; i++ )
    {
        sum += rta::test::nikon_d200_sensitivity[i][1] *
               rta::test::illuminant_ISO7589_SPD[i];
    }

    double scale = 1 / sum;
    OIIO_CHECK_EQUAL_THRESH( scale, 0.13655488147514236, 1e-30 );

    const auto &si = illuminant.data.at( "main" )[0];
    for ( size_t i = 0; i < 81; i++ )
        OIIO_CHECK_EQUAL_THRESH(
            si.values[i], rta::test::illuminant_ISO7589_SPD[i] * scale, 1e-16 );
}

void test_daylight()
{
    std::cerr << "test_daylight" << std::endl;
    core::SpectralSolver solver;
    core::SpectralData   illuminant;

    util::SpectralLoader spectral_loader;
    spectral_loader.verbosity = verbosity;

    auto camera_path =
        spectral_loader.find_file( "misc/CIE_illum_Dxx_comp.json" );
    auto components = spectral_loader.load_spectral_file( camera_path, true );

    solver.create_daylight_illuminant( 65, components, illuminant );
    auto values = illuminant.data["main"][0].values;

    OIIO_CHECK_EQUAL( values.size(), 81 );
    for ( size_t i = 0; i < 81; i++ )
    {
        OIIO_CHECK_EQUAL_THRESH( values[i], rta::test::D65[i], 0.015 );
    }
}

void test_best_illuminant()
{
    std::cerr << "test_best_illuminant" << std::endl;
    const struct
    {
        const char   *illum_name;
        const double *illum_data;
    } tests[] = { { "D65", rta::test::D65 },
                  { "3200K", rta::test::blackbody_3200K },
                  { "iso7589", rta::test::illuminant_ISO7589_SPD } };

    util::SpectralLoader spectral_loader;
    spectral_loader.verbosity = verbosity;

    auto camera_path =
        spectral_loader.find_file( "camera/nikon_d200_380_780_5.json" );
    auto camera = spectral_loader.load_spectral_file( camera_path );
    check_camera( camera );

    for ( size_t i = 0; i < sizeof( tests ) / sizeof( tests[0] ); i++ )
    {
        const double       *illuminant_data = tests[i].illum_data;
        std::vector<double> wb( 3, 0 );

        for ( size_t j = 0; j < 81; j++ )
        {
            for ( size_t k = 0; k < 3; k++ )
            {
                wb[k] += illuminant_data[j] *
                         rta::test::nikon_d200_sensitivity[j][k];
            }
        }

        wb[0] = wb[1] / wb[0];
        wb[2] = wb[1] / wb[2];
        wb[1] = 1.0;

        auto illuminant =
            spectral_loader.find_best_illuminant( camera, wb, 100 );
        OIIO_CHECK_EQUAL( illuminant.illuminant, tests[i].illum_name );
    }
}

void test_CAT()
{
    // clang-format off
    const double src_wb[3] = {
        1.06789561546674380,
        1.00000000000000000,
        0.41004691570276436
    };
    
    const double dst_wb[3] = {
        0.95264607456984595,
        1.00000000000000000,
        1.00882518435159010
    };
    
    const double cat[3][3] = {
        {  0.8908949453633443500, -0.113493669233048970, 0.27986294595392236},
        { -0.0822708087452159040,  1.042220088327919800, 0.11129591728250078},
        {  0.0057446406440006493,  0.019137946949782242, 2.39863421277044120}
    };
    // clang-format on

    std::vector<double> src( src_wb, src_wb + 3 );
    auto               &dst = core::ACES_white_XYZ;

    auto result = core::calculate_CAT( src, dst );

    OIIO_CHECK_EQUAL( result.size(), 3 );
    OIIO_CHECK_EQUAL( result[0].size(), 3 );

    for ( size_t i = 0; i < 3; i++ )
        for ( size_t j = 0; j < 3; j++ )
            OIIO_CHECK_EQUAL_THRESH( result[i][j], cat[i][j], 1e-6 );
}

void test_solver()
{
    util::SpectralLoader spectral_loader;
    spectral_loader.verbosity = verbosity;

    auto camera_path =
        spectral_loader.find_file( "camera/nikon_d200_380_780_5.json" );
    auto camera = spectral_loader.load_spectral_file( camera_path );
    check_camera( camera );

    auto illuminant_path = spectral_loader.find_file(
        "illuminant/iso7589_stutung_380_780_5.json" );
    auto illuminant = spectral_loader.load_spectral_file( illuminant_path );
    check_illuminant( illuminant );

    auto observer_path = spectral_loader.find_file( "cmf/cmf_1931.json" );
    auto observer = spectral_loader.load_spectral_file( observer_path, true );
    check_cmf( observer );

    auto training_path =
        spectral_loader.find_file( "training/training_spectral.json" );
    auto training_data = spectral_loader.load_spectral_file( training_path );

    rta::core::SpectralSolver solver;
    solver.camera = camera;
    solver.scale_illuminant( illuminant );
    check_scale_illuminant( illuminant );

    std::cerr << "solver.calculate_white_balance" << std::endl;

    auto wb = solver.calculate_white_balance( illuminant, camera );
    OIIO_CHECK_EQUAL( wb.size(), 3 );
    OIIO_CHECK_EQUAL_THRESH( wb[0], 1.1397265403538983, 1e-30 );
    OIIO_CHECK_EQUAL( wb[1], 1 );
    OIIO_CHECK_EQUAL_THRESH( wb[2], 2.3240151642033657, 1e-30 );

    solver.verbosity     = verbosity;
    solver.illuminant    = illuminant;
    solver.white_balance = wb;
    solver.observer      = observer;
    solver.training_data = training_data;

    //    solver.scale_illuminant(solver.illuminant);

    std::cerr << "solver.prepare_training_data" << std::endl;
    solver.prepare_training_data();

    std::cerr << "solver.calculate_training_XYZ" << std::endl;
    auto XYZ = solver.calculate_training_XYZ();

    std::cerr << "solver.calculate_training_RGB" << std::endl;
    auto RGB = solver.calculate_training_RGB();

    std::cerr << "solver.calculate_IDT" << std::endl;
    std::vector<std::vector<double>> IDT;
    std::cerr << "solver.calculate_IDT1" << std::endl;
    bool result = solver.calculate_IDT( IDT );

    std::cerr << "solver.calculate_IDT2" << std::endl;
    double true_IDT[3][3] = { { 0.744769240, 0.143420257, 0.111810501 },
                              { 0.045175883, 1.008261897, -0.053437780 },
                              { 0.024714133, -0.124552527, 1.099838394 } };

    OIIO_CHECK_ASSERT( result );
    OIIO_CHECK_EQUAL( IDT.size(), 3 );
    for ( size_t i = 0; i < 3; i++ )
    {
        OIIO_CHECK_EQUAL( IDT[i].size(), 3 );
        for ( size_t j = 0; j < 3; j++ )
            OIIO_CHECK_EQUAL_THRESH( IDT[i][j], true_IDT[i][j], 1.3e-9 );
    }
}

int main( int, char ** )
{
    setenv( "RAWTOACES_DATA_PATH", data_path.c_str(), true );

    test_daylight();
    test_best_illuminant();

    for ( auto use_eigen: { false, true } )
    {
        rta::core::use_eigen = use_eigen;

        for ( auto use_ceres: { false, true } )
        {
            rta::core::use_ceres = use_ceres;

            test_solver();
        }
    }

    return unit_test_failures;
}
