// Copyright Contributors to the rawtoaces project.
// SPDX-License-Identifier: Apache-2.0
// https://github.com/AcademySoftwareFoundation/rawtoaces

#include <OpenImageIO/unittest.h>

#include "../src/core/spectral_solver.hpp"
#include "../src/core/matrix_solver.hpp"
#include <rawtoaces/define.h>
#include <rawtoaces/rta.h>

//#include "../src/rawtoaces_util2/data_adapters/spectral_data.hpp"

using namespace std;
using namespace rta;

const std::string data_path = "../../../data/";

void test_linear_solver()
{
    double A[3][3] = { { 2, 1, -1 }, { -3, -1, 2 }, { -2, 1, 2 } };

    double B[3] = { 8, -11, -3 };

    std::vector<std::vector<double>> a( 3 );
    std::vector<double>              b( 3 );

    for ( size_t i = 0; i < 3; i++ )
    {
        a[i] = std::vector<double>( A[i], A[i] + 3 );
        b[i] = B[i];
    }

    core::solve_linear( a, b );
    OIIO_CHECK_EQUAL_THRESH( b[0], 2, 1e-30 );
    OIIO_CHECK_EQUAL_THRESH( b[1], 3, 1e-30 );
    OIIO_CHECK_EQUAL_THRESH( b[2], -1, 1e-30 );
}

void check_camera(
    const std::string &camera_path,
    const std::string &camera_make,
    const std::string &camera_model,

    const std::string &illuminant_path,
    const std::string &illuminant_name )
{
    std::string observer_path = data_path + "cmf/cmf_1931.json";
    std::string training_path = data_path + "training/training_spectral.json";

    std::vector<double> old_wb;
    std::vector<double> new_wb;

    std::vector<std::vector<double>> old_mat;
    std::vector<std::vector<double>> new_mat;

    Illum illum;

    {
        Idt idt;
        idt.loadCameraSpst(
            camera_path, camera_make.c_str(), camera_model.c_str() );

        vector<string> illumPaths;
        illumPaths.push_back( illuminant_path );
        idt.loadIlluminant( illumPaths, illuminant_name );

        idt.loadTrainingData( training_path );
        idt.loadCMF( observer_path );

        idt.chooseIllumType( illuminant_name.c_str(), 0 );

        illum = idt.getBestIllum();

        idt.calIDT();
        old_wb  = idt.getWB();
        old_mat = idt.getIDT();
    }

    {
        rta::core::Spectrum illuminant_spectrum;
        illuminant_spectrum.values = illum.getIllumData();

        rta::core::SpectralData illuminant;
        illuminant.data["main"].push_back( illuminant_spectrum );

        util::SpectralLoader spectral_loader;

        auto camera = spectral_loader.load_spectral_file( camera_path );

        rta::core::SpectralSolver solver;
        solver.camera = camera;
        solver.scale_illuminant( illuminant.data["main"][0] );

        new_wb = solver.calculate_white_balance(
            illuminant.data["main"][0], solver.current_camera );

        rta::core::SpectralData cmf( observer_path, true );

        solver.current_illuminant = illuminant.data["main"][0];
        solver.white_balance      = new_wb;
        solver.colour_matching    = cmf.data["main"];

        rta::core::SpectralData training_data(
            data_path + "training/training_spectral.json" );

        solver.prepare_training_data( training_data.data["main"] );
        new_mat = solver.calculate_IDT();
    }

    for ( size_t i = 0; i < 3; i++ )
        OIIO_CHECK_EQUAL_THRESH( old_wb[i], new_wb[i], 1e-14 );

    for ( size_t i = 0; i < 3; i++ )
        for ( size_t j = 0; j < 3; j++ )
            OIIO_CHECK_EQUAL_THRESH( old_mat[i][j], new_mat[i][j], 1.2e-6 );
}

// The 2 commented out cameras don't work (the matrices produced by the two
// solvers don't match), due to a (possibly) incorrect scaling in the old
// solver. See the "TODO" comment in rta.cpp. With the proposed change in place
// all tests here pass successfully.
const char *cameras[][3] = {
    //    { "camera/arri_d21_380_780_5.json", "arri", "d21" },
    { "camera/canon_eos_5d_mark_ii_380_780_5.json", "canon", "eos 5d mark ii" },
    { "camera/canon_powershot_s90_380_780_5.json", "canon", "powershot s90" },
    { "camera/canon_xti_380_780_5.json", "canon", "xti" },
    { "camera/nikon_d70_380_780_5.json", "nikon", "d70" },
    { "camera/nikon_d200_380_780_5.json", "nikon", "d200" },
    //    { "camera/nikon_d5100_380_780_5.json", "nikon", "d5100" },
    { "camera/nikon_d700_380_780_5.json", "nikon", "d700" },
    { "camera/nikon_d7000_380_780_5.json", "nikon", "d7000" },
    { "camera/sony_ilce-7rm2_380_780_5.json", "sony", "ilce-7rm2" },
    { "camera/sony_ilce-7sm2_380_780_5.json", "sony", "ilce-7sm2" }
};

const char *illuminants[][2] = {
    { "illuminant/iso7589_stutung_380_780_5.json", "iso7589" },
};

void check_illuminant(
    const std::string &illuminant_path, const std::string &illuminant_name )
{
    std::cerr << "Checking illuminant " << illuminant_name << std::endl;

    for ( size_t c = 0; c < sizeof( cameras ) / sizeof( cameras[0] ); c++ )
    {

        const std::string camera_path  = data_path + cameras[c][0];
        const std::string camera_make  = cameras[c][1];
        const std::string camera_model = cameras[c][2];
        std::cerr << "   Checking camera " << camera_make << " " << camera_model
                  << std::endl;

        check_camera(
            camera_path,
            camera_make,
            camera_model,
            illuminant_path,
            illuminant_name );
    }
}

void verify_matrix()
{

    for ( size_t i = 0; i < sizeof( illuminants ) / sizeof( illuminants[0] );
          i++ )
    {
        const std::string illuminant_path = data_path + illuminants[i][0];
        const std::string illuminant_name = illuminants[i][1];

        check_illuminant( illuminant_path, illuminant_name );
    }

    // Blackbody - pre-calculate
    for ( int i = 1500; i < 4000; i += 500 )
    {
        const std::string illuminant_path = "";
        const std::string illuminant_name = to_string( i ) + "k";

        check_illuminant( illuminant_path, illuminant_name );
    }

    // Daylight - pre-calculate
    for ( int i = 4000; i <= 25000; i += 500 )
    {
        const std::string illuminant_path = "";
        const std::string illuminant_name = "d" + ( to_string( i / 100 ) );

        check_illuminant( illuminant_path, illuminant_name );
    }
}

int main( int, char ** )
{
    test_linear_solver();
    //    verify_matrix();
    return unit_test_failures;
}
