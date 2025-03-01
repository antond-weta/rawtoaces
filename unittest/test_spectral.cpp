// Copyright Contributors to the rawtoaces project.
// SPDX-License-Identifier: Apache-2.0
// https://github.com/AcademySoftwareFoundation/rawtoaces

#include <OpenImageIO/unittest.h>

#include <rawtoaces/spectral_solver.hpp>
#include <rawtoaces/spectral_data.hpp>

using namespace std;
using namespace rta;

const std::string data_path = "../../data/";

// clang-format off

const double nikon_d200_sensitivity[81][3] = {
    { 0.0243, 0.0296, 0.0477 }, { 0.0160, 0.0201, 0.0336 },
    { 0.0078, 0.0106, 0.0195 }, { 0.0062, 0.0086, 0.0152 },
    { 0.0045, 0.0065, 0.0108 }, { 0.0043, 0.0059, 0.0086 },
    { 0.0040, 0.0052, 0.0064 }, { 0.0056, 0.0073, 0.0536 },
    { 0.0073, 0.0095, 0.1008 }, { 0.0175, 0.0234, 0.2714 },
    { 0.0277, 0.0373, 0.4419 }, { 0.0304, 0.0572, 0.6097 },
    { 0.0332, 0.0771, 0.7776 }, { 0.0300, 0.0906, 0.8411 },
    { 0.0268, 0.1040, 0.9046 }, { 0.0255, 0.1311, 0.9340 },
    { 0.0243, 0.1583, 0.9634 }, { 0.0257, 0.2160, 0.9511 },
    { 0.0270, 0.2737, 0.9389 }, { 0.0277, 0.3003, 0.8601 },
    { 0.0285, 0.3269, 0.7812 }, { 0.0280, 0.3384, 0.7242 },
    { 0.0275, 0.3498, 0.6672 }, { 0.0287, 0.4458, 0.5531 },
    { 0.0298, 0.5417, 0.4391 }, { 0.0343, 0.6403, 0.3571 },
    { 0.0388, 0.7388, 0.2751 }, { 0.0523, 0.8374, 0.1973 },
    { 0.0659, 0.9359, 0.1195 }, { 0.0715, 0.9680, 0.0807 },
    { 0.0771, 1.0000, 0.0419 }, { 0.0589, 0.9779, 0.0306 },
    { 0.0406, 0.9559, 0.0194 }, { 0.0301, 0.9072, 0.0137 },
    { 0.0196, 0.8584, 0.0081 }, { 0.0196, 0.7891, 0.0058 },
    { 0.0197, 0.7198, 0.0035 }, { 0.0427, 0.6480, 0.0027 },
    { 0.0658, 0.5762, 0.0020 }, { 0.2576, 0.5059, 0.0021 },
    { 0.4494, 0.4355, 0.0022 }, { 0.5820, 0.3596, 0.0026 },
    { 0.7146, 0.2836, 0.0030 }, { 0.7150, 0.2226, 0.0028 },
    { 0.7155, 0.1617, 0.0026 }, { 0.6795, 0.1179, 0.0025 },
    { 0.6436, 0.0741, 0.0023 }, { 0.6070, 0.0583, 0.0021 },
    { 0.5703, 0.0426, 0.0019 }, { 0.5215, 0.0340, 0.0027 },
    { 0.4727, 0.0254, 0.0036 }, { 0.4253, 0.0218, 0.0033 },
    { 0.3779, 0.0182, 0.0031 }, { 0.3227, 0.0148, 0.0034 },
    { 0.2674, 0.0113, 0.0037 }, { 0.2406, 0.0105, 0.0033 },
    { 0.2137, 0.0097, 0.0029 }, { 0.1695, 0.0087, 0.0029 },
    { 0.1253, 0.0078, 0.0028 }, { 0.0929, 0.0065, 0.0023 },
    { 0.0604, 0.0053, 0.0019 }, { 0.0368, 0.0042, 0.0022 },
    { 0.0131, 0.0032, 0.0024 }, { 0.0088, 0.0029, 0.0031 },
    { 0.0044, 0.0026, 0.0038 }, { 0.0040, 0.0024, 0.0033 },
    { 0.0037, 0.0022, 0.0027 }, { 0.0035, 0.0023, 0.0033 },
    { 0.0033, 0.0023, 0.0039 }, { 0.0028, 0.0020, 0.0034 },
    { 0.0023, 0.0017, 0.0029 }, { 0.0018, 0.0014, 0.0024 },
    { 0.0014, 0.0011, 0.0019 }, { 0.0009, 0.0008, 0.0014 },
    { 0.0004, 0.0005, 0.0009 }, { 0.0005, 0.0006, 0.0009 },
    { 0.0005, 0.0006, 0.0010 }, { 0.0006, 0.0007, 0.0010 },
    { 0.0007, 0.0008, 0.0011 }, { 0.0007, 0.0008, 0.0011 },
    { 0.0008, 0.0009, 0.0012 }
};

const double illuminant_ISO7589_SPD[81] = {
    0.0400, 0.0500, 0.0600, 0.0700, 0.0800, 0.0900, 0.1000, 0.1100,
    0.1200, 0.1325, 0.1450, 0.1575, 0.1700, 0.1800, 0.1900, 0.2025,
    0.2150, 0.2275, 0.2400, 0.2525, 0.2650, 0.2800, 0.2950, 0.3075,
    0.3200, 0.3350, 0.3500, 0.3650, 0.3800, 0.3925, 0.4050, 0.4225,
    0.4400, 0.4550, 0.4700, 0.4850, 0.5000, 0.5125, 0.5250, 0.5400,
    0.5550, 0.5675, 0.5800, 0.5950, 0.6100, 0.6225, 0.6350, 0.6475,
    0.6600, 0.6750, 0.6900, 0.7025, 0.7150, 0.7275, 0.7400, 0.7525,
    0.7650, 0.7750, 0.7850, 0.7975, 0.8100, 0.8225, 0.8350, 0.8475,
    0.8600, 0.8700, 0.8800, 0.8900, 0.9000, 0.9100, 0.9200, 0.9275,
    0.9350, 0.9450, 0.9550, 0.9650, 0.9750, 0.9800, 0.9850, 0.9925,
    1.0000
};

// clang-format on

void check_dimensions( const rta::core::SpectralData &data, size_t c )
{
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
    check_dimensions( camera, 3 );

    for ( size_t i = 0; i < 3; i++ )
    {
        const auto &s = camera.data.at( "main" )[i];
        for ( size_t j = 0; j < 81; j++ )
            OIIO_CHECK_EQUAL_THRESH(
                s.values[j], nikon_d200_sensitivity[j][i], 1e-30 );
    }
}

// Check if the illuminant contains the ISO7589 studio tungsten light spectrum.
void check_illuminant( const rta::core::SpectralData &illuminant )
{
    check_dimensions( illuminant, 1 );

    const auto &s = illuminant.data.at( "main" )[0];

    for ( size_t j = 0; j < 81; j++ )
        OIIO_CHECK_EQUAL_THRESH(
            s.values[j], illuminant_ISO7589_SPD[j], 1e-30 );
}

// Check if the illuminant contains the ISO7589 studio tungsten light spectrum
// scaled to the Nikon D200 curves.
void check_scale_illuminant( const rta::core::SpectralData &illuminant )
{
    check_dimensions( illuminant, 1 );

    double sum = 0;
    for ( size_t i = 0; i < 81; i++ )
    {
        sum += nikon_d200_sensitivity[i][1] * illuminant_ISO7589_SPD[i];
    }

    double scale = 1 / sum;
    OIIO_CHECK_EQUAL_THRESH( scale, 0.13655488147514236, 1e-30 );

    const auto &si = illuminant.data.at( "main" )[0];
    for ( size_t i = 0; i < 81; i++ )
        OIIO_CHECK_EQUAL_THRESH(
            si.values[i], illuminant_ISO7589_SPD[i] * scale, 1e-16 );
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

    std::vector<double> src( 3 );
    std::vector<double> dst( 3 );

    src.assign( src_wb, src_wb + 3 );
    dst.assign( dst_wb, dst_wb + 3 );

    rta::core::SpectralSolver solver;
    auto                      result = solver.calculate_CAT( src, dst );

    OIIO_CHECK_EQUAL( result.size(), 3 );
    OIIO_CHECK_EQUAL( result[0].size(), 3 );

    for ( size_t i = 0; i < 3; i++ )
        for ( size_t j = 0; j < 3; j++ )
            OIIO_CHECK_EQUAL_THRESH( result[i][j], cat[i][j], 1e-6 );
}

void test_solver()
{
    rta::core::SpectralData camera(
        data_path + "camera/nikon_d200_380_780_5.json" );
    check_camera( camera );

    rta::core::SpectralData illuminant(
        data_path + "illuminant/iso7589_stutung_380_780_5.json" );
    check_illuminant( illuminant );

    rta::core::SpectralSolver solver;
    solver.current_camera = camera;
    solver.scale_illuminant( illuminant );
    check_scale_illuminant( illuminant );

    auto wb = solver.calculate_white_balance( illuminant );
    OIIO_CHECK_EQUAL_THRESH( wb[0], 1.1397265403538983, 1e-30 );
    OIIO_CHECK_EQUAL( wb[1], 1 );
    OIIO_CHECK_EQUAL_THRESH( wb[2], 2.3240151642033657, 1e-30 );
}

int main( int, char ** )
{
    test_CAT();
    test_solver();
    return unit_test_failures;
}
