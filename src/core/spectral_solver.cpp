// Copyright Contributors to the rawtoaces project.
// SPDX-License-Identifier: Apache-2.0
// https://github.com/AcademySoftwareFoundation/rawtoaces

#include <assert.h>

#include <rawtoaces/spectral_solver.hpp>
#include <rawtoaces/define.h>
#include <rawtoaces/mathOps.h>

namespace rta
{
namespace core
{

// clang-format off
static const double cat02_inv[3][3] = {
    {  1.096124, -0.278869,  0.182745 },
    {  0.454369,  0.473533,  0.072098 },
    { -0.009628, -0.005698,  1.015326 }
};
// clang-format on

void SpectralSolver::scale_illuminant( SpectralData &illuminant ) const
{
    assert( current_camera.data.contains( "main" ) );
    const auto &camera_data = current_camera.data.at( "main" );
    assert( camera_data.size() == 3 );
    assert( camera_data[0].shape.step > 0 );
    assert( camera_data[0].values.size() > 0 );

    assert( illuminant.data.contains( "main" ) );
    const auto &illuminant_data = illuminant.data.at( "main" );
    assert( illuminant_data.size() == 1 );
    assert( illuminant_data[0].shape.step > 0 );
    assert( illuminant_data[0].values.size() > 0 );

    assert( camera_data[0].values.size() == illuminant_data[0].values.size() );

    const double max_r = camera_data[0].max();
    const double max_g = camera_data[1].max();
    const double max_b = camera_data[2].max();

    size_t max_ch = 0;
    if ( max_g >= max_r && max_g >= max_b )
        max_ch = 1;
    else if ( max_b >= max_r && max_b >= max_g )
        max_ch = 2;

    double scale =
        1.0 / ( illuminant_data[0] * camera_data[max_ch] ).integrate();
    illuminant.data["main"][0] *= scale;
}

std::vector<double>
SpectralSolver::calculate_white_balance( SpectralData &illuminant ) const
{
    assert( current_camera.data.contains( "main" ) );
    const auto &camera_data = current_camera.data.at( "main" );
    assert( camera_data.size() == 3 );
    assert( camera_data[0].shape.step > 0 );
    assert( camera_data[0].values.size() > 0 );

    assert( illuminant.data.contains( "main" ) );
    const auto &illuminant_data = illuminant.data.at( "main" );
    assert( illuminant_data.size() == 1 );
    assert( illuminant_data[0].shape.step > 0 );
    assert( illuminant_data[0].values.size() > 0 );

    assert( camera_data[0].values.size() == illuminant_data[0].values.size() );

    double red   = ( camera_data[0] * illuminant_data[0] ).integrate();
    double green = ( camera_data[1] * illuminant_data[0] ).integrate();
    double blue  = ( camera_data[2] * illuminant_data[0] ).integrate();

    std::vector<double> result( 3, 0.0 );
    result[0] = green / red;
    result[1] = 1.0;
    result[2] = green / blue;
    return result;
}

std::vector<std::vector<double>> SpectralSolver::calculate_CAT(
    const std::vector<double> &src_white_XYZ,
    const std::vector<double> &dst_white_XYZ ) const
{
    assert( src_white_XYZ.size() == 3 );
    assert( dst_white_XYZ.size() == 3 );

    std::vector<std::vector<double>> v1( 3, std::vector<double>( 3 ) );
    std::vector<std::vector<double>> v2( 3, std::vector<double>( 3 ) );
    std::vector<std::vector<double>> v3( 3, std::vector<double>( 3 ) );

    for ( size_t i = 0; i < 3; i++ )
    {
        for ( size_t j = 0; j < 3; j++ )
        {
            v1[i][j] = cat02[i][j];
            v3[i][j] = cat02_inv[i][j];
        }
    }

    std::vector<double> src_white_LMS = mulVector( src_white_XYZ, v1 );
    std::vector<double> dst_white_LMS = mulVector( dst_white_XYZ, v1 );

    for ( size_t i = 0; i < 3; i++ )
        for ( size_t j = 0; j < 3; j++ )
            v2[i][j] = i == j ? dst_white_LMS[i] / src_white_LMS[i] : 0;

    auto x1 = mulVector( v2, transposeVec( v1 ) );
    auto x2 = mulVector( v3, transposeVec( x1 ) );
    return x2;
}

} // namespace core
} // namespace rta
