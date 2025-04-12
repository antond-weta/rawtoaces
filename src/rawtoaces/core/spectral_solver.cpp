// Copyright Contributors to the rawtoaces project.
// SPDX-License-Identifier: Apache-2.0
// https://github.com/AcademySoftwareFoundation/rawtoaces

#include <assert.h>
#include <iostream>
#include <string>

// This is needed for Visual Studio to expose M_PI
#define _USE_MATH_DEFINES
#include <math.h>
#include <cmath>

#include "constants.hpp"
#include "spectral_solver.hpp"
#include "matrix_solver.hpp"

namespace rta
{
namespace core
{

void SpectralSolver::create_blackbody_illuminant(
    float temp, SpectralData &out ) const
{
    // Planck's law
    // L = c1 / lambda^5 / (e^(c2 / lambda / T) - 1)
    // L - spectral radiant exitance
    // c1 = 2 * pi * h * c
    // c2 = h * c / k
    // h = 6.62607015×10−34 - Planck constant
    // v - frequency
    // c = 299792458 - speed of light in vacuum
    // k = 1.380649×10−23 - Bolzmann constant
    // T - temperature

    out.illuminant = "Blackbody " + std::to_string( int( temp ) ) + "K";

    if ( !out.data.contains( "main" ) )
    {
        out.data["main"].emplace_back();
    }

    auto  &shape  = out.data["main"][0].shape;
    auto  &values = out.data["main"][0].values;
    size_t i      = 0;
    for ( float wl = shape.first; wl <= shape.last; wl += shape.step )
    {
        constexpr double h = 6.62607015e-34;
        constexpr double c = 2.99792458e8;
        constexpr double k = 1.380649e-23;

        double lambda = wl / 1e9;
        double c1     = M_PI * 2.0 * h * c * c;
        double c2     = h * c / k;
        double L =
            c1 / std::pow( lambda, 5 ) / ( std::exp( c2 / lambda / temp ) - 1 );
        values[i++] = L;
    }
}

void SpectralSolver::create_daylight_illuminant(
    int temp, const SpectralData &components, SpectralData &out ) const
{
    out.illuminant = "Daylight D" + std::to_string( temp );

    if ( temp < 40 )
    {
        std::cerr
            << "Daylight illuminant should not be used with CCT below 4000K. "
            << temp * 100 << "K has been requested. "
            << "Please use the Blackbody illuminant instead." << std::endl;
    }

    if ( temp > 250 )
    {
        std::cerr << "Daylight illuminant is only defined up to 25000K. "
                  << temp * 100 << "K has been requested. "
                  << "The calculations may be inaccurate." << std::endl;
    }

    double cct = (double)temp * 100.0 * 1.438776877 / 1.4380;

    const double factors[2][4] = { { 0.244063, 0.09911, 2.9678, -4.6070 },
                                   { 0.237040, 0.24748, 1.9018, -2.0064 } };

    size_t ind = cct <= 7000.0 ? 0 : 1;

    double xd;
    double tmp = 1.0;
    for ( size_t i = 0; i < 4; i++ )
    {
        xd += factors[ind][i] * tmp;
        tmp *= 1000.0 / cct;
    }

    double yd = -3.0 * xd * xd + 2.870 * xd - 0.275;

    double M  = 0.0241 + 0.2562 * xd - 0.7341 * yd;
    double M1 = ( -1.3515 - 1.7703 * xd + 5.9114 * yd ) / M;
    double M2 = ( 0.03 - 31.4424 * xd + 30.0717 * yd ) / M;

    auto &cmp      = components.data.at( "main" );
    auto &spectrum = out.data["main"].emplace_back();
    spectrum       = cmp[0] + cmp[1] * M1 + cmp[2] * M2;
}

void SpectralSolver::scale_illuminant( SpectralData &illuminant ) const
{
    auto &cd = camera.data.at( "main" );
    assert( cd.size() == 3 );
    assert( cd[0].shape.step > 0 );
    assert( cd[0].values.size() > 0 );

    auto &id = illuminant.data.at( "main" );
    assert( id.size() == 1 );
    assert( id[0].shape.step > 0 );
    assert( id[0].values.size() > 0 );

    assert( cd[0].values.size() == id[0].values.size() );

    const double max_r = cd[0].max();
    const double max_g = cd[1].max();
    const double max_b = cd[2].max();

    size_t max_ch = 0;
    if ( max_g >= max_r && max_g >= max_b )
        max_ch = 1;
    else if ( max_b >= max_r && max_b >= max_g )
        max_ch = 2;

    double scale = 1.0 / ( id[0] * cd[max_ch] ).integrate();
    id[0] *= scale;
}

std::vector<double> SpectralSolver::calculate_white_balance(
    const SpectralData &illuminant, const SpectralData &observer ) const
{
    auto &id = illuminant.data.at( "main" );
    assert( id.size() == 1 );
    assert( id[0].shape.step > 0 );
    assert( id[0].values.size() > 0 );

    auto &od = observer.data.at( "main" );
    assert( od.size() == 3 );
    assert( od[0].shape.step > 0 );
    assert( od[0].values.size() > 0 );

    assert( od[0].values.size() == id[0].values.size() );

    double red   = ( od[0] * id[0] ).integrate();
    double green = ( od[1] * id[0] ).integrate();
    double blue  = ( od[2] * id[0] ).integrate();

    std::vector<double> result( 3, 0.0 );
    result[0] = green / red;
    result[1] = 1.0;
    result[2] = green / blue;
    return result;
}

void SpectralSolver::prepare_training_data()
{
    // Preparing training data requires a valid illuminant.
    auto &id = illuminant.data.at( "main" );
    assert( id.size() == 1 );
    assert( id[0].shape.step > 0 );
    assert( id[0].values.size() > 0 );

    auto &td = training_data.data.at( "main" );
    assert( td.size() > 0 );
    assert( td[0].shape.step > 0 );
    assert( td[0].values.size() > 0 );

    assert( td[0].values.size() == id[0].values.size() );

    // Multiply each patch spectral transittance curve with the illuminant SPD
    // to get spectral radiance.
    for ( auto &i: td )
    {
        i *= id[0];
    }
}

std::vector<std::vector<double>> SpectralSolver::calculate_training_XYZ() const
{
    auto &id = illuminant.data.at( "main" );
    assert( id.size() == 1 );
    assert( id[0].shape.step > 0 );
    assert( id[0].values.size() > 0 );

    auto &od = observer.data.at( "main" );
    assert( od.size() == 3 );
    assert( od[0].shape.step > 0 );
    assert( od[0].values.size() > 0 );

    auto &td = training_data.data.at( "main" );
    assert( td.size() > 0 );
    assert( td[0].shape.step > 0 );
    assert( td[0].values.size() > 0 );

    assert( td[0].values.size() == id[0].values.size() );
    assert( td[0].values.size() == od[0].values.size() );

    // Source white point
    auto src_white = calculate_white_balance( illuminant, observer );
    src_white[0]   = 1.0 / src_white[0];
    src_white[2]   = 1.0 / src_white[2];

    // Destination white point
    const std::vector<double> &dst_white = ACES_white_XYZ;

    // Chromatic adaptation matrix
    auto cat = calculate_CAT( src_white, dst_white );

    // Scale the colour matching functions to the illuminant
    double scale = 1.0 / ( id[0] * od[1] ).integrate();

    std::vector<std::vector<double>> result(
        td.size(), std::vector<double>( 3 ) );
    for ( size_t i = 0; i < td.size(); i++ )
    {
        result[i][0] = ( td[i] * od[0] ).integrate() * scale;
        result[i][1] = ( td[i] * od[1] ).integrate() * scale;
        result[i][2] = ( td[i] * od[2] ).integrate() * scale;
        result[i]    = product( cat, result[i] );
    }

    return result;
}

std::vector<std::vector<double>> SpectralSolver::calculate_training_RGB() const
{
    auto &cd = camera.data.at( "main" );
    assert( cd.size() == 3 );
    assert( cd[0].shape.step > 0 );
    assert( cd[0].values.size() > 0 );

    auto &td = training_data.data.at( "main" );
    assert( td.size() > 0 );
    assert( td[0].shape.step > 0 );
    assert( td[0].values.size() > 0 );

    assert( td[0].values.size() == cd[0].values.size() );

    std::vector<std::vector<double>> result(
        td.size(), std::vector<double>( 3 ) );

    for ( size_t i = 0; i < td.size(); i++ )
    {
        result[i][0] = ( td[i] * cd[0] ).integrate() * white_balance[0];
        result[i][1] = ( td[i] * cd[1] ).integrate() * white_balance[1];
        result[i][2] = ( td[i] * cd[2] ).integrate() * white_balance[2];
    }

    return result;
}

bool SpectralSolver::calculate_IDT(
    std::vector<std::vector<double>> &out ) const
{
    auto                      src_points = calculate_training_RGB();
    auto                      xyz_points = calculate_training_XYZ();
    const std::vector<double> white_rgb  = { 1, 1, 1 };
    return solve_IDT_LAB(
        src_points, xyz_points, white_rgb, ACES_to_XYZ, out, verbosity );
}

} // namespace core
} // namespace rta
