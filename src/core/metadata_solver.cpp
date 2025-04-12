// Copyright Contributors to the rawtoaces project.
// SPDX-License-Identifier: Apache-2.0
// https://github.com/AcademySoftwareFoundation/rawtoaces

#include <assert.h>
#include <iostream>
#include <map>

#include "metadata_solver.hpp"

#include "constants.hpp"
#include "matrix_solver.hpp"

// Needs to be below matrix_solver
#include <cmath>

namespace rta
{
namespace core
{

MetadataSolver::MetadataSolver( const Metadata &metadata )
    : _metadata( metadata )
{}

// Adapted from DNGIdt::lightSourceToColorTemp()
double CCT_from_lightsource_tag( unsigned short tag )
{
    if ( tag >= 32768 )
        return tag - 32768;

    const std::map<unsigned short, double> mapping = {
        { 2, 3500 },  { 3, 3400 },  { 10, 5550 }, { 17, 2856 },
        { 18, 4874 }, { 19, 6774 }, { 21, 6500 }, { 22, 7500 }
    };

    if ( mapping.contains( tag ) )
        return mapping.at( tag );
    return 5500.0;
}

// Adapted from DNGIdt::XYZtoCameraWeightedMatrix()
std::vector<std::vector<double>> XYZtoCameraWeightedMatrix(
    const double              &mir0,
    const double              &mir1,
    const double              &mir2,
    const std::vector<double> &matrix1,
    const std::vector<double> &matrix2,
    size_t                     rows )
{
    assert( matrix1.size() == matrix2.size() );

    size_t cols = matrix1.size() / rows;
    assert( matrix1.size() == rows * cols );

    double weight =
        std::max( 0.0, std::min( 1.0, ( mir1 - mir0 ) / ( mir1 - mir2 ) ) );

    std::vector<std::vector<double>> result(
        rows, std::vector<double>( cols ) );

    for ( size_t r = 0; r < rows; r++ )
        for ( size_t c = 0; c < cols; c++ )
            result[r][c] = matrix1[r * cols + c] * ( 1.0 - weight ) +
                           matrix2[r * cols + c] * weight;
    return result;
}

// Adapted from mathOps::XYZTouv
template <typename T>
std::vector<T> XYZ_to_CIE1960uv( const std::vector<T> &XYZ )
{
    T              scale  = XYZ[0] + 15 * XYZ[1] + 3 * XYZ[2];
    std::vector<T> result = { ( 4.0 * XYZ[0] ) / scale,
                              ( 6.0 * XYZ[1] ) / scale };
    return result;
};

template <typename T>
std::vector<T> CIE1960uv_to_XYZ( const std::vector<T> &uv )
{
    T scale = 2 * uv[0] - 8 * uv[1] + 4;
    T x     = 3.0 * uv[0] / scale;
    T y     = 2.0 * uv[0] / scale;
    T z     = 1.0 - x - y;

    std::vector<T> result = { x, y, z };
    return result;
};

// Adapted from DNGIdt::robertsonLength()
double
robertsonLength( const std::vector<double> &uv, const std::vector<double> &uvt )
{
    double              t    = uvt[2];
    double              sign = t < 0 ? -1.0 : t > 0 ? 1.0 : 0.0;
    std::vector<double> slope( 2 );
    slope[0] = -sign / std::sqrt( 1 + t * t );
    slope[1] = t * slope[0];

    std::vector<double> uvr( uvt.begin(), uvt.begin() + 2 );
    return slope[0] * ( uv[1] - uvr[1] ) - slope[1] * ( uv[0] - uvr[0] );
}

const std::vector<std::pair<double, std::tuple<double, double, double>>>
    RobertsonTable = {
        // clang-format off
    { 1e-10, { 0.18006,  0.26352,   -0.24341 }},
    { 10.0,  { 0.18066,  0.26589,   -0.25479 }},
    { 20.0,  { 0.18133,  0.26846,   -0.26876 }},
    { 30.0,  { 0.18208,  0.27119,   -0.28539 }},
    { 40.0,  { 0.18293,  0.27407,   -0.3047  }},
    { 50.0,  { 0.18388,  0.27709,   -0.32675 }},
    { 60.0,  { 0.18494,  0.28021,   -0.35156 }},
    { 70.0,  { 0.18611,  0.28342,   -0.37915 }},
    { 80.0,  { 0.18740,  0.28668,   -0.40955 }},
    { 90.0,  { 0.18880,  0.28997,   -0.44278 }},
    {100.0,  { 0.19032,  0.29326,   -0.47888 }},
    {125.0,  { 0.19462,  0.30141,   -0.58204 }},
    {150.0,  { 0.19962,  0.30921,   -0.70471 }},
    {175.0,  { 0.20525,  0.31647,   -0.84901 }},
    {200.0,  { 0.21142,  0.32312,   -1.0182  }},
    {225.0,  { 0.21807,  0.32909,   -1.2168  }},
    {250.0,  { 0.22511,  0.33439,   -1.4512  }},
    {275.0,  { 0.23247,  0.33904,   -1.7298  }},
    {300.0,  { 0.24010,  0.34308,   -2.0637  }},
    {325.0,  { 0.24792,  0.34655,   -2.4681  }},
    {350.0,  { 0.25591,  0.34951,   -2.9641  }},
    {375.0,  { 0.26400,  0.35200,   -3.5814  }},
    {400.0,  { 0.27218,  0.35407,   -4.3633  }},
    {425.0,  { 0.28039,  0.35577,   -5.3762  }},
    {450.0,  { 0.28863,  0.35714,   -6.7262  }},
    {475.0,  { 0.29685,  0.35823,   -8.5955  }},
    {500.0,  { 0.30505,  0.35907,  -11.324   }},
    {525.0,  { 0.31320,  0.35968,  -15.628   }},
    {550.0,  { 0.32129,  0.36011,  -23.325   }},
    {575.0,  { 0.32931,  0.36038,  -40.77    }},
    {600.0,  { 0.33724,  0.36051, -116.45    }}
        // clang-format on
    };

// Adapted from DNGIdt::XYZToColorTemperature()
double XYZToColorTemperature( const std::vector<double> &XYZ )
{

    std::vector<double> uv      = XYZ_to_CIE1960uv( XYZ );
    size_t              Nrobert = RobertsonTable.size();
    size_t              i;

    double mired;
    double RDthis = 0.0, RDprevious = 0.0;

    for ( i = 0; i < Nrobert; i++ )
    {
        auto                uvt       = RobertsonTable[i];
        std::vector<double> robertson = { std::get<0>( uvt.second ),
                                          std::get<1>( uvt.second ),
                                          std::get<2>( uvt.second ) };
        if ( ( RDthis = robertsonLength( uv, robertson ) ) <= 0.0 )
            break;
        RDprevious = RDthis;
    }
    if ( i <= 0 )
        mired = RobertsonTable[0].first;
    else if ( i >= Nrobert )
        mired = RobertsonTable[Nrobert - 1].first;
    else
        mired = RobertsonTable[i - 1].first +
                RDprevious *
                    ( RobertsonTable[i].first - RobertsonTable[i - 1].first ) /
                    ( RDprevious - RDthis );

    double cct = 1.0e06 / mired;
    cct        = std::max( 2000.0, std::min( 50000.0, cct ) );

    return cct;
}

// Adapted from DNGIdt::colorTemperatureToXYZ()
std::vector<double> colorTemperatureToXYZ( const double &cct )
{

    double              mired = 1.0e06 / cct;
    std::vector<double> uv( 2, 1.0 );

    size_t Nrobert = RobertsonTable.size();
    size_t i;

    for ( i = 0; i < Nrobert; i++ )
    {
        if ( RobertsonTable[i].first >= mired )
            break;
    }

    if ( i <= 0 )
    {
        auto [u, v, t] = RobertsonTable[0].second;
        uv             = { u, v, t };
    }
    else if ( i >= Nrobert )
    {
        auto [u, v, t] = RobertsonTable[Nrobert - 1].second;
        uv             = { u, v, t };
    }
    else
    {
        double weight =
            ( mired - RobertsonTable[i - 1].first ) /
            ( RobertsonTable[i].first - RobertsonTable[i - 1].first );

        auto [u1, v1, t1] = RobertsonTable[i].second;
        auto [u2, v2, t2] = RobertsonTable[i + 1].second;
        uv                = { u1 * weight + u2 * ( 1.0 - weight ),
                              v1 * weight + v2 * ( 1.0 - weight ),
                              t1 * weight + t2 * ( 1.0 - weight ) };
    }

    return CIE1960uv_to_XYZ( uv );
}

inline double kelvin_to_mired( double kelvin )
{
    return 1e06 / kelvin;
}

std::vector<std::vector<double>>
vec_to_mat( const std::vector<double> vec, size_t rows )
{
    size_t cols = vec.size() / rows;
    assert( vec.size() == rows * cols );

    std::vector<std::vector<double>> result(
        rows, std::vector<double>( cols ) );

    for ( size_t r = 0; r < rows; r++ )
        for ( size_t c = 0; c < cols; c++ )
            result[r][c] = vec[r * cols + c];
    return result;
}

std::vector<std::vector<double>> MetadataSolver::calculate_IDT() const
{
    std::vector<std::vector<double>> xyz_to_camera;

    // Adapted from DNGIdt::findXYZtoCameraMtx()
    {
        assert( _metadata.xyz2rgbMatrix1.size() > 0 );

        if ( _metadata.xyz2rgbMatrix2.size() == 0 )
        {
            std::cerr << " Only one calibration matrix found." << std::endl;
            xyz_to_camera = vec_to_mat(
                _metadata.xyz2rgbMatrix1, _metadata.xyz2rgbMatrix1.size() / 3 );
        }

        if ( _metadata.neutralRGB.size() == 0 )
        {
            std::cerr << " no neutral RGB values were found." << std::endl;
            xyz_to_camera = vec_to_mat(
                _metadata.xyz2rgbMatrix1, _metadata.xyz2rgbMatrix1.size() / 3 );
        }

        double cct1 =
            CCT_from_lightsource_tag( _metadata.calibrationIlluminant1 );
        double cct2 =
            CCT_from_lightsource_tag( _metadata.calibrationIlluminant2 );

        double mir1 = kelvin_to_mired( cct1 );
        double mir2 = kelvin_to_mired( cct2 );

        double maxMir = kelvin_to_mired( 2000.0 );
        double minMir = kelvin_to_mired( 50000.0 );

        double lomir =
            std::max( minMir, std::min( maxMir, std::min( mir1, mir2 ) ) );
        double himir =
            std::max( minMir, std::min( maxMir, std::max( mir1, mir2 ) ) );
        double mirStep = std::max( 5.0, ( himir - lomir ) / 50.0 );

        double mir = 0.0, lastMired = 0.0, estimatedMired = 0.0, lerror = 0.0,
               lastError = 0.0, smallestError = 0.0;

        for ( mir = lomir; mir < himir; mir += mirStep )
        {
            auto curr_xyz_to_camera = XYZtoCameraWeightedMatrix(
                mir,
                mir1,
                mir2,
                _metadata.xyz2rgbMatrix1,
                _metadata.xyz2rgbMatrix2,
                _metadata.xyz2rgbMatrix1.size() / 3 );
            std::vector<std::vector<double>> curr_camera_to_xyz;
            bool result = inverse( curr_xyz_to_camera, curr_camera_to_xyz );
            if ( !result )
                continue;

            auto white_xyz =
                product( curr_camera_to_xyz, _metadata.neutralRGB );
            auto cct = XYZToColorTemperature( white_xyz );

            lerror = mir - kelvin_to_mired( cct );

            if ( std::fabs( lerror - 0.0 ) <= 1e-09 )
            {
                estimatedMired = mir;
                break;
            }
            if ( std::fabs( mir - lomir - 0.0 ) > 1e-09 &&
                 lerror * lastError <= 0.0 )
            {
                estimatedMired = mir + ( lerror / ( lerror - lastError ) *
                                         ( mir - lastMired ) );
                break;
            }
            if ( std::fabs( mir - lomir ) <= 1e-09 ||
                 std::fabs( lerror ) < std::fabs( smallestError ) )
            {
                estimatedMired = mir;
                smallestError  = lerror;
            }

            lastError = lerror;
            lastMired = mir;
        }

        xyz_to_camera = XYZtoCameraWeightedMatrix(
            estimatedMired,
            mir1,
            mir2,
            _metadata.xyz2rgbMatrix1,
            _metadata.xyz2rgbMatrix2,
            _metadata.xyz2rgbMatrix1.size() / 3 );
    }

    std::vector<std::vector<double>> camera_to_xyz;
    std::vector<double>              camera_white_point_XYZ;

    // From DNGIdt::getCameraXYZMtxAndWhitePoint()
    {
        bool result = inverse( xyz_to_camera, camera_to_xyz );
        assert( result );

        double scale = std::pow( 2.0, _metadata.baselineExposure );

        for ( auto &i: camera_to_xyz )
            for ( auto &j: i )
                j *= scale;

        if ( _metadata.neutralRGB.size() > 0 )
        {
            camera_white_point_XYZ =
                product( camera_to_xyz, _metadata.neutralRGB );
        }
        else
        {
            camera_white_point_XYZ = colorTemperatureToXYZ(
                CCT_from_lightsource_tag( _metadata.calibrationIlluminant1 ) );
        }

        camera_white_point_XYZ[0] /= camera_white_point_XYZ[1];
        camera_white_point_XYZ[2] /= camera_white_point_XYZ[1];
        camera_white_point_XYZ[1] = 1.0;
    }

    auto CAT_matrix = calculate_CAT( camera_white_point_XYZ, ACES_white_XYZ );
    auto IDT_matrix = product( XYZ_D65_to_ACES, CAT_matrix );
    return IDT_matrix;
}

} // namespace core
} // namespace rta
