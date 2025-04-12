// Copyright Contributors to the rawtoaces project.
// SPDX-License-Identifier: Apache-2.0
// https://github.com/AcademySoftwareFoundation/rawtoaces

#include <OpenImageIO/unittest.h>
#include <vector>

#include "../src/core/constants.hpp"
#include "../src/core/matrix_solver.hpp"
#include "test_data.hpp"

void test_transpose()
{
    const auto &m1 = rta::core::XYZ_to_ACES_transposed;
    auto        m2 = rta::core::transposed( m1 );
    const auto &m3 = rta::core::XYZ_to_ACES;

    OIIO_CHECK_EQUAL( m2.size(), 4 );
    OIIO_CHECK_EQUAL( m2[0].size(), 4 );
    for ( size_t i = 0; i < 3; i++ )
        for ( size_t j = 0; j < 3; j++ )
            OIIO_CHECK_EQUAL_THRESH( m2[i][j], (float)m3[i][j], 1e-10 );
}

void test_product_mat_mat()
{
    const auto &m1 = rta::core::ACES_to_XYZ;
    const auto &m2 = rta::core::XYZ_to_ACES;
    auto        m3 = rta::core::product( m1, m2 );

    OIIO_CHECK_EQUAL( m3.size(), 3 );
    OIIO_CHECK_EQUAL( m3[0].size(), 3 );
    for ( size_t i = 0; i < 3; i++ )
        for ( size_t j = 0; j < 3; j++ )
            OIIO_CHECK_EQUAL_THRESH( m3[i][j], ( i == j ? 1.0 : 0.0 ), 1e-10 );
}

void test_product_mat_vec()
{
    const auto               &m1 = rta::core::ACES_to_XYZ;
    const std::vector<double> v( 3, 1 );
    auto                      v2 = rta::core::product( m1, v );

    OIIO_CHECK_EQUAL( v2.size(), 3 );
    for ( size_t i = 0; i < 3; i++ )
        OIIO_CHECK_EQUAL_THRESH( v2[i], rta::core::ACES_white_XYZ[i], 1e-15 );
}

void test_inverse()
{
    const auto                      &m1 = rta::core::ACES_to_XYZ;
    std::vector<std::vector<double>> m2;
    bool                             result = rta::core::inverse( m1, m2 );
    OIIO_CHECK_ASSERT( result );

    for ( size_t i = 0; i < 3; i++ )
        for ( size_t j = 0; j < 3; j++ )
            OIIO_CHECK_EQUAL_THRESH(
                m2[i][j], rta::core::XYZ_to_ACES[i][j], 1e-10 );
}

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

    rta::core::solve_linear( a, b );
    OIIO_CHECK_EQUAL_THRESH( b[0], 2, 1e-30 );
    OIIO_CHECK_EQUAL_THRESH( b[1], 3, 1e-30 );
    OIIO_CHECK_EQUAL_THRESH( b[2], -1, 1e-30 );
}

void test_CAT()
{
    // clang-format off
    const std::vector<double> src = {
        1.06789561546674380,
        1.00000000000000000,
        0.41004691570276436
    };
    
    const auto & dst = rta::core::ACES_white_XYZ;
    
    const std::vector<std::vector<double>> cat = {
        {  0.8908949453633443500, -0.113493669233048970, 0.27986294595392236},
        { -0.0822708087452159040,  1.042220088327919800, 0.11129591728250078},
        {  0.0057446406440006493,  0.019137946949782242, 2.39863421277044120}
    };
    // clang-format on

    auto result = rta::core::calculate_CAT( src, dst );

    OIIO_CHECK_EQUAL( result.size(), 3 );
    OIIO_CHECK_EQUAL( result[0].size(), 3 );

    for ( size_t i = 0; i < 3; i++ )
        for ( size_t j = 0; j < 3; j++ )
            OIIO_CHECK_EQUAL_THRESH( result[i][j], cat[i][j], 1e-6 );
}

int main( int, char ** )
{
    for ( auto use_eigen: { false, true } )
    {
        rta::core::use_eigen = use_eigen;

        test_transpose();
        test_product_mat_mat();
        test_product_mat_vec();
        test_inverse();
        test_linear_solver();

        test_CAT();
    }

    return unit_test_failures;
}
