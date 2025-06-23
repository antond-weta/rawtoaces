// Copyright Contributors to the rawtoaces project.
// SPDX-License-Identifier: Apache-2.0
// https://github.com/AcademySoftwareFoundation/rawtoaces

#include <OpenImageIO/unittest.h>
#include <vector>

#include <rawtoaces/core/matrix_solver.hpp>
#include <rawtoaces/core/constants.hpp>

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

int main( int, char ** )
{
    rta::core::use_eigen = false;
    test_transpose();
    return unit_test_failures;
}
