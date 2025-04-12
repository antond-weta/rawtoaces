// Copyright Contributors to the rawtoaces project.
// SPDX-License-Identifier: Apache-2.0
// https://github.com/AcademySoftwareFoundation/rawtoaces

#include <OpenImageIO/unittest.h>
#include <vector>

#include "../src/core/constants.hpp"
#include "../src/core/matrix_solver.hpp"
#include "../src/core/spectral_solver.hpp"
#include "test_data.hpp"

void test_blackbody()
{
    rta::core::SpectralSolver solver;
    rta::core::SpectralData   illuminant;

    solver.create_blackbody_illuminant( 3200, illuminant );
    auto values = illuminant.data["main"][0].values;

    for ( size_t i = 0; i < 81; i++ )
    {
        OIIO_CHECK_EQUAL_THRESH(
            values[i], rta::test::blackbody_3200K[i] * 1e12, 10 );
    }
}

int main( int, char ** )
{
    for ( auto use_eigen: { false, true } )
    {
        rta::core::use_eigen = use_eigen;

        test_blackbody();
    }

    return unit_test_failures;
}
