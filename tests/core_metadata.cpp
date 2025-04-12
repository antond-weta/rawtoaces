// Copyright Contributors to the rawtoaces project.
// SPDX-License-Identifier: Apache-2.0
// https://github.com/AcademySoftwareFoundation/rawtoaces

#include <OpenImageIO/unittest.h>
#include <vector>

#include "../src/core/constants.hpp"
#include "../src/core/metadata.hpp"
#include "../src/core/matrix_solver.hpp"
#include "../src/core/metadata_solver.hpp"
#include "test_data.hpp"

void test_metadata_IDT()
{
    rta::core::Metadata metadata;

    metadata.baselineExposure       = 2.4;
    metadata.calibrationIlluminant1 = 17;
    metadata.calibrationIlluminant2 = 21;

    metadata.neutralRGB = { 0.6289999865031245, 1, 0.79040003045288199 };

    metadata.xyz2rgbMatrix1 = { 1.3119699954986572,   -0.49678999185562134,
                                0.011559999547898769, -0.41723001003265381,
                                1.4423700571060181,   0.045279998332262039,
                                0.067230001091957092, 0.21709999442100525,
                                0.72650998830795288 };

    metadata.xyz2rgbMatrix2 = { 1.0088499784469604,    -0.27351000905036926,
                                -0.082580000162124634, -0.48996999859809875,
                                1.3444099426269531,    0.11174000054597855,
                                -0.064060002565383911, 0.32997000217437744,
                                0.5391700267791748 };

    rta::core::MetadataSolver solver( metadata );
    auto                      idt = solver.calculate_IDT();

    const double true_idt[3][3] = {
        { 1.0536467571155164, 0.0039044966623301416, 0.0049082062575832638 },
        { -0.48995620076727442, 1.3614783753698549, 0.10208492479454145 },
        { -0.00245012513064022, 0.0060496121419582978, 1.0139163722126245 }
    };

    OIIO_CHECK_EQUAL( idt.size(), 3 );
    OIIO_CHECK_EQUAL( idt[0].size(), 3 );
    for ( size_t i = 0; i < 3; i++ )
        for ( size_t j = 0; j < 3; j++ )
            OIIO_CHECK_EQUAL_THRESH( idt[i][j], true_idt[i][j], 1e-15 );
}

int main( int, char ** )
{
    for ( auto use_eigen: { false, true } )
    {
        rta::core::use_eigen = use_eigen;

        test_metadata_IDT();
    }

    return unit_test_failures;
}
