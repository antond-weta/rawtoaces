// Copyright Contributors to the rawtoaces project.
// SPDX-License-Identifier: Apache-2.0
// https://github.com/AcademySoftwareFoundation/rawtoaces

#include <OpenImageIO/unittest.h>

#include <rawtoaces/util/rawtoaces_util.hpp>

void test_util()
{
    rta::util::ImageConverter converter;
}

int main( int, char ** )
{
    test_util();
    return unit_test_failures;
}
