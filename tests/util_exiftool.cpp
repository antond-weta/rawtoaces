// Copyright Contributors to the rawtoaces project.
// SPDX-License-Identifier: Apache-2.0
// https://github.com/AcademySoftwareFoundation/rawtoaces

#include <OpenImageIO/unittest.h>

#include "../src/rawtoaces/util/data_adapters/exiftool.hpp"

void check_exiftool()
{

#ifdef RTA_XCODE
    const std::string path = "../../../tests/materials/BatteryPark.NEF";
#else
    const std::string path = "../../tests/materials/BatteryPark.NEF";
#endif

    const std::vector<std::string> keys = {
        "cameraMake", "cameraModel", "lensModel", "aperture", "focalLength",
    };

    OIIO::ImageSpec spec;
    rta::util::exiftool::fetch_metadata( spec, path, keys );

    OIIO_CHECK_EQUAL(
        spec.get_string_attribute( "cameraMake" ), "NIKON CORPORATION" );
    OIIO_CHECK_EQUAL(
        spec.get_string_attribute( "cameraModel" ), "NIKON D200" );
    OIIO_CHECK_EQUAL(
        spec.get_string_attribute( "lensModel" ),
        "AF Zoom-Nikkor 28-70mm f/3.5-4.5D" );
    OIIO_CHECK_EQUAL( spec.get_float_attribute( "aperture" ), 8 );
    OIIO_CHECK_EQUAL( spec.get_float_attribute( "focalLength" ), 28 );
}

int main( int, char ** )
{
    check_exiftool();

    return unit_test_failures;
}
