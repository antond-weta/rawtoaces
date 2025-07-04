// Copyright Contributors to the rawtoaces project.
// SPDX-License-Identifier: Apache-2.0
// https://github.com/AcademySoftwareFoundation/rawtoaces

#include <OpenImageIO/unittest.h>

#include "../src/rawtoaces/util/data_adapters/exiftool.hpp"

#include "test_data.hpp"

using namespace std;
using namespace rta;

#if defined( WIN32 ) || defined( WIN64 )
const std::string separator = ";";

int setenv( const char *name, const char *value, int overwrite )
{
    return _putenv_s( name, value );
}

#else
const std::string separator = ":";
#endif

const std::string data_path = "../../../data/" + separator + "../../data/";

const int verbosity = 2;

void check_exiftool()
{
    OIIO::ImageSpec spec;
    
    const std::string path = "../../tests/materials/BatteryPark.NEF";
    
    const std::vector<std::string> keys = {
        "cameraMake",
        "cameraModel",
        "lensModel",
        "aperture",
        "focalLength",
    };
    
    rta::util::exiftool::fetch_metadata(spec, path, keys);
    
    OIIO_CHECK_EQUAL( spec.get_string_attribute("cameraMake"), "NIKON CORPORATION");
    OIIO_CHECK_EQUAL( spec.get_string_attribute("cameraModel"), "NIKON D200");
    OIIO_CHECK_EQUAL( spec.get_string_attribute("lensModel"), "AF Zoom-Nikkor 28-70mm f/3.5-4.5D");
    OIIO_CHECK_EQUAL( spec.get_float_attribute("aperture"), 8);
    OIIO_CHECK_EQUAL( spec.get_float_attribute("focalLength"), 28);
}



int main( int, char ** )
{
    check_exiftool();

    return unit_test_failures;
}
