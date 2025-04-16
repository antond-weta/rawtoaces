// Copyright Contributors to the rawtoaces project.
// SPDX-License-Identifier: Apache-2.0
// https://github.com/AcademySoftwareFoundation/rawtoaces

#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <filesystem>

#ifdef _MSC_VER
#    define strncasecmp _strnicmp
#    define strcasecmp _stricmp
#endif

#include "spectral_loader.hpp"

#include <nlohmann/json.hpp>

#include "../core/spectral_data.hpp"
#include "../core/spectral_solver.hpp"

namespace rta
{
namespace util
{

// Adapted from define.h pathsFinder()
const std::vector<std::string> &data_paths()
{
    static std::vector<std::string> data_paths;
    static bool                     first_time = 1;

    if ( first_time )
    {
#if defined( WIN32 ) || defined( WIN64 )
        const std::string separator    = ";";
        const std::string default_path = ".";
#else
        char              separator    = ':';
        const std::string default_path = "/usr/local/share/rawtoaces/data";
#endif

        std::string path;

        {
            const char *env = getenv( "AMPAS_DATA_PATH" );
            if ( env != nullptr )
            {
                if ( !path.empty() )
                    path += separator;
                path += env;
            }
        }

        {
            const char *env = getenv( "RAWTOACES_DATA_PATH" );
            if ( env != nullptr )
            {
                if ( !path.empty() )
                    path += separator;
                path += env;
            }
        }

        if ( path.empty() )
        {
            path = default_path;
        }

        size_t pos = 0;

        while ( pos < path.size() )
        {
            size_t end = path.find( separator, pos );

            if ( end == std::string::npos )
                end = path.size();

            std::string pathItem = path.substr( pos, end - pos );

            if ( find( data_paths.begin(), data_paths.end(), pathItem ) ==
                 data_paths.end() )
                data_paths.push_back( pathItem );

            pos = end + 1;
        }
    }
    return data_paths;
};

// Adapted from define.h isValidCT()
// Function to check if a string is a valid input
// to represent color temperature(s) (e.g., D60, 3200K)
bool is_valid_CT( std::string str )
{
    int    i      = 0;
    size_t length = str.length();

    if ( length == 0 )
        return false;

    bool first_is_d = std::tolower( str[0] ) == 'd';
    bool last_id_k  = std::tolower( str[length - 1] ) == 'k';

    // other light sources
    if ( !first_is_d && !last_id_k )
    {
        while ( str[i] )
        {
            if ( !isalnum( str[i] ) && str[i] != '-' )
                return false;
            i++;
        }
    }
    // daylight ("D" + Numeric values)
    else if ( first_is_d && !last_id_k )
    {
        i = 1;
        while ( str[i] )
        {
            if ( !isdigit( str[i] ) )
                return false;
            i++;
        }
    }
    // daylight (Numeric values + "K")
    else if ( !first_is_d && last_id_k )
    {
        while ( str[i] && i < str.length() - 1 )
        {
            if ( !isdigit( str[i] ) )
                return false;
            i++;
        }
    }
    // It is rarely for "D" and "K" appear
    // together as the first letter and
    // the last letter of the string
    else if ( first_is_d && last_id_k )
        return false;

    return true;
};

std::string SpectralLoader::find_file( const std::string &filename )
{
    if ( verbosity > 1 )
    {
        std::cerr << "Looking for a file: " << filename << std::endl;
    }

    for ( const auto &path: data_paths() )
    {
        std::string full_path = path + "/" + filename;
        if ( std::filesystem::exists( full_path ) )
        {
            if ( verbosity > 1 )
            {
                std::cerr << "Found file: " << full_path << std::endl;
            }
            return full_path;
        }
    }

    if ( verbosity > 0 )
    {
        std::cerr << "File " << filename << "not found." << std::endl;
    }
    return "";
}

std::vector<std::string>
SpectralLoader::collect_data_files( const std::string &type )
{
    std::vector<std::string> result;

    for ( const auto &path: data_paths() )
    {
        if ( std::filesystem::is_directory( path ) )
        {
            auto type_path = path + "/" + type;
            if ( std::filesystem::exists( type_path ) )
            {
                auto it = std::filesystem::directory_iterator( type_path );
                for ( auto filename2: it )
                {
                    auto p = filename2.path();
                    if ( filename2.path().extension() == ".json" )
                    {
                        result.push_back( filename2.path().string() );
                    }
                }
            }
        }
    }
    return result;
}

inline void
parse_string( nlohmann::json &j, std::string &dst, const std::string &key )
{
    auto &v = j[key];
    if ( v.is_null() )
        dst = "";
    else
        dst = v;
}

core::SpectralData
SpectralLoader::load_spectral_file( std::string path, bool reshape )
{
    core::SpectralData    result;
    core::Spectrum::Shape shape;

    std::ifstream  i( path );
    nlohmann::json file_data = nlohmann::json::parse( i );

    nlohmann::json &h = file_data["header"];
    parse_string( h, result.manufacturer, "manufacturer" );
    parse_string( h, result.model, "model" );
    parse_string( h, result.illuminant, "illuminant" );
    parse_string( h, result.catalog_number, "catalog_number" );
    parse_string( h, result.description, "description" );
    parse_string( h, result.document_creator, "document_creator" );
    parse_string( h, result.unique_identifier, "unique_identifier" );
    parse_string( h, result.measurement_equipment, "measurement_equipment" );
    parse_string( h, result.laboratory, "laboratory" );
    parse_string( h, result.creation_date, "document_creation_date" );
    parse_string( h, result.comments, "comments" );
    parse_string( h, result.license, "license" );

    nlohmann::json &d = file_data["spectral_data"];
    parse_string( d, result.units, "units" );
    parse_string( d, result.reflection_geometry, "reflection_geometry" );
    parse_string( d, result.transmission_geometry, "transmission_geometry" );
    parse_string( d, result.bandwidth_FWHM, "bandwidth_FWHM" );
    parse_string( d, result.bandwidth_corrected, "bandwidth_corrected" );

    float prev_wavelength = -1;

    std::vector<float>       wavelengths;
    std::vector<std::string> keys;

    nlohmann::json spectral_index = file_data["spectral_data"]["index"];

    for ( auto &[k, v]: spectral_index.items() )
    {
        std::string key = k;
        auto        i   = result.index.emplace( key, 0 );

        for ( auto kk: v )
        {
            i.first->second.push_back( kk );
        }
    }

    nlohmann::json spectral_data = file_data["spectral_data"]["data"];

    for ( auto &[k, v]: spectral_data.items() )
    {
        std::string key   = k;
        size_t      count = result.index[key].size();

        auto &vec = result.data[key];
        for ( size_t c = 0; c < count; c++ )
        {
            vec.emplace_back( 0, shape );
        }

        for ( auto &[sub_key, value]: spectral_data[key].items() )
        {
            float this_wavelength = std::stof( sub_key );

            if ( prev_wavelength != -1 )
            {
                float new_step = this_wavelength - prev_wavelength;

                if ( shape.step != 0 && new_step != shape.step )
                {
                    return result;
                }

                shape.step = new_step;
            }
            else
            {
                shape.first = this_wavelength;
            }

            prev_wavelength = this_wavelength;

            for ( int j = 0; j < count; j++ )
            {
                vec[j].values.push_back( value[j] );
            }
        }
    }

    shape.last = prev_wavelength;

    for ( auto &[k, v]: result.data )
    {
        for ( auto &vv: v )
        {
            vv.shape = shape;
            if ( reshape )
            {
                vv.reshape();
            }
        }
    }

    return result;
}

std::vector<core::SpectralData> SpectralLoader::load_all_cameras()
{
    std::vector<core::SpectralData> result;

    auto paths = collect_data_files( "camera" );
    for ( auto &path: paths )
    {
        result.push_back( load_spectral_file( path, false ) );
    }

    return result;
}

bool SpectralLoader::find_camera(
    const std::string &make, const std::string &model, core::SpectralData &out )
{
    auto paths = collect_data_files( "camera" );

    for ( auto &path: paths )
    {
        out = load_spectral_file( path, false );

        if ( strcasecmp( make.c_str(), out.manufacturer.c_str() ) == 0 &&
             strcasecmp( model.c_str(), out.model.c_str() ) == 0 )
        {
            return true;
        }
    }

    return false;
}

std::vector<core::SpectralData> SpectralLoader::load_all_illuminants()
{
    std::vector<core::SpectralData> result;

    auto paths = collect_data_files( "illuminant" );
    for ( auto &path: paths )
    {
        auto illuminant = load_spectral_file( path, false );
        if ( illuminant.data["main"].size() == 1 )
        {
            result.push_back( illuminant );
        }
        else
        {
            if ( verbosity > 0 )
                std::cerr
                    << "The file " << path
                    << " contains more than one spectral component, skipping."
                    << std::endl;
        }
    }

    return result;
}

bool SpectralLoader::load_illuminant(
    const std::string &type, core::SpectralData &out )
{
    // First check if the requested type is one of the built-in illuminants.
    if ( is_valid_CT( type ) )
    {
        out.illuminant = type;

        if ( std::tolower( type[0] ) == 'd' )
        {
            auto path = find_file( "misc/CIE_illum_Dxx_comp.json" );
            if ( path.length() == 0 )
            {
                std::cerr << "Failed to find the illuminant D components "
                          << "data file \"misc/CIE_illum_Dxx_comp.json\". "
                          << "Please check the RAWTOACES_DATA_PATH environment "
                          << "variable." << std::endl;
                return false;
            }
            auto components = load_spectral_file( path, true );

            float cct = std::stoi( type.substr( 1 ) );

            core::SpectralSolver solver;
            solver.verbosity = verbosity;
            solver.create_daylight_illuminant( cct, components, out );
        }
        else if ( std::tolower( type.back() ) == 'k' )
        {
            float cct = std::stoi( type.substr( 0, type.length() - 2 ) );

            core::SpectralSolver solver;
            solver.verbosity = verbosity;
            solver.create_blackbody_illuminant( cct, out );
        }
        return true;
    }

    // Now try to find in the supplied illuminant files.
    auto paths = collect_data_files( "illuminant" );
    for ( auto &path: paths )
    {
        auto temp = load_spectral_file( path, true );
        if ( strcasecmp( temp.model.c_str(), type.c_str() ) == 0 )
        {
            out = temp;
            return true;
        }
    }

    return false;
}

core::SpectralData SpectralLoader::find_best_illuminant(
    const core::SpectralData  &camera,
    const std::vector<double> &wb,
    int                        kelvin_step )
{
    assert( wb.size() == 3 );

    if ( verbosity > 0 )
        std::cerr << "Finding best illuminant for camera ("
                  << camera.manufacturer << ", " << camera.model
                  << ") and white balancing weights (" << wb[0] << ", " << wb[1]
                  << ", " << wb[2] << ")." << std::endl;

    std::vector<core::SpectralData> illuminants = load_all_illuminants();

    core::SpectralSolver solver;
    solver.camera = camera;

    // Blackbody illuminants
    for ( int i = 1500; i < 4000; i += kelvin_step )
    {
        core::SpectralData illuminant;
        solver.create_blackbody_illuminant( i, illuminant );
        illuminants.push_back( illuminant );
    }

    auto components_path = find_file( "misc/CIE_illum_Dxx_comp.json" );
    auto components      = load_spectral_file( components_path, true );

    // Daylight illuminants
    for ( int i = 4000; i <= 25000; i += kelvin_step )
    {
        core::SpectralData illuminant;
        solver.create_daylight_illuminant( i / 100, components, illuminant );
        illuminants.push_back( illuminant );
    }

    double              best_error      = std::numeric_limits<double>::max();
    core::SpectralData *best_illuminant = nullptr;
    for ( auto &illuminant: illuminants )
    {
        solver.scale_illuminant( illuminant );
        auto curr_wb = solver.calculate_white_balance( illuminant, camera );
        assert( curr_wb.size() == 3 );

        double error = 0;
        for ( size_t i = 0; i < 3; i++ )
        {
            error += pow( curr_wb[i] / wb[i] - 1.0, 2 );
        }

        if ( verbosity > 1 )
            std::cerr << "Trying illuminant: " << illuminant.illuminant
                      << ". Error = " << error << ".";

        if ( error < best_error )
        {
            best_error      = error;
            best_illuminant = &illuminant;
            if ( verbosity > 1 )
                std::cerr << " <- new best.";
        }

        if ( verbosity > 1 )
            std::cerr << std::endl;
    }

    if ( verbosity > 0 )
        std::cerr << "Found best illuminant: " << best_illuminant->illuminant
                  << "." << std::endl;

    return *best_illuminant;
}

} // namespace util
} // namespace rta
