// Copyright Contributors to the rawtoaces project.
// SPDX-License-Identifier: Apache-2.0
// https://github.com/AcademySoftwareFoundation/rawtoaces

#pragma once

#include <vector>
#include <map>
#include <iostream>
#include <fstream>

#include <nlohmann/json.hpp>

#include "spectrum.hpp"

namespace rta
{
namespace core
{

inline void
parse_string( nlohmann::json &j, std::string &dst, const std::string &key )
{
    auto &v = j[key];
    if ( v.is_null() )
        dst = "";
    else
        dst = v;
}

struct SpectralData
{
    // Header data
    std::string manufacturer;
    std::string model;
    std::string catalog_number;
    std::string description;
    std::string document_creator;
    std::string unique_identifier;
    std::string measurement_equipment;
    std::string laboratory;
    std::string creation_date;
    std::string comments;
    std::string license;

    // Spectral data
    std::string units;
    std::string reflection_geometry;
    std::string transmission_geometry;
    std::string bandwidth_FWHM;
    std::string bandwidth_corrected;

    std::map<std::string, std::vector<std::string>> index;
    std::map<std::string, std::vector<Spectrum>>    data;

    SpectralData() {}

    SpectralData( const std::string &path, bool reshape = false )
    {
        Spectrum::Shape shape;

        std::ifstream  i( path );
        nlohmann::json file_data = nlohmann::json::parse( i );

        nlohmann::json &h = file_data["header"];
        parse_string( h, manufacturer, "manufacturer" );
        parse_string( h, model, "model" );
        parse_string( h, catalog_number, "catalog_number" );
        parse_string( h, description, "description" );
        parse_string( h, document_creator, "document_creator" );
        parse_string( h, unique_identifier, "unique_identifier" );
        parse_string( h, measurement_equipment, "measurement_equipment" );
        parse_string( h, laboratory, "laboratory" );
        parse_string( h, creation_date, "document_creation_date" );
        parse_string( h, comments, "comments" );
        parse_string( h, license, "license" );

        nlohmann::json &d = file_data["spectral_data"];
        parse_string( d, units, "units" );
        parse_string( d, reflection_geometry, "reflection_geometry" );
        parse_string( d, transmission_geometry, "transmission_geometry" );
        parse_string( d, bandwidth_FWHM, "bandwidth_FWHM" );
        parse_string( d, bandwidth_corrected, "bandwidth_corrected" );

        float prev_wavelength = -1;

        std::vector<float>       wavelengths;
        std::vector<std::string> keys;

        nlohmann::json spectral_index = file_data["spectral_data"]["index"];

        for ( auto &[k, v]: spectral_index.items() )
        {
            std::cerr << k << std::endl;
            std::string key = k;
            auto        i   = index.emplace( key, 0 );

            for ( auto kk: v )
            {
                i.first->second.push_back( kk );
            }
        }

        nlohmann::json spectral_data = file_data["spectral_data"]["data"];

        for ( auto &[k, v]: spectral_data.items() )
        {
            std::string key   = k;
            size_t      count = index[key].size();

            auto &vec = data[key];
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
                        return;
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

        for ( auto &[k, v]: data )
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
    }
};

} // namespace core
} // namespace rta
