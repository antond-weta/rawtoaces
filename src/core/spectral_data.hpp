// Copyright Contributors to the rawtoaces project.
// SPDX-License-Identifier: Apache-2.0
// https://github.com/AcademySoftwareFoundation/rawtoaces

#pragma once

#include <vector>
#include <map>

#include "spectrum.hpp"

namespace rta
{
namespace core
{

struct SpectralData
{
    // Header data
    std::string manufacturer;
    std::string model;
    std::string illuminant;
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
};

} // namespace core
} // namespace rta
