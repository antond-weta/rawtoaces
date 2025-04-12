// Copyright Contributors to the rawtoaces project.
// SPDX-License-Identifier: Apache-2.0
// https://github.com/AcademySoftwareFoundation/rawtoaces

#pragma once

#include "spectral_data.hpp"

namespace rta
{
namespace core
{

class SpectralSolver
{
public:
    SpectralData illuminant;
    SpectralData camera;
    SpectralData observer;
    SpectralData training_data;

    std::vector<double> white_balance;

    int verbosity = 0;

    void create_blackbody_illuminant( float temp, SpectralData &out ) const;

    void create_daylight_illuminant(
        int temp, const SpectralData &components, SpectralData &out ) const;

    void scale_illuminant( SpectralData &illuminant ) const;

    std::vector<double> calculate_white_balance(
        const SpectralData &illuminant, const SpectralData &observer ) const;

    void prepare_training_data();

    std::vector<std::vector<double>> calculate_training_XYZ() const;
    std::vector<std::vector<double>> calculate_training_RGB() const;
    bool calculate_IDT( std::vector<std::vector<double>> &out ) const;
};

} // namespace core
} // namespace rta
