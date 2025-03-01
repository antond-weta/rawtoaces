// Copyright Contributors to the rawtoaces project.
// SPDX-License-Identifier: Apache-2.0
// https://github.com/AcademySoftwareFoundation/rawtoaces

#pragma once

#include <rawtoaces/spectral_data.hpp>

namespace rta
{
namespace core
{

class SpectralSolver
{
public:
    SpectralData current_camera;
    SpectralData current_illuminant;
    SpectralData colour_matching;
    SpectralData training_data;

    std::vector<double> white_balance;

    void scale_illuminant( SpectralData &illuminant ) const;

    std::vector<double> calculate_white_balance(
        const SpectralData &illuminant, const SpectralData &observer ) const;

    std::vector<std::vector<double>> calculate_CAT(
        const std::vector<double> &src_white_XYZ,
        const std::vector<double> &dst_white_XYZ ) const;

    void prepare_training_data( const std::string &path );

    std::vector<std::vector<double>> calculate_training_XYZ() const;
    std::vector<std::vector<double>> calculate_training_RGB() const;
};

} // namespace core
} // namespace rta
