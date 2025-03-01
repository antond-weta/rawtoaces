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

    void scale_illuminant( SpectralData &illuminant ) const;
    std::vector<double>
    calculate_white_balance( SpectralData &illuminant ) const;

    std::vector<std::vector<double>> calculate_CAT(
        const std::vector<double> &src_white_XYZ,
        const std::vector<double> &dst_white_XYZ ) const;
};

} // namespace core
} // namespace rta
