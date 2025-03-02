// Copyright Contributors to the rawtoaces project.
// SPDX-License-Identifier: Apache-2.0
// https://github.com/AcademySoftwareFoundation/rawtoaces

#pragma once

#include <rawtoaces/spectral_data.hpp>

namespace rta
{
namespace core
{

class MatrixSolver
{
private:
    const std::vector<std::vector<double>> &_camera_points;
    const std::vector<std::vector<double>> &_xyz_points;
    const double ( &_white_xyz )[3];
    const double ( &_white_rgb )[3];
    const double ( &_out_to_xyz_mat )[3][3];

    double _x[6];

    std::vector<std::vector<double>> _lab_points;

public:
    size_t max_steps       = 100;
    double error_threshold = 1e-17;

    MatrixSolver(
        const std::vector<std::vector<double>> &camera_points,
        const std::vector<std::vector<double>> &xyz_points,
        const double ( &white_xyz )[3],
        const double ( &white_rgb )[3],
        const double ( &out_to_xyz_mat )[3][3] )
        : _camera_points( camera_points )
        , _xyz_points( xyz_points )
        , _white_xyz( white_xyz )
        , _white_rgb( white_rgb )
        , _out_to_xyz_mat( out_to_xyz_mat ){};

    std::vector<std::vector<double>> solve();

private:
    void objective( size_t i, int j, double ( &result )[3] );
};

} // namespace core
} // namespace rta
