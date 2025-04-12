// Copyright Contributors to the rawtoaces project.
// SPDX-License-Identifier: Apache-2.0
// https://github.com/AcademySoftwareFoundation/rawtoaces

#pragma once

namespace rta
{
namespace core
{

extern bool use_eigen;
extern bool use_ceres;

// clang-format off

static const std::vector<std::vector<double>> CAT02 = {
    {  0.7328, 0.4296, -0.1624 },
    { -0.7036, 1.6975,  0.0061 },
    {  0.0030, 0.0136,  0.9834 }
};

static const std::vector<std::vector<double>> CAT02_inv = {
    {  1.096124, -0.278869,  0.182745 },
    {  0.454369,  0.473533,  0.072098 },
    { -0.009628, -0.005698,  1.015326 }
};

// clang-format on

template <typename T>
std::vector<std::vector<T>> transposed( std::vector<std::vector<T>> mat );

template <typename T>
std::vector<std::vector<T>>
product( std::vector<std::vector<T>> mat1, std::vector<std::vector<T>> mat2 );

template <typename T>
std::vector<T> product( std::vector<std::vector<T>> mat, std::vector<T> vec );

template <typename T>
bool inverse(
    const std::vector<std::vector<T>> &in, std::vector<std::vector<T>> &out );

/// Solve a system of linear equations using Gauss elimination.
bool solve_linear(
    std::vector<std::vector<double>> &a, std::vector<double> &b );

template <typename Objective>
bool solve_nonlinear(
    double ( &x )[],
    const Objective &objective,
    double           error_threshold = 1e-17,
    size_t           max_steps       = 100,
    size_t           verbosity       = 0 );

bool solve_IDT_LAB(
    const std::vector<std::vector<double>> &src_points,
    const std::vector<std::vector<double>> &xyz_points,
    const std::vector<double>              &white_rgb,
    const std::vector<std::vector<double>> &out_to_xyz_mat,
    std::vector<std::vector<double>>       &out,
    int                                     verbosity );

std::vector<std::vector<double>> calculate_CAT(
    const std::vector<double> &src_white_XYZ,
    const std::vector<double> &dst_white_XYZ );

} // namespace core
} // namespace rta
