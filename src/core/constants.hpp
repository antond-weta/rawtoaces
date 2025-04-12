// Copyright Contributors to the rawtoaces project.
// SPDX-License-Identifier: Apache-2.0
// https://github.com/AcademySoftwareFoundation/rawtoaces

#pragma once

#include <vector>

namespace rta
{
namespace core
{

// clang-format off
static const std::vector<double> ACES_white_XYZ = {
    0.952646074569846, 1.0,    1.00882518435159
};

static const std::vector<std::vector<double>> ACES_to_XYZ = {
    { 0.952552395938186, 0.0,                9.36786316604686e-05 },
    { 0.343966449765075, 0.728166096613485, -0.0721325463785608   },
    { 0.0,               0.0,                1.00882518435159     }
};

static const std::vector<std::vector<double>> XYZ_to_ACES = {
    {  1.0498110175, 0.0000000000, -0.0000974845 },
    { -0.4959030231, 1.3733130458,  0.0982400361 },
    {  0.0000000000, 0.0000000000,  0.9912520182 }
};

static const std::vector<std::vector<double>> XYZ_D65_to_ACES = {
    {  1.0634731317028,    0.00639793641966071, -0.0157891874506841 },
    { -0.492082784686793,  1.36823709310019,     0.0913444629573544 },
    { -0.0028137154424595, 0.00463991165243123,  0.91649468506889   }
};

static const std::vector<std::vector<float>> XYZ_to_ACES_transposed = {
    {  1.0498110175, -0.4959030231, 0.0000000000, 0.0 },
    {  0.0000000000,  1.3733130458, 0.0000000000, 0.0 },
    { -0.0000974845,  0.0982400361, 0.9912520182, 0.0 },
    {  0.0000000000,  0.0000000000, 0.0000000000, 1.0 }
};

} // namespace core
} // namespace rta
