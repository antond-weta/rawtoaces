// Copyright Contributors to the rawtoaces project.
// SPDX-License-Identifier: Apache-2.0
// https://github.com/AcademySoftwareFoundation/rawtoaces

#pragma once

#include "metadata.hpp"

namespace rta
{
namespace core
{

class MetadataSolver
{
public:
    MetadataSolver( const Metadata &metadata );
    std::vector<std::vector<double>> calculate_IDT() const;

private:
    Metadata _metadata;
};

} // namespace core
} // namespace rta
