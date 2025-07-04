// Copyright Contributors to the rawtoaces project.
// SPDX-License-Identifier: Apache-2.0
// https://github.com/AcademySoftwareFoundation/rawtoaces

#pragma once

#include <OpenImageIO/imagebuf.h>
#include <OpenImageIO/imagebufalgo.h>
#include <list>

#include "../cache_base.hpp"

namespace rta
{
namespace util
{
namespace lens
{

enum LensCorrectionsCacheEntryType
{
    Aberration,
    Distortion,
    Vignette
};

class LensCorrectionsDescriptor
    : public rta::cache::CacheEntryDescriptor<OIIO::ImageBuf>
{
public:
    LensCorrectionsCacheEntryType type;
    OIIO::ImageSpec               spec;
    OIIO::ParamValueList          options;
    int                           nthreads;

    bool   operator==( const CacheEntryDescriptor &other ) const override;
    bool   operator==( const LensCorrectionsDescriptor &other ) const;
    size_t map_index() const override;
    bool   fetch( OIIO::ImageBuf &data, int verbosity = 0 ) const override;

    friend std::ostream &operator<<(
        std::ostream &stream, const LensCorrectionsDescriptor &descriptor );

    std::tuple<const OIIO::ImageSpec &, OIIO::InitializePixels &&>
    construct_entry() const
    {
        auto xxx = std::forward_as_tuple( spec, OIIO::InitializePixels::No );
        return xxx;
    };
};

class LensCorrectionsCache
    : public rta::cache::CacheBase<LensCorrectionsDescriptor, OIIO::ImageBuf, 3>
{
public:
    LensCorrectionsCache() : CacheBase()
    {
        capacity = 5;
        name     = "Lens Corrections cache";
    }
};

} // namespace lens
} // namespace util
} // namespace rta
