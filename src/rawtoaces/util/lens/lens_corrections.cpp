// Copyright Contributors to the rawtoaces project.
// SPDX-License-Identifier: Apache-2.0
// https://github.com/AcademySoftwareFoundation/rawtoaces

#include "lens_corrections.hpp"

#if defined( RTA_ENABLE_LENSFUN )
#    include "lens_corrections_lensfun.hpp"
namespace ns = rta::util::lens::lensfun;
#endif

namespace rta
{
namespace util
{
namespace lens
{

bool make_vignette_map(
    OIIO::ImageBuf             &dst,
    OIIO::ROI                   roi,
    const OIIO::ImageSpec      &spec,
    const OIIO::ParamValueList &options,
    bool                        inverse,
    int                         nthreads )
{
    return ns::make_vignette_map( dst, roi, spec, options, inverse, nthreads );
}

OIIO::ImageBuf make_vignette_map(
    const OIIO::ImageSpec      &spec,
    const OIIO::ParamValueList &options,
    bool                        inverse,
    int                         nthreads )
{
    OIIO::ImageBuf result( spec, OIIO::InitializePixels::No );
    bool ok = make_vignette_map( result, {}, spec, options, nthreads );
    if ( !ok )
    {
        result.errorfmt( "Failed to calculate vignette map" );
    }
    return result;
}

bool make_distortion_map(
    OIIO::ImageBuf             &dst,
    OIIO::ROI                   roi,
    const OIIO::ImageSpec      &spec,
    const OIIO::ParamValueList &options,
    int                         nthreads )
{
    return ns::make_distortion_map( dst, roi, spec, options, nthreads );
}

OIIO::ImageBuf make_distortion_map(
    const OIIO::ImageSpec      &spec,
    const OIIO::ParamValueList &options,
    int                         nthreads )
{
    OIIO::ImageBuf result( spec, OIIO::InitializePixels::No );
    OIIO::ROI      roi;
    bool ok = make_distortion_map( result, roi, spec, options, nthreads );
    if ( !ok )
    {
        result.errorfmt( "Failed to calculate distortion map" );
    }
    return result;
}

bool make_aberration_map(
    OIIO::ImageBuf             &dst,
    OIIO::ROI                   roi,
    const OIIO::ImageSpec      &spec,
    const OIIO::ParamValueList &options,
    bool                        inverse,
    int                         nthreads )
{
    return ns::make_aberration_map(
        dst, roi, spec, options, inverse, nthreads );
}

OIIO::ImageBuf make_aberration_map(
    const OIIO::ImageSpec      &spec,
    const OIIO::ParamValueList &options,
    bool                        inverse,
    int                         nthreads )
{
    OIIO::ImageBuf result( spec, OIIO::InitializePixels::No );
    bool           ok =
        make_aberration_map( result, {}, spec, options, inverse, nthreads );
    if ( !ok )
    {
        result.errorfmt( "Failed to calculate chromatic aberration map" );
    }
    return result;
}

} // namespace lens
} // namespace util
} // namespace rta
