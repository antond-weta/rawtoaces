// Copyright Contributors to the rawtoaces project.
// SPDX-License-Identifier: Apache-2.0
// https://github.com/AcademySoftwareFoundation/rawtoaces

#pragma once

#include <OpenImageIO/imagebuf.h>
#include <OpenImageIO/imagebufalgo.h>

namespace rta
{
namespace util
{
namespace lens
{
namespace lensfun
{

bool make_vignette_map(
    OIIO::ImageBuf             &dst,
    OIIO::ROI                   roi,
    const OIIO::ImageSpec      &spec,
    const OIIO::ParamValueList &options,
    bool                        inverse  = false,
    int                         nthreads = 0 );

bool make_distortion_map(
    OIIO::ImageBuf             &dst,
    OIIO::ROI                   roi,
    const OIIO::ImageSpec      &spec,
    const OIIO::ParamValueList &options,
    int                         nthreads = 0 );

bool make_aberration_map(
    OIIO::ImageBuf             &dst,
    OIIO::ROI                   roi,
    const OIIO::ImageSpec      &spec,
    const OIIO::ParamValueList &options,
    bool                        inverse  = false,
    int                         nthreads = 0 );

} // namespace lensfun
} // namespace lens
} // namespace util
} // namespace rta
