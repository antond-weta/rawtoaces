// SPDX-License-Identifier: Apache-2.0
// Copyright Contributors to the rawtoaces Project.

#pragma once

#include <vector>

#include <OpenImageIO/imagebuf.h>

class ExposureStacking
{
public:
    bool
    precheck( const OIIO::ImageSpec &spec, float scale, bool is_reference );
    bool process(
        const OIIO::ImageBuf &image,
        float                 black_level    = 0.0f,
        float                 low_threshold  = 0.8f,
        float                 high_threshold = 0.9f,
        float                 scale          = 1.0f,
        int                   nthreads       = 0 );
    bool finalise( int nthreads = 0 );

    const OIIO::ImageBuf &stacked_image();
    const OIIO::ImageBuf &clipping_map();

    float white_level     = 1.0f;
    float reference_scale = 0.0f;

private:
    OIIO::ImageBuf _stacked_image;
    OIIO::ImageBuf _clipping_map;

    struct Accumulator
    {
        /// The quality of the accumulated samples.
        /// 0 - Either an empty pixel (if count == 0), or a fully saturated pixel (count > 0).
        /// 1 - A pixel approaching clipping.
        /// 2 - A good linear pixel.
        size_t quality;
        /// The number of accumulated samples.
        size_t count;
        /// Accumulated value.
        float value;
        /// Total weight of the accumulated value.
        float weight;
        /// Total clipping.
        float clipping;
    };

    std::vector<Accumulator> _accumulator;
};
