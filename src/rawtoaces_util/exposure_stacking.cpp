// SPDX-License-Identifier: Apache-2.0
// Copyright Contributors to the rawtoaces Project.

#include "exposure_stacking.h"
#include <OpenImageIO/imagebufalgo.h>
#include <OpenImageIO/imagebufalgo_util.h>

bool ExposureStacking::precheck(
    const OIIO::ImageSpec &spec, float scale, bool is_reference )
{
    if ( ( spec.width > 0 && _stacked_image.spec().width == 0 ) ||
         ( spec.height > 0 && _stacked_image.spec().height == 0 ) )
    {
        OIIO::ImageSpec new_spec = spec;
        new_spec.format          = OIIO::TypeDesc::FLOAT;

        _stacked_image = OIIO::ImageBuf( new_spec );
        _clipping_map  = OIIO::ImageBuf( new_spec );

        _accumulator.resize( spec.width * spec.height );
    }

    if ( is_reference )
    {
        reference_scale = scale;
    }

    return true;
}

bool ExposureStacking::process(
    const OIIO::ImageBuf &image,
    float                 black_level,
    float                 low_threshold,
    float                 high_threshold,
    float                 scale )
{
    auto spec           = _stacked_image.spec();
    auto image_iterator = OIIO::ImageBuf::ConstIterator<uint16_t, uint16_t>(
        image, image.spec().x, image.spec().y );

    if ( reference_scale > 0.0 )
        scale /= reference_scale;

    OIIO::ROI roi_full = image.roi();
    OIIO::ImageBufAlgo::parallel_image( roi_full, 0, [&]( OIIO::ROI roi ) {
        OIIO::ImageBuf::ConstIterator<uint16_t, uint16_t> image_iterator(
            image, roi );

        for ( int y = roi.ybegin; y < roi.yend; y++ )
        {
            auto accum_iterator =
                _accumulator.begin() +
                ( y - roi_full.ybegin ) * roi_full.width() +
                ( roi.xbegin - roi_full.xbegin );

            for ( int x = roi.xbegin; x < roi.xend; x++ )
            {
                float value = image_iterator[0];
                auto &accum = *accum_iterator;

                size_t new_quality = 0;
                if ( value < low_threshold )
                    new_quality = 2;
                else if ( value < high_threshold )
                    new_quality = 1;

                if ( accum.quality < new_quality )
                {
                    // This sample is of better quality than the ones in the
                    // accumulator. Reset the accumulator state.
                    accum.quality  = new_quality;
                    accum.count    = 0;
                    accum.value    = 0;
                    accum.weight   = 0;
                    accum.clipping = 0;
                }

                if ( accum.quality == new_quality )
                {
                    if ( new_quality == 2 )
                    {
                        // For all samples below the lower threshold, add them
                        // up together, with the weight equal the inverse of the
                        // scale. The actual values don't need to be scaled.
                        // This gives the values already scaled according to
                        // their exposure, e.g. the samples having lower
                        // exposure (and higher noise) would have lower weight.
                        accum.value += value - black_level;
                        accum.weight += 1.0 / scale;
                        accum.clipping += 0.0;
                        accum.count++;
                    }
                    else if ( new_quality == 1 )
                    {
                        // A sample approaching clipping. Set the weight to
                        // the position of the value between the low and high
                        // threshold. Also scale to account for exposure.
                        float weight = ( value - low_threshold ) /
                                       ( high_threshold - low_threshold );
                        accum.value +=
                            ( value - black_level ) * scale * weight;
                        accum.weight += weight;
                        accum.clipping += weight;
                        accum.count++;
                    }
                    else
                    {
                        // Hard clipping. The sample with the highest exposure
                        // wins.
                        accum.value =
                            std::max( accum.value, high_threshold * scale );
                        accum.weight   = 1.0;
                        accum.clipping = 1.0;
                        accum.count    = 1;
                    }
                }

                image_iterator++;
                accum_iterator++;
            }
        }
    } );

    return true;
}

bool ExposureStacking::finalise()
{
    OIIO::ROI roi_full = _stacked_image.roi();
    OIIO::ImageBufAlgo::parallel_image(
        roi_full, 0, [&]( OIIO::ROI roi ) {
            OIIO::ImageBuf::Iterator<float> image_iterator(
                _stacked_image, roi );
            OIIO::ImageBuf::Iterator<float> clipping_iterator(
                _clipping_map, roi );

            for ( int y = roi.ybegin; y < roi.yend; y++ )
            {
                auto accum_iterator =
                    _accumulator.begin() +
                    ( y - roi_full.ybegin ) * roi_full.width() +
                    ( roi.xbegin - roi_full.xbegin );

                for ( int x = roi.xbegin; x < roi.xend; x++ )
                {
                    auto &accum = *accum_iterator;

                    if ( accum.count == 0 )
                    {
                        // An empty pixel with no samples.
                        // This should never happen.
                        // Set to zero and mark as saturated?
                        accum.count    = 1;
                        accum.weight   = 1.0;
                        accum.clipping = 1.0;
                    }

                    float value = ( accum.value / accum.weight ) / white_level;
                    image_iterator[0]    = value;
                    clipping_iterator[0] = accum.clipping / accum.count;

                    image_iterator++;
                    clipping_iterator++;
                    accum_iterator++;
                }
            }
        } );

    return true;
}

const OIIO::ImageBuf &ExposureStacking::stacked_image()
{
    return _stacked_image;
}

const OIIO::ImageBuf &ExposureStacking::clipping_map()
{
    return _clipping_map;
}
