// Copyright Contributors to the rawtoaces project.
// SPDX-License-Identifier: Apache-2.0
// https://github.com/AcademySoftwareFoundation/rawtoaces

#include "lens_corrections_lensfun.hpp"

#include <OpenImageIO/imagebufalgo_util.h>
#include <OpenImageIO/imagebufalgo.h>
#include <lensfun/lensfun.h>

namespace rta
{
namespace util
{
namespace lens
{
namespace lensfun
{

static lfDatabase *Database;

const lfModifier *modifierFromSpec(
    const OIIO::ImageSpec      &spec,
    const OIIO::ParamValueList &options,
    bool                        inverse,
    bool                        enable_vignetting,
    bool                        enable_distortion,
    bool                        enable_aberration )
{
    std::string camera_make    = options.get_string( "CameraMake" );
    std::string camera_model   = options.get_string( "CameraModel" );
    std::string lens_make      = options.get_string( "LensMake" );
    std::string lens_model     = options.get_string( "LensModel" );
    float       focal_length   = options.get_float( "FocalLength" );
    float       aperture       = options.get_float( "Aperture" );
    float       focus_distance = options.get_float( "FocusDistance" );

    if ( Database == nullptr )
    {
        char *path = getenv( "RAWTOACES_LENSFUNDB_PATH" );
        if ( path == nullptr )
        {
            std::cerr
                << "Lensfun DB not found, please provide the path to "
                << "the database directory via the RAWTOACES_LENSFUNDB_PATH "
                << "environment variable." << std::endl;
            exit( -1 );
        }

        Database = new lfDatabase();
#if ( LF_VERSION <= ( 0 << 24 ) | ( 3 << 16 ) | ( 4 << 8 ) )
        Database->LoadDirectory( path );
#else
        Database->Load( path );
#endif // LF_VERSION
    }

    const lfCamera **cameras =
        Database->FindCamerasExt( camera_make.c_str(), camera_model.c_str() );
    auto cam = cameras[0];

    const lfLens **lenses =
        Database->FindLenses( cam, NULL, lens_model.c_str() );
    auto lens = lenses[0];

#if ( LF_VERSION <= ( 0 << 24 ) | ( 3 << 16 ) | ( 4 << 8 ) )

    int flags = 0;

    if ( enable_vignetting )
        flags += LF_MODIFY_VIGNETTING;

    if ( enable_distortion )
        flags += LF_MODIFY_DISTORTION;

    if ( enable_aberration )
        flags += LF_MODIFY_TCA;

    lfModifier *mod = new lfModifier(
        lens, cam->CropFactor, spec.full_width, spec.full_height );
    mod->Initialize(
        lens,
        LF_PF_F32,
        focal_length,
        aperture,
        focus_distance,
        1.0,
        LF_UNKNOWN,
        flags,
        inverse );

    return mod;
#else
    lfModifier *mod = new lfModifier(
        lens,
        focal_length,
        cam->CropFactor,
        spec.full_width,
        spec.full_height,
        LF_PF_F32,
        inverse );

    if ( enable_vignetting )
        mod->EnableVignettingCorrection( aperture, focus_distance );

    if ( enable_distortion )
        mod->EnableDistortionCorrection();

    if ( enable_aberration )
        mod->EnableTCACorrection();

    return mod;
#endif // LF_VERSION
}

bool make_vignette_map(
    OIIO::ImageBuf             &dst,
    OIIO::ROI                   roi,
    const OIIO::ImageSpec      &spec,
    const OIIO::ParamValueList &options,
    bool                        inverse,
    int                         nthreads )
{
    const lfModifier *modifier =
        modifierFromSpec( spec, options, inverse, true, false, false );

    if ( !roi.defined() )
        roi = dst.roi();

    OIIO::ImageBufAlgo::parallel_image( roi, nthreads, [&]( OIIO::ROI roi ) {
        OIIO::ImageBuf::Iterator<float> iterator( dst, roi );
        size_t                          count = roi.width();

        std::vector<float> buff( count );

        for ( int y = roi.ybegin; y < roi.yend; y++ )
        {
            for ( size_t x = 0; x < count; x++ )
            {
                buff[x] = 1.0;
            }

            modifier->ApplyColorModification(
                buff.data(),
                roi.xbegin - spec.full_x,
                y - spec.full_y,
                roi.width(),
                1,
                LF_CR_1( INTENSITY ),
                0 );

            size_t i = 0;
            for ( int x = roi.xbegin; x < roi.xend; x++ )
            {
                for ( int c = roi.chbegin; c < roi.chend; c++ )
                {
                    iterator[c] = buff[i];
                }
                iterator++;
                i++;
            }
        }
    } );

    delete modifier;

    return true;
}

bool make_distortion_map(
    OIIO::ImageBuf             &dst,
    OIIO::ROI                   roi,
    const OIIO::ImageSpec      &spec,
    const OIIO::ParamValueList &options,
    int                         nthreads )
{
    const lfModifier *modifier =
        modifierFromSpec( spec, options, false, false, true, false );

    if ( !roi.defined() )
        roi = dst.roi();

    const float offsets[2] = { (float)spec.full_x, (float)spec.full_y };

    const float scales[2] = { 1.0f / spec.full_width, 1.0f / spec.full_height };

    OIIO::ImageBufAlgo::parallel_image( roi, nthreads, [&]( OIIO::ROI roi ) {
        OIIO::ImageBuf::Iterator<float> iterator( dst, roi );
        size_t                          count = roi.width() * 2 * 3;

        std::vector<float> buff( count );

        for ( int y = roi.ybegin; y < roi.yend; y++ )
        {
            modifier->ApplyGeometryDistortion(
                (float)roi.xbegin - spec.full_x + 0.5,
                (float)y - spec.full_y + 0.5,
                roi.width(),
                1,
                buff.data() );

            size_t i = 0;
            for ( int x = roi.xbegin; x < roi.xend; x++ )
            {
                for ( int j = 0; j < 2; j++ )
                {
                    float val = buff[i++];
                    val       = ( val + offsets[j] ) * scales[j];
                    iterator[roi.chbegin + j] = val;
                }
                iterator++;
            }
        }
    } );

    return true;
}

bool make_aberration_map(
    OIIO::ImageBuf             &dst,
    OIIO::ROI                   roi,
    const OIIO::ImageSpec      &spec,
    const OIIO::ParamValueList &options,
    bool                        inverse,
    int                         nthreads )
{
    const lfModifier *modifier =
        modifierFromSpec( spec, options, false, false, false, true );

    if ( !roi.defined() )
        roi = dst.roi();

    const float offsets[6] = { (float)spec.full_x, (float)spec.full_y,
                               (float)spec.full_x, (float)spec.full_y,
                               (float)spec.full_x, (float)spec.full_y };

    const float scales[6] = {
        1.0f / spec.full_width, 1.0f / spec.full_height,
        1.0f / spec.full_width, 1.0f / spec.full_height,
        1.0f / spec.full_width, 1.0f / spec.full_height,
    };

    OIIO::ImageBufAlgo::parallel_image( roi, nthreads, [&]( OIIO::ROI roi ) {
        OIIO::ImageBuf::Iterator<float> iterator( dst, roi );
        size_t                          count = roi.width() * 2 * 3;

        std::vector<float> buff( count );

        for ( int y = roi.ybegin; y < roi.yend; y++ )
        {

            modifier->ApplySubpixelDistortion(
                (float)roi.xbegin - spec.full_x + 0.5,
                (float)y - spec.full_y + 0.5,
                roi.width(),
                1,
                buff.data() );

            int i = 0;
            for ( int x = roi.xbegin; x < roi.xend; x++ )
            {
                for ( int j = 0; j < 6; j++ )
                {
                    float val = buff[i++];
                    val       = ( val + offsets[j] ) * scales[j];
                    iterator[roi.chbegin + j] = val;
                }
                iterator++;
            }
        }
    } );

    return true;
}

} // namespace lensfun
} // namespace lens
} // namespace util
} // namespace rta
