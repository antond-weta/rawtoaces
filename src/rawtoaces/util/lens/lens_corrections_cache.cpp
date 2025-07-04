// Copyright Contributors to the rawtoaces project.
// SPDX-License-Identifier: Apache-2.0
// https://github.com/AcademySoftwareFoundation/rawtoaces

#include "lens_corrections.hpp"
#include "lens_corrections_cache.hpp"

namespace rta
{
namespace util
{
namespace lens
{

std::ostream &
operator<<( std::ostream &stream, const LensCorrectionsDescriptor &descriptor )
{
    switch ( descriptor.type )
    {
        case LensCorrectionsCacheEntryType::Aberration:
            stream << "chromatic aberration:";
            break;
        case LensCorrectionsCacheEntryType::Distortion:
            stream << "distortion:";
            break;
        case LensCorrectionsCacheEntryType::Vignette:
            stream << "vignette:";
            break;
    }

    stream << std::endl;

    const char *keys[] = {
        "CameraMake", "CameraModel", "LensMake", "LensModel", "FocalLength"
    };

    for ( size_t i = 0; i < sizeof( keys ) / sizeof( keys[0] ); i++ )
    {
        stream << "  " << keys[i] << ": "
               << descriptor.options.get_string( keys[i] ) << std::endl;
    }

    if ( descriptor.type == LensCorrectionsCacheEntryType::Vignette )
    {
        const char *keys[] = { "Aperture", "FocusDistance" };

        for ( size_t i = 0; i < sizeof( keys ) / sizeof( keys[0] ); i++ )
        {
            stream << "  " << keys[i] << ": "
                   << descriptor.options.get_string( keys[i] ) << std::endl;
        }
    }

    //    stream << std::endl;
    return stream;
}

bool LensCorrectionsDescriptor::operator==(
    const CacheEntryDescriptor &other ) const
{
    return false;
}

bool LensCorrectionsDescriptor::operator==(
    const LensCorrectionsDescriptor &other ) const
{
    if ( type != other.type )
        return false;

    const char *keys[] = {
        "CameraMake", "CameraModel", "LensMake", "LensModel", "FocalLength"
    };

    for ( size_t i = 0; i < sizeof( keys ) / sizeof( keys[0] ); i++ )
    {
        auto s1 = options.get_string( keys[i] );
        auto s2 = other.options.get_string( keys[i] );
        if ( s1 != s2 )
            return false;
    }

    if ( type == LensCorrectionsCacheEntryType::Vignette )
    {
        const char *keys[] = { "Aperture", "FocusDistance" };

        for ( size_t i = 0; i < sizeof( keys ) / sizeof( keys[0] ); i++ )
        {
            auto s1 = options.get_string( keys[i] );
            auto s2 = other.options.get_string( keys[i] );
            if ( s1 != s2 )
                return false;
        }
    }

    return true;
}

size_t LensCorrectionsDescriptor::map_index() const
{
    return type;
}

bool LensCorrectionsDescriptor::fetch(
    OIIO::ImageBuf &data, int verbosity ) const
{
    bool result;

    switch ( type )
    {
        case LensCorrectionsCacheEntryType::Aberration:
            result =
                make_aberration_map( data, {}, spec, options, false, nthreads );
            break;
        case LensCorrectionsCacheEntryType::Distortion:
            result = make_distortion_map( data, {}, spec, options, nthreads );
            break;
        case LensCorrectionsCacheEntryType::Vignette:
            result =
                make_vignette_map( data, {}, spec, options, false, nthreads );
            break;
    }

    return result;
}

} // namespace lens
} // namespace util
} // namespace rta

template class rta::cache::
    CacheBase<rta::util::lens::LensCorrectionsDescriptor, OIIO::ImageBuf, 3>;
