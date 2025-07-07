// Copyright Contributors to the rawtoaces project.
// SPDX-License-Identifier: Apache-2.0
// https://github.com/AcademySoftwareFoundation/rawtoaces

#include "rawtoaces_util.hpp"

#include <rawtoaces/core/constants.hpp>
#include <rawtoaces/core/metadata.hpp>

#include "transform_cache.hpp"
#include "data_adapters/spectral_loader.hpp"

#include <OpenImageIO/imageio.h>
#include <OpenImageIO/imagebuf.h>
#include <OpenImageIO/imagebufalgo.h>

#include <filesystem>

namespace rta
{
namespace util
{

cache::TransformCache transform_cache;

const char *HelpString =
    "Rawtoaces converts raw image files from a digital camera to "
    "the Academy Colour Encoding System (ACES) compliant images.\n"
    "The process consists of two parts:\n"
    "- the colour values get converted from the camera native colour "
    "space to the ACES AP0 (see \"SMPTE ST 2065-1\"), and \n"
    "- the image file gets converted from the camera native raw "
    "file format to the ACES Image Container file format "
    "(see \"SMPTE ST 2065-4\").\n"
    "\n"
    "Rawtoaces supports the following white-balancing modes:\n"
    "- \"metadata\" uses the white-balancing coefficients from the raw "
    "image file, provided by the camera.\n"
    "- \"illuminant\" performs white balancing to the illuminant, "
    "provided in the \"--illuminant\" parameter. The list of the "
    "supported illuminants can be seen using the "
    "\"--list-illuminants\" parameter. This mode requires spectral "
    "sensitivity data for the camera model the image comes from. "
    "The list of cameras such data is available for, can be "
    "seen using the \"--list-cameras\" parameter. In addition to the named "
    "illuminants, which are stored under ${RAWTOACES_DATA_PATH}/illuminant, "
    "blackbody illuminants of a given colour temperature can me used (use 'K' "
    "suffix, i.e. '3200K'), as well as daylight illuminants (use the 'D' "
    "prefix, i.e. 'D65').\n"
    "- \"box\" performs white-balancing to make the given region of "
    "the image appear neutral gray. The box position (origin and size) "
    "can be specified using the \"--wb-box\" parameter. In case no such "
    "parameter provided, the whole image is used for white-balancing.\n"
    "- \"custom\" uses the custom white balancing coefficients "
    "provided using the -\"custom-wb\" parameter.\n"
    "\n"
    "Rawtoaces supports the following methods of color matrix "
    "computation:\n"
    "- \"spectral\" uses the camera sensor's spectral sensitivity data "
    "to compute the optimal matrix. This mode requires spectral "
    "sensitivity data for the camera model the image comes from. "
    "The list of cameras such data is available for, can be "
    "seen using the \"--list-cameras\" parameter.\n"
    "- \"metadata\" uses the matrix (matrices) contained in the raw "
    "image file metadata. This mode works best with the images using "
    "the DNG format, as the DNG standard mandates the presense of "
    "such matrices.\n"
    "- \"Adobe\" uses the Adobe coefficients provided by LibRaw. \n"
    "- \"custom\" uses a user-provided color conversion matrix. "
    "A matrix can be specified using the \"--custom-mat\" parameter.\n"
    "\n"
    "The paths rawtoaces uses to search for the spectral sensitivity "
    "data can be specified in the RAWTOACES_DATA_PATH environment "
    "variable.\n";

const char *UsageString =
    "\n"
    "    rawtoaces --wb-method METHOD --mat-method METHOD [PARAMS] "
    "path/to/dir/or/file ...\n"
    "Examples: \n"
    "    rawtoaces --wb-method metadata --mat-method metadata raw_file.dng\n"
    "    rawtoaces --wb-method illuminant --illuminant 3200K --mat-method "
    "spectral raw_file.cr3\n";

template <typename T, typename F1, typename F2>
bool check_param(
    const std::string    &mode_name,
    const std::string    &mode_value,
    const std::string    &param_name,
    const std::vector<T> &param_value,
    size_t                correct_size,
    const std::string    &default_value_message,
    bool                  is_correct_mode,
    F1                    on_success,
    F2                    on_failure )
{
    if ( is_correct_mode )
    {
        if ( param_value.size() == correct_size )
        {
            on_success();
            return true;
        }
        else
        {
            if ( ( param_value.size() == 0 ) ||
                 ( ( param_value.size() == 1 ) && ( param_value[0] == 0 ) ) )
            {
                std::cerr << "Warning: " << mode_name << " was set to \""
                          << mode_value << "\", but no \"--" << param_name
                          << "\" parameter provided. " << default_value_message
                          << std::endl;

                on_failure();
                return false;
            }

            std::cerr << "Warning: The parameter \"" << param_name
                      << "\" must have " << correct_size << " values. "
                      << default_value_message << std::endl;

            on_failure();
            return false;
        }
    }
    else
    {
        if ( ( param_value.size() > 1 ) ||
             ( ( param_value.size() == 1 ) && ( param_value[0] != 0 ) ) )
        {
            std::cerr << "Warning: the \"--" << param_name
                      << "\" parameter provided, but the " << mode_name
                      << " is different from \"" << mode_value << "\". "
                      << default_value_message << std::endl;

            on_failure();
            return false;
        }
        else
        {
            return true;
        }
    }
}

void ImageConverter::init_parser( OIIO::ArgParse &argParse )
{
    argParse.intro( HelpString );
    argParse.usage( UsageString );
    argParse.print_defaults( true );
    argParse.add_help( true );
    argParse.add_version( RTA_VERSION );

    argParse.arg( "--wb-method" )
        .help(
            "White balance method. Supported options: metadata, illuminant, "
            "box, custom." )
        .metavar( "STR" )
        .defaultval( "metadata" )
        .action( OIIO::ArgParse::store() );

    argParse.arg( "--mat-method" )
        .help(
            "IDT matrix calculation method. Supported options: spectral, "
            "metadata, Adobe, custom." )
        .metavar( "STR" )
        .defaultval( "spectral" )
        .action( OIIO::ArgParse::store() );

    argParse.arg( "--illuminant" )
        .help( "Illuminant for white balancing. (default = D55)" )
        .metavar( "STR" )
        .action( OIIO::ArgParse::store() );

    argParse.arg( "--wb-box" )
        .help(
            "Box to use for white balancing. (default = (0,0,0,0) - full "
            "image)" )
        .nargs( 4 )
        .metavar( "X Y W H" )
        .action( OIIO::ArgParse::store<int>() );

    argParse.arg( "--custom-wb" )
        .help( "Custom white balance multipliers." )
        .nargs( 4 )
        .metavar( "R G B G" )
        .action( OIIO::ArgParse::store<float>() );

    argParse.arg( "--custom-mat" )
        .help( "Custom camera RGB to XYZ matrix." )
        .nargs( 9 )
        .metavar( "Rr Rg Rb Gr Gg Gb Br Bg Bb" )
        .action( OIIO::ArgParse::store<float>() );

    argParse.arg( "--custom-camera-make" )
        .help(
            "Camera manufacturer name to be used for spectral sensitivity "
            "curves lookup. "
            "If present, overrides the value stored in the file metadata." )
        .metavar( "STR" )
        .action( OIIO::ArgParse::store() );

    argParse.arg( "--custom-camera-model" )
        .help(
            "Camera model name to be used for spectral sensitivity "
            "curves lookup. "
            "If present, overrides the value stored in the file metadata." )
        .metavar( "STR" )
        .action( OIIO::ArgParse::store() );

    argParse.arg( "--headroom" )
        .help( "Highlight headroom factor." )
        .metavar( "VAL" )
        .defaultval( 6.0f )
        .action( OIIO::ArgParse::store<float>() );

    argParse.separator( "General options:" );

    argParse.arg( "--overwrite" )
        .help(
            "Allows overwriting existing files. If not set, trying to write "
            "to an existing file will generate an error." )
        .action( OIIO::ArgParse::store_true() );

    argParse.arg( "--output-dir" )
        .help(
            "The directory to write the output files to. "
            "This gets applied to every input directory, so it is better to "
            "be used with a single input directory." )
        .metavar( "STR" )
        .action( OIIO::ArgParse::store() );

    argParse.arg( "--create-dirs" )
        .help( "Create output directories if they don't exist." )
        .action( OIIO::ArgParse::store_true() );

    argParse.arg( "--disable-cache" )
        .help( "Disable the colour space transform cache." )
        .action( OIIO::ArgParse::store_true() );

    argParse.separator( "Raw conversion options:" );

    argParse.arg( "--auto-bright" )
        .help( "Enable automatic exposure adjustment." )
        .action( OIIO::ArgParse::store_true() );

    argParse.arg( "--adjust-maximum-threshold" )
        .help(
            "Automatically lower the linearity threshold provided in the "
            "metadata by this scaling factor." )
        .metavar( "VAL" )
        .defaultval( 0.75f )
        .action( OIIO::ArgParse::store<float>() );

    argParse.arg( "--black-level" )
        .help( "If >= 0, override the black level." )
        .metavar( "VAL" )
        .defaultval( -1 )
        .action( OIIO::ArgParse::store<int>() );

    argParse.arg( "--saturation-level" )
        .help(
            "If not 0, override the level which appears to be saturated "
            "after normalisation." )
        .metavar( "VAL" )
        .defaultval( 0 )
        .action( OIIO::ArgParse::store<int>() );

    argParse.arg( "--chromatic-aberration" )
        .help(
            "Red and blue scale factors for chromatic aberration correction. "
            "The value of 1 means no correction." )
        .metavar( "R B" )
        .nargs( 2 )
        .defaultval( 1.0f )
        .action( OIIO::ArgParse::store<float>() );

    argParse.arg( "--half-size" )
        .help( "If present, decode image at half size resolution." )
        .action( OIIO::ArgParse::store_true() );

    argParse.arg( "--highlight-mode" )
        .help( "0 = clip, 1 = unclip, 2 = blend, 3..9 = rebuild." )
        .metavar( "VAL" )
        .defaultval( 0 )
        .action( OIIO::ArgParse::store<int>() );

    argParse.arg( "--crop-box" )
        .help(
            "Apply custom crop. If not present, the default crop is applied, "
            "which should match the crop of the in-camera JPEG." )
        .nargs( 4 )
        .metavar( "X Y W H" )
        .action( OIIO::ArgParse::store<int>() );

    argParse.arg( "--crop-mode" )
        .help(
            "Cropping mode. Supported options: 'none' (write out the full "
            "sensor area), 'soft' (write out full image, mark the crop as the "
            "display window), 'hard' (write out only the crop area)." )
        .metavar( "STR" )
        .defaultval( "soft" )
        .action( OIIO::ArgParse::store() );

    argParse.arg( "--flip" )
        .help(
            "If not 0, override the orientation specified in the metadata. "
            "1..8 correspond to EXIF orientation codes "
            "(3 = 180 deg, 6 = 90 deg CCW, 8 = 90 deg CW.)" )
        .metavar( "VAL" )
        .defaultval( 0 )
        .action( OIIO::ArgParse::store<int>() );

    argParse.arg( "--denoise_threshold" )
        .help( "Wavelet denoising threshold" )
        .metavar( "VAL" )
        .defaultval( 0 )
        .action( OIIO::ArgParse::store<int>() );

    argParse.separator( "Benchmarking and debugging:" );

    argParse.arg( "--list-cameras" )
        .help( "Shows the list of cameras supported in spectral mode." )
        .action( OIIO::ArgParse::store_true() );

    argParse.arg( "--list-illuminants" )
        .help( "Shows the list of illuminants supported in spectral mode." )
        .action( OIIO::ArgParse::store_true() );

    argParse.arg( "--verbose" )
        .help(
            "(-v) Print progress messages. "
            "Repeated -v will increase verbosity." )
        .action( [&]( OIIO::cspan<const char *> argv ) { verbosity++; } );

    argParse.arg( "-v" ).hidden().action(
        [&]( OIIO::cspan<const char *> argv ) { verbosity++; } );
}

bool ImageConverter::parse_params( const OIIO::ArgParse &argParse )
{
    if ( argParse["list-cameras"].get<int>() )
    {
        std::cout
            << "Spectral sensitivity data are available for the following "
            << "cameras:" << std::endl;

        std::vector<std::string> camera_list;
        util::SpectralLoader     spectral_loader;

        auto cameras = spectral_loader.load_all_cameras();
        for ( auto &camera: cameras )
        {
            camera_list.push_back( camera.manufacturer + " " + camera.model );
        }

        std::sort( camera_list.begin(), camera_list.end() );

        for ( auto s: camera_list )
            std::cout << "- " << s << std::endl;
        std::cout << std::endl;
        exit( 0 );
    }

    if ( argParse["list-illuminants"].get<int>() )
    {
        std::cout << "The following illuminants are supported:" << std::endl;
        std::cout << "- The standard illuminant series D (e.g., D60, D6025)"
                  << std::endl;
        std::cout << "- Black-body radiation (e.g., 3200K)" << std::endl;

        std::vector<std::string> illuminant_list;
        util::SpectralLoader     spectral_loader;

        auto illuminants = spectral_loader.load_all_illuminants();
        for ( auto &illuminant: illuminants )
        {
            illuminant_list.push_back( illuminant.illuminant );
        }

        std::sort( illuminant_list.begin(), illuminant_list.end() );

        for ( auto s: illuminant_list )
            std::cout << "- " << s << std::endl;
        std::cout << std::endl;
        exit( 0 );
    }

    std::string wb_method = argParse["wb-method"].get();

    if ( wb_method == "metadata" )
    {
        wbMethod = WBMethod::Metadata;
    }
    else if ( wb_method == "illuminant" )
    {
        wbMethod = WBMethod::Illuminant;
    }
    else if ( wb_method == "box" )
    {
        wbMethod = WBMethod::Box;
    }
    else if ( wb_method == "custom" )
    {
        wbMethod = WBMethod::Custom;
    }
    else
    {
        std::cerr << std::endl
                  << "Unsupported white balancing method: \"" << wb_method
                  << "\"." << std::endl;

        return false;
    }

    std::string mat_method = argParse["mat-method"].get();

    if ( mat_method == "spectral" )
    {
        matrixMethod = MatrixMethod::Spectral;
    }
    else if ( mat_method == "metadata" )
    {
        matrixMethod = MatrixMethod::Metadata;
    }
    else if ( mat_method == "Adobe" )
    {
        matrixMethod = MatrixMethod::Adobe;
    }
    else if ( mat_method == "custom" )
    {
        matrixMethod = MatrixMethod::Custom;
    }
    else
    {
        std::cerr << std::endl
                  << "Unsupported matrix method: \"" << mat_method << "\"."
                  << std::endl;

        return false;
    }

    headroom = argParse["headroom"].get<float>();

    illuminant = argParse["illuminant"].get();

    if ( wbMethod == WBMethod::Illuminant )
    {
        if ( illuminant.empty() )
        {
            std::cerr << "Warning: the white balancing method was set to "
                      << "\"illuminant\", but no \"--illuminant\" parameter "
                      << "provided. " << illuminant << " will be used."
                      << std::endl;
        }
    }
    else
    {
        if ( !illuminant.empty() )
        {
            std::cerr << "Warning: the \"--illuminant\" parameter provided "
                      << "but the white balancing mode different from "
                      << "\"illuminant\" "
                      << "requested. The custom illuminant will be ignored."
                      << std::endl;
        }
    }

    auto box = argParse["wb-box"].as_vec<int>();
    check_param(
        "white balancing mode",
        "box",
        "wb-box",
        box,
        4,
        "The box will be ignored.",
        wbMethod == WBMethod::Box,
        [&]() {
            for ( int i = 0; i < 4; i++ )
                wbBox[i] = box[i];
        },
        [&]() {
            for ( int i = 0; i < 4; i++ )
                wbBox[i] = 0;
        } );

    auto custom_wb = argParse["custom-wb"].as_vec<float>();
    check_param(
        "white balancing mode",
        "custom",
        "custom-wb",
        custom_wb,
        4,
        "The scalers will be ignored. The default values of (1, 1, 1, 1) will be used",
        wbMethod == WBMethod::Custom,
        [&]() {
            for ( int i = 0; i < 4; i++ )
                customWB[i] = custom_wb[i];
        },
        [&]() {
            for ( int i = 0; i < 4; i++ )
                customWB[i] = 1.0;
        } );

    auto custom_mat = argParse["custom-mat"].as_vec<float>();
    check_param(
        "matrix mode",
        "custom",
        "custom-mat",
        custom_mat,
        9,
        "Identity matrix will be used",
        matrixMethod == MatrixMethod::Custom,
        [&]() {
            for ( int i = 0; i < 3; i++ )
                for ( int j = 0; j < 3; j++ )
                    customMatrix[i][j] = custom_mat[i * 3 + j];
        },
        [&]() {
            for ( int i = 0; i < 3; i++ )
                for ( int j = 0; j < 3; j++ )
                    customMatrix[i][j] = i == j ? 1.0 : 0.0;
        } );

    auto crop = argParse["crop-box"].as_vec<int>();
    if ( crop.size() == 4 )
    {
        for ( size_t i = 0; i < 4; i++ )
            cropbox[i] = crop[i];
    }

    std::string cropping_mode = argParse["crop-mode"].get();

    if ( cropping_mode == "off" )
    {
        crop_mode = CropMode::Off;
    }
    else if ( cropping_mode == "soft" )
    {
        crop_mode = CropMode::Soft;
    }
    else if ( cropping_mode == "hard" )
    {
        crop_mode = CropMode::Hard;
    }
    else
    {
        std::cerr << std::endl
                  << "Unsupported cropping mode: \"" << cropping_mode << "\"."
                  << std::endl;

        return false;
    }

    custom_camera_make  = argParse["custom-camera-make"].get();
    custom_camera_model = argParse["custom-camera-model"].get();

    auto_bright              = argParse["auto-bright"].get<int>();
    adjust_maximum_threshold = argParse["adjust-maximum-threshold"].get<int>();
    black_level              = argParse["black-level"].get<int>();
    saturation_level         = argParse["saturation-level"].get<int>();
    half_size                = argParse["half-size"].get<int>();
    highlight_mode           = argParse["highlight-mode"].get<int>();
    flip                     = argParse["flip"].get<int>();

    disable_cache = argParse["disable-cache"].get<int>();
    overwrite     = argParse["overwrite"].get<int>();
    create_dirs   = argParse["create-dirs"].get<int>();
    output_dir    = argParse["output-dir"].get();

    return true;
}

/// Normalise the metadata in the cases where the OIIO attribute name
/// doesn't match the standard OpenEXR and/or ACES Container attribute name.
/// We only check the attribute names which are set by the raw input plugin.
void fix_metadata( OIIO::ImageSpec &spec )
{
    const std::map<std::string, std::string> standard_mapping = {
        { "Make", "cameraMake" }, { "Model", "cameraModel" }
    };

    for ( auto i: standard_mapping )
    {
        auto &src_name = i.first;
        auto &dst_name = i.second;

        auto src_attribute = spec.find_attribute( src_name );
        auto dst_attribute = spec.find_attribute( dst_name );

        if ( dst_attribute == nullptr && src_attribute != nullptr )
        {
            auto type = src_attribute->type();
            if ( type.arraylen == 0 )
            {
                if ( type.basetype == OIIO::TypeDesc::STRING )
                    spec[dst_name] = src_attribute->get_string();
                else if ( type.basetype == OIIO::TypeDesc::FLOAT )
                    spec[dst_name] = src_attribute->get_float();
            }
            spec.erase_attribute( src_name );
        }
    }
}

bool ImageConverter::configure(
    const OIIO::ImageSpec &imageSpec, OIIO::ParamValueList &options )
{
    options["raw:use_camera_wb"] = 0;
    options["raw:use_auto_wb"]   = 0;

    options["raw:auto_bright"]        = (int)auto_bright;
    options["raw:adjust_maximum_thr"] = adjust_maximum_threshold;
    options["raw:user_black"]         = black_level;
    options["raw:user_sat"]           = saturation_level;
    options["raw:half_size"]          = (int)half_size;
    options["raw:user_flip"]          = flip;
    options["raw:HighlightMode"]      = highlight_mode;

    if ( cropbox[2] != 0 && cropbox[3] != 0 )
    {
        options.attribute(
            "raw:cropbox", OIIO::TypeDesc( OIIO::TypeDesc::INT, 4 ), cropbox );
    }

    _is_DNG = imageSpec.extra_attribs.find( "raw:dng:version" )->get_int() > 0;

    switch ( wbMethod )
    {
        case WBMethod::Metadata: {
            float user_mul[4];

            for ( int i = 0; i < 4; i++ )
            {
                user_mul[i] = imageSpec.find_attribute( "raw:cam_mul" )
                                  ->get_float_indexed( i );
            }

            options.attribute(
                "raw:user_mul",
                OIIO::TypeDesc( OIIO::TypeDesc::FLOAT, 4 ),
                user_mul );

            //            if ( _read_raw )
            {
                _WB_mults.resize( 4 );
                for ( size_t i = 0; i < 4; i++ )
                    _WB_mults[i] = user_mul[i];
            }
            break;
        }
        case WBMethod::Illuminant: {
            std::string lower_illuminant = OIIO::Strutil::lower( illuminant );
            break;
        }
        case WBMethod::Box:

            if ( wbBox[2] == 0 || wbBox[3] == 0 )
            {
                // Empty box, use whole image.
                options["raw:use_auto_wb"] = 1;
            }
            else
            {
                int32_t box[4];
                for ( int i = 0; i < 4; i++ )
                {
                    box[i] = wbBox[i];
                }
                options.attribute(
                    "raw:greybox",
                    OIIO::TypeDesc( OIIO::TypeDesc::INT, 4 ),
                    box );
            }
            break;

        case WBMethod::Custom:
            options.attribute(
                "raw:user_mul",
                OIIO::TypeDesc( OIIO::TypeDesc::FLOAT, 4 ),
                customWB );

            //            if ( _read_raw )
            {
                _WB_mults.resize( 4 );
                for ( size_t i = 0; i < 4; i++ )
                    _WB_mults[i] = customWB[i];
            }
            break;

        default:
            std::cerr
                << "ERROR: This white balancing method has not been configured "
                << "properly." << std::endl;
            return false;
    }

    switch ( matrixMethod )
    {
        case MatrixMethod::Spectral:
            options["raw:ColorSpace"]        = "raw";
            options["raw:use_camera_matrix"] = 0;
            break;
        case MatrixMethod::Metadata:
            options["raw:ColorSpace"]        = "XYZ";
            options["raw:use_camera_matrix"] = _is_DNG ? 1 : 3;
            break;
        case MatrixMethod::Adobe:
            options["raw:ColorSpace"]        = "XYZ";
            options["raw:use_camera_matrix"] = 1;
            break;
        case MatrixMethod::Custom:
            options["raw:ColorSpace"]        = "raw";
            options["raw:use_camera_matrix"] = 0;

            _IDT_matrix.resize( 3 );
            for ( int i = 0; i < 3; i++ )
            {
                _IDT_matrix[i].resize( 3 );
                for ( int j = 0; j < 3; j++ )
                {
                    _IDT_matrix[i][j] = customMatrix[i][j];
                }
            }
            break;
        default:
            std::cerr
                << "ERROR: This matrix method has not been configured properly."
                << std::endl;
            return false;
    }

    bool spectral_white_balance = wbMethod == WBMethod::Illuminant;
    bool spectral_matrix        = matrixMethod == MatrixMethod::Spectral;

    if ( spectral_white_balance || spectral_matrix )
    {
        if ( !prepareIDT_spectral(
                 imageSpec, spectral_white_balance, spectral_matrix ) )
        {
            std::cerr << "ERROR: the colour space transfor has not been "
                      << "configured properly (spectral mode)." << std::endl;
            return false;
        }

        if ( spectral_white_balance )
        {
            float user_mul[4];

            for ( int i = 0; i < _WB_mults.size(); i++ )
            {
                user_mul[i] = _WB_mults[i];
            }
            if ( _WB_mults.size() == 3 )
                user_mul[3] = _WB_mults[1];

            options.attribute(
                "raw:user_mul",
                OIIO::TypeDesc( OIIO::TypeDesc::FLOAT, 4 ),
                user_mul );
        }
    }

    if ( matrixMethod == MatrixMethod::Metadata )
    {
        if ( _is_DNG )
        {
            options["raw:use_camera_matrix"] = 1;
            options["raw:use_camera_wb"]     = 1;

            if ( !prepareIDT_DNG( imageSpec ) )
            {
                std::cerr << "ERROR: the colour space transfor has not been "
                          << "configured properly (metadata mode)."
                          << std::endl;
                return false;
            }
        }
        else
        {
            std::cerr
                << "ERROR: This matrix method \"metadata\" only works with "
                << "DNG files." << std::endl;
            return false;
        }
    }
    else if ( matrixMethod == MatrixMethod::Adobe )
    {
        if ( !prepareIDT_nonDNG( imageSpec ) )
        {
            std::cerr << "ERROR: the colour space transfor has not been "
                      << "configured properly (Adobe mode)." << std::endl;
            return false;
        }
    }

    return true;
}

bool ImageConverter::configure(
    const std::string &input_filename, OIIO::ParamValueList &options )
{
    _read_raw = options.get_string( "raw:Demosaic" ) == "none";

    OIIO::ImageSpec imageSpec;

    options["raw:ColorSpace"]    = "XYZ";
    options["raw:use_camera_wb"] = 0;
    options["raw:use_auto_wb"]   = 0;

    OIIO::ImageSpec temp_spec;
    temp_spec.extra_attribs = options;

    auto imageInput = OIIO::ImageInput::create( "raw", false, &temp_spec );
    bool result     = imageInput->open( input_filename, imageSpec, temp_spec );
    if ( !result )
    {
        return false;
    }

    fix_metadata( imageSpec );
    return configure( imageSpec, options );
}

bool ImageConverter::applyMatrix(
    const std::vector<std::vector<double>> &matrix,
    OIIO::ImageBuf                         &dst,
    const OIIO::ImageBuf                   &src,
    OIIO::ROI                               roi )
{
    float M[4][4];

    size_t n = matrix.size();

    if ( n )
    {
        size_t m = matrix[0].size();

        for ( size_t i = 0; i < n; i++ )
        {
            for ( size_t j = 0; j < m; j++ )
            {
                M[j][i] = matrix[i][j];
            }

            for ( size_t j = m; j < 4; j++ )
                M[j][i] = 0;
        }

        for ( size_t i = n; i < 4; i++ )
        {
            for ( size_t j = 0; j < m; j++ )
                M[j][i] = 0;
            for ( size_t j = m; j < 4; j++ )
                M[j][i] = 1;
        }
    }

    return OIIO::ImageBufAlgo::colormatrixtransform( dst, src, M, false, roi );
}

bool ImageConverter::apply_matrix(
    OIIO::ImageBuf &dst, const OIIO::ImageBuf &src, OIIO::ROI roi )
{
    bool success = true;

    if ( !roi.defined() )
        roi = dst.roi();

    if ( _IDT_matrix.size() )
    {
        success = applyMatrix( _IDT_matrix, dst, src, roi );
        if ( !success )
            return false;
    }

    if ( _CAT_matrix.size() )
    {
        success = applyMatrix( _CAT_matrix, dst, dst, roi );
        if ( !success )
            return false;

        success = applyMatrix( core::XYZ_to_ACES, dst, dst, roi );
        if ( !success )
            return false;
    }

    return success;
}

bool ImageConverter::apply_scale(
    OIIO::ImageBuf &dst, const OIIO::ImageBuf &src, OIIO::ROI roi )
{
    return OIIO::ImageBufAlgo::mul( dst, src, headroom );
}

bool ImageConverter::apply_crop(
    OIIO::ImageBuf &dst, const OIIO::ImageBuf &src, OIIO::ROI roi )
{
    if ( crop_mode == CropMode::Off )
    {
        if ( !OIIO::ImageBufAlgo::copy( dst, src ) )
        {
            return false;
        }
        dst.specmod().full_x      = dst.specmod().x;
        dst.specmod().full_y      = dst.specmod().y;
        dst.specmod().full_width  = dst.specmod().width;
        dst.specmod().full_height = dst.specmod().height;
    }
    else if ( crop_mode == CropMode::Hard )
    {
        // OIIO can not currently crop in place.
        if ( &dst == &src )
        {
            OIIO::ImageBuf temp;
            if ( !OIIO::ImageBufAlgo::copy( temp, src ) )
            {
                return false;
            }

            if ( !OIIO::ImageBufAlgo::crop( dst, temp, temp.roi_full() ) )
            {
                return false;
            }
        }
        else
        {
            if ( !OIIO::ImageBufAlgo::crop( dst, src, src.roi_full() ) )
            {
                return false;
            }
        }
        dst.specmod().x      = 0;
        dst.specmod().y      = 0;
        dst.specmod().full_x = 0;
        dst.specmod().full_y = 0;
    }

    return true;
}

bool ImageConverter::make_output_path(
    std::string &path, const std::string &suffix )
{
    std::filesystem::path temp_path( path );

    temp_path.replace_extension();
    temp_path += suffix + ".exr";

    if ( !output_dir.empty() )
    {
        auto new_directory = std::filesystem::path( output_dir );

        auto filename      = temp_path.filename();
        auto old_directory = temp_path.remove_filename();

        new_directory = old_directory / new_directory;

        if ( !std::filesystem::exists( new_directory ) )
        {
            if ( create_dirs )
            {
                if ( !std::filesystem::create_directory( new_directory ) )
                {
                    std::cerr << "ERROR: Failed to create directory "
                              << new_directory << "." << std::endl;
                    return false;
                }
            }
            else
            {
                std::cerr << "ERROR: The output directory " << new_directory
                          << " does not exist." << std::endl;
                return false;
            }
        }

        temp_path = std::filesystem::absolute( new_directory / filename );
    }

    if ( !overwrite && std::filesystem::exists( temp_path ) )
    {
        std::cerr
            << "ERROR: file " << temp_path << " already exists. Use "
            << "--overwrite to allow overwriting existing files. Skipping "
            << "this file." << std::endl;
        return false;
    }

    path = temp_path.string();
    return true;
}

bool ImageConverter::save(
    const std::string &output_filename, const OIIO::ImageBuf &buf )
{
    const float chromaticities[] = { 0.7347, 0.2653, 0,       1,
                                     0.0001, -0.077, 0.32168, 0.33767 };

    OIIO::ImageSpec imageSpec = buf.spec();
    imageSpec.set_format( OIIO::TypeDesc::HALF );
    imageSpec["acesImageContainerFlag"] = 1;
    imageSpec["compression"]            = "none";
    imageSpec.attribute(
        "chromaticities",
        OIIO::TypeDesc( OIIO::TypeDesc::FLOAT, 8 ),
        chromaticities );

    auto imageOutput = OIIO::ImageOutput::create( "exr" );
    bool result      = imageOutput->open( output_filename, imageSpec );
    if ( result )
    {
        result = buf.write( imageOutput.get() );
    }
    return result;
}

bool ImageConverter::process_image( const std::string &input_filename )
{
    std::string output_filename = input_filename;
    if ( !make_output_path( output_filename ) )
    {
        return ( false );
    }

    OIIO::ParamValueList options;

    if ( !configure( input_filename, options ) )
    {
        std::cerr << "Failed to configure the reader for the file: "
                  << input_filename << std::endl;
        return ( false );
    }

    OIIO::ImageSpec imageSpec;
    imageSpec.extra_attribs = options;

    OIIO::ImageBuf buffer =
        OIIO::ImageBuf( input_filename, 0, 0, nullptr, &imageSpec, nullptr );

    if ( !buffer.read(
             0, 0, 0, buffer.nchannels(), true, OIIO::TypeDesc::FLOAT ) )
    {
        std::cerr << "Failed to read the file: " << input_filename << std::endl;
        return ( false );
    }

    if ( !apply_matrix( buffer, buffer ) )
    {
        std::cerr << "Failed to apply colour space conversion to the file: "
                  << input_filename << std::endl;
        return ( false );
    }

    if ( !apply_scale( buffer, buffer ) )
    {
        std::cerr << "Failed to apply scale to the file: " << input_filename
                  << std::endl;
        return ( false );
    }

    if ( !apply_crop( buffer, buffer ) )
    {
        std::cerr << "Failed to apply crop to the file: " << input_filename
                  << std::endl;
        return ( false );
    }

    if ( !save( output_filename, buffer ) )
    {
        std::cerr << "Failed to save the file: " << output_filename
                  << std::endl;
        return ( false );
    }

    return ( true );
}

bool ImageConverter::fetch_camera_make_and_model(
    const OIIO::ImageSpec &spec,
    std::string           &camera_make,
    std::string           &camera_model )
{
    camera_make = custom_camera_make;
    if ( camera_make.empty() )
        camera_make = spec["cameraMake"];
    if ( camera_make.empty() )
    {
        std::cerr << "Missing the camera manufacturer name in the file "
                  << "metadata. You can provide a camera make using the "
                  << "--custom-camera-make parameter" << std::endl;
        return false;
    }

    camera_model = custom_camera_model;
    if ( camera_model.empty() )
        camera_model = spec["cameraModel"];
    if ( camera_model.empty() )
    {
        std::cerr << "Missing the camera model name in the file metadata. "
                  << "You can provide a camera make using the "
                  << "--custom-camera-make parameter" << std::endl;
        return false;
    }

    return true;
}

bool ImageConverter::prepareIDT_spectral(
    const OIIO::ImageSpec &imageSpec,
    bool                   calc_white_balance,
    bool                   calc_matrix )
{
    std::string lower_illuminant = OIIO::Strutil::lower( illuminant );
    if ( lower_illuminant.empty() )
        lower_illuminant = "na";

    std::string camera_make, camera_model;
    if ( !fetch_camera_make_and_model( imageSpec, camera_make, camera_model ) )
        return false;

    cache::TransformDescriptor descriptor;
    descriptor.camera_make  = camera_make;
    descriptor.camera_model = camera_model;

    if ( lower_illuminant == "na" )
    {
        std::vector<double> wb_multipliers( 4 );

        if ( _WB_mults.size() == 4 )
        {
            for ( int i = 0; i < 3; i++ )
                wb_multipliers[i] = _WB_mults[i];
        }
        else
        {
            auto attr = imageSpec.find_attribute( "raw:pre_mul" );
            for ( int i = 0; i < 4; i++ )
                wb_multipliers[i] = attr->get_float_indexed( i );
        }

        if ( wb_multipliers[3] != 0 )
            wb_multipliers[1] = ( wb_multipliers[1] + wb_multipliers[3] ) / 2.0;
        wb_multipliers.resize( 3 );

        float min_val = std::numeric_limits<float>::max();
        for ( int i = 0; i < 3; i++ )
            if ( min_val > wb_multipliers[i] )
                min_val = wb_multipliers[i];

        if ( min_val > 0 && min_val != 1 )
            for ( int i = 0; i < 3; i++ )
                wb_multipliers[i] /= min_val;

        descriptor.type = cache::TransformEntryType::Illum_from_WB;
        cache::WB &wb   = descriptor.value.emplace<cache::WB>();
        wb.value[0]     = wb_multipliers[0];
        wb.value[1]     = wb_multipliers[1];
        wb.value[2]     = wb_multipliers[2];
    }
    else
    {
        descriptor.type  = cache::TransformEntryType::WB_from_Illum;
        descriptor.value = lower_illuminant;
    }

    transform_cache.disabled  = disable_cache;
    transform_cache.verbosity = verbosity;
    const cache::TransformCacheEntryData &wb_data =
        transform_cache.fetch( descriptor );

    if ( !wb_data.initialised )
    {
        std::cerr << "Failed to calculate the white balancing weights."
                  << std::endl;
        return false;
    }

    if ( calc_white_balance )
    {
        auto p = std::get<std::pair<cache::WB, std::string>>( wb_data.value );
        const cache::WB &wb = p.first;
        _WB_mults.resize( 3 );
        _WB_mults[0] = wb.value[0];
        _WB_mults[1] = wb.value[1];
        _WB_mults[2] = wb.value[2];
    }

    if ( calc_matrix )
    {
        auto &p = std::get<std::pair<cache::WB, std::string>>( wb_data.value );
        std::string illum = p.second;

        cache::TransformDescriptor descriptor;
        descriptor.camera_make  = camera_make;
        descriptor.camera_model = camera_model;
        descriptor.type         = cache::TransformEntryType::Mat_from_Illum;
        descriptor.value        = illum;

        const cache::TransformCacheEntryData &mat_data =
            transform_cache.fetch( descriptor );

        if ( !mat_data.initialised )
        {
            std::cerr << "Failed to calculate the colour transform matrix."
                      << std::endl;
            return false;
        }

        const cache::Matrix33 &mat =
            std::get<cache::Matrix33>( mat_data.value );

        _IDT_matrix.resize( 3 );
        for ( size_t i = 0; i < 3; i++ )
        {
            _IDT_matrix[i].resize( 3 );
            for ( size_t j = 0; j < 3; j++ )
            {
                _IDT_matrix[i][j] = mat.value[i][j];
            }
        }

        _CAT_matrix.resize( 0 );
    }

    return true;
}

bool fetch_matrix(
    cache::TransformDescriptor       &descriptor,
    std::vector<std::vector<double>> &matrix )
{
    const cache::TransformCacheEntryData &mat_data =
        transform_cache.fetch( descriptor );

    if ( !mat_data.initialised )
        return false;

    const cache::Matrix33 &mat = std::get<cache::Matrix33>( mat_data.value );

    matrix.resize( 3 );
    for ( size_t i = 0; i < 3; i++ )
    {
        matrix[i].resize( 3 );
        for ( size_t j = 0; j < 3; j++ )
        {
            matrix[i][j] = mat.value[i][j];
        }
    }

    return true;
}

bool ImageConverter::prepareIDT_DNG( const OIIO::ImageSpec &imageSpec )
{
    std::string camera_make, camera_model;
    if ( !fetch_camera_make_and_model( imageSpec, camera_make, camera_model ) )
        return false;

    cache::TransformDescriptor descriptor;
    descriptor.type          = cache::TransformEntryType::Mat_from_DNG;
    descriptor.camera_make   = camera_make;
    descriptor.camera_model  = camera_model;
    core::Metadata &metadata = descriptor.value.emplace<core::Metadata>();

    metadata.neutralRGB.resize( 3 );
    metadata.xyz2rgbMatrix1.resize( 9 );
    metadata.xyz2rgbMatrix2.resize( 9 );
    metadata.cameraCalibration1.resize( 9 );
    metadata.cameraCalibration2.resize( 9 );

    metadata.baselineExposure =
        imageSpec.get_float_attribute( "raw:dng:baseline_exposure" );
    metadata.calibrationIlluminant1 =
        imageSpec.get_int_attribute( "raw:dng:calibration_illuminant1" );
    metadata.calibrationIlluminant2 =
        imageSpec.get_int_attribute( "raw:dng:calibration_illuminant2" );

    for ( int i = 0; i < 3; i++ )
    {
        metadata.neutralRGB[i] =
            1.0 /
            imageSpec.find_attribute( "raw:cam_mul" )->get_float_indexed( i );

        for ( int j = 0; j < 3; j++ )
        {
            metadata.xyz2rgbMatrix1[i * 3 + j] =
                imageSpec.find_attribute( "raw:dng:color_matrix1" )
                    ->get_float_indexed( i * 3 + j );
            metadata.xyz2rgbMatrix2[i * 3 + j] =
                imageSpec.find_attribute( "raw:dng:color_matrix2" )
                    ->get_float_indexed( i * 3 + j );
            metadata.cameraCalibration1[i * 3 + j] =
                imageSpec.find_attribute( "raw:dng:camera_calibration1" )
                    ->get_float_indexed( i * 4 + j );
            metadata.cameraCalibration2[i * 3 + j] =
                imageSpec.find_attribute( "raw:dng:camera_calibration2" )
                    ->get_float_indexed( i * 4 + j );
        }
    }

    if ( !fetch_matrix( descriptor, _IDT_matrix ) )
    {
        std::cerr << "Failed to calculate the colour transform matrix."
                  << std::endl;
        return false;
    }

    // Do not apply CAT for DNG
    _CAT_matrix.resize( 0 );
    return true;
}

bool ImageConverter::prepareIDT_nonDNG( const OIIO::ImageSpec &imageSpec )
{
    std::string camera_make, camera_model;
    if ( !fetch_camera_make_and_model( imageSpec, camera_make, camera_model ) )
        return false;

    if ( _read_raw )
    {
        auto   mat  = imageSpec.find_attribute( "raw:cam_xyz" );
        size_t size = mat->type().arraylen;

        if ( size == 12 )
        {
            cache::TransformDescriptor descriptor;
            descriptor.type        = cache::TransformEntryType::Mat_from_nonDNG;
            descriptor.camera_make = camera_make;
            descriptor.camera_model = camera_model;

            auto &p =
                descriptor.value
                    .emplace<std::pair<cache::Vector3, cache::Matrix33>>();

            cache::Vector3  &wb   = p.first;
            cache::Matrix33 &mat2 = p.second;

            int idx = 0;
            for ( size_t row = 0; row < 3; row++ )
            {
                for ( size_t col = 0; col < 3; col++ )
                {
                    mat2.value[row][col] = mat->get_float_indexed( idx++ );
                }
            }

            if ( !fetch_matrix( descriptor, _IDT_matrix ) )
            {
                std::cerr << "Failed to calculate the colour transform matrix."
                          << std::endl;
                return false;
            }
        }
    }
    else
    {
        _IDT_matrix.resize( 0 );
    }

    cache::Vector3 src, dst;
    for ( size_t i = 0; i < 3; i++ )
    {
        src.value[i] = rta::core::D65_white_XYZ[i];
        dst.value[i] = rta::core::ACES_white_XYZ[i];
    }

    cache::TransformDescriptor descriptor;
    descriptor.camera_make  = camera_make;
    descriptor.camera_model = camera_model;
    descriptor.type         = cache::TransformEntryType::Mat_from_Illum;
    descriptor.value        = std::pair<cache::WB, cache::WB>( src, dst );

    if ( !fetch_matrix( descriptor, _CAT_matrix ) )
    {
        std::cerr << "Failed to calculate the colour transform matrix."
                  << std::endl;
        return false;
    }

    return true;
}

} // namespace util
} // namespace rta
