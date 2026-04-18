// SPDX-License-Identifier: Apache-2.0
// Copyright Contributors to the rawtoaces Project.

#include <rawtoaces/image_converter.h>

/// Amend the error message by adding info on the relevant command line parameter.
void update_error_message( std::string &error_message )
{
    static const std::map<std::string, std::string> message_map = {
        { "Missing the camera manufacturer name in the file metadata.",
          " You can provide a custom value using the --custom-camera-make "
          "parameter." },
        { "Missing the camera model name in the file metadata.",
          " You can provide a custom value using the --custom-camera-model "
          "parameter." },
        { "Missing the lens model name in the file metadata.",
          " You can provide a custom value using the --custom-lens-model "
          "parameter." },
        { "Missing the aperture value in the file metadata.",
          " You can provide a custom value using the --custom-aperture "
          "parameter." },
        { "Missing the focal length value in the file metadata.",
          " You can provide a custom value using the --custom-focal-length "
          "parameter." },
        { "Please provide the database path using the RAWTOACES_DATA_PATH "
          "environment variable",
          " or the --data-dir parameter" }
    };

    for ( auto [key, value]: message_map )
    {
        auto found = error_message.find( key );
        if ( found != std::string::npos )
        {
            error_message.insert( found + key.length(), value );
        }
    }
}

bool stack_files(
    const std::vector<std::string> &files,
    rta::util::ImageConverter      &converter )
{
    converter.process_stack( files );
    return true;
}

bool stack_batch(
    const std::vector<std::string> &batch,
    int                             stack_size,
    rta::util::ImageConverter      &converter )
{
    if ( batch.size() % stack_size != 0 )
    {
        std::cerr << "ERROR: the number of files in a batch must be a multiple "
                  << "of the stack size." << std::endl;
        return false;
    }

    size_t stacks = batch.size() / stack_size;

    for ( size_t i = 0; i < stacks; i++ )
    {
        size_t begin_offset = i * stack_size;
        size_t end_offset   = begin_offset + stack_size;

        std::vector<std::string> files_to_stack = std::vector<std::string>(
            batch.begin() + begin_offset, batch.begin() + end_offset );
        if ( !stack_files( files_to_stack, converter ) )
        {
            return false;
        }
    }

    return true;
}

int main( int argc, const char *argv[] )
{
#ifndef WIN32
    putenv( (char *)"TZ=UTC" );
#else
    _putenv( (char *)"TZ=UTC" );
#endif

    rta::util::ImageConverter converter;

    OIIO::ArgParse arg_parser;
    arg_parser.arg( "filename" ).action( OIIO::ArgParse::append() ).hidden();

    arg_parser.separator( "Stacking options:" );

    arg_parser.arg( "--stack-size" )
        .help( "Number of images to stack." )
        .metavar( "VAL" )
        .defaultval( 5 )
        .action( OIIO::ArgParse::store<int>() );

    arg_parser.separator( "Options inherited from rawtoaces:" );

    converter.init_parser( arg_parser );

    arg_parser.parse_args( argc, argv );

    if ( !converter.parse_parameters( arg_parser ) )
    {
        std::string error_message = converter.last_error_message;
        update_error_message( error_message );
        std::cerr << "Error: " << error_message << std::endl;
        return 1;
    }

    int stack_size = arg_parser["stack-size"].get<int>();

    auto files = arg_parser["filename"].as_vec<std::string>();
    if ( files.empty() || ( files.size() == 1 && files[0] == "" ) )
    {
        arg_parser.print_help();
        return 1;
    }

    // Gather all the raw images from arg list
    std::vector<std::vector<std::string>> batches =
        rta::util::collect_image_files( files ); // LCOV_EXCL_LINE

    // Process raw files
    //    bool empty  = true;
    bool result = true;

    for ( auto const &batch: batches )
    {
        if ( !batch.empty() )
            stack_batch( batch, stack_size, converter );
    }

    //    if ( empty )
    //        arg_parser.print_help();

    return result ? 0 : 1;
}
