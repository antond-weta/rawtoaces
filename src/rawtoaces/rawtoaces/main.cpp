// Copyright Contributors to the rawtoaces project.
// SPDX-License-Identifier: Apache-2.0
// https://github.com/AcademySoftwareFoundation/rawtoaces

#include <rawtoaces/util/rawtoaces_util.hpp>

#include <filesystem>

bool check_and_add(
    const std::filesystem::path &path, std::vector<std::string> &batch )
{
    if ( std::filesystem::is_regular_file( path ) ||
         std::filesystem::is_symlink( path ) )
    {
        auto filename = path.filename();

        if ( filename == ".DS_Store" )
            return false;

        auto extension = path.extension();

        if ( extension == ".exr" || extension == ".EXR" )
            return false;
        if ( extension == ".jpg" || extension == ".JPG" )
            return false;
        if ( extension == ".jpeg" || extension == ".JPEG" )
            return false;

        std::string str = path.string();
        batch.push_back( str );
    }
    else
    {
        std::cerr << "Not a regular file: " << path << std::endl;
    }
    return true;
}

int main( int argc, const char *argv[] )
{
    OIIO::ArgParse argParse;
    argParse.arg( "filename" ).action( OIIO::ArgParse::append() ).hidden();

    rta::util::ImageConverter converter;
    converter.init_parser( argParse );

    if ( argParse.parse_args( argc, argv ) < 0 )
    {
        return 1;
    }

    if ( !converter.parse_params( argParse ) )
    {
        return 1;
    }

    // Create a separate batch for each input directory.
    // Reserve the first batch for the individual input files.
    std::vector<std::vector<std::string>> batches( 1 );

    auto files = argParse["filename"].as_vec<std::string>();

    if ( files.empty() || ( files.size() == 1 && files[0] == "" ) )
    {
        argParse.print_help();
        return 1;
    }

    for ( auto filename: files )
    {
        if ( !std::filesystem::exists( filename ) )
        {
            std::cerr << "File or directory not found: " << filename
                      << std::endl;
            return 1;
        }

        auto canonical_filename = std::filesystem::canonical( filename );

        if ( std::filesystem::is_directory( filename ) )
        {
            std::vector<std::string> &curr_batch = batches.emplace_back();
            auto it = std::filesystem::directory_iterator( filename );

            for ( auto filename2: it )
            {
                if ( !check_and_add( filename2, curr_batch ) )
                    continue;
            }
        }
        else
        {
            if ( !check_and_add( filename, batches[0] ) )
                continue;
        }
    }

    bool result = true;
    for ( auto const &batch: batches )
    {
        for ( auto const &input_filename: batch )
        {
            result = converter.process_image( input_filename );
            if ( !result )
                break;
        }
        if ( !result )
            break;
    }

    return result ? 0 : 1;
}
