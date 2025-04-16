// Copyright Contributors to the rawtoaces project.
// SPDX-License-Identifier: Apache-2.0
// https://github.com/AcademySoftwareFoundation/rawtoaces

#include "transform_cache.hpp"

#include "../core/spectral_solver.hpp"
#include "../core/metadata_solver.hpp"
#include "../core/matrix_solver.hpp"
#include "data_adapters/spectral_loader.hpp"

#include <OpenImageIO/strutil.h>

namespace rta
{
namespace cache
{

std::ostream &
operator<<( std::ostream &stream, const TransformDescriptor &descriptor )
{
    switch ( descriptor.type )
    {
        case TransformEntryType::Illum_from_WB:
            stream << "illuminant from WB weights:";
            break;
        case TransformEntryType::WB_from_Illum:
            stream << "WB weights from illuminant:";
            break;
        case TransformEntryType::Mat_from_Illum:
            stream << "IDT matrix from illuminant:";
            break;
        case TransformEntryType::Mat_from_DNG:
            stream << "IDT matrix from DNG metadata:";
            break;
        case TransformEntryType::Mat_from_nonDNG:
            stream << "IDT matrix from non-DNG metadata:";
            break;
        case TransformEntryType::CAT_from_WB:
            stream << "CAT matrix from source and destination white point:";
            break;
    }

    stream << std::endl;
    stream << "  camera make: " << descriptor.camera_make << std::endl;
    stream << "  camera model: " << descriptor.camera_model << std::endl;

    switch ( descriptor.type )
    {
        case TransformEntryType::Illum_from_WB: {
            auto &wb = std::get<WB>( descriptor.value );

            stream << "  white balance: (";
            stream << wb.value[0] << ", ";
            stream << wb.value[1] << ", ";
            stream << wb.value[2] << ")";
            break;
        }
        case TransformEntryType::WB_from_Illum:
        case TransformEntryType::Mat_from_Illum: {
            auto &illum = std::get<Illuminant>( descriptor.value );
            stream << "  illuminant: " << illum;
            break;
        }
        case TransformEntryType::Mat_from_DNG:
        case TransformEntryType::Mat_from_nonDNG: break;
        case TransformEntryType::CAT_from_WB: {
            auto wb = std::get<std::pair<WB, WB>>( descriptor.value );
            stream << " source white point ( " << wb.first.value[0] << ", "
                   << wb.first.value[1] << ", " << wb.first.value[2] << "); "
                   << " destination white point ( " << wb.second.value[0]
                   << ", " << wb.second.value[1] << ", " << wb.second.value[2]
                   << ")" << std::endl;
            break;
        }
    }

    stream << std::endl;
    return stream;
}

bool TransformDescriptor::operator==( const TransformDescriptor &other ) const
{
    if ( camera_make != other.camera_make )
        return false;
    if ( camera_model != other.camera_model )
        return false;
    if ( value != other.value )
        return false;
    return true;
}

bool TransformDescriptor::operator==( const CacheEntryDescriptor &other ) const
{
    return false;
}

size_t TransformDescriptor::map_index() const
{
    return (size_t)type;
}

//void prepare_solver(
//    core::SpectralSolver &solver,
//    const std::string    camera_make,
//    const std::string    camera_model,
//    const std::string    illuminant_name,
//    int                  verbosity )
//{
//    solver.verbosity = verbosity;
//
//    bool found_camera = false;
//    auto camera_paths = util::collectDataFiles( "camera" );
//
//    if (!spectral_l::find_camera(camera_make, camera_model, solver.camera))
//    {
//        std::cerr << "Camera spectral sensitivity data not found for "
//                  << camera_make << " " << camera_model << ". "
//                  << "Please check that the data is available "
//                  << "at the location(s) specified in RAWTOACES_DATA_PATH"
//                  << std::endl;
//        exit( 1 );
//    }
//
//    auto training = findFile( "training/training_spectral.json" );
//    if ( training.length() )
//    {
//        solver.training_data = util::load_spectral_file(training);
//    }
//
//    auto cmf = findFile( "cmf/cmf_1931.json" );
//    if ( cmf.length() )
//    {
//        solver.observer = util::load_spectral_file(cmf);
//    }
//
//    if (!util::load_illuminant(illuminant_name, solver.illuminant))
//    {
//        //TODO
//        exit( 1 );
//    }
//}

bool TransformDescriptor::fetch(
    TransformCacheEntryData &data, int verbosity ) const
{
    switch ( type )
    {
        case TransformEntryType::Illum_from_WB: {
            auto               &wb1 = std::get<WB>( value );
            std::vector<double> wb2( 3 );
            wb2[0] = wb1.value[0];
            wb2[1] = wb1.value[1];
            wb2[2] = wb1.value[2];

            util::SpectralLoader spectral_loader;

            core::SpectralSolver solver;
            solver.verbosity = verbosity;

            if ( !spectral_loader.find_camera(
                     camera_make, camera_model, solver.camera ) )
            {
                std::cerr << "Failed to find spectral data for camera "
                          << camera_make << " " << camera_model << std::endl;
                return false;
            }

            solver.illuminant =
                spectral_loader.find_best_illuminant( solver.camera, wb2 );

            auto wb3 = solver.calculate_white_balance(
                solver.illuminant, solver.camera );

            WB wb4;
            wb4.value[0] = wb3[0];
            wb4.value[1] = wb3[1];
            wb4.value[2] = wb3[2];

            data.value =
                std::pair<WB, std::string>( wb4, solver.illuminant.illuminant );
            break;
        }
        case TransformEntryType::WB_from_Illum: {
            auto &illuminant_type = std::get<Illuminant>( value );

            core::SpectralSolver solver;
            solver.verbosity = verbosity;

            util::SpectralLoader spectral_loader;

            if ( !spectral_loader.find_camera(
                     camera_make, camera_model, solver.camera ) )
            {
                std::cerr << "Failed to find spectral data for camera "
                          << camera_make << " " << camera_model << std::endl;
                return false;
            }

            if ( !spectral_loader.load_illuminant(
                     illuminant_type, solver.illuminant ) )
            {
                std::cerr << "Failed to find illuminant " << illuminant_type
                          << "." << std::endl;
                return false;
            }

            auto wb1 = solver.calculate_white_balance(
                solver.illuminant, solver.camera );

            WB wb2;
            wb2.value[0] = wb1[0];
            wb2.value[1] = wb1[1];
            wb2.value[2] = wb1[2];

            std::pair<WB, std::string>( wb2, illuminant_type );
            break;
        }
        case TransformEntryType::Mat_from_Illum: {
            auto       &illuminant_type = std::get<Illuminant>( value );
            std::string lower_illuminant =
                OIIO::Strutil::lower( illuminant_type );

            core::SpectralSolver solver;
            solver.verbosity = verbosity;

            util::SpectralLoader spectral_loader;

            if ( !spectral_loader.find_camera(
                     camera_make, camera_model, solver.camera ) )
            {
                std::cerr << "Failed to find spectral data for camera "
                          << camera_make << " " << camera_model << std::endl;
                return false;
            }

            if ( !spectral_loader.load_illuminant(
                     illuminant_type, solver.illuminant ) )
            {
                std::cerr << "Failed to find illuminant " << illuminant_type
                          << ". "
                          << "Please check the RAWTOACES_DATA_PATH environment "
                          << "variable." << std::endl;
                return false;
            }

            solver.white_balance = solver.calculate_white_balance(
                solver.illuminant, solver.camera );

            auto training_data_path =
                spectral_loader.find_file( "training/training_spectral.json" );
            if ( training_data_path.length() == 0 )
            {
                std::cerr << "Failed to find training data file "
                          << "\"training/training_spectral.json\". "
                          << "Please check the RAWTOACES_DATA_PATH environment "
                          << "variable." << std::endl;
                return false;
            }
            solver.training_data =
                spectral_loader.load_spectral_file( training_data_path );

            auto observer_path =
                spectral_loader.find_file( "cmf/cmf_1931.json" );
            if ( observer_path.length() == 0 )
            {
                std::cerr << "Failed to find the standard observer data file "
                          << "\"cmf/cmf_1931.json\". "
                          << "Please check the RAWTOACES_DATA_PATH environment "
                          << "variable." << std::endl;
                return false;
            }

            solver.observer =
                spectral_loader.load_spectral_file( observer_path, true );

            solver.scale_illuminant( solver.illuminant );

            std::vector<std::vector<double>> mat1;
            bool result = solver.calculate_IDT( mat1 );
            if ( !result )
                return false;

            Matrix33 &mat2 = data.value.emplace<Matrix33>();

            for ( size_t i = 0; i < 3; i++ )
                for ( size_t j = 0; j < 3; j++ )
                    mat2.value[i][j] = mat1[i][j];

            break;
        }
        case TransformEntryType::Mat_from_DNG: {
            const core::Metadata &metadata = std::get<core::Metadata>( value );

            core::MetadataSolver solver( metadata );

            auto      mat1 = solver.calculate_IDT();
            Matrix33 &mat2 = data.value.emplace<Matrix33>();

            for ( size_t i = 0; i < 3; i++ )
                for ( size_t j = 0; j < 3; j++ )
                    mat2.value[i][j] = mat1[i][j];
            break;
        }
        case TransformEntryType::Mat_from_nonDNG: {
            auto p = std::get<std::pair<WB, Matrix33>>( value );

            Vector3                         &wb   = p.first;
            Matrix33                        &mat1 = p.second;
            std::vector<std::vector<double>> xyz2cam;

            xyz2cam.resize( 3 );
            for ( size_t row = 0; row < 3; row++ )
            {
                xyz2cam[row].resize( 3 );

                for ( size_t col = 0; col < 3; col++ )
                {
                    xyz2cam[row][col] = mat1.value[row][col];
                }
            }

            std::vector<std::vector<double>> cam2xyz;
            bool result = core::inverse( xyz2cam, cam2xyz );
            if ( !result )
                return false;

            for ( size_t row = 0; row < 3; row++ )
            {
                for ( size_t col = 0; col < 3; col++ )
                {
                    cam2xyz[row][col] /= wb.value[row];
                }
            }

            Matrix33 &mat2 = data.value.emplace<Matrix33>();

            for ( size_t i = 0; i < 3; i++ )
                for ( size_t j = 0; j < 3; j++ )
                    mat2.value[i][j] = cam2xyz[i][j];

            break;
        }
        case TransformEntryType::CAT_from_WB: {

            auto                wb = std::get<std::pair<WB, WB>>( value );
            std::vector<double> src_XYZ( wb.first.value, wb.first.value + 3 );
            std::vector<double> dst_XYZ( wb.second.value, wb.second.value + 3 );

            auto mat1 = core::calculate_CAT( src_XYZ, dst_XYZ );

            Matrix33 &mat2 = data.value.emplace<Matrix33>();

            for ( size_t i = 0; i < 3; i++ )
                for ( size_t j = 0; j < 3; j++ )
                    mat2.value[i][j] = mat1[i][j];

            break;
        }
    }

    return true;
}

template class cache::
    CacheBase<TransformDescriptor, TransformCacheEntryData, 5>;

} // namespace cache
} // namespace rta
