#include "ZipLoader.h"
#include "Paths.h"
#include <stdexcept>
#include <iostream>
#include <fstream>


ZipLoader::ZipLoader( const std::string& pathToArchive ) : 
    pathToZip_( pathToArchive ) 
{
    if ( !std::filesystem::exists( pathToZip_ ) || 
         pathToZip_.extension() != ZIP_EXTENSION )
    {
        throw "Incorrect path to archive";
    }

    int err = 0;
    zipFile_ = zip_open( pathToArchive.c_str(), 0, &err );
    if ( !zipFile_ ) // err != 0
    {
        throw std::invalid_argument( "Invalid archive" );
    }

    if ( zip_get_num_files( zipFile_ ) == 0 )
    {
        throw std::invalid_argument( "Empty archive" );
    }
}

ZipLoader::~ZipLoader()
{
    zip_close( zipFile_ );
}

std::filesystem::path ZipLoader::unarchive() const
{
    const std::filesystem::path& pathToFolder = pathToZip_.parent_path() / pathToZip_.stem();
    if ( std::filesystem::exists( pathToFolder ) )
    {
        std::filesystem::remove_all( pathToFolder );
    }

    std::filesystem::create_directory( pathToFolder );

    const int filesNum = zip_get_num_files( zipFile_ );
    for ( int i = 0; i < filesNum; ++i ) {
        struct zip_stat fileInfo;
        zip_stat_init( &fileInfo );
        zip_stat_index( zipFile_, i, 0, &fileInfo );
        std::cout << fileInfo.name << std::endl;

        std::ofstream newFile( pathToFolder / fileInfo.name ); // TODO: create the File class

        struct zip_file* f = zip_fopen_index( zipFile_, i, 0 );
        if ( f ) {
            char* contents = new char[fileInfo.size];
            zip_fread( f, contents, fileInfo.size );
            zip_fclose( f );
            newFile << contents;
            // TODO: check the end of files
            delete[] contents;
        }

        newFile.close();
    };

    return pathToFolder;
}
