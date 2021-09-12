#include "ZipLoader.h"
#include "WaveReader.h"
#include <stdexcept>
#include <iostream>


ZipLoader::ZipLoader( const std::string& pathToArchive ) : 
    path_( pathToArchive ) 
{
    if ( !std::filesystem::exists( path_ ) || 
         path_.extension() != ".zip" )
    {
        throw "Incorrect path to archive";
    }

    int err = 0;
    zipFile_ = zip_open( pathToArchive.c_str(), 0, &err );
    if ( !zipFile_ ) // err != 0
    {
        throw std::invalid_argument( "Invalid archive" );
    }

    const int filesNum = zip_get_num_files( zipFile_ );
    if ( filesNum == 0 )
    {
        throw std::invalid_argument( "Empty archive" );
    }

    // TODO: make unzip

/*    for ( int i = 0; i < filesNum; ++i ) {
        struct zip_stat file_info;
        zip_stat_init( &file_info );
        zip_stat_index( zipFile_, i, 0, &file_info );

        const std::string& fileName = file_info.name;
        // std::cout << fileName << std::endl;
    };*/
}

ZipLoader::~ZipLoader()
{
    zip_close( zipFile_ );
    // TODO: delete unzipped folder
}

zip* ZipLoader::getZipFile() const
{
    return zipFile_;
}

bool ZipLoader::loadWaves()
{
    return true;
}
