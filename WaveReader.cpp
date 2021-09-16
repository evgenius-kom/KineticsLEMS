#include "WaveReader.h"
#include <filesystem>
#include <stdexcept>
#include <sstream>
#include <iostream>
#include <fstream>


WaveReader::WaveReader( const std::filesystem::path& pathToFile ) :
	pathToFile_( pathToFile )
{
    if ( !std::filesystem::exists( pathToFile_ ) )
    {
        throw std::invalid_argument( "Invalid file " + pathToFile_.u8string() );
    }
}

bool WaveReader::read()
{
    std::ifstream file;
    file.open( pathToFile_ );
    if ( !file.good() )
    {
        return false;
    }

    std::string line;
    while ( std::getline( file, line ) )
    {
        std::istringstream iss( line );
        float x, y;
        if ( !( iss >> x >> y ) )
        {
            break;
        }

        Point point { x, y };
        wave_.addPoint( point );
    }
    return !wave_.empty();
}
