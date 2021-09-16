#include "CaseLoader.h"
#include <iostream>


bool checkParams( const int argc, char* argv[], std::string& pathToArchive );

int main( int argc, char* argv[] )
{
    std::string pathToArchive;

    if ( !checkParams( argc, argv, pathToArchive ) )
    { 
        return false;
    }

    try
    {
        CaseData caseData;
        if ( !CaseLoader( caseData, pathToArchive ).load() )
        {
            std::cerr << "BAD CASE" << std::endl;
        }

        for ( const auto& [value, wave] : caseData.waves )
        {
            std::cout << "Wave values for a heating rate: " << value << " C/min :\n";
            for ( const Point& point : wave.points() )
            {
                std::cout << "x = " << point.x << " -> y = " << point.y << std::endl;
            }
        }
    }
    catch ( const char* e )
    {
        std::cout << e << std::endl;
    }
    catch ( const std::exception& e )
    {
        std::cout << e.what() << std::endl;
    }

    return 0;
}

bool checkParams( const int argc, char* argv[], std::string& pathToArchive )
{
    if ( argc == 1 )
    {
        std::cerr << "Very few number of input args\n";
        return false;
    }

    // TODO: parse & check args
    pathToArchive = argv[1];
    return true;
}
