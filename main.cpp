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
        CaseLoader caseLoader( caseData, pathToArchive );
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
