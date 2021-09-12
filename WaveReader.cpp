#include "WaveReader.h"


WaveReader::WaveReader( const std::string& fileName ) :
	fileName_( fileName )
{
    //Alloc memory for its uncompressed contents
/*    char *contents = new char[st.size];

    //Read the compressed file
    zip_file* f = zip_fopen( z, fileName_, 0 );
    zip_fread( f, contents, st.size );
    zip_fclose( f );

    std::cout << contents << std::endl;
    delete[] contents;
*/}
