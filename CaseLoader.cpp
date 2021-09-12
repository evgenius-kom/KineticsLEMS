#include "CaseLoader.h"


CaseLoader::CaseLoader( const std::string& pathToCase, CaseData& caseData ) :
	zipLoader_( std::make_unique<ZipLoader>( pathToCase ) ),
	caseData_( caseData )
{

}
