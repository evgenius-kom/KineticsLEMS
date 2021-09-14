#include "CaseLoader.h"
#include "Settings.h"
#include "Paths.h"
#include <iostream>


CaseLoader::CaseLoader( CaseData& caseData, const std::string& pathToArchive ) :
	caseData_( caseData ),
	zipLoader_( std::make_unique<ZipLoader>( pathToArchive ) ),
	pathToFolder_( zipLoader_->unarchive() ),
	caseParamsParser_( std::make_unique<CaseParamsParser>( caseData.params, 
														   Settings( pathToFolder_ / CASE_SETTINGS_FILE ).getJson() ) )
{
	if ( !caseParamsParser_->parse() )
	{
		throw std::invalid_argument( "Invalid case settings file" );
	}

	// TODO: read waves
}
