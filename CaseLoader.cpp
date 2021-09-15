#include "CaseLoader.h"
#include "Settings.h"
#include "Paths.h"
#include <iostream>


CaseLoader::CaseLoader( CaseData& caseData, const std::string& pathToArchive ) :
	caseData_( caseData ),
	zipLoader_( std::make_unique<ZipLoader>( pathToArchive ) )
{
}

bool CaseLoader::load()
{
	const std::filesystem::path& pathToFolder = zipLoader_->unarchive();

	const Settings settings( pathToFolder / CASE_SETTINGS_FILE );
	if ( !CaseParamsParser( caseData_.params, settings.getJson() ).parse() )
	{
		return false;
	}

	// TODO: read waves
	return true;
}
