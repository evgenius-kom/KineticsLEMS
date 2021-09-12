#pragma once
#include "CaseData.h"
#include "ZipLoader.h"
#include <memory>

/*
    1. Read the archive
    2. Read .txt files inside the archive
    3. Make waves from these .txt files
    4. Parse case settings.json file -> fill CaseData
*/

class CaseLoader
{
public:
	CaseLoader( const std::string& pathToCase, CaseData& caseData );

private:
	const std::unique_ptr<ZipLoader> zipLoader_;
	CaseData& caseData_;

};
