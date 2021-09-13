#pragma once
#include "CaseData.h"
#include "ZipLoader.h"
#include <memory>

/*
    1. Read the archive +
    2. Unzip files +
    3. Read .txt files inside the archive +
    4. Make waves from these .txt files
    5. Parse case settings.json file -> fill CaseData
*/

class CaseLoader
{
public:
	CaseLoader( const std::string& pathToCase, CaseData& caseData );

private:
	const std::unique_ptr<ZipLoader> zipLoader_;
	CaseData& caseData_;

};
