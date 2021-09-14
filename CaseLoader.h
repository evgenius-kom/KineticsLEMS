#pragma once
#include "ZipLoader.h"
#include "CaseParamsParser.h"
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
    CaseLoader( CaseData& caseData, const std::string& pathToArchive );

private:
    CaseData& caseData_;
    const std::unique_ptr<const ZipLoader> zipLoader_;
    const std::filesystem::path pathToFolder_;
    const std::unique_ptr<CaseParamsParser> caseParamsParser_;

};
