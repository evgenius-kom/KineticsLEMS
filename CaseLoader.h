#pragma once
#include "ZipLoader.h"
#include "CaseParamsParser.h"
#include <memory>


class CaseLoader
{
public:
    CaseLoader( CaseData& caseData, const std::string& pathToArchive );

    bool load();

private:
    CaseData& caseData_;
    const std::unique_ptr<const ZipLoader> zipLoader_;

};
