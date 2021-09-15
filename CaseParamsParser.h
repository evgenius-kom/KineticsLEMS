#pragma once
#include "CaseData.h"
#include "thirdparty/json.hpp"


class CaseParamsParser
{
public:
	CaseParamsParser( CaseParams& params, const nlohmann::json& json );

	bool parse();

private:
	CaseParams& params_;
	const nlohmann::json& json_;

};
