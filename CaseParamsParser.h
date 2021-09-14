#pragma once
#include "CaseData.h"
#include "thirdparty/json.hpp"


class CaseParamsParser
{
public:
	CaseParamsParser( CaseParams& params, const nlohmann::json& json );

	bool parse() noexcept;

private:
	CaseParams& params_;
	const nlohmann::json& json_;	

    bool getMaterial_();
    bool getExpType_();
    bool getMethod_();
    bool getFileToConditionMap_();

};
