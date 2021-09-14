#include "CaseParamsParser.h"


CaseParamsParser::CaseParamsParser( CaseParams& params, const nlohmann::json& json ) :
	params_( params ),
	json_( json )
{
	if ( !json_.contains( SETTINGS_FIELD ) )
	{
		throw "No " + SETTINGS_FIELD + " field in the settings file";
	}
}

bool CaseParamsParser::parse() noexcept
{
	return getMaterial_() && getExpType_() && getMethod_() && getFileToConditionMap_();
}

bool CaseParamsParser::getMaterial_()
{
	if ( json_.contains( MATERIAL_FIELD ) )
	{
		return true;	
	}
	return false;
}

bool CaseParamsParser::getExpType_()
{
	if ( json_.contains( EXP_TYPE_FIELD ) )
	{
		return true;	
	}
	return false;
}

bool CaseParamsParser::getMethod_()
{
	if ( json_.contains( METHOD_FIELD ) )
	{
		return true;	
	}
	return false;
}

bool CaseParamsParser::getFileToConditionMap_()
{
	if ( json_.contains( CONDITIONS_FIELD ) )
	{
		return true;	
	}
	return false;
}
