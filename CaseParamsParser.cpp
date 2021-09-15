#include "CaseParamsParser.h"
#include <iostream>


CaseParamsParser::CaseParamsParser( CaseParams& params, const nlohmann::json& json ) :
	params_( params ),
	json_( json )
{
	if ( !json_.contains( SETTINGS_FIELD ) )
	{
		throw "No " + SETTINGS_FIELD + " field in the settings file";
	}

	if ( !json_[SETTINGS_FIELD].is_array() )
	{
		throw SETTINGS_FIELD + " : is not an array";
	}
}

bool CaseParamsParser::parse()
{
	bool isMaterialFound    = false;
	bool isExpTypeFound     = false; 
	bool isMethodFound      = false;
	bool areConditionsFound = false; 

	for ( const auto& item : json_[SETTINGS_FIELD] )
	{
		if ( item["name"] == MATERIAL_FIELD )
		{
			params_.material = item["value"];
			isMaterialFound = true;
		}
		else if ( item["name"] == EXP_TYPE_FIELD )
		{
			if      ( item["value"] == "isothermal" ) { params_.expType = ExpType::ISOTHERMAL; }
			else if ( item["value"] == "heating" )    { params_.expType = ExpType::HEATING; }
			else if ( item["value"] == "cooling" )    { params_.expType = ExpType::COOLING; }
			else { std::cout << "Invalid experiment type found\n"; return false; }
			isExpTypeFound = true;
		}
		else if ( item["name"] == METHOD_FIELD )
		{
			if      ( item["value"] == "DSC" )  { params_.method = Method::DSC; }
			else if ( item["value"] == "TGA" )  { params_.method = Method::TGA; }
			else if ( item["value"] == "UFSC" ) { params_.method = Method::UFSC; }
			else if ( item["value"] == "POM" )  { params_.method = Method::POM; }
			else { std::cout << "Invalid method found\n"; return false; }
			isMethodFound = true;
		}
		else if ( item["name"] == CONDITIONS_FIELD )
		{
			// TODO: parse
			areConditionsFound = true;
		}
		else
		{
			std::cout << "Invalid case setting type found\n";
		}
	}
	return isMaterialFound && isExpTypeFound && isMethodFound && areConditionsFound;
}
