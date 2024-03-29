#include "Settings.h"
#include <iostream>
#include <fstream>


Settings::Settings( const std::filesystem::path& path )
{
	std::ifstream settingsFile;
	settingsFile.open( path );
	if ( settingsFile.good() )
	{
		settingsFile >> json_;
	}
	else
	{
		throw std::invalid_argument( "Invalid settings file" );
	}
}

nlohmann::json Settings::getJson() const
{
	return json_;
}
