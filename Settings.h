#pragma once
#include "filesystem"
#include "thirdparty/json.hpp"


class Settings
{
public:
	explicit Settings( const std::filesystem::path& path );
	Settings( const Settings& settings );

	nlohmann::json getJson() const;

private:
	nlohmann::json json_;

};
