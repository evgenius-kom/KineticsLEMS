#pragma once
#include <string>
#include <map>


enum class Extension
{
	TXT,
	ZIP,
	JSON
};

static const std::map<Extension, std::string> extToStrMap = {
	{ Extension::TXT,  ".txt"  },
	{ Extension::ZIP,  ".zip"  },
	{ Extension::JSON, ".json" }
};

static const std::string OUTPUT_FOLDER = "Output";

static const std::string CASE_SETTINGS_FILENAME = "settings";
static const std::string CASE_SETTINGS_FILE     = CASE_SETTINGS_FILENAME + extToStrMap.at( Extension::JSON );

static const std::string APP_SETTINGS_FILENAME = "app_settings";
static const std::string APP_SETTINGS_FILE     = APP_SETTINGS_FILENAME + extToStrMap.at( Extension::JSON );
