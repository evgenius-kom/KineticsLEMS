#pragma once
#include <string>
#include <map>


enum class Extension
{
	TXT,
	ZIP,
	JSON
};

static const std::string TXT_EXTENSION  = ".txt";
static const std::string ZIP_EXTENSION  = ".zip";
static const std::string JSON_EXTENSION = ".json";

static const std::string OUTPUT_FOLDER = "Output";

static const std::string CASE_SETTINGS_FILENAME = "settings";
static const std::string CASE_SETTINGS_FILE     = CASE_SETTINGS_FILENAME + JSON_EXTENSION;

static const std::string APP_SETTINGS_FILENAME = "app_settings";
static const std::string APP_SETTINGS_FILE     = APP_SETTINGS_FILENAME + JSON_EXTENSION;
