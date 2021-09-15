#include "Wave.h"
#include <filesystem>
#include <string>
#include <map>


static const std::string VERSION_FIELD  = "Version";
static const std::string SETTINGS_FIELD = "Settings";

// fields in the app_settings.json file
static const std::string CASE_PATH_FIELD = "PathToArchive";

// fields in the case settings.json file
static const std::string MATERIAL_FIELD   = "Material";
static const std::string EXP_TYPE_FIELD   = "ExperimentType";
static const std::string METHOD_FIELD     = "Method";
static const std::string CONDITIONS_FIELD = "Conditions";


enum class X
{
	TEMPERATURE,
	TIME
};

enum class Y
{
	HEAT_FLOW,
	MASS,
	OPT_DENSITY
};

enum class ExpType
{
	ISOTHERMAL,
	HEATING,
	COOLING,
	INVALID
};

enum class Method
{
	DSC,    // differential scanning calorimetry
	TGA,    // thermogravimetry
	UFSC,   // ultrafast scanning calorimetry on a chip, or nanocalorimetry
	POM,    // polarized optical microscopy
	INVALID
};

struct CaseParams
{
	std::string material;
	ExpType expType;
	Method method;
	std::map<std::string, float> fileToConditionMap;
};

struct CaseData
{
	CaseParams params;
	std::filesystem::path path;
	std::map<float, Wave> waves;
};
