#include "Wave.h"
#include <filesystem>
#include <string>
#include <map>


enum class ExpType
{
	ISO,
	NON_ISO_HEAT,
	NON_ISO_COOL,
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
	ExpType expType;
	Method method;
	std::string material;
	std::map<std::string, float> fileToConditionMap;
};

struct CaseData
{
	CaseParams params;
	std::filesystem::path path;
	std::map<float, Wave> waves;
};
