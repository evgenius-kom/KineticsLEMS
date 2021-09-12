#include <filesystem>
#include <string>
#include <map>


enum class ExpType
{
	ISO,
	NON_ISO_HEAT,
	NON_ISO_COOL
};

enum class Method
{
	DSC,   // differential scanning calorimetry
	TGA,   // thermogravimetry
	UFSC,  // ultrafast scanning calorimetry on a chip, or nanocalorimetry
	POM    // polarized optical microscopy
};

struct CaseParams
{
	ExpType expType;
	Method method;
	std::string material;
	std::map<string, float> fileToConditionMap;
};

struct CaseData
{
	CaseParams params;
	std::filesystem::path path;
};
