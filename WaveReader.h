#pragma once
#include "Wave.h"
#include <filesystem>


class WaveReader
{
public:
	explicit WaveReader( const std::filesystem::path& pathToFile );

	bool read();
	Wave getWave() const { return wave_; };

private:
	const std::filesystem::path pathToFile_;
	Wave wave_;

};
