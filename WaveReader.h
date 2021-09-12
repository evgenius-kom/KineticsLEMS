#pragma once

#include "Wave.h"
#include <string>


class WaveReader
{
public:
	explicit WaveReader( const std::string& fileName );

private:
	const std::string fileName_;

};
