#pragma once

#include "Point.h"
#include <vector>


class Wave
{
public:
	Wave();
	Wave( const std::vector<Point>& points ) : points_( points ) {}

    // TODO: overload operators to work with waves

private:
	std::vector<Point> points_;

};

