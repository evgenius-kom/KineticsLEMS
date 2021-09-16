#pragma once
#include "Point.h"
#include <vector>


class Wave
{
public:
	bool empty() const { return points_.empty(); };
	void addPoint( const Point& point );

	std::vector<Point> points() const { return points_; };
    
    // TODO: overload operators to work with waves

private:
	std::vector<Point> points_;

};

