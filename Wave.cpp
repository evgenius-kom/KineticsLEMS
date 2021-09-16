#include "Wave.h"


void Wave::addPoint( const Point& point )
{
	points_.emplace_back( point );
}
