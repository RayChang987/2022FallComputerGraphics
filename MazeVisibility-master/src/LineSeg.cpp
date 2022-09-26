/************************************************************************
     File:        LineSeg.cpp

     Author:     
                  Stephen Chenney, schenney@cs.wisc.edu
     Modifier
                  Yu-Chi Lai, yu-chi@cs.wisc.edu

     Comment:    
						(c) 2001-2002 Stephen Chenney, University of Wisconsin at Madison
						Class header file for LineSeg class.
		

     Platform:    Visio Studio.Net 2003 (converted to 2005)

*************************************************************************/

#include "LineSeg.h"

#include <iostream>
using namespace std;
//**********************************************************************
//
// * Constructor from an edge
//======================================================================
Point2D::Point2D(float x, float y) {
	this->x = x; this->y = y;
}
LineSeg::
LineSeg(Edge *e)
//======================================================================
{
	start[0] = e->endpoints[Edge::START]->posn[Vertex::X];
	start[1] = e->endpoints[Edge::START]->posn[Vertex::Y];

	end[0] = e->endpoints[Edge::END]->posn[Vertex::X];
	end[1] = e->endpoints[Edge::END]->posn[Vertex::Y];
	if (end[0] == start[0]) {m = 1e20+9;}
	else {m = (end[1] - start[1]) / (end[0] - start[0]);}
	b = start[1] - m * start[0];
}

//**********************************************************************
//
// * Constructor for specifyng start and end point
//======================================================================
LineSeg::
LineSeg(float xs, float ys, float xe, float ye)
//======================================================================
{
	start[0] = xs;
	start[1] = ys;
	end[0] = xe;
	end[1] = ye;
	if (end[0] == start[0]) { m = 1e20 + 9; }
	else { m = (end[1] - start[1]) / (end[0] - start[0]); }
	b = start[1] - m * start[0];
}


//**********************************************************************
//
// * Return the parameter value at which this segment crosses the given
//   segment. This will return parameter values outside the range 0,1
//   THIS FUNCTION IS EXTREMELY USEFUL FOR CLIPPING, but it 
//   DOES NOT tell you whether the edge is "entering" or "leaving".
//   But you can use tests like Edge::Point_Side() to figure that out.
//======================================================================
Point2D LineSeg::find_intersection(LineSeg l) {
	//parametric form
	float a1, a2, b1, b2, c1, c2;
	a1 = l.m, b1 = -1, c1 = l.b;
	a2 = m, b2 = -1, c2 = b;
	float det = a1 * b2 - a2 * b1;
	return Point2D((b1 * c2 - b2 * c1) / det, (c1*a2-a1*c2) / det);
}
LineSeg::LineSeg(Point2D s, Point2D e) {
	start[0] = s.x; start[1] = s.y;
	end[0] = e.x; end[1] = e.y;
	if (end[0] == start[0]) { m = 1e20 + 9; }
	else { m = (end[1] - start[1]) / (end[0] - start[0]); }
	b = start[1] - m * start[0];
}
float LineSeg::
Cross_Param(LineSeg e)
//======================================================================
{
	float   dx1, dy1, dx2, dy2;
	float   denom, s;

	// This computation comes from writing each segment in parametric form,
	// and solving a simulataneous equation to determine the parameter
	// value of the intersection point on this LineSeg.
	dx1 = e.end[0] - e.start[0];
	dy1 = e.end[1] - e.start[1];
	dx2 = end[0] - start[0];
	dy2 = end[1] - start[1];

	if ( ( denom = dx2 * dy1 - dy2 * dx1 ) == 0.0 )
		// Parallel segments.
		return 1.0e20f;

	s = ( e.start[0] - start[0] ) * dy1 - ( e.start[1] - start[1] ) * dx1;

	return s / denom;
}
char LineSeg::Point_Side(float x, float y) {
	// Compute the determinant: | xs ys 1 |
	//                          | xe ye 1 |
	//                          | x  y  1 |
// Use its sign to get the answer.

	float   det;

	det = start[0] *
		(end[1] - y) -
		start[1] *
		(end[0] - x) +
		end[0] * y -
		end[1] * x;

	if (det == 0.0)
		return Edge::ON;
	else if (det > 0.0)
		return Edge::LEFT;
	else
		return Edge::RIGHT;
}
LineSeg::LineSeg(float* s, float* e) {
	start[0] = s[0];
	start[1] = s[2];
	end[0] = e[0];
	end[1] = e[2];
	if (end[0] == start[0]) { m = 1e20 + 9; }
	else { m = (end[1] - start[1]) / (end[0] - start[0]); }
	b = start[1] - m * start[0];
}
bool LineSeg::onSeg(Point2D p) {
	bool x_mono = (start[0] >= p.x && p.x >= end[0]) || (end[0] >= p.x && p.x >= start[0]);
	bool y_mono = (start[1] >= p.y && p.y >= end[1]) || (end[1] >= p.y && p.y >= start[1]);
	return x_mono && y_mono;
}
void LineSeg::check() {
	cout << "{" << start[0] << " " << start[1] << "}" << "{" << end[0] << " " << end[1] << "}" << endl;
}