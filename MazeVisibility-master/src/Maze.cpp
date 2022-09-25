/************************************************************************
     File:        Maze.cpp

     Author:     
                  Stephen Chenney, schenney@cs.wisc.edu
     Modifier
                  Yu-Chi Lai, yu-chi@cs.wisc.edu

     Comment:    
						(c) 2001-2002 Stephen Chenney, University of Wisconsin at Madison

						Class header file for Maze class. Manages the maze.
		

     Platform:    Visio Studio.Net 2003 (converted to 2005)

*************************************************************************/

#include "Maze.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <FL/Fl.h>
#include <FL/fl_draw.h>
#include<iostream> //for testing

using namespace std;

const char Maze::X = 0;
const char Maze::Y = 1;
const char Maze::Z = 2;

const float Maze::BUFFER = 0.1f;


//**********************************************************************
//
// * Constructor for the maze exception
//======================================================================
MazeException::
MazeException(const char *m)
//======================================================================
{
	message = new char[strlen(m) + 4];
	strcpy(message, m);
}


//**********************************************************************
//
// * Constructor to create the default maze
//======================================================================
Maze::
Maze(const int nx, const int ny, const float sx, const float sy)
//======================================================================
{
	// Build the connectivity structure.
	Build_Connectivity(nx, ny, sx, sy);

	// Make edges transparent to create a maze.
	Build_Maze();

	// Set the extents of the maze
	Set_Extents();

	// Default values for the viewer.
	viewer_posn[X] = viewer_posn[Y] = viewer_posn[Z] = 0.0;
	viewer_dir = 0.0;
	viewer_fov = 45.0;

	// Always start on the 0th frame.
	frame_num = 0;
}


//**********************************************************************
//
// * Construtor to read in precreated maze
//======================================================================
Maze::
Maze(const char *filename)
//======================================================================
{
	char    err_string[128];
	FILE    *f;
	int	    i;

	// Open the file
	if ( ! ( f = fopen(filename, "r") ) )
		throw new MazeException("Maze: Couldn't open file");

	// Get the total number of vertices
	if ( fscanf(f, "%d", &num_vertices) != 1 )
		throw new MazeException("Maze: Couldn't read number of vertices");

	// Read in each vertices
	vertices = new Vertex*[num_vertices];
	for ( i = 0 ; i < num_vertices ; i++ ) {
		float x, y;
		if ( fscanf(f, "%g %g", &x, &y) != 2 )	{
			sprintf(err_string, "Maze: Couldn't read vertex number %d", i);
			throw new MazeException(err_string);
		}
		vertices[i] = new Vertex(i, x, y);
	}

	// Get the number of edges
	if ( fscanf(f, "%d", &num_edges) != 1 )
		throw new MazeException("Maze: Couldn't read number of edges");

	// read in all edges
	edges = new Edge*[num_edges];
	for ( i = 0 ; i < num_edges ; i++ ){
		int     vs, ve, cl, cr, o;
		float	r, g, b;
		if ( fscanf(f, "%d %d %d %d %d %g %g %g",
						&vs, &ve, &cl, &cr, &o, &r, &g, &b) != 8) {
			sprintf(err_string, "Maze: Couldn't read edge number %d", i);
			throw new MazeException(err_string);
		}
		edges[i] = new Edge(i, vertices[vs], vertices[ve], r, g, b);
		edges[i]->Add_Cell((Cell*)cl, Edge::LEFT);
		edges[i]->Add_Cell((Cell*)cr, Edge::RIGHT);
		edges[i]->opaque = o ? true : false;
	}

	// Read in the number of cells
	if ( fscanf(f, "%d", &num_cells) != 1 )
		throw new MazeException("Maze: Couldn't read number of cells");


	// Read in all cells
	cells = new Cell*[num_cells];
	for ( i = 0 ; i < num_cells ; i++ )	{
		int epx, epy, emx, emy;
		if ( fscanf(f, "%d %d %d %d", &epx, &epy, &emx, &emy) != 4 ){
			sprintf(err_string, "Maze: Couldn't read cell number %d", i);
			throw new MazeException(err_string);
		}
		cells[i] = new Cell(i, epx >= 0 ? edges[epx] : NULL,
									epy >= 0 ? edges[epy] : NULL,
									emx >= 0 ? edges[emx] : NULL,
									emy >= 0 ? edges[emy] : NULL);
		if ( cells[i]->edges[0] ) {
			if ( cells[i]->edges[0]->neighbors[0] == (Cell*)i )
				cells[i]->edges[0]->neighbors[0] = cells[i];
			else if ( cells[i]->edges[0]->neighbors[1] == (Cell*)i )
				cells[i]->edges[0]->neighbors[1] = cells[i];
			else	{
				sprintf(err_string,
						  "Maze: Cell %d not one of edge %d's neighbors",
							i, cells[i]->edges[0]->index);
				throw new MazeException(err_string);
			}
		}

		if ( cells[i]->edges[1] )	{
			if ( cells[i]->edges[1]->neighbors[0] == (Cell*)i )
				cells[i]->edges[1]->neighbors[0] = cells[i];
			else if ( cells[i]->edges[1]->neighbors[1] == (Cell*)i )
				cells[i]->edges[1]->neighbors[1] = cells[i];
			else {
				sprintf(err_string,
							"Maze: Cell %d not one of edge %d's neighbors",
							i, cells[i]->edges[1]->index);
				throw new MazeException(err_string);
			}
		}
		if ( cells[i]->edges[2] ) {
			if ( cells[i]->edges[2]->neighbors[0] == (Cell*)i )
				cells[i]->edges[2]->neighbors[0] = cells[i];
			else if ( cells[i]->edges[2]->neighbors[1] == (Cell*)i )
				cells[i]->edges[2]->neighbors[1] = cells[i];
			else	{
				sprintf(err_string,
							"Maze: Cell %d not one of edge %d's neighbors",
							i, cells[i]->edges[2]->index);
				throw new MazeException(err_string);
			}
		}
		if ( cells[i]->edges[3] ) {
			if ( cells[i]->edges[3]->neighbors[0] == (Cell*)i )
				cells[i]->edges[3]->neighbors[0] = cells[i];
			else if ( cells[i]->edges[3]->neighbors[1] == (Cell*)i )
				cells[i]->edges[3]->neighbors[1] = cells[i];
			else	{
				sprintf(err_string,
							"Maze: Cell %d not one of edge %d's neighbors",
							i, cells[i]->edges[3]->index);
				throw new MazeException(err_string);
			}
		}
	}

	if ( fscanf(f, "%g %g %g %g %g",
					 &(viewer_posn[X]), &(viewer_posn[Y]), &(viewer_posn[Z]),
					 &(viewer_dir), &(viewer_fov)) != 5 )
		throw new MazeException("Maze: Error reading view information.");

	// Some edges have no neighbor on one side, so be sure to set their
	// pointers to NULL. (They were set at -1 by the save/load process.)
	for ( i = 0 ; i < num_edges ; i++ )	{
		if ( edges[i]->neighbors[0] == (Cell*)-1 )
			edges[i]->neighbors[0] = NULL;
		if ( edges[i]->neighbors[1] == (Cell*)-1 )
			edges[i]->neighbors[1] = NULL;
	}

	fclose(f);

	Set_Extents();

	// Figure out which cell the viewer is in, starting off by guessing the
	// 0th cell.
	Find_View_Cell(cells[0]);

	frame_num = 0;
}


//**********************************************************************
//
// * Destructor must free all the memory allocated.
//======================================================================
Maze::
~Maze(void)
//======================================================================
{
	int i;

	for ( i = 0 ; i < num_vertices ; i++ )
		delete vertices[i];
	delete[] vertices;

	for ( i = 0 ; i < num_edges ; i++ )
		delete edges[i];
	delete[] edges;

	for ( i = 0 ; i < num_cells ; i++ )
		delete cells[i];
	delete[] cells;
}


//**********************************************************************
//
// * Randomly generate the edge's opaque and transparency for an empty maze
//======================================================================
void Maze::
Build_Connectivity(const int num_x, const int num_y,
                   const float sx, const float sy)
//======================================================================
{
	int	i, j, k;
	int edge_i;

	// Ugly code to allocate all the memory for a new maze and to associate
	// edges with vertices and faces with edges.

	// Allocate and position the vertices.
	num_vertices = ( num_x + 1 ) * ( num_y + 1 );
	vertices = new Vertex*[num_vertices];
	k = 0;
	for ( i = 0 ; i < num_y + 1 ; i++ ) {
		for ( j = 0 ; j < num_x + 1 ; j++ )	{
			vertices[k] = new Vertex(k, j * sx, i * sy);
			k++;
		}
	}

	// Allocate the edges, and associate them with their vertices.
	// Edges in the x direction get the first num_x * ( num_y + 1 ) indices,
	// edges in the y direction get the rest.
	num_edges = (num_x+1)*num_y + (num_y+1)*num_x;
	edges = new Edge*[num_edges];
	k = 0;
	for ( i = 0 ; i < num_y + 1 ; i++ ) {
		int row = i * ( num_x + 1 );
		for ( j = 0 ; j < num_x ; j++ ) {
			int vs = row + j;
			int ve = row + j + 1;
			edges[k] = new Edge(k, vertices[vs], vertices[ve],
			rand() / (float)RAND_MAX * 0.5f + 0.25f,
			rand() / (float)RAND_MAX * 0.5f + 0.25f,
			rand() / (float)RAND_MAX * 0.5f + 0.25f);
			k++;
		}
	}

	edge_i = k;
	for ( i = 0 ; i < num_y ; i++ ) {
		int row = i * ( num_x + 1 );
		for ( j = 0 ; j < num_x + 1 ; j++ )	{
			int vs = row + j;
			int ve = row + j + num_x + 1;
			edges[k] = new Edge(k, vertices[vs], vertices[ve],
			rand() / (float)RAND_MAX * 0.5f + 0.25f,
			rand() / (float)RAND_MAX * 0.5f + 0.25f,
			rand() / (float)RAND_MAX * 0.5f + 0.25f);
			k++;
		}
	}

	// Allocate the cells and associate them with their edges.
	num_cells = num_x * num_y;
	cells = new Cell*[num_cells];
	k = 0;
	for ( i = 0 ; i < num_y ; i++ ) {
		int row_x = i * ( num_x + 1 );
		int row_y = i * num_x;
		for ( j = 0 ; j < num_x ; j++ )	{
			int px = edge_i + row_x + 1 + j;
			int py = row_y + j + num_x;
			int mx = edge_i + row_x + j;
			int my = row_y + j;
			cells[k] = new Cell(k, edges[px], edges[py], edges[mx], edges[my]);
			edges[px]->Add_Cell(cells[k], Edge::LEFT);
			edges[py]->Add_Cell(cells[k], Edge::RIGHT);
			edges[mx]->Add_Cell(cells[k], Edge::RIGHT);
			edges[my]->Add_Cell(cells[k], Edge::LEFT);
			k++;
		}
	}
}


//**********************************************************************
//
// * Add edges from cell to the set that are available for removal to
//   grow the maze.
//======================================================================
static void
Add_To_Available(Cell *cell, int *available, int &num_available)
//======================================================================
{
	int i, j;

	// Add edges from cell to the set that are available for removal to
	// grow the maze.

	for ( i = 0 ; i < 4 ; i++ ){
		Cell    *neighbor = cell->edges[i]->Neighbor(cell);

		if ( neighbor && ! neighbor->counter )	{
			int candidate = cell->edges[i]->index;
			for ( j = 0 ; j < num_available ; j++ )
				if ( candidate == available[j] ) {
					printf("Breaking early\n");
					break;
			}
			if ( j == num_available )  {
				available[num_available] = candidate;
				num_available++;
			}
		}
	}

	cell->counter = 1;
}


//**********************************************************************
//
// * Grow a maze by removing candidate edges until all the cells are
//   connected. The edges are not actually removed, they are just made
//   transparent.
//======================================================================
void Maze::
Build_Maze()
//======================================================================
{
	Cell    *to_expand;
	int     index;
	int     *available = new int[num_edges];
	int     num_available = 0;
	int	    num_visited;
	int	    i;

	srand(time(NULL));

	// Choose a random starting cell.
	index = (int)floor((rand() / (float)RAND_MAX) * num_cells);
	to_expand = cells[index];
	Add_To_Available(to_expand, available, num_available);
	num_visited = 1;

	// Join cells up by making edges opaque.
	while ( num_visited < num_cells && num_available > 0 ) {
		int ei;

		index = (int)floor((rand() / (float)RAND_MAX) * num_available);
		to_expand = NULL;

		ei = available[index];

		if ( edges[ei]->neighbors[0] && 
			 !edges[ei]->neighbors[0]->counter )
			to_expand = edges[ei]->neighbors[0];
		else if ( edges[ei]->neighbors[1] && 
			 !edges[ei]->neighbors[1]->counter )
			to_expand = edges[ei]->neighbors[1];

		if ( to_expand ) {
			edges[ei]->opaque = false;
			Add_To_Available(to_expand, available, num_available);
			num_visited++;
		}

		available[index] = available[num_available-1];
		num_available--;
	}

	for ( i = 0 ; i < num_cells ; i++ )
		cells[i]->counter = 0;
}


//**********************************************************************
//
// * Go through all the vertices looking for the minimum and maximum
//   extents of the maze.
//======================================================================
void Maze::
Set_Extents(void)
//======================================================================
{
	int i;

	min_xp = vertices[0]->posn[Vertex::X];
	max_xp = vertices[0]->posn[Vertex::X];
	min_yp = vertices[0]->posn[Vertex::Y];
	max_yp = vertices[0]->posn[Vertex::Y];
	for ( i = 1 ; i < num_vertices ; i++ ) {
		if ( vertices[i]->posn[Vertex::X] > max_xp )
			 max_xp = vertices[i]->posn[Vertex::X];
		if ( vertices[i]->posn[Vertex::X] < min_xp )
			 min_xp = vertices[i]->posn[Vertex::X];
		if ( vertices[i]->posn[Vertex::Y] > max_yp )
			 max_yp = vertices[i]->posn[Vertex::Y];
		if ( vertices[i]->posn[Vertex::Y] < min_yp )
			 min_yp = vertices[i]->posn[Vertex::Y];
    }
}


//**********************************************************************
//
// * Figure out which cell the view is in, using seed_cell as an
//   initial guess. This procedure works by repeatedly checking
//   whether the viewpoint is in the current cell. If it is, we're
//   done. If not, Point_In_Cell returns in new_cell the next cell
//   to test. The new cell is the one on the other side of an edge
//   that the point is "outside" (meaning that it might be inside the
//   new cell).
//======================================================================
void Maze::
Find_View_Cell(Cell *seed_cell)
//======================================================================
{
	Cell    *new_cell;

	// 
	while ( ! ( seed_cell->Point_In_Cell(viewer_posn[X], viewer_posn[Y],
													 viewer_posn[Z], new_cell) ) ) {
		if ( new_cell == 0 ) {
			// The viewer is outside the top or bottom of the maze.
			throw new MazeException("Maze: View not in maze\n");
		}

		seed_cell = new_cell;
    }
    
    view_cell = seed_cell;
}


//**********************************************************************
//
// * Move the viewer's position. This method will do collision detection
//   between the viewer's location and the walls of the maze and prevent
//   the viewer from passing through walls.
//======================================================================
void Maze::
Move_View_Posn(const float dx, const float dy, const float dz)
//======================================================================
{
	Cell    *new_cell;
	float   xs, ys, zs, xe, ye, ze;

	// Move the viewer by the given amount. This does collision testing to
	// prevent walking through walls. It also keeps track of which cells the
	// viewer is in.

	// Set up a line segment from the start to end points of the motion.
	xs = viewer_posn[X];
	ys = viewer_posn[Y];
	zs = viewer_posn[Z];
	xe = xs + dx;
	ye = ys + dy;
	ze = zs + dz;

	// Fix the z to keep it in the maze.
	if ( ze > 1.0f - BUFFER )
		ze = 1.0f - BUFFER;
	if ( ze < BUFFER - 1.0f )
		ze = BUFFER - 1.0f;

	// Clip_To_Cell clips the motion segment to the view_cell if the
	// segment intersects an opaque edge. If the segment intersects
	// a transparent edge (through which it can pass), then it clips
	// the motion segment so that it _starts_ at the transparent edge,
	// and it returns the cell the viewer is entering. We keep going
	// until Clip_To_Cell returns NULL, meaning we've done as much of
	// the motion as is possible without passing through walls.
	while ( ( new_cell = view_cell->Clip_To_Cell(xs, ys, xe, ye, BUFFER) ) )
		view_cell = new_cell;

	// The viewer is at the end of the motion segment, which may have
	// been clipped.
	viewer_posn[X] = xe;
	viewer_posn[Y] = ye;
	viewer_posn[Z] = ze;
}

//**********************************************************************
//
// * Set the viewer's location 
//======================================================================
void Maze::
Set_View_Posn(float x, float y, float z)
//======================================================================
{
	// First make sure it's in some cell.
	// This assumes that the maze is rectangular.
	if ( x < min_xp + BUFFER )
		x = min_xp + BUFFER;
	if ( x > max_xp - BUFFER )
		x = max_xp - BUFFER;
	if ( y < min_yp + BUFFER )
		y = min_yp + BUFFER;
	if ( y > max_yp - BUFFER )
		y = max_yp - BUFFER;
	if ( z < -1.0f + BUFFER )
		z = -1.0f + BUFFER;
	if ( z > 1.0f - BUFFER )
		z = 1.0f - BUFFER;

	viewer_posn[X] = x;
	viewer_posn[Y] = y;
	viewer_posn[Z] = z;

	// Figure out which cell we're in.
	Find_View_Cell(cells[0]);
}


//**********************************************************************
//
// * Set the angle in which the viewer is looking.
//======================================================================
void Maze::
Set_View_Dir(const float d)
//======================================================================
{
	viewer_dir = d;
}


//**********************************************************************
//
// * Set the horizontal field of view.
//======================================================================
void Maze::
Set_View_FOV(const float f)
//======================================================================
{
	viewer_fov = f;
}


//**********************************************************************
//
// * Draws the map view of the maze. It is passed the minimum and maximum
//   corners of the window in which to draw.
//======================================================================
void Maze::
Draw_Map(int min_x, int min_y, int max_x, int max_y)
//======================================================================
{
	int	    height;
	float   scale_x, scale_y, scale;
	int	    i;

	// Figure out scaling factors and the effective height of the window.
	scale_x = ( max_x - min_x - 10 ) / ( max_xp - min_xp );
	scale_y = ( max_y - min_y - 10 ) / ( max_yp - min_yp );
	scale = scale_x > scale_y ? scale_y : scale_x;
	height = (int)ceil(scale * ( max_yp - min_yp ));

	min_x += 5;
	min_y += 5;

	// Draw all the opaque edges.
	for ( i = 0 ; i < num_edges ; i++ )
		if ( edges[i]->opaque )	{
			float   x1, y1, x2, y2;

			x1 = edges[i]->endpoints[Edge::START]->posn[Vertex::X];
			y1 = edges[i]->endpoints[Edge::START]->posn[Vertex::Y];
			x2 = edges[i]->endpoints[Edge::END]->posn[Vertex::X];
			y2 = edges[i]->endpoints[Edge::END]->posn[Vertex::Y];

			fl_color((unsigned char)floor(edges[i]->color[0] * 255.0),
					 (unsigned char)floor(edges[i]->color[1] * 255.0),
					 (unsigned char)floor(edges[i]->color[2] * 255.0));
			fl_line_style(FL_SOLID);
			fl_line(min_x + (int)floor((x1 - min_xp) * scale),
					  min_y + height - (int)floor((y1 - min_yp) * scale),
					  min_x + (int)floor((x2 - min_xp) * scale),
					  min_y + height - (int)floor((y2 - min_yp) * scale));
		}
}


//**********************************************************************
//
// * Draws the first-person view of the maze. It is passed the focal distance.
//   THIS IS THE FUINCTION YOU SHOULD MODIFY.
//======================================================================

void matrix_mul_4x4(float* a, float* b, float* c) {
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			c[i * 4 + j] = 0;
			for (int k = 0; k < 4; ++k) {
				c[i * 4 + j] += a[i * 4 + k] * b[k * 4 + j];
			}
			//cout << T[i * 4 + j] << " ";
		}
		//cout << endl;
	}
}
void input4x4(float* mat, float* x, float* y) {
	for (int i = 0; i < 4; ++i) {
		y[i] = 0;
		for (int j = 0; j < 4; ++j) {
			y[i] += mat[i * 4 + j] * x[j];
		}
	}
}
float t_modelview[16];
float t_viewerview[16] = { 2.4142, 0, 0, 0, 0, 2.4142, 0, 0,0,0,-1.0001, -0.0200, 0, 0, -1.0000, 0 }; //相機內部算出來的

bool clip(LineSeg frustum_side, float* start, float* end) {
	char s_side = frustum_side.Point_Side(start[0], start[2]);
	char e_side = frustum_side.Point_Side(end[0], end[2]);

	if (s_side == Edge::RIGHT) {
		if (e_side == Edge::LEFT) {
			float percent = frustum_side.Cross_Param(LineSeg(start, end));
			end[0] = frustum_side.start[0] + (frustum_side.end[0] - frustum_side.start[0]) * percent;
			end[2] = frustum_side.start[1] + (frustum_side.end[1] - frustum_side.start[1]) * percent;
		}
	}
	else if (e_side == Edge::RIGHT) {
		float percent = frustum_side.Cross_Param(LineSeg(start, end));
		start[0] = frustum_side.start[0] + (frustum_side.end[0] - frustum_side.start[0]) * percent;
		start[2] = frustum_side.start[1] + (frustum_side.end[1] - frustum_side.start[1]) * percent;
	}
	else {
		return false;
	}
	return true;
}
#define square(x) ((x)*(x))
void Normalize3x1(float* x) {
	float Norm = sqrt(square(x[0]) + square(x[1]) + square(x[2]));
	x[0] /= Norm;
	x[1] /= Norm;
	x[2] /= Norm;
}

void ComputeNormalOfPlane(float* result, float* A, float* B) {
	result[0] = A[1] * B[2] - A[2] * B[1];
	result[1] = A[2] * B[0] - A[0] * B[2];
	result[2] = A[0] * B[1] - A[1] * B[0];
}

void check44(float* a) {
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			cout << a[i * 4 + j] << " ";
		}
		cout << endl;
	}
}
void transpose44(float* a) {
	float t[16];
	for (int i = 0; i < 16; ++i) {t[i] = a[i];}
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			a[i * 4 + j] = t[j * 4 + i];
		}
	}
}

void Maze::LookAt(float posX, float posY, float posZ, float centerX, float centerY, float centerZ, float upX, float upY, float upZ) {
	float center[3] = { centerX, centerY, centerZ };
	float pos[3] = { posX, posY, posZ };
	float up[3] = { upX, upY, upZ };
	float L[3], S[3], U[3];
	float M[16];
	float T[16];
	for (int i = 0; i < 3; ++i) {
		L[i] = center[i] - pos[i];
	}
	Normalize3x1(L);
	
	//for (int i = 0; i < 3; ++i) { cout << L[i] << " "; }cout << endl;
	ComputeNormalOfPlane(S, L, up);
	//for (int i = 0; i < 3; ++i) { cout << S[i] << ","; }cout << endl;
	//for (int i = 0; i < 3; ++i) { cout << L[i] << ","; }cout << endl;
	//for (int i = 0; i < 3; ++i) { cout << up[i] << ","; }cout << endl;
	Normalize3x1(S);
	ComputeNormalOfPlane(U, S, L);
	//for (int i = 0; i < 3; ++i) { cout << U[i] << ","; }cout << endl;
	
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			if (i == j) { T[i * 4 + j] = 1.0; }
			else { T[i * 4 + j] = 0.0; }
		}
	}
	T[3] = -posX, T[7] = -posY, T[11] = -posZ;
	M[0] = S[0]; M[4] = S[1]; M[8] = S[2]; M[12] = 0.0;
	M[1] = up[0]; M[5] = up[1]; M[9] = up[2]; M[13] = 0.0;
	M[2] = -L[0]; M[6] = -L[1]; M[10] = -L[2]; M[14] = 0.0; //negative z
	M[3] = M[7] = M[11] = 0.0;     M[15] = 1.0;
	transpose44(M);
	
	matrix_mul_4x4(M, T, t_modelview);


	//check44(T); 
	//check44(M); 
	//check44(t_modelview);
	//system("PAUSE");
}
void Maze::Perspective(float fov, float aspect, float zNear,float zFar) {
	for (int i = 0; i < 16; ++i) { t_viewerview[i] = 0; }
	float f = 1.0 / tan(To_Radians(fov/2));
	t_viewerview[0] = f/ aspect;
	t_viewerview[5] = f;
	t_viewerview[10] = (zFar + zNear) / (zNear - zFar);
	t_viewerview[11] = (2 * zFar * zNear) / (zNear - zFar);
	t_viewerview[14] = -1;

}
void Maze::
Draw_View(const float aspect)
//======================================================================
{

	frame_num++;
	//###################################################################
	// TODO
	// The rest is up to you!
	//###################################################################
	

	//glEnable(GL_DEPTH_TEST);
	//cout << num_edges << endl;
	//float t_modelviewGL[16], t_viewerviewGL[16];
	

	glGetFloatv(GL_MODELVIEW_MATRIX, t_modelviewGL);
	glGetFloatv(GL_PROJECTION_MATRIX, t_viewerviewGL);
	transpose44(t_modelviewGL);
	transpose44(t_viewerviewGL);
	//for (int i = 0; i < 4; ++i) { cout << t_modelviewGL[i] << " "; }cout << endl;
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	float viewer_pos[3] = { viewer_posn[Maze::Y], 0.0f, viewer_posn[Maze::X] };
	LookAt(viewer_pos[Maze::X], viewer_pos[Maze::Y], viewer_pos[Maze::Z],
		viewer_pos[Maze::X] + sin(Maze::To_Radians(viewer_dir)), viewer_pos[Maze::Y], viewer_pos[Maze::Z] + cos(Maze::To_Radians(viewer_dir)),
		0.0,1.0, 0.0
	);
	Perspective(viewer_fov, aspect, 0.01, 200);
	for (int i = 0; i < num_cells; i++) { cells[i]->vis = 0; }
	draw_cell(view_cell,
		LineSeg(z_near * tan(To_Radians(viewer_fov * 0.5f)), -z_near, z_far * tan(To_Radians(viewer_fov * 0.5f)), -z_far),
		LineSeg(-z_far * tan(To_Radians(viewer_fov * 0.5f)), -z_far, -z_near * tan(To_Radians(viewer_fov * 0.5f)), -z_near));
}


//**********************************************************************
//
// * Draws the frustum on the map view of the maze. It is passed the
//   minimum and maximum corners of the window in which to draw.
//======================================================================
void Maze::
Draw_Frustum(int min_x, int min_y, int max_x, int max_y)
//======================================================================
{
	int	  height;
	float   scale_x, scale_y, scale;
	float   view_x, view_y;

	// Draws the view frustum in the map. Sets up all the same viewing
	// parameters as draw().
	scale_x	= ( max_x - min_x - 10 ) / ( max_xp - min_xp );
	scale_y	= ( max_y - min_y - 10 ) / ( max_yp - min_yp );
	scale		= scale_x > scale_y ? scale_y : scale_x;
	height	= (int)ceil(scale * ( max_yp - min_yp ));

	min_x += 5;
	min_y += 5;

	view_x = ( viewer_posn[X] - min_xp ) * scale;
	view_y = ( viewer_posn[Y] - min_yp ) * scale;
	fl_line(min_x + (int)floor(view_x + 
			  cos(To_Radians(viewer_dir+viewer_fov / 2.0)) * scale),
			  min_y + height- 
			  (int)floor(view_y + 
							 sin(To_Radians(viewer_dir+viewer_fov / 2.0)) * 
							 scale),
				min_x + (int)floor(view_x),
				min_y + height - (int)floor(view_y));
	fl_line(min_x + (int)floor(view_x + 
										cos(To_Radians(viewer_dir-viewer_fov / 2.0))	* 
										scale),
				min_y + height- 
				(int)floor(view_y + sin(To_Radians(viewer_dir-viewer_fov / 2.0)) *
				scale),
				min_x + (int)floor(view_x),
				min_y + height - (int)floor(view_y));
	}


//**********************************************************************
//
// * Draws the viewer's cell and its neighbors in the map view of the maze.
//   It is passed the minimum and maximum corners of the window in which
//   to draw.
//======================================================================
void Maze::
Draw_Neighbors(int min_x, int min_y, int max_x, int max_y)
//======================================================================
{
	int	    height;
	float   scale_x, scale_y, scale;
	int	    i, j;

	// Draws the view cell and its neighbors in the map. This works
	// by drawing just the neighbor's edges if there is a neighbor,
	// otherwise drawing the edge. Every edge is shared, so drawing the
	// neighbors' edges also draws the view cell's edges.

	scale_x = ( max_x - min_x - 10 ) / ( max_xp - min_xp );
	scale_y = ( max_y - min_y - 10 ) / ( max_yp - min_yp );
	scale = scale_x > scale_y ? scale_y : scale_x;
	height = (int)ceil(scale * ( max_yp - min_yp ));

	min_x += 5;
	min_y += 5;

	for ( i = 0 ; i < 4 ; i++ )   {
		Cell	*neighbor = view_cell->edges[i]->Neighbor(view_cell);

		if ( neighbor ){
			for ( j = 0 ; j < 4 ; j++ ){
				Edge    *e = neighbor->edges[j];

				if ( e->opaque )	{
					float   x1, y1, x2, y2;

					x1 = e->endpoints[Edge::START]->posn[Vertex::X];
					y1 = e->endpoints[Edge::START]->posn[Vertex::Y];
					x2 = e->endpoints[Edge::END]->posn[Vertex::X];
					y2 = e->endpoints[Edge::END]->posn[Vertex::Y];

					fl_color((unsigned char)floor(e->color[0] * 255.0),
							  (unsigned char)floor(e->color[1] * 255.0),
							  (unsigned char)floor(e->color[2] * 255.0));
					fl_line_style(FL_SOLID);
					fl_line( min_x + (int)floor((x1 - min_xp) * scale),
							 min_y + height - (int)floor((y1 - min_yp) * scale),
							 min_x + (int)floor((x2 - min_xp) * scale),
							 min_y + height - (int)floor((y2 - min_yp) * scale));
				}
			}
		}
		else {
			Edge    *e = view_cell->edges[i];

			if ( e->opaque ){
				float   x1, y1, x2, y2;

				x1 = e->endpoints[Edge::START]->posn[Vertex::X];
				y1 = e->endpoints[Edge::START]->posn[Vertex::Y];
				x2 = e->endpoints[Edge::END]->posn[Vertex::X];
				y2 = e->endpoints[Edge::END]->posn[Vertex::Y];

				fl_color((unsigned char)floor(e->color[0] * 255.0),
							 (unsigned char)floor(e->color[1] * 255.0),
							 (unsigned char)floor(e->color[2] * 255.0));
				fl_line_style(FL_SOLID);
				fl_line(min_x + (int)floor((x1 - min_xp) * scale),
							min_y + height - (int)floor((y1 - min_yp) * scale),
							min_x + (int)floor((x2 - min_xp) * scale),
							min_y + height - (int)floor((y2 - min_yp) * scale));
			 }
		}
	}
}

void	Maze::draw_cell(Cell* current_cell, LineSeg L, LineSeg R){
	current_cell->vis = true;
	LineSeg hori(R.end[0], R.end[1], L.start[0], L.start[1]);
	//cout << "Compare the model view matrix" << endl;
	//for (int i = 0; i < 16; ++i) { cout << t_modelview[i]<< " "; }cout << endl;
	//for (int i = 0; i < 16; ++i) {cout << t_modelviewGL[i] << " "; }cout << endl;
	//cout << endl;
	//cout << "Compare the projection matrix" << endl;
	//for (int i = 0; i < 16; ++i) { cout << t_viewerview[i] << " "; }cout << endl;
	//for (int i = 0; i < 16; ++i) { cout << t_viewerviewGL[i] << " "; }cout << endl;
	//cout << endl;
	for (int i = 0; i < 4; ++i) {
		LineSeg wall(current_cell->edges[i]);
		float start[4] = { wall.start[1], 1.0f, wall.start[0], 1.0f }, end[4] = { wall.end[1], 1.0f, wall.end[0], 1.0f };
		float start_modelview[4], end_modelview[4], start_projview[4], end_projview[4];
		input4x4(t_modelview, start, start_modelview);
		input4x4(t_modelview, end, end_modelview);
		//for (int i = 0; i < 16; ++i) { cout << t_modelview[i]<< " "; }cout << endl;
		//for (int i = 0; i < 16; ++i) {cout << t_modelviewGL[i] << " "; }cout << endl;
		//cout << endl;
		if (!clip(L, start_modelview, end_modelview) || !clip(R, start_modelview, end_modelview)) continue;
		wall.start[0] = start_modelview[0];
		wall.start[1] = start_modelview[2];
		wall.end[0] = end_modelview[0]; //update the clipped wall
		wall.end[1] = end_modelview[2];
		if (current_cell->edges[i]->opaque) {
			if (!clip(hori, start_modelview, end_modelview)) continue;
			
			input4x4(t_viewerview, start_modelview, start_projview);
			input4x4(t_viewerview, end_modelview, end_projview);
			//if (start_projview[3] < z_near && end_projview[3] < z_near) continue;
			for (int i = 0; i < 3; ++i) {
				start_projview[i] /= start_projview[3];
				end_projview[i] /= end_projview[3];
			}
			//system("PAUSE");
			glBegin(GL_POLYGON);
			glColor3fv(current_cell->edges[i]->color);
			glVertex2f(start_projview[0], start_projview[1]);
			glVertex2f(end_projview[0], end_projview[1]);
			glVertex2f(end_projview[0], -end_projview[1]);
			glVertex2f(start_projview[0], -start_projview[1]);
			glEnd();
		}
		else {
			if (current_cell->edges[i]->Neighbor(current_cell) == NULL) continue;
			static float pre_Lx = 0.0f, pre_Rx = 0.0f, pre_Ly = 0.0f, pre_Ry = 0.0f; //因下面的if else if 可能都不會進去，所以要讓之前的值維持，才不會沒初始化變數就呼叫function
			float Lx = pre_Lx, Rx = pre_Rx, Ly = pre_Ly, Ry = pre_Ry;

			LineSeg midline(0.0f, 0.0f, (wall.start[0] + wall.end[0]) * 0.5, (wall.start[1] + wall.end[1]) * 0.5);
			if (midline.Point_Side(wall.start[0], wall.start[1]) == Edge::LEFT && midline.Point_Side(wall.end[0], wall.end[1]) == Edge::RIGHT) {
				Lx = wall.start[0];
				Ly = wall.start[1];
				Rx = wall.end[0];
				Ry = wall.end[1];
			}
			else if (midline.Point_Side(wall.start[0], wall.start[1]) == Edge::RIGHT && midline.Point_Side(wall.end[0], wall.end[1]) == Edge::LEFT) {
				Lx = wall.end[0];
				Ly = wall.end[1];
				Rx = wall.start[0];
				Ry = wall.start[1];
			}
			pre_Lx = Lx;
			pre_Rx = Rx;
			pre_Ly = Ly;
			pre_Ry = Ry;
			LineSeg newL(Lx, Ly, Lx / Ly * -z_far, -z_far);
			LineSeg newR(Rx / Ry * -z_far, -z_far, Rx, Ry);
			if (!current_cell->edges[i]->Neighbor(current_cell)->vis && fabs((Lx / Ly * -z_far) - (Rx / Ry * -z_far)) > z_near) {
				draw_cell(current_cell->edges[i]->Neighbor(current_cell), newL, newR);
			}
		}
	}
}
//**********************************************************************
//
// * Save the maze to a file of the given name.
//======================================================================
bool Maze::
Save(const char *filename)
//======================================================================
{
	FILE    *f = fopen(filename, "w");
	int	    i;

	// Dump everything to a file of the given name. Returns false if it
	// couldn't open the file. True otherwise.

	if ( ! f )  {
		return false;
   }

	fprintf(f, "%d\n", num_vertices);
	for ( i = 0 ; i < num_vertices ; i++ )
		fprintf(f, "%g %g\n", vertices[i]->posn[Vertex::X],
			      vertices[i]->posn[Vertex::Y]);

		fprintf(f, "%d\n", num_edges);
	for ( i = 0 ; i < num_edges ; i++ )
	fprintf(f, "%d %d %d %d %d %g %g %g\n",
				edges[i]->endpoints[Edge::START]->index,
				edges[i]->endpoints[Edge::END]->index,
				edges[i]->neighbors[Edge::LEFT] ?
				edges[i]->neighbors[Edge::LEFT]->index : -1,
				edges[i]->neighbors[Edge::RIGHT] ?
				edges[i]->neighbors[Edge::RIGHT]->index : -1,
				edges[i]->opaque ? 1 : 0,
				edges[i]->color[0], edges[i]->color[1], edges[i]->color[2]);

	fprintf(f, "%d\n", num_cells);
	for ( i = 0 ; i < num_cells ; i++ )
		fprintf(f, "%d %d %d %d\n",
					cells[i]->edges[0] ? cells[i]->edges[0]->index : -1,
					cells[i]->edges[1] ? cells[i]->edges[1]->index : -1,
					cells[i]->edges[2] ? cells[i]->edges[2]->index : -1,
					cells[i]->edges[3] ? cells[i]->edges[3]->index : -1);

	   fprintf(f, "%g %g %g %g %g\n",
					viewer_posn[X], viewer_posn[Y], viewer_posn[Z],
					viewer_dir, viewer_fov);

	fclose(f);

	return true;
}
