#ifndef GRID_H
#define GRID_H
#include<vector>

class Grid
{

// n_x, n_y = total number of intervals
// h_x, h_y = node spacing in x and y
public:

unsigned int n_x,n_y;
double h_x,h_y;

///u approximation vector
double *vec;
unsigned int size=0; //vec is of zero length 

///Constructor Declarations
Grid(unsigned int size_x,unsigned int size_y);
double& operator ()(unsigned int x,unsigned int y);

void release(){delete [] vec;}
};

#endif
