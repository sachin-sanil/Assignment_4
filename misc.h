#include "grid.h"
#define pi 3.1415926535897932
#include <math.h>
#include <fstream>
#include <vector>

void set_boundary(Grid &u){
  unsigned int i=0,j=u.n_y;;
  //loop over y= n_y (ngp_y-1)
  for (i=0; i != (u.n_x)+1; ++i){
    //set u
    u(i,j)=sin(2*pi*i*u.h_x)*sinh(2*pi);
    //std::cout<< i*u.h_x << " " << u.n_y << " " << u(i,j) << std::endl;
  }      
}

void set_rhs(Grid &f){
 unsigned int i=0,j=0;
  //std::cout << "setting rhs" << std::endl;
  for (j=0; j != f.n_y+1; ++j)
     for(i=0; i!= f.n_x+1; ++i)
        {    
        f(i,j)=4*pi*pi*sin(2*pi*i*f.h_x)*sinh(2*pi*j*f.h_y);
        }   
}

//////////////SUM////////////////////////////////////////
// SUM(PARENT_ARG, SUM_ARG, CONST_WEIGHT)

void sum(Grid &a, Grid &b, Grid &c, const double & alpha) // performs a = b + alpha *c
{
	


	double * x, *y, *z;
	x = a.vec;
	y = b.vec;
	z = c.vec;

	for (size_t j = 0; j< (a.size); ++j)
		x[j] = y[j] + alpha * z[j];


}

////////////////////////////////////////////////////////



