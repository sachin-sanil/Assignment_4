#include "grid.h"
#define pi 3.1415926535897932
#include <math.h>
#include <fstream>
#include <vector>


void set_init_bc(Grid &u){
  unsigned int x=0,y=0;
  double h_x=u.h_x, h_y=u.h_y;
  double n_x=u.n_x, n_y=u.n_y;


  //loop over y= n_y (ngp_y-1)
	for(y=0;y!=n_y+1;++y)
	  for (x=0; x != n_x+1; ++x)
	{
    //set u
    u(x,y)=sin(pi*x*h_x)*sin(pi*y*h_y);
    //std::cout<< i*u.h_x << " " << u.n_y << " " << u(i,j) << std::endl;
  }      
}

void set_rhs(Grid &u_old, double tau,double k, double alpha){
 unsigned int i=0,j=0;
 double n_x=u_old.n_x, n_y=u_old.n_y;
 double h_x=u_old.h_x, h_y=u_old.h_y;

  //std::cout << "setting rhs" << std::endl;
  for (j=1; j != n_y; ++j)
     for(i=1; i!= n_x; ++i)
        {    
        u_old(i,j) += tau*(1-alpha)*k*((u_old(i-1,j) + u_old(i+1,j))/(h_x*h_x) + (u_old(i,j-1) + u_old(i,j+1))/(h_y*h_y) - 2*u_old(i,j)/(h_x*h_x) - 2*u_old(i,j)/(h_y*h_y));
        }   

	for (j=0; j != n_y+1; j=j+n_y)
     for(i=0; i!= n_x+1; ++i)
		{
			u_old(i,j) += tau*(1-alpha)*k*((- 2*u_old(i,j)/(h_x*h_x) - 2*u_old(i,j)/(h_y*h_y)));
		}
	
	// update the lateral boundary
	for(j=1; j!=n_y;++j)
		for(i=0;i!=n_x+1; i=i+n_x)
		{
		 u_old(i,j) += tau*(1-alpha)*k*((- 2*u_old(i,j)/(h_x*h_x) - 2*u_old(i,j)/(h_y*h_y)));
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



