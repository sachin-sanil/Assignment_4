#include "grid.h"
#define pi 3.1415926535897932
#include <math.h>
#include <fstream>
#include <vector>
void write1(Grid &u,std::string output_file);
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
  //write1(u,"u_init.txt");
}

void set_rhs(Grid &f,Grid &u,double &alpha, double &tau, double &k){
unsigned int x=0,y=0;
  //double h_x=u.h_x, h_y=u.h_y;
  double n_x=u.n_x, n_y=u.n_y;
  //std::cout << "setting rhs" << std::endl;
  for(y=1;y<n_y;++y)
	  for (x=1; x < n_x ; ++x)
        {    
        f(x,y) = u(x,y) + k*(1-alpha)*tau*(((u(x+1,y) + u(x-1,y)-2*u(x,y))*n_x*n_x) + ((u(x,y+1)+u(x,y-1)-2*u(x,y))*n_y*n_y));
        }   
		//std::cout << x << " "<< y << std::endl;
		//std::cout << alpha << " "<< tau << " "<< k <<std::endl;
		//write1(f,"u_rhs.txt");
		//write1(u,"u_init.txt");
}

//////////////////////////////////////////////////////////////////////////

void write1(Grid &u,std::string output_file){
unsigned int i=0,j=0;
 //initialize output_file
 std::ofstream file(output_file);
 //loop over j ... n_y
 for (j=0; j != u.n_y+1; ++j){
    //loop over i ... n_x
    for (i=0; i != u.n_x+1; ++i)
    {    
      //print x y value
     file << i*u.h_x << " " << j*u.h_y << " " << u(i,j) << std::endl;     
    }         
 } 
 file.close();
}


////////////////////////////////////////////////////////



