#include"grid.h"
#include<iostream>


Grid::Grid(unsigned int size_x,unsigned int size_y): n_x(size_x),n_y(size_y)
{
    h_x = 1.0/n_x;
    h_y = 1.0/n_y;
    unsigned int n=(n_x+1)*(n_y+1);

    vec = new double[n];
	size = n;
    for (unsigned int i=0; i!=(n_x+1)*(n_y+1); ++i)
      vec[i]=0;
} 

double& Grid::operator ()(unsigned int x,unsigned int y){
	//std::cout << x << " " << y << " " << (y*(n_x+1) + x) << std::endl;
	return vec[(y*(n_x+1) + x)];
}
