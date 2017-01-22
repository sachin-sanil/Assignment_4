#include <vector>
#include <iostream>
#include "Timer.h"
#include <sstream>
#include "grid.h"
#include <math.h>
#include "misc.h"
#include <assert.h>
#include <mpi.h>


#define pi 3.1415926535897932

int rank, num_proc;
  
void solver(Grid &, Grid &, int, double);

int main(int argc, char *argv[]){
  
  
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  std::stringstream str;  
  unsigned int n_x=0, n_y=0, c=0;
  double eps = 0.0;
  str << argc << " "<< argv[1] << " " << argv[2] <<" " << argv[3] << " " << argv[4];
  str >> n_x >> n_x >> n_y >> c>>eps; //n_x is wrongly written first and rewritten later
  
//default initialization   
  Grid u(n_x,n_y);
  Grid f(n_x,n_y);  

//setboundaryconditions function
  set_boundary(u); 
  set_rhs(f);
  
//RBGS
  solver(u,f,c,eps);
  
//write grid
  /*if (rank == 0)
  write(u, "solution.txt");*/
  
  //std::cout<<"This is processor" <<rank<<std::endl;
  u.release(); f.release(); 
  MPI_Finalize(); 
  return 0;

}












































