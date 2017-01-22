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
  unsigned int n_x=0, n_y=0, c=0, timesteps,vtk_spacing;
  double eps = 0.0,tau=0.0, k=0.0,alpha=0.0;
  str << argc << " "<< argv[1] << " " << argv[2] <<" " << argv[3] << " " << argv[4] << " " << argv[5] << " " << argv[6] << " " << argv[7] << " " << argv[8]<< " " << argv[9];
  str >> n_x >> n_x >> n_y >> c>>eps>>timesteps>>tau>>k>>alpha>>vtk_spacing; //n_x is wrongly written first and rewritten later
  
	//std::cout<< n_x << " " << n_y << " " << c << " " << eps << " "<< timesteps<< " "<< tau<< " " << k<< " " << alpha<< " "<< vtk_spacing << std::endl;
//default initialization   
  Grid u_old(n_x,n_y);
  Grid u_new(n_x,n_y);  

//setboundaryconditions function
  set_init_bc(u_old);
  
  //Timestep iteration starts from here
	// for (i =0 ; i <= timesteps; i =i+tau){ 
  set_rhs(u_old,tau, k, alpha);
  
//RBGS
  //solver(u_new,u_old,c,eps);
	//}
  
//write grid
  /*if (rank == 0)
  write(u, "solution.txt");*/
  
  //std::cout<<"This is processor" <<rank<<std::endl;
  u_old.release(); u_new.release(); 
  MPI_Finalize(); 
  return 0;

}












































