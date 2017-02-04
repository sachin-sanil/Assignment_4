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
  
void solver(Grid &u,Grid &f,int c,double eps,double &alpha, double &tau, double &k);
void write(Grid &u,std::string output_file);

int main(int argc, char *argv[]){
  
  
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  std::stringstream str;  
  
  unsigned int n_x=0, n_y=0, c=0,vtk_spacing=0;
  double eps = 0.0,tau=0.0, k=0.0,alpha=0.0,timesteps=0.0;
  str << argc << " "<< argv[1] << " " << argv[2] <<" " << argv[3] << " " << argv[4] << " " << argv[5] << " " << argv[6] << " " << argv[7] << " " << argv[8]<< " " << argv[9];
  str >> n_x >> n_x >> n_y >> c>>eps>>timesteps>>tau>>k>>alpha>>vtk_spacing; //n_x is wrongly written first and rewritten later
  
	//std::cout<< "n_x" << " " << "n_y" << " " << "c" << " " << "eps" << " "<< "timesteps"<< " "<< "tau"<< " " << "k"<< " " << "alpha"<< " "<< "vtk_spacing" << std::endl;
	//std::cout<< n_x << " " << n_y << " " << c << " " << eps << " "<< timesteps<< " "<< tau<< " " << k<< " " << alpha<< " "<< vtk_spacing << std::endl;
	
//default initialization   
  Grid u_old(n_x,n_y);
  Grid u_new(n_x,n_y);
  Grid f(n_x,n_y);  

//setboundaryconditions function
  set_init_bc(u_old); 
  
  for (double i =0.0; i<=timesteps; i = i + tau){
  set_rhs(f,u_old,alpha,tau,k);
  
//cg solver
  solver(u_new,f,c,eps,alpha,tau,k);
	for(unsigned int l =0; l<u_new.size;++l)
		u_old.vec[l]=u_new.vec[l];
  }
 
//write grid
  if (rank == 0)
	write(u_old, "solution1.txt");

  
  //std::cout<<"This is processor" <<rank<<std::endl;
  u_old.release();u_new.release(); f.release(); 
  MPI_Finalize(); 
  return 0;

}

void write(Grid &u,std::string output_file){
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











































