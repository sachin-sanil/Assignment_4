#include <vector>
#include <iostream>
#include "Timer.h"
#include <sstream>
#include "grid.h"
#include <math.h>
#include <omp.h>
#include <assert.h>
#include<algorithm>
#include <mpi.h>
#include <fstream>
#define pi 3.1415926535897932

double dot(Grid &u, Grid &d);
//double norm_squared(Grid &u);
void mult_A(Grid &u, Grid &v,double &alpha, double &tau, double &k) ;
void sum(Grid& u, double alpha, Grid& d);
void sum(Grid&d, Grid& r, double beta);
void write(Grid &u,std::string output_file);

int coord[2];
unsigned int init;
unsigned int fin;
int cart_rank;
int ndims=2 ,dim[2],periods[2]={0,0} ,reorder=0;
int q,rem;
MPI_Comm comm_cart;
extern int num_proc;
MPI_Request reqs[4];
MPI_Status status[4];

int u_r, d_r;
int down_send, up_send,down_rec, up_rec;

//c here is maximum number of cg iterations.

void solver(Grid &u,Grid &f,int c,double eps,double &alpha, double &tau, double &k){
	
	
	dim[0] = num_proc;
	dim[1] = 1;
	
	//std::cout << "something"<< std::endl;
	//std::cout << ndims << " " << dim[0] << " " << dim[1] << " " << periods[0] << " " << periods[1] << " " << reorder << std::endl;
	MPI_Cart_create(MPI_COMM_WORLD, ndims, dim, periods, reorder, &comm_cart);
	MPI_Comm_rank (comm_cart , &cart_rank);
	MPI_Cart_shift(comm_cart, 0, 1, &d_r,&u_r);
	
	//std::cout<< cart_rank << " "  << u_r << " " << d_r << std::endl;
	
	q = (u.n_y-1)/num_proc;	
	rem = ((int)u.n_y-1)%num_proc;
	
	int elem_proc, disp_proc;
	
	int recvcounts[num_proc], disps[num_proc];
	
	MPI_Cart_coords(comm_cart, cart_rank, 2, coord);
//	std::cout<<coord[0]<<" "<<cart_rank<<std::endl;
	
		 if (coord[0]!=num_proc-1)
			{
				init = 1+coord[0]*q;
				fin = (coord[0]+1)*q;
			}
		else
			{
				init = 1+coord[0]*q;
				fin = (coord[0]+1)*q+rem;					
			}
		//std:: cout << coord[0] << " " << coord[1] << " " << init << " " << fin << std::endl;
	elem_proc=(fin-init+1)*(u.n_x+1);
	disp_proc=init*(u.n_x+1);
	//std::cout << init << " " << fin << " " << elem_proc << " for rank " << cart_rank << std::endl;
	
	unsigned int n_x = u.n_x;
	unsigned int n_y = u.n_y;
	//int ngp_x = n_x+1; //number of grid points
	//int ngp_y = n_y+1;

	//Grid res(n_x,n_y);
	double numer = ((double)(n_x*n_x)*2 + (double)(n_y*n_y)*2);  
	//double denom = 1.0/numer;           // 1/(2/h_x*h_x  + 2/h_y*h_y + k*k)
	siwir::Timer timer;
	double time = 100.0;
	unsigned int x=1, y=1;
	double _nx2 = (double)n_x*(double)n_x; 
	double _ny2 = (double)n_y*(double)n_y;

	Grid r(n_x,n_y), z(n_x,n_y),d(n_x,n_y); // Grids initialized Just check if these are resized properly
	double delta_0=0, delta_1=0, alpha1=0, beta=0; // Variables used in CG
	double res = 0; // norm of residual ||r||2
	double dtz =0 ; //stores inner product d'z
	//only solver to be timed
	timer.reset();
	double temp = 0; 
	//CG Algorithm starts from here
	for(y=1; y <= n_y-1 ; ++y){
	for(x=1; x <= n_x-1 ; ++x){
			temp = f(x,y)+ ((alpha*k*tau)*((_nx2 * (u(x-1,y)+u(x+1,y))) +  (_ny2 * (u(x,y-1)+u(x,y+1))) - (numer*u(x,y))))- u(x,y); // r = f-Au
			r(x,y) = temp;
			d(x,y) = temp; //d = r	
	}
	}

	delta_0 = dot(r,r);
	res = delta_0; 
	
	if (sqrt(res) > eps) //exit criterion check if residual is less than some epsilon
	{
	for(int i=0; i<c; ++i){
	MPI_Barrier(comm_cart);
	mult_A(z,d,alpha,tau,k);MPI_Barrier(comm_cart);
	dtz = dot(d,z); MPI_Barrier(comm_cart);
	alpha1 = delta_0/dtz; MPI_Barrier(comm_cart);
	sum(u,alpha1,d); MPI_Barrier(comm_cart);
	sum(r,-alpha1,z); MPI_Barrier(comm_cart);
	delta_1 = dot(r,r); MPI_Barrier(comm_cart);
	res = delta_1 ; MPI_Barrier(comm_cart);
	
	if (sqrt(res) <= eps)
	  break;
	else
		beta = delta_1/delta_0;
		sum(d,r,beta);
		delta_0=delta_1;
	}
	} 
	//CG Algorithm ends here
	
	
	//gather u for all processors
	MPI_Barrier(comm_cart);
	MPI_Allgather(&elem_proc, 1, MPI_INT, &(recvcounts), 1, MPI_INT, comm_cart);
	MPI_Allgather(&disp_proc, 1, MPI_INT, &(disps), 1, MPI_INT, comm_cart);
/* 	if (cart_rank==0){
	//print disp_proc
		for (auto i=0; i!=num_proc; i++)
			std::cout<<disps[i]<<" "<<recvcounts[i] << " ";
		std::cout << std::endl;
	} */

	
	time = std::min(time, timer.elapsed());
	
	//std::cout << cart_rank << " " << disp_proc+elem_proc-1 << " " << u.vec[0] << std::endl;
	MPI_Barrier(comm_cart);
	MPI_Allgatherv(&(u.vec[disp_proc]),elem_proc, MPI_DOUBLE, &(u.vec[0]), recvcounts, disps, MPI_DOUBLE, comm_cart);

if (cart_rank==0)
	std::cout << "Time Taken " << time << std::endl;
	
if (cart_rank==0)
	write(u, "solution1.txt");
 
if(cart_rank ==0)
		{
		std::ofstream bout;
		bout.open("solution.txt");
		
		bout<<"# x y u(x,y)"<<std::endl;
		int j=0;
					for (unsigned int i =0;i<=u.n_x;i++)
					{
						bout<<i*u.h_x<<" "<<j*u.h_y<<" "<<u(i,j)<<std::endl;
					}
					bout<<std::endl;
				bout.close();
		
		}
		for(int p =0;p<num_proc;p++)
		{
			MPI_Barrier(comm_cart);
			if(coord[0] ==p)
			{
				//cout<<"Writing the results to solution"<<rank<<".txt :"<<endl;
				std::ofstream fout;
				fout.open("solution.txt",std::ios_base::app);
				
				//fout<<"# x y u(x,y)"<<endl;
				for(unsigned int j = init;j<=fin;j++)
				{
					for (unsigned int i =0;i<=u.n_x;i++)
					{
						fout<<i*u.h_x<<" "<<j*u.h_y<<" "<<u(i,j)<<std::endl;
					}
					fout<<std::endl;
				}
				fout.close();
			}
		}
		MPI_Barrier(comm_cart);
		if(cart_rank ==0)
		{
			std::ofstream cout;
			cout.open("solution.txt",std::ios_base::app);
			int j=u.n_y;
						for (unsigned int i =0;i<=u.n_x;i++)
						{
							cout<<i*u.h_x<<" "<<j*u.h_y<<" "<<u(i,j)<<std::endl;
						}
						cout<<std::endl;
					cout.close();
		
		}

	r.release(); z.release(); d.release();
}

double dot(Grid &u, Grid &v) {
	assert(u.size == v.size);
	double part_res = 0, res=0;
	for (unsigned int i = init*(u.n_x +1); i<= (fin * (u.n_x+1)) + u.n_x; i++) {
		part_res += u.vec[i] * v.vec[i];
	}
	MPI_Allreduce(&part_res, &res, 1, MPI_DOUBLE, MPI_SUM,comm_cart);
	return res;
}



void mult_A(Grid &u, Grid &v,double &alpha, double &tau, double &k) { //does. u = A*d

	//Functioning variables
	assert(u.size == v.size);
	double numer = ((double)(u.n_x*u.n_x)*2 + (double)(u.n_y*u.n_y) * 2) ;
	//double denom = 1.0 / numer;
	double d_nx2 = u.n_x*u.n_x*0.25;
	double d_ny2 = u.n_y*u.n_y;
	unsigned int n_x=u.n_x;
	double aux1=0.0, aux2=0.0, aux3=0.0,aux4=0.0;
	
	//std::cout<<&(v.vec[up_send])<<"error adress"<<std::endl;
	
	down_send = (u.n_x+1)*init + 1;
	up_send = (u.n_x+1)*fin+1;
	down_rec = (u.n_x+1)*(init-1)+1;
	up_rec = (u.n_x+1)*(fin+1)+1;

	MPI_Isend( &(v.vec[up_send]), (u.n_x-1), MPI_DOUBLE, u_r, 0, comm_cart, &reqs[0]);
	MPI_Isend( &(v.vec[down_send]), (u.n_x-1), MPI_DOUBLE, d_r, 1, comm_cart, &reqs[1]);
	MPI_Irecv( &(v.vec[down_rec]), (u.n_x-1), MPI_DOUBLE, d_r, 0, comm_cart, &reqs[2] );
	MPI_Irecv( &(v.vec[up_rec]), (u.n_x-1), MPI_DOUBLE, u_r, 1, comm_cart, &reqs[3] ); 
	
	MPI_Waitall(4, reqs, status);
	
	//only inner grid points are touched
	for (unsigned int y = init; y <= fin; y = y + 1)
		for (unsigned int x = 1; x < v.n_x; x = x + 1) {
			aux1 = v.vec[(y*(n_x+1) + x)];
			aux2 = (numer*v.vec[(y*(n_x+1) + x)]);
			aux3 = (d_nx2 * (v.vec[(y*(n_x+1)+x-1)] + v.vec[(y*(n_x+1) + x+1)]));
			aux4 = (d_ny2 * (v.vec[((y-1)*(n_x+1) + x)] + v.vec[((y+1)*(n_x+1) + x)]));
			
			u.vec[(y*(n_x+1) + x)] = aux1 -(alpha*k*tau)* (aux3+aux4-aux2);
			
			
			//u.vec[(y*(n_x+1) + x)] = v.vec[(y*(n_x+1) + x)]-((alpha*k*tau)((-numer*v.vec[(y*(n_x+1) + x)]) + (d_nx2 * (v.vec[(y*(n_x+1)+x-1)] + v.vec[(y*(n_x+1) + x+1)])) + (d_ny2 * (v.vec[((y-1)*(n_x+1) + x)] + v.vec[((y+1)*(n_x+1) + x)]))));
		}

}

void sum(Grid& u, double alpha, Grid& d) {
	for (unsigned int i = init*(u.n_x +1); i<= (fin * (u.n_x+1)) + u.n_x; i++)
		u.vec[i] += alpha*d.vec[i];
}

void sum(Grid&d, Grid& r, double beta) {
	for (unsigned int i = init*(d.n_x +1); i<= (fin * (d.n_x+1)) + d.n_x; i++)
		d.vec[i] = r.vec[i] + beta*d.vec[i];
}

