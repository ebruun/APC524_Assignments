#include <stdio.h> 
#include <stdlib.h>
#include <math.h>
#include "integrator.h"


int main(int argc, char* argv[])
{
	if (argc != 3){
		printf("enter a timestep and step number you silly goose!\n");
		exit(0);
	} // so dont seg fault if forget to pass both args

	double dt = atof(argv[1]); //timestep 
	int nsteps = atoi(argv[2]); //number of timesteps

	//printf("integrate with timestep %f for %d steps \n", dt, nsteps);

	const int n=2; //dimension of system
	double t=0; // initialize time
	double x[n]; // initialize x, x(0)=q x(1)=y
	// want to keep track of prev timesteps for multi step methods

	x[0]=0; // initial conditions 
	x[1]=0;


	int Duffode(int n, double t, const double *x, double *Dx); 
	// will compute Dx of duffing ode
	FuncPtr rhs=&Duffode; // rhs =address of fn to compute duffing ode


	Integrator *ints; // pointer to integrator object
    ints = integrator_new(n, dt, rhs); // create integrator structure

	printf("%15.8f %15.8f %15.8f \n",t, x[0], x[1]);
	// print initial time and state
	for (int ii = 1; ii <= nsteps; ++ii)
	{
		t=ii*dt; // increment time
		// preform integrator step using linked scheme
		integrator_step(ints,  t,  x);
		// passes pointer x so can update
		// print output in desired format
		printf("%15.8f %15.8f %15.8f \n",t, x[0], x[1]);
	}
	integrator_free(ints); // clear integrator from heap
	return 0;
}

// passing a  pointer to a fn allows to change value
int Duffode(int n, double t, const double *x, double *Dx)
{
	double delta=0.2;
	double gamma=0.3;
	double omega=1;
	// duffing ode below
	Dx[0] = x[1];
	Dx[1] =  x[0] - pow(x[0], 3) - delta*x[1] + gamma*cos(omega*t);
	return 0;
}

