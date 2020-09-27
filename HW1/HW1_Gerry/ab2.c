#include <stdio.h> 
#include <stdlib.h>
#include "integrator.h"

struct integrator_t
{
    int n;
    double dt;
    FuncPtr rhs;
    int Eulerflag; // use to check if first step or not
    double *prevDx; //use to store Dx from previous step
}; // define structure type Integrator
// a lot more comments are in euler.c
// most of this is repeated

Integrator *integrator_new(int n, double dt, FuncPtr rhs)
{
    Integrator* ints = malloc(sizeof(Integrator));
    // pointer to integrator structure
   	ints->n = n; // same as (*ints).n=n, pointing to ints.n
    ints->dt = dt;
    ints->rhs = rhs;
    ints->Eulerflag =0;
    ints->prevDx=malloc(n*sizeof(double)); //space for appropriate sized array
    return ints; //goal was to create pointer ints
}


int integrator_step(Integrator *integrator, double t, double *x)
{
	double Dx[integrator->n]; // Dx is array length n
	//double Dpast[integrator->n]; // use to store Dx at pref timestep
	for (int ii = 0; ii < integrator->n; ++ii)
	{
		Dx[ii] = 0; //initialize Dx or get segfaults
	//	Dpast[ii]=0;
	}
	// calculate Dx for timestep
	int errflag = ((integrator->rhs))(integrator->n, t, x, Dx);

	if (integrator->Eulerflag==0) //first timestep do explicit euler
	{
		for (int ii = 0; ii < integrator->n; ++ii)
		{
			x[ii] = x[ii] + integrator->dt * Dx[ii]; //Explicit Euler step
			integrator->prevDx[ii] = Dx[ii]; // save val of Dx for next step
		}
		integrator->Eulerflag=1; //dont do exp euler again
		//printf(" set eulerflag to 1, this should only happen once\n");
		return 0; // exit fn early since already updated x
	}

	for (int ii = 0; ii < integrator->n; ++ii)
	{	// use prevDx saved in past step to integrator to compute ab2 step
		x[ii] = x[ii] +  3*integrator->dt * Dx[ii]/2 - integrator->dt * integrator->prevDx[ii]/2; 
		integrator->prevDx[ii]= Dx[ii]; // save val of Dx for next step
	}
	
	return 0;
}


void integrator_free(Integrator *integrator)
{
    free(integrator);
    integrator=0;  //the C tutorial said to do this
}