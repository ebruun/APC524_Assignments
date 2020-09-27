#include <stdio.h> 
#include <stdlib.h>
#include "integrator.h"

struct integrator_t
{
    int n;
    double dt;
    FuncPtr rhs;
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
    return ints; //goal was to create pointer ints
}


int integrator_step(Integrator *integrator, double t, double *x)
{
	double Dx[integrator->n]; // Dx is array length n
	double k1[integrator->n]; // rk4 parameters
	double k2[integrator->n];
	double k3[integrator->n];
	double k4[integrator->n];
	double xtemp[integrator->n]; // use to compute rk4 params
	for (int ii = 0; ii < integrator->n; ++ii)
	{
		Dx[ii] = 0; //initialize
	}
	int errflag = ((integrator->rhs))(integrator->n, t, x, Dx);
	// now Dx is k1, use to calc new xtemp, pass to rhs to calc k2
	// seperate for loops for each k so works for arb n
	for (int ii = 0; ii < integrator->n; ++ii)
	{
		k1[ii]=Dx[ii];
		xtemp[ii]=x[ii]+ (integrator->dt)*k1[ii]/2;
	}
	// use xtemp to update Dx to k2
	errflag = ((integrator->rhs))(integrator->n, t+(integrator->dt)/2, xtemp, Dx);
	for (int ii = 0; ii < integrator->n; ++ii) 	// now calc k2 using xtemp
	{
		k2[ii]=Dx[ii]; //updated Dx so its now k2, now set
		xtemp[ii]=x[ii]+ (integrator->dt)*k2[ii]/2;  //pass to rhs to calc k3
	}
	errflag = ((integrator->rhs))(integrator->n, t+(integrator->dt)/2, xtemp, Dx);
	for (int ii = 0; ii < integrator->n; ++ii)
	{
		k3[ii]=Dx[ii];  //Dx is k3 after update
		xtemp[ii]=x[ii]+ (integrator->dt)*k3[ii]; //pass to rhs to calc k4
	}
	errflag = ((integrator->rhs))(integrator->n, t+integrator->dt, xtemp, Dx);
	for (int ii = 0; ii < integrator->n; ++ii)
	{
		k4[ii]=Dx[ii];  //Dx is k4
		// now have all k, can calc next step x immediatly
		 x[ii] = x[ii] + integrator->dt * (k1[ii]/6 + k2[ii]/3 + k3[ii]/3 + k4[ii]/6);
	}
	return 0;
}


void integrator_free(Integrator *integrator)
{
    free(integrator);
    integrator=0;  //the C tutorial said to do this
}