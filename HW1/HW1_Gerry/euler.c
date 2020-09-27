#include <stdio.h> 
#include <stdlib.h>
#include "integrator.h"

struct integrator_t
{
    int n;
    double dt;
    FuncPtr rhs;
}; // define structure type Integrator


Integrator *integrator_new(int n, double dt, FuncPtr rhs)
{
    Integrator* ints = malloc(sizeof(Integrator));
    // ints is a pointer to integrator
    // need to allocate it memory in heap for Integrator struct type
   	ints->n = n; // same as (*ints).n=n, pointing to ints.n
    ints->dt = dt;
    ints->rhs = rhs;
    // ints is a pointer to a struct, it has a field rhs
    // that field rhs is a pointer to a function->C is silly

    return ints; //goal was to create pointer ints
}


int integrator_step(Integrator *integrator, double t, double *x)
{
	double Dx[integrator->n]; // Dx is array length n
	for (int ii = 0; ii < integrator->n; ++ii)
	{
		Dx[ii] = 0; //initialize Dx to be safe
	}
	// need to dereference integrator.n
	//ok next step is a doozie: update Dx via duffing ode
	int errflag = (*(integrator->rhs))(integrator->n, t, x, Dx);
	// first need to dereference rhs from integrator pointer
	// rhs is pointer to duffingode, so derefrence it
	// now have fn, need to pass integrator.n t x Dx
	// x not *x cuz dont want to dereference and pass value
	// function outputs zero if sucessfull, so save this as errorflag
	
	// now have updated Dx so use to calculate Euler step
	// loop so can handle arbitrary dim n
	for (int ii = 0; ii < integrator->n; ++ii)
	{
		x[ii] = x[ii] + integrator->dt * Dx[ii]; //Explicit Euler step
	}
	return 0;
}


void integrator_free(Integrator *integrator)
{
    free(integrator); // frees up memory in heap that was allocated
    integrator=0;  //the C tutorial said to do this
}