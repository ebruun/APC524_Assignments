/*
 * euler.c
 */

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "integrator.h"


/* Integrator object for Forward Euler method */
struct integrator_t
{
  int n;		   	/* dimension of state vector */
  double dt;		  	/* time step */
  FuncPtr rhs;		      	/* right-hand-side of \dot x = f(x,t) */
};
  
/* "Constructor" */
Integrator *integrator_new(const int n, const double dt, FuncPtr rhs)
{
  Integrator *integrator = (Integrator *) malloc(sizeof(*integrator));
  assert(integrator);
  
  integrator->n = n;
  integrator->dt = dt;
  integrator->rhs = rhs;

  return integrator;
}

/* "Destructor" */
void integrator_free(Integrator *integrator)
{
  assert(integrator);
  free(integrator);
}

/* Stepper */
int integrator_step(Integrator *integrator, double t, double *x_dx)
{
  assert(integrator);


  // Dimension of the problem (2 in this case)
  const int n     = integrator->n;
  const double dt = integrator->dt;

  // This is a VLA = dynamically changing size
  double Dx[n];		
  
  /* 
  *1. Bubble up any errors from RHS function
  *2. Call ddo_ode function
  *    called with pointer to the function in the integrator struct
  *3. Goal is to update the value dx = f(x,t), slope of function
  */
  int rhserr = integrator->rhs(n, t, x_dx, Dx);
  if (rhserr != 0)
  {
    return rhserr;
  }

  /*
  * 1. Forward Euler: Update the values of x and dx
  *     x[1]   = x[0]  + dt*dx[0] = x[0]  + dt*Dx[0]
  *     dx[1]  = dx[0] + dt*f(t,x[0],dx[0]) = dx[0]  + dt*Dx[1]
  * 2. Scale values by timestep size
  *      x is related to the first term of Dx
  *      dx is related to the second term in Dx
  * 3. The problem is of n = 2 dimensions
  */
  for (int i = 0; i < n; ++i)
  {
    x_dx[i] += dt * Dx[i];
  }
  
  /* Successful exit */
  return 0;
}




