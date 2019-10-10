/*
 * ab2.c
 */

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "integrator.h"


/* Integrator object for Forward Euler method */
struct integrator_t
{
  int n;		   	    /* dimension of state vector */
  double dt;		  	/* time step */
  double *Dx_prev;  // NEW: Pointer to array storing Dx from previous step
  FuncPtr rhs;	    /* right-hand-side of \dot x = f(x,t) */
};
  
/* "Constructor" */
Integrator *integrator_new(const int n, const double dt, FuncPtr rhs)
{
  Integrator *integrator = (Integrator *) malloc(sizeof(*integrator));
  assert(integrator);
  
  integrator->n = n;
  integrator->dt = dt;
  integrator->rhs = rhs;

  //Initialize and point to memory for array
  integrator->Dx_prev = malloc(n*sizeof(double)); 

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

  //Pointer to the 2-element array, update throughout integration
  // (Just to simplify syntax later)
  double *Dx_prev = integrator->Dx_prev; 

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


  // First time step do a forward Euler
  if (t == 0)
  {
    for (int i = 0; i < n; ++i)
    {
      x_dx[i]   += dt * Dx[i];
      Dx_prev[i] = Dx[i];
    } 
  }
  // For all other steps use Dx from previous timestep
  else
  {
    for (int i = 0; i < n; ++i)
    {
      x_dx[i]   += dt * (3.0*Dx[i]/2.0 - Dx_prev[i]/2.0);
      Dx_prev[i] = Dx[i];
    }
  }
  
  /* Successful exit */
  return 0;
}