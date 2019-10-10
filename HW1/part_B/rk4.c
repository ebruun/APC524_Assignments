/*
 * rk4.c
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

  int rhserr; //error flag

  // Dimension of the problem (2 in this case)
  const int n     = integrator->n;
  const double dt = integrator->dt;

  // This is a VLA = dynamically changing size
  double Dx[n];

  //rk4 terms
  double k1[n];
  double k2[n];
  double k3[n];
  double k4[n];

  double t_rk4 = t + dt/2; //rk4 modified time
  double x_dx_temp[n]; //vector used as input to calculate new Dx depending on k,i

  /* 
  *1. Bubble up any errors from RHS function
  *2. Call ddo_ode function
  *    called with pointer to the function in the integrator struct
  *    Update the value going into the rhs function with calculated k-values (as per rk4)
  *    Update the t-value going into the rhs function (as per rk4)
  *3. Once k1,k2,k3,k4 are calculated, find new x_dx values
  */

  //x_dx(previous step) --> Dx --> k1, t = t
  rhserr = integrator->rhs(n, t, x_dx, Dx);
  if (rhserr != 0)
  {
    return rhserr;
  }

  for (int i = 0; i < n; ++i)
  {
    k1[i]         = dt*Dx[i];
    x_dx_temp[i]  = x_dx[i] + k1[i]/2;
  }


  //(x_dx + k1/2) --> Dx --> k2, t = t + dt/2
  rhserr = integrator->rhs(n, t_rk4, x_dx_temp, Dx);
  if (rhserr != 0)
  {
    return rhserr;
  }

  for (int i = 0; i < n; ++i)
  {
    k2[i]         = dt*Dx[i];
    x_dx_temp[i]  = x_dx[i] + k2[i]/2;
  }


  //(x_dx + k2/2) --> Dx --> k3, t = t + dt/2
  rhserr = integrator->rhs(n, t_rk4, x_dx_temp, Dx);
  if (rhserr != 0)
  {
    return rhserr;
  }  

  for (int i = 0; i < n; ++i)
  {
    k3[i]         = dt*Dx[i];
    x_dx_temp[i]  = x_dx[i] + k3[i];
  }


  //(x_dx + k3) --> Dx --> k4, t = t + dt
  rhserr = integrator->rhs(n, (t+dt), x_dx_temp, Dx);
  if (rhserr != 0)
  {
    return rhserr;
  }  

  for (int i = 0; i < n; ++i)
  {
    k4[i]         = dt*Dx[i];
  }


  // Update x and dx, with all the new rk4 terms
  for (int i = 0; i < n; ++i)
  {
    
    x_dx[i] += (k1[i]/6 + k2[i]/3 + k3[i]/3 + k4[i]/6);
  }
  

  /* Successful exit */
  return 0;
}