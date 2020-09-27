#ifndef INTEGRATOR_H
#define INTEGRATOR_H

/*
 * Pointer to a function that computes the time derivative of x
 * Dx = f(x,t) ,
 * where n is the dimension of the solution x.
 *
 * Should return 0 if successful, nonzero if error.
 */
typedef int (*FuncPtr)(int n, double t, const double *x, double *Dx);

typedef struct integrator_t Integrator;

/*
 
 * Return a new Integrator object
 *
 * n:   dimension of state vector 
 * dt:  timestep
 * rhs: pointer to a function to compute right-hand side of
 *      Dx = f(x,t).
 
 Can't change the value (const)
 
 Typedef - defining a new type of variable
 
 why make a pointer to a function, we want to provide a pointer to a function to another function
 
 */
Integrator *integrator_new(int n, double dt, FuncPtr rhs);
/* we want to pass just memory address, so no *before FuncPtr
 spits out a pointer to an integrator
 */

/* Free the Integrator object and any associated memory */
void integrator_free(Integrator *integrator);

/*
 * Advance x by one timestep
 *
 * Should return 0 if successful, nonzero if error.
*/
int integrator_step(Integrator *integrator, double t, double *x);


#endif /* INTEGRATOR_H */
