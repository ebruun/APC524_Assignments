/* 
 * Part B
 * ddo.c
 */

#include <stdlib.h> //atof function doesn't work without it
#include <stdio.h>
#include <math.h>
#include "integrator.h"

//Declaring a global variable as per Gabe's correction to assignment
float w;

int main(int argc, char* argv[])
{

  //USER INPUT
   if( argc == 3 ) {
      printf("The Driving Frequency: %s\n", argv[1]);
      printf("The Number of Timesteps: %s\n", argv[2]);
   }
   else if( argc > 3 ) {
      printf("Too many arguments specified!!!\n");
      exit(0);
   }
   else {
      printf("Not enough arguments specified!!!\n");
      exit(0);
   }

   //VARIABLE DEFINITION
  w    = atof(argv[1]);  //frequency
  int n_step = atoi(argv[2]);  //number of timesteps

  double t  = 0; //Initialize time
  double dt = 18.84955592153876/n_step; //Size of timestep based on T = 6*pi
  int    n  = 2; //Problem Size

  /*
  *Initialize the x[0] & dx[0] variables to 0. Recall form:
  *   x[1]   = x[0]  + dt*dx[0] = x[0]  + dt*Dx[0]
  *   dx[1]  = dx[0] + dt*f(t,x[0],dx[0]) = dx[0]  + dt*Dx[1]
  */
  double x_dx[2] = {0.0,0.0};

  
  //RUNNING PROGRAM

  /*
  *1. Initialize the function that computes Dx = slope of x and dx
  *2. This format is declared in the header file
  *3. Function defined after main()
  */
  int ddo_ode(int n, double t, const double *x_dx, double *Dx);
  FuncPtr rhs=&ddo_ode; //the address of function to calculate slopes


  /*
  *1. Create a pointer to a structure of type Integrator
  *    Integrator is defined as 'struct integrator_t'
  *    The typedef occurs in the header file
  *    'struct integrator_t' is defined in .c file for integration method
  *2. The pointer variable created here is called integrator
  *    confusing naming, change if I have time
  *3. Create new Integrator structure, and then point to it
  */
  Integrator *integrator;
  integrator = integrator_new(n,dt,rhs);


  //Print initial conditions
  printf("%15.8f %15.8f %15.8f\n", t, x_dx[0], x_dx[1]);

  //Enter solution loop
  for (int i = 1; i <= n_step; ++i)
  {
    /*
    *Step forward based on the linked integrator
    * euler, rk4, ab2
    * Pointer to x and dx updated in previous step
    */
    integrator_step(integrator, t, x_dx);

    // Increment Time
    t += dt;

    // print output in desired format
    printf("%15.8f %15.8f %15.8f \n",t, x_dx[0], x_dx[1]);


  }


  //Free the memory assigned to integrator
  integrator_free(integrator);
  
  return 0;
}


// Function that computes Dx (slope of x and dx)
int ddo_ode(int n, double t, const double *x_dx, double *Dx)
{
  //Integration function constants (provided)
  double beta=0.25;
  double omega_0=1;
  

  /*
  * 1. Update the right-hand side for integration scheme
  *     Dx[0] = rhs for x
  *     Dx[1] = rhs for dx
  * 2. This is not an efficient way to do this
  *     Using the 'n' variable so the compiler doesn't throw an error
  *     Could just do the Dx[i] equations explicitly since we know there are 2 states
  */

  for (int i=0; i<n; i++)
  {
    if (i==0)
    {
      Dx[i] = x_dx[1];
    }
    else
    {
      Dx[i] = (cos(w*t) - omega_0*x_dx[0] - 2*beta*x_dx[1]);
    }
  }

  //Could do this as well...
  //Dx[0] = x_dx[1];
  //Dx[1] = (cos(w*t) - omega_0*x_dx[0] - 2*beta*x_dx[1]);
  return 0;
}