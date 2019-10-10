/*
 * lessnasty.c
 *
 * To compile this on adroit, give gcc the -lm flag
 * gcc -std=c11 -o lessnasty.out lessnasty.c -lm
 */

  /*Note: I have re-written Part A to more closely match my final structure for Part B */

#include <stdlib.h> //atof function doesn't work without it
#include <stdio.h>
#include <math.h>

// Declaring global variables
// Not good practice, but this allows Euler function to conform closer to RHS function from Part B
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

  //VARIABLE INITIALIZATION & DEFINITION
  w           = atof(argv[1]);  //frequency
  int n_step  = atoi(argv[2]);  //number of timesteps

  double t  = 0;                        //Initialize time
  double dt = 18.84955592153876/n_step; //Size of timestep based on T = 6*pi
  int    n  = 2;                        //Problem Size

  /*
  *Initialize the x[0] & dx[0] variables to 0. Recall form:
  *   x[1]   = x[0]  + dt*dx[0] = x[0]  + dt*Dx[0]
  *   dx[1]  = dx[0] + dt*f(t,x[0],dx[0]) = dx[0]  + dt*Dx[1]
  */
  double x_dx[2] = {0.0,0.0};
  double Dx[2]   = {0.0,0.0};

  void euler(int n, double t, const double *x_dx, double *Dx);


  //RUNNING PROGRAM
  //Print initial conditions
  printf("%15.8f %15.8f %15.8f\n", t, x_dx[0], x_dx[1]);

  //Enter solution loop
  for (int i = 1; i <= n_step; ++i)
  {
    /*
    *Step forward based on the Euler scheme
    *Values updated in 'euler' function, with pointers
    */
    euler(n, t, x_dx, Dx);

    for (int i = 0; i < n; ++i)
    {
      x_dx[i] += dt * Dx[i];
    }    

    //Print output in desired format
    printf("%15.8f %15.8f %15.8f \n",t, x_dx[0], x_dx[1]);

    // Increment Time
    t += dt;
  }
  
  return 0;
}


//F3: Euler forward integration function
void euler(int n, double t, const double *x_dx, double *Dx)
{
  int i = 0;

  //Integration function constants (provided)
  double beta=0.25;
  double omega_0=1;

  /*
  * 1. Update the right-hand side for integration scheme
  *     Dx[0] = rhs for x
  *     Dx[1] = rhs for dx
  * 2. This is not an efficient way to do this
  *     Could just do the Dx[i] equations explicitly since we know there are 2 states
  *     This keeps it "flexible", for situations where n > 2
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
  }
