//  file: orbital_path_solver.cpp
//
//  C++ Program to calculate the orbital path of an object given initial conditions
//   
//
//  Programmer:  Dick Furnstahl  furnstahl.1@osu.edu
//							 Cole Howell 	cwhowell42@tntech.edu
//
//  Revision history:
//      4/25/2022:	Original version, adapted from ode_test.cpp by Dick Furnstahl
//      
//      
//
//  Notes:  
//   * Example taken from the GNU Scientific Library Reference Manual
//      Edition 1.1, for GSL Version 1.1 9 January 2002
//      URL: gsl/ref/gsl-ref_23.html#SEC364
//   * Compile and link with:
//       g++ -Wall -o orbital_path_solver orbital_path_solver.cpp -lgsl -lgslcblas
//   * gsl routines have built-in 
//       extern "C" {
//          <header stuff>
//       }
//      so they can be called from C++ programs without modification
//
//********************************************************************
//The following section is leftover from the original program, the equation this program solves is the 2D gravitational force equations
// The following details are taken from the GSL documentation
// 
// The following program solves the second-order nonlinear 
//  Van der Pol oscillator equation (see background notes),
//
//     x"(t) + \mu x'(t) (x(t)^2 - 1) + x(t) = 0
//
// This can be converted into a first order system suitable for 
//  use with the library by introducing a separate variable for 
//  the velocity, v = x'(t).  We assign x --> y[0] and v --> y[1].
//  So the equations are:
// x' = v                  ==>  dy[0]/dt = f[0] = y[1]
// v' = -x + \mu v (1-x^2) ==>  dy[1]/dt = f[1] = -y[0] + mu*y[1]*(1-y[0]*y[0])
//
//*********************************************************************

// include files
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>		// C++ stringstream class (can omit iostream)
#include <cmath>
using namespace std;

#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_const_mksa.h>
#include "GnuplotPipe.h"

struct Params{

	double mass;		//mass of large object in system
	double G;				//universal gravitation constant

};

// function prototypes 
int rhs (double t, const double y[], double f[], void *params_ptr);
//int rhsy (double t, const double y[], double f[], void *params_ptr);
//int xjacobian (double t, const double y[], double *dfdy,
	//      double dfdt[], void *params_ptr);
//int yjacobian (double t, const double y[], double *dfdy,
	//      double dfdt[], void *params_ptr);

//*************************** main program ****************************

int
main (void)
{
  const int dimension = 4;	// number of differential equations 

  const double eps_abs = 1.e-8;	        // absolute error requested 
  const double eps_rel = 1.e-10;	// relative error requested 

  // Define the type of routine for making steps: 
  const gsl_odeiv_step_type *type_ptr = gsl_odeiv_step_rk4;
  // some other possibilities (see GSL manual):          
  //   = gsl_odeiv_step_rk4;
  //   = gsl_odeiv_step_rkck;
  //   = gsl_odeiv_step_rk8pd;
  //   = gsl_odeiv_step_rk4imp;
  //   = gsl_odeiv_step_bsimp;  
  //   = gsl_odeiv_step_gear1;
  //   = gsl_odeiv_step_gear2;

  // Allocate/initialize the stepper, the control function, and the
  //  evolution function.
  gsl_odeiv_step *step_ptr = gsl_odeiv_step_alloc (type_ptr, dimension);
  gsl_odeiv_control *control_ptr = gsl_odeiv_control_y_new (eps_abs, eps_rel);
  gsl_odeiv_evolve *evolve_ptr = gsl_odeiv_evolve_alloc (dimension);

  gsl_odeiv_system my_system;	// structure with the rhs function, etc. 
	//gsl_odeiv_system my_y_system;	// structure with the rhsy function, etc.

  Params p; //declare a params sruct

	p.mass = GSL_CONST_MKSA_SOLAR_MASS;		//input the solar mass into the struct
	p.G = GSL_CONST_MKSA_GRAVITATIONAL_CONSTANT;	//input the Gravitational constant 

  // Load values into the my_system structure 
  my_system.function = rhs;	// the right-hand-side functions dy[i]/dt 
  my_system.jacobian = NULL;	// the Jacobian df[i]/dy[j] 
  my_system.dimension = dimension;	// number of diffeq's 
  my_system.params = &p;	// parameters to pass to rhs and jacobian 

/*
	my_y_system.function = rhsy;	// the right-hand-side functions dy[i]/dt 
  my_y_system.jacobian = yjacobian;	// the Jacobian df[i]/dy[j] 
  my_y_system.dimension = dimension;	// number of diffeq's 
  my_y_system.params = &p;	// parameters to pass to rhs and jacobian 
*/
  double tmin = 0.;		// starting t value 
  double tmax = 3.154e9;		// final t value 
  double delta_t = 8.64e4;        // step size in time
  double t = tmin;		// initialize t 

	int plot_skip = 10;	//plot every plot_skip pointsw
	int plot_delay = 10;	//wait plot_delay msec between points

	GnuplotPipe myPipe;
	myPipe.set_title ("Orbital Path of Object Around the Sun");
	myPipe.set_xlabel("x (m)");
	myPipe.set_ylabel("y (m)");
	myPipe.set_delay(1000*plot_delay);

	myPipe.init();

  double y[4];			// current solution vector for all variables
  y[0] = 1.52e11;			// initial x value, aphelion of the earth

	cout << "Enter the initial x velocity (m/s): ";
	cin >> y[1];

  //y[1] = -7200;			// initial vx value
	
	//double y[2];
	y[2] = 0;	//initial y value

	cout << "Enter the initial y velocity (m/s): ";
	cin >> y[3];

	//y[3] = 1500; //initial vy value

	

  // Set up a file names with the initial values
  ostringstream my_stringstream;	// declare a stringstream object
  my_stringstream << "orbital_path_solver" << "_vx0_" << setprecision(2) << y[1]
                  << "_vy0_" << setprecision(2) << y[3] << ".dat";

  ofstream ode_out;		// now open a stream to a file for output
  // use .str() to convert to a string, then .c_str() to convert to a char *  
  ode_out.open(my_stringstream.str().c_str());	

  // print initial values
  ode_out << "# Running ode_test with x0 = " << setprecision(2) << y[0]
          << ", vx0 = " << setprecision(2) << y[1] << ", y0 = " << setprecision(2) << y[2] << ", and vy0 = " << setprecision(2) << y[3] << endl;
  ode_out << scientific << setprecision (5) << setw (12) << t << " " 
          << setw (12) << y[0] << " " << setw (12) << y[1] << setw(12) << y[2] << setw(12) << y[3] << endl;


	myPipe.plot(y[0], y[2]);
	myPipe.plot2(y[0], y[2]);

	int point_count = 0; //initialize point counter

  // step to tmax from tmin 
  double h = 1e-6;		// starting step size for ode solver 
  for (double t_next = tmin + delta_t; t_next <= tmax; t_next += delta_t)
  {
		point_count++;
    while (t < t_next)	// evolve from t to t_next 
	  {
	    gsl_odeiv_evolve_apply (evolve_ptr, control_ptr, step_ptr,
			                    	  &my_system, &t, t_next, &h, y);

			if((point_count%plot_skip)==0)
			{
				myPipe.plot(y[0], y[2]);
				myPipe.plot2(y[0], y[2]);
			}
	  }

    // print at t = t_next
    ode_out << scientific << setprecision (5) << setw (12) << t << " " 
            << setw (12) << y[0] << " " << setw (12) << y[1] << " " << setw(12) << y[2] << " " << setw(12) << y[3] <<  endl;
  }
  ode_out.close();

	cout << "Wait while the pipe is closed..." << flush;
	myPipe.finish(); //close the pipe to gnuplot
	cout << "All done!" << endl;

  // all done; free up the gsl_odeiv stuff 
  gsl_odeiv_evolve_free (evolve_ptr);
  gsl_odeiv_control_free (control_ptr);
  gsl_odeiv_step_free (step_ptr);

  return 0;
}

//*************************** rhs ****************************
// 
// Define the array of right-hand-side functions y[i] to be integrated.
//  The equations are:
// x' = vx                  ==>  dy[0]/dt = f[0] = y[1]
// vx' = -GMx/(x^2 + y^2)^3/2 ==>  dy[1]/dt = f[1] = -GMy[0]/((y[0])^2 + (y[2])^2)^3/2
//
// y' = vy                  ==>  dy[2]/dt = f[2] = y[3]
// vy' = -GMy/(x^2 + y^2)^3/2 ==>  dy[3]/dt = f[1] = -GMy[2]/((y[0])^2 + (y[2])^2)^3/2
//
//  * params is a void pointer that is used in many GSL routines
//     to pass parameters to a function
//
int
rhs (double , const double y[], double f[], void *params_ptr)
{
  // get parameter(s) from params_ptr

	struct Params* passed_ptr;
	passed_ptr = (struct Params *) params_ptr;
 
  double G = passed_ptr->G;
	double M = passed_ptr->mass;

  // evaluate the right-hand-side functions at t 
  f[0] = y[1];
  f[1] = -G*M*y[0]/pow(y[0]*y[0] + y[2]*y[2], 1.5);
	f[2] = y[3];
  f[3] = -G*M*y[2]/pow(y[0]*y[0] + y[2]*y[2], 1.5);

  return GSL_SUCCESS;		// GSL_SUCCESS defined in gsl/errno.h as 0 
}

/*	Rk4 stepping does not need the Jacobian matrix
int
rhsy (double , const double y[], double f[], void *params_ptr)
{
  // get parameter(s) from params_ptr; here, just a double 
  double G = *(double *) params_ptr->G;
	double M = *(double *) params_ptr->mass;

  // evaluate the right-hand-side functions at t 
  f[0] = y[1];
  f[1] = -G*M*y[1]/pow(x[0]*x[0] + y[0]*y[0], 1.5);

  return GSL_SUCCESS;		// GSL_SUCCESS defined in gsl/errno.h as 0 
}

*************************** Jacobian ****************************
//
// Define the Jacobian matrix using GSL matrix routines.
//  (see the GSL manual under "Ordinary Differential Equations") 
//
//  * params is a void pointer that is used in many GSL routines
//     to pass parameters to a function
//
int
xjacobian (double , const double x[], const double y[], double *dfdy,
	  double dfdt[], void *params_ptr)
{
  // get parameter(s) from params_ptr 
  double G = *(double *) params_ptr->G;
	double M = *(double *) params_ptr->mass;

  gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdy, 2, 2);

  gsl_matrix *m_ptr = &dfdy_mat.matrix;	// m_ptr points to the matrix 

  // fill the Jacobian matrix as shown 
  gsl_matrix_set (m_ptr, 0, 0, 0.0);	// df[0]/dy[0] = 0 
  gsl_matrix_set (m_ptr, 0, 1, 1.0);	// df[0]/dy[1] = 1 
  gsl_matrix_set (m_ptr, 1, 0, G*M*(2*x[0]*x[0]-y[0]*y[0])/(pow(x[0]*x[0]+y[0]*y[0],5./2.)));	// df[1]/dy[0] 
  gsl_matrix_set (m_ptr, 1, 1, 0);	    // df[1]/dy[1] 

  // set explicit t dependence of f[i] 
  dfdt[0] = 0.0;
  dfdt[1] = 0.0;

  return GSL_SUCCESS;		// GSL_SUCCESS defined in gsl/errno.h as 0 
}

int
yjacobian (double , const double y[], double *dfdy,
	  double dfdt[], void *params_ptr)
{
  // get parameter(s) from params_ptr 
  double G = *(double *) params_ptr->G;
	double M = *(double *) params_ptr->mass;

  gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdy, 2, 2);

  gsl_matrix *m_ptr = &dfdy_mat.matrix;	// m_ptr points to the matrix 

  // fill the Jacobian matrix as shown 
  gsl_matrix_set (m_ptr, 0, 0, 0.0);	// df[0]/dy[0] = 0 
  gsl_matrix_set (m_ptr, 0, 1, 1.0);	// df[0]/dy[1] = 1 
  gsl_matrix_set (m_ptr, 1, 0, G*M*(2*y[0]*y[0]-x[0]*x[0])/(pow(x[0]*x[0]+y[0]*y[0],5./2.)));	// df[1]/dy[0] 
  gsl_matrix_set (m_ptr, 1, 1, 0);	    // df[1]/dy[1] 

  // set explicit t dependence of f[i] 
  dfdt[0] = 0.0;
  dfdt[1] = 0.0;

  return GSL_SUCCESS;		// GSL_SUCCESS defined in gsl/errno.h as 0 
}
*/





