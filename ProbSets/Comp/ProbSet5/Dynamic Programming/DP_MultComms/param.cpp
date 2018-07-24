/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
 * param.cpp
 * 
 * This file defines the parameters of the problem. 
 * It defines the preference parameters, 
 * the state space, the level of discretization, 
 * and the number of iterations. 
 * It also defines parameters useful for formatting the output files of the program.
 * 
 * Carlos Rangel & Simon Scheidegger -- 06/17.
 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/ 
#include "econ.hpp"

using namespace Eigen;

// Preference Parameters
double alpha=0.25;
double gamm=5.0;
double beta=0.9;


// Capital Stock Interval

// Minimum and Maximum Capital Stock
double kmin=0.85;
double kmax=1.15;

// Number of capital stocks and discretization level
int nk=36000;//31;
double kappa=(kmax-kmin)/(nk-1.0);

// Number of theta states and the minimum and maximum theta values
int ntheta=5;
double theta_min=1-10.0/100;
double theta_step=5.0/100;

// Capital and Theta grids
ArrayXd kgrid=fill_kgrid();
ArrayXd thetagrid=fill_thetagrid();


// Number of Iterations
int numstart=0; // iterations begin at this timestep
int Numits=2;  // iterations end at this timestep

// Maximum difference between succesive Value Function Iterates (initialized to 0)
double errmax=0;

// Number of Columns for Output files
int nout=3;

// Frequency with which data will be printed
int datafreq=1;

// Hybrid Switch ( 0 for Comm Splitting; 1 to use regular hybrid)
int hybrid_switch=0;




