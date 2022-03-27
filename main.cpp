/* -------------------------------------------------------------------*/
/*                                                                    */
/*  File Description : This file contains the input data required for */ 
/*          running the Euler 2-D simulation and the main function.   */
/*                                                                    */
/*  Programer: Aaron Ng                                               */
/*                                                                    */
/*  Last Revision: 01 February 2022                                   */
/*                                                                    */
/* -------------------------------------------------------------------*/

#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <fstream>
#include <array>

#include "system_Euler.H"
#include "numerical_method.H"
#include "array.H"

typedef std::vector<double> Vectorofdouble;
typedef std::vector<std::vector<double> > Vector_vectorofdouble;

double NVAR = 4;
int XNCELLS = 100;
double X0 = 0;
double X1 = 2;
double X_MIDPOINT = 1;
double C2 = 0.85;
double C1 = 0.2*C2;
double GAMMA = 1.4;
double TSTART = 0; 
double TEND = 0.25;
int DIMENSION = 2;

int YNCELLS = 100;
double Y0 = 0;
double Y1 = 2;
double Y_MIDPOINT = 1;

double RHO_L_I = 1.0;
double VX_L_I = 0.0;
double VY_L_I = 0;
double P_L_I = 1.0;
double RHO_R_I = 0.125;
double VX_R_I = 0.0;
double VY_R_I = 0;
double P_R_I = 0.1;

int main()
{
    class system* S = new class system();
    S->setparameters(NVAR, XNCELLS, X0, X1, X_MIDPOINT, C1, C2, GAMMA, TSTART, TEND, DIMENSION,
    YNCELLS, Y0, Y1, Y_MIDPOINT);
    S->computeInitialCondition(DIMENSION, RHO_L_I, VX_L_I, VY_L_I, P_L_I, RHO_R_I, VX_R_I, VY_R_I, P_R_I);
    S->iteration();
    S->outputData();
}






