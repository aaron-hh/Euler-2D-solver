/* -------------------------------------------------------------------*/
/*                                                                    */
/*  File Description : This file contains the functions for computing */
/*      time step, defining initial conditions, defining boundary     */
/*      condition and iteration loops.                                */
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

#include "system_Euler.H"
#include "eos_Euler.H"
#include "numerical_method.H"
#include "array.H"

typedef std::vector<double> Vectorofdouble;
typedef std::vector<std::vector<double> >VVectorofdouble;

using namespace std;

//initiating classes
numerical_method* NM = new numerical_method();
eos* E = new eos();
m_array U;
m_array PRIM;
m_array UPLUS1;
m_array UIL_NPLUSHALF;
m_array UIR_NPLUSHALF;
m_array PRIML_NPLUSHALF;
m_array PRIMR_NPLUSHALF;
m_array WAVESPEED;
m_array FLUX;

//extra data storage for 2D
m_array UBAR;
m_array PRIMBAR;
m_array UIL_NPLUSHALF_BAR;
m_array UIR_NPLUSHALF_BAR;
m_array PRIML_NPLUSHALF_BAR;
m_array PRIMR_NPLUSHALF_BAR;
m_array WAVESPEED_BAR;
m_array YFLUX;

//allocating memory space for static
double system::m_nVar;
int system::m_xnCells;
double system::m_x0;
double system::m_x1;
double system::m_dx;
double system::m_x_midpoint;
double system::m_c1;
double system::m_c2;
double system::m_gamma;
double system::m_tstart;
double system::m_tend;
int system::m_d;

//2D
int system::m_ynCells;
double system::m_y0;
double system::m_y1;
double system::m_dy;
double system::m_y_midpoint;

//constructor
system::system()
{}

//define a function for parameter settings
void system::setparameters(double NVAR, int XNCELLS, double X0, double X1, double X_MIDPOINT, 
double C1, double C2, double GAMMA, double TSTART, double TEND, int DIMENSION, int YNCELLS,
double Y0, double Y1, double Y_MIDPOINT)
{
    m_nVar = NVAR;
    m_xnCells = XNCELLS;
    m_x0 = X0;
    m_x1 = X1;
    m_dx = (X1-X0)/XNCELLS;
    m_x_midpoint = X_MIDPOINT;
    m_c1 = C1;
    m_c2 = C2;
    m_gamma = GAMMA;
    m_tstart = TSTART;
    m_tend = TEND; 
    m_d = DIMENSION;

    //2D
    m_ynCells = YNCELLS;
    m_y0 = Y0;
    m_y1 = Y1;
    m_dy = (Y1-Y0)/YNCELLS;
    m_y_midpoint = Y_MIDPOINT;

    U.setSize(m_xnCells+4, m_ynCells+4, 1, m_nVar);
    PRIM.setSize(m_xnCells+4,  m_ynCells+4, 1, m_nVar);
    UPLUS1.setSize(m_xnCells+4,  m_ynCells+4, 1, m_nVar);
    UIL_NPLUSHALF.setSize(m_xnCells+2,  m_ynCells+2, 1, m_nVar);
    UIR_NPLUSHALF.setSize(m_xnCells+2,  m_ynCells+2, 1, m_nVar);
    PRIML_NPLUSHALF.setSize(m_xnCells+2,  m_ynCells+2, 1, m_nVar);
    PRIMR_NPLUSHALF.setSize(m_xnCells+2,  m_ynCells+2, 1, m_nVar);
    WAVESPEED.setSize(m_xnCells+1,  m_ynCells+1, 1, 3);
    FLUX.setSize(m_xnCells+1,  m_ynCells+1, 1, m_nVar);

    //2D
    UBAR.setSize(m_xnCells+4,  m_ynCells+4, 1, m_nVar);
    PRIMBAR.setSize(m_xnCells+4,  m_ynCells+4, 1, m_nVar);
    UIL_NPLUSHALF_BAR.setSize(m_xnCells+2,  m_ynCells+2, 1, m_nVar);
    UIR_NPLUSHALF_BAR.setSize(m_xnCells+2,  m_ynCells+2, 1, m_nVar);
    PRIML_NPLUSHALF_BAR.setSize(m_xnCells+2,  m_ynCells+2, 1, m_nVar);
    PRIMR_NPLUSHALF_BAR.setSize(m_xnCells+2,  m_ynCells+2, 1, m_nVar);
    WAVESPEED_BAR.setSize(m_xnCells+1,  m_ynCells+1, 1, 3);
    YFLUX.setSize(m_xnCells+1,  m_ynCells+1, 1, m_nVar);
}

//compute time step
void system::computeTimeStep(double m_dx, double m_c, double m_xnCells, double m_ynCells, double m_gamma, double& dt)
{

  Vector_vectorofdouble Cs(m_xnCells+4,Vectorofdouble(m_ynCells+4));
  Vector_vectorofdouble div(m_xnCells+4,Vectorofdouble(m_ynCells+4));
  Vector_vectorofdouble v_m(m_xnCells+4,Vectorofdouble(m_ynCells+4));

  //computing Cs 
  for(int i=0; i<m_xnCells+4; i++)
  {
    for(int j=0; j<m_ynCells+4; j++)
    {
    div[i][j] = PRIM(i,j,0,3) / PRIM(i,j,0,0);
    Cs[i][j] = sqrt(m_gamma*div[i][j]);
    v_m[i][j] = sqrt(pow(PRIM(i,j,0,1),2) + pow(PRIM(i,j,0,2),2));
    }
  }

  //finding a_max
  double a_max;
  a_max = 0;
    for(int i=0; i<m_xnCells+4; i++)
    {
      for(int j=0; j<m_ynCells+4; j++)
      {
        if((Cs[i][j]+v_m[i][j])>a_max)
        {
        a_max = Cs[i][j]+v_m[i][j];
        }
      }
    }

  dt = m_c*m_dx/a_max;
}

//setting initial data in primitive form
void system::computeInitialCondition(int m_d, double RHO_L_I, double VX_L_I, double VY_L_I, double P_L_I, double RHO_R_I, double VX_R_I, double VY_R_I, double P_R_I)
{
  Vectorofdouble IC_L(m_nVar), IC_R(m_nVar);
  IC_L[0] = RHO_L_I;
  IC_L[1] = VX_L_I;  
  IC_L[m_d+1] = P_L_I;   

  IC_R[0] = RHO_R_I;
  IC_R[1] = VX_R_I;  
  IC_R[m_d+1] = P_R_I;

  if(m_d == 2)
  {
  IC_L[2] = VY_L_I;
  IC_R[2] = VY_R_I;
  }    

  for(int i=0; i<m_xnCells+4; i++)
  {
    for(int j=0; j<m_ynCells+4; j++)
    {
      double x = m_x0 + (i-1)*m_dx;
      double y = m_y0 + (j-1)*m_dy;

      if(sqrt(pow((x-1),2)+pow((y-1),2))<0.4)
      {
        PRIM(i,j,0) = IC_L;
      }
      else
      {
        PRIM(i,j,0) = IC_R;
      }
    }
  }
}

//defining function for applying boundary condition (transmission)
void system::applyBoundaryCondition(double m_xnCells, double m_ynCells, int D)
{   
  if(D == 2)
  {
  //transmissive top bottom
  for(int i=0; i<m_xnCells+4; i++)
  {
  U(i,0,0) = U(i,2,0);
  U(i,1,0) = U(i,2,0);
  U(i,m_ynCells+2,0) = U(i,m_ynCells+1,0);
  U(i,m_ynCells+3,0) = U(i,m_ynCells+1,0);

  PRIM(i,0,0) = PRIM(i,2,0);
  PRIM(i,1,0) = PRIM(i,2,0);
  PRIM(i,m_ynCells+2,0) = PRIM(i,m_ynCells+1,0);
  PRIM(i,m_ynCells+3,0) = PRIM(i,m_ynCells+1,0);
  }
  }

  //transmissive left right
  for(int j=0; j<m_ynCells+4; j++) 
  {
  U(0,j,0) = U(2,j,0);
  U(1,j,0) = U(2,j,0);
  U(m_xnCells+2,j,0) = U(m_xnCells+1,j,0);
  U(m_xnCells+3,j,0) = U(m_xnCells+1,j,0);

  PRIM(0,j,0) = PRIM(2,j,0);
  PRIM(1,j,0) = PRIM(2,j,0);
  PRIM(m_xnCells+2,j,0) = PRIM(m_xnCells+1,j,0);
  PRIM(m_xnCells+3,j,0) = PRIM(m_xnCells+1,j,0);
  }
}

//defining function for iteration
void system::iteration()
{
double t = m_tstart;
double dt;
int counter;

//convert from primitive to conservative
for(int i=0; i<m_xnCells+4; i++)
{
  for(int j=0; j<m_ynCells+4; j++)
  {
    E->set_gamma(m_gamma);
    U(i,j,0) = E->prim_to_u(PRIM(i,j,0), m_gamma, m_nVar, m_d);
  }
}

do
{
  //use lower CFL number for the first few iteration for stability
  counter +=1;
  double C;
  if(counter<20)
  {
    C = m_c1;
  }
  else
  {
    C = m_c2;
  }

  //compute the stable time step for this iteration
  system* S = new system;
  S->computeTimeStep(m_dx, C, m_xnCells, m_ynCells, m_gamma, dt);
  t = t + dt;
  std::cout<<t<<std::endl;

  //applying boundary condition
  S->applyBoundaryCondition(m_xnCells, m_ynCells, m_d);

  //calculate uiL_nplushalf, uiR_nplushalf for x cells
   for(int i = 0; i<m_xnCells+2; i++)
   {
     for(int j = 0; j<m_ynCells+2; j++)
     {
       Vectorofdouble a_x = U(i,j+2,0);
       Vectorofdouble b_x = U(i+1,j+2,0);
       Vectorofdouble c_x = U(i+2,j+2,0);
      
       NM->compute_nplushalf_variables(a_x, b_x, c_x, m_dx, dt, m_gamma, UIL_NPLUSHALF(i,j,0), UIR_NPLUSHALF(i,j,0), m_nVar, m_d);

       E->set_gamma(m_gamma);
       PRIML_NPLUSHALF(i,j,0) = E->u_to_prim(UIL_NPLUSHALF(i,j,0),m_gamma,m_nVar,m_d);
       PRIMR_NPLUSHALF(i,j,0) = E->u_to_prim(UIR_NPLUSHALF(i,j,0),m_gamma,m_nVar,m_d);
     }
  }

  //calculate xHLLCflux
  for(int i = 0; i<m_xnCells+1; i++)
  {
    for(int j = 0; j<m_ynCells+1; j++)
    {
    double x = m_x0 + (i-1)*m_dx;
    NM->wavespeedestimate(PRIMR_NPLUSHALF(i,j,0), PRIML_NPLUSHALF(i+1,j,0), m_xnCells, m_gamma, WAVESPEED(i,j,0), m_d);
    NM->compute_HLLCflux(UIR_NPLUSHALF(i,j,0), PRIMR_NPLUSHALF(i,j,0), UIL_NPLUSHALF(i+1,j,0), 
    PRIML_NPLUSHALF(i+1,j,0), WAVESPEED(i,j,0), FLUX(i,j,0), m_nVar, m_gamma, m_d, x, t);

    // FLUX(i,j,0) = NM->getFORCEflux_bar(UIR_NPLUSHALF(i,j,0), UIL_NPLUSHALF(i+1,j,0), PRIMR_NPLUSHALF(i,j,0), 
    // PRIML_NPLUSHALF(i+1,j,0), m_dx, dt, m_gamma, 4, m_d);
    }
  }

  // x update the data using MUSCL Hancock from 2 to 102
  for(int i = 2; i<m_xnCells+2; i++)
  {
    for(int j = 2; j<m_ynCells+2; j++)
    {
      for(int k = 0; k<m_nVar; k++)
      {
      UBAR(i,j,0,k) = U(i,j,0,k) - (dt/m_dx) * (FLUX(i-1,j-2,0,k) - FLUX(i-2,j-2,0,k));
      }
    }
  }

  //transmissive top bottom
  if(m_d == 2)
  {
  for(int i=0; i<m_xnCells+4; i++)
  {
  UBAR(i,0,0) = UBAR(i,2,0);
  UBAR(i,1,0) = UBAR(i,2,0);
  UBAR(i,m_ynCells+2,0) = UBAR(i,m_ynCells+1,0);
  UBAR(i,m_ynCells+3,0) = UBAR(i,m_ynCells+1,0);
  }
  }

  //transmissive left right
  for(int j=0; j<m_ynCells+4; j++)
  {
  UBAR(0,j,0) = UBAR(2,j,0);
  UBAR(1,j,0) = UBAR(2,j,0);
  UBAR(m_xnCells+2,j,0) = UBAR(m_xnCells+1,j,0);
  UBAR(m_xnCells+3,j,0) = UBAR(m_xnCells+1,j,0);
  }

  //convert ubar to primbar
  for(int i=0; i<m_xnCells+4; i++)
  {
    for(int j=0; j<m_ynCells+4; j++)
    {
    E->set_gamma(m_gamma);
    PRIMBAR(i,j,0) = E->u_to_prim(UBAR(i,j,0),m_gamma,m_nVar,m_d);
    }
  }

  //calculate uiL_nplushalf, uiR_nplushalf for y cells
   for(int i = 0; i<m_xnCells+2; i++)
   {
     for(int j = 0; j<m_ynCells+2; j++)
     {
       Vectorofdouble a_y = UBAR(i+2,j,0);
       Vectorofdouble b_y = UBAR(i+2,j+1,0);
       Vectorofdouble c_y = UBAR(i+2,j+2,0);

       E->set_gamma(m_gamma);
       NM->ycompute_nplushalf_variables(a_y, b_y, c_y, m_dy, dt, m_gamma, UIL_NPLUSHALF_BAR(i,j,0), UIR_NPLUSHALF_BAR(i,j,0), m_nVar, m_d);

       PRIML_NPLUSHALF_BAR(i,j,0) = E->u_to_prim(UIL_NPLUSHALF_BAR(i,j,0),m_gamma,m_nVar,m_d);
       PRIMR_NPLUSHALF_BAR(i,j,0) = E->u_to_prim(UIR_NPLUSHALF_BAR(i,j,0),m_gamma,m_nVar,m_d);
     }
  }

  //calculate yHLLCflux
  for(int i = 0; i<m_xnCells+1; i++)
  {
    for(int j = 0; j<m_ynCells+1; j++)
    {
    double y = m_y0 + (j-1)*m_dy;
    NM->ywavespeedestimate(PRIMR_NPLUSHALF_BAR(i,j,0), PRIML_NPLUSHALF_BAR(i,j+1,0), m_ynCells, m_gamma, WAVESPEED_BAR(i,j,0), m_d);
    NM->ycompute_HLLCflux(UIR_NPLUSHALF_BAR(i,j,0), PRIMR_NPLUSHALF_BAR(i,j,0), UIL_NPLUSHALF_BAR(i,j+1,0), 
    PRIML_NPLUSHALF_BAR(i,j+1,0), WAVESPEED_BAR(i,j,0), YFLUX(i,j,0), m_nVar, m_gamma, m_d, y, t);

    // YFLUX(i,j,0) = NM->ygetFORCEflux_bar(UIR_NPLUSHALF_BAR(i,j,0), UIL_NPLUSHALF_BAR(i,j+1,0), PRIMR_NPLUSHALF_BAR(i,j,0), 
    // PRIML_NPLUSHALF_BAR(i,j+1,0), m_dy, dt, m_gamma, 4, m_d);
    }
  }

  // y update the data using MUSCL Hancock
  for(int i = 2; i<m_xnCells+2; i++)
  {
    for(int j = 2; j<m_ynCells+2; j++)
    {
      for(int k = 0; k<m_nVar; k++)
      {
      U(i,j,0,k) = UBAR(i,j,0,k) - (dt/m_dy) * (YFLUX(i-2,j-1,0,k) - YFLUX(i-2,j-2,0,k));
      }
    }
  }

  //convert u to prim
  for(int i=0; i<m_xnCells+4; i++)
  {
    for(int j=0; j<m_ynCells+4; j++)
    {
    E->set_gamma(m_gamma);
    PRIM(i,j,0) = E->u_to_prim(U(i,j,0),m_gamma,m_nVar,m_d);
    }
  }

} while (t<m_tend);
}

void system::outputData()
{
//delete file from previous run
std::remove("ttEuler_2D.dat"); 

//outputing data
std::ofstream output ("ttEuler_2D.dat");
for(int i=2; i<m_xnCells+2; i++)
{
  for(int j=2; j<m_ynCells+2; j++)
  {
  double x = m_x0 + (i-1)*m_dx;
  double y = m_y0 + (j-1)*m_dy;
  output<<x<<" "<<y<<" "<<PRIM(i,j,0,0)<<" "<<PRIM(i,j,0,1)<<" "<<PRIM(i,j,0,2)<<" "<<PRIM(i,j,0,3)<<std::endl;
  }
  output<<std::endl;
}
}
