/* -------------------------------------------------------------------*/
/*                                                                    */
/*  File Description : This file contains the ideal gas equation of   */
/*        state functions.                                            */
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
#include "eos_Euler.H"

typedef std::vector<double> Vectorofdouble;

//constructor 
eos::eos()
{}

//setter for gamma 
void eos::set_gamma(double GAMMA)
{
    m_gamma = GAMMA;
}

//defining compute pressure function
double eos::compute_pressure(double density, double v_x, double v_y, double energy, double m_gamma)
{
    double pressure = (m_gamma-1) * (energy - 0.5*density*((v_x*v_x) + (v_y*v_y)));

    return pressure;
}

//defining compute energy function
double eos::compute_energy(double pressure, double density, double v_x, double v_y, double m_gamma)
{   
    double energy = pressure / (m_gamma-1) + 0.5*density*((v_x*v_x) + (v_y*v_y));

    return energy;
}

//defining compute wave speed function
void eos::compute_wavespeed(double cs, double pressure, double density, double m_gamma)
{
    cs = sqrt(pressure * m_gamma / density);
}

//defining primitive to conservative function
Vectorofdouble eos::prim_to_u(Vectorofdouble prim, double m_gamma, int nVar, int D) 
{
    Vectorofdouble u(nVar);

    u[0] = prim[0];
    u[1] = prim[1] * prim[0];
    eos* E = new eos();

    if(D==2)
    {
    u[2] = prim[2] * prim[0];
    }
    else if (D==1)
    {}

    u[D+1] = E->compute_energy(prim[D+1], prim[0], prim[1], prim[D], m_gamma);

    return u;
}

//defining conservative to primitive function
Vectorofdouble eos::u_to_prim(Vectorofdouble u, double m_gamma, int nVar, int D)
{  
    Vectorofdouble prim(nVar);

    prim[0] = u[0];
    prim[1] = u[1]/u[0];

    if(D==2)
    {
    prim[2] = u[2]/u[0];
    }
    else if (D==1)
    {}

    eos* E = new eos();
    prim[D+1] = E->compute_pressure(prim[0], prim[1], prim[D], u[D+1], m_gamma);


    return prim;
}

//defining flux function
Vectorofdouble eos::fluxf(Vectorofdouble u, Vectorofdouble prim, int nVar, int D)
{
    Vectorofdouble x_fluxfunction(nVar);

    x_fluxfunction[0] = u[1];
    x_fluxfunction[1] = u[1] * prim[1] + prim[D+1];
    x_fluxfunction[D+1] = (u[D+1] + prim[D+1]) * prim[1]; 

    if(D==2)
    {
    x_fluxfunction[2] = prim[0] * prim[1] * prim[2];
    }
    else if (D==1)
    {}

    return x_fluxfunction;
}

//defining flux function
Vectorofdouble eos::y_fluxf(Vectorofdouble u, Vectorofdouble prim, int nVar, int D)
{
  Vectorofdouble y_fluxfunction(nVar);

  y_fluxfunction[0] = u[2];
  y_fluxfunction[1] = prim[0] * prim[1] * prim[2];
  y_fluxfunction[2] = u[2] * prim[2] + prim[3]; 
  y_fluxfunction[3] = (u[3] + prim[3]) * prim[2]; 

  return y_fluxfunction;
}

