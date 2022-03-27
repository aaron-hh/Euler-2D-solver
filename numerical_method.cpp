/* -------------------------------------------------------------------*/
/*                                                                    */
/*  File Description : This file contains the functions for SLIC and  */
/*      MUSCL Hancock solver (FORCE & HLLC).                          */
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
#include "numerical_method.H"

eos* E1 = new eos();

//constructor
numerical_method::numerical_method()
{}

//function to find Minbee constant
Vectorofdouble numerical_method::Minbee(Vectorofdouble ui, Vectorofdouble uiMinus1, Vectorofdouble uiPlus1, int nVar)
{
  Vectorofdouble r(nVar);
  for(int i=0; i<4; i++)
  {
    if (fabs(uiPlus1[i] - ui[i])>0)
    {
      r[i] = (ui[i] - uiMinus1[i]) / (uiPlus1[i] - ui[i]);
    }
    else 
    {
      r[i] = 0;
    }
  }

  Vectorofdouble Mb(nVar);
  for(int i=0; i<4; i++)
  {
    if(r[i]<=0)
    {
      Mb[i] = 0;
    }
    else if(r[i]>0 && r[i]<=1)
    {
      Mb[i] = r[i];
    }
    else
    {
      if((2/(1+r[i])>1))
      {
        Mb[i] = 2/(1+r[i]);
      }
      else
      {
        Mb[i] = 1;
      }
    }
  }
  return Mb;
}

//function to find Vanleer constant
Vectorofdouble numerical_method::Vanleer(Vectorofdouble ui, Vectorofdouble uiMinus1, Vectorofdouble uiPlus1, int nVar)
{
  Vectorofdouble r(nVar);
  Vectorofdouble xi(nVar);
  Vectorofdouble Vl(nVar);

  for(int i=0; i<nVar; i++)
  {
    if (fabs(uiPlus1[i] - ui[i])>0)
    {
      r[i] = (ui[i] - uiMinus1[i]) / (uiPlus1[i] - ui[i]);
    }
    else 
    {
      r[i] = 0;
    }
    xi[i] = 2/(1+r[i]);
    //r[i] = 0;
  }

  for(int i=0; i<nVar; i++)
  {
    if(r[i]<=0)
    {
      Vl[i] = 0;
    }
    else if(r[i]>0)
    {
      if(xi[i]>(2*r[i]/(1+r[i])))
      {
        Vl[i] = 2*r[i]/(1+r[i]);
      }
      else 
      {
        Vl[i] = xi[i];
      }
    }
  }
  return Vl;
}

//defining wavespeedestimate
//p represents primR_nplushalf[i], q represents primL_nplushalf[i+1]
//Sl(minimum wave speed from minimum vx - maximum Cf), Sr(maximum wave speed from maximum vx + maximum Cf), S_star (intermediate wave speed)
void numerical_method::wavespeedestimate(Vectorofdouble p, Vectorofdouble q, double n, double gamma, Vectorofdouble& wavespeed, int D)
{
  double rho_l = p[0];
  double rho_r = q[0];
  double vx_l = p[1];
  double vx_r = q[1];
  double p_l = p[D+1];
  double p_r = q[D+1];
  double a_l = sqrt(gamma*p_l/rho_l);
  double a_r = sqrt(gamma*p_r/rho_r);
  double vy_l;
  double vy_r;

  double S_lr = std::max(abs(vx_l)+a_l, abs(vx_r)+a_r);
  double Sl = - S_lr;
  double Sr = + S_lr;

  double S_star = (p_r-p_l+((rho_l*vx_l)*(Sl-vx_l))-((rho_r*vx_r)*(Sr-vx_r)))/((rho_l*(Sl-vx_l))-(rho_r*(Sr-vx_r)));

  wavespeed[0] = Sl;
  wavespeed[1] = Sr;
  wavespeed[2] = S_star;
}

void numerical_method::ywavespeedestimate(Vectorofdouble p, Vectorofdouble q, double n, double gamma, Vectorofdouble& wavespeed, int D)
{
  double rho_l = p[0];
  double rho_r = q[0];
  double vx_l = p[1];
  double vx_r = q[1];
  double p_l = p[D+1];
  double p_r = q[D+1];
  double a_l = sqrt(gamma*p_l/rho_l);
  double a_r = sqrt(gamma*p_r/rho_r);
  double vy_l;
  double vy_r;

  if(D == 1)
  {
  vy_l = 0;
  vy_r = 0;
  }
  else if(D == 2)
  {
  vy_l = p[D];
  vy_r = q[D];
  }

  double S_lr = std::max(abs(vy_l)+a_l, abs(vy_r)+a_r);
  double Sl = - S_lr;
  double Sr = + S_lr;

  double S_star = (p_r-p_l + ((rho_l*vy_l)*(Sl-vy_l)) - ((rho_r*vy_r)*(Sr-vy_r))) / ((rho_l*(Sl-vy_l))-(rho_r*(Sr-vy_r)));

  wavespeed[0] = Sl;
  wavespeed[1] = Sr;
  wavespeed[2] = S_star;
}

//computing nplushalf
void numerical_method::compute_nplushalf_variables(Vectorofdouble a, Vectorofdouble b, Vectorofdouble c, double dx, double dt, double gamma, Vectorofdouble& uiL_nplushalf, Vectorofdouble& uiR_nplushalf, int nVar, int D)
{
  Vectorofdouble ui(nVar), uiMinus1(nVar), uiPlus1(nVar), Vl(nVar), Mb(nVar), delta(nVar);
  Vectorofdouble uiL(nVar), uiR(nVar), primiL(nVar), primiR(nVar), fuiL(nVar), fuiR(nVar);

  //defining ui, uiMinus1, uiPlus1
  for(int i=0; i<nVar; i++)
  {
  uiMinus1[i] = a[i];
  ui[i] = b[i];
  uiPlus1[i] = c[i];
  }

  //finding Vl
  numerical_method *Vanleer_ui = new numerical_method();
  Vl = Vanleer_ui->Vanleer(ui, uiMinus1, uiPlus1, nVar);
  double Vl_min;
  Vl_min = *std::min_element(Vl.begin(), Vl.end());

  //finding Mb
  numerical_method *Minbee_ui = new numerical_method();
  Mb = Minbee_ui->Minbee(ui, uiMinus1, uiPlus1, nVar);

  //define delta
  for(int i=0; i<nVar; i++)
  {
    delta[i] = 0.5*(uiPlus1[i]-uiMinus1[i]);
  }

  //define uiL, uiR
  for(int i=0; i<nVar; i++)
  {
    uiL[i] = ui[i] - 0.5*Vl_min*delta[i];
    uiR[i] = ui[i] + 0.5*Vl_min*delta[i];
  }

  //convert uiL, uiR to primiL, primiR
  E1->set_gamma(gamma);
  primiL = E1->u_to_prim(uiL, gamma, nVar, D);
  primiR = E1->u_to_prim(uiR, gamma, nVar, D);

  fuiL = E1->fluxf(uiL, primiL, nVar, D);
  fuiR = E1->fluxf(uiR, primiL, nVar, D);

  //calculate uiL_nplushalf, uiR_nplushalf
  for(int i=0; i<nVar; i++)
  {
    uiL_nplushalf[i] = uiL[i] - 0.5*dt/dx*(fuiR[i]-fuiL[i]);
    uiR_nplushalf[i] = uiR[i] - 0.5*dt/dx*(fuiR[i]-fuiL[i]);
  }
}

//computing nplushalf
void numerical_method::ycompute_nplushalf_variables(Vectorofdouble a, Vectorofdouble b, Vectorofdouble c, double dx, double dt, double gamma, Vectorofdouble& uiL_nplushalf, Vectorofdouble& uiR_nplushalf, int nVar, int D)
{
  Vectorofdouble ui(nVar), uiMinus1(nVar), uiPlus1(nVar), Vl(nVar), Mb(nVar), delta(nVar);
  Vectorofdouble uiL(nVar), uiR(nVar), primiL(nVar), primiR(nVar), fuiL(nVar), fuiR(nVar);

  //defining ui, uiMinus1, uiPlus1
  for(int i=0; i<nVar; i++)
  {
  uiMinus1[i] = a[i];
  ui[i] = b[i];
  uiPlus1[i] = c[i];
  }

  //finding Vl
  numerical_method *Vanleer_ui = new numerical_method();
  Vl = Vanleer_ui->Vanleer(ui, uiMinus1, uiPlus1, nVar);
  double Vl_min;
  Vl_min = *std::min_element(Vl.begin(), Vl.end());

  //finding Mb
  numerical_method *Minbee_ui = new numerical_method();
  Mb = Minbee_ui->Minbee(ui, uiMinus1, uiPlus1, nVar);

  //define delta
  for(int i=0; i<nVar; i++)
  {
    delta[i] = 0.5*(uiPlus1[i]-uiMinus1[i]);
  }

  //define uiL, uiR
  for(int i=0; i<nVar; i++)
  {
    uiL[i] = ui[i] - 0.5*Vl_min*delta[i];
    uiR[i] = ui[i] + 0.5*Vl_min*delta[i];
  }

  //convert uiL, uiR to primiL, primiR
  E1->set_gamma(gamma);
  primiL = E1->u_to_prim(uiL, gamma, nVar, D);
  primiR = E1->u_to_prim(uiR, gamma, nVar, D);

  fuiL = E1->y_fluxf(uiL, primiL, nVar, D);
  fuiR = E1->y_fluxf(uiR, primiL, nVar, D);

  //calculate uiL_nplushalf, uiR_nplushalf
  for(int i=0; i<nVar; i++)
  {
    uiL_nplushalf[i] = uiL[i] - 0.5*dt/dx*(fuiR[i]-fuiL[i]);
    uiR_nplushalf[i] = uiR[i] - 0.5*dt/dx*(fuiR[i]-fuiL[i]);
  }
}

void numerical_method::compute_HLLCflux(Vectorofdouble a, Vectorofdouble p, Vectorofdouble b, Vectorofdouble q, Vectorofdouble wavespeed, Vectorofdouble& Fhllc, int nVar, double gamma, int D, double x, double t)
{
  //define input variables 
  double Sl = wavespeed[0];
  double Sr = wavespeed[1];
  double S_star = wavespeed[2];

  double rho_l = p[0];
  double vx_l = p[1];
  double p_l = p[D+1];
  double U_l = a[D+1];

  double rho_r = q[0];
  double vx_r = q[1];
  double p_r = q[D+1];
  double U_r = b[D+1];

  Vectorofdouble Fl(nVar),Fr(nVar);
  E1->set_gamma(gamma);
  Fl = E1->fluxf(a, p, nVar, D);
  Fr = E1->fluxf(b, q, nVar, D);

  //compute Ql* & Fl*, Qr* & Fr*
  Vectorofdouble Ql_star(nVar),Qr_star(nVar),Q_hllc(nVar),prim_hllc(nVar);
  Vectorofdouble Fl_star(nVar),Fr_star(nVar);

  Ql_star[0] = rho_l*((Sl-vx_l)/(Sl-S_star));
  Ql_star[1] = rho_l*((Sl-vx_l)/(Sl-S_star))*(S_star);
  Ql_star[D+1] = rho_l*((Sl-vx_l)/(Sl-S_star))*((U_l/rho_l)+(S_star-vx_l)*(S_star+(p_l/(rho_l*(Sl-vx_l)))));

  Qr_star[0] = rho_r*((Sr-vx_r)/(Sr-S_star));
  Qr_star[1] = rho_r*((Sr-vx_r)/(Sr-S_star))*(S_star);
  Qr_star[D+1] = rho_r*((Sr-vx_r)/(Sr-S_star))*((U_r/rho_r)+(S_star-vx_r)*(S_star+(p_r/(rho_r*(Sr-vx_r)))));

  if(D == 2)
  {
  double vy_l = p[2];
  double vy_r = q[2];

  Ql_star[2] = rho_l*((Sl-vx_l)/(Sl-S_star))*(vy_l);
  Qr_star[2] = rho_r*((Sr-vx_r)/(Sr-S_star))*(vy_r);
  }
  else 
  {}

  for(int i=0; i<nVar; i++)
  {
  Fl_star[i] = Fl[i]  + Sl*(Ql_star[i] - a[i]);
  Fr_star[i] = Fr[i]  + Sr*(Qr_star[i] - b[i]);
  }
	
	if(Sl >= 0)
	{	
		Fhllc = Fl;
		return; 
	}

	if(Sl < 0 && S_star >= 0)
	{
		Fhllc = Fl_star;
		return;
	}
	
	if(S_star < 0 && Sr > 0)
	{
		Fhllc = Fr_star;
		return;
	}

	if(Sr <= 0)
	{	Fhllc = Fr;
		return; 
	}		
}

void numerical_method::ycompute_HLLCflux(Vectorofdouble a, Vectorofdouble p, Vectorofdouble b, Vectorofdouble q, Vectorofdouble wavespeed, Vectorofdouble& Fhllc, int nVar, double gamma, int D, double x, double t)
{
  //define input variables 
  double Sl = wavespeed[0];
  double Sr = wavespeed[1];
  double S_star = wavespeed[2];

  double rho_l = p[0];
  double vx_l = p[1];
  double p_l = p[D+1];
  double U_l = a[D+1];

  double rho_r = q[0];
  double vx_r = q[1];
  double p_r = q[D+1];
  double U_r = b[D+1];

  Vectorofdouble Fl(nVar),Fr(nVar);
  E1->set_gamma(gamma);
  Fl = E1->y_fluxf(a, p, nVar, D);
  Fr = E1->y_fluxf(b, q, nVar, D);

  //compute Ql* & Fl*, Qr* & Fr*
  Vectorofdouble Ql_star(nVar),Qr_star(nVar),Q_hllc(nVar),prim_hllc(nVar);
  Vectorofdouble Fl_star(nVar),Fr_star(nVar);

  double vy_l = p[2];
  double vy_r = q[2];

  Ql_star[2] = rho_l*((Sl-vy_l)/(Sl-S_star))*(S_star);
  Qr_star[2] = rho_r*((Sr-vy_r)/(Sr-S_star))*(S_star);

  Ql_star[0] = rho_l*((Sl-vy_l)/(Sl-S_star));
  Ql_star[1] = rho_l*((Sl-vy_l)/(Sl-S_star))*(vx_l);
  Ql_star[D+1] = rho_l*((Sl-vy_l)/(Sl-S_star))*((U_l/rho_l)+(S_star-vy_l)*(S_star+(p_l/(rho_l*(Sl-vy_l)))));

  Qr_star[0] = rho_r*((Sr-vy_r)/(Sr-S_star));
  Qr_star[1] = rho_r*((Sr-vy_r)/(Sr-S_star))*(vx_r);
  Qr_star[D+1] = rho_r*((Sr-vy_r)/(Sr-S_star))*((U_r/rho_r)+(S_star-vy_r)*(S_star+(p_r/(rho_r*(Sr-vy_r)))));

  for(int i=0; i<nVar; i++)
  {
  Fl_star[i] = Fl[i]  + Sl*(Ql_star[i] - a[i]);
  Fr_star[i] = Fr[i]  + Sr*(Qr_star[i] - b[i]);
  }

  if(Sl >= 0)
  {	
    Fhllc = Fl;
    return; 
  }

  if(Sl <0 && S_star >= 0)
  {
    Fhllc = Fl_star;
    return;
  }
  
  if(S_star < 0 && Sr > 0)
  {
    Fhllc = Fr_star;
    return;
  }

  if(Sr <= 0)
  {	Fhllc = Fr;
    return; 
  }		
}

//defining getFORCEflux for single cell (q=prim; p=prim(i+1))
Vectorofdouble numerical_method::getFORCEflux_bar(Vectorofdouble a, Vectorofdouble b,Vectorofdouble p,Vectorofdouble q, 
double dx, double dt, double gamma, double nVar, int D)
{
  Vectorofdouble fluxfunc(nVar);
  Vectorofdouble fluxfuncplus1(nVar);
  Vectorofdouble LFflux(nVar);
  Vectorofdouble RIflux(nVar);
  Vectorofdouble RI(nVar);
  Vectorofdouble FORCEflux(nVar);

  E1->set_gamma(gamma);
  fluxfunc = E1->fluxf(a, p, nVar, D);
  fluxfuncplus1 = E1->fluxf(b, q, nVar, D);

  //FORCE method
  for(int k=0; k<nVar; k++)
  {
  LFflux[k] = (0.5*(dx/dt)*(a[k]-b[k]))+ (0.5*(fluxfunc[k]+fluxfuncplus1[k]));
  }

  // //RI in u(i+0.5)
  for(int i=0; i<nVar; i++)
  {
  RI[i] = (0.5*(a[i]+b[i])) - (0.5*(dt/dx)*(fluxfuncplus1[i]-fluxfunc[i]));
  }

  //converting RI using the same operations of u(i+0.5) to prim(i+0.5) conversion
  Vectorofdouble RI_prim;
  RI_prim = E1->u_to_prim(RI, gamma, nVar, D);

  //RI*prim(i+0.5)
  RIflux = E1->fluxf(RI, RI_prim, nVar, D);

  for(int i=0; i<4; i++)
  {
  FORCEflux[i] = 0.5*(LFflux[i]+RIflux[i]);
  }
  return FORCEflux;
}

//defining getFORCEflux for single cell (q=prim; p=prim(i+1))
Vectorofdouble numerical_method::ygetFORCEflux_bar(Vectorofdouble a, Vectorofdouble b,Vectorofdouble p,Vectorofdouble q, 
double dx, double dt, double gamma, double nVar, int D)
{
  Vectorofdouble fluxfunc(nVar);
  Vectorofdouble fluxfuncplus1(nVar);
  Vectorofdouble LFflux(nVar);
  Vectorofdouble RIflux(nVar);
  Vectorofdouble RI(nVar);
  Vectorofdouble FORCEflux(nVar);

  E1->set_gamma(gamma);
  fluxfunc = E1->y_fluxf(a, p, nVar, D);
  fluxfuncplus1 = E1->y_fluxf(b, q, nVar, D);

  //FORCE method
  for(int k=0; k<nVar; k++)
  {
  LFflux[k] = (0.5*(dx/dt)*(a[k]-b[k]))+ (0.5*(fluxfunc[k]+fluxfuncplus1[k]));
  }

  // //RI in u(i+0.5)
  for(int i=0; i<nVar; i++)
  {
  RI[i] = (0.5*(a[i]+b[i])) - (0.5*(dt/dx)*(fluxfuncplus1[i]-fluxfunc[i]));
  }

  //converting RI using the same operations of u(i+0.5) to prim(i+0.5) conversion
  Vectorofdouble RI_prim;
  RI_prim = E1->u_to_prim(RI, gamma, nVar, D);

  //RI*prim(i+0.5)
  RIflux = E1->y_fluxf(RI, RI_prim, nVar, D);

  for(int i=0; i<4; i++)
  {
  FORCEflux[i] = 0.5*(LFflux[i]+RIflux[i]);
  }
  return FORCEflux;
}
