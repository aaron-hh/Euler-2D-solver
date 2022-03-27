#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>

typedef std::vector<double> Vectorofdouble;

class eos
{
private:
    double m_gamma;
    double compute_pressure(double density, double v_x, double v_y, double energy, double m_gamma);
    double compute_energy(double pressure, double density, double v_x, double v_y, double m_gamma);

public:
    eos();
    //~eos();

    void set_gamma(double GAMMA);
    void compute_wavespeed(double cs, double pressure, double density, double m_gamma);
    virtual Vectorofdouble prim_to_u(Vectorofdouble prim, double gamma, int nVar, int D);
    virtual Vectorofdouble u_to_prim(Vectorofdouble u, double gamma, int nVar, int D);
    Vectorofdouble fluxf(Vectorofdouble u, Vectorofdouble prim, int nVar, int D);
    Vectorofdouble y_fluxf(Vectorofdouble u, Vectorofdouble prim, int nVar, int D);
};
