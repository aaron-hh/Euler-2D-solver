#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>

typedef std::vector<double> Vectorofdouble;
typedef std::vector<std::vector<double> > Vector_vectorofdouble;

class system
{
private:
    static double m_nVar;
    static int m_xnCells;
    static double m_x0;
    static double m_x1;
    static double m_dx;
    static double m_x_midpoint;
    static double m_c1;
    static double m_c2;
    static double m_gamma;
    static double m_tstart;
    static double m_tend;
    static int m_d;

    static int m_ynCells;
    static double m_y0;
    static double m_y1;
    static double m_dy;
    static double m_y_midpoint;
    //static int m_total_cells;

    void computeTimeStep(double m_dx, double m_c, double m_xnCells, double m_ynCells, double m_gamma, double& dt);
    void applyBoundaryCondition(double m_xnCells, double m_ynCells, int m_d);

public:

    system(); 
    void setparameters(double NVAR, int XNCELLS, double X0, double X1, double X_MIDPOINT, 
    double C1, double C2, double GAMMA, double TSTART, double TEND, int DIMENSION, int YNCELLS,
    double Y0, double Y1, double Y_MIDPOINT);
    void computeInitialCondition(int m_d, double RHO_L_I, double VX_L_I, double VY_L_I, double P_L_I, double RHO_R_I, double VX_R_I, double VY_R_I, double P_R_I);
    void iteration();
    void outputData();
};