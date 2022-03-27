typedef std::vector<double> Vectorofdouble;

class numerical_method
{
private:
    Vectorofdouble Vanleer(Vectorofdouble ui, Vectorofdouble uiMinus1, Vectorofdouble uiPlus1, int nVar);
    Vectorofdouble Minbee(Vectorofdouble ui, Vectorofdouble uiMinus1, Vectorofdouble uiPlus1, int nVar);

public:
    numerical_method();
    
    void wavespeedestimate(Vectorofdouble p, Vectorofdouble q, double n, double gamma, Vectorofdouble& wavespeed, int D);
    void ywavespeedestimate(Vectorofdouble p, Vectorofdouble q, double n, double gamma, Vectorofdouble& wavespeed, int D);
    void compute_nplushalf_variables(Vectorofdouble a, Vectorofdouble b, Vectorofdouble c, 
    double dx, double dt, double gamma, Vectorofdouble& uiL_nplushalf, Vectorofdouble& uiR_nplushalf, int nVar, int D);
    void compute_HLLCflux(Vectorofdouble a, Vectorofdouble p, Vectorofdouble b, 
    Vectorofdouble q, Vectorofdouble wavespeed, Vectorofdouble& Fhllc, int nVar, double gamma, int D, double x, double t);

    void ycompute_nplushalf_variables(Vectorofdouble a, Vectorofdouble b, Vectorofdouble c, 
    double dx, double dt, double gamma, Vectorofdouble& uiL_nplushalf, Vectorofdouble& uiR_nplushalf, int nVar, int D);
    void ycompute_HLLCflux(Vectorofdouble a, Vectorofdouble p, Vectorofdouble b, 
    Vectorofdouble q, Vectorofdouble wavespeed, Vectorofdouble& Fhllc, int nVar, double gamma, int D, double x, double t);
    Vectorofdouble getFORCEflux_bar(Vectorofdouble a, Vectorofdouble b,Vectorofdouble p,Vectorofdouble q, 
    double dx, double dt, double gamma, double nVar, int D);
    Vectorofdouble ygetFORCEflux_bar(Vectorofdouble a, Vectorofdouble b,Vectorofdouble p,Vectorofdouble q, 
    double dx, double dt, double gamma, double nVar, int D);
};


