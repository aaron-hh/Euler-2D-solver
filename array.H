#include <cmath>
#include <vector>
#include <algorithm>
#include <fstream>

class m_array
{
private:
    size_t m_xSize;
    size_t m_ySize;
    size_t m_zSize;
    size_t m_TSize;

    std::vector<std::vector<double> >m_data; //2 vectors as one for index and another for data“

public:
    m_array();

    void setSize(double XSIZE, double YSIZE, double ZSIZE, double n_Var);
    std::vector<double>& operator()(int i, int j, int k);
    double& operator()(int i, int j, int k, int l);
};


