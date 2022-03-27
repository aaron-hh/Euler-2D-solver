/* -------------------------------------------------------------------*/
/*                                                                    */
/*  File Description : This file contains the necessary functions for */
/*        the data structure used in the simulation.                  */
/*                                                                    */
/*  Programer: Aaron Ng                                               */
/*                                                                    */
/*  Last Revision: 01 February 2022                                   */
/*                                                                    */
/* -------------------------------------------------------------------*/

#include"array.H"

m_array::m_array()
{}

void m_array::setSize(double XSIZE, double YSIZE, double ZSIZE, double n_Var)
{
    m_xSize = XSIZE;
    m_ySize = YSIZE;
    m_zSize = ZSIZE;
    m_TSize = XSIZE * YSIZE * ZSIZE;
    m_data.resize(m_TSize, std::vector<double> (n_Var));
}

std::vector<double>& m_array::operator()(int i, int j, int k) 
{
    return m_data[i + j*m_xSize + k*m_xSize*m_ySize];
}

double& m_array::operator()(int i, int j, int k, int l) 
{
    return m_data[i + j*m_xSize + k*m_xSize*m_ySize][l];
}