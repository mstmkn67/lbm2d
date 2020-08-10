#ifndef _LBM_2D_LATTICE_H_
#define _LBM_2D_LATTICE_H_

#include "../SimpleArray.h"
#include <iostream>
using namespace std;

class LBM2dLattice{
public:
	LBM2dLattice(int nx,int ny,double dx);
	virtual ~LBM2dLattice();
	int nx,ny;
	double lx,ly;
	double dx;
	Array2d<double> u,v;
	Array2d<double> f[9];//distribution function
	Array2d<double> f0[9];//distribution function previous time
	Array2d<double> bfx,bfy;//body force
	Array2d<double> rho;//density
	Array2d<double> vis;//viscosity field
private:
};

#endif // _LBM_2D_LATTICE_H_
