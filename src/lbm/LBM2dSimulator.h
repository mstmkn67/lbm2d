#ifndef _LBM_2D_SIMULATOR_H_
#define _LBM_2D_SIMULATOR_H_

#include "LBM2dLattice.h"
#include <cmath>
#include <stdlib.h>
using namespace std;

class LBM2dSimulator{
public:
	LBM2dSimulator(LBM2dLattice* lattice,
	               double viscosity,double density,double dt,double delp=0.0);
	virtual ~LBM2dSimulator();
	virtual void initial();
	virtual void update();
	//
	virtual void set_bulk_viscosity();
	virtual void reset_body_force();
	//
	LBM2dLattice* lattice;
	//
	double delp;
protected:
	virtual void collide()=0;
	virtual void calc_density_velocity();
	virtual void calc_collide(int n,int i,int j);
	virtual void time_update();
	//virtual void set_boundary_distribution()=0;
	virtual void add_boundary_force()=0;
	virtual void set_boundary_viscosity()=0;
	//
	double viscosity;
	double density;
	double dt;
	//
	double c,cs2;
	double w[9];
	int cx[9];
	int cy[9];
};

#endif // _LBM_2D_SIMULATOR_H_
