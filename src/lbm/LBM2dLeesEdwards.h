#ifndef _LBM_2D_LEES_EDWARDS_H_
#define _LBM_2D_LEES_EDWARDS_H_

#include "LBM2dSimulator.h"

class LBM2dLeesEdwards:public LBM2dSimulator{
public:
	LBM2dLeesEdwards(LBM2dLattice* lattice,
                  double viscosity,double density,double dt,
                  int* step,
                  double shear_rate=0.0);
	virtual ~LBM2dLeesEdwards();
	virtual void collide();
protected:
	virtual void collide_up_down_periodic(int i);
	virtual void collide_sides_periodic(int j);
	virtual void collide_corner();
	//
	virtual void periodic_down(int n,int i);//n=2,5,6
	virtual void periodic_up(int n,int i);//n=4,7,8
	//virtual void set_boundary_distribution();
	virtual void add_boundary_force();
	virtual void set_boundary_viscosity();
private:
	int* step;
	double gdot;
	//
	double Du;
	double a0,a1;
	int Dn;
};


#endif // _LBM_2D_LEES_EDWARDS_H_
