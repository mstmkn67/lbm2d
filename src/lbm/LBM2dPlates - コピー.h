#ifndef _LBM_2D_PLATES_H_
#define _LBM_2D_PLATES_H_

#include "LBM2dSimulator.h"

class LBM2dPlates:public LBM2dSimulator{
public:
	LBM2dPlates(LBM2dLattice* lattice,
	           double viscosity,double density,double dt,
	           double delta_p=0.0,
	           double u_down=0.0,double v_down=0.0,
	           double u_up=0.0,double v_up=0.0);
	virtual ~LBM2dPlates();
	virtual void collide();
protected:
	virtual void collide_down_wall(int i);
	virtual void collide_up_wall(int i);
	virtual void collide_sides_periodic(int j);
	virtual void collide_corner0();
	virtual void collide_corner1();
	virtual void collide_corner2();
	//virtual void set_boundary_distribution();
	virtual void add_boundary_force();
	virtual void set_boundary_viscosity();
private:
	double u0,v0;
	double u2,v2;
	double delp;
	//
	//double rhoo,rhoi;/////////////
};


#endif // _LBM_2D_PLATES_H_
