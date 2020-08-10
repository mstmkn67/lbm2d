#include "LBM2dPlates.h"

LBM2dPlates::LBM2dPlates(
	LBM2dLattice* lattice,
	double viscosity,double density,double dt,
	double delta_p,
	double u_down,double v_down,double u_up,double v_up)
	:LBM2dSimulator(lattice,viscosity,density,dt,delta_p),
	 u0(u_down),v0(v_down),u2(u_up),v2(v_up){
}

LBM2dPlates::~LBM2dPlates(){}

void LBM2dPlates::collide(){
	for(int j=1;j<lattice->ny;j++){
		for(int i=1;i<lattice->nx;i++){
			for(int n=0;n<9;n++){
				calc_collide(n,i,j);
			}
		}
	}
	for(int i=1;i<lattice->nx;i++){
		collide_down_wall(i);//down
		collide_up_wall(i);//up
	}
	for(int j=1;j<lattice->ny;j++){
		collide_sides_periodic(j);
	}
	collide_corner();
}

void LBM2dPlates::collide_down_wall(int i){
	int j=0;
	//0,1,3,4,7,8
	calc_collide(0,i,j);
	calc_collide(1,i,j);
	calc_collide(3,i,j);
	calc_collide(4,i,j);
	calc_collide(7,i,j);
	calc_collide(8,i,j);
	//2,5,6
	nonslip_down_wall(i);
}

void LBM2dPlates::nonslip_down_wall(int i){
	int j=0;
	double rhow=lattice->f[0][i][j]+lattice->f[1][i][j]+lattice->f[3][i][j]
		         +2.0*(lattice->f[4][i][j]+lattice->f[7][i][j]+lattice->f[8][i][j]);
	rhow/=1.0-v0/c;
	lattice->f[2][i][j]=lattice->f[4][i][j]+2./3*rhow*v0/c;
	lattice->f[5][i][j]=lattice->f[7][i][j]+0.5*(lattice->f[3][i][j]-lattice->f[1][i][j])
	                   +1./6*rhow*v0/c+0.5*rhow*u0/c;
	lattice->f[6][i][j]=lattice->f[8][i][j]+0.5*(lattice->f[1][i][j]-lattice->f[3][i][j])
	                   +1./6*rhow*v0/c-0.5*rhow*u0/c;
}

void LBM2dPlates::collide_up_wall(int i){
	int j=lattice->ny;
	//0,1,2,3,5,6
	calc_collide(0,i,j);
	calc_collide(1,i,j);
	calc_collide(2,i,j);
	calc_collide(3,i,j);
	calc_collide(5,i,j);
	calc_collide(6,i,j);
	//4,7,8
	nonslip_up_wall(i);
}

void LBM2dPlates::nonslip_up_wall(int i){
	int j=lattice->ny;
	double rhow=lattice->f[0][i][j]+lattice->f[1][i][j]+lattice->f[3][i][j]
		         +2.0*(lattice->f[2][i][j]+lattice->f[5][i][j]+lattice->f[6][i][j]);
	rhow/=1.0+v2/c;
	lattice->f[4][i][j]=lattice->f[2][i][j]-2./3*rhow*v2/c;
	lattice->f[7][i][j]=lattice->f[5][i][j]+0.5*(lattice->f[1][i][j]-lattice->f[3][i][j])
	                   -1./6*rhow*v2/c-0.5*rhow*u2/c;
	lattice->f[8][i][j]=lattice->f[6][i][j]+0.5*(lattice->f[3][i][j]-lattice->f[1][i][j])
	                   -1./6*rhow*v2/c+0.5*rhow*u2/c;
}

void LBM2dPlates::collide_sides_periodic(int j){
	int i0=0,i1=lattice->nx;
	//left side
	calc_collide(0,i0,j);
	calc_collide(2,i0,j);
	calc_collide(3,i0,j);
	calc_collide(4,i0,j);
	calc_collide(6,i0,j);
	calc_collide(7,i0,j);
	//right side
	calc_collide(0,i1,j);
	calc_collide(1,i1,j);
	calc_collide(2,i1,j);
	calc_collide(4,i1,j);
	calc_collide(5,i1,j);
	calc_collide(8,i1,j);
	//
	//left side
	lattice->f[1][i0][j]=lattice->f[1][i1][j];
	lattice->f[5][i0][j]=lattice->f[5][i1][j];
	lattice->f[8][i0][j]=lattice->f[8][i1][j];
	//right side
	lattice->f[3][i1][j]=lattice->f[3][i0][j];
	lattice->f[6][i1][j]=lattice->f[6][i0][j];
	lattice->f[7][i1][j]=lattice->f[7][i0][j];
}

void LBM2dPlates::collide_corner(){
	////standard collisions
	//left down
	int i=0,j=0;
	calc_collide(0,i,j);
	calc_collide(3,i,j);
	calc_collide(4,i,j);
	calc_collide(7,i,j);
	//right down
	i=lattice->nx;j=0;
	calc_collide(0,i,j);
	calc_collide(1,i,j);
	calc_collide(4,i,j);
	calc_collide(8,i,j);
	//left up
	i=0;j=lattice->ny;
	calc_collide(0,i,j);
	calc_collide(2,i,j);
	calc_collide(3,i,j);
	calc_collide(6,i,j);	
	//right up
	i=lattice->nx;j=lattice->ny;
	calc_collide(0,i,j);
	calc_collide(1,i,j);
	calc_collide(2,i,j);
	calc_collide(5,i,j);
	////periodic boundary conditions
	//down
	lattice->f[1][0          ][0]=lattice->f[1][lattice->nx][0];
	lattice->f[8][0          ][0]=lattice->f[8][lattice->nx][0];
	lattice->f[3][lattice->nx][0]=lattice->f[3][0          ][0];
	lattice->f[7][lattice->nx][0]=lattice->f[7][0          ][0];
	//up
	lattice->f[1][0          ][lattice->ny]=lattice->f[1][lattice->nx][lattice->ny];
	lattice->f[5][0          ][lattice->ny]=lattice->f[5][lattice->nx][lattice->ny];
	lattice->f[3][lattice->nx][lattice->ny]=lattice->f[3][0          ][lattice->ny];
	lattice->f[6][lattice->nx][lattice->ny]=lattice->f[6][0          ][lattice->ny];
	////non slip conditions
	//2,5,6
	nonslip_down_wall(0);
	nonslip_down_wall(lattice->nx);
	//4,7,8
	nonslip_up_wall(0);
	nonslip_up_wall(lattice->nx);
}

void LBM2dPlates::add_boundary_force(){
}

void LBM2dPlates::set_boundary_viscosity(){
}
