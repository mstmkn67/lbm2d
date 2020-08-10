#include "LBM2dLeesEdwards.h"

LBM2dLeesEdwards::LBM2dLeesEdwards(
	LBM2dLattice* lattice,
	double viscosity,double density,double dt,
	int* st,double sr)
	:LBM2dSimulator(lattice,viscosity,density,dt),
	 step(st),gdot(sr){
	Du=gdot*lattice->ly;
}

LBM2dLeesEdwards::~LBM2dLeesEdwards(){}

void LBM2dLeesEdwards::collide(){
	double Dx=(*step)*dt*gdot*lattice->ly;
	Dx-=floor(Dx/lattice->lx)*lattice->lx;
	Dn=(int)floor(Dx/lattice->dx);
	double DDx=Dx-Dn*lattice->dx;
	a0=DDx/lattice->dx;a1=1.0-a0;
	//
	for(int j=1;j<lattice->ny;j++){
		for(int i=1;i<lattice->nx;i++){
			for(int n=0;n<9;n++){
				calc_collide(n,i,j);
			}
		}
	}
	for(int i=1;i<lattice->nx;i++){
		collide_up_down_periodic(i);
	}
	for(int j=1;j<lattice->ny;j++){
		collide_sides_periodic(j);
	}
	collide_corner();
}

void LBM2dLeesEdwards::collide_up_down_periodic(int i){
	//// standard collisions of down side
	//0,1,3,4,7,8
	calc_collide(0,i,0);
	calc_collide(1,i,0);
	calc_collide(3,i,0);
	calc_collide(4,i,0);
	calc_collide(7,i,0);
	calc_collide(8,i,0);
	//// standard collisions of up side
	//0,1,2,3,5,6
	calc_collide(0,i,lattice->ny);
	calc_collide(1,i,lattice->ny);
	calc_collide(2,i,lattice->ny);
	calc_collide(3,i,lattice->ny);
	calc_collide(5,i,lattice->ny);
	calc_collide(6,i,lattice->ny);
	//// Lees Edwards condisions of down side
	//2,5,6
	periodic_down(2,i);
	periodic_down(5,i);
	periodic_down(6,i);
	//// Lees Edwards condisions of up side
	//4,7,8
	periodic_up(4,i);
	periodic_up(7,i);
	periodic_up(8,i);
}

void LBM2dLeesEdwards::periodic_down(int n,int i){
	int i0=i+Dn,i1=i+Dn+1;
	if(i0<0)i0+=lattice->nx;          if(i1<0)i1+=lattice->nx;
	if(i0>lattice->nx)i0-=lattice->nx;if(i1>lattice->nx)i1-=lattice->nx;
	double u=a1*lattice->u[i0][lattice->ny]+a0*lattice->u[i1][lattice->ny];
	double v=a1*lattice->v[i0][lattice->ny]+a0*lattice->v[i1][lattice->ny];
	double vis=a1*lattice->vis[i0][lattice->ny]+a0*lattice->vis[i1][lattice->ny];
	double rho=a1*lattice->rho[i0][lattice->ny]+a0*lattice->rho[i1][lattice->ny];
	double bfx=a1*lattice->bfx[i0][lattice->ny]+a0*lattice->bfx[i1][lattice->ny];
	double bfy=a1*lattice->bfy[i0][lattice->ny]+a0*lattice->bfy[i1][lattice->ny];
	double tau=vis/rho/cs2/dt+0.5;
	double ueq=u+tau*bfx*dt/rho;
	double ueq2=u-Du+tau*bfx*dt/rho;
	double veq=v+tau*bfy*dt/rho;
	double uc=ueq*(cx[n]*c)+veq*(cy[n]*c);
	double uc2=ueq2*(cx[n]*c)+veq*(cy[n]*c);
	double uu=ueq*ueq+veq*veq;
	double uu2=ueq2*ueq2+veq*veq;
	double feq=rho*w[n]*(1.0+uc/cs2+0.5*uc*uc/cs2/cs2-0.5*uu/cs2);
	double feq2=rho*w[n]*(1.0+uc2/cs2+0.5*uc2*uc2/cs2/cs2-0.5*uu2/cs2);
	lattice->f[n][i][0]=a1*lattice->f[n][i0][lattice->ny]+a0*lattice->f[n][i1][lattice->ny]+feq2-feq;
}

void LBM2dLeesEdwards::periodic_up(int n,int i){
	int i0=i-Dn-1,i1=i-Dn;
	if(i0<0)i0+=lattice->nx;          if(i1<0)i1+=lattice->nx;
	if(i0>lattice->nx)i0-=lattice->nx;if(i1>lattice->nx)i1-=lattice->nx;
	double u=a0*lattice->u[i0][0]+a1*lattice->u[i1][0];
	double v=a0*lattice->v[i0][0]+a1*lattice->v[i1][0];
	double vis=a0*lattice->vis[i0][0]+a1*lattice->vis[i1][0];
	double rho=a0*lattice->rho[i0][0]+a1*lattice->rho[i1][0];
	double bfx=a0*lattice->bfx[i0][0]+a1*lattice->bfx[i1][0];
	double bfy=a0*lattice->bfy[i0][0]+a1*lattice->bfy[i1][0];
	double tau=vis/rho/cs2/dt+0.5;
	double ueq=u+tau*bfx*dt/rho;
	double ueq2=u+Du+tau*bfx*dt/rho;
	double veq=v+tau*bfy*dt/rho;
	double uc=ueq*(cx[n]*c)+veq*(cy[n]*c);
	double uc2=ueq2*(cx[n]*c)+veq*(cy[n]*c);
	double uu=ueq*ueq+veq*veq;
	double uu2=ueq2*ueq2+veq*veq;
	double feq=rho*w[n]*(1.0+uc/cs2+0.5*uc*uc/cs2/cs2-0.5*uu/cs2);
	double feq2=rho*w[n]*(1.0+uc2/cs2+0.5*uc2*uc2/cs2/cs2-0.5*uu2/cs2);
	lattice->f[n][i][lattice->ny]=a0*lattice->f[n][i0][0]+a1*lattice->f[n][i1][0]+feq2-feq;
}

void LBM2dLeesEdwards::collide_sides_periodic(int j){
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

void LBM2dLeesEdwards::collide_corner(){
	////standard collisions
	//down left
	calc_collide(0,0,0);
	calc_collide(3,0,0);
	calc_collide(7,0,0);
	calc_collide(4,0,0);
	//down right
	calc_collide(0,lattice->nx,0);
	calc_collide(1,lattice->nx,0);
	calc_collide(4,lattice->nx,0);
	calc_collide(8,lattice->nx,0);
	//up left
	calc_collide(0,0,lattice->ny);
	calc_collide(2,0,lattice->ny);
	calc_collide(3,0,lattice->ny);
	calc_collide(6,0,lattice->ny);
	//up right
	calc_collide(0,lattice->nx,lattice->ny);
	calc_collide(1,lattice->nx,lattice->ny);
	calc_collide(2,lattice->nx,lattice->ny);
	calc_collide(5,lattice->nx,lattice->ny);
	////periodic condition
	//down left
	lattice->f[1][0][0]=lattice->f[1][lattice->nx][0];
	lattice->f[8][0][0]=lattice->f[8][lattice->nx][0];
	//down right
	lattice->f[3][lattice->nx][0]=lattice->f[3][0][0];
	lattice->f[7][lattice->nx][0]=lattice->f[7][0][0];
	//up left
	lattice->f[1][0][lattice->ny]=lattice->f[1][lattice->nx][lattice->ny];
	lattice->f[5][0][lattice->ny]=lattice->f[5][lattice->nx][lattice->ny];
	//up right
	lattice->f[3][lattice->nx][lattice->ny]=lattice->f[3][0][lattice->ny];
	lattice->f[6][lattice->nx][lattice->ny]=lattice->f[6][0][lattice->ny];
	////Lees Edwards condition
	//down left
	periodic_down(2,0);
	periodic_down(5,0);
	periodic_down(6,0);
	//down right
	periodic_down(2,lattice->nx);
	periodic_down(5,lattice->nx);
	periodic_down(6,lattice->nx);
	//up left
	periodic_up(4,0);
	periodic_up(7,0);
	periodic_up(8,0);
	//up right
	periodic_up(4,lattice->ny);
	periodic_up(7,lattice->ny);
	periodic_up(8,lattice->ny);
}

void LBM2dLeesEdwards::add_boundary_force(){
}

void LBM2dLeesEdwards::set_boundary_viscosity(){
}
