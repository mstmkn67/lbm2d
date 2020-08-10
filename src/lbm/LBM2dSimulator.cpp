#include "LBM2dSimulator.h"

#define PI 3.1415926535

LBM2dSimulator::LBM2dSimulator(LBM2dLattice* lat,double vis,double den,double t,double dp)
:lattice(lat),delp(dp),viscosity(vis),density(den),dt(t){
	c=lattice->dx/dt;
	cs2=c*c/3.0;
	cx[0]= 0;cy[0]= 0;w[0]=4.0/9.0;
	cx[1]= 1;cy[1]= 0;w[1]=1.0/9.0;
	cx[2]= 0;cy[2]= 1;w[2]=1.0/9.0;
	cx[3]=-1;cy[3]= 0;w[3]=1.0/9.0;
	cx[4]= 0;cy[4]=-1;w[4]=1.0/9.0;
	cx[5]= 1;cy[5]= 1;w[5]=1.0/36.0;
	cx[6]=-1;cy[6]= 1;w[6]=1.0/36.0;
	cx[7]=-1;cy[7]=-1;w[7]=1.0/36.0;
	cx[8]= 1;cy[8]=-1;w[8]=1.0/36.0;
}

LBM2dSimulator::~LBM2dSimulator(){
}

void LBM2dSimulator::initial(){
	for(int j=0;j<=lattice->ny;j++){
		for(int i=0;i<=lattice->nx;i++){
			for(int n=0;n<9;n++){
				lattice->f[n][i][j]=density*w[n];
				lattice->f0[n][i][j]=density*w[n];
			}
		}
	}
	reset_body_force();
	calc_density_velocity();
}

void LBM2dSimulator::update(){
	//set_boundary_distribution();
	add_boundary_force();
	set_boundary_viscosity();
	collide();
	calc_density_velocity();
	time_update();
}

void LBM2dSimulator::set_bulk_viscosity(){
	for(int j=0;j<=lattice->ny;j++){
		for(int i=0;i<=lattice->nx;i++){
			lattice->vis[i][j]=viscosity;
			//test
			//if(j>=lattice->ny/2){
			//	lattice->vis[i][j]*=2;
			//}
			//
		}
	}
}

void LBM2dSimulator::reset_body_force(){
	for(int j=0;j<=lattice->ny;j++){
		for(int i=0;i<=lattice->nx;i++){
			lattice->bfx[i][j]=delp/lattice->lx;
			lattice->bfy[i][j]=0.0;
		}
	}
}

void LBM2dSimulator::calc_density_velocity(){
	for(int j=0;j<=lattice->ny;j++){
		for(int i=0;i<=lattice->nx;i++){
			lattice->rho[i][j]=lattice->f[0][i][j];
			for(int n=1;n<9;n++){
				lattice->rho[i][j]+=lattice->f[n][i][j];
			}
		}
	}
	for(int j=0;j<=lattice->ny;j++){
		for(int i=0;i<=lattice->nx;i++){
			lattice->u[i][j]=lattice->f[0][i][j]*(cx[0]*c);
			lattice->v[i][j]=lattice->f[0][i][j]*(cy[0]*c);
			for(int n=1;n<9;n++){
				lattice->u[i][j]+=lattice->f[n][i][j]*(cx[n]*c);
				lattice->v[i][j]+=lattice->f[n][i][j]*(cy[n]*c);
			}
			lattice->u[i][j]/=lattice->rho[i][j];
			lattice->v[i][j]/=lattice->rho[i][j];
		}
	}
}

void LBM2dSimulator::calc_collide(int n,int i,int j){
	int ii=i-cx[n],jj=j-cy[n];
	double tau=lattice->vis[ii][jj]/lattice->rho[ii][jj]/cs2/dt+0.5;
	double ueq=lattice->u[ii][jj]+tau*lattice->bfx[ii][jj]*dt/lattice->rho[ii][jj];
	double veq=lattice->v[ii][jj]+tau*lattice->bfy[ii][jj]*dt/lattice->rho[ii][jj];
	double uc=ueq*(cx[n]*c)+veq*(cy[n]*c);
	double uu=ueq*ueq+veq*veq;
	double feq=lattice->rho[ii][jj]*w[n]*(1.0+uc/cs2+0.5*uc*uc/cs2/cs2-0.5*uu/cs2);
	double a=lattice->f0[n][ii][jj]-(lattice->f0[n][ii][jj]-feq)/tau;
	//a+=w[n]*dt/cs2*(lattice->bfx[ii][jj]*cx[n]*c+lattice->bfy[ii][jj]*cy[n]*c);
	lattice->f[n][i][j]=a;
}

void LBM2dSimulator::time_update(){
	for(int j=0;j<=lattice->ny;j++){
		for(int i=0;i<=lattice->nx;i++){
			for(int n=0;n<9;n++){
				lattice->f0[n][i][j]=lattice->f[n][i][j];
			}
		}
	}
}

#undef PI
