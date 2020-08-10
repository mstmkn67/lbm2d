#include "LBM2dPlates.h"

LBM2dPlates::LBM2dPlates(
	LBM2dLattice* lattice,
	double viscosity,double density,double dt,
	double delta_p,
	double u_down,double v_down,double u_up,double v_up)
	:LBM2dSimulator(lattice,viscosity,density,dt),
	 u0(u_down),v0(v_down),u2(u_up),v2(v_up),delp(delta_p){
}

LBM2dPlates::~LBM2dPlates(){}

void LBM2dPlates::collide(){
	for(int j=1;j<=lattice->ny;j++){
		for(int i=1;i<=lattice->nx;i++){
			for(int n=0;n<9;n++){
				int ii=i-cx[n],jj=j-cy[n];
				lattice->f[n][i][j]=lattice->f0[n][ii][jj]+calc_rhs_collide(n,i,j);
			}
		}
	}
	for(int i=1;i<=lattice->nx;i++){
		collide_down_wall(i);//down
		collide_up_wall(i);//up
	}
	////////////////////////////////
	/*rhoo=0.0;rhoi=0.0;
	for(int j=0;j<=lattice->ny+1;j++){
		rhoi+=lattice->rho[0][j];
		rhoo+=lattice->rho[lattice->nx+1][j];
	}
	rhoi/=lattice->ny+2;rhoo/=lattice->ny+2;*/
	////////////////////////////////////////
	for(int j=1;j<=lattice->ny;j++){
		collide_sides_periodic(j);
	}
	collide_corner0();
	collide_corner1();
	collide_corner2();
}

void LBM2dPlates::collide_down_wall(int i){
	int j=0;
	//0,1,3,4,7,8
	int ii=i-cx[0],jj=j-cy[0];
	lattice->f[0][i][j]=lattice->f0[0][ii][jj]+calc_rhs_collide(0,i,j);
	ii=i-cx[1];jj=j-cy[1];
	lattice->f[1][i][j]=lattice->f0[1][ii][jj]+calc_rhs_collide(1,i,j);
	ii=i-cx[3];jj=j-cy[3];
	lattice->f[3][i][j]=lattice->f0[3][ii][jj]+calc_rhs_collide(3,i,j);
	ii=i-cx[4];jj=j-cy[4];
	lattice->f[4][i][j]=lattice->f0[4][ii][jj]+calc_rhs_collide(4,i,j);
	ii=i-cx[7];jj=j-cy[7];
	lattice->f[7][i][j]=lattice->f0[7][ii][jj]+calc_rhs_collide(7,i,j);
	ii=i-cx[8];jj=j-cy[8];
	lattice->f[8][i][j]=lattice->f0[8][ii][jj]+calc_rhs_collide(8,i,j);
	//2,5,6
	double rhow=lattice->f[0][i][j]+lattice->f[1][i][j]+lattice->f[3][i][j]
		         +2.0*(lattice->f[4][i][j]+lattice->f[7][i][j]+lattice->f[8][i][j]);
	rhow/=1.0-v0;
	double rhop=rhow*v0+lattice->f[4][i][j]+lattice->f[7][i][j]+lattice->f[8][i][j];
	rhop*=6./(1.+3.*v0+3.*v0*v0);
	double upx=rhow*u0-(lattice->f[1][i][j]-lattice->f[3][i][j]
	                   +lattice->f[8][i][j]-lattice->f[7][i][j]);
	upx*=6./rhop;
	upx-=u0+3.*u0*v0;
	upx/=1.+3.*v0;
	lattice->f[2][i][j]=rhop/9.*(1.+3.*v0+4.5*v0*v0
	                            -1.5*((u0+upx)*(u0+upx)+v0*v0));
	lattice->f[5][i][j]=rhop/36.*(1.+3.*(u0+upx+v0)
	                               +4.5*(u0+upx+v0)*(u0+upx+v0)
	                               -1.5*((u0+upx)*(u0+upx)+v0*v0));
	lattice->f[6][i][j]=rhop/36.*(1.+3.*(-u0-upx+v0)
	                               +4.5*(-u0-upx+v0)*(-u0-upx+v0)
	                               -1.5*((u0+upx)*(u0+upx)+v0*v0));
	/////////////////
	//lattice->f[2][i][j]=lattice->f[4][i][j];
	//lattice->f[5][i][j]=lattice->f[7][i][j];
	//lattice->f[6][i][j]=lattice->f[8][i][j];
	/////////////////
}

void LBM2dPlates::collide_up_wall(int i){
	int j=lattice->ny+1;
	//0,1,2,3,5,6
	int ii=i-cx[0],jj=j-cy[0];
	lattice->f[0][i][j]=lattice->f0[0][ii][jj]+calc_rhs_collide(0,i,j);
	ii=i-cx[1];jj=j-cy[1];
	lattice->f[1][i][j]=lattice->f0[1][ii][jj]+calc_rhs_collide(1,i,j);
	ii=i-cx[2];jj=j-cy[2];
	lattice->f[2][i][j]=lattice->f0[2][ii][jj]+calc_rhs_collide(2,i,j);
	ii=i-cx[3];jj=j-cy[3];
	lattice->f[3][i][j]=lattice->f0[3][ii][jj]+calc_rhs_collide(3,i,j);
	ii=i-cx[5];jj=j-cy[5];
	lattice->f[5][i][j]=lattice->f0[5][ii][jj]+calc_rhs_collide(5,i,j);
	ii=i-cx[6];jj=j-cy[6];
	lattice->f[6][i][j]=lattice->f0[6][ii][jj]+calc_rhs_collide(6,i,j);
	//4,7,8
	double rhow=lattice->f[0][i][j]+lattice->f[1][i][j]+lattice->f[3][i][j]
		         +2.0*(lattice->f[2][i][j]+lattice->f[5][i][j]+lattice->f[6][i][j]);
	rhow/=1.0-v2;
	double rhop=rhow*v2+lattice->f[2][i][j]+lattice->f[5][i][j]+lattice->f[6][i][j];
	rhop*=6./(1.+3.*v2+3.*v2*v2);
	double upx=rhow*u2-(lattice->f[1][i][j]-lattice->f[3][i][j]
	                   +lattice->f[5][i][j]-lattice->f[6][i][j]);
	upx*=6./rhop;
	upx-=u2+3.*u2*v2;
	upx/=1.+3.*v2;
	lattice->f[4][i][j]=rhop/9.*(1.+3.*v2+4.5*v2*v2
	                            -1.5*((u2+upx)*(u2+upx)+v2*v2));
	lattice->f[8][i][j]=rhop/36.*(1.+3.*(u2+upx+v2)
	                               +4.5*(u2+upx+v2)*(u2+upx+v2)
	                               -1.5*((u2+upx)*(u2+upx)+v2*v2));
	lattice->f[7][i][j]=rhop/36.*(1.+3.*(-u2-upx+v2)
	                               +4.5*(-u2-upx+v2)*(-u2-upx+v2)
	                               -1.5*((u2+upx)*(u2+upx)+v2*v2));
	/////////////////
	//lattice->f[4][i][j]=lattice->f[2][i][j];
	//lattice->f[8][i][j]=lattice->f[6][i][j];
	//lattice->f[7][i][j]=lattice->f[5][i][j];
	/////////////////
}

void LBM2dPlates::collide_sides_periodic(int j){
	int i0=0,i1=lattice->nx+1;
	//left side
	int ii=i0-cx[0],jj=j-cy[0];
	lattice->f[0][i0][j]=lattice->f0[0][ii][jj]+calc_rhs_collide(0,i0,j);
	ii=i0-cx[2];jj=j-cy[2];
	lattice->f[2][i0][j]=lattice->f0[2][ii][jj]+calc_rhs_collide(2,i0,j);
	ii=i0-cx[3];jj=j-cy[3];
	lattice->f[3][i0][j]=lattice->f0[3][ii][jj]+calc_rhs_collide(3,i0,j);
	ii=i0-cx[4];jj=j-cy[4];
	lattice->f[4][i0][j]=lattice->f0[4][ii][jj]+calc_rhs_collide(4,i0,j);
	ii=i0-cx[6];jj=j-cy[6];
	lattice->f[6][i0][j]=lattice->f0[6][ii][jj]+calc_rhs_collide(6,i0,j);
	ii=i0-cx[7];jj=j-cy[7];
	lattice->f[7][i0][j]=lattice->f0[7][ii][jj]+calc_rhs_collide(7,i0,j);
	//right side
	ii=i1-cx[0];jj=j-cy[0];
	lattice->f[0][i1][j]=lattice->f0[0][ii][jj]+calc_rhs_collide(0,i1,j);
	ii=i1-cx[1];jj=j-cy[1];
	lattice->f[1][i1][j]=lattice->f0[1][ii][jj]+calc_rhs_collide(1,i1,j);
	ii=i1-cx[2];jj=j-cy[2];
	lattice->f[2][i1][j]=lattice->f0[2][ii][jj]+calc_rhs_collide(2,i1,j);
	ii=i1-cx[4];jj=j-cy[4];
	lattice->f[4][i1][j]=lattice->f0[4][ii][jj]+calc_rhs_collide(4,i1,j);
	ii=i1-cx[5];jj=j-cy[5];
	lattice->f[5][i1][j]=lattice->f0[5][ii][jj]+calc_rhs_collide(5,i1,j);
	ii=i1-cx[8];jj=j-cy[8];
	lattice->f[8][i1][j]=lattice->f0[8][ii][jj]+calc_rhs_collide(8,i1,j);
	//
	double c=delp-(lattice->f[0][i0][j]-lattice->f[0][i1][j]
	              +lattice->f[2][i0][j]-lattice->f[2][i1][j]
	              +lattice->f[4][i0][j]-lattice->f[4][i1][j])/3.0;
	//left side
	lattice->f[1][i0][j]=lattice->f[1][i1][j]+c;
	lattice->f[5][i0][j]=lattice->f[5][i1][j]+0.25*c;
	lattice->f[8][i0][j]=lattice->f[8][i1][j]+0.25*c;
	//right side
	lattice->f[3][i1][j]=lattice->f[3][i0][j]-c;
	lattice->f[6][i1][j]=lattice->f[6][i0][j]-0.25*c;
	lattice->f[7][i1][j]=lattice->f[7][i0][j]-0.25*c;
		/////////////////////////////////////
	//left side
	/*
	lattice->f[1][i0][j]=lattice->f[1][i1][j]*(density+delp/cs2)/rhoo;
	lattice->f[5][i0][j]=lattice->f[5][i1][j]*(density+delp/cs2)/rhoo;
	lattice->f[8][i0][j]=lattice->f[8][i1][j]*(density+delp/cs2)/rhoo;
	//right side
	lattice->f[3][i1][j]=lattice->f[3][i0][j]*(density-delp/cs2)/rhoi;
	lattice->f[6][i1][j]=lattice->f[6][i0][j]*(density-delp/cs2)/rhoi;
	lattice->f[7][i1][j]=lattice->f[7][i0][j]*(density-delp/cs2)/rhoi;*/
	////////////////////////////////////////

}

void LBM2dPlates::collide_corner0(){
	//left down
	int i=0,j=0;
	int ii=i-cx[0],jj=j-cy[0];
	lattice->f[0][i][j]=lattice->f0[0][ii][jj]+calc_rhs_collide(0,i,j);
	ii=i-cx[3];jj=j-cy[3];
	lattice->f[3][i][j]=lattice->f0[3][ii][jj]+calc_rhs_collide(3,i,j);
	ii=i-cx[4];jj=j-cy[4];
	lattice->f[4][i][j]=lattice->f0[4][ii][jj]+calc_rhs_collide(4,i,j);
	ii=i-cx[7];jj=j-cy[7];
	lattice->f[7][i][j]=lattice->f0[7][ii][jj]+calc_rhs_collide(7,i,j);
	//right down
	i=lattice->nx+1;j=0;
	ii=i-cx[0];jj=j-cy[0];
	lattice->f[0][i][j]=lattice->f0[0][ii][jj]+calc_rhs_collide(0,i,j);
	ii=i-cx[1];jj=j-cy[1];
	lattice->f[1][i][j]=lattice->f0[1][ii][jj]+calc_rhs_collide(1,i,j);
	ii=i-cx[4];jj=j-cy[4];
	lattice->f[4][i][j]=lattice->f0[4][ii][jj]+calc_rhs_collide(4,i,j);
	ii=i-cx[8];jj=j-cy[8];
	lattice->f[8][i][j]=lattice->f0[8][ii][jj]+calc_rhs_collide(8,i,j);
	
	//left up
	i=0;j=lattice->ny+1;
	ii=i-cx[0];jj=j-cy[0];
	lattice->f[0][i][j]=lattice->f0[0][ii][jj]+calc_rhs_collide(0,i,j);
	ii=i-cx[2];jj=j-cy[2];
	lattice->f[2][i][j]=lattice->f0[2][ii][jj]+calc_rhs_collide(2,i,j);
	ii=i-cx[3];jj=j-cy[3];
	lattice->f[3][i][j]=lattice->f0[3][ii][jj]+calc_rhs_collide(3,i,j);
	ii=i-cx[6];jj=j-cy[6];
	lattice->f[6][i][j]=lattice->f0[6][ii][jj]+calc_rhs_collide(6,i,j);	
	//right up
	i=lattice->nx+1;j=lattice->ny+1;
	ii=i-cx[0];jj=j-cy[0];
	lattice->f[0][i][j]=lattice->f0[0][ii][jj]+calc_rhs_collide(0,i,j);
	ii=i-cx[1];jj=j-cy[1];
	lattice->f[1][i][j]=lattice->f0[1][ii][jj]+calc_rhs_collide(1,i,j);
	ii=i-cx[2];jj=j-cy[2];
	lattice->f[2][i][j]=lattice->f0[2][ii][jj]+calc_rhs_collide(2,i,j);
	ii=i-cx[5];jj=j-cy[5];
	lattice->f[5][i][j]=lattice->f0[5][ii][jj]+calc_rhs_collide(5,i,j);	
}

void LBM2dPlates::collide_corner1(){
	//down
	int i0=0,i1=lattice->nx+1,j=0;
	double c=delp-(lattice->f[0][i0][j]-lattice->f[0][i1][j]
	              +lattice->f[2][i0][j]-lattice->f[2][i1][j]
	              +lattice->f[4][i0][j]-lattice->f[4][i1][j])/3.0;/**/
	/*double c=1.2*delp
		      -0.4*(lattice->f[0][i0][j]-lattice->f[0][i1][j]
	             +lattice->f[2][i0][j]-lattice->f[2][i1][j]
	             +lattice->f[4][i0][j]-lattice->f[4][i1][j]
	             +lattice->f[5][i0][j]-lattice->f[5][i1][j]
	             +lattice->f[6][i0][j]-lattice->f[6][i1][j]);*/
	lattice->f[1][i0][j]=lattice->f[1][i1][j]+c;
	lattice->f[8][i0][j]=lattice->f[8][i1][j]+0.25*c;
	//
	lattice->f[5][i0][j]=lattice->f[7][i0][j];//
	lattice->f[2][i0][j]=lattice->f[4][i0][j];
	lattice->f[6][i0][j]=lattice->f[8][i0][j];
	///////////////////////////
	lattice->f[3][i1][j]=lattice->f[3][i0][j]-c;
	lattice->f[7][i1][j]=lattice->f[7][i0][j]-0.25*c;
	//
	lattice->f[6][i1][j]=lattice->f[8][i1][j];//
	lattice->f[2][i1][j]=lattice->f[4][i1][j];
	lattice->f[5][i1][j]=lattice->f[7][i1][j];
	//up
	j=lattice->ny+1;
	/*c=delp-(lattice->f[0][i0][j]-lattice->f[0][i1][j]
	       +lattice->f[2][i0][j]-lattice->f[2][i1][j]
	       +lattice->f[2][i0][j]-lattice->f[2][i1][j])/3.0;*/
	c=delp-(lattice->f[0][i0][j]-lattice->f[0][i1][j]
	       +lattice->f[2][i0][j]-lattice->f[2][i1][j]
	       +lattice->f[4][i0][j]-lattice->f[4][i1][j])/3.0;
	lattice->f[1][i0][j]=lattice->f[1][i1][j]+c;
	lattice->f[5][i0][j]=lattice->f[5][i1][j]+0.25*c;
	//
	lattice->f[8][i0][j]=lattice->f[6][i0][j];//
	lattice->f[4][i0][j]=lattice->f[2][i0][j];
	lattice->f[7][i0][j]=lattice->f[5][i0][j];
	////////////////////////
	lattice->f[3][i1][j]=lattice->f[3][i0][j]-c;
	lattice->f[6][i1][j]=lattice->f[6][i0][j]-0.25*c;
	//
	lattice->f[7][i1][j]=lattice->f[5][i1][j];//
	lattice->f[4][i1][j]=lattice->f[2][i1][j];
	lattice->f[8][i1][j]=lattice->f[6][i1][j];
}

void LBM2dPlates::collide_corner2(){
	int i0=0,j0=0,i1=lattice->nx+1,j1=0;
	//2,5,6
	double rhow=lattice->f[0][i0][j0]+lattice->f[1][i0][j0]+lattice->f[3][i0][j0]
		         +2.0*(lattice->f[4][i0][j0]+lattice->f[7][i0][j0]+lattice->f[8][i0][j0]);
	rhow/=1.0-v0;
	double rhop=rhow*v0+lattice->f[4][i0][j0]+lattice->f[7][i0][j0]+lattice->f[8][i0][j0];
	rhop*=6./(1.+3.*v0+3.*v0*v0);
	double upx=rhow*u0-(lattice->f[1][i0][j0]-lattice->f[3][i0][j0]
	                   +lattice->f[8][i0][j0]-lattice->f[7][i0][j0]);
	upx*=6./rhop;
	upx-=u0+3.*u0*v0;
	upx/=1.+3.*v0;
	lattice->f[2][i0][j0]=rhop/9.*(1.+3.*v0+4.5*v0*v0
	                            -1.5*((u0+upx)*(u0+upx)+v0*v0));
	lattice->f[5][i0][j0]=rhop/36.*(1.+3.*(u0+upx+v0)
	                               +4.5*(u0+upx+v0)*(u0+upx+v0)
	                               -1.5*((u0+upx)*(u0+upx)+v0*v0));
	lattice->f[6][i0][j0]=rhop/36.*(1.+3.*(-u0-upx+v0)
	                               +4.5*(-u0-upx+v0)*(-u0-upx+v0)
	                               -1.5*((u0+upx)*(u0+upx)+v0*v0));
	//
	//i=lattice->nx+1;j=0;
	rhow=lattice->f[0][i1][j1]+lattice->f[1][i1][j1]+lattice->f[3][i1][j1]
	    +2.0*(lattice->f[4][i1][j1]+lattice->f[7][i1][j1]+lattice->f[8][i1][j1]);
	rhow/=1.0-v0;
	rhop=rhow*v0+lattice->f[4][i1][j1]+lattice->f[7][i1][j1]+lattice->f[8][i1][j1];
	rhop*=6./(1.+3.*v0+3.*v0*v0);
	upx=rhow*u0-(lattice->f[1][i1][j1]-lattice->f[3][i1][j1]
	            +lattice->f[8][i1][j1]-lattice->f[7][i1][j1]);
	upx*=6./rhop;
	upx-=u0+3.*u0*v0;
	upx/=1.+3.*v0;
	lattice->f[2][i1][j1]=rhop/9.*(1.+3.*v0+4.5*v0*v0
	                            -1.5*((u0+upx)*(u0+upx)+v0*v0));
	lattice->f[5][i1][j1]=rhop/36.*(1.+3.*(u0+upx+v0)
	                               +4.5*(u0+upx+v0)*(u0+upx+v0)
	                               -1.5*((u0+upx)*(u0+upx)+v0*v0));
	lattice->f[6][i1][j1]=rhop/36.*(1.+3.*(-u0-upx+v0)
	                               +4.5*(-u0-upx+v0)*(-u0-upx+v0)
	                               -1.5*((u0+upx)*(u0+upx)+v0*v0));
	//4,7,8
	i0=0;j0=lattice->ny+1;i1=lattice->nx+1;j1=lattice->ny+1;
	rhow=lattice->f[0][i0][j0]+lattice->f[1][i0][j0]+lattice->f[3][i0][j0]
	    +2.0*(lattice->f[2][i0][j0]+lattice->f[5][i0][j0]+lattice->f[6][i0][j0]);
	rhow/=1.0-v2;
	rhop=rhow*v2+lattice->f[2][i0][j0]+lattice->f[5][i0][j0]+lattice->f[6][i0][j0];
	rhop*=6./(1.+3.*v2+3.*v2*v2);
	upx=rhow*u2-(lattice->f[1][i0][j0]-lattice->f[3][i0][j0]
	            +lattice->f[5][i0][j0]-lattice->f[6][i0][j0]);
	upx*=6./rhop;
	upx-=u2+3.*u2*v2;
	upx/=1.+3.*v2;
	lattice->f[4][i0][j0]=rhop/9.*(1.+3.*v2+4.5*v2*v2
	                            -1.5*((u2+upx)*(u2+upx)+v2*v2));
	lattice->f[8][i0][j0]=rhop/36.*(1.+3.*(u2+upx+v2)
	                               +4.5*(u2+upx+v2)*(u2+upx+v2)
	                               -1.5*((u2+upx)*(u2+upx)+v2*v2));
	lattice->f[7][i0][j0]=rhop/36.*(1.+3.*(-u2-upx+v2)
	                               +4.5*(-u2-upx+v2)*(-u2-upx+v2)
	                               -1.5*((u2+upx)*(u2+upx)+v2*v2));
	//
	rhow=lattice->f[0][i1][j1]+lattice->f[1][i1][j1]+lattice->f[3][i1][j1]
	    +2.0*(lattice->f[2][i1][j1]+lattice->f[5][i1][j1]+lattice->f[6][i1][j1]);
	rhow/=1.0-v2;
	rhop=rhow*v2+lattice->f[2][i1][j1]+lattice->f[5][i1][j1]+lattice->f[6][i1][j1];
	rhop*=6./(1.+3.*v2+3.*v2*v2);
	upx=rhow*u2-(lattice->f[1][i1][j1]-lattice->f[3][i1][j1]
	            +lattice->f[5][i1][j1]-lattice->f[6][i1][j1]);
	upx*=6./rhop;
	upx-=u2+3.*u2*v2;
	upx/=1.+3.*v2;
	lattice->f[4][i1][j1]=rhop/9.*(1.+3.*v2+4.5*v2*v2
	                            -1.5*((u2+upx)*(u2+upx)+v2*v2));
	lattice->f[8][i1][j1]=rhop/36.*(1.+3.*(u2+upx+v2)
	                               +4.5*(u2+upx+v2)*(u2+upx+v2)
	                               -1.5*((u2+upx)*(u2+upx)+v2*v2));
	lattice->f[7][i1][j1]=rhop/36.*(1.+3.*(-u2-upx+v2)
	                               +4.5*(-u2-upx+v2)*(-u2-upx+v2)
	                               -1.5*((u2+upx)*(u2+upx)+v2*v2));
}

void LBM2dPlates::add_boundary_force(){
}

void LBM2dPlates::set_boundary_viscosity(){
}
