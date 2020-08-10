#include "LBM2dLattice.h"

LBM2dLattice::LBM2dLattice(int _nx,int _ny,double _dx)
:nx(_nx),ny(_ny),lx(_nx*_dx),ly(_ny*_dx),dx(_dx),
	u(0,nx,0,ny),v(0,nx,0,ny),
	//p(0,nx+1,0,ny+1),
	bfx(0,nx,0,ny),bfy(0,nx,0,ny),
	rho(0,nx,0,ny),vis(0,nx,0,ny)
{
	for(int i=0;i<9;i++){
		f[i].setBounds(0,nx,0,ny);
		f0[i].setBounds(0,nx,0,ny);
	}
}

LBM2dLattice::~LBM2dLattice(){
}
