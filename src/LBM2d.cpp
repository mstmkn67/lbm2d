#include "LBM2d.h"
//
#include "lbm/LBM2dPlates.h"
#include "lbm/LBM2dLeesEdwards.h"

LBM2d::LBM2d(UDFManager* _in,UDFManager* _out)
:in(_in),out(_out),iteration(0),
 fou(false),fop(false),forho(false),fof(false){
	timer=new Timer;
}

LBM2d::~LBM2d(){
	delete lbm;
	delete lattice;
	delete timer;
}

void LBM2d::update(){
	input_lbm();
	input_other();
	lbm->initial();
	lbm->set_bulk_viscosity();
	output();
	//
	int total=in->i("simulation.time.simulation_steps");
	int report=in->i("simulation.time.record_steps");
	for(iteration=1;iteration<=total;iteration++){
		lbm->update(); //
		if(iteration%report==0){
			cout << "report step at " << iteration << endl;
			output();
		}
	}
}

void LBM2d::input_lbm(){
	double dt=in->d("simulation.time.dt");
	double viscosity=in->d("fluid.viscosity");
	double density=in->d("fluid.mass_density");
	//
	int nx=in->i("system_size.nx");
	int ny=in->i("system_size.ny");
	double dx=in->d("system_size.dx");
	lattice=new LBM2dLattice(nx,ny,dx);
	if(in->s("simulation.boundary_condition.type")=="periodic"){
		lbm=new LBM2dLeesEdwards(lattice,viscosity,density,dt,
			                       &iteration,0.0);
	}else if(in->s("simulation.boundary_condition.type")=="Couette"){
		double vxu=in->d("simulation.boundary_condition.Couette.velocity_x_up");
		double vxd=in->d("simulation.boundary_condition.Couette.velocity_x_down");
		lbm=new LBM2dPlates(lattice,viscosity,density,dt,
		                    0.0, vxd,0.0, vxu,0.0);
	}else if(in->s("simulation.boundary_condition.type")=="Poiseuille"){
		double delp=in->d("simulation.boundary_condition.Poiseuille.delta_pressure");
		lbm=new LBM2dPlates(lattice,viscosity,density,dt,
		                    delp, 0.0,0.0, 0.0,0.0);
	}else if(in->s("simulation.boundary_condition.type")=="LeesEdwards"){
		double gd=in->d("simulation.boundary_condition.LeesEdwards.shear_rate");
		lbm=new LBM2dLeesEdwards(lattice,viscosity,density,dt,
			                       &iteration, gd);
	}
}

void LBM2d::input_other(){
	if(in->s("simulation.field_output.u")=="true")fou=true;
	if(in->s("simulation.field_output.p")=="true")fop=true;
	if(in->s("simulation.field_output.rho")=="true")forho=true;
	if(in->s("simulation.field_output.f")=="true")fof=true;
}

void LBM2d::output(){
	out->newRecord();
	for(int i=0;i<=lattice->nx;i++){
		for(int j=0;j<=lattice->ny;j++){
			if(fou){
				double u=lattice->u[i][j];
				nan_check(u,'u');//////////////////////////////////
				out->put("simulation_result.lattice.u[][]",u,INDEX(i,j));
				double v=lattice->v[i][j];
				nan_check(v,'v');/////////////////////////
				out->put("simulation_result.lattice.v[][]",v,INDEX(i,j));
			}
			if(fop){
				double gp=-lbm->delp/lattice->lx;
				double x=lattice->dx*i;
				nan_check(lattice->rho[i][j],'p');/////////////////////////
				out->put("simulation_result.lattice.p[][]",lattice->rho[i][j]/3.0+gp*x,INDEX(i,j));
			}
			if(forho){
				out->put("simulation_result.lattice.rho[][]",lattice->rho[i][j],INDEX(i,j));
			}
			if(fof){
				for(int n=0;n<9;n++){
					double f=lattice->f[n][i][j];
					nan_check(f,'f');//////////////////////////////////
					out->put("simulation_result.lattice.f[][][]",f,INDEX(n,i,j));
				}
			}
		}
	}
	//
	out->put("simulation_result.cpu_time",timer->get());
	out->write();
}

void LBM2d::nan_check(double x,char c){
  if( isnan(x) == true ){ 
    cout  << "Not a Number is found: " << c << endl;
    exit(1);
  }
}
