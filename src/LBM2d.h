#ifndef _LBM2D_H_
#define _LBM2D_H_

#include "udfmanager.h"
#include "lbm/LBM2dSimulator.h"
#include "Timer.h"
#include <cmath>
#include <sstream>
using namespace std;

class LBM2d{
public:
	LBM2d(UDFManager* in,UDFManager* out);
	virtual ~LBM2d();
	virtual void update();
protected:
	virtual void input_lbm();
	virtual void input_other();
	virtual void output();
	virtual void nan_check(double x,char c='0');
private:
	//virtual void field_test();
	UDFManager *in,*out;
	Timer* timer;
	int iteration;
	//flow
	LBM2dSimulator *lbm;
	LBM2dLattice *lattice;
	//output
	bool fou,fop,forho,fof;
};

#endif // _MIST_SIMULATOR2D_H_
