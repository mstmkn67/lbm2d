#include "udf/gourmain.h"
#include "LBM2d.h"

void udfHeaderCheck()
{
	string version("1.0"),engine("lbm2d");
	cout << "**************************************************************" << endl;
	cout <<  "              " <<  engine << "  " << version << "           " << endl;
	cout << "                                        Masato MAKINO         " << endl;
	cout << "**************************************************************" << endl;
	cout <<  endl;
}

void error_massage(){
	cout << "usage: lbm2d -I input_udf [-O output_udf] " << endl;
}


int gourmain(UDFManager* in,UDFManager* out){
	udfHeaderCheck();
	LBM2d* sim=new LBM2d(in,out);
	sim->update();
	delete sim;
	return 0;
}
