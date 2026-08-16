//proton global OMP
//from A.J. Koning and J.P. Delaroche, Nucl. Phys. A 713 (2003) 231-310
//entitled Local and global nucleon optical models from 1 keV to 200 MeV
#include "/home/caijw/prjs/include/cpp/outColor.h"
#include "p_gOMP.C"
#include "n_gOMP.C"

void koning_2003(double A,double Z, double Kd){
	koning_2003_p(A,Z,Kd/2);
	koning_2003_n(A,Z,Kd/2);
}
void varner_1991(double A,double Z, double Kd){
	varner_1991_p(A,Z,Kd/2);
	varner_1991_n(A,Z,Kd/2);
}
void d_ADWA_gOMP(double A,double Z, double Kd){
	if(A<Z){
		cout<<BOLDRED<<"A must NOT be smaller than Z!!!"<<RESET<<"\n";
		exit(1);
	}
	cout<<"incident "<<BOLDRED<<"deuteron"<<RESET<<": "<<Kd<<" MeV; Target A="<<A<<", Z="<<Z<<"/"<<endl;
	koning_2003(A,Z,Kd);
	varner_1991(A,Z,Kd);
}
