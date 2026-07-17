//triton global OMP
//Xiaohua Li, Chuntian Liang and Chonghai Cai, Nucl. Phys. A 789 (2007) 103-113
//Global triton optical model potential
#include "/home/caijw/prjs/include/cpp/outColor.h"
void Li_2007(const double &A,const double &Z, const double &Kt){
	double rc = 1.4219;
	double V=137.6-0.1456*Kt+0.0436*pow(Kt,2)+4.3751*(A-Z-Z)/A+1.0474*Z/cbrt(A);
	double rr=1.1201-0.1504/cbrt(A);
	double ar=0.6833+0.0191*cbrt(A);
	double Wv=7.383+0.5025*Kt-pow(0.0097,2);
	double rv=1.3202-0.1776/cbrt(A);
	double av=1.119+0.01913*cbrt(A);
	double Ws=37.06-0.6451*Kt-47.19*(A-Z-Z)/A;
	double rs=1.251-0.4622/cbrt(A);
	double as=0.8114+0.01159*cbrt(A);
	const double Vso=1.9029;
	double rso=0.46991+0.1294/cbrt(A);
	double aso=0.3545-0.0522*cbrt(A);

	cout<<BOLDGREEN<<"Xiaohua Li, Chuntian Liang & Chonghai Cai, NPA 2007"<<RESET<<"\n";
	cout<<GREEN<<"\ttype=0 p(1:3)="<<A<<" 0 "<<rc<<"/\n";
	cout<<"\ttype=1 p(1:6)="<<V<<" "<<rr<<" "<<ar<<" "<<Wv<<" "<<rv<<" "<<av<<"/\n";
	cout<<"\ttype=2 p(4:6)="<<Ws<<" "<<rs<<" "<<as<<"/\n";
	cout<<"\ttype=3 p(1:3)="<<Vso<<" "<<rso<<" "<<aso<<"/"<<RESET<<endl;
}
//D. Y. Pang, P. Roussel-Chomaz and et al, Phys. Rev. C 79, 024615 (2009)
//entitled Global optical model potential for A = 3 projectiles
void Pang_2009(const double &A,const double &Z,const double &Kt){
	const double rc=1.24, rc0=0.12;
	double rc_act=(rc*cbrt(A)+rc0)/cbrt(A);
	double EC=6.*1.*Z*1.44/(5.*(rc*cbrt(A)+rc0));
	double V=118.3-0.13*(Kt-EC);
	double rr=(1.3*cbrt(A)-0.48)/cbrt(A);
	const double ar=0.82;
	double Wv=38.5/(1.+exp((156.1-(Kt-EC))/52.4));
	double Ws_t=(35.0-34.2*(A-Z-Z)/A)/(1.+exp(((Kt-EC)-30.8)/106.4));
	double Ws_He3=(35.0+34.2*(A-Z-Z)/A)/(1.+exp(((Kt-EC)-30.8)/106.4));
	double rw=(1.31*cbrt(A)-0.13)/cbrt(A);
	double rs=rw;
	const double aw=0.84;
	const double as=aw;

	cout<<BOLDMAGENTA<<"D. Y. Pang, P. Roussel-Chomaz & et al, PRC 2009"<<RESET<<"\n";
	cout<<MAGENTA<<"\ttype=0 p(1:3)="<<A<<" 0 "<<rc_act<<"/\n";
	cout<<"\ttype=1 p(1:6)="<<V<<" "<<rr<<" "<<ar<<" "<<Wv<<" "<<rw<<" "<<aw<<"/\n";
	cout<<"\ttype=2 p(4:6)="<<Ws_t<<" "<<rs<<" "<<as<<"/"<<RESET<<endl;
}
//D. Y. Pang, W. M. Dean, and A. M. Mukhamedzhanov, Phys. Rev. C 91, 024611 (2015)
//entitled Optical model potential of A=3 projectiles for 1p-shell nuclei
void Pang_2015(const double &A,const double &Z,const double &Kt){
	const double rc=1.30;
	double EC=6.*1.*Z*1.44/(5.*rc*cbrt(A));
	double Vr=155.1-0.678*(Kt-EC);
	double rr=(0.920*cbrt(A)+0.108+0.0031*(Kt-EC))/cbrt(A);
	const double ar=0.792;
	double Wv=33.1/(1.+exp((156.1-(Kt-EC))/52.4));
	double Ws_t=(21.8-13.1*(A-Z-Z)/A)/(1.+exp(((Kt-EC)-30.8)/106.4));
	double Ws_He3=(21.8+13.1*(A-Z-Z)/A)/(1.+exp(((Kt-EC)-30.8)/106.4));
	double rw=(1.43*cbrt(A)-0.16)/cbrt(A);
	double rs=rw;
	const double aw=0.801;
	const double as=aw;

	cout<<BOLDCYAN<<"D. Y. Pang, W. M. Dean & et al, PRC 2015"<<RESET<<"\n";
	cout<<CYAN<<"\ttype=0 p(1:3)="<<A<<" 0 "<<rc<<"/\n";
	cout<<"\ttype=1 p(1:6)="<<Vr<<" "<<rr<<" "<<ar<<" "<<Wv<<" "<<rw<<" "<<aw<<"/\n";
	cout<<"\ttype=2 p(4:6)="<<Ws_t<<" "<<rs<<" "<<as<<"/"<<RESET<<endl;
}
void t_gOMP(double A,double Z, double Kt){
	cout<<"incident "<<BOLDRED<<"triton"<<RESET<<": "<<Kt<<" MeV; Target A="<<A<<", Z="<<Z<<"\n";
	Li_2007(A,Z,Kt);
	Pang_2009(A,Z,Kt);
	Pang_2015(A,Z,Kt);
}
