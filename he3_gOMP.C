//He3 global OMP
//D.Y.Pang et al, Phys. Rev. C 79, 024615 (2009)
//entitled Global optical model potential for A = 3 projectiles
#include "/home/caijw/prjs/include/cpp/outColor.h"
void Pang_2009(const double &A,const double &Z,const double &K){
	const double rc=1.24, rc0=0.12;
	double rc_act=(rc*cbrt(A)+rc0)/cbrt(A);
	double EC=6.*1.*Z*1.44/(5.*(rc*cbrt(A)+rc0));
	double V=118.3-0.13*(K-EC);
	double rr=(1.3*cbrt(A)-0.48)/cbrt(A);
	const double ar=0.82;
	double Wv=38.5/(1.+exp((156.1-(K-EC))/52.4));
	double Ws_t=(35.0-34.2*(A-Z-Z)/A)/(1.+exp(((K-EC)-30.8)/106.4));
	double Ws_He3=(35.0+34.2*(A-Z-Z)/A)/(1.+exp(((K-EC)-30.8)/106.4));
	double rw=(1.31*cbrt(A)-0.13)/cbrt(A);
	double rs=rw;
	const double aw=0.84;
	const double as=aw;

	cout<<BOLDGREEN<<"D. Y. Pang, P. Roussel-Chomaz & et al"<<RESET<<"\n";
	cout<<GREEN<<"type=0 p(1:3)="<<A<<" 0 "<<rc_act<<"/\n";
	cout<<"type=1 p(1:6)="<<V<<" "<<rr<<" "<<ar<<" "<<Wv<<" "<<rw<<" "<<aw<<"/\n";
	cout<<"type=2 p(4:6)="<<Ws_He3<<" "<<rs<<" "<<as<<"/"<<RESET<<endl;
}
//D. Y. Pang, W. M. Dean, and A. M. Mukhamedzhanov, Phys. Rev. C 91, 024611 (2015)
//entitled Optical model potential of A=3 projectiles for 1p-shell nuclei
void Pang_2015(const double &A,const double &Z,const double &K){
	const double rc=1.30;
	double EC=6.*1.*Z*1.44/(5.*rc*cbrt(A));
	double Vr=155.1-0.678*(K-EC);
	double rr=(0.920*cbrt(A)+0.108+0.0031*(K-EC))/cbrt(A);
	const double ar=0.792;
	double Wv=33.1/(1.+exp((156.1-(K-EC))/52.4));
	double Ws_t=(21.8-13.1*(A-Z-Z)/A)/(1.+exp(((K-EC)-30.8)/106.4));
	double Ws_He3=(21.8+13.1*(A-Z-Z)/A)/(1.+exp(((K-EC)-30.8)/106.4));
	double rw=(1.43*cbrt(A)-0.16)/cbrt(A);
	double rs=rw;
	const double aw=0.801;
	const double as=aw;

	cout<<BOLDCYAN<<"D. Y. Pang, W. M. Dean & et al"<<RESET<<"\n";
	cout<<CYAN<<"\ttype=0 p(1:3)="<<A<<" 0 "<<rc<<"/\n";
	cout<<"\ttype=1 p(1:6)="<<Vr<<" "<<rr<<" "<<ar<<" "<<Wv<<" "<<rw<<" "<<aw<<"/\n";
	cout<<"\ttype=2 p(4:6)="<<Ws_He3<<" "<<rs<<" "<<as<<"/"<<RESET<<endl;
}
void he3_gOMP(double A,double Z, double KHe3){
	cout<<"incident 3He: "<<KHe3<<" MeV; Target A="<<A<<", Z="<<Z<<"/"<<endl;
	Pang_2009(A,Z,KHe3);
	Pang_2015(A,Z,KHe3);
}
