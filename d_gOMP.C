//deuteron global OMP
//from Haixia An and Chonghai Cai, Phys. Rev. C 73, 054605 (2006)
//entitled Global deuteron optical model potential for the energy range up to 183 MeV
#include "/home/caijw/prjs/include/cpp/outColor.h"
void An_Cai_2006(const double &A,const double &Z,const double &Kd){
	double V=91.85-0.249*Kd+0.000116*Kd*Kd+0.642*Z/cbrt(A);
	double rr=1.152-0.00776/cbrt(A);
	double ar=0.719+0.0126*cbrt(A);
	double Wv=1.104+0.0622*Kd;
	double rv=1.305+0.0997/cbrt(A);
	double av=0.855-0.10*cbrt(A);
	double Ws=10.83-0.0306*Kd;
	double rs=1.334+0.152/cbrt(A);
	double as=0.531+0.062*cbrt(A);
	const double Vso=3.557;
	const double aso=1.011,rso=0.972,rc=1.303;

	cout<<BOLDGREEN<<"Haixia An & Chonghai Cai, 2006"<<RESET<<"\n";
	cout<<GREEN<<"\ttype=0 p(1:3)="<<A<<" 0 "<<rc<<"/\n";
	cout<<"\ttype=1 p(1:6)="<<V<<" "<<rr<<" "<<ar<<" "<<Wv<<" "<<rv<<" "<<av<<"/\n";
	cout<<"\ttype=2 p(4:6)="<<Ws<<" "<<rs<<" "<<as<<"/\n";
	cout<<"\ttype=3 p(1:3)="<<Vso<<" "<<rso<<" "<<aso<<"/"<<RESET<<endl;
}
//Y. Zhang, D. Y. Pang and J. L. Lou, Phys. Rev. C 94, 014619 (2016)
//entitled Optical model potential for deuteron elastic scattering with 1p-shell nuclei
void Zhang_2016(const double &A,const double &Z,const double &Kd){
	const double rc=1.3;
	double EC=6.*1.*Z*1.44/(5.*rc*cbrt(A));
	double rr=(1.11*cbrt(A)-0.167+0.00117*(Kd-EC))/cbrt(A);
	double rw=(0.561*cbrt(A)+3.07-0.00449*(Kd-EC))/cbrt(A);
	double V=98.9-0.279*(Kd-EC);
	double Wv=11.5/(1.+exp((18.1-(Kd-EC))/5.97));
	double Ws=7.56/(1.+exp(((Kd-EC)-14.3)/4.55));
	const double ar=0.775;
	const double aw=0.744;

	cout<<BOLDMAGENTA<<"Y. Zhang, D. Y. Pang & J. L. Lou, 2016"<<RESET<<"\n";
	cout<<MAGENTA<<"\ttype=0 p(1:3)="<<A<<" 0 "<<rc<<"/\n";
	cout<<"\ttype=1 p(1:6)="<<V<<" "<<rr<<" "<<ar<<" "<<Wv<<" "<<rw<<" "<<aw<<"/\n";
	cout<<"\ttype=2 p(4:6)="<<Ws<<" "<<rw<<" "<<aw<<"/"<<RESET<<endl;
}
//W. W. Daehnick, J. D. Childs and Z. Vrcelj, Phys. Rev. C 21, 2253 (1980)
//Global optical model potential for elastic deuteron scattering from 12 to 90 MeV
void Daehnick_1980_L(const double &A,const double &Z,const double &Kd){//L for nonrelativistic kinematics
	double Vr=88.5-0.26*Kd+0.88*Z/cbrt(A);
	const double r0=1.17;
	double a0=0.709+0.0017*Kd;
	double beta=-pow(Kd/100.,2);
	double Ms[]={8,20,28,50,82,126};//magic numbers
	short nMs=sizeof(Ms)/sizeof(Ms[0]);
	double nn=A-Z;//number of neutrons
	double shell_effect=0.;
	for(short i=0;i<nMs;i++){
		double mu=pow((Ms[i]-nn)*0.5,2);
		shell_effect=shell_effect+exp(-mu);
	}
	double Ws=(12.2+0.026*Kd)*(1-exp(beta));
	double Wd=(12.2+0.026*Kd)*exp(beta);
	double rI=1.325;
	double aI=0.53+0.07*cbrt(A)-0.04*shell_effect;
	const double rC=1.30;
	double Vls=7.33-0.029*Kd;
	const double rls=1.07,als=0.66;

	cout<<BOLDCYAN<<"W. W. Daehnick, J. D. Childs & Z. Vrcelj, 1980, for nonrelativistic kinematics"<<RESET<<"\n";
	cout<<CYAN<<"\ttype=0 p(1:3)="<<A<<" 0 "<<rC<<"/\n";
	cout<<"\ttype=1 p(1:6)="<<Vr<<" "<<r0<<" "<<a0<<" "<<Ws<<" "<<rI<<" "<<aI<<"/\n";
	cout<<"\ttype=2 p(4:6)="<<Wd<<" "<<rI<<" "<<aI<<"/\n";
	cout<<"\ttype=3 p(1:3)="<<Vls<<" "<<rls<<" "<<als<<"/"<<RESET<<endl;
}
void Daehnick_1980_F(const double &A,const double &Z,const double &Kd){//F for lativistic kinematics
	double Vr=88.0-0.283*Kd+0.88*Z/cbrt(A);
	const double r0=1.17;
	double a0=0.717+0.0012*Kd;
	double beta=-pow(Kd/100.,2);
	double Ms[]={8,20,28,50,82,126};//magic numbers
	short nMs=sizeof(Ms)/sizeof(Ms[0]);
	double nn=A-Z;//number of neutrons
	double shell_effect=0.;
	for(short i=0;i<nMs;i++){
		double mu=pow((Ms[i]-nn)*0.5,2);
		shell_effect=shell_effect+exp(-mu);
	}
	double Ws=(12.0+0.031*Kd)*(1-exp(beta));
	double Wd=(12.0+0.031*Kd)*exp(beta);
	double rI=1.376-0.01*sqrt(Kd);
	double aI=0.52+0.07*cbrt(A)-0.04*shell_effect;
	const double rC=1.30;
	double Vls=7.20-0.032*Kd;
	const double rls=1.07,als=0.66;

	cout<<BOLDCYAN<<"W. W. Daehnick, J. D. Childs & Z. Vrcelj, 1980, for relativistic kinematics"<<RESET<<"\n";
	cout<<CYAN<<"\ttype=0 p(1:3)="<<A<<" 0 "<<rC<<"/\n";
	cout<<"\ttype=1 p(1:6)="<<Vr<<" "<<r0<<" "<<a0<<" "<<Ws<<" "<<rI<<" "<<aI<<"/\n";
	cout<<"\ttype=2 p(4:6)="<<Wd<<" "<<rI<<" "<<aI<<"/\n";
	cout<<"\ttype=3 p(1:3)="<<Vls<<" "<<rls<<" "<<als<<"/"<<RESET<<endl;
}

void d_gOMP(double A,double Z,double Kd){
	cout<<"incident "<<BOLDRED<<"deuteron"<<RESET<<": "<<Kd<<" MeV; Target A="<<A<<", Z="<<Z<<"\n";
	An_Cai_2006(A,Z,Kd);
	Zhang_2016(A,Z,Kd);
	Daehnick_1980_L(A,Z,Kd);
	Daehnick_1980_F(A,Z,Kd);
}
