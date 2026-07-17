//proton global OMP
//from A.J. Koning and J.P. Delaroche, Nucl. Phys. A 713 (2003) 231-310
//entitled Local and global nucleon optical models from 1 keV to 200 MeV
#include "/home/caijw/prjs/include/cpp/outColor.h"
void koning_2003(const double &A,const double &Z,const double &Kp){
	double vps[4]={59.30+21.0*(A-Z-Z)/A-0.024*A, 0.007067+4.23e-6*A, 1.729e-5+1.136e-8*A, 7e-9};
	double wps[2]={14.667+0.009629*A, 73.55+0.0795*A};
	double dps[3]={16.-16.*(A-Z-Z)/A, 0.018+0.003802/(1+exp((A-156.)/8.)), 11.5};
	double vsops[2]={5.922+0.003*A, 0.004};
	double wsops[2]={-3.1, 160.};
	double Epf=-8.4075+0.01378*A;
	double rC=1.198+0.697/pow(A,2./3.)+12.994/pow(A,5./3.);
	double VC_bar=1.73/(rC*Z/cbrt(A));
	double delta_E=Kp-Epf;

	double Vv=vps[0]*(1- vps[1]*delta_E+ vps[2]*pow(delta_E,2)- vps[3]*pow(delta_E,3))+
			VC_bar*vps[0]*(vps[1]-2*vps[2]*delta_E+ 3*vps[3]*pow(delta_E,2));
	double Wv=wps[0]*pow(delta_E,2)/(pow(delta_E,2) +wps[1]*wps[1]);
	double rV=1.3039-0.4054/cbrt(A);
	double aV=0.6778-1.487e-4*A;
	double WD=dps[0]*pow(delta_E,2)*exp(-dps[1]*delta_E)/(pow(delta_E,2) +dps[2]*dps[2]);
	double rD=1.3424-0.01585*cbrt(A);
	double aD=0.5187+5.205e-4*A;
	double VSO=vsops[0]*exp(-vsops[1]*delta_E);
	double WSO=wsops[0]*pow(delta_E,2)/(pow(delta_E,2) +wsops[1]*wsops[1]);
	double rSO=1.1854-0.647/cbrt(A);
	const double aSO=0.59;

	cout<<BOLDGREEN<<"A.J. Koning and J.P. Delaroche"<<RESET<<"\n";
	cout<<GREEN<<"type=0 p(1:3)="<<A<<" 0 "<<rC<<"/"<<endl;
	cout<<"type=1 p(1:6)="<<Vv<<" "<<rV<<" "<<aV<<" "<<Wv<<" "<<rV<<" "<<aV<<"/"<<endl;
	cout<<"type=2 p(4:6)="<<WD<<" "<<rD<<" "<<aD<<"/"<<endl;
	cout<<"type=3 p(1:6)="<<max(VSO,0.0)<<" "<<rSO<<" "<<aSO<<" "<<max(WSO,0.0)<<" "<<rSO<<" "<<aSO<<"/"<<RESET<<endl;
}
//from R.L. Varner, W.J. Thompson and et al, Phys. Rep. 201, 57 (1991)
//entitled A global nucleon optical model potential
void varner_1991(const double &A,const double &Z,const double &Kp){
	const double V0=52.9, Vt=13.1,Ve=-0.299;
	const double r0=1.25, r00=-0.225,a0=0.690;
	double Rc=1.238*cbrt(A)+0.116;
	double EC=1.73*Z/Rc;
	double rc=Rc/cbrt(A);
	const double rso=1.34,rso0=-1.2;
	const double Wv0=7.8, Wve0=35,Wvew=16;
	const double Ws0=10., Wst=18,Wse0=36,Wsew=37;
	const double rw=1.33, rw0=-0.42,aw=0.69;

	double Vr=V0+Vt*(A-Z-Z)/A+(Kp-EC)*Ve;
	double rVr=(r0*cbrt(A)+r00)/cbrt(A);
	double Wv=Wv0/(1.+exp((Wve0-(Kp-EC))/Wvew));
	double Ws=(Ws0+Wst*(A-Z-Z)/A)/(1.+exp(((Kp-EC)-Wse0)/Wsew));
	double rVso=(rso*cbrt(A)+rso0)/cbrt(A);
	double rW=(rw*cbrt(A)+rw0)/cbrt(A);
	const double Vso=5.9, aso=0.63;

	cout<<BOLDYELLOW<<"R.L. Varner, W.J. Thompson and et al"<<RESET<<"\n";
	cout<<YELLOW<<"type=0 p(1:3)="<<A<<" 0 "<<rc<<"/"<<endl;
	cout<<"type=1 p(1:6)="<<Vr<<" "<<rVr<<" "<<a0<<" "<<Wv<<" "<<rW<<" "<<aw<<"/"<<endl;
	cout<<"type=2 p(4:6)="<<Ws<<" "<<rW<<" "<<aw<<"/"<<endl;
	cout<<"type=3 p(1:3)="<<Vso<<" "<<rVso<<" "<<aso<<"/"<<RESET<<endl;
}
void p_gOMP(double A,double Z, double Kp){
	if(A<Z){
		cout<<BOLDRED<<"A must NOT be smaller than Z!!!"<<RESET<<"\n";
		exit(1);
	}
	cout<<"incident "<<BOLDRED<<"proton"<<RESET<<": "<<Kp<<" MeV; Target A="<<A<<", Z="<<Z<<"/"<<endl;
	koning_2003(A,Z,Kp);
	varner_1991(A,Z,Kp);
}
