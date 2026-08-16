//neutron global OMP
//from A.J. Koning and J.P. Delaroche, Nucl. Phys. A 713 (2003) 231-310
//entitled Local and global nucleon optical models from 1 keV to 200 MeV
#include "/home/caijw/prjs/include/cpp/outColor.h"
void koning_2003_n(const double &A,const double &Z,const double &Kn){
	double vns[4]={59.30-21.0*(A-Z-Z)/A-0.024*A, 0.007228-1.48e-6*A, 1.994e-5-2e-8*A, 7e-9};
	double wns[2]={12.195+0.0167*A, 73.55+0.0795*A};
	double dns[3]={16.-16.*(A-Z-Z)/A, 0.018+0.003802/(1+exp((A-156.)/8.)), 11.5};
	double vsons[2]={5.922+0.003*A, 0.004};
	double wsons[2]={-3.1, 160.};
	double Enf=-11.2814+0.02646*A;
	double delta_E=Kn-Enf;

	double Vv=vns[0]*(1- vns[1]*delta_E+ vns[2]*pow(delta_E,2)- vns[3]*pow(delta_E,3));
	double Wv=wns[0]*pow(delta_E,2)/(pow(delta_E,2) +wns[1]*wns[1]);
	double rV=1.3039-0.4054/cbrt(A);
	double aV=0.6778-1.487e-4*A;
	double WD=dns[0]*pow(delta_E,2)*exp(-dns[1]*delta_E)/(pow(delta_E,2) +dns[2]*dns[2]);
	double rD=1.3424-0.01585*cbrt(A);
	double aD=0.5446-1.656e-4*A;
	double VSO=vsons[0]*exp(-vsons[1]*delta_E);
	double WSO=wsons[0]*pow(delta_E,2)/(pow(delta_E,2) +wsons[1]*wsons[1]);
	double rSO=1.1854-0.647/cbrt(A);
	const double aSO=0.59;

	cout<<BOLDGREEN<<"A.J. Koning and J.P. Delaroche"<<RESET<<"\n";
	cout<<GREEN<<"type=1 p(1:7)="<<Vv<<" "<<rV<<" "<<aV<<" "<<Wv<<" "<<rV<<" "<<aV<<" "<<A<<"/"<<endl;
	cout<<"type=2 p(4:6)="<<WD<<" "<<rD<<" "<<aD<<"/"<<endl;
	cout<<"type=3 p(1:6)="<<max(VSO,0.0)<<" "<<rSO<<" "<<aSO<<" "<<max(WSO,0.0)<<" "<<rSO<<" "<<aSO<<"/"<<RESET<<endl;
}
//from R.L. Varner, W.J. Thompson and et al, Phys. Rep. 201, 57 (1991)
//entitled A global nucleon optical model potential
void varner_1991_n(const double &A,const double &Z,const double &Kn){
	const double V0=52.9, Vt=13.1,Ve=-0.299;
	const double r0=1.25, r00=-0.225,a0=0.690;
	double Rc=1.238*cbrt(A)+0.116;
	double EC=1.73*Z/Rc;
	const double rso=1.34,rso0=-1.2;
	const double Wv0=7.8, Wve0=35,Wvew=16;
	const double Ws0=10., Wst=18,Wse0=36,Wsew=37;
	const double rw=1.33, rw0=-0.42,aw=0.69;

	double Vr=V0-Vt*(A-Z-Z)/A+(Kn-EC)*Ve;
	double rVr=(r0*cbrt(A)+r00)/cbrt(A);
	double Wv=Wv0/(1.+exp((Wve0-(Kn-EC))/Wvew));
	double Ws=(Ws0+Wst*(A-Z-Z)/A)/(1.+exp(((Kn-EC)-Wse0)/Wsew));
	double rVso=(rso*cbrt(A)+rso0)/cbrt(A);
	double rW=(rw*cbrt(A)+rw0)/cbrt(A);
	const double Vso=5.9, aso=0.63;

	cout<<BOLDYELLOW<<"R.L. Varner, W.J. Thompson and et al"<<RESET<<"\n";
	cout<<YELLOW<<"type=1 p(1:7)="<<Vr<<" "<<rVr<<" "<<a0<<" "<<Wv<<" "<<rW<<" "<<aw<<" "<<A<<endl;
	cout<<"type=2 p(4:6)="<<Ws<<" "<<rW<<" "<<aw<<"/"<<endl;
	cout<<"type=3 p(1:3)="<<Vso<<" "<<rVso<<" "<<aso<<"/"<<RESET<<endl;
}
void n_gOMP(double A,double Z, double Kn){
	if(A<Z){
		cout<<BOLDRED<<"A must NOT be smaller than Z!!!"<<RESET<<"\n";
		exit(1);
	}
	cout<<"incident "<<BOLDRED<<"neutron"<<RESET<<": "<<Kn<<" MeV; Target A="<<BOLDRED<<A<<RESET<<", Z="<<BOLDRED<<Z<<RESET<<"/"<<endl;
	koning_2003_n(A,Z,Kn);
	varner_1991_n(A,Z,Kn);
}
