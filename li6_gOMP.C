//Li6 global OMP
//Yongli Xu1, Yinlu Han and et al, Phys. Rev. C 98, 024619 (2018)
//6Li global phenomenological optical model potential
#include "/home/caijw/prjs/include/cpp/outColor.h"
void li6_gOMP(double A, double Z, double Kli6) {
	double v0= 265.736,v1 = -0.183;
	double w0= 28.850,w1 = -0.0989;
	double u0= -5.226,u1 = 0.118,u2 = 0.000379;
	double rr= 1.120,ar = 0.814;
	double rs= 1.311,as = 0.939;
	double rv= 1.537,av = 0.726;
	double rc= 1.674;

	double vr = v0 + v1 * Kli6;
	double ws = TMath::Max(0.,w0 + w1 * Kli6);
	double wv = TMath::Max(0.,u0 + u1 * Kli6 + u2 * Kli6 * Kli6);

	cout<<"incident "<<BOLDRED<<"6Li"<<RESET<<": "<<Kli6<<" MeV; Target A="<<A<<", Z="<<Z<<"/"<<endl;
	printf("type=1 p1=%-7.3f p2=%-6.4f p3=%-6.4f p4=%-8.3f p5=%-6.4f p6=%-6.4f /\n",
          vr, rr, ar, wv, rv, av);
	printf("type=2 p4=%-7.3f p5=%-6.4f p6=%-6.4f /\n",ws, rs, as);
}
