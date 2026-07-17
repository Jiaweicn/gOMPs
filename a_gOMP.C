//alpha global OMP
//M. Nolte, H. Machner, and J. Bojowald, Phys. Rev. C 36, 1312 (1987)
//Global optical potential for alpha particles with energies above 80 MeV
#include "/home/caijw/prjs/include/cpp/outColor.h"
void a_gOMP(double A, double Z, double Ka) {
	const double a0= 101.1, a1= 6.051, a2=-0.248;
	const double b0= 26.82, b1=-1.706, b2=0.006;
	const double c0= 0.817, c1=-0.0085;
	const double d0= 0.692, d1=-0.020;
	double a3rd = cbrt(A);
   
	double v = a0 + a1 * Z / a3rd + a2 *Ka;
	double av = c0 + c1 * a3rd;
	double rv = 1.245;
	double w = b0 + b1 * a3rd + b2 *Ka;
	double aw = d0 + d1 * a3rd;
	double rw = 1.570;

	cout<<"incident "<<BOLDRED<<"alpha"<<RESET<<": "<<Ka<<" MeV; Target A="<<A<<", Z="<<Z<<"/"<<endl;
	printf("type=1 p1=%-7.3f p2=%-6.4f p3=%-6.4f p4=%-7.3f p5=%-6.4f p6=%-6.4f /\n",
          v, rv, av, w, rw, aw);
}
