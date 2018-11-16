#include <math.h>
#include "mp.h"
#include <stdlib.h>
#define EPS 1.0e-16
#define FPMIN 1.0e-30
#define MAXIT 10000
#define XMIN 2.0
#define NUSE1 5
#define NUSE2 5

static char vcid[] = "$Id: nrerror.c,v 3.2 1999/08/23 23:59:06 vicadmin Exp $";

/*Returns the modified Bessel functions ri = Iν , rk = Kν and their derivatives rip = Iν ,
rkp = Kν , for positive x and for xnu = ν ≥ 0. The relative accuracy is within one or two
significant digits of MYEPS. FPMIN is a number close to the machine’s smallest floating-point
number. All internal arithmetic is in double precision. To convert the entire routine to double
precision, change the float declarations above to double and decrease MYEPS to 10−16 . Also
convert the function beschb.*/

void bessik(real x, real xnu, real *ri, real *rk, real *rip, real *rkp)
{
	void beschb(real x, double *gam1, double *gam2, double *gampl,
		double *gammi);
	void nrerror(char error_text[]);
	int i,l,nl;
	double a,a1,b,c,d,del,del1,delh,dels,e,f,fact,fact2,ff,gam1,gam2,
		gammi,gampl,h,p,pimu,q,q1,q2,qnew,ril,ril1,rimu,rip1,ripl,
		ritemp,rk1,rkmu,rkmup,rktemp,s,sum,sum1,x2,xi,xi2,xmu,xmu2;
    
    if (x <= 0.0 || xnu < 0.0) nrerror("bad arguments in bessik");
	nl=(int)(xnu+0.5);
	xmu=xnu-nl;
	xmu2=xmu*xmu;
	xi=1.0/x;
	xi2=2.0*xi;
	h=xnu*xi;
	if (h < FPMIN) h=FPMIN;
	b=xi2*xnu;
	d=0.0;
	c=h;
	for (i=1;i<=MAXIT;i++) {
		b += xi2;
		d=1.0/(b+d);
		c=b+1.0/c;
		del=c*d;
		h=del*h;
		if (fabs(del-1.0) < EPS) break;
	}
	if (i > MAXIT) nrerror("x too large in bessik; try asymptotic expansion");
	ril=FPMIN;
	ripl=h*ril;
	ril1=ril;
	rip1=ripl;
	fact=xnu*xi;
	for (l=nl;l>=1;l--) {
		ritemp=fact*ril+ripl;
		fact -= xi;
		ripl=fact*ritemp+ril;
		ril=ritemp;
	}
	f=ripl/ril;
	if (x < XMIN) {
		x2=0.5*x;
		pimu=PI*xmu;
		fact = (fabs(pimu) < EPS ? 1.0 : pimu/sin(pimu));
		d = -log(x2);
		e=xmu*d;
		fact2 = (fabs(e) < EPS ? 1.0 : sinh(e)/e);
		beschb(xmu,&gam1,&gam2,&gampl,&gammi);
		ff=fact*(gam1*cosh(e)+gam2*fact2*d);
		sum=ff;
		e=exp(e);
		p=0.5*e/gampl;
		q=0.5/(e*gammi);
		c=1.0;
		d=x2*x2;
		sum1=p;
		for (i=1;i<=MAXIT;i++) {
			ff=(i*ff+p+q)/(i*i-xmu2);
			c *= (d/i);
			p /= (i-xmu);
			q /= (i+xmu);
			del=c*ff;
			sum += del;
			del1=c*(p-i*ff);
			sum1 += del1;
			if (fabs(del) < fabs(sum)*EPS) break;
		}
		if (i > MAXIT) nrerror("bessk series failed to converge");
		rkmu=sum;
		rk1=sum1*xi2;
	} else {
		b=2.0*(1.0+x);
		d=1.0/b;
		h=delh=d;
		q1=0.0;
		q2=1.0;
		a1=0.25-xmu2;
		q=c=a1;
		a = -a1;
		s=1.0+q*delh;
		for (i=2;i<=MAXIT;i++) {
			a -= 2*(i-1);
			c = -a*c/i;
			qnew=(q1-b*q2)/a;
			q1=q2;
			q2=qnew;
			q += c*qnew;
			b += 2.0;
			d=1.0/(b+a*d);
			delh=(b*d-1.0)*delh;
			h += delh;
			dels=q*delh;
			s += dels;
			if (fabs(dels/s) < EPS) break;
		}
		if (i > MAXIT) nrerror("bessik: failure to converge in cf2");
		h=a1*h;
		rkmu=sqrt(PI/(2.0*x))*exp(-x)/s;
		rk1=rkmu*(xmu+x+0.5-h)*xi;
	}
	rkmup=xmu*xi*rkmu-rk1;
	rimu=xi/(f*rkmu-rkmup);
	*ri=(rimu*ril1)/ril;
	*rip=(rimu*rip1)/ril;
	for (i=1;i<=nl;i++) {
		rktemp=(xmu+i)*xi2*rk1+rkmu;
		rkmu=rk1;
		rk1=rktemp;
	}
	*rk=rkmu;
	*rkp=xnu*xi*rkmu-rk1;
}


void beschb(real x, double *gam1, double *gam2, double *gampl, double *gammi)
{
	float chebev(real a, real b, real c[], int m, real x);
	real xx;
	static real c1[] = {
		-1.142022680371172e0,6.516511267076e-3,
		3.08709017308e-4,-3.470626964e-6,6.943764e-9,
		3.6780e-11,-1.36e-13};
	static real c2[] = {
		1.843740587300906e0,-0.076852840844786e0,
		1.271927136655e-3,-4.971736704e-6,-3.3126120e-8,
		2.42310e-10,-1.70e-13,-1.0e-15};

	xx=8.0*x*x-1.0;
	*gam1=chebev(-1.0,1.0,c1,NUSE1,xx);
	*gam2=chebev(-1.0,1.0,c2,NUSE2,xx);
	*gampl= *gam2-x*(*gam1);
	*gammi= *gam2+x*(*gam1);
}

float chebev(real a, real b, real c[], int m, real x)
{
	void nrerror(char error_text[]);
	float d=0.0,dd=0.0,sv,y,y2;
	int j;

	if ((x-a)*(x-b) > 0.0) nrerror("x not in range in routine chebev");
	y2=2.0*(y=(2.0*x-a-b)/(b-a));
	for (j=m-1;j>=1;j--) {
		sv=d;
		d=y2*d-dd+c[j];
		dd=sv;
	}
	return y*d-dd+0.5*c[0];
}

void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
//	void _exit();

	mastererr("Model run-time error...\n");
	mastererr("%s\n",error_text);
	mastererr("...now exiting to system...\n");
	exit(1);
}

#undef NUSE1
#undef NUSE2
#undef EPS
#undef FPMIN
#undef MAXIT
#undef XMIN
/* (C) Copr. 1986-92 Numerical Recipes Software ?421.1-9. */

