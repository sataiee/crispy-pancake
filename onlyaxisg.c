/** \file onlyaxisg.c
 * 
 * This file contains the functions for calculating the axi-symetric part of 
 * the disc self-gravity as in Kley1996. An array that contains the elliptical
 * integrals are calculated here and passed to the FillForcesArrays() to apply the 
 * sg-acceleration of the azimuthally averaged density on each cell */
 
#include "mp.h"
#include <stdarg.h>

real rf(real x, real y, real z)
{
  real ERRTOL= 0.08, TINY = 1.5e-38, BIG = 3.0e37, THIRD = 1.0/3.0;
  real C1 = (1.0/24.0), C2  = 0.1, C3 = (3.0/44.0), C4 = (1.0/14.0);
  
	real alamb,ave,delx,dely,delz,e2,e3,sqrtx,sqrty,sqrtz,xt,yt,zt;

	if (fmin(fmin(x,y),z) < 0.0 || fmin(fmin(x+y,x+z),y+z) < TINY || fmax(fmax(x,y),z) > BIG)
			nrerror("invalid arguments in rf");
	xt=x;
	yt=y;
	zt=z;
	do {
		sqrtx=sqrt(xt);
		sqrty=sqrt(yt);
		sqrtz=sqrt(zt);
		alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz;
		xt=0.25*(xt+alamb);
		yt=0.25*(yt+alamb);
		zt=0.25*(zt+alamb);
		ave=THIRD*(xt+yt+zt);
		delx=(ave-xt)/ave;
		dely=(ave-yt)/ave;
		delz=(ave-zt)/ave;
	} while (fmax(fmax(fabs(delx),fabs(dely)),fabs(delz)) > ERRTOL);
	e2=delx*dely-delz*delz;
	e3=delx*dely*delz;
	return (1.0+(C1*e2-C2-C3*e3)*e2+C4*e3)/sqrt(ave);
}

real rd(real x, real y, real z)
{
  real ERRTOL = 0.05, TINY = 1.0e-25, BIG = 4.5e21;
  real C1 = (3.0/14.0), C2 = (1.0/6.0), C3 = (9.0/22.0), C4 = (3.0/26.0);
  real C5 = (0.25*C3), C6 = (1.5*C4);
  
	real alamb,ave,delx,dely,delz,ea,eb,ec,ed,ee,fac,sqrtx,sqrty,
		sqrtz,sum,xt,yt,zt;

	if (fmin(x,y) < 0.0 || fmin(x+y,z) < TINY || fmax(fmax(x,y),z) > BIG)
		nrerror("invalid arguments in rd");
	xt=x;
	yt=y;
	zt=z;
	sum=0.0;
	fac=1.0;
	do {
		sqrtx=sqrt(xt);
		sqrty=sqrt(yt);
		sqrtz=sqrt(zt);
		alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz;
		sum += fac/(sqrtz*(zt+alamb));
		fac=0.25*fac;
		xt=0.25*(xt+alamb);
		yt=0.25*(yt+alamb);
		zt=0.25*(zt+alamb);
		ave=0.2*(xt+yt+3.0*zt);
		delx=(ave-xt)/ave;
		dely=(ave-yt)/ave;
		delz=(ave-zt)/ave;
	} while (fmax(fmax(fabs(delx),fabs(dely)),fabs(delz)) > ERRTOL);
	ea=delx*dely;
	eb=delz*delz;
	ec=ea-eb;
	ed=ea-6.0*eb;
	ee=ed+ec+ec;
	return 3.0*sum+fac*(1.0+ed*(-C1+C5*ed-C6*delz*ee)
		+delz*(C2*ee+delz*(-C3*ec+delz*C4*ea)))/(ave*sqrt(ave));
}

real ellf(real phi, real ak)
{
	real s;

	s=sin(phi);
	return s*rf(cos(phi)*cos(phi),(1.0-s*ak)*(1.0+s*ak),1.0);
}

real elle(real phi, real ak)
{
	real rd(real x, real y, real z);
	real rf(real x, real y, real z);
	real cc,q,s;

	s=sin(phi);
	cc=cos(phi) * cos(phi);
	q=(1.0-s*ak)*(1.0+s*ak);
	return s*(rf(cc,q,1.0)-(s*ak*s*ak)*rd(cc,q,1.0)/3.0);
} 

/* This function calculates the constant part of axi-symmetric 
 * SG acceleration of the disc on its potential (Kley 1996).
 * The name of the arrays are chosen consistenly with the paper*/

void CalculateAxiSGDiskPotentialTools(SGAarray)
  real *SGAarray;
{
  int i, j, k, nr, l;
  MPI_Request req1;
  masterprint("Calculating the tools for axi-symmetric SG acceleration , please be patient...\n");
  // local variables
  nr = NRAD+1;
  real *Aarray;
  Aarray = (real *)malloc(sizeof(real) * NRAD*(GLOBALNRAD+1));
  
  // Calculate A array for each processor
  CalculateAs(Aarray);
  
  // Make the global A arrays 
  // (because the size of arrays are not the same for all processors, I cannot use Allgather)
  if ( CPU_Number == 1 ) {
    for ( i = 0; i < GLOBALNRAD; i++){
      for (j = 0; j <= GLOBALNRAD; j++){
        l = i*(GLOBALNRAD+1) + j;
        SGAarray[l] = Aarray[l];
      }
    }
    printf("\n");
  }
  if ( CPU_Number > 1 ) {
    if ( CPU_Rank == 0 ) {
      for (i = 0; i < Max_or_active; i++){
        for (j = 0; j <= GLOBALNRAD; j++){
          l = i*(GLOBALNRAD+1) + j;
          SGAarray[l] = Aarray[l];
        }
      }
      MPI_Isend (SGAarray, (GLOBALNRAD)*(GLOBALNRAD+1), MPI_DOUBLE, CPU_Next, 0, MPI_COMM_WORLD, &req1);
      MPI_Wait (&req1, &fargostat);
    }
    if ( CPU_Rank != 0 ) {
      MPI_Irecv (SGAarray, (GLOBALNRAD)*(GLOBALNRAD+1), MPI_DOUBLE, CPU_Prev, 0, MPI_COMM_WORLD, &req1);
      MPI_Wait (&req1, &fargostat);
      for (i = Zero_or_active; i < Max_or_active; i++){
        for (j = 0; j <= GLOBALNRAD; j++){
          l = (i+IMIN)*(GLOBALNRAD+1) + j ;
          k = i*(GLOBALNRAD+1) + j;
          SGAarray[l] = Aarray[k];
        }
      }
      if ( CPU_Rank != CPU_Highest ) {
        MPI_Isend (SGAarray, (GLOBALNRAD)*(GLOBALNRAD+1), MPI_DOUBLE, CPU_Next, 0, MPI_COMM_WORLD, &req1);
        MPI_Wait (&req1, &fargostat);
      } 
    }
    MPI_Barrier (MPI_COMM_WORLD);
    MPI_Bcast (SGAarray, (GLOBALNRAD)*(GLOBALNRAD+1), MPI_DOUBLE, CPU_Highest, MPI_COMM_WORLD);
  }
  free(Aarray);
}

void CalculateAs(Aarray)
  real *Aarray;
{
  int i,j,l;
  real k;
  for (i = 0; i < NRAD; i++){
    for (j = 0; j <= (i+IMIN); j++){
      l = i*(GLOBALNRAD+1) + j;
      k = Rmed[i]/Radii[j];
      Aarray[l] = k * elle(PI/2, 1./k) - (k*k-1)/k * ellf(PI/2., 1./k);
    }
    for (j = (i+1+IMIN); j <= GLOBALNRAD; j++){
      l = i*(GLOBALNRAD+1) + j;
      k = Rmed[i]/Radii[j];
      Aarray[l] = elle(PI/2., k);
    }
  }
}

/* These are the elliptical function the Sareh wrote the functions
 * one can compare the results with the ones from numerical recipe
 * that comes later 
 */
/*real elle(phi, k)
  real k, phi;
{
  int j, ns=0;
  real E = 0;
  while (azimuth[ns] < phi)
    ns++;
  for (j = 0; j < ns; j++)
    E += sqrt(1-k*k * sin(azimuth[j])*sin(azimuth[j]));
  return E*(azimuth[1]-azimuth[0]);
}

real ellf(phi, k)
  real k, phi;
{
  int j, ns=0;
  real F = 0;
  while (azimuth[ns] < phi)
    ns++;
  for (j = 0; j < ns; j++)
    F += 1./sqrt(1-k*k * sin(azimuth[j]*sin(azimuth[j])));
  return F*(azimuth[1]-azimuth[0]);
} */ 

