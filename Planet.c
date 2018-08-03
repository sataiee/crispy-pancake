/** \file Planet.c

Accretion of disk material onto the planets, and solver of planetary
orbital elements.  The prescription used for the accretion is the one
designed by W. Kley.

*/

#include "mp.h"

void AccreteOntoPlanets (Rho, Vrad, Vtheta, dt, sys)
     real dt;
     PolarGrid *Rho, *Vrad, *Vtheta;
     PlanetarySystem *sys;
{
  real RRoche, Rplanet, distance, dx, dy, deltaM, angle, temp;
  int i_min,i_max, j_min, j_max, i, j, l, jf, ns, nr, lip, ljp, k;
  real Xplanet, Yplanet, Mplanet, VXplanet, VYplanet;
  real facc, facc1, facc2, frac1, frac2; /* We adopt the same notations as W. Kley */
  real *dens, *abs, *ord, *vrad, *vtheta;
  real PxPlanet, PyPlanet, vrcell, vtcell, vxcell, vycell, xc, yc;
  real dPxPlanet, dPyPlanet, dMplanet=0;
  real mass45, mass4575, mtot, meandr, dM, massfloor;
  int ip;
  extern boolean FargoPlanete, OneDRun;
  nr     = Rho->Nrad;
  ns     = Rho->Nsec;
  dens   = Rho->Field;
  abs    = CellAbscissa->Field;
  ord    = CellOrdinate->Field;
  vrad   = Vrad->Field;
  vtheta = Vtheta->Field;
  if (!OneDRun){
  for (k=0; k < sys->nb; k++) {
    Xplanet = sys->x[k];
    Yplanet = sys->y[k];
    Mplanet = sys->mass[k];
    Rplanet = sqrt(Xplanet*Xplanet+Yplanet*Yplanet);
    RRoche = pow(1.0/3.0*Mplanet,1.0/3.0)*Rplanet; 
    /* Central mass is 1.0 */
    ip = ReturnIndex(Rplanet);
    meandr = (Radii[ip+2]-Radii[ip-2])/4.;
    if ((RRoche*0.45/meandr) >= 2.0 ){
        if (FargoPlanete){
            dM = CalculateFracMass(Rho,sys,k,0.45,1);
            MPI_Allreduce (&dM, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            mass45 = temp;
            dM = CalculateFracMass(Rho,sys,k,0.75,1);
            MPI_Allreduce (&dM, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            mass4575 = temp;
            mass4575 -= mass45;
            dM = CalculateFracMass(Rho,sys,k,0.75,0);
            MPI_Allreduce (&dM, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            massfloor = temp;
            sys->acc[k] = MdotEnvelope[k] / (mass45+1./3.*mass4575-massfloor);
            if (((sys->acc[k]*dt) > 1.0) || (MenvCount[k] > 100)){
               if (MenvCount[k] > 100)
                 mastererr("Planet %d is detached.\n", k);
               else {
                 mastererr("Accretion rate on the planet %d is very high, it should be about %e. we set it to the maximum value 1/dt\n",k, (mass45+1./3.*mass4575)/dt);
                 MenvCount[k] +=1;
               }
               if (dt > (0.8*pow(Rplanet,1.5)))
                 sys->acc[k] = 1./dt;
               else 
                 sys->acc[k] = 1./0.8/pow(Rplanet,1.5);
             }
        }
      if (sys->acc[k] > 1e-10) {
          dMplanet = dPxPlanet = dPyPlanet = 0.0;
          /* Hereafter : initialization of W. Kley's parameters */
          facc = dt*(sys->acc[k]);
          facc1 = 1.0/3.0*facc;
          facc2 = 2.0/3.0*facc;
          frac1 = 0.75;
          frac2 = 0.45;
          /* W. Kley's parameters initialization finished */
          VXplanet = sys->vx[k];
          VYplanet = sys->vy[k];
          
          i_min=0;
          i_max=nr-1;
          while ((Rsup[i_min] < Rplanet-RRoche) && (i_min < nr)) i_min++;
          while ((Rinf[i_max] > Rplanet+RRoche) && (i_max > 0)) i_max--;
          angle = atan2 (Yplanet, Xplanet);
          j_min =(int)((real)ns/(PMAX-PMIN)*(angle - 2.0*RRoche/Rplanet));
          j_max =(int)((real)ns/(PMAX-PMIN)*(angle + 2.0*RRoche/Rplanet));
          PxPlanet = Mplanet*VXplanet;
          PyPlanet = Mplanet*VYplanet;
    #pragma omp parallel for private(j,jf,vrcell,vtcell,vxcell,vycell,l,lip,ljp,xc,yc,dx,dy,distance,deltaM) shared(dPxPlanet, dPyPlanet, dMplanet)
          for (i = i_min; i <= i_max; i++) {
        for (j = j_min; j <= j_max; j++) {
          jf = j;
          while (jf <  0)  jf += ns;
          while (jf >= ns) jf -= ns;
          l   = jf+i*ns;
          lip = l+ns;
          ljp = l+1;
          if (jf == ns-1) ljp = i*ns;
          xc = abs[l];
          yc = ord[l];
          dx = Xplanet-xc;
          dy = Yplanet-yc;
          distance = sqrt(dx*dx+dy*dy);
          vtcell=0.5*(vtheta[l]+vtheta[ljp])+Rmed[i]*OmegaFrame;
          vrcell=0.5*(vrad[l]+vrad[lip]);
          vxcell=(vrcell*xc-vtcell*yc)/Rmed[i];
          vycell=(vrcell*yc+vtcell*xc)/Rmed[i];
            if (distance < frac1*RRoche) {
//Sareh modified this part and added floordens to exclude the floor density mass
              if (dens[l]>floordens)
                deltaM = facc1*(dens[l]-floordens)*Surf[i]; 
              else
                deltaM = 0.0;
              if (i < Zero_or_active) deltaM = 0.0;
              if (i >= Max_or_active) deltaM = 0.0;
              dens[l] *= (1.0 - facc1);
    #pragma omp atomic
              dPxPlanet    += deltaM*vxcell;
    #pragma omp atomic
              dPyPlanet    += deltaM*vycell;
    #pragma omp atomic
              dMplanet     += deltaM;
            }
            if (distance < frac2*RRoche) {
//Sareh modified this part and added floordens to exclude the floor density mass
              if (dens[l]>floordens)
                deltaM = facc2*(dens[l]-floordens)*Surf[i]; 
              else
                deltaM = 0.0;
              if (i < Zero_or_active) deltaM = 0.0;
              if (i >= Max_or_active) deltaM = 0.0;
              dens[l] *= (1.0 - facc2);
    #pragma omp atomic
              dPxPlanet    += deltaM*vxcell;
    #pragma omp atomic
              dPyPlanet    += deltaM*vycell;
    #pragma omp atomic
              dMplanet     += deltaM;
            }
           }
         }
          MPI_Allreduce (&dMplanet, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
          dMplanet = temp;
          MPI_Allreduce (&dPxPlanet, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
          dPxPlanet = temp;
          MPI_Allreduce (&dPyPlanet, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
          dPyPlanet = temp;
          PxPlanet += dPxPlanet;
          PyPlanet += dPyPlanet;
          Mplanet  += dMplanet;
          Menvelope[k] += dMplanet;
          MenvRemoved[k] += dMplanet;
          MenvAccreted[k] += (mass45+1./3.*mass4575)- MenvRemained[k];
          MenvRemained[k] = (mass45+1./3.*mass4575)-dMplanet;
          if (sys->FeelDisk[k] == YES) {
        sys->vx[k] = PxPlanet/Mplanet;
        sys->vy[k] = PyPlanet/Mplanet;
          }
          sys->mass[k] = Mplanet;
          AccMassPls += dMplanet;
      }
    } 
  }
 }
}


void FindOrbitalElements (sys, x, y, vx, vy, m, n, l)
     PlanetarySystem *sys;
     real x, y, vx, vy, m; // m is mstar + mplanet
     int n, l;
{
  int nbplanets;
  real Ax, Ay, e, d, h, a, E, M, V, lambda, varpi;
  real PerihelionPA, MenvelopeMax;
  FILE *output, *outputplanete;
  char name[256], nameplanete[256];
  extern boolean FargoPlanete;
  if (CPU_Rank != CPU_Highest) return;
  sprintf (name, "%sorbit%d.dat", OUTPUTDIR, n);
  output = fopen (name, "a");
  if (output == NULL) {
    message ("Can't open 'orbit%d.dat'. Exited.\n",n);
    prs_exit (1);
  }
  nbplanets = sys->nb;
  h = x*vy-y*vx;
  d = sqrt(x*x+y*y);
  Ax = x*vy*vy-y*vx*vy-G*m*x/d;
  Ay = y*vx*vx-x*vx*vy-G*m*y/d;
  e = sqrt(Ax*Ax+Ay*Ay)/G/m;
  a = h*h/G/m/(1-e*e);
  sys->a[n] = a;
  sys->e[n] = e;
  // E: eccentric anomaly
  if (e != 0.0) {
    E = acos((1.0-d/a)/e);
  } else {
    E = 0.0;
  }
  if ((x*y*(vy*vy-vx*vx)+vx*vy*(x*x-y*y)) < 0) E= -E;
  // M: mean anomaly
  M = E-e*sin(E);
  if (e != 0.0) {
    V = acos ((a*(1.0-e*e)/d-1.0)/e);
  } else {
    V = 0.0;
  }
  if (E < 0.0) V = -V;
  // Mean longitude
  lambda = atan2(y,x);
  // Periastron position angle (definition 1)
  if (e != 0.0) {
    PerihelionPA=atan2(Ay,Ax);
  } else {
    PerihelionPA=atan2(y,x);
  }
  // Periastron longitude, definition 2
  if (e != 0.0)
    varpi = 2.0*atan2(Ay,1.+Ax);
  else
    varpi = lambda;
  fprintf (output, "%.12g\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g\n", PhysicalTime, e, a, M, V, PerihelionPA, lambda, varpi);
  fclose (output);
  /* NEW (Jan 2014): option to write orbit_lightn.dat files at every
     end of Fargo runs. Project Fargo/planet coupling */
  if ( (l==1) && (FargoPlanete) ) {
    sprintf (nameplanete, "%sfargo_to_planete%d.dat", OUTPUTDIR, n);
    outputplanete = fopen (nameplanete, "w");
    if (outputplanete == NULL) {
      message ("Can't open 'fargo_to_planete%d.dat'. Exited.\n",n);
      prs_exit (1);
    }
    /* Time counter, Physical time, Number of Planets, Planet's mass,
       semi-major axis, eccentricity, mean longitude, and pericentre
       position angle, maximum available mass for the planet for accretion 
       (include the whole removed mass+the accreted mass+ what is left the Hill radius from start of Runtime)
     */
    MenvelopeMax = MenvRemoved[n]+MenvAccreted[n]+MenvRemained[n];
    fprintf (outputplanete, "%d\t%g\t%d\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%d\n", TimeStep,PhysicalTime,nbplanets,m-1.0,a,e,lambda,PerihelionPA,MenvelopeMax, Menvelope[n], MenvCount[n]);
    fclose (outputplanete);
  }
}

real CalculateFracMass(Rho,sys,k,frac, index)
       PolarGrid *Rho;
       PlanetarySystem *sys;
       int k, index;
       real frac;
{
  real RRoche, Rplanet, distance, dx, dy, deltaM, angle, temp;
  int i_min,i_max, j_min, j_max, i, j, l, jf, ns, nr, lip, ljp, ip;
  real Xplanet, Yplanet, Mplanet, VXplanet, VYplanet;
  real *dens, *abs, *ord, *vrad, *vtheta;
  real xc, yc, dM;
  nr     = Rho->Nrad;
  ns     = Rho->Nsec;
  dens   = Rho->Field;
  abs    = CellAbscissa->Field;
  ord    = CellOrdinate->Field;
  dM= 0.0;
  Xplanet = sys->x[k];
  Yplanet = sys->y[k];
  Mplanet = sys->mass[k];
  Rplanet = sqrt(Xplanet*Xplanet+Yplanet*Yplanet);
  RRoche = pow((1.0/3.0*Mplanet),(1.0/3.0))*Rplanet; 
  ip = ReturnIndex(Rplanet);
  /* Central mass is 1.0 */
  i_min=0;
  i_max=nr-1;
  while ((Rsup[i_min] < Rplanet-RRoche) && (i_min < nr)) i_min++;
  while ((Rinf[i_max] > Rplanet+RRoche) && (i_max > 0)) i_max--;
  angle = atan2 (Yplanet, Xplanet);
  j_min =(int)((real)ns/(PMAX-PMIN)*(angle - 2.0*RRoche/Rplanet));
  j_max =(int)((real)ns/(PMAX-PMIN)*(angle + 2.0*RRoche/Rplanet));
#pragma omp parallel for private(j,jf,vrcell,vtcell,vxcell,vycell,l,lip,ljp,xc,yc,dx,dy,distance,deltaM) shared(dPxPlanet, dPyPlanet, dMplanet)
  for (i = i_min; i <= i_max; i++) {
    for (j = j_min; j <= j_max; j++) {
       jf = j;
       while (jf <  0)  jf += ns;
       while (jf >= ns) jf -= ns;
       l   = jf+i*ns;
       lip = l+ns;
       ljp = l+1;
       if (jf == ns-1) ljp = i*ns;
       xc = abs[l];
       yc = ord[l];
       dx = Xplanet-xc;
       dy = Yplanet-yc;
       distance = sqrt(dx*dx+dy*dy);
       if (distance < frac*RRoche) {
         if (index == 1)
           deltaM = dens[l]*Surf[i];
         else
           deltaM = floordens*Surf[i];
         if (i < Zero_or_active) deltaM = 0.0;
         if (i >= Max_or_active) deltaM = 0.0;
         dM += deltaM;
       }
    }
  }
  return dM;
}
