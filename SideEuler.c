/** \file SideEuler.c

Total mass and angular momentum monitoring, and  conditions.
In addition, this file contains a few low-level functions that
manipulate PolarGrid 's or initialize the forces evaluation.

*/

#include "mp.h"

extern boolean OpenInner, KNOpen, NonReflecting, OuterSourceMass, Evanescent, Alexboundary;
extern boolean SelfGravity, SGZeroMode, EnergyEquation, MixedBC;
extern Pair DiskOnPrimaryAcceleration;
extern int dimfxy;
real Hp0, Hg0, Ht0;
extern boolean AccBoundary, ThermalDiffusion, RadiativeDiffusion;

real GasTotalMass (array)
     PolarGrid *array;
{
  int i, j, ns;
  real *density, total = 0.0, fulltotal=0.0;
  ns = array->Nsec;
  density = array->Field;
  if (FakeSequential && (CPU_Rank > 0)) 
    MPI_Recv (&total, 1, MPI_DOUBLE, CPU_Rank-1, 0, MPI_COMM_WORLD, &fargostat);
  for (i = Zero_or_active; i < Max_or_active; i++) {
    for (j = 0; j < ns; j++) {
      total += Surf[i]*density[j+i*ns];
    }
  }
  if (FakeSequential) {
    if (CPU_Rank < CPU_Number-1)
      MPI_Send (&total, 1, MPI_DOUBLE, CPU_Rank+1, 0, MPI_COMM_WORLD);
  } else {
    MPI_Allreduce (&total, &fulltotal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  }
  if (FakeSequential) {
    MPI_Bcast (&total, 1, MPI_DOUBLE, CPU_Number-1, MPI_COMM_WORLD);
    fulltotal = total;
  }
  return fulltotal;
}

real GasMomentum (Density, Vtheta)
     PolarGrid *Density, *Vtheta;
{
  int i,j,l,ns;
  real vt_cent;
  real *density, *vtheta, total = 0.0, fulltotal=0.0;
  ns = Density->Nsec;
  density = Density->Field;
  vtheta = Vtheta->Field;
  if (FakeSequential && (CPU_Rank > 0)) 
    MPI_Recv (&total, 1, MPI_DOUBLE, CPU_Rank-1, 2, MPI_COMM_WORLD, &fargostat);
  for (i = Zero_or_active; i < Max_or_active; i++) {
    for (j = 0; j < ns; j++) {
      l = i*ns + j;
      /* centered-in-cell azimuthal velocity */
      if (j < ns-1)
       vt_cent = 0.5*(vtheta[l]+vtheta[l+1]) + Rmed[i]*OmegaFrame;
      else
       vt_cent = 0.5*(vtheta[l]+vtheta[i*ns]) + Rmed[i]*OmegaFrame;
      total += Surf[i]*density[l]*Rmed[i]*vt_cent;
    }
  }
  if (FakeSequential) {
    if (CPU_Rank < CPU_Number-1)
      MPI_Send (&total, 1, MPI_DOUBLE, CPU_Rank+1, 2, MPI_COMM_WORLD);
  }
  else
    MPI_Allreduce (&total, &fulltotal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  if (FakeSequential) {
    MPI_Bcast (&total, 1, MPI_DOUBLE, CPU_Number-1, MPI_COMM_WORLD);
    fulltotal = total;
  }
  return fulltotal;
}

real GasTotalEnergy (Density, Vrad, Vtheta, Energy)
     PolarGrid *Density, *Vrad, *Vtheta, *Energy;
{
  int i, j, l, ns;
  real *density, *vrad, *vtheta, *energy, *pot;
  real vr_cent, vt_cent;
  real total = 0.0, fulltotal=0.0;
  ns = Density->Nsec;
  density = Density->Field;
  vrad = Vrad->Field;
  vtheta = Vtheta->Field;
  energy = Energy->Field;
  pot = Potential->Field;
  if (FakeSequential && (CPU_Rank > 0)) 
    MPI_Recv (&total, 1, MPI_DOUBLE, CPU_Rank-1, 2, MPI_COMM_WORLD, &fargostat);
  for (i = Zero_or_active; i < Max_or_active; i++) {
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      /* centered-in-cell radial velocity */
      vr_cent = (Rmed[i]-Rinf[i])*vrad[l+ns] + (Rsup[i]-Rmed[i])*vrad[l];
      vr_cent /= (Rsup[i]-Rinf[i]);
      /* centered-in-cell azimuthal velocity */
      if (j < ns-1)
       vt_cent = 0.5*(vtheta[l]+vtheta[l+1]) + Rmed[i]*OmegaFrame;
      else
       vt_cent = 0.5*(vtheta[l]+vtheta[i*ns]) + Rmed[i]*OmegaFrame;
      total += 0.5*Surf[i]*density[l]*(vr_cent*vr_cent + vt_cent*vt_cent) + \
       Surf[i]*energy[l] -                                          \
       Surf[i]*density[l]*pot[l];
      /* Gas total energy is the sum of its kinematic energy, internal energy */
      /* and gravitational potential energy, including self-gravity */
    }
  }
  if (FakeSequential) {
    if (CPU_Rank < CPU_Number-1)
      MPI_Send (&total, 1, MPI_DOUBLE, CPU_Rank+1, 2, MPI_COMM_WORLD);
  } else {
    MPI_Allreduce (&total, &fulltotal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  }
  if (FakeSequential) {
    MPI_Bcast (&total, 1, MPI_DOUBLE, CPU_Number-1, MPI_COMM_WORLD);
    fulltotal = total;
  }
  return fulltotal;
}


void CheckMomentumConservation (Density, Vtheta, sys)
     PolarGrid *Density, *Vtheta;
     PlanetarySystem *sys;
{
  FILE *fichmom;
  char name[256];
  int k;
  real totalmomentum, plmom;
  real xplanet, yplanet, vxplanet, vyplanet;
  real rpl, thetapl, vazimpl, masspl;
  real gasmom, planetsmom;
  gasmom = GasMomentum (Density, Vtheta);
  planetsmom = 0.;
  
  for ( k = 0; k < sys->nb; k++ ) {
    xplanet     = sys->x[k];
    yplanet     = sys->y[k];
    rpl         = sqrt( xplanet*xplanet + yplanet*yplanet );
    thetapl     = atan2 (yplanet, xplanet);
    vxplanet    = sys->vx[k];
    vyplanet    = sys->vy[k];
    vazimpl     = -vxplanet*sin(thetapl) + vyplanet*cos(thetapl);
    masspl      = sys->mass[k];
    plmom       = masspl*rpl*vazimpl;
    planetsmom += plmom;
  }
  totalmomentum = gasmom + planetsmom;
  if ( PhysicalTime < 1e-10 ) {
    Hp0 = plmom;
    Hg0 = gasmom;
    Ht0 = totalmomentum;
    printf("time = %lg, Hp0 = %lg, Hg0 = %lg et Ht0 = %lg\n", PhysicalTime, Hp0, Hg0, Ht0);
  }
  if (!CPU_Master) return;
  sprintf (name, "%s%s.dat", OUTPUTDIR, "Momentum");
  fichmom = fopen(name, "a");
  if (fichmom == NULL) {
    fprintf (stderr, "Can't write 'Momentum.dat' file. Aborting.\n");
    prs_exit (1);
  }
  plmom = fabs (plmom - Hp0);
  gasmom = fabs (gasmom - Hg0);
  totalmomentum = fabs (totalmomentum - Ht0);
  fprintf (fichmom, "%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\n", PhysicalTime, plmom, gasmom, totalmomentum, totalmomentum / Ht0);
  fclose (fichmom);
}


void DivisePolarGrid (Num, Denom, Res)
     PolarGrid *Num, *Denom, *Res;
{
  int i,j,l,nr,ns;
  real *num, *denom, *res;
  num = Num->Field;
  denom=Denom->Field;
  res = Res->Field;
  ns = Res->Nrad;
  nr = Res->Nsec;
#pragma omp parallel for private(j,l)
  for (i = 0; i <= nr; i++) {
    for (j = 0; j < ns; j++) {
      l = j+ns*i;
      res[l] = num[l]/(denom[l]+1e-20);
    }
  }
}

void InitComputeAccel ()
{
  int i, j, l, nr, ns;
  real *abs, *ord;
  CellAbscissa = CreatePolarGrid (NRAD,NSEC,"abscissa");
  CellOrdinate = CreatePolarGrid (NRAD,NSEC,"ordinate");
  nr = CellAbscissa->Nrad;
  ns = CellAbscissa->Nsec;
  abs = CellAbscissa->Field;
  ord = CellOrdinate->Field;
  for (i = 0; i < nr; i++) {
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      abs[l] = Rmed[i] * cos(azimuth[j]);
      ord[l] = Rmed[i] * sin(azimuth[j]);
    }
  }
}
  
Pair ComputeAccel (force, Rho, x, y, rsmoothing, mass, psys, index)
     Force *force;
     PolarGrid *Rho;
     real x, y, rsmoothing, mass;
     PlanetarySystem *psys;
     int index;
{
  Pair acceleration;
  ComputeForce (force, Rho, x, y, rsmoothing, mass, dimfxy, psys, index);
  if (ExcludeHill) {
    acceleration.x = force->fx_ex_inner+force->fx_ex_outer;
    acceleration.y = force->fy_ex_inner+force->fy_ex_outer;
  } else {
    acceleration.x = force->fx_inner+force->fx_outer;
    acceleration.y = force->fy_inner+force->fy_outer;
  }
  return acceleration;
}

/* This function calculates the velocities using the torque formula
 * given in Paardekooper et. al 2011 and Paardekooper 2010. This is used only when Planete is
 * running with Fargo and mass of some planets are smaller than the critical values.
 * Another important point is applying the timescales on the torque calculation
 * before T_uturn, the torque is in linear regime,
 * before liberation time, the torque is unsaturated,
 * and after that it is equal to the value in the paper.
 * Sareh applied a smooth transition function on torques based on the time scale comparison */
Pair AccelFromFormula (force, Rho, x, y, smoothing, mass, sys, index, flag)
    Force *force;
    PolarGrid *Rho;
    real x,y, smoothing, mass;
    PlanetarySystem *sys;
    int index;
    int flag;
{
    real theta, r2, r, rdot, Gamma, GammaL, GammaHSB, GammaHSE, Gamma0, Gammaeff;
    real GammaSum, GammaLB, GammaLE, GammaCTotal, GammaUnsat, GammaLin, GammaPBK;
    int ip,i, n1, n2, lip;
    real lnr[GLOBALNRAD], lnf[GLOBALNRAD];
    real alphaf, betaf, dlogr, h, chi;
    real Q, denominator, gamq, gamq2, xs, K, Pnu=0, Pk=0, zeta, Rh, dist;
    real coeff[2];
    Pair accel;
    real TauMig, TauEcc, vdotr;
    real omega, new_r, denom, dtheta, vx, vy, gradT, Rl, lambda;
    real Tlib, TUturn, torbit;
    r2 = x*x+y*y;
    theta  = atan2(y,x);
    r = sqrt(r2);
    omega = sqrt((1.+mass)/r/r/r);
    ip = ReturnIndex(r);
    lip = ip+1;
    if (lip == GLOBALNRAD) lip = ip; 
    Rh = r*pow(mass/3.,1./3);
    dist = (Rh > (Radii[ip+1]-Radii[ip]) ? Rh : (Radii[ip+1]-Radii[ip]));
    n1 = ReturnIndex(r-dist);
    n2 = ReturnIndex(r+dist);
    if ((ip == -1) || (ip >= GLOBALNRAD) ){
        mastererr("Wrong index for planet position \n");
        mastererr("ip=%d, r=%lg \n", ip, r);
        prs_exit(1);
    }        
    if (EnergyEquation)
        h = axics[ip]/omega/r/sqrt(ADIABATICINDEX);
    else
        h = AspectRatio(r)* pow(r,FLARINGINDEX);
    /* torque scaling coeeficients, see section 3 and 5.6 of the paper */
    Gamma0 = mass/h/h * r * axidens[ip];
    for (i=0; i<GLOBALNRAD; i++){
        lnr[i] = log(Radii[i]);
        lnf[i] = log(axidens[i]);
    }
    LinearInterpole(lnf, lnr , coeff , n1, n2);
    alphaf = -coeff[0];       
    for (i=0; i<GLOBALNRAD; i++){
       lnf[i] = log(axitemp[i]);
    }
    LinearInterpole(lnf, lnr , coeff, n1, n2);
    betaf = -coeff[0];              
    K = sqrt(r) /(2 * PI * FViscosity(r, axics[i])) ;
    zeta = betaf - (ADIABATICINDEX-1)*alphaf;

    if (EnergyEquation){
       if (ThermalDiffusion) {
            Q = 2*DIFFUSIVITY/3/(h*h*h)/sqrt(r);
            gamq = ADIABATICINDEX*Q;
            gamq2 = gamq*gamq;
            denominator = sqrt(pow(gamq2+1,2)-16*(Q*Q)*sqrt(ADIABATICINDEX-1)) + gamq2 - 1;
            denominator *= 2;
            denominator = 0.5*sqrt(denominator);
            denominator += gamq;
            Gammaeff = 2*gamq/denominator;
            xs = 1.1/pow(Gammaeff,0.25)*pow(0.4/THICKNESSSMOOTHING,0.25)*sqrt(mass/h);
            Pk = sqrt(sqrt(r)*xs*xs*xs/2/PI/DIFFUSIVITY);
        } else if (RadiativeDiffusion) {
            gradT = sqrt(pow((axitemp[lip]-axitemp[ip])/(GlobalRmed[lip]-GlobalRmed[ip]),2));
            Rl = 8*h/axidens[ip]/opaaxi[ip]*gradT/axitemp[ip];
            if (Rl <= 2){
               lambda = 2./(3.*sqrt(9+10*Rl*Rl));
            } else {
               lambda = 10./sqrt(10.*Rl+9+sqrt(180.*Rl+81.));
            }
            chi  = lambda*64.*ADIABATICINDEX*(ADIABATICINDEX-1)*sigma_SB*pow(axitemp[ip],4)\
                     / (opaaxi[ip]*pow(axidens[ip],2)*pow(GlobalRmed[ip],-3));
            Q = 2*chi/3/(h*h*h)/sqrt(r);
            gamq = ADIABATICINDEX*Q;
            gamq2 = gamq*gamq;
            denominator = sqrt(pow(gamq2+1,2)-16*(Q*Q)*(ADIABATICINDEX-1)) + gamq2 - 1;
            denominator *= 2;
            denominator = 0.5*sqrt(denominator);
            denominator += gamq;
            Gammaeff = 2*gamq/denominator;
            xs = 1.1/pow(Gammaeff,0.25)*pow(0.4/THICKNESSSMOOTHING,0.25)*sqrt(mass/h);
            Pk = sqrt(sqrt(r)*xs*xs*xs/2/PI/chi);
        } else {
            Gammaeff = ADIABATICINDEX;
            xs = 1.1/pow(Gammaeff,0.25)*pow(0.4/THICKNESSSMOOTHING,0.25)*sqrt(mass/h);       
        }
    } else {  /*if locally isothermal */
        Gammaeff = 1;
        xs = 1.1/pow(Gammaeff,0.25)*pow(0.4/THICKNESSSMOOTHING,0.25)*sqrt(mass/h);
        Pk = 0;
    }
    Pnu = 2./3. * sqrt(xs*xs*xs*K);
    /* Barotropic vortensity related horseshoe drag Eq. 4 with modification from 2010 paper*/
    GammaHSB = 1.1 * (1.5-alphaf) * (0.4/THICKNESSSMOOTHING);
    /* Entropy related linear corotation torque Eq. 7 with modification from 2010 paper*/
    GammaLE = 2.2 * zeta * pow(0.4/THICKNESSSMOOTHING,0.71) - 1.4 * zeta/Gammaeff * pow(0.4/THICKNESSSMOOTHING,1.26);
    /* Barotropic vortensity related linear corotaion torque */
    GammaLB = 0.7 * (1.5-alphaf) * pow(0.4/THICKNESSSMOOTHING,1.26);
    /* Entropy related horseshoe drag, Eq. 5 with modifications from Paardekooper 2010 paper*/
    GammaHSE = zeta/Gammaeff *(0.4/THICKNESSSMOOTHING) * (10.1*sqrt(0.4/THICKNESSSMOOTHING)-2.2);
    /* Lindblad torque relation 47 */
    GammaL = (-2.5 - 1.7*betaf + 0.1*alphaf)* pow(0.4/THICKNESSSMOOTHING,0.71);
    if (EnergyEquation){
       if (ThermalDiffusion || RadiativeDiffusion){
           GammaCTotal = GammaHSE * Ffunc(Pnu)*Ffunc(Pk) * sqrt(Gfunc(Pnu)*Gfunc(Pk));
           GammaCTotal += GammaLE * sqrt((1-Kfunc(Pnu)) * (1-Kfunc(Pk)));
       } else {
           GammaCTotal = 0.0;
       }
    } else {  /*if locally isothermal */
       GammaCTotal = GammaHSE * Ffunc(Pnu)*Ffunc(Pk) * sqrt(Gfunc(Pnu)*Gfunc(Pk));
       GammaCTotal += GammaLE * sqrt((1-Kfunc(Pnu)) * (1-Kfunc(Pk)));
    }
    GammaCTotal += GammaHSB * Gfunc(Pnu) * Ffunc(Pnu);
    GammaCTotal += GammaLB * (1-Kfunc(Pnu));

    /* Applying the eccentricity correction on the corotation torque
     * using the relation 8 and 10 of Fendyke 2013 */
    GammaCTotal *= exp(-1 * sys->e[index] / (0.5*h+0.01)); 

    /* Comparing the timescales */
/*        Tlib = 8*PI*r / (3*omega*xs);
        TUturn = 1.14*h*Tlib / sqrt(Gammaeff);
        torbit = PhysicalTime/2/PI;
        GammaLin = GammaL + GammaLE + GammaLB;*/
    GammaPBK = GammaL + GammaCTotal; 
    if (EnergyEquation){
        GammaUnsat = GammaL + GammaHSE + GammaHSB;
    } else {
        GammaUnsat = GammaL + GammaHSB + GammaLE;
    }
/*        if (torbit <= TUturn){
          Gamma = GammaLin * pow(1-torbit/TUturn,2) + GammaUnsat * pow(torbit/TUturn,2);
        } else if ((torbit > TUturn) && (torbit <= Tlib)){
          Gamma = GammaUnsat * exp(10*(1-torbit/Tlib)) + GammaPBK * exp(10*(1-Tlib/torbit));
          Gamma /= (exp(10*(1-torbit/Tlib)) + exp(10*(1-Tlib/torbit)));
        } else {*/
    Gamma = GammaPBK;
//        }

    Gamma *= Gamma0/Gammaeff;
  /* In fargo torque is special torque that means it is devided 
       by the mass of the planet */  
               
//       rdot = 2* pow(r,0.5) * Gamma; //Be careful: devision to 2 is only for testing purposes.
        
    if (flag == 1){
        if (sys->FeelDisk[index] == YES){
            TauMig = - r*r * omega / Gamma; //Note that migration timescale is twice of orbital decay timescale.
            TauEcc =  2.6 * h*h * r*r * omega / Gamma0;
            vdotr = sys->vx[index]*sys->x[index]+sys->vy[index]*sys->y[index];
            accel.x =  -1. * sys->vx[index]/TauMig;
            accel.x +=  -2. /r/r * vdotr * sys->x[index]/TauEcc;
            accel.y =  -1. * sys->vy[index]/TauMig;
            accel.y +=  -2. /r/r * vdotr * sys->y[index]/TauEcc;
        }
    }       
    force->fy_inner = Gamma;
    force->fy_outer = Gamma0;
    force->fy_ex_inner = GammaL;
    force->fy_ex_outer = GammaCTotal;
    force->fx_inner = GammaLB;
    force->fx_outer = GammaLE;
    force->fx_ex_inner = GammaHSB;
    force->fx_ex_outer = GammaHSE;
      
    return accel;
}

int ReturnIndex(r)
    real r;
{
    int index=0;
    while ( index < GLOBALNRAD &&  Radii[index] < r) index++;
    return (index > GLOBALNRAD ? -1 : index);
}

real Ffunc(p)
    real p;
{
    return 1./(1+p*p/1.3/1.3);
}

real Gfunc(p)
    real p;
{
    real crit;
    crit = 8./45/PI;
    if (p < sqrt(crit))
        return 16./25.*pow(crit,-3./4)*pow(p,1.5);
    else 
        return 1-9./25.*pow(crit,4./3)*pow(p,-8./3);
}

real Kfunc(p)
    real p;
{
    real crit;
    crit = 28./45./PI;
    if (p < sqrt(crit))
        return 16./25.*pow(crit,-3./4)*pow(p,1.5);
    else 
        return 1-9./25.*pow(crit,4./3)*pow(p,-8./3);
}

void OpenBoundary (Vrad, Vtheta, Rho, Energy)
     PolarGrid *Vrad, *Vtheta, *Rho, *Energy;
{
  int i,j,l,ns,nr;
  extern boolean DontApplySubKeplerian;
  real *rho, *vr, *vt, *energy, *cs, visc;
  cs = SoundSpeed->Field;
  ns = Rho->Nsec;
  nr = Rho->Nrad;
  rho = Rho->Field;
  vr  = Vrad->Field;
  vt  = Vtheta->Field;
  energy = Energy->Field;
  /* -------------------------------- */
  /* Inner Boundary Condition         */
  /* -------------------------------- */
  if (CPU_Rank == 0) {
    i = 1;
#pragma omp parallel for private(l)
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      if ( KNOpen ) {
        /* Kley and Nelson (2008) prescription */
        rho[l-ns] = rho[l] ;      // zero gradient for surface density
        energy[l-ns] = energy[l]; // zero gradient for thermal energy
        visc = FViscosity(Radii[i], axics[i]); // azimuthally averaged viscosity
        vr[l] = -1.5*visc/Radii[i];
        visc = FViscosity(Radii[i-1], axics[i-1]);
        vr[l-ns] = -1.5*visc/Radii[i-1];
      } else {
        /* Standard outflow prescription */
        rho[l-ns] = rho[l] ;      // zero gradient for surface density
        energy[l-ns] = energy[l]; // zero gradient for thermal energy
        /* NEW (Sept 28 2011): if subkeplerian BC is not applied, then
        impose zero gradient on vtheta */
        if (DontApplySubKeplerian)
          vt[l-ns] = vt[l];
        if (vr[l+ns] >= 0.0){
          vr[l] = 0.0;            // vr set to zero when directed outward
        } else {
          vr[l-ns] = vr[l];       // vr extrapolated otherwise 
        }
      }
    }
  }
  /* -------------------------------- */
  /* Outer boundary condition         */
  /* -------------------------------- */
  if ( (CPU_Rank == CPU_Highest) && (!MixedBC) ) {
    i = nr-1;
#pragma omp parallel for private(l)
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      rho[l] = rho[l-ns];       // zero gradient for surface density
      energy[l] = energy[l-ns]; // zero gradient for thermal energy
      /* NEW (Sept 28 2011): if subkeplerian BC is not applied, then
       impose zero gradient on vtheta */
     if ( KNOpen ) {
      /* Kley and Nelson (2008) prescription */
      rho[l] = rho[l-ns] ;      // zero gradient for surface density
      energy[l] = energy[l-ns]; // zero gradient for thermal energy
      visc = FViscosity(Radii[i+IMIN], axics[i+IMIN]); // azimuthally averaged viscosity
      vr[l] = -1.5*visc/Radii[i+IMIN];
      rho[l+ns] = rho[l] ;      // zero gradient for surface density
      energy[l+ns] = energy[l]; // zero gradient for thermal energy
      vr[l+ns] = vr[l];
      } else {
        if (DontApplySubKeplerian)
          vt[l] = vt[l-ns];
        if (vr[l-ns] < 0.0 )
          vr[l] = 0.0;            // vr set to zero when directed inward
        else
          vr[l] = vr[l-ns];       // vr extrapolated otherwise 
        vr[l+ns] = vr[l];
        rho[l+ns] = rho[l] ;      // zero gradient for surface density
        energy[l+ns] = energy[l]; // zero gradient for thermal energy
     }
    }
  }
}

/* This boundary is written for the code when we use it with Planete
 * It has an imposed time dependent accretion at the outer and inner boundaries
 * in ordet to establish a constant accretion rate in the disc  */
void AccretingBoundary (Vrad, Vtheta, Rho, Energy, step)
     PolarGrid *Vrad, *Vtheta, *Rho, *Energy;
     real step;
{
  int i,j,l,ns,nr;
  real *dens, *vr, *vt, *energy, *cs, visco, visci, *temperature;
  real mdot, h, omega, vr_med, axi[GLOBALNRAD], timeyear;
  extern boolean MdotHartmann, DecInner, OpInner;
  cs = SoundSpeed->Field;
  ns = Rho->Nsec;
  nr = Rho->Nrad;
  dens = Rho->Field;
  vr  = Vrad->Field;
  vt  = Vtheta->Field;
  energy = Energy->Field;
  temperature = Temperature->Field;
  if (MdotHartmann){
    timeyear = (PhysicalTime*unit_time)/31556926.0; //time in year
    mdot = 1e-8 * pow((timeyear + THARTMANN)/1e6, -1.4); //Hartmann1998, modified to give a smaller Mdot because our disc is evolved
    mdot *= -(1.9891e30/31556926.0 / unit_mass*unit_time); //convert to code unit
  } else {
    if ((PhysicalTime/2./PI) <= MDOTTIME){
      mdot = MDOTINIT + (MDOTFINAL - MDOTINIT) * (PhysicalTime/2./PI) / MDOTTIME;
      mdot /= -1;
    } else {
      mdot = -MDOTFINAL;
    }
  }
  /*----------------------------------------------*/
  /* Inner boundary: three options exist:
       1) open outflow boundary (OpInner)
       2) decrease the inner surface density progressively (DecInner)
       3) damp to viscous profiles */
  /*----------------------------------------------*/
  if (CPU_Rank == 0) {
    i = 1;
#pragma omp parallel for private(l)
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      /*This boundary has a discontinuity in T and is very slow if used with radiative diffusion */
      if (DecInner){
          if (PhysicalTime != 0.0)
            dens[l] *= 0.999;
         dens[l-ns] = dens[l];
         vr[l-ns] = vr[l];
      } else if (OpInner){
      /* Standard outflow prescription */
          if (vr[l+ns] > 0)
            vr[l] = 0;
          else
            vr[l] = vr[l+ns];
          dens[l-ns] = dens[l];
             vr[l-ns] = vr[l];
          energy[l-ns] = energy[l];
     } else {
      /* Accreting value is imposed at the boundaries
         if energy equation is used, all terms must be included and 
         the temperature profiles is calculated using Bitsch+2014  */
           if (ViscosityAlpha){
              visco = ALPHAVISCOSITY*axics[i+IMIN] \
                  * axics[i+IMIN] / sqrt(ADIABATICINDEX) * pow(Rmed[i], 1.5);
              visci = ALPHAVISCOSITY*ALPHAVISCOSITY*axics[i-1+IMIN] \
                  * axics[i-1+IMIN]/ sqrt(ADIABATICINDEX) * pow(Rmed[i-1], 1.5);
           } else {
              visco = VISCOSITY;
              visci = VISCOSITY;
           }
         vr[l] = -3*(visco+visci)*0.5/2/Rinf[i];
         dens[l] = mdot/2./PI/Rmed[i]/vr[l];
         energy[l] = temperature[l]*dens[l]/(ADIABATICINDEX-1);
         vr[l-ns] = -3*visci/2/Rinf[i-1];
         dens[l-ns] = mdot/2./PI/Rmed[i-1]/vr[l-ns];
         energy[l] = temperature[l-ns]*dens[l-ns]/(ADIABATICINDEX-1);
     }
     vt[l-ns] = vt[l]*sqrt(Rmed[i]/Rmed[i-1]);
   }
  }
  /* -------------------------------- */
  /* Outer boundary condition: set to stable accretion values */
  /* -------------------------------- */
  if (CPU_Rank == CPU_Highest) {
    i = nr-1;
#pragma omp parallel for private(l)
    for (j = 0; j < ns; j++) {
      l = j+i*ns;   
      vr[l] = VradMed[i+IMIN];
      dens[l] = SigmaMed[i];
      vr[l-ns] = VradMed[i-1+IMIN];
      dens[l] = SigmaMed[i-1];
    }
  }
}

void NonReflectingBoundary (Vrad, Rho, Energy, Vtheta)
     PolarGrid *Vrad, *Rho, *Energy, *Vtheta;
{
  int i, j, l, ns, nr, jp, lp, i_angle;
  real *rho, *vr, *vt, *cs, *energy;
  real dangle, mean;
  real cs0, cs1, vr_med, csnrm1, csnrm2;
  real csinit0, csinit1;
  real omega;
  cs = SoundSpeed->Field;
  energy = Energy->Field;
  ns = Rho->Nsec;
  nr = Rho->Nrad;
  rho = Rho->Field;
  vr  = Vrad->Field;
  vt = Vtheta->Field;
  /* ======================================= */
  /* INNER NON-REFLECTING BOUNDARY CONDITION */
  /* ======================================= */
  if (CPU_Rank == 0) {
    i=1;
    /* ------------------------ */
    /* STANDARD BAROTROPIC CASE */
    /* ------------------------ */
    if (!EnergyEquation) {
      /* (a) sound speed in first two rings */
      cs0 = 0.0;
      cs1 = 0.0;
      for (j = 0; j < ns; j++) {
       cs0 += cs[j];
       cs1 += cs[ns+j];
      }
      cs0 /= (real)ns;
      cs1 /= (real)ns;
      /* (b) pitch angle calculation */
      /* The expression below should be refined as we need to know the
       orbital frequency of the nearest planet */
      omega = 1.0;
      dangle = (pow(Rinf[1],-1.5)-omega)/(.5*(cs0+cs1));
      dangle *= (Rmed[1]-Rmed[0]);
      i_angle = (int)(dangle/(PMAX-PMIN)*(real)NSEC+.5);
      /* (c) boundary on vrad and dens fields */
#pragma omp parallel for private(l,jp,lp,vr_med)
      for (j = 0; j < ns; j++) {
       l = j+i*ns;   // recall i=1
       jp = j+i_angle;
       if (jp >= ns) jp -= ns;
       if (jp < 0) jp += ns;
       lp = jp;
       rho[lp] = rho[l];       /* copy first ring into ghost ring */
       vr_med = -cs1*(rho[l]-SigmaMed[1])/SigmaMed[1];
       vr[l] = 2.*vr_med-vr[l+ns];
      }
      /* (d) density adjustment */
      mean = 0.0;
      for (j = 0; j < ns; j++) {
       mean += rho[j];
      }
      mean /= (real)ns;
      for (j = 0; j < ns; j++) {
       rho[j] += (SigmaMed[0]-mean);
      }
    }
    /* ------------------------ */
    /*      ADIABATIC CASE      */
    /* ------------------------ */
    if (EnergyEquation) {
      /* (a) sound speed in first two rings */
      cs0 = 0.0;
      cs1 = 0.0;
      for (j = 0; j < ns; j++) {
       cs0 += cs[j];
       cs1 += cs[ns+j];
      }
      cs0 /= (real)ns;
      cs1 /= (real)ns;
      csinit0 = sqrt(ADIABATICINDEX)*ASPECTRATIO*pow(Rmed[0],-0.5+FLARINGINDEX);
      csinit1 = sqrt(ADIABATICINDEX)*ASPECTRATIO*pow(Rmed[1],-0.5+FLARINGINDEX);
      dangle = (pow(Rinf[1],-1.5)-omega)/(.5*(csinit0+csinit1));
      dangle *= (Rmed[1]-Rmed[0]);
      i_angle = (int)(dangle/(PMAX-PMIN)*(real)NSEC+.5);
      /* (a) shift on energy and vrad */
#pragma omp parallel for private(l,jp,lp,vr_med)
      for (j = 0; j < ns; j++) {
       /* The expression below should be refined as we need to know the
          orbital frequency of the nearest planet */
       l = j+i*ns;   // recall i=1
       jp = j+i_angle;
       if (jp >= ns) jp -= ns;
       if (jp < 0) jp += ns;
       lp = jp;
       energy[lp] = energy[l];
       vr_med = -csinit1*(rho[l]-SigmaMed[1])/SigmaMed[1];
       vr[l] = 2.*vr_med-vr[l+ns];
      }
      /* (b) shift on acoustic part of density */
#pragma omp parallel for private(l,jp,lp,vr_med)
      for (j = 0; j < ns; j++) {
       /* The expression below should be refined as we need to know the
          orbital frequency of the nearest planet */
       rho[j] = SigmaMed[0] + (ADIABATICINDEX-1.0)*(energy[j]-EnergyMed[0])/csinit0/csinit0;
      }
      /* (c) density and energy adjustments */
      mean = 0.0;
      for (j = 0; j < ns; j++) {
       mean += rho[j];
      }
      mean /= (real)ns;
      for (j = 0; j < ns; j++) {
       rho[j] += (SigmaMed[0]-mean);
      }
      mean = 0.0;
      for (j = 0; j < ns; j++) {
       mean += energy[j];
      }
      mean /= (real)ns;
      for (j = 0; j < ns; j++) {
       energy[j] += (EnergyMed[0]-mean);
      }
    }
  }
  /* ======================================= */
  /* OUTER NON-REFLECTING BOUNDARY CONDITION */
  /* ======================================= */
  if (CPU_Rank == CPU_Highest) {
    csnrm2 = 0.0;
    csnrm1 = 0.0;
    for (j=0; j<ns; j++) {
      csnrm2 += cs[(nr-2)*ns+j];
      csnrm1 += cs[(nr-1)*ns+j];
    }
    csnrm1 /= (real)ns;
    csnrm2 /= (real)ns;
    i = nr-1;               /* The expression below should be refined */
    /* We need to know the orbital frequency of the nearest planet */
    dangle = (pow(Rinf[nr-2],-1.5)-1.0)/(.5*(csnrm1+csnrm2));
    dangle *= (Rmed[nr-1]-Rmed[nr-2]);
    i_angle = (int)(dangle/(PMAX-PMIN)*(real)NSEC+.5);
#pragma omp parallel for private(l,jp,lp,vr_med)
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      jp = j-i_angle;
      if (jp >= ns) jp -= ns;
      if (jp < 0) jp += ns;
      lp = jp+(i-1)*ns;
      rho[l] = rho[lp];              /* copy first ring into ghost ring */
      energy[l] = energy[lp];       /* copy first ring into ghost ring */
      if (!EnergyEquation) 
       vr_med = csnrm1*(rho[l-ns]-SigmaMed[nr-2])/SigmaMed[nr-2];
      else
       vr_med = cs[l]*(rho[l-ns]-SigmaMed[nr-2])/SigmaMed[nr-2];
      vr[l] = 2.*vr_med-vr[l-ns];
    }
    mean = 0.0;
    for (j = 0; j < ns; j++) {
      mean += rho[j+ns*(nr-1)];
    }
    mean /= (real)ns;
    for (j = 0; j < ns; j++) {
      rho[j+(nr-1)*ns] += SigmaMed[nr-1]-mean;
    }
    mean = 0.0;
    for (j = 0; j < ns; j++) {
      mean += energy[j+ns*(nr-1)];
    }
    mean /= (real)ns;
    for (j = 0; j < ns; j++) {
      energy[j+(nr-1)*ns] += EnergyMed[nr-1]-mean;
    }
  }
}


/* Sareh changed the inner boundary to have open in the planet shock heating
 * project */
void EvanescentBoundary (Vrad, Vtheta, Rho, Energy, step, mdot)
     PolarGrid *Vrad, *Vtheta, *Rho, *Energy;
     real step, mdot;
{
  int i, j, l, nr, ns;
  real *vrad, *vtheta, *dens, *energ;
  real vrad0, vtheta0, dens0, energ0;
  real damping, Tin, Tout, lambda;
  real temp0, temp1, visc;
  extern boolean DampToIni, DampToAxi, DecInner, OpInner;
  vrad = Vrad->Field;
  vtheta = Vtheta->Field;
  dens = Rho->Field;
  energ = Energy->Field;
  nr = Rho->Nrad;
  ns = Rho->Nsec;
  /* Orbital period at inner and outer boundary */
  Tin = 2.0*PI*pow(GlobalRmed[0],3./2);
  Tout = 2.0*PI*pow(GlobalRmed[GLOBALNRAD-1],3./2);
  /* WKZRMIN AND WKZRMAX are global Radii boundaries of killing wave zones */

  lambda = 0.0;
  for (i = 0; i < nr; i++) {
    if ( (Rmed[i] < WKZRMIN) || (Rmed[i] > WKZRMAX) ) {
      /* Damping operates only inside the wave killing zones */
      if (Rmed[i] < WKZRMIN) {
       damping = (Rmed[i]-WKZRMIN)/(GlobalRmed[0]-WKZRMIN);
       lambda = damping*damping*WKZTIN*step/Tin;
      }
      if (Rmed[i] > WKZRMAX) {
       damping = (Rmed[i]-WKZRMAX)/(GlobalRmed[GLOBALNRAD-1]-WKZRMAX);
       lambda = damping*damping*WKZTOUT*step/Tout;
      }
      // OLD Evanescent BC with damping wrt initial profiles
      if (!SelfGravity) {
       vtheta0 = sqrt ( G*1.0/Rmed[i] *                            \
                      ( 1.0 - (1.0+SIGMASLOPE-2.0*FLARINGINDEX)*       \
                        pow(AspectRatio(Rmed[i]),2.0)*pow(Rmed[i],2.0*FLARINGINDEX) ) );
      }
      if (SelfGravity) {
       vtheta0 = sqrt (  G*1.0/Rmed[i] *                            \
                       ( 1.0 - (1.0+SIGMASLOPE-2.0*FLARINGINDEX)*       \
                         pow(AspectRatio(Rmed[i]),2.0)*pow(Rmed[i],2.0*FLARINGINDEX) ) - \
                       Rmed[i]*GLOBAL_AxiSGAccr[i+IMIN] );
      }
      // this could be refined if CentrifugalBalance is used...
      /* NEW July 2012: two options -> damping wrt initial fields
       (default case) or damping wrt instantaneous axisymmetric
       fields */
      if (DampToIni) {
        if (Rmed[i] < WKZRMIN){
          if (OpInner){
          /* Standard outflow prescription */
            for (j = 0; j < ns; j++){
              l = i*ns+j;
              if (vrad[l+ns] > 0)
                vrad[l] = 0;
              else
                vrad[l] = vrad[l+ns];
              dens[l-ns] = dens[l];
              vrad[l-ns] = vrad[l];
              energ[l-ns] = energ[l];
            }
          } else {
           vtheta0 = VthetaMed[i];
           vrad0 = VradMed[i];
           dens0 = SigmaMed[i];
           energ0 = EnergyMed[i];
          }
        }
        if (Rmed[i] > WKZRMAX){
           vtheta0 = VthetaMed[i+IMIN];
           vrad0 = VradMed[i+IMIN];
           dens0 = SigmaMed[i];
           energ0 = EnergyMed[i];
        }
      }
      if (DampToAxi) {
        if (Rmed[i] < WKZRMIN){
            vtheta0 = 0.0;
            vrad0 = 0.0;
            dens0 = 0.0;
            energ0 = 0.0;
            for (j = 0; j < ns; j++) {
             l = i*ns + j;
             vrad0   += vrad[l];
             vtheta0 += vtheta[l];
             dens0   += dens[l];
             energ0  += energ[l];
           }
           vrad0   /= (real)ns;
           vtheta0 /= (real)ns;
           dens0   /= (real)ns;
           energ0  /= (real)ns;
       } else if (Rmed[i] > WKZRMAX) {
          vrad0   = 0.0;
         vtheta0 = 0.0;
         dens0   = 0.0;
         energ0  = 0.0;
         for (j = 0; j < ns; j++) {
          l = i*ns + j;
         vrad0   += vrad[l];
         vtheta0 += vtheta[l];
         dens0   += dens[l];
         energ0  += energ[l];
       }
       vrad0   /= (real)ns;
       vtheta0 /= (real)ns;
       dens0   /= (real)ns;
       energ0  /= (real)ns;
      }
      }
      /* Do not modify lines below */
      for (j = 0; j < ns; j++) {
       l = i*ns + j;
       vrad[l]   = (vrad[l]+lambda*vrad0)/(1.0+lambda);
       vtheta[l] = (vtheta[l]+lambda*vtheta0)/(1.0+lambda);
       dens[l]   = (dens[l]+lambda*dens0)/(1.0+lambda);
       if (EnergyEquation)
         energ[l]  = (energ[l]+lambda*energ0)/(1.0+lambda);
      }
    }
  }
}

/* The damping boundary for the radiative disc simulations as Alex does, 
 * only radial velocities are damped to zero*/
void AlexBoundary (Vrad, step)
     PolarGrid *Vrad;
     real step;
{
  int i, j, l, nr, ns;
  real *vrad, vrad0;
  real damping, Tau, lambda, A, B, C, dr;
  vrad = Vrad->Field;
  nr = Vrad->Nrad;
  ns = Vrad->Nsec;
  /* WKZRMIN AND WKZRMAX are global Radii boundaries of killing wave zones */
  lambda = 0.0;
  for (i = 0; i < nr; i++) {
    if ( (Rmed[i] < WKZRMIN) || (Rmed[i] > WKZRMAX) ) {
      /* Damping operates only inside the wave killing zones */
      if (Rmed[i] < WKZRMIN) {
       Tau = 2.0*PI*pow(GlobalRmed[0],3./2);
       dr = WKZRMIN-Rmed[0];
       A = 1./dr/dr;
       B = -2. * A * WKZRMIN;
       C = A * WKZRMIN * WKZRMIN;
      }
      if (Rmed[i] > WKZRMAX) {
       Tau = 2.0*PI*pow(GlobalRmed[GLOBALNRAD-1],3./2);
       dr = WKZRMAX-Rmed[GLOBALNRAD-1];
       A = 1./dr/dr;
       B = -2. * A * WKZRMAX;
       C = A * WKZRMAX * WKZRMAX;
      }
      damping = A*Rmed[i]*Rmed[i] + B*Rmed[i] + C;
      lambda = step/Tau*damping;
      vrad0 = 0.0;
      /* Do not modify lines below */
      for (j = 0; j < ns; j++) {
       l = i*ns + j;
       vrad[l]   -= (vrad[l]-vrad0)*lambda;
      }
    }
  }
}


void ApplyOuterSourceMass (Rho, Vrad)
     PolarGrid *Rho, *Vrad;
{
  int i, j, l, nr, ns;
  real *rho, average_rho = 0.0, *vr, penul_vr;
  nr = Rho->Nrad;
  ns = Rho->Nsec;
  rho= Rho->Field;
  vr = Vrad->Field;
  if (CPU_Rank == CPU_Highest) {
    i = nr-1;
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      average_rho += rho[l];
    }
    average_rho /= (real)ns;
    average_rho = SigmaMed[nr-1]-average_rho;
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      rho[l] += average_rho;
    }
    penul_vr = IMPOSEDDISKDRIFT*pow((Rinf[nr-1]/1.0),-SIGMASLOPE);
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      vr[l] = penul_vr;
    }
  }
}

void ApplySubKeplerianBoundary (Vtheta, sys)
     PolarGrid *Vtheta;
     PlanetarySystem *sys;
{
  int i, j, l, k, nr, ns;
  real VKepIn, VKepOut, totalmass, xp, yp, mp, rp;
  real *vt;
  vt = Vtheta->Field;
  nr = Vtheta->Nrad;
  ns = Vtheta->Nsec;
#pragma omp single
  {
    totalmass = 1.0;
    /* NEW (28 Sept 2011), but needs to be refined (the distance from
       the center of mass of all objects inside the inner edge should
       be calculated... */
    for (k = 0; k < sys->nb; k++) {
      xp = sys->x[k];
      yp = sys->y[k];
      mp = sys->mass[k];
      rp = sqrt( xp*xp + yp*yp );
      if (rp < GlobalRmed[0]) totalmass += mp;
    }
    if ( !SelfGravity ) {
      VKepIn = sqrt (  G*totalmass/Rmed[0] *                                   \
                     ( 1.0 - (1.0+SIGMASLOPE-2.0*FLARINGINDEX)/totalmass* \
                      pow(AspectRatio(Rmed[0]),2.0)*pow(Rmed[0],2.0*FLARINGINDEX) ) );
      VKepOut = sqrt (  G*totalmass/Rmed[nr-1] *                     \
                     ( 1.0 - (1.0+SIGMASLOPE-2.0*FLARINGINDEX)/totalmass* \
                       pow(AspectRatio(Rmed[nr-1]),2.0)*pow(Rmed[nr-1],2.0*FLARINGINDEX) ) );
    } else {
      if ( !SGZeroMode )
       mpi_make1Dprofile (SG_Accr, GLOBAL_AxiSGAccr);
      else
       GLOBAL_AxiSGAccr = SG_Accr;
      VKepIn = sqrt (  G*totalmass/Rmed[0] *                                   \
                     ( 1.0 - (1.0+SIGMASLOPE-2.0*FLARINGINDEX)/totalmass* \
                      pow(AspectRatio(Rmed[0]),2.0)*pow(Rmed[0],2.0*FLARINGINDEX) ) - \
                     Rmed[0]*GLOBAL_AxiSGAccr[0]/totalmass );
      VKepOut = sqrt (  G*totalmass/Rmed[nr-1] *                            \
                     ( 1.0 - (1.0+SIGMASLOPE-2.0*FLARINGINDEX)/totalmass* \
                       pow(AspectRatio(Rmed[nr-1]),2.0)*pow(Rmed[nr-1],2.0*FLARINGINDEX)/totalmass ) - \
                     Rmed[nr-1]*GLOBAL_AxiSGAccr[nr-1+IMIN] );
    }
    /* ----- */
    /* i = 0 */
    /* ----- */
    if ( CPU_Rank == 0 ) {
      i = 0;
      for (j = 0; j < ns; j++) {
       l = i*ns + j;
       vt[l] = VKepIn-Rmed[i]*OmegaFrame;
      }
    }
    /* ---------- */
    /* i = nr - 1 */
    /* ---------- */
    if ( CPU_Rank == CPU_Highest ) {
      i = nr - 1;
      for (j = 0; j < ns; j++) {
       l = i*ns + j;
       vt[l] = VKepOut-Rmed[i]*OmegaFrame;
      }
    }
  }
}

void ApplyBoundaryCondition (Vrad, Vtheta, Rho, Energy, step, sys)
     PolarGrid *Vrad, *Vtheta, *Rho, *Energy;
     real step;
     PlanetarySystem *sys;
{
  if (OpenInner == YES) OpenBoundary (Vrad, Vtheta, Rho, Energy);
  if (NonReflecting == YES) {
    if (EnergyEquation)
      ComputeSoundSpeed (Rho, Energy, sys);
    NonReflectingBoundary (Vrad, Rho, Energy, Vtheta);
  }
  if (Evanescent == YES) EvanescentBoundary (Vrad, Vtheta, Rho, Energy, step, 0);
  if (Alexboundary == YES) AlexBoundary (Vrad, step);
  /* New 'mixed' boundary condition, where an open BC is applied at
the grid's inner edge, and an evanescent BC at the outer edge (WKRMAX
needs to be specified) */
  if (MixedBC == YES) {
    OpenBoundary (Vrad, Vtheta, Rho, Energy);
    EvanescentBoundary (Vrad, Vtheta, Rho, Energy, step, 0);
  }
  if (AccBoundary == YES) AccretingBoundary(Vrad, Vtheta, Rho, Energy, step);
  if (OuterSourceMass == YES) ApplyOuterSourceMass (Rho, Vrad);
}

void CorrectVtheta (vtheta, domega)
     PolarGrid *vtheta;
     real domega;
{
  int i, j, l, nr, ns;
  real *vt;
  nr = vtheta->Nrad;
  ns = vtheta->Nsec;
  vt = vtheta->Field;
  for (i = 0; i < nr; i++) {
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      vt[l] -= domega*Rmed[i];
    }
  }
}

void DampDensity(Vrad, Vtheta, Rho, Energy, step, sys)
     PolarGrid *Vrad, *Vtheta, *Rho, *Energy;
     real step;
     PlanetarySystem *sys;
{
  int i, j, l, nr, ns;
  real *vrad, *vtheta, *dens, *energ;
  real Tadd, axidens;
  real lambda, x0, y0, rp0, x1, y1, rp1;
  real cutrmin, cutrmax, cutrmed, damping, normdampdist;
  real vrad0, vtheta0, energ0, dens0;
  vrad = Vrad->Field;
  vtheta = Vtheta->Field;
  dens = Rho->Field;
  energ = Energy->Field;
  nr = Rho->Nrad;
  ns = Rho->Nsec;
  /* OLD Damping used for adding mass in turbulent runs */
  /*
  for (i = 0; i < nr; i++) {
    // Damping timescale is 20 local orbital periods
    Tadd = 20.0*2.0*M_PI*pow(Rmed[i],1.5);
    axidens = 0.0;
    for (j = 0; j < ns; j++) {
      l = i*ns + j;
      axidens += dens[l];
    }
    // axisymmetric density at ring i 
    axidens /= (real)ns;
    for (j = 0; j < ns; j++) {
      l = i*ns + j;
      dens[l] = dens[l] - (axidens-SigmaMed[i])*dt/Tadd;
    }
  }
  */
  lambda = 0.0;
  x0 = sys->x[0];
  y0 = sys->y[0];
  rp0 = sqrt(x0*x0 + y0*y0);
  x1 = sys->x[1];
  y1 = sys->y[1];
  rp1 = sqrt(x1*x1 + y1*y1);
  cutrmin = rp0*(1.0+DENSDAMPRAD*ASPECTRATIO);
  cutrmax = rp1*(1.0-DENSDAMPRAD*ASPECTRATIO);
  cutrmed = 0.5*(cutrmin+cutrmax);
  for (i = 0; i < nr; i++) {
    if ( (Rmed[i] > cutrmin) && (Rmed[i] < cutrmax) ) {
      if ( (Rmed[i] > cutrmin) && (Rmed[i] < cutrmed) ) 
       normdampdist = (Rmed[i]-cutrmin)/(cutrmed-cutrmin);
      if ( (Rmed[i] >= cutrmed) && (Rmed[i] < cutrmax) )
       normdampdist = (Rmed[i]-cutrmax)/(cutrmed-cutrmax);
      damping = pow(sin((normdampdist)*0.5*M_PI),2.);
      lambda = damping*step;  // to be checked / customized...
      /* damping towards instantaneous axisymmetric fields... */
      vrad0   = 0.0;
      vtheta0 = 0.0;
      dens0   = 0.0;
      energ0  = 0.0;
      for (j = 0; j < ns; j++) {
       l = i*ns + j;
       vrad0   += vrad[l];
       vtheta0 += vtheta[l];
       dens0   += dens[l];
       energ0  += energ[l];
      }
      vrad0   /= (real)ns;
      vtheta0 /= (real)ns;
      dens0   /= (real)ns;
      energ0  /= (real)ns;
      /* Do not modify lines below */
      for (j = 0; j < ns; j++) {
       l = i*ns + j;
       vrad[l]   = (vrad[l]+lambda*vrad0)/(1.0+lambda);
       vtheta[l] = (vtheta[l]+lambda*vtheta0)/(1.0+lambda);
       dens[l]   = (dens[l]+lambda*dens0)/(1.0+lambda);
       if (EnergyEquation)
         energ[l]  = (energ[l]+lambda*energ0)/(1.0+lambda);
      }
    }
  }
}

void Evaporation(Rho, dt)
     PolarGrid *Rho;
     real dt;
{
  int i, j, l, nr, ns;
  real *dens;
  real axidens, floor;
  dens = Rho->Field;
  nr = Rho->Nrad;
  ns = Rho->Nsec;
  if (ievaporation == 0) {
    /* Calculate axisymmetric density field */
    GLOBAL_Axidens = (real *) malloc(sizeof(real) * GLOBALNRAD);
    if ( GLOBAL_Axidens==NULL ) {
      fprintf (stderr, "Not enough memory for allocation of GLOBAL_Axidens in SideEuler.c \n");
      prs_exit (1);
    }
    mpi_make1Dprofile (dens, axidens);
    ievaporation = 1;
  }

  /* Target density = minimum density allowed*/
  for (i = 0; i < nr; i++) {
    floor = 1e-3*GLOBAL_Axidens[i+IMIN];
    axidens = 0.0;
    for (j = 0; j < ns; j++) {
      l = i*ns + j;
      axidens += dens[l];
    }
    axidens /= (real)ns;
    for (j = 0; j < ns; j++) {
      l = i*ns + j;
      dens[l] = dens[l] - (axidens-floor)*dt/TEVAP/2.0/M_PI;
    }
  }
}

/* Calculate the scaling factor for photevaporation 
 * calculation using formula by Owen2012. */
real PhotoEvaporation(Vrad, Rho, dt)
     PolarGrid *Rho, *Vrad;
     real dt;
{
  int i, j, l, nr, ns;
  real *dens, *vrad;
  real x, sigmareduc;
  double scale;
  double mdot, summ;
  double world_summ;
  double checkmass;
  double world_mass, world_rhole;
  real l0, lx, l10, coeff;
  real a1,b1,c1,d1,e1,f1,g1;
  real sigd, axidens[GLOBALNRAD];
  real a2,b2,c2,d2,e2,f2;
  checkmass = 0.0;
  a2 = -0.438226;
  b2 = -0.10658387;
  c2 = 0.5699464;
  d2 = 0.010732277;
  e2 = -0.131809597;
  f2 = -1.32285709;
  l10 = log(10);
  a1 = 0.15138;
  b1 = -1.2182;
  c1 = 3.4046;
  d1 = -3.5717;
  e1 = -0.32762;
  f1 = 3.6064;
  g1 = -2.4918;
  dens = Rho->Field;
  vrad = Vrad->Field;
  nr = Rho->Nrad;
  ns = Rho->Nsec;
  /* In order to avoid getting stuck before an artificial bump close to 
   * the inner boundary, we make the critreria a little fluffy */
  if (sigcrit == floordens){
     if ((axidens[0] <= floordens)|| (axidens[1] <= floordens)) 
        sigcrit = 1.1 * floordens;
     else
        sigcrit = floordens;
  }
  while ((axidens[ihole] <= sigcrit) && (ihole <GLOBALNRAD)){
     Rhole = GlobalRmed[ihole];
     ihole++;
  }

  if (Rhole > GlobalRmed[1]) {
    for (i = 0; i <= nr; i++){
      if (Rsup[i] <= Rhole){
        for (j = 0; j <= ns; j++){
          l = i*ns+j;
          dens[l] = floordens;
        }
      }
    }
  summ=0;
  for (i = Zero_or_active; i < Max_or_active; i++){
     x = 0.95 *((Rmed[i]-Rhole)*unit_length/1.49598e11)/(unit_mass/1.9891e30);
     if (x>0){
        coeff = a2*b2*exp(b2*x) + c2*d2*exp(d2*x)+ e2*f2*exp(f2*x);
        coeff /= (Rmed[i]*unit_length/1.49598e11);
        sigd = coeff * exp(-pow((x/57),10));
     } else if (x<0){
        sigd = 0;
     }
     summ += ns*Surf[i]*sigd;
  }
  MPI_Allreduce (&summ, &world_summ, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  summ = world_summ;    
  summ *= (unit_length*1e2)*(unit_length*1e2);
  mdot = 4.8e-9 * pow((unit_mass/1.9891e30),(-0.148)) * pow((LX/1e30),1.14);
  mdot *= 1.9891e33/31556926;
  scale = mdot/summ; 
  scale /= (unit_mass*1e3); //scale is calculated in cgs
  scale *= (unit_length*1e2)*(unit_length*1e2);
  scale *= unit_time; //converting sigmadot to code units       
  for (i = 0; i < nr; i++) {
      x = 0.95 *((Rmed[i]-Rhole)*unit_length/1.49598e11)/(unit_mass/1.9891e30);
      if (x>=0){
         coeff = a2*b2*exp(b2*x) + c2*d2*exp(d2*x)+ e2*f2*exp(f2*x);
         coeff /= (Rmed[i]*unit_length/1.49598e11);
         sigd = coeff * exp(-pow((x/57),10));
      } else if (x<0){
         sigd = 0;
      }
      sigmareduc = scale * sigd ;
      for (j = 0; j < ns; j++) {
         l = i*ns + j;
         dens[l] = dens[l] - sigmareduc*dt;
         if ((dens[l] <= floordens)){
             dens[l] = floordens;
         } else {
             checkmass += sigmareduc * Surf[i] * dt;
         }
      }
    }
  } else {
    summ = 0;
    for (i = Zero_or_active; i < Max_or_active; i++){
        x = 0.85 *(Rmed[i]*unit_length/1.49598e11)/(unit_mass/1.9891e30);
        l0 = log10(x);
        lx = log(x);
        if ( x>0.7 ){
           sigd = pow(10,(a1*pow(l0,6) + b1 * pow(l0,5)+ c1 * pow(l0,4)+ d1* pow(l0,3)\
                  + e1* pow(l0,2) + f1* l0+ g1));
           coeff = 6*a1*pow(lx,5)/pow(x,2)/pow(l10,7) + 5*b1*pow(lx,4)/pow(x,2)/pow(l10,6)\
                  + 4*c1*pow(lx,3)/pow(x,2)/pow(l10,5);
           coeff += 3*d1*pow(lx,2)/pow(x,2)/pow(l10,4) + 2*e1*lx/pow(x,2)/pow(l10,3) + f1/pow(x,2)/pow(l10,2);
           sigd *= coeff * exp(-pow((x/100.),10));
        } else if (x<0.7){
           sigd = 0;
        }
        summ += ns*Surf[i]*sigd;
    }
    MPI_Allreduce (&summ, &world_summ, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    summ = world_summ;
    summ *= (unit_length*1e2)*(unit_length*1e2);
    mdot = 6.25e-9 * pow((unit_mass/1.9891e30),-0.068) * pow((LX/1e30),1.14);
    mdot *= 1.9891e33/31556926;
    scale = mdot/summ;
    scale /= (unit_mass*1e3); //scale is calculated in sgs
    scale *= (unit_length*1e2)*(unit_length*1e2);
    scale *= unit_time; //converting sigmadot to code units   
    for (i = 0; i < nr; i++) {
        x = 0.85 *(Rmed[i]*unit_length/1.49598e11)/(unit_mass/1.9891e30);
        l0 = log10(x);
        lx = log(x);
        if (x>0.7){
           sigd = pow(10,(a1*pow(l0,6) + b1 * pow(l0,5)+ c1 * pow(l0,4)+ d1* pow(l0,3)\
                 + e1* pow(l0,2) + f1* l0+ g1));
           coeff = 6*a1*pow(lx,5)/pow(x,2)/pow(l10,7) + 5*b1*pow(lx,4)/pow(x,2)/pow(l10,6)\
                 + 4*c1*pow(lx,3)/pow(x,2)/pow(l10,5);
           coeff += 3*d1*pow(lx,2)/pow(x,2)/pow(l10,4) + 2*e1*lx/pow(x,2)/pow(l10,3) + f1/pow(x,2)/pow(l10,2);
           sigd *= coeff * exp(-pow((x/100.),10));
        } else if (x<0.7){
           sigd = 0;
        }       
        sigmareduc = scale * sigd ;
        SigmaDotW[i] = sigmareduc;
        SigmaDotV[i] = 1.5*ALPHAVISCOSITY*(ASPECTRATIO*ASPECTRATIO)*dens[i]*pow(Rmed[i],-1.5);
        for (j = 0; j < ns; j++) {
           l = i*ns + j;
           dens[l] = dens[l] - sigmareduc*dt;
           if (dens[l] <= floordens){
                dens[l] = floordens;
           } else {
                checkmass += sigmareduc * Surf[i] * dt;
           }
        }
    }
  }
  MPI_Allreduce (&checkmass, &world_mass, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  checkmass = world_mass;
  return checkmass;
}

void WriteSigmaDotFile(timestep)
       int timestep;
{
  FILE *output;
  char name[256];
  int i;
  if (!CPU_Master) return;
  printf ("Updating 'SigmaDotFile.dat'...");
  fflush (stdout);
  sprintf (name, "%sSigmaDot%d.dat", OUTPUTDIR, timestep);
  output = fopen (name, "a");
  if (output == NULL) {
    fprintf (stderr, "Can't write 'SigmaDot%d.dat' file. Aborting.\n", timestep);
    prs_exit (1);
  }
  for (i = 0; i < GLOBALNRAD; i++){
       fprintf (output, "%#.18g\t%#.18g\t%#.18g\n", Rmed[i], SigmaDotW[i], SigmaDotV[i]);
  }
  fclose (output);
  printf ("done\n");
  fflush (stdout);
}

void SetRhoFloor(Rho)
  PolarGrid *Rho;
{
  int nr, ns, i, j, l;
  real Sfloor;
  real *dens;
  nr = Rho->Nrad;
  ns =  Rho->Nsec;
  dens = Rho->Field;
  
  Sfloor = 1e-7 * SIGMA0;
  for (i = 0; i < nr; i++){
    for (j = 0; j < ns; j++){
       l = i*ns + j;
         if (dens[l] < Sfloor)
              dens[l] = Sfloor;
    }
  }
}
    
void SetEnergyFloor(Rho, Energy)
  PolarGrid *Energy;
  PolarGrid *Rho;
{
  int nr, ns, i, j, l;
  real Efloor, Tfloor;
  real *energy, *dens, *temp;
  nr = Energy->Nrad;
  ns =  Energy->Nsec;
  energy = Energy->Field;
  dens = Rho->Field;
  temp = Temperature->Field;
  Tfloor = 3.0/unit_temperature;
  for (i = 0; i < nr; i++){
    for (j = 0; j < ns; j++){
       l = i*ns + j;
       Efloor = Tfloor * dens[l]/(ADIABATICINDEX-1);
       if (energy[l] < Efloor){
           energy[l] = Efloor;
           temp[l] = energy[l] / dens[l] * (ADIABATICINDEX-1);
       }
    }
  }
}

real BitschTemperature(Mdot,r)
  real Mdot;
  real r;
{
  real Btemp, chi, rau;
  Mdot *= (unit_mass/unit_time) / (1.9891e30/31556926.0); //convert to Msun/year
  Mdot *= -1;  //Making it positive
  rau = r * FACTORUNITLENGTH;
  chi = log10(ZMETAL/0.5)/log10(2.);
  if (ZMETAL <0.005)
    Btemp = LowMetal(Mdot,rau, chi);
  else
    Btemp = HighMetal(Mdot,rau,chi);
  return Btemp/unit_temperature; 
}

real LowMetal(Mdot,r,chi)
   real Mdot;
   real r;
   real chi;
{
  real lowtemp, Mdotref, muR;
  real A1, B1, si12;
  real A2, B2, si21, si23;
  real A3, B3, si32;
  real R12, R21, R23, R32;
  if (Mdot > 3.5e-8){
     Mdotref = 3.5e-8;
     muR = log10(Mdot/Mdotref);
     A1 = 410.*pow(1.4,chi);
     B1 = 1.3 * pow(0.7,chi);
     si12 = -(atan(r-3.*(pow(1.832,muR)+(0.15+0.05*muR*muR)*chi))-PI/2.)/PI;
     A2 = 650 * pow(1.3,chi);
     B2 = 1.8 * pow(0.75,chi);
     si21 = (atan(r-2.6*(pow(1.832,muR)+(0.1666+0.05*muR*muR)*chi))+PI/2.)/PI;
     si23 = (-atan(r-8.3*(pow(1.372,muR)+(0.1+muR*muR)*chi))+PI/2.)/PI;
     A3 = 185 * pow(0.995,chi);
     B3 = 1.4 * pow(1.1,chi);
     si32 = (atan(r-8.3*(pow(1.372,muR)+(0.1+muR*muR)*chi))+PI/2.)/PI;
     lowtemp = A1* pow(B1,muR)* pow(r,-6./7-log10(ZMETAL+0.5)/4.)* si12 \
             + A2* pow(B2,muR)* pow(r,-8./7)* si21*si23 \
             + A3* pow(B3,muR)* pow(r,-4./7)* si32;
  } else {
     if (Mdot > 1.75e-8){
       Mdotref = 3.5e-8;
       muR = log10(Mdot/Mdotref);
       A1 = 410.*pow(1.4,chi);
       B1 = 3. * pow(0.8,chi);
       si12 = -(atan(r-3.*(pow(1.832,muR)+(0.15+0.2*muR*muR)*chi))-PI/2.)/PI;
       A2 = 650 * pow(1.3,chi);
       B2 = 2.1 * pow(1.1,chi);
       si21 = (atan(r-2.6*(pow(1.832,muR)+(0.1666+0.2*muR*muR)*chi))+PI/2.)/PI;
       si23 = (-atan(r-8.3*(pow(1.372,muR)+(0.15-2.5*muR*muR)*chi))+PI/2.)/PI;
       A3 = 185 * pow(0.995,chi);
       B3 = 1.176 * pow(1.01,chi);
       si32 = (atan(r-8.3*(pow(1.372,muR)+(0.1-2.5*muR*muR)*chi))+PI/2.)/PI;
       lowtemp = A1* pow(B1,muR)* pow(r,-6./7-log10(ZMETAL+0.5)/4.)* si12 \
               + A2* pow(B2,muR)* pow(r,-8./7)* si21*si23 \
               + A3* pow(B3,muR)* pow(r,-4./7)* si32;
     } else {
        if (Mdot > 8.75e-9){
         Mdotref = 1.75e-8;
         muR = log10(Mdot/Mdotref);
         R12 = pow(1.832,log10(0.5))+(0.15+0.2*(log10(0.5)*log10(0.5)))*chi;
         R21 = pow(1.832,log10(0.5))+(0.1666+0.2*(log10(0.5)*log10(0.5)))*chi;
         R23 = pow(1.372,log10(0.5))+(0.1-2.5*(log10(0.5)*log10(0.5)))*chi;
         R32 = pow(1.372,log10(0.5))+(0.1-2.5*(log10(0.5)*log10(0.5)))*chi;
         A1 = 410.* pow(3*pow(0.8,chi),log10(0.5))* pow(1.4,chi);
         B1 = 8.* pow(1.5,chi);
         si12 = -(atan(r-3.*R12*pow(3.-0.2*chi*chi,muR))-PI/2.)/PI;
         A2 = 650* pow(2.1*pow(1.1,chi),log10(0.5))* pow(1.3,chi);
         B2 = 30.0 * pow(1.1,chi);
         si21 = (atan(r-2.6*R21*pow(5-0.2*chi*chi,muR))+PI/2.)/PI;
         si23 = (-atan(r-8.3*R23*pow(3-4.5*chi,muR))+PI/2.)/PI;
         A3 = 185* pow(1.176*pow(1.01,chi),log10(0.5))* pow(0.995,chi);
         B3 = 2.2 * pow(0.125,chi);
         si32 = (atan(r-8.3* R32* pow(3.5-4.5*chi,muR))+PI/2.)/PI;
         lowtemp = A1* pow(B1,muR)* pow(r,-6./7-log10(ZMETAL+0.5)/4.*(1-muR/log10(0.5)))* si12 \
                 + A2* pow(B2,muR)* pow(r,-8./7)* si21*si23 \
                 + A3* pow(B3,muR)* pow(r,-4./7-muR*(1./ZMETAL)/8.)* si32;
        } else {
           Mdotref = 8.75e-9;
           muR = log10(Mdot/Mdotref);
           R12 = pow(1.832,log10(0.5))+(0.15+0.2*(log10(0.5)*log10(0.5)))*chi \
                 *pow(3-0.2*chi*chi,log10(0.5));
           R21 = pow(1.832,log10(0.5))+(0.1666+0.2*(log10(0.5)*log10(0.5)))*chi \
                 *pow(5-0.2*chi*chi,log10(0.5));
           R23 = pow(1.372,log10(0.5))+(0.1-2.5*(log10(0.5)*log10(0.5)))*chi \
                 *pow(3-4.5*chi,log10(0.5));
           R32 = pow(1.372,log10(0.5))+(0.1-2.5*(log10(0.5)*log10(0.5)))*chi \
                 *pow(3.5-4.5*chi,log10(0.5));
           A1 = 410.* pow(3*pow(0.8,chi),log10(0.5))* pow(8*pow(1.5,chi),log10(0.5))* pow(1.4,chi);
           B1 = 2.* pow(1.4,chi);
           si12 = -(atan(r-3.*R12*pow(3.,muR))-PI/2.)/PI;
           A2 = 650* pow(2.1*pow(1.1,chi),log10(0.5))* pow(30.*pow(1.1,chi),log10(0.5))* pow(1.3,chi);
           B2 = 1.5 * pow(0.7,chi);
           si21 = (atan(r-2.6*R21*pow(5,muR))+PI/2.)/PI;
           si23 = (-atan(r-8.3*R23*pow(-.5*chi-2*muR,muR))+PI/2.)/PI;
           A3 = 185* pow(1.176*pow(1.01,chi),log10(0.5))* pow(2.2*pow(0.125,chi),log10(0.5))* pow(0.995,chi);
           B3 = 1.9 * pow(0.9,chi);
           si32 = (atan(r-8.3* R32* pow(-.5*chi-2.*muR,muR))+PI/2.)/PI;
           lowtemp = A1* pow(B1,muR)* pow(r,-6./7)* si12 \
                   + A2* pow(B2,muR)* pow(r,-8./7)* si21*si23 \
                   + A3* pow(B3,muR)* pow(r,-4./7-log10(0.5)*(1./ZMETAL)/8)* si32;
        }
     }
  }
  return lowtemp;
}

real HighMetal(Mdot,r,chi)
   real Mdot;
   real r;
   real chi;
{
  real hightemp, Mdotref, muR;
  real A1, B1, si12;
  real A2, B2, si21, si23;
  real A3, B3, si32;
  real R12, R21, R23, R32;
  if (Mdot > 3.5e-8){
     Mdotref = 3.5e-8;
     muR = log10(Mdot/Mdotref);
     A1 = 410.*pow(1.175,chi);
     B1 = 1.3 * pow(1.2,chi);
     si12 = -(atan(r-3.*(pow(1.832,muR)+(0.15+0.55*muR*muR)*chi))-PI/2.)/PI;
     A2 = 650 * pow(1.15,chi);
     B2 = 1.8 * pow(1.15,chi);
     si21 = (atan(r-2.6*(pow(1.832,muR)+(0.1666+0.666*muR*muR)*chi))+PI/2.)/PI;
     si23 = (-atan(r-8.3*(pow(1.372,muR)+(0.2+muR*muR)*chi))+PI/2.)/PI;
     A3 = 185 * pow(1.025,chi);
     B3 = 1.4 * pow(1.05,chi);
     si32 = (atan(r-8.3*(pow(1.372,muR)+(0.2+muR*muR)*chi))+PI/2.)/PI;
     hightemp = A1* pow(B1,muR)* pow(r,-6./7-log10(ZMETAL+0.5)/4.)* si12 \
             + A2* pow(B2,muR)* pow(r,-8./7)* si21*si23 \
             + A3* pow(B3,muR)* pow(r,-4./7)* si32;
  } else {
     if (Mdot > 1.75e-8){
       Mdotref = 3.5e-8;
       muR = log10(Mdot/Mdotref);
       A1 = 410.*pow(1.175,chi);
       B1 = 3. * pow(0.9,chi);
       si12 = -(atan(r-3.*(pow(1.832,muR)+(0.15-0.4*muR*muR)*chi))-PI/2.)/PI;
       A2 = 650 * pow(1.15,chi);
       B2 = 2.1 * pow(1.15,chi);
       si21 = (atan(r-2.6*(pow(1.832,muR)+(0.1666-0.5*muR*muR)*chi))+PI/2.)/PI;
       si23 = (-atan(r-8.3*(pow(1.372,muR)+(0.2-1.5*muR*muR)*chi))+PI/2.)/PI;
       A3 = 185 * pow(1.025,chi);
       B3 = 1.176 * pow(1.05,chi);
       si32 = (atan(r-8.3*(pow(1.372,muR)+(0.2-1.5*muR*muR)*chi))+PI/2.)/PI;
       hightemp = A1* pow(B1,muR)* pow(r,-6./7-log10(ZMETAL+0.5)/4.)* si12 \
               + A2* pow(B2,muR)* pow(r,-8./7)* si21*si23 \
               + A3* pow(B3,muR)* pow(r,-4./7)* si32;
     } else {
        if (Mdot > 8.75e-9){
         Mdotref = 1.75e-8;
         muR = log10(Mdot/Mdotref);
         R12 = pow(1.832,log10(0.5))+(0.15-0.4*(log10(0.5)*log10(0.5)))*chi;
         R21 = pow(1.832,log10(0.5))+(0.1666-0.5*(log10(0.5)*log10(0.5)))*chi;
         R23 = pow(1.372,log10(0.5))+(0.2-1.5*(log10(0.5)*log10(0.5)))*chi;
         R32 = pow(1.372,log10(0.5))+(0.2-1.5*(log10(0.5)*log10(0.5)))*chi;
         A1 = 410.* pow(3*pow(0.9,chi),log10(0.5))* pow(1.175,chi);
         B1 = 8.* pow(1.2,chi);
         si12 = -(atan(r-3.*R12*pow(3.-0.2*chi*chi,muR))-PI/2.)/PI;
         A2 = 650* pow(2.1*pow(1.15,chi),log10(0.5))* pow(1.15,chi);
         B2 = 30.0 * pow(1.2-chi/(chi+3),chi);
         si21 = (atan(r-2.6*R21*pow(5-0.2*chi*chi,muR))+PI/2.)/PI;
         si23 = (-atan(r-8.3*R23*pow(3+chi,muR))+PI/2.)/PI;
         A3 = 185* pow(1.176*pow(1.05,chi),log10(0.5))* pow(1.025,chi);
         B3 = 2.2 * pow(0.9,chi);
         si32 = (atan(r-8.3* R32* pow(3.5+chi,muR))+PI/2.)/PI;
         hightemp = A1* pow(B1,muR)* pow(r,-6./7-log10(ZMETAL+0.5)/4.*(1-muR/log10(0.5)))* si12 \
                 + A2* pow(B2,muR)* pow(r,-8./7)* si21*si23 \
                 + A3* pow(B3,muR)* pow(r,-4./7-muR*(1./ZMETAL)/8.)* si32;
        } else {
           Mdotref = 8.75e-9;
           muR = log10(Mdot/Mdotref);
           R12 = pow(1.832,log10(0.5))+(0.15-0.4*(log10(0.5)*log10(0.5)))*chi \
                 *pow(3-0.2*chi*chi,log10(0.5));
           R21 = pow(1.832,log10(0.5))+(0.1666-0.5*(log10(0.5)*log10(0.5)))*chi \
                 *pow(5-0.2*chi*chi,log10(0.5));
           R23 = pow(1.372,log10(0.5))+(0.2-1.5*(log10(0.5)*log10(0.5)))*chi \
                 *pow(3+chi,log10(0.5));
           R32 = pow(1.372,log10(0.5))+(0.2-1.5*(log10(0.5)*log10(0.5)))*chi \
                 *pow(3.5+chi,log10(0.5));
           A1 = 410.* pow(3*pow(0.9,chi),log10(0.5))* pow(8*pow(1.2,chi),log10(0.5))* pow(1.175,chi);
           B1 = 2.* pow(1.4,chi);
           si12 = -(atan(r-3.*R12*pow(3.,muR))-PI/2.)/PI;
           A2 = 650* pow(2.1*pow(1.15,chi),log10(0.5))* pow(30.*pow(1.2-chi/(chi+3),chi),log10(0.5))* pow(1.15,chi);
           B2 = 1.5 * pow(1.8,chi);
           si21 = (atan(r-2.6*R21*pow(5,muR))+PI/2.)/PI;
           si23 = (-atan(r-8.3*R23*pow(0.1*chi-2*muR,muR))+PI/2.)/PI;
           A3 = 185* pow(1.176*pow(1.05,chi),log10(0.5))* pow(2.2*pow(0.9,chi),log10(0.5))* pow(1.025,chi);
           B3 = 1.9 * pow(0.9,chi);
           si32 = (atan(r-8.3* R32* pow(0.2*chi*chi-2.*muR,muR))+PI/2.)/PI;
           hightemp = A1* pow(B1,muR)* pow(r,-6./7)* si12 \
                   + A2* pow(B2,muR)* pow(r,-8./7)* si21*si23 \
                   + A3* pow(B3,muR)* pow(r,-4./7-log10(0.5)*(1./ZMETAL)/8)* si32;
        }
     }
  }
  return hightemp;
}

