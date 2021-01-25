/** \file Theo.c

A few functions that manipulate the surface density, internal energy
and cooling time profiles.

*/

#include "mp.h"

extern real ScalingFactor;
extern boolean AccBoundary;

/* Surface density */
real Sigma(r)
     real r;
{
  real cavity = 1.0;
  real sigma, sigref, rdecay, rmin, rmax;
  real sigmabg, deltz, rtarget, deltmod, fedge;
  extern boolean ExponentialDecay;
  //if (r < CAVITYRADIUS) cavity = 1.0/CAVITYRATIO; 
  rmin = CAVITYRADIUS-CAVITYWIDTH*ASPECTRATIO;
  rmax = CAVITYRADIUS+CAVITYWIDTH*ASPECTRATIO;
  if (r < rmin) cavity /= CAVITYRATIO;
  if ((r >= rmin) && (r <= rmax)) {
    cavity /= exp((rmax-r)/(rmax-rmin)*log(CAVITYRATIO));
  }
  /* This is *not* a steady state */
  /* profile, if a cavity is defined. It first needs */
  /* to relax towards steady state, on a viscous time scale */
  sigma = cavity*ScalingFactor*SIGMA0*pow(r,-SIGMASLOPE);
  /* Smooth transition as the inner edge of a disc */
  if (EdgeTransition){
    fedge = 1./(exp(-(r-EDGERMID)/EDGEDELTA)+1);
    sigma *= fedge + EDGESIGMADROP*(1-fedge);
  }
  if (ExponentialDecay) {
     rdecay = 0.5*RMAX;
     sigma *= exp(-r/rdecay);
  }
  return sigma;
}

void FillSigma()
{
  int i, ind1;
  real visc, dsmoothin=0.1, mdot0;
  extern boolean MdotHartmann, DecInner, EnergyEquation;
  if (AccBoundary){
    if (ViscosityAlpha){
      if (EnergyEquation) 
      visc = ADIABATICINDEX *ALPHAVISCOSITY * (ASPECTRATIO*ASPECTRATIO);
      else
      visc = ALPHAVISCOSITY * (ASPECTRATIO*ASPECTRATIO);
    } else { 
      visc = VISCOSITY;
    }
    SIGMA0 = MDOTINIT / 3./PI/visc;
    if (MdotHartmann){
      mdot0 = 1e-8 * pow(THARTMANN/1e6, -1.4);
      mdot0 *= (1.9891e30/31556926.0 / unit_mass*unit_time);
      SIGMA0 = mdot0/ 3./PI/visc;
    }
  }
  for (i = 0; i < NRAD; i++) {
    SigmaMed[i] = Sigma(Rmed[i]);
    SigmaInf[i] = Sigma(Rinf[i]);
    if ((AccBoundary) && (DecInner)) {
      if (Rmed[i] < (GlobalRmed[0]+dsmoothin))
        SigmaMed[i] = (SigmaMed[i]-floordens) * exp(-pow(Rmed[i]-(GlobalRmed[0]+dsmoothin),2)/2./(dsmoothin*dsmoothin)) + floordens;
      if (Rinf[i] < (Radii[0]+dsmoothin))
        SigmaInf[i] = (SigmaInf[i]-floordens) * exp(-pow(Radii[i]-(Radii[0]+dsmoothin),2)/2./(dsmoothin*dsmoothin)) + floordens;
      SigmaMed[0] = SigmaMed[1];
    }
  }
}

void RefillSigma (Surfdens)
     PolarGrid *Surfdens;
{
  int i, j, nr, ns, l;
  real *field;
  real moy;
  nr = Surfdens->Nrad;
  ns = Surfdens->Nsec;
  field = Surfdens->Field;
  for (i = 0; i < nr; i++) {
    moy = 0.0;
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      moy += field[l];
    }
    moy /= (real)ns;
    SigmaMed[i] = moy;
  }
  SigmaInf[0] = SigmaMed[0];
  for (i = 1; i < nr; i++) {
    SigmaInf[i] = (SigmaMed[i-1]*(Rmed[i]-Rinf[i])+\
    SigmaMed[i]*(Rinf[i]-Rmed[i-1]))/(Rmed[i]-Rmed[i-1]);
  }
}

/* Thermal energy */
real Energy(r)
     real r;
{
  int i;
  real energy0;
  real cavity = 1.0;
  extern boolean DecInner;
  if (ADIABATICINDEX == 1.0) {
    fprintf (stderr, "The adiabatic index must differ from unity to initialize the gas internal energy. I must exit.\n");
    prs_exit (1);
  } else {
    if ((AccBoundary) && (DecInner)){
      i = ReturnIndex(r);
      energy0 = R/MU/(ADIABATICINDEX-1.0)*SigmaMed[i-1-IMIN]*(ASPECTRATIO*ASPECTRATIO)*pow(r,-1.0+2.0*FLARINGINDEX);
    } else { 
      energy0 = R/MU/(ADIABATICINDEX-1.0)*SIGMA0*(ASPECTRATIO*ASPECTRATIO)*pow(r,-SIGMASLOPE-1.0+2.0*FLARINGINDEX);
    }
  }
  if (r < CAVITYRADIUS) cavity = 1.0/CAVITYRATIO; 
  return cavity*ScalingFactor*energy0;
}

void FillEnergy() {
  int i;
  for (i = 0; i < NRAD; i++)
    EnergyMed[i] = Energy(Rmed[i]);
}


void RefillEnergy (energy)
     PolarGrid *energy;
{
  int i, j, nr, ns, l;
  real *field;
  real moy;
  nr = energy->Nrad;
  ns = energy->Nsec;
  field = energy->Field;
  for (i = 0; i < nr; i++) {
    moy = 0.0;
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      moy += field[l];
    }
    moy /= (real)ns;
    EnergyMed[i] = moy;
  }
}

/* Temperature prescription time */
real PrescTime(r)
     real r;
{
  real pt0;
  pt0 = PRESCTIME0*pow(r,2.0+2.0*FLARINGINDEX);
  return pt0;
}

void FillPrescTime() {
  int i;
  for (i = 0; i < NRAD; i++)
    PrescTimeMed[i] = PrescTime(Rmed[i]);
}


/* Beta cooling prescription time */
real ComputeBetaCooling(r)
     real r;
{
  real pt0;
  pt0 = BETACOOLINGTIME*pow(r,-BETACOOLINGSLOPE);
  return pt0;
}
