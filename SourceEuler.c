/** \file SourceEuler.c 

Contains routines used by the hydrodynamical loop. More specifically,
it contains the main loop itself and all the source term substeps
(with the exception of the evaluation of the viscous force). The
transport substep is treated elsewhere. */

#include "mp.h"

#define CFLSECURITY 0.5              /* Maximum fraction of zone size */
                            /* swept in one timestep */

#define CVNR 1.41              /* Shocks are spread over CVNR zones:       */
                                /* von Neumann-Richtmyer viscosity constant */
                            /* Beware of misprint in Stone and Norman's */
                            /* paper : use C2^2 instead of C2           */

#define BINARYSECURITY 100.0

static PolarGrid *TemperInt;
static PolarGrid *VradNew, *VradInt;
static PolarGrid *VthetaNew, *VthetaInt;
static PolarGrid *EnergyNew, *EnergyInt, *TempInt;
static real timeCRASH;  
extern boolean Corotating, ModifiedSoundSpeed;
extern boolean EnergyEquation, ThermalDiffusion, ThermalCooling, ViscousHeating, RadiativeDiffusion;
extern boolean SelfGravity, ZMPlus;
int FirstGasStepFLAG=1;
static int AlreadyCrashed = 0, GasTimeStepsCFL;
extern boolean FastTransport, IsDisk;
extern boolean AddMass, DontApplySubKeplerian, DiscEvaporation;
Pair DiskOnPrimaryAcceleration;
extern boolean PhotoEvapor, StellarIrradiation;

boolean DetectCrash (array)
     PolarGrid *array;
{
  int i, j, l, nr, ns;
  real *ptr;
  boolean bool = NO;
  nr = array->Nrad;
  ns = array->Nsec;
  ptr= array->Field;
#pragma omp parallel for private(j,l) shared(bool)
  for (i = 0; i < nr; i++) {
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      if (ptr[l] < 0.0) 
       bool = YES;
    }
  }
  return bool;
}
 
void FillPolar1DArrays ()
{
  extern boolean Restart;
  FILE *input, *output;
  int i, ii, j;
  real drrsep;
  float temporary;
  char InputName[256], OutputName[256];
  drrsep = (RMAX-RMIN)/(real)GLOBALNRAD;
  sprintf (InputName, "%s%s", OUTPUTDIR, "radii.dat");
  if (!Restart) 
    sprintf (OutputName, "%s%s", OUTPUTDIR, "used_rad.dat");
  else
    sprintf (OutputName, "%s%s", OUTPUTDIR, "newused_rad.dat");
  input = fopen (InputName, "r");
  if (input == NULL) {
    mastererr ("Warning : no `radii.dat' file found. Using default.\n");
    if (LogGrid == YES) {
      for (i = 0; i <= GLOBALNRAD; i++)
       Radii[i] = RMIN*exp((real)i/(real)GLOBALNRAD*log(RMAX/RMIN));
    } else {
      for (i = 0; i <= GLOBALNRAD; i++)
       Radii[i] = RMIN+drrsep*(real)(i);
    }
  } else {
    mastererr ("Reading 'radii.dat' file.\n");
    for (i = 0; i <= GLOBALNRAD; i++) {
      fscanf (input, "%f", &temporary);
      Radii[i] = (real)temporary;
    }
  }
  for (i = 0; i < GLOBALNRAD; i++) {
    GlobalRmed[i] = 2.0/3.0*(Radii[i+1]*Radii[i+1]*Radii[i+1]-Radii[i]*Radii[i]*Radii[i]);
    GlobalRmed[i] = GlobalRmed[i] / (Radii[i+1]*Radii[i+1]-Radii[i]*Radii[i]);
  }
  for (i = 0; i < NRAD; i++) {
    ii = i+IMIN;
    Rinf[i] = Radii[ii];
    Rsup[i] = Radii[ii+1];
    Rmed[i] = 2.0/3.0*(Rsup[i]*Rsup[i]*Rsup[i]-Rinf[i]*Rinf[i]*Rinf[i]);
    Rmed[i] = Rmed[i] / (Rsup[i]*Rsup[i]-Rinf[i]*Rinf[i]);
    Surf[i] = 0.5*(PMAX-PMIN)*(Rsup[i]*Rsup[i]-Rinf[i]*Rinf[i])/(real)NSEC;
    InvRmed[i] = 1.0/Rmed[i];
    InvSurf[i] = 1.0/Surf[i];
    InvDiffRsup[i] = 1.0/(Rsup[i]-Rinf[i]);
    InvRinf[i] = 1.0/Rinf[i];
  }
  Rinf[NRAD]=Radii[NRAD+IMIN];
  for (i = 1; i < NRAD; i++) {
    InvDiffRmed[i] = 1.0/(Rmed[i]-Rmed[i-1]);
  }
  if (CPU_Master) {
    output = fopen (OutputName, "w");
    if (output == NULL) {
      mastererr ("Can't write %s.\nProgram stopped.\n", OutputName);
      prs_exit (1);
    }
    for (i = 0; i <= GLOBALNRAD; i++)
      fprintf (output, "%.18g\n", Radii[i]);
    fclose (output);
  }
  if (input != NULL) fclose (input);
  /* ---------------------------------- */
  /* New (Apr 12 2010): definition of a global array with cells azimuths */
  /* ---------------------------------- */
  for (j = 0; j < NSEC; j++) {
    azimuth[j] = PMIN + (PMAX-PMIN)*(real)j/(real)NSEC;
    /* case where azimuthal extent smaller than 2pi */
    if ( fabs(PMAX-PMIN-2.*M_PI) > 0.2 ) 
      azimuth[j] = PMIN + (PMAX-PMIN)*(real)j/(real)(NSEC-1);
  }
  if (CPU_Master) {
    sprintf (OutputName, "%s%s", OUTPUTDIR, "used_azi.dat");
    output = fopen (OutputName, "w");
    if (output == NULL) {
      mastererr ("Can't write %s.\nProgram stopped.\n", OutputName);
      prs_exit (1);
    }
    for (j = 0; j < NSEC; j++) {
      fprintf (output, "%.18g\n", azimuth[j]);
    }
    fclose (output);
  }
}

void InitEuler (Vr, Vt, Rho, Energy, sys)
     PolarGrid *Vr, *Vt, *Rho, *Energy;
     PlanetarySystem *sys;
{
  extern boolean ImposedAlpha, InitEquilibrium;;
  InitTransport ();
  InitViscosity ();
  RhoStar      = CreatePolarGrid(NRAD, NSEC, "RhoStar");
  RhoInt       = CreatePolarGrid(NRAD, NSEC, "RhoInt");
  VradNew      = CreatePolarGrid(NRAD, NSEC, "VradNew");
  VradInt      = CreatePolarGrid(NRAD, NSEC, "VradInt");
  VthetaNew    = CreatePolarGrid(NRAD, NSEC, "VthetaNew");
  VthetaInt    = CreatePolarGrid(NRAD, NSEC, "VthetaInt");
  EnergyNew    = CreatePolarGrid(NRAD, NSEC, "EnergyNew");
  EnergyInt    = CreatePolarGrid(NRAD, NSEC, "EnergyInt");
  TemperInt    = CreatePolarGrid(NRAD, NSEC, "TemperInt");
  TempInt      = CreatePolarGrid(NRAD, NSEC, "TempInt");
  Potential    = CreatePolarGrid(NRAD, NSEC, "Potential");
  TurbPotential= CreatePolarGrid(NRAD, NSEC, "TurbPotential");
  Pressure     = CreatePolarGrid(NRAD, NSEC, "Pressure");
  SoundSpeed   = CreatePolarGrid(NRAD, NSEC, "SoundSpeed");
  Temperature  = CreatePolarGrid(NRAD, NSEC, "Temperature");
  Opacity      = CreatePolarGrid(NRAD, NSEC, "Opacity");
  ViscHeat     = CreatePolarGrid(NRAD, NSEC, "ViscousHeating");
  ThermHeat    = CreatePolarGrid(NRAD, NSEC, "ThermalHeating");
  ThermCool    = CreatePolarGrid(NRAD, NSEC, "ThermalCooling");
  RadDiffusion    = CreatePolarGrid(NRAD, NSEC, "RadiativeDiffusion");
  StarIrradiation = CreatePolarGrid(NRAD, NSEC, "StarIrradiation");
  ArtViscHeat     = CreatePolarGrid(NRAD, NSEC, "ArtViscousHeating");
  pdvEnergy       = CreatePolarGrid(NRAD, NSEC, "pdvEnergy");
  Test    = CreatePolarGrid(NRAD, NSEC, "Test");
  InitComputeAccel ();
  /* Rho and Energy are already initialized: cf main.c */
  if (ImposedAlpha)
    InitImposedAlpha ();
  if (InitEquilibrium)
    EqInitialize(Rho, Energy);
  ComputeSoundSpeed (Rho, Energy, sys);
  ComputePressureField (Rho, Energy);
  ComputeTemperatureField (Rho, Energy);
  InitGasVelocities (Vr, Vt, Rho);
  /* To output heating source terms at t=0... */
  if (EnergyEquation)
     ComputeOpacities (Rho, Energy);
  ComputeViscousTerms (Vr, Vt, Rho);
  if (ViscousHeating)
    ComputeViscousHeating (Rho);
  if (ThermalDiffusion)
    ComputeThermalDiffusion (Rho, Energy);
  if (ThermalCooling)
    ComputeThermalCooling (Rho, Energy);
  if (IrradStar)
    ComputeStarIrrad (Rho);
  if (RadiativeDiffusion)
    ComputeRadiavtiveDiffusion(Rho, Energy);
}

real min2 (a,b)
     real a,b;
{
  if (b < a) return b;
  return a;
}

real max2 (a,b)
     real a,b;
{
  if (b > a) return b;
  return a;
}


void ActualiseGas (array, newarray)
     PolarGrid *array, *newarray;
{
  int i,j,l,ns,nr;
  real *old, *new;
  nr = array->Nrad;
  ns = array->Nsec;
  old= array->Field;
  new= newarray->Field;
#pragma omp parallel for private(j,l)
  for (i = 0; i <= nr; i++) {
    for (j = 0; j < ns; j++) {
      l = j+ns*i;
      old[l] = new[l];
    }
  }
}


void AlgoGas (force, Rho, Vrad, Vtheta, Energy, Label, sys)
     Force *force;
     PolarGrid *Rho, *Vrad, *Vtheta, *Energy, *Label;
     PlanetarySystem *sys;
{
  real dthydro, dtnbody, dt, buf, dtemp=0.0;
  real xk, xj, yk, yj, mk, mj, dist;
  real OmegaNew, domega;
  int gastimestepcfl, k, j, NbPlanets;
  int ip;
  boolean Crashed=NO;
  extern boolean FargoPlanete, AccBoundary;
  extern real Runtime;
  extern real Mdisc0;
  real Mdiscnow, masscorenew, *Mswitch,rp;
  extern int dimfxy;
  extern boolean       Write_Sigdot;
  FirstGasStepFLAG=1;
  gastimestepcfl = 1;
  NbPlanets = sys->nb;
  Mswitch = (real *)malloc(NbPlanets*sizeof(real));
  if (EnergyEquation || ModifiedSoundSpeed) {
    ComputeSoundSpeed (Rho, Energy, sys);
    /* it is necessary to update calculation of soundspeed if one uses
       alphaviscosity in FViscosity. It is not necessary in locally
       isothermal runs since sounsspeed is constant. It is computed
       here for the needs of ConditionCFL. */
  }
  mpi_make1Dprofile (SoundSpeed->Field, axics);
  mpi_make1Dprofile (Rho->Field, axidens);
  mpi_make1Dprofile (Temperature->Field, axitemp);
  mpi_make1Dprofile (Opacity->Field, opaaxi);
  if (IsDisk == YES) {
    CommunicateBoundaries (Rho, Vrad, Vtheta, Energy, Label);
    if (SloppyCFL == YES)
      gastimestepcfl = ConditionCFL (Rho,Vrad, Vtheta, DT-dtemp);
  }
  MPI_Allreduce (&gastimestepcfl, &GasTimeStepsCFL, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  /* dthydro is the hydrodynamic timestep */
  dthydro = DT / (real)GasTimeStepsCFL;
  dtnbody = dthydro;
  if (NbPlanets > 1) {
    /* dtnbody is the n-body timestep */
    for (k = 0; k < NbPlanets; k++) {
      sys->mass[k] = FinalPlanetMass[k];
      for (j = k+1; j < NbPlanets; j++) {
       sys->mass[j] = FinalPlanetMass[j];
       xk = sys->x[k];
       xj = sys->x[j];
       yk = sys->y[k];
       yj = sys->y[j];
       mk = sys->mass[k];
       mj = sys->mass[j];
       dist = sqrt( (xk-xj)*(xk-xj) + (yk-yj)*(yk-yj) );
       buf = 2.0*M_PI*sqrt( dist*dist*dist / (mk+mj+1e-8) )/BINARYSECURITY;
       if ((mk > 0) && (mj > 0))
         dtnbody = min2(dtnbody,buf);
      }
    }
  }
  /* dt is the minimum between dthydro and dtnbody */
  dt = min2(dthydro,dtnbody);
  while (dtemp < 0.999999999*DT) {
    if (!FargoPlanete) {
      MassTaper = (PhysicalTime-PhysicalTimeInitial)/(MASSTAPER*2.0*M_PI);
      MassTaper = (MassTaper > 1.0 ? 1.0 : pow(sin(MassTaper*M_PI/2.0),2.0));
    } else {
      MassTaper = (PhysicalTime-PhysicalTimeInitial)/Runtime;
    }
    for (k = 0; k < NbPlanets; k++){
      if (!FargoPlanete){
        sys->mass[k] = PlanetMassAtRestart[k] + (FinalPlanetMass[k]-PlanetMassAtRestart[k])*MassTaper;
      } else {
        masscorenew = PlanetMassAtRestart[k] + (FinalPlanetMass[k]-PlanetMassAtRestart[k])*MassTaper;
        sys->mass[k] = masscorenew + Menvelope[k];
        if (FinalPlanetMass[k] == 0.0)
          sys->mass[k] = 0.0;
        rp = sqrt(pow(sys->x[k],2)+pow(sys->y[k],2));
        ip = ReturnIndex(rp);
        Mswitch[k] = MCRIFACTOR * pow(axics[ip]*pow(GlobalRmed[ip],0.5),3);
        if ((sys->mass[k] < Mswitch[k]) && !(sys->TorqueFlag[k]))
          sys->TorqueFlag[k] = YES;
      }
    }
    /* The current planet masses at t=PhysicalTime */
    dtnbody = dthydro;
    if (NbPlanets > 1) {
      for (k = 0; k < NbPlanets; k++) {
        for (j = k+1; j < NbPlanets; j++) {
          xk = sys->x[k];
          xj = sys->x[j];
          yk = sys->y[k];
          yj = sys->y[j];
          mk = sys->mass[k];
          mj = sys->mass[j];
          dist = sqrt( (xk-xj)*(xk-xj) + (yk-yj)*(yk-yj) );
          buf = 2.0*M_PI*sqrt( dist*dist*dist / (mk+mj+1e-8) )/BINARYSECURITY;
          if ((mk > 0) && (mj > 0))
            dtnbody = min2(dtnbody,buf);
            /* dtnbody is the n-body timestep */
        }
      }
    }
    if (IsDisk == YES) {
      CommunicateBoundaries (Rho, Vrad, Vtheta, Energy, Label);
      if (SloppyCFL == NO) {
        gastimestepcfl = 1;
        gastimestepcfl = ConditionCFL (Rho,Vrad, Vtheta, DT-dtemp);
        MPI_Allreduce (&gastimestepcfl, &GasTimeStepsCFL, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        /* dthydro is the hydrodynamic timestep */
        dthydro = (DT-dtemp)/(real)GasTimeStepsCFL;
      }
      AccreteOntoPlanets (Rho, Vrad, Vtheta, dt, sys);
    }
    /* dt is the minimum between dthydro and dtnbody */
    dt = min2(dthydro,dtnbody);
    dtemp += dt;
    DiskOnPrimaryAcceleration.x = 0.0;
    DiskOnPrimaryAcceleration.y = 0.0;
    if (Corotating == YES)
     GetPsysInfo (sys, MARK);
    if (IsDisk == YES) {
      /* Indirect term of star potential */
      DiskOnPrimaryAcceleration   = ComputeAccel (force, Rho, 0.0, 0.0, 0.0, 0.0, sys, 2);
      /* Gravitational potential from star and planet(s) */
      FillForcesArrays (sys, Rho, Energy, Vtheta, dt);
      /* Planets' velocities are updated with gravitationnal
      interaction with disk */
      AdvanceSystemFromDisk (force, Rho, Energy, sys, dt);
    }
    /* Planets' positions and velocities are updated with
    gravitational interaction with star and other planets */
    AdvanceSystemRK5 (sys, dt);
    /* Below we correct vtheta, the planet's position and velocities
    if we work in a frame non-centered on the primary */
    if (Corotating == YES) {
      OmegaNew = GetPsysInfo(sys, GET) / dt;
      domega = OmegaNew-OmegaFrame;
      if (IsDisk == YES)
        CorrectVtheta (Vtheta, domega);
      OmegaFrame = OmegaNew;
    }
    RotatePsys (sys, OmegaFrame*dt);
    /* Now we update gas fields */
    if (IsDisk == YES) {
      ApplyBoundaryCondition (Vrad, Vtheta, Rho, Energy, dt, sys);
      Crashed = DetectCrash (Rho);    /* test for negative density values */
      if (Crashed == YES) {
        if (AlreadyCrashed == 0) {
          timeCRASH=PhysicalTime;   /* if it appears to be the first crash */
          fprintf (stdout,"\nCrash in density! at time %.12g\n", timeCRASH);
          WriteDiskPolar (Rho, 999);    /* We write the HD arrays */
          WriteDiskPolar (Vrad, 999);   /* in order to keep a track */
          WriteDiskPolar (Vtheta, 999); /* of what happened */
          WriteDiskPolar (Energy, 999);
        }
        AlreadyCrashed++;
        masterprint ("c");
      }
      Crashed = DetectCrash (Energy);  /* test for negative energy values */
      if (Crashed == YES) {
        if (AlreadyCrashed == 0) {
          timeCRASH=PhysicalTime;   /* if it appears to be the first crash */
          fprintf (stdout,"\nCrash in energy! at time %.12g\n", timeCRASH);
          WriteDiskPolar (Rho, 999);    /* We write the HD arrays */
          WriteDiskPolar (Vrad, 999);   /* in order to keep a track */
          WriteDiskPolar (Vtheta, 999); /* of what happened */
          WriteDiskPolar (Energy, 999);
        }
        AlreadyCrashed++;
        masterprint ("c");
      } else {
        masterprint (".");
      }
      fflush (stdout);
      if (ZMPlus) {
        /* To model the non-axisymmetric component of the gas
        self-gravity with an anisotropic pressure (see aniso.c) */
        compute_anisotropic_pressurecoeff (sys);
      }
      /* Thermal diffusion needs to be applied first */
      if (EnergyEquation) {
        ComputeOpacities (Rho, Energy);
        if (ThermalDiffusion) {
          ComputeThermalDiffusion (Rho, Energy);
          SubStep0(Rho, Energy, dt);
        } else if (RadiativeDiffusion){
          ComputeRadiavtiveDiffusion(Rho, Energy);
          SubStep0(Rho, Energy, dt);
        }
      }
      ComputeTemperatureField (Rho, Energy);
      ComputeSoundSpeed (Rho, Energy, sys);
      ComputePressureField (Rho, Energy);
      /* Update vrad and vtheta with pressure, gravity and curvature
      source terms */
      SubStep1 (Vrad, Vtheta, Rho, sys, dt);
      /* Add some artifical viscosity */
      SubStep2 (Rho, Energy, dt);
      ActualiseGas (Vrad, VradNew);
      ActualiseGas (Vtheta, VthetaNew);
      if (EnergyEquation) {
        /* Update thermal energy with heating, cooling source terms */
        if (ViscousHeating) {
          ComputeViscousTerms (Vrad, Vtheta, Rho);
          ComputeViscousHeating (Rho);
        }
        if (ThermalCooling)
          ComputeThermalCooling (Rho, Energy);
        if (IrradStar)
          ComputeStarIrrad (Rho);
        SubStep3 (Rho, Vtheta, dt);
        ActualiseGas (Energy, EnergyNew);
      }
      ApplyBoundaryCondition (Vrad, Vtheta, Rho, Energy, dt, sys);
      /* Update velocities, surface density and thermal energy with
      advective terms */
      Transport (Rho, Vrad, Vtheta, Energy, Label, dt);
      ApplyBoundaryCondition (Vrad, Vtheta, Rho, Energy, dt, sys);
      /* Call to routine that reestablishes initial surface density on
      fixed timescale */
      if (AddMass) {
        DampDensity(Vrad, Vtheta, Rho, Energy, dt, sys);
        ApplyBoundaryCondition (Vrad, Vtheta, Rho, Energy, dt, sys);
      }
      if (DiscEvaporation) {
        Evaporation(Rho, dt);
        ApplyBoundaryCondition (Vrad, Vtheta, Rho, Energy, dt, sys);
      }
      if (PhotoEvapor){
        EvapMass += PhotoEvaporation(Vrad,Rho, dt);
        ApplyBoundaryCondition (Vrad, Vtheta, Rho, Energy, dt, sys);
      }
      /* Update just for outputs... */
      ComputeTemperatureField (Rho, Energy);
      /* Calculate mass (and mass excess) inside the planet's
      circumplanetary disk */
      mdcp = CircumPlanetaryMass (Rho, sys);
      exces_mdcp = mdcp - mdcp0;
      SetRhoFloor(Rho);
      if (EnergyEquation)
        SetEnergyFloor(Rho,Energy);
    }
    /* Check if the disc is gone */
    Mdiscnow =   GasTotalMass(Rho);
    if ((Mdiscnow/Mdisc0 <= 0.01) || (Rhole/GlobalRmed[GLOBALNRAD-1] >= 0.8)){
      SendOutput (TimeStep, Rho, Vrad, Vtheta, Energy, Label,sys);
      WritePlanetSystemFile (sys, TimeStep+1);
      WriteMassTrack (TimeStep+1, DiskMass, EvapMass, AccMassPls);
      if (Write_Sigdot)
        WriteSigmaDotFile(TimeStep+1);
      UpdateLog (force, sys, Rho, Energy, TimeStep+1, PhysicalTime, dimfxy);
      masterprint("Disc is gone\n");
      prs_exit(0);              
    }
    mpi_make1Dprofile (SoundSpeed->Field, axics);
    mpi_make1Dprofile (Rho->Field, axidens);
    mpi_make1Dprofile (Temperature->Field, axitemp);
    mpi_make1Dprofile (Opacity->Field, opaaxi);
    PhysicalTime += dt;
  }
  masterprint ("\n");
  free(Mswitch);
}


void SubStep0 (Rho,Energy,dt)
     PolarGrid *Rho;
     PolarGrid *Energy;
     real dt;
{
  int i, j, l, nr, ns;
  real *energy, *therheat, *dens, *raddiff;
  nr = Rho->Nrad;
  ns = Rho->Nsec;
  dens = Rho->Field;
  energy = Energy->Field;
  therheat = ThermHeat->Field;
  raddiff = RadDiffusion->Field;
    /* Update (explicite) with thermal diffusion */
  for (i = 0; i < nr; i++) {
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      if (ThermalDiffusion) {
        energy[l] += therheat[l]*dt;
      } else if (RadiativeDiffusion){
        energy[l] += raddiff[l]*dt;
      }
    }
  }
}


void SubStep1 (Vrad, Vtheta, Rho, sys, dt)
     PolarGrid *Vrad, *Vtheta, *Rho;
     PlanetarySystem *sys;
     real dt;
{
  int i, j, l, lim, ljm, ljp, nr, ns;
  boolean selfgravityupdate;
  real *vrad, *vtheta, *rho;
  real *Pot, *Press;
  real *vradint, *vthetaint;
  real gradp, gradphi, vt2, dxtheta;
  real invdxtheta;
  real supp_torque=0.0;              /* for imposed disk drift */
  nr = Vrad->Nrad;
  ns = Vrad->Nsec;
  rho = Rho->Field;
  vrad = Vrad->Field;
  vtheta = Vtheta->Field;
  vradint = VradInt->Field;
  vthetaint = VthetaInt->Field;
  Pot = Potential->Field;
  Press = Pressure->Field;
  /* In this substep we take into account the source terms of Euler
   equations. We update velocities with pressure gradients,
   gravitational forces and curvature terms */
#pragma omp parallel private(j,l,lim,ljm,ljp,dxtheta,vt2,gradp,gradphi,invdxtheta,supp_torque)
  {
#pragma omp for nowait
    for (i = 1; i < nr; i++) {
      for (j = 0; j < ns; j++) {
       l = j+i*ns;
       lim = l-ns;
       ljp = l+1;
       if (j == ns-1) ljp = i*ns;
       gradp = (Press[l]-Press[lim])*2.0/(rho[l]+rho[lim])*InvDiffRmed[i];
       gradphi = (Pot[l]-Pot[lim])*InvDiffRmed[i];
       vt2 = vtheta[l]+vtheta[ljp]+vtheta[lim]+vtheta[ljp-ns];
       vt2 = vt2/4.0+Rinf[i]*OmegaFrame;
       vt2 = vt2*vt2;
       vradint[l] = vrad[l]+dt*(-gradp-gradphi+vt2*InvRinf[i]);
      }
    }
#pragma omp for
    for (i = 0; i < nr; i++) {
      supp_torque = IMPOSEDDISKDRIFT*0.5*pow(Rmed[i],-2.5+SIGMASLOPE);
      dxtheta = (PMAX-PMIN)/(real)ns*Rmed[i];
      invdxtheta = 1.0/dxtheta;
      for (j = 0; j < ns; j++) {
       l = j+i*ns;
       ljm = l-1;
       if (j == 0) ljm = i*ns+ns-1;
       gradp = (Press[l]-Press[ljm])*2.0/(rho[l]+rho[ljm])*invdxtheta;
       if ( ZMPlus ) {
         /* To model the non-axisymmetric component of the gas
            self-gravity with an anisotropic pressure (see aniso.c) */
         gradp *= SG_aniso_coeff;
       }
       gradphi = (Pot[l]-Pot[ljm])*invdxtheta;
       vthetaint[l] = vtheta[l]-dt*(gradp+gradphi);
       vthetaint[l] += dt*supp_torque;
      }
    }
  }
   /* Update velocities with gas self-gravity */
  if ( SelfGravity ) {
    selfgravityupdate = YES;
    compute_selfgravity(Rho, VradInt, VthetaInt, dt, selfgravityupdate);
  }
  
  /* Update velocities with gas viscosity */
  ComputeViscousTerms (VradInt, VthetaInt, Rho);
  UpdateVelocitiesWithViscosity (VradInt, VthetaInt, Rho, dt);
  if ( !DontApplySubKeplerian ) {
    ApplySubKeplerianBoundary (VthetaInt, sys);
  }
}

void SubStep2 (Rho, Energy, dt)
     PolarGrid *Rho, *Energy;
     real dt;
{
  int i, j, l, lim, lip, ljm, ljp, nr, ns;
  real *vrad, *vtheta, *rho, *energy, *artvischeat;
  real *vradnew, *vthetanew, *qt, *qr, *energyint;
  real dxtheta, invdxtheta;
  real dv;
  nr = Rho->Nrad;
  ns = Rho->Nsec;
  rho = Rho->Field;
  vrad = VradInt->Field;
  vtheta = VthetaInt->Field;
  qr = RhoInt->Field;
  vradnew = VradNew->Field;
  vthetanew = VthetaNew->Field;
  qt = TemperInt->Field;
  energy = Energy->Field;
  energyint = EnergyInt->Field;
  artvischeat = ArtViscHeat->Field;
#pragma omp parallel for private(j,dxtheta,l,lim,lip,ljm,ljp,dv)
  for (i = 0; i < nr; i++) {
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      lip = l+ns;
      ljp = l+1;
      if (j == ns-1) ljp = i*ns;
      dv = vrad[lip]-vrad[l];
      if (dv < 0.0)
        qr[l] = CVNR*CVNR*rho[l]*dv*dv;
      else 
        qr[l] = 0.0; 
      dv = vtheta[ljp]-vtheta[l];
      if (dv < 0.0)
        qt[l] = CVNR*CVNR*rho[l]*dv*dv;
      else
       qt[l] = 0.0;
    }
  }
#pragma omp parallel private(l,lim,lip,ljm,ljp,j,dxtheta,invdxtheta)
  {
#pragma omp for nowait
    for (i = 1; i < nr; i++) {
      for (j = 0; j < ns; j++) {
       l = j+i*ns;
       lim = l-ns;
       vradnew[l] = vrad[l]-dt*2.0/(rho[l]+rho[lim])*(qr[l]-qr[lim])*InvDiffRmed[i];
      }
    }
#pragma omp for
    for (i = 0; i < nr; i++) {
      dxtheta = (PMAX-PMIN)/(real)ns*Rmed[i];
      invdxtheta = 1.0/dxtheta;
      for (j = 0; j < ns; j++) {
       l = j+i*ns;
       ljm = l-1;
       if (j == 0) ljm = i*ns+ns-1;
       vthetanew[l] = vtheta[l]-dt*2.0/(rho[l]+rho[ljm])*(qt[l]-qt[ljm])*invdxtheta;
      }
    }
    if (EnergyEquation) {
#pragma omp for nowait
      for (i = 0; i < nr; i++) {
       dxtheta = (PMAX-PMIN)/(real)ns*Rmed[i];
       invdxtheta = 1.0/dxtheta;
       for (j = 0; j < ns; j++) {
         l = j+i*ns;
         lip = l+ns;
         ljp = l+1;
         if (j == ns-1) ljp = i*ns;
         energyint[l] = energy[l] -                            \
           dt*qr[l]*(vrad[lip]-vrad[l])*InvDiffRsup[i] -       \
           dt*qt[l]*(vtheta[ljp]-vtheta[l])*invdxtheta;
         artvischeat[l] = -qr[l]*(vrad[lip]-vrad[l])*InvDiffRsup[i] -       \
           qt[l]*(vtheta[ljp]-vtheta[l])*invdxtheta;
       }
      }
    }
  }
}
              
void SubStep3 (Rho, Vtheta, dt)
     PolarGrid *Rho, *Vtheta;
     real dt;
{
  extern boolean TempPresc, BetaCooling, IrradStar;
  int i, j, l, nr, ns;
  real *energy, *energynew, *dens, *divergence, *vischeat, *therheat, *thercool, *vtheta, *opacity, *soundspeed, *stirrad, *pdvenergy;
  real num, den, omega, beta, coolingtime, horizontaltau, tau, h, epsilon = 0.5, taueff;
  nr = Rho->Nrad;
  ns = Rho->Nsec;
  dens = Rho->Field;
  vtheta = Vtheta->Field;
  energy = EnergyInt->Field;
  energynew = EnergyNew->Field;
  divergence = DivergenceVelocity->Field;
  vischeat = ViscHeat->Field;
  therheat = ThermHeat->Field;
  thercool = ThermCool->Field;
  stirrad = StarIrradiation->Field;
  opacity = Opacity->Field;
  soundspeed = SoundSpeed->Field;
  pdvenergy = pdvEnergy->Field;
  /* In this substep, we update the gas thermal energy with source
     terms (compression/dilatation, viscous heating, thermal
     diffusion, temperature prescription) */
#pragma omp parallel private(j,l,num,den)
  {
#pragma omp for
    /* Update (implicite) with pdV only */
    for (i = 0; i < nr; i++) {
      for (j = 0; j < ns; j++) {
       l = j+i*ns;
       num = energy[l];
       den = 1.0+(ADIABATICINDEX-1.0)*dt*divergence[l];
       energynew[l] = num/den;
       pdvenergy[l] = energynew[l]-energy[l];
      }
    }
    /* Update (explicite) with viscous heating (Q+ term, see
       e.g. d'Angelo et al. 03). Viscous heating array already 
       calculated before call to substep3 routine. */
    if (ViscousHeating) {
      for (i = 0; i < nr; i++) {
       for (j = 0; j < ns; j++) {
         l = j+i*ns;
         energynew[l] += vischeat[l]*dt;
       }
      }
    }
    /* Stellar irradiation as in Muller&Kley2013 but devided by opacity to be consistent with Pierens2016*/
    if (IrradStar){  
      for (i = 0; i < nr; i++) {
        for (j = 0; j < ns; j++) {
         l = j+i*ns;
          energynew[l] += stirrad[l]*dt; 
        }
      }
    }
    if (ThermalCooling) {
      /* Update (explicite) with thermal cooling (beta bersion). Note
        that this array is defined as positive. Thermal cooling array
        already calculated before call to substep3 routine. */
      for (i = 0; i < nr; i++) {
       for (j = 0; j < ns; j++) {
         l = j+i*ns;
         energynew[l] -= thercool[l]*dt;
       }
      }
    }
    /* Here cooling is modeled through a temperature prescription,
       with characteristic time is PrescTimeMed(r). Energy update is
       implicite. */
    if (TempPresc) {
      for (i = 0; i < nr; i++) {
       for (j = 0; j < ns; j++) {
         l = j+i*ns;
         num = energynew[l] + EnergyMed[i]*(dens[l]/SigmaMed[i])*(dt/PrescTimeMed[i]);
         den = 1.0+dt/PrescTimeMed[i]; 
         energynew[l] = num/den;
       }
      }
    }
    /* Case of a simple beta-cooling, where the cooling source term
       reads -e/tau, with tau = beta/Omega */
    if (BetaCooling) {
      for (i = 0; i < nr; i++) {
       omega = 0.0;
       for (j = 0; j < ns; j++) {
         l = j+i*ns;
         omega += (vtheta[l] + Rmed[i]*OmegaFrame);
       }
       omega /= (real)ns;
       omega /= Rmed[i];
       beta  = ComputeBetaCooling(Rmed[i]);
       for (j = 0; j < ns; j++) {
         l = j+i*ns;
         num = energynew[l];
         coolingtime = beta / omega;
         den = 1.0+dt/coolingtime; 
         energynew[l] = num/den;
       }
      }
    }
    /* --- */
  }
}

int ConditionCFL (Rho,Vrad, Vtheta, deltaT)
     PolarGrid *Rho, *Vrad, *Vtheta;
     real deltaT;
{
  static real Vresidual[MAX1D], Vmoy[MAX1D];
  int i, j, l, ns, nr, lip, ljp;
  real invdt1, invdt2, invdt3, invdt4, invdt5, invdt6, cs, newdt, dt;
  int ideb, jdeb;
  real itdbg1, itdbg2, itdbg3, itdbg4, itdbg5, itdbg6, mdtdbg; /* debugging variables */
  real *vt, *vr, dxrad, dxtheta, dvr, dvt, viscr, visct;
  real *soundspeed, *opacity, *dens, *temperature, chi, h, Rl, lambda, gradT, invdphi;
  soundspeed = SoundSpeed->Field;
  dens = Rho->Field;
  temperature = Temperature->Field;
  opacity = Opacity->Field;
  ns = Vtheta->Nsec;
  nr = Vtheta->Nrad;
  vt = Vtheta->Field;
  vr = Vrad->Field;
  newdt = 1e30;
  invdphi = (real)ns/(PMAX-PMIN);
  for (i = 0; i < nr; i++) {
    dxrad = Rsup[i]-Rinf[i];
    dxtheta=Rmed[i]*(PMAX-PMIN)/(real)ns;
    Vmoy[i] = 0.0;
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      Vmoy[i] += vt[l];
    }
    Vmoy[i] /= (real)ns;
  }
  for (i = One_or_active; i < Max_or_active; i++) {
    dxrad = Rsup[i]-Rinf[i];
    dxtheta=Rmed[i]*(PMAX-PMIN)/(real)ns;
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      if (FastTransport == YES)
       Vresidual[j] = vt[l]-Vmoy[i];  /* FARGO algorithm */
      else
       Vresidual[j] = vt[l];              /* Standard algorithm */
    }
    Vresidual[ns]=Vresidual[0];
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      lip = l+ns;
      ljp = l+1;
      if (j == ns-1) ljp = i*ns;
      cs = soundspeed[l];
      invdt1 = cs/(min2(dxrad,dxtheta));
      invdt2 = fabs(vr[l])/dxrad;
      invdt3 = fabs(Vresidual[j])/dxtheta;
      dvr = vr[lip]-vr[l];
      dvt = vt[ljp]-vt[l];
      if (dvr >= 0.0) dvr = 1e-10;
      else dvr = -dvr;
      if (dvt >= 0.0) dvt = 1e-10;
      else dvt = -dvt;
      invdt4 = max2(dvr/dxrad,dvt/dxtheta);
      invdt4*= 4.0*CVNR*CVNR;
      if ( ViscosityAlpha || (VISCOSITY != 0.0) ) 
       invdt5 = FViscosity(Rmed[i],cs)*4.0/pow(min2(dxrad,dxtheta),2.0);
      else 
       invdt5 = 1e-10;
      if (ThermalDiffusion){ 
       invdt6 = DIFFUSIVITY*4.0/pow(min2(dxrad,dxtheta),2.0);
      } else if (RadiativeDiffusion){
        h = cs/pow(Rmed[i],-0.5);
  /* calculating lambda for cell l*/
      gradT = sqrt(pow(InvDiffRsup[i]*(temperature[lip]-temperature[l]),2)+ pow(InvRmed[i]*invdphi*(temperature[ljp]-temperature[l]),2));
      Rl = 8*h*Rmed[i]/dens[l]/opacity[l]*gradT/temperature[l];
      if (Rl <= 2){
        lambda = 2./(3.*sqrt(9+10*Rl*Rl));
      } else {
        lambda = 10./sqrt(10.*Rl+9+sqrt(180.*Rl+81.));
      }
        chi  = lambda*64.*ADIABATICINDEX*(ADIABATICINDEX-1)*sigma_SB*pow(temperature[l],4)/ (opacity[l]*pow(dens[l],2)*pow(Rmed[i],-3));
        invdt6 = chi*4.0/pow(min2(dxrad,dxtheta),2.0);
      } else {
       invdt6 = 1e-10;
      }
      dt = CFLSECURITY/sqrt(invdt1*invdt1+invdt2*invdt2+invdt3*invdt3+invdt4*invdt4+invdt5*invdt5+invdt6*invdt6);
      if (dt < newdt) {
       newdt = dt;
       if (debug == YES) {
         ideb = i;
         jdeb = j;
         itdbg1 = 1.0/invdt1; itdbg2=1.0/invdt2; itdbg3=1.0/invdt3; itdbg4=1.0/invdt4; itdbg5=1.0/invdt5; itdbg6=1.0/invdt6;
         mdtdbg = newdt;
         viscr = dxrad/dvr/4.0/CVNR/CVNR;     
         visct = dxtheta/dvt/4.0/CVNR/CVNR;
       }
      }  
    }
  }
  for (i = Zero_or_active; i < MaxMO_or_active; i++) {
    dt = (PMAX-PMIN)*CFLSECURITY/(real)NSEC/fabs(Vmoy[i]*InvRmed[i]-Vmoy[i+1]*InvRmed[i+1]);
    if (dt < newdt) newdt = dt;
  }
  if (deltaT < newdt) newdt = deltaT;
  if (debug == YES) {
    printf ("Timestep control information for CPU %d: \n", CPU_Rank);
    printf ("Most restrictive cell at i=%d and j=%d\n", ideb, jdeb);
    printf ("located at radius Rmed         : %g\n", Rmed[ideb]);
    printf ("Sound speed limit              : %g\n", itdbg1);
    printf ("Radial motion limit            : %g\n", itdbg2);
    printf ("Residual circular motion limit : %g\n", itdbg3);
    printf ("Artificial viscosity limit     : %g\n", itdbg4);
    printf ("Viscosity limit                : %g\n", itdbg5);
    printf ("Thermal Diffusivity limit      : %g\n", itdbg6);
    printf ("   Arise from r with limit     : %g\n", viscr);
    printf ("   and from theta with limit   : %g\n", visct);
    printf ("Limit time step for this cell  : %g\n", mdtdbg);
    printf ("Limit time step adopted        : %g\n", newdt);
    if (newdt < mdtdbg) {
      printf ("Discrepancy arise either from shear.\n");
      printf ("or from the imposed DT interval.\n");
    }
  }
  return (int)(ceil(deltaT/newdt));
}

void ComputeViscousHeating (Rho)
     PolarGrid *Rho;
{
  int i, j, l, nr, ns;
  int lip, li2p;
  real r, rip, ri2p, qpip, qpi2p, viscosity;
  real *dens, *divergence, *Trr, *Trp, *Tpp, *vischeat, *cs;
  nr = Rho->Nrad;
  ns = Rho->Nsec;
  dens = Rho->Field;
  divergence = DivergenceVelocity->Field;
  vischeat = ViscHeat->Field;
  Trr = TAURR->Field;
  Trp = TAURP->Field;
  Tpp = TAUPP->Field;
  cs  = SoundSpeed->Field; 
  /* We calculate the heating source term from i=1 */
  for (i = 1; i < nr; i++) {     /* Trp defined from i=1 */
    if (!ViscosityAlpha && (VISCOSITY == 0.0) )
      viscosity = 0.0;
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      if (ViscosityAlpha || (VISCOSITY != 0.0) )
        viscosity = FViscosity (Rmed[i], cs[l]);
      if (viscosity != 0.0) {
       vischeat[l] = 0.5/viscosity/dens[l]*( Trr[l]*Trr[l] +              \
                                         2.0*Trp[l]*Trp[l] +       \
                                         Tpp[l]*Tpp[l] );
       vischeat[l] += (2.0/9.0)*viscosity*dens[l]*divergence[l]*divergence[l];
      }
      else
       vischeat[l] = 0.0;
    }
  }
  /* We calculate the heating source term Vischeat for i=0 */
  i = 0;
  r    = Rmed[i];
  rip  = Rmed[i+1];
  ri2p = Rmed[i+2];
  for (j = 0; j < ns; j++) {
    l = j+i*ns;
    lip = l+ns;
    li2p = lip+ns;
    qpip = vischeat[lip];   // vischeat(i=1,j)
    qpi2p = vischeat[li2p]; // vischeat(i=2,j)
    if (viscosity != 0.0) {
      // power-law extrapolation
    vischeat[l] = qpip*exp( log(qpip/qpi2p) * log(r/rip) / log(rip/ri2p) );
    }
    else
      vischeat[l] = 0.0;
  }
}

void ComputeOpacities (Rho, Energy)
     PolarGrid *Rho, *Energy;
{
  int i, j, l, nr, ns;
   real *energ;
   real *dens, rho3D, phys_dens;
   real *opacity;
   real *test;
   real temp, phys_temp;
   real roversigma, buf;
   real temp_transition_34, temp_transition_45, temp_transition_56, temp_transition_67, temp_transition_78;
   dens = Rho->Field;
   energ = Energy->Field;
   opacity = Opacity->Field;
   test = Test->Field;
   nr = Rho->Nrad;
   ns = Rho->Nsec;
   for (i = 0; i < nr; i++) {
     for (j = 0; j < ns; j++) {
       l = i*ns + j;
       /* Convert code temperature into Kelvins */
       temp = (ADIABATICINDEX-1.0)*energ[l]*pow(dens[l],-1.0);  // temperature in code units
       phys_temp = temp * unit_temperature;
       /* Convert 3D volume density into g.cm^-3 */
       roversigma = Rmed[i] / dens[l];
       buf = (ADIABATICINDEX-1.0)*energ[l]*pow(roversigma,3.);
       rho3D = 1./sqrt(buf*2*M_PI);  // 3D density = sigma / sqrt(2*pi) H, 
       phys_dens = rho3D * unit_mass * pow(unit_length, -3.);  // in kg.m^(-3)
       phys_dens *= 1e-3;  // in g.cm^(-3)
       opacity[l] = opLBL94 (phys_dens, phys_temp);
//       opacity[l] =  OpacityOriginal(phys_dens, phys_temp);
       test[l] = opacity[l];   // in cm^2 / g
       /* We convert opacities them in m^2 / kg, before translating
         the result into code units */
       opacity[l] *= (0.1 * pow(unit_length,-2.0) * pow(unit_mass,1.0));
       opacity[l] *= ZMETAL; 
     }
   }
}

real OpacityOriginal(phys_dens, phys_temp)
   real phys_dens, phys_temp;
{
  real opac, temp_transition_34, temp_transition_45, temp_transition_56;
  real temp_transition_67, temp_transition_78;
  /* Opacities are calculated in cm^2 / g from Bell and Lin (94)
    tables  */
       if ( phys_temp < 167.0 )
        opac = 2e-4*pow(phys_temp,2.0);
       else {
        if ( phys_temp < 203.0 )
          opac = 2e16*pow(phys_temp,-7.0);
        else {
          temp_transition_34 = pow(2e82*phys_dens,2./49);
          if ( phys_temp < temp_transition_34 ) 
            opac = 0.1*pow(phys_temp,0.5);
          else {
            temp_transition_45 = pow(2e89*pow(phys_dens,1./3),1./27);
            if ( phys_temp < temp_transition_45 )
              opac = 2e81*pow(phys_dens,1.0)*pow(phys_temp,-24.);
            else {
              temp_transition_56 = pow(1e28*pow(phys_dens,1./3),1./7);
              if ( phys_temp < temp_transition_56 )
                     opac = 1e-8*pow(phys_dens,2./3)*pow(phys_temp,3.);
              else {
           temp_transition_67 = pow(1.5e56*pow(phys_dens,2./3),0.08);
                     if ( phys_temp < temp_transition_67 )
                       opac = 1e-36*pow(phys_dens,1./3)*pow(phys_temp,10.);
                       else {
                         temp_transition_78 = pow(4.31e20*phys_dens,2./5);
                         if ( phys_temp < temp_transition_78 )
                           opac = 1.5e20*pow(phys_dens,1.)*pow(phys_temp,-2.5);
                         else
                           opac = 0.348;
                       }
               }
             }
           }
         }
       }
  return opac;
}

real opLBL94 (rho3d, temp)
   real rho3d, temp;
{
   real power[3] = {2.8369e-2, 1.1464e-2, 2.2667e-1};
   real trans[3] = {1.46e3, 4.51e3, 2.37e6};
//   coefficients for opacity laws 1, 2, and 3 in cgs units.
   real ak[3] = {2.e-4, 2.e16, 0.1};
//   coefficients for opacity laws 3, 4, 5, 6, 7, and 8 in T_4 units.
   real bk[6] = {10.,2.e-15,1e4,1.e4,1.5e10,0.348};
   real ts4, rho13, rho23, ts42, ts44, ts48, t2, t4, t8, t10;
   real o1, o2, o3, o4, o5, o6, o7, o8, o1an, o2an, o3an, o4an, o6an, o7an, o8an ;
   real opacity;
/* Opacities are calculated in cm^2 / g from Bell and Lin (94)
         tables with smoothing as Lin and  Papaloizou 1985 */
   if(temp > trans[0]*pow(rho3d,power[0])){
      ts4=1.e-4*temp;
      rho13=pow(rho3d,0.333333333);
      rho23=rho13*rho13;
      ts42=ts4*ts4;
      ts44=ts42*ts42;
      ts48=ts44*ts44;
      if(temp > trans[1]*pow(rho3d,power[1])){
          if((temp < trans[2]*pow(rho3d,power[2])) || ((rho3d < 1e10) && (temp < 1e4))){
             o5=bk[2]*rho23*ts42*ts4;
             o6=bk[3]*rho13*ts44*ts44*ts42;
             o7=bk[4]*rho3d/(ts42*sqrt(ts4));
             o6an=o6*o6;
             o7an=o7*o7;
             opacity = pow(o6an*o7an/(o6an+o7an),2);
             opacity += pow(o5/(1+pow(ts4/1.1/pow(rho3d,0.04762),10)),4);
             opacity = pow(opacity,0.25);
            } else {
             o7=bk[4]*rho3d/(ts42*sqrt(ts4));
             o8=bk[5];       
             o7an=o7*o7;
             o8an=o8*o8;
             opacity = pow(o7an*o7an+o8an*o8an,0.25);
            }
          return opacity;
      }
      o3=bk[0]*sqrt(ts4);
      o4=bk[1]*rho3d/(ts48*ts48*ts48);
      o5=bk[2]*rho23*ts42*ts4;
      o4an=pow(o4,4);
      o3an=pow(o3,4);
      opacity = pow((o4an*o3an/(o4an+o3an))+pow(o5/(1+6.561e-5/ts48*1e2*rho23),4),0.25);
    } else {
       t2=temp*temp;
       t4=t2*t2;
       t8=t4*t4;
       t10=t8*t2;
       o1=ak[0]*t2;
       o2=ak[1]*temp/t8;
       o3=ak[2]*sqrt(temp);
       o1an=o1*o1;
       o2an=o2*o2;
       opacity = pow(pow(o1an*o2an/(o1an+o2an),2)+pow(o3/(1+1.e22/t10),4),0.25);
    }
   return opacity;
}


void ComputeThermalCooling (Rho, Energy)
     PolarGrid *Rho, *Energy;
{
  int i, j, l, nr, ns;
  real *dens, *energ, *thercool, *opacity;
  real temp, tau, tau_eff, temp_irr;
  dens = Rho->Field;
  energ = Energy->Field;
  thercool = ThermCool->Field;
  opacity = Opacity->Field;
  nr = Rho->Nrad;
  ns = Rho->Nsec;
  for (i = 0; i < nr; i++) {
    for (j = 0; j < ns; j++) {
      l = i*ns + j;
      tau = opacity[l]*dens[l] / sqrt(2*M_PI); 
      tau_eff = 0.375*tau + 0.25*sqrt(3.0) + 0.25/(tau+1e-20); // effective optical depth
      temp = (ADIABATICINDEX-1.0)*energ[l]/dens[l];  // temperature
      if (!StellarIrradiation){
       thercool[l] = 2.0*sigma_SB*pow(temp,4.)*pow(tau_eff,-1.0);
      } else {
       temp_irr = (160/unit_temperature) * pow(Rmed[i]*unit_length/1.49598e11,-0.5); // assuming Tirr = 160K x (R/1AU)^-1/2
       thercool[l] = 2.0*sigma_SB*(pow(temp,4.)-pow(temp_irr,4.))*pow(tau_eff,-1.0);
      }
    }
  }
}

void ComputeRadiavtiveDiffusion (Rho, Energy)
     PolarGrid *Rho, *Energy;
{
  int i, j, l, nr, ns;
  int lip, lim, ljp, ljm, liip, ljjp, ljjm;
  real *energy, *dens, *raddiff, *temp, *cs, *opacity;
  real dphi, invdphi;
  real Omega, Kip, Kl, Kjp, Kim, Kjm, Hscale, laplacienT, lambda, Rl, gradT;
  nr = Rho->Nrad;
  ns = Rho->Nsec;
  energy = Energy->Field;
  dens = Rho->Field;
  temp =  Temperature->Field;
  raddiff = RadDiffusion->Field;
  cs = SoundSpeed->Field;
  opacity = Opacity->Field;
  dphi = (PMAX-PMIN)/(real)ns;
  invdphi = 1.0/dphi;
  /* Radiative diffusion implemented as in Muller & Kley 2013
     and is -2H div(K gradT) */
  /* calculating K gradT */
  for (i = 1; i < nr-2; i++) {
    Omega = pow(Rmed[i],-1.5);
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      lip = l+ns;
      lim = l-ns;
      ljp = l+1;
      if (j == ns-1) ljp = i*ns;
      ljm = l-1;
      if (j == 0) ljm = i*ns+ns-1;
      liip = lip+ns;
      ljjp = ljp+ns;
      if (j == ns-2) ljjp = i*ns;
      ljjm = ljm-1;
      if (j == 1) ljjm = i*ns+ns-1;
      Hscale = cs[l]/Omega;
  /* calculating lambda for cell l and Kl*/
      gradT = sqrt(pow(InvDiffRsup[i]*(temp[lip]-temp[l]),2)+ pow(InvRmed[i]*invdphi*(temp[ljp]-temp[l]),2));
      Rl = 8*cs[l]/pow(Rmed[i],-1.5)/dens[l]/opacity[l]*gradT/temp[l];
      if (Rl <= 2){
        lambda = 2./(3.*sqrt(9+10*Rl*Rl));
      } else {
        lambda = 10./sqrt(10.*Rl+9+sqrt(180.*Rl+81.));
      }
      Kl = -lambda*16.*sigma_SB*2*Hscale*pow(temp[l],3) / (dens[l]*opacity[l]);
  /* calculating lambda for cell lip*/
      gradT = sqrt(pow(InvDiffRsup[i+1]*(temp[liip]-temp[lip]),2)+ pow(InvRmed[i+1]*invdphi*(temp[ljp+ns]-temp[lip]),2));
      Rl = 8*cs[lip]/pow(Rmed[i+1],-1.5)/dens[lip]/opacity[lip]*gradT/temp[lip];
      if (Rl <= 2){
        lambda = 2./(3.*sqrt(9+10*Rl*Rl));
      } else {
        lambda = 10./sqrt(10.*Rl+9+sqrt(180.*Rl+81.));
      }
      Kip = -lambda*16.*sigma_SB*2*cs[lip]*pow(temp[lip],3) / pow(Rmed[i+1],-1.5)/ (dens[lip]*opacity[lip]);
  /* calculating lambda for cell ljp*/
      gradT = sqrt(pow(InvDiffRsup[i]*(temp[ljp+ns]-temp[ljp]),2)+ pow(InvRmed[i]*invdphi*(temp[ljjp]-temp[ljp]),2));
      Rl = 8*cs[ljp]/pow(Rmed[i],-1.5)/dens[ljp]/opacity[ljp]*gradT/temp[ljp];
      if (Rl <= 2){
        lambda = 2./(3.*sqrt(9+10*Rl*Rl));
      } else {
        lambda = 10./sqrt(10.*Rl+9+sqrt(180.*Rl+81.));
      }
      Kjp = -lambda*16.*sigma_SB*2*cs[ljp]*pow(temp[ljp],3) / Omega / (dens[ljp]*opacity[ljp]);
  /* calculating lambda for cell lim*/
      gradT = sqrt(pow(InvDiffRsup[i-1]*(temp[l]-temp[lim]),2)+ pow(InvRmed[i-1]*invdphi*(temp[lim]-temp[ljm-ns]),2));
      Rl = 8*cs[lim]/pow(Rmed[i-1],-1.5)/dens[lim]/opacity[lim]*gradT/temp[lim];
      if (Rl <= 2){
        lambda = 2./(3.*sqrt(9+10*Rl*Rl));
      } else {
        lambda = 10./sqrt(10.*Rl+9+sqrt(180.*Rl+81.));
      }
      Kim = -lambda*16.*sigma_SB*2*cs[lim]*pow(temp[lim],3) / Omega / (dens[lim]*opacity[lim]);
  /* calculating lambda for cell ljm*/
      gradT = sqrt(pow(InvDiffRsup[i]*(temp[ljm+ns]-temp[ljm]),2)+ pow(InvRmed[i]*invdphi*(temp[ljm]-temp[ljjm]),2));
      Rl = 8*cs[ljm]/pow(Rmed[i],-1.5)/dens[ljm]/opacity[ljm]*gradT/temp[ljm];
      if (Rl <= 2){
        lambda = 2./(3.*sqrt(9+10*Rl*Rl));
      } else {
        lambda = 10./sqrt(10.*Rl+9+sqrt(180.*Rl+81.));
      }
      Kjm = -lambda*16.*sigma_SB*2*cs[ljm]*pow(temp[ljm],3) / Omega / (dens[ljm]*opacity[ljm]);
      raddiff[l] = InvRmed[i]*InvDiffRsup[i]*(Rsup[i]*0.5*(Kl+Kip)*(temp[lip]-temp[l])*InvDiffRmed[i+1]\
-Rsup[i-1]*0.5*(Kl+Kim)*(temp[l]-temp[lim])*InvDiffRmed[i]);
      raddiff[l] += InvRmed[i]*InvRmed[i]*invdphi*invdphi*(0.5*(Kjp+Kl)*(temp[ljp]-temp[l])\
-0.5*(Kl+Kjm)*(temp[l]-temp[ljm]));
      raddiff[l] *= -2.*Hscale;
    }
  }
}

void ComputeThermalDiffusion (Rho, Energy)
     PolarGrid *Rho, *Energy;
{
  int i, j, l, nr, ns;
  int lip, lim, ljp, ljm;
  real *energy, *tempint, *dens, *therheat;
  real dphi, invdphi, laplacien, buf;
  nr = Rho->Nrad;
  ns = Rho->Nsec;
  /* tempint is an intermediary array: energy / density */
  tempint = TempInt->Field;
  energy = Energy->Field;
  dens = Rho->Field;
  therheat = ThermHeat->Field;
  dphi = (PMAX-PMIN)/(real)ns;
  invdphi = 1.0/dphi;
  /* Thermal diffusion implemented as in Paardekooper, Baruteau & Kley
     2011 (see their Eq. 1) with actually entropy diffusion. This
     avoids the background disc structure to evolve
     significantly... */
  for (i = 0; i < nr; i++) {
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      buf = (ADIABATICINDEX-1.0)*energy[l]/pow(dens[l],ADIABATICINDEX);
      if (buf <= 0.0) {
       tempint[l] = 0.0;
      } 
      if (buf > 0.0) {
       tempint[l] = log( buf );
      }
      //tempint[l] = energy[l]/dens[l];
    }
  }
  for (i = 1; i < nr-1; i++) {
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      lip = l+ns;
      lim = l-ns;
      ljp = l+1;
      if (j == ns-1) ljp = i*ns;
      ljm = l-1;
      if (j == 0) ljm = i*ns+ns-1;
      laplacien = InvDiffRsup[i]*InvRmed[i]*(
                                        Rsup[i]*InvDiffRmed[i]*(tempint[lip]-tempint[l]) - \
                                        Rinf[i]*InvDiffRmed[i]*(tempint[l]-tempint[lim]) 
                                        ) +                     \
       InvRmed[i]*InvRmed[i]*invdphi*invdphi*(tempint[ljp]+tempint[ljm]-2.0*tempint[l]);
      therheat[l] = energy[l]*DIFFUSIVITY*laplacien;
    }
  }
}

void ComputeStarIrrad (Rho)
     PolarGrid *Rho;
{
  int i, j, l, nr, ns;
  real *dens, *stirrad, *opacity, *soundspeed;
  real h, tau, taueff, qirrad, epsilon=0.5;
  dens = Rho->Field;
  stirrad = StarIrradiation->Field;
  opacity = Opacity->Field;
  soundspeed = SoundSpeed ->Field;
  nr = Rho->Nrad;
  ns = Rho->Nsec; 
  for (i = 0; i < nr; i++) {
  /* a factor 2 is multiplied to count for both sides of the disc*/
     qirrad = sigma_SB * (tstar*tstar*tstar*tstar) * (rstar*rstar/Rmed[i]/Rmed[i]);
     for (j = 0; j < ns; j++) {
       l = j+i*ns;
       h = soundspeed[l] * sqrt(Rmed[i]);
       tau = opacity[l]*dens[l]/sqrt(2*M_PI);
       taueff = 3./8*tau+sqrt(3.)/4.+1./(4*tau+1e-20);
       /* the qirrad is calculated considering grazing angle as in Pierens2016 */
       stirrad[l] = 2.*(qirrad*(1-epsilon) * h * 2./7.)/taueff; //factor of 2 is because of disc illumination in both sides
     }
   }
}

void ComputeSoundSpeed (Rho, Energy, sys)
     PolarGrid *Rho;
     PolarGrid *Energy;
     PlanetarySystem *sys;
{
  int i, j, k, l, nr, ns, NbPlanets;
  real *dens, *energ, *cs, *abs, *ord;
  real num, den, h, rsmoothing;
  real xp, yp, mp, xc, yc, dx, dy, d2, dp;
  real omsq, omplasq, omeff;
  real csiso, cstarget;
  nr = Rho->Nrad;
  ns = Rho->Nsec;
  dens = Rho->Field;
  energ = Energy->Field;
  cs = SoundSpeed->Field;
  abs = CellAbscissa->Field;
  ord = CellOrdinate->Field;
  NbPlanets = sys->nb;
  for ( i = 0; i < nr; i++ ) {
    if (ModifiedSoundSpeed) {
      csiso = AspectRatio(Rmed[i])*sqrt(G*1.0/Rmed[i])*pow(Rmed[i], FLARINGINDEX);
      rsmoothing = compute_smoothing(Rmed[i]);
    }
    for ( j = 0; j < ns; j++ ) {
      l = i*ns + j;
      if (!EnergyEquation) {
       if (!ModifiedSoundSpeed) {
         cs[l] = AspectRatio(Rmed[i])*sqrt(G*1.0/Rmed[i])*pow(Rmed[i], FLARINGINDEX); 
       } else {
         /* NEW: sound speed expression in Peplinski et al. 08 and in
            Lin & Papaloizou 2011/12. We slowly increase the sound
            speed from its locally isothermal expression (csiso(r)) to
            its 'target' expression given e.g. by Eq. (5) in Lin &
            Papaloizou 2012. Tapering timescale and procedure
            identical as that of the planet's potential */
         // effective sound speed with planet index 0
         for (k = 0; k < NbPlanets; k++) {
           xp = sys->x[k];
           yp = sys->y[k];
           mp = FinalPlanetMass[k];
           xc = abs[l];
           yc = ord[l];
           dx = xc-xp;
           dy = yc-yp;
           d2 = dx*dx + dy*dy + rsmoothing*rsmoothing;
           // smoothed distance from planet
           dp = sqrt(d2);
           // disc aspect ratio at r=Rmed[i]
           h = AspectRatio(Rmed[i])*pow(Rmed[i], FLARINGINDEX);
           num = h*Rmed[i] * PLANETASPECTRATIO*dp;
           den = pow(pow(h*Rmed[i],7./2) + pow(PLANETASPECTRATIO*dp,7./2),2./7);
           omsq = G*1.0*pow(Rmed[i],-3.0);
           omplasq = G*mp*pow(dp,-3.0);
           // Effective omega squared
           omeff = sqrt(omsq + omplasq);
           if (k==0) 
             cstarget = num*omeff/den;
           else
             cstarget += num*omeff/den;
         }
         cstarget /= (NbPlanets+0.0);
         cs[l] = csiso + (cstarget-csiso)*MassTaper;
       }
      } else {
       cs[l] = sqrt( ADIABATICINDEX*(ADIABATICINDEX-1.0)*energ[l]/dens[l] );
     }
    }
  }
}

void ComputePressureField (Rho, Energy)
     PolarGrid *Rho;
     PolarGrid *Energy;
{
  int i, j, l, nr, ns;
  real *dens, *pres, *energ, *cs;
  nr = Rho->Nrad;
  ns = Rho->Nsec;
  dens = Rho->Field;
  pres = Pressure->Field;
  energ = Energy->Field;
  cs = SoundSpeed->Field;
  for ( i = 0; i < nr; i++ ) {
    for ( j = 0; j < ns; j++ ) {
      l = i*ns + j;
      if (!EnergyEquation) {
       pres[l] = dens[l]*cs[l]*cs[l]; /* since SoundSpeed is not updated */
                                       /* from initialization, cs remains */ 
                                       /* axisymmetric */
      }
      else
       pres[l] = (ADIABATICINDEX-1.0)*energ[l];
    }
  }
}

void ComputeTemperatureField (Rho, Energy)
     PolarGrid *Rho, *Energy;
{
  int i, j, l, nr, ns;
  real *dens, *pres, *energ, *temp;
  nr = Rho->Nrad;
  ns = Rho->Nsec;
  dens = Rho->Field;
  pres = Pressure->Field;
  energ = Energy->Field;
  temp = Temperature->Field;
  for ( i = 0; i < nr; i++ ) {
    for ( j = 0; j < ns; j++ ) {
      l = i*ns + j;
      if (!EnergyEquation)
       temp[l] = MU/R* pres[l]/dens[l];
      else
       temp[l] = MU/R*(ADIABATICINDEX-1.0)*energ[l]/dens[l];
    }
  }
}

real CircumPlanetaryMass (Rho, sys)
     PolarGrid *Rho;
     PlanetarySystem *sys;
{
  int i, j, l, ns;
  real xpl, ypl, rpl;
  real dist, mdcplocal, mdcptotal, MyHillRadius;
  real *dens, *abs, *ord;
  ns = Rho->Nsec;
  dens = Rho->Field;
  abs = CellAbscissa->Field;
  ord = CellOrdinate->Field;
  xpl = sys->x[0];
  ypl = sys->y[0];
  rpl = sqrt( xpl*xpl + ypl*ypl );
  mdcplocal = 0.0;
  mdcptotal = 0.0;
  MyHillRadius = rpl * pow( sys->mass[0]/3., 1./3. );
  if (FakeSequential && (CPU_Rank > 0)) 
    MPI_Recv (&mdcplocal, 1, MPI_DOUBLE, CPU_Rank-1, 0, MPI_COMM_WORLD, &fargostat);
  for ( i = Zero_or_active; i < Max_or_active; i++ ) {
    for ( j = 0; j < ns; j++ ) {
      l = i*ns + j;
      dist = sqrt ( (abs[l]-xpl)*(abs[l]-xpl) +              \
                  (ord[l]-ypl)*(ord[l]-ypl) );
      if ( dist < MyHillRadius ) {
       mdcplocal += Surf[i] * dens[l];
      }
    }
  }
  if (FakeSequential) {
     if (CPU_Rank < CPU_Number-1)
       MPI_Send (&mdcplocal, 1, MPI_DOUBLE, CPU_Rank+1, 0, MPI_COMM_WORLD);
  } else {
    MPI_Allreduce (&mdcplocal, &mdcptotal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  }
  if (FakeSequential) {
    MPI_Bcast (&mdcplocal, 1, MPI_DOUBLE, CPU_Number-1, MPI_COMM_WORLD);
    mdcptotal = mdcplocal;
  }
  return mdcptotal;
}
