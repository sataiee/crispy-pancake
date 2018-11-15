/** \file Pframeforce.c

Functions that evaluate the %force between the planets and the disk.
The FillForcesArrays() function is ill-named: it rather fill an array
of the potential, that is later use to derive the force acting on the
disk at every zone.  The name of this file is due to the fact that we
work in the frame centered on the primary (which is therefore not
inertial). Old versions of fargo also featured the file Gframeforce.c,
in which we worked in the frame centered on the center of gravity of
the system.  The present file also contains the functions necessary to
update the planets' position and velocities (taking into account, or
not, the indirect term, ie the term arising from the fact that the
frame is not inertial), as well as the function that initializes the
hydrodynamics fields with analytic prescription.

*/

#include "mp.h"

extern boolean AllowAccretion, Indirect_Term;
extern Pair DiskOnPrimaryAcceleration;
static Pair IndirectTerm;
static real q0[MAX1D], q1[MAX1D], PlanetMasses[MAX1D];
static real vt_int[MAX1D], vt_cent[MAX1D];
extern boolean OneDRun;

void ComputeIndirectTerm () {
  IndirectTerm.x = -DiskOnPrimaryAcceleration.x;
  IndirectTerm.y = -DiskOnPrimaryAcceleration.y; 
  if (Indirect_Term == NO) {
    IndirectTerm.x = 0.0;
    IndirectTerm.y = 0.0;
  }
}
/* Below : work in non-rotating frame */
/* centered on the primary */
void FillForcesArrays (sys, Rho, Energy, Vtheta, dt)
     PlanetarySystem *sys;
     PolarGrid *Rho, *Energy, *Vtheta;
     real dt;
{
  int i, j, l, nr, ns, k, NbPlanets;
  real x, y, angle, distance, distancesmooth;
  real xplanet, yplanet, RRoche,smooth, mplanet;
  real PlanetDistance, *Pot, pot=0, smoothing;
  //real *test;
  real InvPlanetDistance3, InvDistance;
  real *cs;
  extern boolean MHDLSA;
  Pot= Potential->Field;
  cs = SoundSpeed->Field;
  nr = Potential->Nrad;
  ns = Potential->Nsec;
  NbPlanets = sys->nb;
  //test = Test->Field;

  /* Indirect term star on gas here */
  ComputeIndirectTerm ();

#pragma omp parallel for
  for (i = 0; i < (nr+1)*ns; i++) Pot[i] = 0.0;
  //for (i = 0; i < (nr+1)*ns; i++) test[i] = 0.0;
  /* -- Gravitational potential from planet on gas -- */
  for (k = 0; k < NbPlanets; k++) {
    xplanet = sys->x[k];
    yplanet = sys->y[k];
    mplanet = sys->mass[k];
    PlanetDistance = sqrt(xplanet*xplanet+yplanet*yplanet);
    InvPlanetDistance3 =  1.0/PlanetDistance/PlanetDistance/PlanetDistance;
    RRoche = PlanetDistance*pow((1.0/3.0*mplanet),1.0/3.0);
    if (RocheSmoothing)
      smoothing = RRoche*ROCHESMOOTHING;
    #pragma omp parallel for private(InvDistance,j,l,angle,x,y,distance,distancesmooth,pot)
    for (i = 0; i < nr; i++) {
      InvDistance = 1.0/Rmed[i];
      for (j = 0; j < ns; j++) {
        l = j+i*ns;
        angle = azimuth[j];
        x = Rmed[i]*cos(angle);
        y = Rmed[i]*sin(angle);
        if (!RocheSmoothing)
          smoothing = cs[l]*pow(x*x+y*y, 1.5)/sqrt(ADIABATICINDEX) * THICKNESSSMOOTHING ;
        smooth = smoothing*smoothing;
        distance = (x-xplanet)*(x-xplanet)+(y-yplanet)*(y-yplanet);
        distancesmooth = sqrt(distance+smooth);
        /* If the planet's mass is smaller than the critical mass and the run 
        * is 1D, the planet's potential is not applied on the disc.
        * Otherwise, the rings close to the planet, feel the planet so strong and 
        * therefore, the density close to the planet increased greatly. */
        if (!OneDRun) {
          pot = -G*mplanet/distancesmooth; /* Direct term from planet */
          if (Indirect_Term == YES)
            pot += G*mplanet*InvPlanetDistance3*(x*xplanet+y*yplanet); /* Indirect term from planet  */
        }
        Pot[l] += pot;
      }
    }
  }
  /* -- Gravitational potential from star on gas -- */
#pragma omp parallel for private(InvDistance,j,l,angle,x,y,pot)
  for (i = 0; i < nr; i++) {
    InvDistance = 1.0/Rmed[i];
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      angle = azimuth[j];
      x = Rmed[i]*cos(angle);
      y = Rmed[i]*sin(angle);
      pot = -G*1.0*InvDistance;  /* Direct term from star */
      /* case where azimuthal extent equals 2pi */
      if ( fabs(PMAX-PMIN-2.*M_PI) < 0.01 )
        pot -= IndirectTerm.x*x + IndirectTerm.y*y; /* Indirect term from star */
      Pot[l] += pot;  
    }
  }
  /* -- Turbulent potential to modelize the MHD turbulence driven by
     MRI, as lin Laughlin et al. (2004, LSA04) -- */
  if (MHDLSA)
    ApplyLSAOnPotential();
}

void AdvanceSystemFromDisk (force, Rho, Energy, sys, dt)
     Force *force;
     PlanetarySystem *sys;
     PolarGrid *Rho, *Energy;
     real dt;           
{
  int NbPlanets, k, ip;
  Pair gamma, accel;
  real x, y, r, m, smoothing;
  NbPlanets = sys->nb;
  for (k = 0; k < NbPlanets; k++) {
    if (sys->FeelDisk[k] == YES) {
      m = sys->mass[k];
      x = sys->x[k];
      y = sys->y[k];
      r = sqrt(x*x + y*y);
      ip = ReturnIndex(r);
      if (RocheSmoothing)
        smoothing = r*pow(m/3.,1./3.)*ROCHESMOOTHING;
      else
        smoothing = compute_smoothing (r);
      if (sys->mass[k] > 0.0) {
        if (sys->TorqueFlag[k] == YES){
          if (CPU_Master){
            accel = AccelFromFormula (force, Rho, x, y, smoothing, m, sys, k,1);             
            sys->vx[k] += dt * accel.x;
            sys->vy[k] += dt * accel.y;
          }
          MPI_Barrier (MPI_COMM_WORLD);
        } else {
          gamma = ComputeAccel (force, Rho, x, y, m, sys);
          sys->vx[k] += dt * gamma.x;
          sys->vy[k] += dt * gamma.y;
        }
      }
      sys->vx[k] += dt * IndirectTerm.x;
      sys->vy[k] += dt * IndirectTerm.y;
    }
  }
}

void AdvanceSystemRK5 (sys, dt)
     PlanetarySystem *sys;
     real dt;
{
  extern boolean ForcedCircular, ForcedInnerCircular;
  int i, n, myimin;
  boolean *feelothers;
  real dtheta, omega, rdot, x, y, r, v, new_r, vx, vy, theta, denom;
  n = sys->nb;
  if (!ForcedCircular) {
    for (i = 0; i < n; i++) {
      q0[i] = sys->x[i];
      q0[i+n] = sys->y[i];
      q0[i+2*n] = sys->vx[i];
      q0[i+3*n] = sys->vy[i];
      PlanetMasses[i] = sys->mass[i];
    }
    feelothers = sys->FeelOthers;
    RungeKunta (q0, dt, PlanetMasses, q1, n, feelothers);
  }
  /* Default case (see below) */
  if (!ForcedInnerCircular) {
    for (i = 1-(PhysicalTime >= RELEASEDATE); i < sys->nb; i++) {
      /* Default case: planets position and velocity updated after 
      Runge Kutta step */
      if (!ForcedCircular) {
        sys->x[i] = q1[i];
        sys->y[i] = q1[i+n];
        sys->vx[i] = q1[i+2*n];
        sys->vy[i] = q1[i+3*n];
      } else {
        /* Case where planets are held on a fixed circular orbit with 
        initial angular frequency omega */
        x = sys->x[i];
        y = sys->y[i];
        theta = atan2(y,x);
        vx = sys->vx[i];
        vy = sys->vy[i];
        r = sqrt(x*x + y*y);
        v = sqrt(vx*vx + vy*vy);
        omega = (-y*vx + x*vy)/r/r;
        dtheta = omega*dt;
        sys->x[i]  = r*cos(theta+dtheta);
        sys->y[i]  = r*sin(theta+dtheta);
        sys->vx[i] = -v*sin(theta+dtheta);
        sys->vy[i] =  v*cos(theta+dtheta);
      }
    }
  } else {
    /* New (july 2012): particular case where inner planet held on a fixed 
    circular orbit */
    for (i = 0; i < n; i++) {
      if (i == 0) {  // inner planet (i=0) fixed -> copy-paste of above
        x = sys->x[i];
        y = sys->y[i];
        theta = atan2(y,x);
        vx = sys->vx[i];
        vy = sys->vy[i];
        r = sqrt(x*x + y*y);
        v = sqrt(vx*vx + vy*vy);
        omega = (-y*vx + x*vy)/r/r;
        dtheta = omega*dt;
        sys->x[i]  = r*cos(theta+dtheta);
        sys->y[i]  = r*sin(theta+dtheta);
        sys->vx[i] = -v*sin(theta+dtheta);
        sys->vy[i] =  v*cos(theta+dtheta);
      } else {  // all planets except that indexed with i=0
        sys->x[i] = q1[i];
        sys->y[i] = q1[i+n];
        sys->vx[i] = q1[i+2*n];
        sys->vy[i] = q1[i+3*n];
      }
    }
  }
  /* Case where the innermost planet (with index 0) is drifted
     manually with a prescribed migration rate tuned by RELEASERADIUS 
     and RELEASETIME in .par file */
  if (PhysicalTime < RELEASEDATE) {
    x = sys->x[0];
    y = sys->y[0];
    r = sqrt(x*x+y*y);
    theta = atan2(y,x);
    rdot = (RELEASERADIUS-r)/(RELEASEDATE-PhysicalTime);
    omega = sqrt((1.+sys->mass[0])/r/r/r);
    new_r = r + rdot*dt;
    denom = r-new_r;
    if (denom != 0.0) {
      dtheta = 2.*dt*r*omega/denom*(sqrt(r/new_r)-1.);
    } else {
      dtheta = omega*dt;
    }
    vx = rdot;
    vy = new_r*sqrt((1.+sys->mass[0])/new_r/new_r/new_r);
    sys->x[0] = new_r*cos(dtheta+theta);
    sys->y[0] = new_r*sin(dtheta+theta);
    sys->vx[0]= vx*cos(dtheta+theta) - vy*sin(dtheta+theta); 
    sys->vy[0]= vx*sin(dtheta+theta) + vy*cos(dtheta+theta); 
  }
}

void SolveOrbits (sys,l)
     PlanetarySystem *sys;
     int l;
{
  int i, n;
  real x, y, vx, vy;
  char name[256];
  FILE *output;
  n = sys->nb;
  for (i = 0; i < n; i++) {
    x = sys->x[i];
    y = sys->y[i];
    vx = sys->vx[i];
    vy = sys->vy[i];
    FindOrbitalElements (sys, x, y, vx, vy, 1.0+sys->mass[i], i, l);
  }
}  

real ConstructSequence (u, v, n)
     real *u, *v;
     int n;
{
  int i;
  real lapl=0.0;
  for (i = 1; i < n; i++)
    u[i] = 2.0*v[i]-u[i-1];
  for (i = 1; i < n-1; i++)
    lapl += fabs(u[i+1]+u[i-1]-2.0*u[i]);
  return lapl;
}

void InitGasDensity (Rho)
     PolarGrid *Rho;
{
  int i, j, l, nr, ns;
  real *dens, randomnb;
  extern boolean ImposedDensity, AddNoise;
  dens = Rho->Field;
  nr = Rho->Nrad;
  ns = Rho->Nsec;
  if (!ImposedDensity) {
    /* Standard case with power-law density profile */
    FillSigma ();
    for (i = 0; i < nr; i++) {
      for (j = 0; j < ns; j++) {
        l = j+i*ns;
        dens[l] = SigmaMed[i];
        /* No random noise is added by default to the initial density
        and velocity profiles. If AddNoise set to yes, white noise
        added to the initial density field with arbitrary 1d-3
        relative amplitude */
        if (AddNoise) {
          randomnb = 2.0*drand48()-1.0;
          dens[l] += 1e-3*SigmaMed[i]*randomnb;
        }
      }
    }
  } else {
    /* NEW Nov 2010: we impose a fixed density profile */
    InitImposedDensity (Rho);
  }
}

void InitImposedDensity (density)
     PolarGrid *density;
{
  int i, ig, j, l, lg, ns, nr;
  real *dens, foo, value;
  FILE *DENSFILE;
  char name_dens[1024];
  real *globaldens;
  sprintf (name_dens, "%saxidens.dat", OUTPUTDIR);
  DENSFILE = fopen (name_dens, "r");
  dens = density->Field;
  nr = density->Nrad;
  ns = density->Nsec;
  if (DENSFILE == NULL)
    erreur ("ERROR: I could not read the file axidens.dat containing the imposed density profile. Please check and run again\n");
    
  globaldens = (real*) malloc(sizeof(real)*ns*GLOBALNRAD);
  for (i = 0; i < GLOBALNRAD; i++) {
    fscanf (DENSFILE, "%lf %lf", &foo, &value);
    for (j = 0; j < ns; j++) {
      l = i*ns + j;
      globaldens[l] = (real)value;
    }
  }
  for (i = 0; i < nr; i++) {
    ig = (i+IMIN)*ns;
    SigmaMed[i] = globaldens[ig];
    if (i > 0)
      SigmaInf[i] = 0.5*(SigmaMed[i]+SigmaMed[i-1]);
    else
      SigmaInf[i] = SigmaMed[i];
    for (j = 0; j < ns; j++) {
      l = i*ns + j;
      lg = (i+IMIN)*ns + j;
      dens[l] = globaldens[lg];
    }
  }
  fclose (DENSFILE);
}

void InitImposedAlpha ()
{
  int i;
  real foo, value;
  FILE *ALPHAFILE;
  char name_alpha[1024];
  sprintf (name_alpha, "%saxialpha.dat", OUTPUTDIR);
  ALPHAFILE = fopen (name_alpha, "r");
  if (ALPHAFILE == NULL) {
    erreur ("ERROR: I could not read the file axialpha.dat containing the imposed alpha viscosity profile. Please check and run again\n");
  }
  for (i = 0; i < GLOBALNRAD; i++) {
    fscanf (ALPHAFILE, "%lf %lf", &foo, &value);
    GLOBAL_ImposedAlpha[i] = (real)value;
  }
  fclose (ALPHAFILE);
}

void InitGasEnergy (energ)
     PolarGrid *energ;
{
  int i, j, l, nr, ns;
  real *energy;
  extern boolean TempPresc;
  energy = energ->Field;
  nr = energ->Nrad;
  ns = energ->Nsec;
  FillEnergy ();
  for (i = 0; i < nr; i++) {
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      energy[l] = EnergyMed[i];
    }
  }
  if (TempPresc)
    FillPrescTime();
}

void InitGasVelocities (Vr, Vt, Rho)
     PolarGrid *Vr, *Vt, *Rho;
{
  extern boolean SGZeroMode, AccBoundary;
  extern boolean SelfGravity, InitEquilibrium, MdotHartmann;
  int i, j, l, nr, ns;
  real *vr, *vt, *pres, *cs, *dens, mdot;
  real  r, omega, ri, moy;
  real viscosity, t1, t2, r1, r2;
  real axipres[GLOBALNRAD];
  vr  = Vr->Field;
  vt  = Vt->Field;
  nr  = Vt->Nrad;
  ns  = Vt->Nsec;
  cs = SoundSpeed->Field;
  pres = Pressure->Field;  /* Pressure is already initialized: cf initeuler in SourceEuler.c ... */
  dens = Rho->Field;
  if (MdotHartmann){
    mdot = 1e-8 * pow(THARTMANN/1e6, -1.4); //Hartmann1998, modified to give a smaller Mdot because our disc is evolved
    mdot *= -(1.9891e30/31556926.0 / unit_mass*unit_time); //convert to code unit
  } else {
    mdot  = -MDOTINIT;
  }
  /* --------- */
  // Initialization of azimutal velocity with exact centrifugal balance
  /* --------- */
  if ( CentrifugalBalance ) {
    /* vt_int \equiv rOmega = grad(P)/sigma +  \partial(phi)/\partial(r)  -  acc_sg_radial */
    mpi_make1Dprofile (pres, axipres);
    /* global axisymmetric pressure field, known by all cpus*/
    for (i = 1; i < GLOBALNRAD; i++) {
      vt_int[i] = ( axipres[i] - axipres[i-1] ) /  \
        (.5*(Sigma(GlobalRmed[i])+Sigma(GlobalRmed[i-1])))/(GlobalRmed[i]-GlobalRmed[i-1]) + \
        G*(1.0/GlobalRmed[i-1]-1.0/GlobalRmed[i])/(GlobalRmed[i]-GlobalRmed[i-1]);
    }
    /* Case of a disk with self-gravity */
    if ( SelfGravity ) { // Better test with CL rigid!
      if ( !SGZeroMode )
        mpi_make1Dprofile (SG_Accr, GLOBAL_AxiSGAccr);
      else
        GLOBAL_AxiSGAccr = SG_Accr;
      for (i = 1; i < GLOBALNRAD; i++)
        vt_int[i] -= ( (Radii[i] - GlobalRmed[i-1])*GLOBAL_AxiSGAccr[i] + \
          (GlobalRmed[i] - Radii[i])*GLOBAL_AxiSGAccr[i-1] ) / (GlobalRmed[i]-GlobalRmed[i-1]);
    }
    for (i = 1; i < GLOBALNRAD; i++)
      vt_int[i] = sqrt(vt_int[i]*Radii[i])-Radii[i]*OmegaFrame;

    t1 = vt_cent[0] = vt_int[1]+.75*(vt_int[1]-vt_int[2]);
    r1 = ConstructSequence (vt_cent, vt_int, GLOBALNRAD);
    vt_cent[0] += .25*(vt_int[1]-vt_int[2]);
    t2 = vt_cent[0];
    r2 = ConstructSequence (vt_cent, vt_int, GLOBALNRAD);
    t1 = t1-r1/(r2-r1)*(t2-t1);
    vt_cent[0] = t1;
    ConstructSequence (vt_cent, vt_int, GLOBALNRAD);
    vt_cent[GLOBALNRAD] = vt_cent[GLOBALNRAD-1];
  }
  /* --------- */
  // Initialization with self-gravity, without exact centrifugal balance
  if (SelfGravity && !CentrifugalBalance)
    init_azimutalvelocity_withSG (Vt);
  /* --------- */

  for (i = 0; i <= nr; i++) {
    if (i == nr) {
      r = Rmed[nr-1];
      ri= Rinf[nr-1];
    } else {
      r = Rmed[i];
      ri= Rinf[i];
    }
    if (!ViscosityAlpha && (VISCOSITY == 0.0) )
      viscosity = 0.0;
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      if (ViscosityAlpha || (VISCOSITY != 0.0) )
        viscosity = FViscosity (r, cs[l], dens[l]);
      /* --------- */
      if (!SelfGravity) {
        omega = sqrt(G*1.0/r/r/r);
        vt[l] = omega*r*sqrt(1.0-pow(ASPECTRATIO,2.0)*      \
            pow(r,2.0*FLARINGINDEX)*      \
            (1.+SIGMASLOPE-2.0*FLARINGINDEX) );
      }
      /* --------- */
      vt[l] -= OmegaFrame*r;
      if (CentrifugalBalance)
        vt[l] = vt_cent[i+IMIN];
      if (i == nr) {
        vr[l] = 0.0;
      } else {
        vr[l] = IMPOSEDDISKDRIFT*SIGMA0/SigmaInf[i]/ri;
        if (AccBoundary){ 
          vr[l] -= 3.0*viscosity/r/2.;
        } else {
          if (ViscosityAlpha) 
            vr[l] -= 3.0*viscosity/r*(-SIGMASLOPE+2.0*FLARINGINDEX+1.0);
            //vr[l] -= 3.0*viscosity/2./r;
          else 
            vr[l] -= 3.0*viscosity/r*(-SIGMASLOPE+.5);
        }
        if (InitEquilibrium)
          vr[l] = mdot/2/PI/dens[l]/ri;  
      }
    }
  }
  if (InitEquilibrium){
    moy = 0.0;
    for (j = 0; j < ns; j++){
      vr[j+ns*nr] = mdot/2/PI/dens[j+ns*(nr-1)]/Rsup[nr-1];  
      moy += vr[j+ns*nr];
    }
    VradMed[GLOBALNRAD] = moy/ns;
  }

   mpi_make1Dprofile (vr, VradMed);
   mpi_make1Dprofile (vt, VthetaMed); 
}
  
