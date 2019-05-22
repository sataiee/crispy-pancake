/** \file Init.c

Contains the functions needed to initialize the hydrodynamics arrays.
These can be initialized by reading a given output (in the case of a
restart) or by calling a function, InitEuler (), which contains
analytic prescription for the different hydrodynamics fields. Note
that this function InitEuler() is located in SourceEuler.c, which
itself calls InitGas(), in the file Pframeforce.c.
Also, note that the present file contains InitLabel(), which sets
the initial value of a passive scalar.
*/

#include "mp.h"

extern boolean Restart;
extern int     NbRestart;

void ReadfromFile (array, fileprefix, filenumber)
     PolarGrid *array;
     char *fileprefix;
     int filenumber;
{
  int nr, ns, c, foo=0;
  real *field;
  char name[256];
  FILE *input;
  /* Simultaneous read access to the same file have been observed to
     give wrong results. */
  /* A sequential reading is imposed below. */
  /* If current CPU has a predecessor, wait for a message from him */
  if (CPU_Rank > 0) MPI_Recv (&foo, 1, MPI_INT, CPU_Prev, 10, MPI_COMM_WORLD, &fargostat);
  sprintf (name, "%s%s%d.dat", OUTPUTDIR, fileprefix, filenumber);
  input = fopen (name, "r");
  if (input == NULL) {
    fprintf (stderr, "WARNING ! Can't read %s. Restarting with t=0 settings.\n", name); 
    if (CPU_Rank < CPU_Highest) MPI_Send (&foo, 1, MPI_INT, CPU_Next, 10, MPI_COMM_WORLD);
    return;
  }
  field = array->Field;
  nr = array->Nrad;
  ns = array->Nsec;
  for (c = 0; c < IMIN; c++) {
    fread (field, sizeof(real), ns, input); 
    /* Can't read at once in order not to overflow 'field' */
  }
  fread (field, sizeof(real), nr*ns, input);
  fclose (input);
  /* Next CPU is waiting. Tell it to start now by sending the message
     that it expects */
  if (CPU_Rank < CPU_Highest) MPI_Send (&foo, 1, MPI_INT, CPU_Next, 10, MPI_COMM_WORLD);
  MPI_Barrier (MPI_COMM_WORLD);       /* previous CPUs do not touch anything
                               meanwhile */
}

void InitLabel (array, sys)
     PolarGrid *array;
     PlanetarySystem *sys;
{
  int nr, ns, i, j, l;
  real xp, yp, rp;
  real x, y, angle, distance, rhill;
  real *field;
  field = array->Field;
  nr = array->Nrad;
  ns = array->Nsec;
  xp = sys->x[0];
  yp = sys->y[0];
  rp = sqrt ( xp*xp + yp*yp );
  rhill = rp * pow( sys->mass[0]/3., 1./3 );
  /* Initialize label as you wish. In this example, label only takes
     into account fluid elements inside the planet's Hill Sphere */
  for (i = 0; i < nr; i++) {
    for (j = 0; j < ns; j++) {
      l = j + i*ns;
      angle = azimuth[j];
      x = Rmed[i] * cos(angle);
      y = Rmed[i] * sin(angle);
      distance = sqrt( (x-xp)*(x-xp) + (y-yp)*(y-yp) );
      if ( distance < rhill )
       field[l] = 1.0;
      else
       field[l] = 0.0;
    }
  }
}

void Initialization (gas_density, gas_v_rad, gas_v_theta, gas_energy, gas_label, pla_sys, SGAarray)
     PolarGrid *gas_density, *gas_v_rad, *gas_v_theta, *gas_energy, *gas_label;
     PlanetarySystem *pla_sys;
     real *SGAarray;
{
  extern boolean EnergyEquation, ThermalDiffusion, ThermalCooling, FargoPlanete;
  extern int NbRestart;
  real *energ, *dens;
  FILE *output;
  char OutputName[256];
  int i, j, l, nr, ns;
  nr = gas_energy->Nrad;
  ns = gas_energy->Nsec;
  ReadPrevDim ();
  InitEuler (gas_v_rad, gas_v_theta, gas_density, gas_energy, pla_sys, SGAarray);
  InitLabel (gas_label, pla_sys);
  if (Restart == YES) {
    CheckRebin (NbRestart);
    /* Now that OldRmed is built, Cpu_master writes again Radii in file
    used_rad.dat */
    if ( CPU_Master ) {
      sprintf (OutputName, "%s%s", OUTPUTDIR, "used_rad.dat");
      output = fopen (OutputName, "w");
      if (output == NULL) {
        mastererr ("Can't write %s.\nProgram stopped.\n", OutputName);
        prs_exit (1);
      }
      for (i = 0; i <= GLOBALNRAD; i++)
        fprintf (output, "%.18g\n", Radii[i]);
      fclose (output);
    }
    MPI_Barrier (MPI_COMM_WORLD);
    
    /* Don't start reading before master has finished rebining... */
    /* It shouldn't be a problem though since a sequential read is */
    /* imposed in the ReadfromFile function below */
    mastererr ("Reading restart files...\n");
    fflush (stderr);
    ReadfromFile (gas_density, "gasdens", NbRestart);
    ReadfromFile (gas_v_rad, "gasvrad", NbRestart);
    ReadfromFile (gas_v_theta, "gasvtheta", NbRestart);
    if (EnergyEquation) {
      ReadfromFile (gas_energy, "gasTemperature", NbRestart);
      /* ! gas_energy accounts for the gas temperature... */
      energ = gas_energy->Field;
      dens = gas_density->Field;
      for (i=0; i<nr; i++) {
        for (j=0; j<ns; j++) {
          l = i*ns + j;
          energ[l] = dens[l]*energ[l]/(ADIABATICINDEX-1.0);
          /* this is e = dens*temp / (gamma-1) */
        }
      }
    }
    if (FargoPlanete){
      mpi_make1Dprofile (gas_v_rad->Field, VradMed);
      mpi_make1Dprofile (gas_v_theta->Field, VthetaMed); 
      RefillSigma (gas_density);
      RefillEnergy (gas_energy);
    }
    /* To output restart fields correctly */
    if (NbRestart == 0){
      mpi_make1Dprofile (gas_v_rad->Field, VradMed);
      mpi_make1Dprofile (gas_v_theta->Field, VthetaMed); 
      RefillSigma (gas_density);
      RefillEnergy (gas_energy);
    }
    ComputeSoundSpeed (gas_density, gas_energy, pla_sys);
    ComputePressureField (gas_density, gas_energy);
    ComputeTemperatureField (gas_density, gas_energy);
    if (ThermalDiffusion)
      ComputeThermalDiffusion (gas_density, gas_energy);
    if (ThermalCooling)
      ComputeThermalCooling (gas_density, gas_energy);
    ComputeViscousTerms (gas_v_rad, gas_v_theta, gas_density);
    ComputeViscousHeating (gas_density);
    ReadfromFile (gas_label, "gaslabel", NbRestart);
    if (EnergyEquation)
      ComputeOpacities (gas_density, gas_energy);
    if (StoreSigma)
     RefillSigma (gas_density);
    if (StoreEnergy)
     RefillEnergy (gas_energy);
    fprintf (stderr, "done\n");
    fflush (stderr);
  }
  make_azi_average_profile (SoundSpeed->Field, axics);
  make_azi_average_profile (gas_density->Field, axidens);
  make_azi_average_profile (Temperature->Field, axitemp);
  WriteDim (); 
}

/* This function calculates the initial temperature and density 
   for a disc with constant accration rate when energy equarion is present
   This function is modified to start from Bitsch2014 profiles (Jan 2017)*/

void EqInitialize(Rho, Energy)
  PolarGrid *Energy, *Rho;
{
  extern boolean DecInner, AccBoundary, MdotHartmann, OpInner;
  int nr, ns, i, j, l, n;
  real *temp, *energ, *dens;
  real Terr, Tnew, Told, dsmoothin=0.1, mdot;
  real epsilon = 1e-3, foo, value;
  nr = Energy->Nrad;
  ns = Energy->Nsec;
  energ = Energy->Field;
  dens = Rho->Field;
  temp = (real*) malloc(sizeof(real)*nr);
  if (MdotHartmann){
    mdot = 1e-8 * pow(THARTMANN/1e6, -1.4); //Hartmann1998, modified to give a smaller Mdot because our disc is evolved
    mdot *= -(1.9891e30/31556926.0 / unit_mass*unit_time); //convert to code unit
  } else {
    mdot  = -MDOTINIT;
  }
  for (i=0; i<nr; i++){
    l = i*ns;
    //    temp[i] = MU/R*(ADIABATICINDEX-1.0)*energ[l]/dens[l];
    temp[i] = BitschTemperature(mdot, Rmed[i]);
    Terr = 100;
    n = 0;
    while(Terr > epsilon){
      if (n < 50){
        Tnew = CalculateTNew(Rmed[i], temp[i], -mdot);
        Terr = fabs(temp[i] - Tnew)/temp[i];
        Told = temp[i];
        n += 1;
        temp[i] = Tnew;
      } else {
        temp[i] = 0.5*(Told+Tnew);
        break;
      }
    }
  }
  for (i=0; i<nr; i++){
    for (j=0; j<ns; j++){
      l = i*ns+j;
      dens[l] = -mdot/3./PI/(ALPHAVISCOSITY * sqrt(ADIABATICINDEX) * temp[i] * pow(Rmed[i], 1.5));
      if ((AccBoundary) && (DecInner)) {
        if (Rmed[i] < (GlobalRmed[0]+dsmoothin))
        dens[l] = (dens[l]-floordens) * exp(-pow(Rmed[i]-(GlobalRmed[0]+dsmoothin),2)/2./pow(dsmoothin,2)) + floordens;
      }
      energ[l] = temp[i] * dens[l] *R/MU/(ADIABATICINDEX-1);
    }
  }
  if ((AccBoundary) && (DecInner)) {
    for (j=0; j<ns; j++)
      dens[j] = dens[j+ns];
  }
  RefillSigma (Rho);
  RefillEnergy (Energy);
}

real CalculateTNew(r, temp, mdot)
  real r, temp, mdot;
{
  real phys_temp, phys_dens, sigma, H, rho3D, opac;
  real tau, taueff, T4, Lstar;
  sigma = mdot/3/M_PI/(ALPHAVISCOSITY*pow(r,1.5) *temp*sqrt(ADIABATICINDEX));
  Lstar = 4*M_PI*sigma_SB*(tstar*tstar*tstar*tstar)*(rstar*rstar);
/* Convert code temperature into Kelvins */
  phys_temp = temp * unit_temperature;
/* Convert 3D volume density into g.cm^-3 */
  H = sqrt(temp * r*r*r);  //H needs isothermal soundspeed
  rho3D = sigma/H/sqrt(2*PI);  // 3D density = sigma / 2H, in code units
  phys_dens = rho3D * unit_mass * pow(unit_length, -3);  // in kg.m^(-3)
  phys_dens *= 1e-3;  // in g.cm^(-3)
  opac = opLBL94 (phys_dens, phys_temp);
  opac *= (0.1 * pow(unit_length,-2) * pow(unit_mass,1)); //make opac dimensionless
  opac *= ZMETAL; //The effect of metalicity as Bitsch+2014
  tau = 0.5 *opac * sigma / sqrt(2*PI); 
  taueff = 3./8*tau + sqrt(3.)/4.+ 0.25/(tau+1e-20);
  T4 = taueff * mdot * pow(r,-3) * 3./8/PI/sigma_SB; //In code unit, viscous heating
  T4 += Lstar/28./M_PI/r/r*(H/r)/sigma_SB;  //Stellar irradiation
  temp = pow(T4,0.25);
  return temp;
}
