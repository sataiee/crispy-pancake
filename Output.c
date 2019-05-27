/** \file Output.c

Contains most of the functions that write the output files.  In
addition to the writing of hydrodynamics files (handled by SendOutput
()), this file also contains the functions that update the planet.dat
and bigplanet.dat files, and the functions that seek information about
the planets at a restart.
*/

#include "mp.h"

static real     Xplanet, Yplanet, VXplanet, VYplanet, MplanetVirtual;
extern real     LostMass;
extern boolean  Write_Density, Write_Velocity, Write_Energy, IsDisk;
extern boolean  Write_Temperature, Write_DivV, Write_TherHeat, Write_TherCool, Write_ViscHeat, ModifiedSoundSpeed, Write_RadDiff, Write_StarIrrad, Write_Opacity;
extern boolean  Write_Potential, Write_Test, Write_OneD_Fields, Write_pdv, Write_ArtVisc;
extern boolean  Write_gr, Write_gtheta, Write_OneD_Viscosity;
extern boolean  AdvecteLabel;

void EmptyPlanetSystemFile (sys)
     PlanetarySystem *sys;
{
  FILE *output;
  char name[256];
  int i, n;
  n = sys->nb;
  if (!CPU_Master) return;
  for (i = 0; i < n; i++) {
    sprintf (name, "%splanet%d.dat", OUTPUTDIR, i);
    output = fopen (name, "w");
    if (output == NULL) {
      fprintf (stderr, "Can't write %s file. Aborting.\n", name);
      prs_exit (1);
    }
    fclose (output);
  }
}

void WritePlanetFile (sys, timestep, n)
     int timestep;
     int n;
     PlanetarySystem *sys;     
{
  FILE *output;
  char name[256];
  real MenvelopeMax;
  if (!CPU_Master) return;
  printf ("Updating 'planet%d.dat'...", n);
  fflush (stdout);
  sprintf (name, "%splanet%d.dat", OUTPUTDIR, n);
  output = fopen (name, "a");
  if (output == NULL) {
    fprintf (stderr, "Can't write 'planet%d.dat' file. Aborting.\n", n);
    prs_exit (1);
  }
   MenvelopeMax = MenvRemoved[n]+MenvAccreted[n]+MenvRemained[n];
fprintf (output, "%d\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%d\t%#.18g\t%#.18g\t%#.18g\t%#.18g\n", timestep, Xplanet, Yplanet, VXplanet, VYplanet, MplanetVirtual, LostMass, PhysicalTime, OmegaFrame, mdcp, exces_mdcp, Menvelope[n],sys->acc[n], MdotEnvelope[n], MenvelopeMax, MenvCount[n],MenvRemained[n],sys->MassTaper[n], MenvRemoved[n], MenvAccreted[n]);
  fclose (output);
  printf ("done\n");
  fflush (stdout);
}

void WriteMassTrack (sys, timestep, Mdisk, Mevap, MaccPlanets)
     int timestep;
     real Mdisk, Mevap, MaccPlanets;  
     PlanetarySystem *sys;   
{
  FILE *output;
  char name[256];
  int i, j;
  real etot = 0.0, ekin=0.0, epot = 0.0, dx, dy;
  if (!CPU_Master) return;
  printf ("Updating 'masstrack.dat'...");
  fflush (stdout);
  sprintf (name, "%smasstrack.dat", OUTPUTDIR);
  output = fopen (name, "a");
  if (output == NULL) {
    fprintf (stderr, "Can't write 'masstrack.dat' file. Aborting.\n");
    prs_exit (1);
  }
  /* We also dump the total planets energy in this file since the other files contain *
   * individual information of planets */
	for (i = 0; i < sys->nb; i++){
		ekin += 0.5 * sys->mass[i] * (sys->vx[i]*sys->vx[i]+sys->vy[i]*sys->vy[i]);
		for (j = i+1; j < sys->nb; j++){
			dx = sys->x[i] - sys->x[j];
			dy = sys->y[i] - sys->y[j];
			epot -= G * sys->mass[i] * sys->mass[j] / sqrt(dx*dx + dy*dy);
		}
	}
	etot = epot + ekin;
  fprintf (output, "%d\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\n", \
                    timestep, PhysicalTime, Mdisk, Mevap, MaccPlanets, Rhole, etot);
  fclose (output);
  printf ("done\n");
  fflush (stdout);
}

void WritePlanetSystemFile (sys, t)
     PlanetarySystem *sys;
     int t;
{
  int i, n;
  n = sys->nb;
  for (i = 0; i < n; i++) {
    Xplanet = sys->x[i];
    Yplanet = sys->y[i];
    VXplanet = sys->vx[i];
    VYplanet = sys->vy[i];
    MplanetVirtual = sys->mass[i];
    WritePlanetFile (sys, t, i);
  }
}
   

void WriteBigPlanetFile (sys, timestep, n)
     int timestep;
     int n;
     PlanetarySystem *sys; 
{
  FILE *output;
  char name[256];
  real MenvelopeMax;
  if (!CPU_Master) return;
  sprintf (name, "%sbigplanet%d.dat", OUTPUTDIR, n);
  output = fopen (name, "a");
  if (output == NULL) {
    fprintf (stderr, "Can't write 'bigplanet.dat' file. Aborting.\n");
    prs_exit (1);
  }
  MenvelopeMax = MenvRemoved[n]+MenvAccreted[n]+MenvRemained[n];
fprintf (output, "%d\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%d\t%#.18g\t%#.18g\t%#.18g\t%#.18g\n", timestep, Xplanet, Yplanet, VXplanet, VYplanet, MplanetVirtual, LostMass, PhysicalTime, OmegaFrame, mdcp, exces_mdcp, Menvelope[n],sys->acc[n], MdotEnvelope[n], MenvelopeMax, MenvCount[n],MenvRemained[n],sys->MassTaper[n], MenvRemoved[n], MenvAccreted[n]);
  fclose (output);
}

void WriteBigPlanetSystemFile (sys, t)
     PlanetarySystem *sys;
     int t;
{
  int i, n;
  n = sys->nb;
  for (i = 0; i < n; i++) {
    Xplanet = sys->x[i];
    Yplanet = sys->y[i];
    VXplanet = sys->vx[i];
    VYplanet = sys->vy[i];
    MplanetVirtual = sys->mass[i];
    WriteBigPlanetFile (sys, t, i);
  }
}

real GetfromPlanetFile (timestep, column, n)
     int timestep, column, n;
{
  FILE *input;
  char name[256];
  char testline[829];
  int time;
  char *pt;
  double value;
  sprintf (name, "%splanet%d.dat", OUTPUTDIR, n);
  input = fopen (name, "r");
  if (input == NULL) {
    mastererr ("Can't read 'planet%d.dat' file. Aborting restart.\n",n);
    prs_exit (1);
  }
  if (column < 2) {
    mastererr ("Invalid column number in 'planet%d.dat'. Aborting restart.\n",n);
    prs_exit (1);
  }
  do {
    pt = fgets (testline, 828, input);
    sscanf (testline, "%d", &time);
  } while ((time != timestep) && (pt != NULL));
  if (pt == NULL) {
    mastererr ("Can't read entry %d in 'planet%d.dat' file. Aborting restart.\n", timestep,n);
    prs_exit (1);
  }
  fclose (input);
  pt = testline;
  while (column > 1) {
    pt += strspn(pt, "eE0123456789-.");
    pt += strspn(pt, "\t :=>_");
    column--;
  }
  sscanf (pt, "%lf", &value);
  return (real)value;
}

real GetfromTrackMassFile (timestep, column)
     int timestep, column;
{
  FILE *input;
  char name[256];
  char testline[256];
  int time;
  char *pt;
  double value;
  char *buftest;
    buftest=(char *)malloc(100*sizeof(char));
    free(buftest);
  sprintf (name, "%smasstrack.dat", OUTPUTDIR);
  input = fopen (name, "r");
  if (input == NULL) {
    printf ("Can't read '%smasstrack.dat' file. Aborting restart.\n",OUTPUTDIR);
    prs_exit (1);
  }
  if (column < 2) {
    mastererr ("Invalid column number in 'masstrack.dat'. Aborting restart.\n");
    prs_exit (1);
  }
  do {
    pt = fgets (testline, 255, input);
    sscanf (testline, "%d", &time);
  } while ((time != timestep) && (pt != NULL));
  if (pt == NULL) {
    mastererr ("Can't read entry %d in 'masstrack.dat' file. Aborting restart.\n", timestep);
    prs_exit (1);
  }
  fclose (input);
  pt = testline;
  while (column > 1) {
    pt += strspn(pt, "eE0123456789-.");
    pt += strspn(pt, "\t :=>_");
    column--;
  }
  sscanf (pt, "%lf", &value);
  return (real)value;
}


void RestartPlanetarySystem (timestep, sys)
     PlanetarySystem *sys;
     int timestep;
{
  int k, n, ip;
  real bufmass, bufvx, bufvy;
  real x,y,vx,vy,h,d,Ax,Ay,a,e,m;
  real mcore, menv, *Mswitch, rp;
  real *cs, axi[GLOBALNRAD], hplanet;
  char name[256];
  FILE *input;
  extern boolean FargoPlanete;
  n = sys->nb;
  cs = SoundSpeed->Field;
  mpi_make1Dprofile (cs, axi);
  Mswitch = (real *)malloc(n*sizeof(real));
  for (k = 0; k < sys->nb; k++) {
  	sprintf (name, "%splanet%d.dat", OUTPUTDIR, k);
		input = fopen (name, "r");
		if (input != NULL){
      fclose(input);
      sys->x[k]  = GetfromPlanetFile (timestep, 2, k);
      sys->y[k]  = GetfromPlanetFile (timestep, 3, k);
      sys->vx[k] = GetfromPlanetFile (timestep, 4, k);
      sys->vy[k] = GetfromPlanetFile (timestep, 5, k);
    }
    rp = sqrt((sys->x[k]*sys->x[k])+(sys->y[k]*sys->y[k]));
    ip = ReturnIndex(rp);
    hplanet = axi[ip]*sqrt(GlobalRmed[ip]);
    Mswitch[k] = MCRIFACTOR * (hplanet*hplanet*hplanet);
    if ((FargoPlanete) && (sys->mass[k] < Mswitch[k])) {
      sys->TorqueFlag[k] = YES;
    } else {
      sys->TorqueFlag[k] = NO;
    }
      /* Below we infer planet's semi-major axis and eccentricity */
      x  = sys->x[k];
      y  = sys->y[k];
      vx = sys->vx[k];
      vy = sys->vy[k];
      h = x*vy-y*vx;
      d = sqrt(x*x+y*y);
      m = 1.0+sys->mass[k];
      Ax = x*vy*vy-y*vx*vy-G*m*x/d;
      Ay = y*vx*vx-x*vx*vy-G*m*y/d;
      e = sqrt(Ax*Ax+Ay*Ay)/G/m;
      a = h*h/G/m/(1-e*e);
      sys->a[k]=a;
      sys->e[k]=e;
  }
  free(Mswitch);
}

void WriteDiskPolar(array, number)
     PolarGrid 	*array;
     int 	number;
{
  int           Nr, Ns;
  FILE          *dump;
  char 		name[256];
  real 		*ptr;
  ptr = array->Field;
  if (CPU_Master)
    sprintf (name, "%s%s%d.dat", OUTPUTDIR, array->Name, number);
  else
    sprintf (name, "%s%s%d.dat.%05d", OUTPUTDIR, array->Name, number, CPU_Rank);
  Nr = array->Nrad;
  Ns = array->Nsec;
  dump = fopen(name, "w");
  if (dump == NULL) {
    fprintf(stderr, "Unable to open '%s'.\n", name);
    prs_exit(1);
  }
  masterprint ("Writing '%s%d.dat'...", array->Name, number);
  fflush (stdout);
  MPI_Barrier (MPI_COMM_WORLD);
/* We strip the first CPUOVERLAP rings if the current CPU is not the 
   innermost one */
  if (CPU_Rank > 0) {
    ptr += CPUOVERLAP*Ns;
    Nr -=CPUOVERLAP ;
  }
/* We strip the last CPUOVERLAP rings if the current CPU is not the outermost
   one, equal to CPU_Highest in all cases */
  if (CPU_Rank != CPU_Highest) {
    Nr -=CPUOVERLAP;
  }
  fwrite (ptr, sizeof(real), Nr*Ns,dump);
  fclose(dump);
  fprintf(stdout, "%d/", CPU_Rank);  
  fflush(stdout);
  MPI_Barrier (MPI_COMM_WORLD);
  masterprint("done\n");
}

void Write1DFields(dens, gasvr, Temperature, number, sys)
     PolarGrid 	*dens, *gasvr, *Temperature;
     int 	number;
     PlanetarySystem *sys;
{
  FILE          *dump;
  char 		name[256];
  real 		*ptr, *OneDarray;
  int i, j, nr, ns, nplanet, n, jp[20];
  real thetp;
  extern boolean FargoPlanete;
  OneDarray = (real *)malloc(GLOBALNRAD*sizeof(real));
  if (FargoPlanete){
    nplanet = sys->nb;
    for (n=0; n<nplanet; n++){
      thetp = atan2(sys->y[n],sys->x[n]);
      if (thetp < 0) thetp += 2*PI;
      j = 0;
      while (j<NSEC && azimuth[j]<thetp) j++;
      if (j == NSEC) j=0;
      jp[n] = j;
    }
  }
  // 1D density -> ascii file ~/gasdens1Dxx.dat
  ptr = dens->Field;
  mpi_make1Dprofile (ptr, axidens);
  if (CPU_Master) {
    sprintf (name, "%s%s1D%d.dat", OUTPUTDIR, dens->Name, number);
    dump = fopen(name, "w");
    if (dump == NULL) {
      fprintf(stderr, "Unable to open '%s'.\n", name);
      prs_exit(1);
    }
    for ( i = 0; i < GLOBALNRAD; i++ ) {
      fprintf (dump,"%-.18e %-.18e\n",GlobalRmed[i],axidens[i]);
    }
    fclose (dump);
  }
   // for using fargo with planete, write the density at the azimuth of each plant    
  if (FargoPlanete){
    for (n=0; n<nplanet; n++){
      mpi_Find1Dprofile(ptr,jp[n],OneDarray);
      if (CPU_Master) {
        sprintf (name, "%s%s%dnp%d.dat", OUTPUTDIR, dens->Name, number, n);
        dump = fopen(name, "w");
        if (dump == NULL) {
          fprintf(stderr, "Unable to open '%s'.\n", name);
           prs_exit(1);
        }
        for ( i = 0; i < GLOBALNRAD; i++ ) {
          fprintf (dump,"%-.18e %-.18e\n",GlobalRmed[i],OneDarray[i]);
        }
        fclose (dump);
      }
    }
  }
  // 1D temperature -> ascii file ~/gasTemperature1Dxx.dat
  ptr = Temperature->Field;
  mpi_make1Dprofile (ptr, axitemp);
  if (CPU_Master) {    
    sprintf (name, "%s%s1D%d.dat", OUTPUTDIR, Temperature->Name, number);
    dump = fopen(name, "w");
    if (dump == NULL) {
      fprintf(stderr, "Unable to open '%s'.\n", name);
      prs_exit(1);
    }
    for ( i = 0; i < GLOBALNRAD; i++ ) {
      fprintf (dump,"%.18g\t%.18g\n",GlobalRmed[i],axitemp[i]);
    }
    fclose (dump);
  }
  if (FargoPlanete){
    for (n=0; n<nplanet; n++){
      mpi_Find1Dprofile(ptr,jp[n],OneDarray);
      if (CPU_Master) {
        sprintf (name, "%s%s%dnp%d.dat", OUTPUTDIR, Temperature->Name, number, n);
        dump = fopen(name, "w");
        if (dump == NULL) {
          fprintf(stderr, "Unable to open '%s'.\n", name);
           prs_exit(1);
        }
        for ( i = 0; i < GLOBALNRAD; i++ ) {
          fprintf (dump,"%-.18e %-.18e\n",GlobalRmed[i],OneDarray[i]);
        }
        fclose (dump);
      }
    }
  }
  // 1D vrad -> ascii file ~/gasvrad1Dxx.dat
  ptr = gasvr->Field;
  mpi_make1Dprofile (ptr, OneDarray);
  if (CPU_Master) {    
    sprintf (name, "%s%s1D%d.dat", OUTPUTDIR, gasvr->Name, number);
    dump = fopen(name, "w");
    if (dump == NULL) {
      fprintf(stderr, "Unable to open '%s'.\n", name);
      prs_exit(1);
    }
    for ( i = 0; i < GLOBALNRAD; i++ ) {
      fprintf (dump,"%.18g\t%.18g\n",GlobalRmed[i],OneDarray[i]);
    }
    fclose (dump);
  }
  free(OneDarray);
}


void Write1DViscosity(Rho, number)
     int 	number;
     PolarGrid 	*Rho;
     
{
  FILE          *dump;
  char 		name[256];
  real 		viscosity;
  int i;
  // 1D viscosity -> ascii file ~/viscosity1DXX.dat
  if (CPU_Master) {
    sprintf (name, "%s%s1D%d.dat", OUTPUTDIR, "Viscosity", number);
    dump = fopen(name, "w");
    if (dump == NULL) {
      fprintf(stderr, "Unable to open '%s'.\n", name);
      prs_exit(1);
    }
    mpi_make1Dprofile (Rho->Field, axidens);
    mpi_make1Dprofile (SoundSpeed->Field, axics);
    for ( i = 0; i < GLOBALNRAD; i++ ) {
      viscosity = FViscosity (GlobalRmed[i], axics[i], axidens[i]);
      fprintf (dump,"%.18e %.18e\n",GlobalRmed[i],viscosity);
    }
    fclose (dump);
  }
}

void WriteDim () {	  
  char filename[256];
  FILE 	*dim;
  if (!CPU_Master) return;
  sprintf (filename, "%sdims.dat", OUTPUTDIR);
  if ((dim = fopen (filename, "w")) == NULL) {
    fprintf (stderr, "Unable to open %s. Program stopped\n", filename);
    prs_exit (1);
  }
  fprintf (dim,"%d\t%d\t\t%d\t%d\t%f\t%d\t%d\t%d\n",0,0,0,0,RMAX, NTOT/NINTERM, GLOBALNRAD, NSEC);
  fclose (dim);
}

void SendOutput (index, dens, gasvr, gasvt, gasenerg, label, sys)
     int          index;
     PolarGrid   *dens, *gasvr, *gasvt, *label, *gasenerg;
     PlanetarySystem *sys;
{
  if (CPU_Master)
    printf ("\n*** OUTPUT %d ***\n", index);
  if (IsDisk == YES) {
    if (AdvecteLabel == YES) WriteDiskPolar (label, index);
    if (Write_Density == YES) WriteDiskPolar (dens, index);
    if (Write_Velocity == YES) {
      WriteDiskPolar (gasvr, index);
      WriteDiskPolar (gasvt, index);
    }
    //if (AdvecteLabel == YES) WriteDiskPolar (label, index);
    if (Write_Energy == YES) WriteDiskPolar (gasenerg, index);
    if (Write_Temperature == YES) {
      if (ModifiedSoundSpeed) 
	WriteDiskPolar (SoundSpeed, index);
      else
	WriteDiskPolar (Temperature, index);
    }
    if (Write_DivV == YES) WriteDiskPolar (DivergenceVelocity, index);
    if (Write_ViscHeat == YES)  WriteDiskPolar (ViscHeat, index);
    if (Write_TherHeat == YES)  WriteDiskPolar (ThermHeat, index);
    if (Write_TherCool == YES)  WriteDiskPolar (ThermCool, index);
    if (Write_RadDiff == YES)  WriteDiskPolar (RadDiffusion, index);
    if (Write_StarIrrad == YES) WriteDiskPolar (StarIrradiation, index);
    if (Write_pdv == YES) WriteDiskPolar (pdvEnergy, index);
    if (Write_ArtVisc == YES) WriteDiskPolar (ArtViscHeat, index);
    if (Write_Opacity == YES)   WriteDiskPolar (Opacity, index);
    if (Write_Potential == YES)  WriteDiskPolar (Potential, index);
    if (Write_Test == YES)  WriteDiskPolar (Test, index);
    if (Write_gr == YES)  WriteDiskPolar (gr, index);
    if (Write_gtheta == YES)  WriteDiskPolar (gtheta, index);
    if (Write_OneD_Fields == YES) Write1DFields (dens, gasvr, Temperature, index, sys);
    if (Write_OneD_Viscosity == YES) Write1DViscosity(dens, index);
    MPI_Barrier (MPI_COMM_WORLD);
    if (Merge && (CPU_Number > 1)) merge (index);
  }
}
