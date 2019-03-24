/** \file Psys.c

Contains the functions that set up the planetary system configuration.
In addition, the last two functions allow to track the first planet
(number 0) of the planetary system, in order to perform a calculation
in the frame corotating either with this planet or with its
guiding-center.

*/

#include "mp.h"

static real Xplanet, Yplanet;
extern boolean GuidingCenter, BinaryCenter, RetrogradeBinary;
extern boolean OneDRun, FargoPlanete;

int FindNumberOfPlanets (filename)
     char *filename;
{
  FILE *input;
  char s[512];
  int Counter=0;
  
  if (FargoPlanete) {
    /* NEW (Feb. 2014): Fargo/Planete coupling: we get the planets
       number from file 'donnees_MP' written by Planete in the output
       directory. The number of planets is the number of lines in that
       file - 1 */
    sprintf (filename, "donnees_MP");
  }
  input = fopen (filename, "r");
  if (input == NULL) {
    fprintf (stderr, "Error : can't find '%s'.\n", filename);
    prs_exit (1);
  }
  while (fgets(s, 510, input) != NULL) {
    if (isalnum(s[0])) // check if first character alphabetical or numerical
      Counter++;
  }
  fclose (input);
  if (FargoPlanete) Counter--;
  return Counter;
}

PlanetarySystem *AllocPlanetSystem (nb)
int nb;
{
  real *mass, *x, *y, *vx, *vy, *acc, *a, *e;
  boolean *feeldisk, *feelothers, *binary, *torqueflag;
  int i, j;
  PlanetarySystem *sys;
  sys  = (PlanetarySystem *)malloc (sizeof(PlanetarySystem));
  if (sys == NULL) {
    fprintf (stderr, "Not enough memory.\n");
    prs_exit (1);
  }
  x    = (real *)malloc (sizeof(real)*(nb+1));
  y    = (real *)malloc (sizeof(real)*(nb+1));
  vy   = (real *)malloc (sizeof(real)*(nb+1));
  vx   = (real *)malloc (sizeof(real)*(nb+1));
  a    = (real *)malloc (sizeof(real)*(nb+1));
  e    = (real *)malloc (sizeof(real)*(nb+1));
  mass = (real *)malloc (sizeof(real)*(nb+1));
  acc  = (real *)malloc (sizeof(real)*(nb+1));
  if ((x == NULL) || (y == NULL) || (vx == NULL) || (vy == NULL) || (acc == NULL) || \
      (mass == NULL) || (a == NULL) || (e == NULL)) {
    fprintf (stderr, "Not enough memory.\n");
    prs_exit (1);
  }
  feeldisk   = (boolean *)malloc (sizeof(real)*(nb+1));
  feelothers = (boolean *)malloc (sizeof(real)*(nb+1));
  binary     = (boolean *)malloc (sizeof(real)*(nb+1));
  torqueflag = (boolean *)malloc (sizeof(real)*(nb+1));
  if ((feeldisk == NULL) || (feelothers == NULL) || (binary == NULL) || (torqueflag == NULL)) {
    fprintf (stderr, "Not enough memory.\n");
    prs_exit (1);
  }
  sys->x = x;
  sys->y = y;
  sys->vx= vx;
  sys->vy= vy;
  sys->a= a;
  sys->e= e;
  sys->acc=acc;
  sys->mass = mass;
  sys->FeelDisk = feeldisk;
  sys->FeelOthers = feelothers;
  sys->Binary = binary;
  sys->TorqueFlag = torqueflag;
  for (i = 0; i < nb; i++) {
    x[i] = y[i] = vx[i] = vy[i] = a[i] = e[i] = mass[i] = acc[i] = 0.0;
    feeldisk[i] = feelothers[i] = YES;
    binary[i] = torqueflag[i] = NO;
  }
  return sys;
}


void FreePlanetary (sys)
     PlanetarySystem *sys;
{
  free (sys->x);
  free (sys->vx);
  free (sys->y);
  free (sys->vy);
  free (sys->a);
  free (sys->e);
  free (sys->mass);
  free (sys->acc);
  free (sys->FeelOthers);
  free (sys->FeelDisk);
  free (sys->Binary);
  free (sys->TorqueFlag);
  free (sys);
}

PlanetarySystem *InitPlanetarySystem (filename,NbRestart, Rho, force)
char *filename;
int NbRestart;
PolarGrid *Rho;
Force *force;
{
  extern boolean CICPlanet;
  FILE *input, *DIM, *tempinput;
  char name_dim[1024];
  char name[256];
  char s[512], nm[512], test1[512], test2[512], test3[512], *s1;
  PlanetarySystem *sys;
  int i=0, j, nb, counter, nbplanets, fixedpls, Foo, ii, nr, ns, l;
  float mass, dist, accret, eccentricity, azimuth, vr, vtheta;
  real MassInside=0, *cs;
  Pair gamma;
  extern real Runtime;
  int OldNSEC;
  real m0, m1, mbin, *Mswitch;
  real a, e, incl, foo, mcore, menv, dcore, stime, deltamenv;
  boolean feeldis, feelothers, binary;
  nb = FindNumberOfPlanets (filename);
  cs = SoundSpeed->Field;
  nr = SoundSpeed->Nrad;
  ns = SoundSpeed->Nsec;
  Mswitch = (real *)malloc(nb*sizeof(real));
  if (CPU_Master)
    printf ("%d planet(s) found.\n", nb);
  sys = AllocPlanetSystem (nb);
  if (FargoPlanete) {
    /* NEW (Feb. 2014): Fargo/Planete coupling: we get the planets
       initial orbital parameters from file 'planete_to_fargo.dat' written by
       Planete in the output directory. */
       sprintf (filename, "%splanete_to_fargo.dat", OUTPUTDIR);
  }
  input = fopen (filename, "r");
  if (input == NULL) {
    fprintf (stderr, "Error : can't find '%s'. Aborting.\n", filename);
    prs_exit (1);
  }
  sys->nb = nb;
  if (!FargoPlanete) {
    for (i = 0; i < nr; i++){
      for (j = 0; j < ns; j++){
        l = i*ns + j;
        cs[l] = AspectRatio(Rmed[i])*sqrt(G*1.0/Rmed[i])*pow(Rmed[i], FLARINGINDEX);
      }
    }
    i = 0;
    j = 0;
    while (fgets(s, 510, input) != NULL) {
      sscanf(s, "%s ", nm);
      if (isalpha(s[0])) {
       s1 = s + strlen(nm);
       sscanf(s1 + strspn(s1, "\t :=>_"), "%f %f %f %s %s %s %f %f", &dist, &mass, &accret, test1, test2, test3, \
                                                                     &eccentricity, &azimuth);
       if ( CICPlanet ) {
         /* initialization puts planet at the interface between two
            cells (with eccentricity = 0 only) */
         j = 0;
         while ( GlobalRmed[j] < dist ) j++;
         dist = Radii[j+1];
       }
       if (NbRestart != 0) {      
        sprintf (name, "%splanet%d.dat", OUTPUTDIR, i);
        tempinput = fopen (name, "r");
        if (tempinput == NULL) {
          masterprint ("Can't read 'planet%d.dat' file. Using configuration information.\n",i);
          PlanetMassAtRestart[i] = (real)mass;
          Menvelope[i] = 0.0;
          MenvRemained[i] = 0.0;
        } else {
          PlanetMassAtRestart[i] = GetfromPlanetFile (NbRestart, 6, i);
          Menvelope[i] = GetfromPlanetFile (NbRestart, 12, i);
          MenvRemained[i] = GetfromPlanetFile (NbRestart, 17, i);
          fclose(tempinput);
        }
       } else { 
         PlanetMassAtRestart[i] = (real)mass;
         Menvelope[i] = 0.0;
         MenvRemained[i] = 0.0;
       }
       MenvRemoved[i] = 0.0;
       MenvAccreted[i] = 0.0;
       if ((MASSTAPER > 1e-3) && (NbRestart == 0))
          PlanetMassAtRestart[i] = 0.0;
       sys->mass[i] =  PlanetMassAtRestart[i];
       /* Mass of planet i at the end of the calculation */
       FinalPlanetMass[i] = (real)mass; 
       feeldis = feelothers = YES;
       binary = NO;
       if (tolower(*test1) == 'n') feeldis = NO;
       if (tolower(*test2) == 'n') feelothers = NO;
       if (tolower(*test3) == 'y') binary = YES;
       MassInside += PlanetMassAtRestart[i];
       sys->x[i] = (real)dist*(1.0+eccentricity)*cos(azimuth);
       sys->y[i] = (real)dist*(1.0+eccentricity)*sin(azimuth);
       sys->a[i] = (real)dist;
       sys->e[i] = eccentricity;
       sys->acc[i] = accret;
       sys->FeelDisk[i] = feeldis;
       sys->FeelOthers[i] = feelothers;
       sys->Binary[i] = binary;
       sys->TorqueFlag[i] = NO;
       if ((feeldis == YES) && (tolower(*CORRECTVPLANET) == 'y')) {
         if (eccentricity > 0.0)
           mastererr("When you are using the correction for the planetary mass,\
                      you cannot have non-zero eccentricity.\n");
         CorrectVgasSG = YES;
         gamma = ComputeAccel (force, Rho, sys->x[i], sys->y[i], PlanetMassAtRestart[i], sys, 2);
         gamma.x -= (1.0+MassInside)/dist/dist;
         vtheta = sqrt(dist*fabs(gamma.x)); // This condition is only for circular planets
       } else {
         vtheta = sqrt((1.0+MassInside)/dist)*sqrt( (1.0-eccentricity)/(1.0+eccentricity) );
       }
       vr = 0.0;
       sys->vy[i] = vr*sin(azimuth) + vtheta*cos(azimuth);
       sys->vx[i] = vr*cos(azimuth) + vtheta*sin(azimuth);
       i++;
      }
    }
  } else {
    /* NEW (Feb. 2014): Fargo/Planete coupling: we get the planets
    initial orbital parameters from file 'donnees_MP' written by
    Planete in the output directory. */
    if (!OneDRun){
      sprintf (name_dim, "%sdims.dat", OUTPUTDIR);
      DIM = fopen (name_dim, "r");
      if (DIM == NULL) 
        masterprint("dims.dat cannot be read in Psys.c\n");
      fscanf (DIM,"%d %d %d %d %lg %d %d %d\n",&Foo,&Foo,&Foo,&Foo,&foo,&Foo,&Foo,&OldNSEC);
      fclose (DIM);
    }
    while (fgets(s, 510, input) != NULL) {
      sscanf(s, "%s ", nm);
      if (isalnum(s[0])) {
        sscanf(s + strspn(s, "\t :=>_"), "%lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %d", &a, &e, \
              &incl, &foo, &foo, &foo, &mcore, &menv, &dcore, &stime, &deltamenv, &fixedpls);
        /* -------------------------- */
        /* Mass of planet i at the end of the calculation */
        FinalPlanetMass[i] = mcore+menv - deltamenv; //deltamcore will be accreted
        if (NbRestart != 0) {
          PlanetMassAtRestart[i] = GetfromPlanetFile (NbRestart, 6, i);
          Menvelope[i] = GetfromPlanetFile (NbRestart, 12, i);
          if (OldNSEC == 1)
            Menvelope[i] = menv - deltamenv;
          MenvRemained[i] = GetfromPlanetFile (NbRestart, 17, i);
          PlanetMassAtRestart[i] -= Menvelope[i];
          if (a <= WKZRMIN)
           FinalPlanetMass[i]=0.0;
          MenvCount[i] = GetfromPlanetFile (NbRestart, 16, i); //Number of timesteps that the planet cound not accrete as much as it should
        } else {
          PlanetMassAtRestart[i] = 0.0;
          Menvelope[i] = 0.0;
          MenvCount[i] = 0;
          MenvRemained[i] = 0.0;
        }
        MenvAccreted[i] = 0.0;
        MenvRemoved[i] = 0.0;
        sys->mass[i] = PlanetMassAtRestart[i];
        MdotEnvelope[i] = deltamenv/Runtime;
        /* This is done only the first time Fargo is called by Planete */
        sys->x[i] = a*(1.0+e);
        sys->y[i] = 0.0;
        sys->a[i] = a;
        sys->e[i] = e;
        sys->vy[i] = (real)sqrt(G*(1.0+PlanetMassAtRestart[i])/a)*       \
        sqrt( (1.0-e)/(1.0+e) );
        sys->vx[i] = -0.0000000001*sys->vy[i];
        sys->acc[i] = 0.0; // since mass accretion is taken care of by Planete
        if (fixedpls == 1) 
          sys->FeelDisk[i] = sys->FeelOthers[i] = NO;
        else 
          sys->FeelDisk[i] = sys->FeelOthers[i] = YES;
        sys->TorqueFlag[i] = NO;
        sys->Binary[i] = NO;  // not implemented...
        Mswitch[i] = MCRIFACTOR * pow(AspectRatio(a)*pow(a, FLARINGINDEX),3); 
        if (sys->mass[i] < Mswitch[i])
          sys->TorqueFlag[i] = YES;
        if (FinalPlanetMass[i] == 0.0) {
          sys->FeelDisk[i] = NO;
          sys->FeelOthers[i] = NO;
        }
        i++;
      }
    }
  }
  fclose(input);
  HillRadius = sys->x[0] * pow( sys->mass[0]/3., 1./3. );
  free(Mswitch);
  return sys; 
}

void ListPlanets (sys)
     PlanetarySystem *sys;
{
  int nb;
  int i;
  nb = sys->nb;
  if (!CPU_Master) return;
  for (i = 0; i < nb; i++) {
    printf ("Planet number %d\n", i);
    printf ("---------------\n");
    printf ("mass = %lg\n", sys->mass[i]);
    printf ("x = %.10f\ty = %.10f\n", sys->x[i],sys->y[i]);
    printf ("vx = %.10f\tvy = %.10f\n", sys->vx[i],sys->vy[i]);
    printf ("a = %.10f\te = %.10f\n", sys->a[i],sys->e[i]);
    if (sys->acc[i] != 0.0) {
              printf ("accretion time = %e\n", 1.0/(sys->acc[i]));
       } else {
              printf("Non-accreting.\n");
       }
    if (sys->FeelDisk[i] == YES) {
      printf ("Feels the disk potential\n");
    } else {
      printf ("Doesn't feel the disk potential\n");
    }
    if (sys->FeelOthers[i] == YES) {
      printf ("Feels the other planets potential\n");
    } else {
      printf ("Doesn't feel the other planets potential\n");
    }
    if (sys->Binary[i] == YES) {
      printf ("Is in a binary system\n");
    } else {
      printf ("Is not in a binary system\n");
    }
    printf ("\n");
  }
}

real GetPsysInfo (sys, action)
     PlanetarySystem *sys;
     boolean action;
{
  extern boolean SelfGravity, SGZeroMode;
  extern boolean CorotateWithOuterPlanet;
  real d1, d2, cross;
  real x,y, vx, vy, m, h, d, Ax, Ay, e, a, E, M;
  real xc, yc, vxc, vyc, omega;
  real x0, x1, y0, y1, vx0, vx1, vy0, vy1, m0, m1;
  real arg, PerihelionPA;
  real ri, rip1, dr, sgacc;
  int ipl;
  if (!CorotateWithOuterPlanet) {
    xc = x = sys->x[0];
    yc = y = sys->y[0];
    vxc = vx= sys->vx[0];
    vyc = vy= sys->vy[0];
    m = sys->mass[0]+1.;
  } else {
    xc = x = sys->x[1];
    yc = y = sys->y[1];
    vxc = vx= sys->vx[1];
    vyc = vy= sys->vy[1];
    m = sys->mass[1]+1.;
  }
  if (BinaryCenter) {
    x0 = sys->x[0];
    x1 = sys->x[1];
    y0 = sys->y[0];
    y1 = sys->y[1];
    vx0 = sys->vx[0];
    vx1 = sys->vx[1];
    vy0 = sys->vy[0];
    vy1 = sys->vy[1];
    m0 = sys->mass[0];
    m1 = sys->mass[1];
    if (m0 == 0.0)
      m0 = FinalPlanetMass[0];
    if (m1 == 0.0)
      m1 = FinalPlanetMass[1];
    xc = x = (m0*x0 + m1*x1) / (m0+m1);
    yc = y = (m0*y0 + m1*y1) / (m0+m1);
    vxc = vx = (m0*vx0 + m1*vx1) / (m0+m1);
    vyc = vy = (m0*vy0 + m1*vy1) / (m0+m1);
    m = m0+m1+1.;
  }
  h = x*vy-y*vx;
  d = sqrt(x*x+y*y);
  Ax = x*vy*vy-y*vx*vy-G*m*x/d;
  Ay = y*vx*vx-x*vx*vy-G*m*y/d;
  e = sqrt(Ax*Ax+Ay*Ay)/m;
  a = h*h/G/m/(1.-e*e);
  if (e == 0.0) {
    arg = 1.0;
  } else {
    arg = (1.0-d/a)/e;
  }
  if (fabs(arg) >= 1.0) 
    E = PI*(1.-arg/fabs(arg))/2.;
  else
    E = acos((1.0-d/a)/e);
  if ((x*y*(vy*vy-vx*vx)+vx*vy*(x*x-y*y)) < 0) E= -E;
  M = E-e*sin(E);
  omega = sqrt(m/a/a/a);
  /* Here omega is modified to include self-gravity */
  if ( SelfGravity ) {
    if ( !SGZeroMode )
      mpi_make1Dprofile (SG_Accr, GLOBAL_AxiSGAccr);
    else
      GLOBAL_AxiSGAccr = SG_Accr;
    ipl = 0;
    while (GlobalRmed[ipl] <= a) ipl++;
    ri = GlobalRmed[ipl];
    rip1 = GlobalRmed[ipl+1];
    dr = rip1 - ri;
    sgacc = (a - ri)*GLOBAL_AxiSGAccr[ipl+1] + (rip1 - a)*GLOBAL_AxiSGAccr[ipl];
    sgacc /= dr;
    omega *= (real)sqrt(1. - a*a*sgacc/m);
  }
  PerihelionPA=atan2(Ay,Ax);
  if (GuidingCenter == YES) {
    xc = a*cos(M+PerihelionPA);
    yc = a*sin(M+PerihelionPA);
    vxc = -a*omega*sin(M+PerihelionPA);
    vyc =  a*omega*cos(M+PerihelionPA);
  } 
  if (e < 1e-8) {
    xc = x;
    yc = y;
    vxc = vx;
    vyc = vy;
  }
  switch (action) {
  case MARK: 
    Xplanet = xc;
    Yplanet = yc;
    return 0.;
    break;
  case GET:
    x = xc;
    y = yc;
    vx = vxc;
    vy = vyc;
    d2 = sqrt(x*x+y*y);
    d1 = sqrt(Xplanet*Xplanet+Yplanet*Yplanet);
    cross = Xplanet*y-x*Yplanet;
    Xplanet = x;
    Yplanet = y;
    return asin(cross/(d1*d2));
    break;
  case FREQUENCY:
    return omega;
    break;
  }
  return 0.0;
}

void RotatePsys (sys, angle)       /* Rotate by angle '-angle' */
     PlanetarySystem *sys;
     real angle;
{
  int nb;
  int i;
  real sint, cost, xt, yt;
  nb = sys->nb;
  sint = sin(angle);
  cost = cos(angle);
  for (i = 0; i < nb; i++) {
    xt = sys->x[i];
    yt = sys->y[i];
    sys->x[i] = xt*cost+yt*sint;
    sys->y[i] = -xt*sint+yt*cost;
    xt = sys->vx[i];
    yt = sys->vy[i];
    sys->vx[i] = xt*cost+yt*sint;
    sys->vy[i] = -xt*sint+yt*cost;
  }
}
