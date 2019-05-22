/** \file Interpret.c

Contains the functions required to read the parameter file, and
functions that provide runtime information. The function var()
associates a string to a global variable. The function ReadVariables()
reads the content of a parameter file.  In addition, this file
contains a function that prints the command line usage to the standard
output, a function that provides verbose information about the setup
(if the -v switch is set on the command line), and functions that act
as a chronometer (if the -t switch is set on the command line).
*/

#include "mp.h"
#define MAXVARIABLES 500

extern int      begin_i;
extern boolean  OpenInner, OneDRun;
static Param    VariableSet[MAXVARIABLES];
static int      VariableIndex = 0;
static int      FirstStep = YES;
static clock_t  First, Preceeding, Current, FirstUser, CurrentUser, PreceedingUser;
static long     Ticks;
boolean         FastTransport = YES, GuidingCenter = NO, BinaryCenter = NO, Indirect_Term = YES, RetrogradeBinary = NO;
boolean         IsDisk = YES, NonReflecting = NO, Corotating = NO, OuterSourceMass = NO, Evanescent = NO, MixedBC = NO;
boolean         Write_Density = YES, Write_Velocity = YES, Write_Energy = NO, Write_pdv=NO, Write_ArtVisc=NO;
boolean         Write_Temperature = NO, Write_DivV = NO, Write_OneD_Fields = NO;
boolean         Write_TherHeat = NO, Write_TherCool = NO, Write_ViscHeat = NO, Write_RadDiff=NO, Write_StarIrrad = NO, Write_Opacity=NO;
boolean         Write_Potential = NO, Write_Test = NO, Write_gr = NO, Write_gtheta = NO;
boolean         SelfGravity = NO, SGZeroMode = NO, ZMPlus = NO, AddNoise = NO;
boolean         EnergyEquation = NO, ModifiedSoundSpeed = NO, ThermalDiffusion = NO, ThermalCooling = NO, ViscousHeating = YES, TempPresc = NO, BetaCooling = NO;
boolean         CICPlanet = NO, ForcedCircular = NO, ForcedInnerCircular = NO, ForcedCircularTemp = NO;;
boolean         MHDLSA = NO, HighMCutoff = NO;
boolean         KNOpen = NO;
boolean         AdvecteLabel = NO;
boolean         SoftWriting = NO;
boolean         RetrogradePlanet = NO, AddMass = NO, ImposedDensity = NO, ImposedAlpha = NO;
boolean         ReadPlanetFileAtRestart = YES, DontApplySubKeplerian = NO;
boolean         CustTqExc = NO;
boolean         DampToIni = YES, DampToAxi = NO, CutForces = NO;
boolean         CorotateWithOuterPlanet = NO;
boolean         DiscEvaporation = NO;
boolean         CustomizedIT = NO;
boolean         FargoPlanete = NO, ExponentialDecay = NO;
boolean         PhotoEvapor = NO, AccBoundary = NO, Write_Sigdot, MdotHartmann = NO, DecInner = NO, OpInner = NO, StellarIrradiation=NO, Alexboundary=NO;
boolean         InitEquilibrium = NO, Write_OneD_Viscosity = NO, SmoothAtPlanet = NO, WriteForces = NO;

void var(name, ptr, type, necessary, deflt)
     char           *name;
     char           *ptr;
     int             type;
     int             necessary;
     char           *deflt;
{
  real            valuer;
  int             valuei;
  double          temp;
  sscanf (deflt, "%lf", &temp);
  valuer = (real) (temp);
  valuei = (int) valuer;
  strcpy(VariableSet[VariableIndex].name, name);
  VariableSet[VariableIndex].variable = ptr;
  VariableSet[VariableIndex].type = type;
  VariableSet[VariableIndex].necessary = necessary;
  VariableSet[VariableIndex].read = NO;
  if (necessary == NO) {
    if (type == INT) {
      *((int *) ptr) = valuei;
    } else if (type == REAL) {
      *((real *) ptr) = valuer;
    } else if (type == STRING) {
      strcpy (ptr, deflt);
    }
  }
  VariableIndex++;
}

void ReadVariables(filename)
     char *filename;
{
  char            nm[300], s[350],stringval[290];
  char           *s1;
  double          temp;
  real            valuer;
  int             i, found, valuei, success, type;
  int            *ptri;
  real           *ptrr;
  FILE           *input;
  struct  stat   stoutput;
  char           CommandLine[1024];
 
  InitVariables();
  input = fopen(filename, "r");
  if (input == NULL) {
    mastererr ("Unable to read '%s'. Program stopped.\n",filename);
    prs_exit(1);
  }
  mastererr ("Reading parameters file '%s'.\n", filename);
  while (fgets(s, 349, input) != NULL) {
    success = sscanf(s, "%s ", nm);
    if ((nm[0] != '#') && (success == 1)) {  /* # begins a comment line */
      s1 = s + strlen(nm);
      sscanf(s1 + strspn(s1, "\t :=>_"), "%lf", &temp);
      sscanf(s1 + strspn(s1, "\t :=>_"), "%289s ", stringval);
      valuer = (real) temp;
      valuei = (int) temp;
      for (i = 0; i < strlen(nm); i++) {
        nm[i] = (char) toupper(nm[i]);
      }
      found = NO;
      for (i = 0; i < VariableIndex; i++) {
        if (strcmp(nm, VariableSet[i].name) == 0) {
          if (VariableSet[i].read == YES)
            mastererr("Warning : %s defined more than once.\n", nm);
          found = YES;
          VariableSet[i].read = YES;
          ptri = (int *) (VariableSet[i].variable);
          ptrr = (real *) (VariableSet[i].variable);
          if (VariableSet[i].type == INT) {
            *ptri = valuei;
          } else if (VariableSet[i].type == REAL) {
            *ptrr = valuer;
          } else if (VariableSet[i].type == STRING) {
            strcpy (VariableSet[i].variable, stringval);
          }
        }
      }
      if (found == NO)
        mastererr("Warning : variable %s defined but non-existent in code.\n", nm);
    }
  }
  /* NEW: project fargo/planet: if -w flag  activated, run with 1D disc */
  if (OneDRun == YES) {
    NSEC = 1;
    Indirect_Term = NO;
    masterprint("Flag -w is active, meaning the disc will be considered as\n");
    masterprint("being 1D, even if NSEC > 1 in parameter file. Check this is OK!\n");
  }
  /* ----------------- */
  found = NO;
  for (i = 0; i < VariableIndex; i++) {
    if ((VariableSet[i].read == NO) && (VariableSet[i].necessary == YES)) {
      if (found == NO) {
        mastererr("Fatal error : undefined mandatory variable(s):\n");
        found = YES;
      }
      mastererr("%s\n", VariableSet[i].name);
    }
    if (found == YES)
    prs_exit(1);
  }
  found = NO;
  for (i = 0; i < VariableIndex; i++) {
    if (VariableSet[i].read == NO) {
      if (found == NO) {
        mastererr("Secondary variables omitted :\n");
        found = YES;
      }
      if ((type = VariableSet[i].type) == REAL)
        mastererr("%s ;\t Default Value : %.5g\n", VariableSet[i].name, *((real *) VariableSet[i].variable));
      if (type == INT)
        mastererr("%s ;\t Default Value : %d\n", VariableSet[i].name, *((int *) VariableSet[i].variable));
      if (type == STRING)
        mastererr("%s ;\t Default Value : %s\n", VariableSet[i].name, VariableSet[i].variable);
    }
  }
  if ((*ADVLABEL == 'y') || (*ADVLABEL == 'Y')) AdvecteLabel = YES;
  if ((*OUTERSOURCEMASS == 'y') || (*OUTERSOURCEMASS == 'Y')) OuterSourceMass = YES;
  if ((*TRANSPORT == 's') || (*TRANSPORT == 'S')) FastTransport = NO;
  if ((*OPENINNERBOUNDARY == 'O') || (*OPENINNERBOUNDARY == 'o')) OpenInner = YES;
  if ((*OPENINNERBOUNDARY == 'N') || (*OPENINNERBOUNDARY == 'n')) NonReflecting = YES;
  if ((*OPENINNERBOUNDARY == 'E') || (*OPENINNERBOUNDARY == 'e')){
    Evanescent = YES;
    if ((*ZEROINNER == 'O') || (*ZEROINNER == 'o')){
      OpInner = YES;
      masterprint("Your inner boundary will be open.\n");
    }
    if ((*ZEROINNER == 'D') || (*ZEROINNER == 'd')){
      DecInner = YES;
      masterprint("You set an innerd disc edge.\n");
     if (DISKINNEREDGE == 0.0)
      mastererr("You forgot to set the radius of the inner edge.\n");
    }
  }
  if ((*OPENINNERBOUNDARY == 'M') || (*OPENINNERBOUNDARY == 'm')) MixedBC = YES;
  if ((*OPENINNERBOUNDARY == 'A') || (*OPENINNERBOUNDARY == 'a')) AccBoundary = YES;
  if ((*OPENINNERBOUNDARY == 'K') || (*OPENINNERBOUNDARY == 'k')) {
    KNOpen = YES;
    OpenInner = YES;
  }
  if ((*OPENINNERBOUNDARY == 'P') || (*OPENINNERBOUNDARY == 'p')) Alexboundary = YES;
  if ((*GRIDSPACING == 'L') || (*GRIDSPACING == 'l')) LogGrid = YES;
  if ((*DISK == 'N') || (*DISK == 'n')) IsDisk = NO;
  if ((*FRAME == 'C') || (*FRAME == 'c')) Corotating = YES;
  if ((*FRAME == 'G') || (*FRAME == 'g')) {
    Corotating = YES;
    GuidingCenter = YES;
  }
  if ((*FRAME == 'B') || (*FRAME == 'b')) {
    Corotating = YES;
    BinaryCenter = YES;
  }
  if ((*RETROGRADEBINARY == 'Y') || (*RETROGRADEBINARY == 'y')) RetrogradeBinary = YES;
  if ((*WRITEVELOCITY == 'N') || (*WRITEVELOCITY == 'n')) Write_Velocity = NO;
  if ((*WRITEDENSITY == 'N') || (*WRITEDENSITY == 'n')) Write_Density = NO;
  if ((*WRITEENERGY == 'Y') || (*WRITEENERGY == 'y')) Write_Energy = YES;
  if ((*WRITETEMPERATURE == 'Y') || (*WRITETEMPERATURE == 'y')) Write_Temperature = YES;
  if ((*WRITEONEDFIELDS == 'Y') || (*WRITEONEDFIELDS == 'y')) Write_OneD_Fields = YES;
  if ((*WRITEDIVV == 'Y') || (*WRITEDIVV == 'y')) Write_DivV = YES;
  if ((*WRITEVISCHEAT == 'Y') || (*WRITEVISCHEAT == 'y')) Write_ViscHeat = YES;
  if ((*WRITETHERHEAT == 'Y') || (*WRITETHERHEAT == 'y')) Write_TherHeat = YES;
  if ((*WRITETHERCOOL == 'Y') || (*WRITETHERCOOL == 'y')) Write_TherCool = YES;
  if ((*WRITERADDIFF == 'Y') || (*WRITERADDIFF == 'y')) Write_RadDiff = YES;
  if ((*WRITESTARIRRAD == 'Y') || (*WRITESTARIRRAD == 'y')) Write_StarIrrad = YES;
  if ((*WRITEPDV == 'Y') || (*WRITEPDV == 'y')) Write_pdv = YES;
  if ((*WRITEARTVISC == 'Y') || (*WRITEARTVISC == 'y')) Write_ArtVisc = YES;
  if ((*WRITEOPACITY =='Y') || (*WRITEOPACITY == 'y')) Write_Opacity = YES;
  if ((*WRITEPOTENTIAL == 'Y') || (*WRITEPOTENTIAL == 'y')) Write_Potential = YES;
  if ((*WRITETEST == 'Y') || (*WRITETEST == 'y')) Write_Test = YES;
  if ((*WRITEGR == 'Y') || (*WRITEGR == 'y')) Write_gr = YES;
  if ((*WRITEGTHETA == 'Y') || (*WRITEGTHETA == 'y')) Write_gtheta = YES;
  if ((*SOFTWRITING == 'y') || (*SOFTWRITING == 'Y')) SoftWriting = YES;
  if ((*RETROGRADEPLANET == 'y') || (*RETROGRADEPLANET == 'Y')) RetrogradePlanet = YES;
  if ((*ADDMASS == 'y') || (*ADDMASS == 'Y')) AddMass = YES;
  if ((*EXPONENTIALDECAY == 'y') || (*EXPONENTIALDECAY == 'Y')) {
    ExponentialDecay = YES;
    CentrifugalBalance = YES;
  }
  if ((*DISCEVAPORATION == 'y') || (*DISCEVAPORATION == 'Y')) {
    DiscEvaporation = YES;
    masterprint("Disc evaporation included, with characteristic timescale = %lg\n", TEVAP);
  }
  if ((*CUSTIT == 'y') || (*CUSTIT == 'Y')) {
    CustomizedIT = YES;
  }
  if ((*DAMPTOINI == 'n') || (*DAMPTOINI == 'N')) DampToIni = NO;
  if ((*DAMPTOAXI == 'y') || (*DAMPTOAXI == 'Y')) {
    DampToAxi = YES;
    DampToIni = NO;
  }
  if ((*CUTFORCES == 'y') || (*CUTFORCES == 'Y')) CutForces = YES;
  if ((*COROTATEWITHOUTERPLANET == 'y') || (*COROTATEWITHOUTERPLANET == 'Y')) CorotateWithOuterPlanet = YES;
  if ((*INDIRECTTERM == 'N') || (*INDIRECTTERM == 'n')) Indirect_Term = NO;
  if ((*SELFGRAVITY == 'Y') || (*SELFGRAVITY == 'y')) SelfGravity = YES;
  if ((*SELFGRAVITY == 'Z') || (*SELFGRAVITY == 'z')) {
    SelfGravity = YES;
    SGZeroMode = YES;
  }
  if ((*ZMPLUS == 'Y') || (*ZMPLUS == 'y')) ZMPlus = YES;
  if ( (ZMPlus) && (!SGZeroMode) ) {
    masterprint ("This is not very meaningfull to involve the anisotropic pressure model (ZMPlus=Yes) without taking into account the axisymmetric component of the disk self-gravity. I decided to put ZMPlus = No. Please check again!");
    ZMPlus = NO;
  }
  if ((*ADDNOISE == 'Y') || (*ADDNOISE == 'y')) {
    AddNoise = YES;
    srand48(time(NULL));
  }
  if ((*ENERGYEQUATION == 'Y') || (*ENERGYEQUATION == 'y')) {
    EnergyEquation = YES;
    Write_Temperature = YES;
  }
  if ((*FARGOPLANETE == 'y') || (*FARGOPLANETE == 'Y')) {
    FargoPlanete = YES;
  }
  if ((*MODIFIEDSOUNDSPEED == 'Y') || (*MODIFIEDSOUNDSPEED == 'y')) {
    ModifiedSoundSpeed = YES;
    Write_Temperature = YES;
  }
  if ((*THERMALDIFFUSION == 'Y') || (*THERMALDIFFUSION == 'y')) ThermalDiffusion = YES; //Entropy diffusion (Clement)
  if ((*THERMALCOOLING == 'Y') || (*THERMALCOOLING == 'y')) ThermalCooling = YES;
  if ((*VISCOUSHEATING == 'N') || (*VISCOUSHEATING == 'n')) ViscousHeating = NO;
  if ((*TEMPPRESC == 'Y') || (*TEMPPRESC == 'y')) TempPresc = YES;
  if ((*BETACOOLING == 'Y') || (*BETACOOLING == 'y')) BetaCooling = YES;
  if ((*STELLARIRRADIATION == 'Y') || (*STELLARIRRADIATION == 'y'))  StellarIrradiation = YES;  //This is the model with a background temperature (Clement)
  if ( (*STARIRRAD == 'Y') || (*STARIRRAD == 'y'))  IrradStar = YES; //This is the model with calculating direct heating (Sareh)
  if (StellarIrradiation && IrradStar){
    mastererr("You cannot have StellarIrradiation and StarIrrad together.\n");
    prs_exit(1);
  }
  if ((EnergyEquation) && (ADIABATICINDEX == 1)) {
    masterprint ("You cannot have EnergyEquation = YES and AdiabatcIndex = 1.");
    masterprint ("I decided to put EnergyEquation = No, to simulate a locally isothermal equation of state.");
    masterprint ("Please check if it is what you really wanted to do!\n");
    EnergyEquation = NO;
  }
  if (!EnergyEquation) {
      ThermalDiffusion = NO;
      ThermalCooling = NO;
      BetaCooling = NO;  
      StellarIrradiation = NO; 
      IrradStar = NO;
      ADIABATICINDEX = 1.0;
  }
  if ((*WRITEENERGY == 'N') || (*WRITEENERGY == 'n')) Write_Energy = NO;
  if ((*MHD == 'L') || (*MHD == 'l')) {
    MHDLSA = YES;
    srand48(time(NULL));
  }
  if ((*HIGHMCUTOFF == 'Y') || (*HIGHMCUTOFF == 'y')) {
    HighMCutoff = YES;
  }
  if ((*IMPOSEDDENSITY == 'y') || (*IMPOSEDDENSITY == 'Y')) ImposedDensity = YES;
  if ((*IMPOSEDALPHA == 'y') || (*IMPOSEDALPHA == 'Y')) {
    ViscosityAlpha = YES;
    masterprint ("Viscosity is of alpha type, and is imposed as a fixed radial profile\n");
    ImposedAlpha = YES;
  }
  if ((*EXCLUDEHILL == 'Y') || (*EXCLUDEHILL == 'y')) ExcludeHill = YES;
  if ((*CICPLANET == 'Y') || (*CICPLANET == 'y')) CICPlanet = YES;
  if ((*FORCEDCIRCULAR == 'Y') || (*FORCEDCIRCULAR == 'y')) ForcedCircular = YES;
  if ((*FORCEDINNERCIRCULAR == 'Y') || (*FORCEDINNERCIRCULAR == 'y')) ForcedInnerCircular = YES;
  if ((*FORCEDCIRCULARTEMPORARY== 'Y') || (*FORCEDCIRCULARTEMPORARY == 'y')){
		 ForcedCircularTemp = YES;
     if (RELEASEDATE == 0.0){
		 		mastererr ("ForcedCircularTemporary needs ReleaseDate larger than 0.\n");
     prs_exit (1);
     }
  }
  if ((*READPLANETFILEATRESTART == 'N') || (*READPLANETFILEATRESTART == 'n')) ReadPlanetFileAtRestart = NO;
  if ((*DONTAPPLYSUBKEPLERIAN == 'Y') || (*DONTAPPLYSUBKEPLERIAN == 'y')) {
    masterprint ("================================================\n");
    masterprint ("I will not apply subkeplerian boundary on vtheta\n");
    masterprint ("================================================\n");
    DontApplySubKeplerian = YES;
  }
  if ((*CUSTTQEXC == 'Y') || (*CUSTTQEXC == 'y')) {
    CustTqExc = YES;
  }
  if (Evanescent) {
    masterprint ("Evanescent wave-killing zones boundary condition is applied,\n");
    masterprint ("I will therefore not apply subKeplerian boundary condition on vtheta.\n");
    DontApplySubKeplerian = YES;
  }
  if ( (EXCLUDEHILLFACTOR < 0.0) || (EXCLUDEHILLFACTOR > 1.0) ) {
    mastererr ("EXCLUDEHILLFACTOR must range between 0 and 1.\n");
    prs_exit (1);
  }
  if ((ALPHAVISCOSITY != 0.0) && (VISCOSITY != 0.0)) {
    mastererr ("You cannot use at the same time\n");
    mastererr ("VISCOSITY and ALPHAVISCOSITY.\n");
    mastererr ("Edit the parameter file so as to remove\n");
    mastererr ("one of these variables and run again.\n");
    prs_exit (1);
  }
  if (ALPHAVISCOSITY != 0.0) {
    ViscosityAlpha = YES;
    masterprint ("Viscosity is of alpha type\n");
  }
  
  if ((THICKNESSSMOOTHING != 0.0) && (ROCHESMOOTHING != 0.0)) {
    mastererr ("You cannot use at the same time\n");
    mastererr ("`ThicknessSmoothing' and `RocheSmoothing'.\n");
    mastererr ("Edit the parameter file so as to remove\n");
    mastererr ("one of these variables and run again.\n");
    prs_exit (1);
  } 
  if ((THICKNESSSMOOTHING <= 0.0) && (ROCHESMOOTHING <= 0.0)) {
    mastererr ("A non-vanishing potential smoothing length is required.\n");
    mastererr ("Please use either of the following variables:\n");
    mastererr ("`ThicknessSmoothing' *or* `RocheSmoothing'.\n");
    mastererr ("before launching the run again.\n");
    prs_exit (1);
  }
  if (SelfGravity && (SGTHICKNESSSMOOTHING == 0.0)) {
    mastererr ("You cannot have a vanishing smoothing length for \n");
    mastererr ("the self-gravitating kernel. Please modify and run again.\n");
    prs_exit (1);
  }
  if (ROCHESMOOTHING != 0.0) {
    RocheSmoothing = YES;
    masterprint ("Planet potential smoothing scales with their Hill sphere.\n");
  } else {
    if ((*SMOOTHINGAT == 'P') || (*SMOOTHINGAT == 'p')){
      SmoothAtPlanet = YES;
      masterprint("Smoothing length will be calculated at the location of the planet.\n");
    }
  }
  if (MCRIFACTOR != 0.0) masterprint("The critical mass factor you set is %e\n", MCRIFACTOR );
  if ((MCRIFACTOR == 0.0) && (FargoPlanete == YES)) {
    mastererr("You should set a critical mass factor when using Planete code with Fargo.\n");
    prs_exit(1);
  }
  if ((*PHOTOEVAPORATION == 'Y') || (*PHOTOEVAPORATION == 'y')) {
    PhotoEvapor = YES;
  }
  if ( (PhotoEvapor == YES) && ( LX == 0.0)) {
    mastererr("You should set the star X-ray luminosity in erg/s for photoevaporation.\n");
    mastererr("Don't also forget FACTORUNITMASS.\n");
    prs_exit(1);
  }
  if (OverridesOutputdir == YES) {
    sprintf (OUTPUTDIR, "%s", NewOutputdir);
  }
  if (Evanescent || Alexboundary){
    if ( (WKZRMIN == 0.0) || (WKZRMAX == 0.0) ) {
      mastererr ("Evanescent 'wave-killing zones' assumed as a boundary condition, \n");
      mastererr ("but you did not specify the radii of the border zones. Please run again.\n");
      mastererr ("by setting the values for WKZRMIN and WKZRMAX");
      prs_exit (1);
    }
  }
  if (MixedBC) {
    WKZRMIN = 0.0;
    if (WKZRMAX == 0.0) {
      mastererr ("Evanescent 'wave-killing zones' assumed as the outer boundary condition, \n");
      mastererr ("but you did not specify the radii of the border zones. Please run again.\n");
      mastererr ("by setting the value of WKZRMAX");
      prs_exit (1);
    }
  }
  if ((*MDOTHARTMANN == 'Y') || (*MDOTHARTMANN == 'y')){
    MdotHartmann = YES;
    masterprint("You want to impose accretion rate by Hartmann 1998 on the boundaries.\n");
  }
  if (AccBoundary) {
    if (!MdotHartmann){
      if ((MDOTINIT == 0.0) || (MDOTFINAL == 0.0) || (MDOTTIME == 1)){
        mastererr ("For having an accreting boundary, you need to set\n");
        mastererr ("MDOTINIT, MDOTFINAL, MDOTTIME in the .par file.\n");
        prs_exit (1);
      }
    } else if (MDOTTIME != 1.0) {
          mastererr("You can not use mdot Hartmann and user a definded value at the same time.\n");
    }
    if ((*ZEROINNER == 'D') || (*ZEROINNER == 'd')){
        DecInner = YES;
        masterprint("Your inner boundary will be decreased progressively.\n");
    } else if ((*ZEROINNER == 'O') || (*ZEROINNER == 'o')){
        OpInner = YES;
        masterprint("Your inner boundary will be open.\n");
    }
    if ((WKZRMIN == 0.0) || (WKZRMAX == 0.0) || (WKZTIN == 0.0) || (WKZTOUT == 0.0)){
       masterprint("some of the damping parameters have not been set.\n");
       prs_exit(1);
    }
  }
  if ((*WRITESIGDOT == 'Y') || (*WRITESIGDOT == 'y')) Write_Sigdot = YES;
  if (DIFFUSIVITY < 1e-13) {
      mastererr("In this version that we used Paardekooper 2011 formula for the");
      mastererr("torque calculation, the DIFFUSIVITY cannot be smaller than");
      mastererr("1e-13 because the bessel function fails. Please give a very small value");
      prs_exit(1);
  }
  /* Add a trailing slash to OUTPUTDIR if needed */
  if (*(OUTPUTDIR+strlen(OUTPUTDIR)-1) != '/')
    strcat (OUTPUTDIR, "/");
  /* Creat the output directory if not exist */
  if (stat(OUTPUTDIR, &stoutput) == -1){
     sprintf (CommandLine, "mkdir -p %s", OUTPUTDIR);
     system (CommandLine);
  }

  if ( (*EQUILIBRIUMINIT == 'Y') || (*EQUILIBRIUMINIT == 'y')){
      if ((MDOTINIT != 0.0) && (EnergyEquation)){
         InitEquilibrium = YES;
      } else {
         mastererr("For starting from equilibrium, you need to set MDOTINIT\n");
         mastererr("and include EneryEquation. Please check this condition\n");
         prs_exit(1);
      }
  }
  if ((*IMPOSEDGAMMA == 'y') || (*IMPOSEDGAMMA == 'Y')){
    if (GAMMAVALUE == 0.0) {
      masterprint("You decided to impose a given torque on the planet but you forgot to give the GAMMAVALUE.\n");
      prs_exit(1);
    }
  }
  if ((*WRITEONEDVISC == 'Y') || (*WRITEONEDVISC == 'y'))
    Write_OneD_Viscosity = YES;
  if ( (*WRITEFORCES == 'Y') || (*WRITEFORCES == 'y'))
    WriteForces = YES;
}

void PrintUsage (execname)
     char *execname;
{
  mastererr("Usage : %s [-abcdeimnptvz] [-(0-9)] [-s number] [-f scaling] parameters file\n", execname);
  mastererr("\n-a : Monitor mass and angular momentum at each timestep\n");
  mastererr("-b : Adjust azimuthal velocity to impose strict centrifugal balance at t=0\n");
  mastererr("-c : Sloppy CFL condition (checked at each DT, not at each timestep)\n");
  mastererr("-d : Print some debugging information on 'stdout' at each timestep\n");
  mastererr("-e : Activate EU test problem torque file output\n");
  mastererr("-f : Scale density array by 'scaling'. Useful to increase/decrease\n");
  mastererr("     disk surface density after a restart, for instance.            \n");
  mastererr("-i : tabulate Sigma profile as given by restart files\n");
  mastererr("-m : Merge output files from different CPUs\n");
  mastererr("-n : Disable simulation. The program just reads parameters file\n");
  mastererr("-o : Overrides output directory of input file.\n");
  mastererr("-p : Give profiling information at each time step\n");
  mastererr("-s : Restart simulation, taking #'number' files as initial conditions\n");
  mastererr("-t : Monitor CPU time usage at each time step\n");
  mastererr("-v : Verbose mode. Tells everything about parameters file\n");
  mastererr("-x : The planets mass is imposed by an external file (project) \n");
  mastererr("-y : Impose duration over which simulation is carried over (project) \n");
  mastererr("-z : fake sequential built when evaluating sums on HD meshes\n");
  mastererr("-(0-9) : only write initial (or restart) HD meshes,\n");
  mastererr("     proceed to the next nth output and exit\n");
  mastererr("     This option must stand alone on one switch (-va -4 is legal, -v4a is not)\n");
  prs_exit (1);
}

real TellNbOrbits (time)
     real time;
{
  return time/(2.0*M_PI)*sqrt(G*1.0/1.0/1.0/1.0);
}

real TellNbOutputs (time)
     real time;
{
  return (time/DT/NINTERM);
}

void TellEverything () {
  real temp, nbfileoutput;
  if (!CPU_Master) return;
  printf ("\nDisc properties:\n");
  printf ("----------------\n");
  printf ("Inner Radius          : %g\n", RMIN);
  printf ("Outer Radius          : %g\n", RMAX);
  printf ("Aspect Ratio          : %g\n", ASPECTRATIO);
  printf ("VKep at inner edge    : %.3g\n", sqrt(G*1.0*(1.-0.0)/RMIN));
  printf ("VKep at outer edge    : %.3g\n", sqrt(G*1.0/RMAX));
  temp=(PMAX-PMIN)*SIGMA0/(2.0-SIGMASLOPE)*(pow(RMAX,2.0-SIGMASLOPE) - pow(RMIN,2.0-SIGMASLOPE));  /* correct this and what follows... */
  printf ("Initial Disk Mass             : %g\n", temp);
  temp=(PMAX-PMIN)*SIGMA0/(2.0-SIGMASLOPE)*(1.0 - pow(RMIN,2.0-SIGMASLOPE));
  printf ("Initial Mass inner to r=1.0  : %g \n", temp);
  temp=(PMAX-PMIN)*SIGMA0/(2.0-SIGMASLOPE)*(pow(RMAX,2.0-SIGMASLOPE) - 1.0);
  printf ("Initial Mass outer to r=1.0  : %g \n", temp);
  printf ("Travelling time for acoustic density waves :\n");
  temp = 2.0/3.0/ASPECTRATIO*(pow(RMAX,1.5)-pow(RMIN,1.5));
  printf (" * From Rmin to Rmax  : %.2g = %.2f orbits ~ %.1f outputs\n", temp, TellNbOrbits(temp), TellNbOutputs(temp));
  temp = 2.0/3.0/ASPECTRATIO*(pow(RMAX,1.5)-pow(1.0,1.5));
  printf (" * From r=1.0 to Rmax: %.2g = %.2f orbits ~ %.1f outputs\n", temp, TellNbOrbits(temp), TellNbOutputs(temp));
  temp = 2.0/3.0/ASPECTRATIO*(pow(1.0,1.5)-pow(RMIN,1.5));
  printf (" * From r=1.0 to Rmin: %.2g = %.2f orbits ~ %.1f outputs\n", temp, TellNbOrbits(temp), TellNbOutputs(temp));
  temp = (PMAX-PMIN)*sqrt(RMIN*RMIN*RMIN/G/1.0);
  printf ("Orbital time at Rmin  : %.3g ~ %.2f outputs\n", temp, TellNbOutputs(temp));
  temp = (PMAX-PMIN)*sqrt(RMAX*RMAX*RMAX/G/1.0);
  printf ("Orbital time at Rmax  : %.3g ~ %.2f outputs\n", temp, TellNbOutputs(temp));
  printf ("Sound speed :\n");
  printf (" * At unit radius     : %.3g\n", ASPECTRATIO*sqrt(G*1.0));
  printf (" * At outer edge      : %.3g\n", ASPECTRATIO*sqrt(G*1.0/RMAX));
  printf (" * At inner edge      : %.3g\n", ASPECTRATIO*sqrt(G*1.0/RMIN));
  printf ("\nGrid properties:\n");
  printf ("----------------\n");
  printf ("Number of rings       : %d\n", NRAD);
  printf ("Number of sectors     : %d\n", NSEC);
  printf ("Total cells           : %d\n", NRAD*NSEC);
  printf ("\nOutputs properties:\n");
  printf ("-------------------\n");
  printf ("Time increment between outputs : %.3f = %.3f orbits\n", NINTERM*DT, TellNbOrbits(NINTERM*DT));
  printf ("At each output #i, the following files are written:\n");
  printf ("gasdens[i].dat : %d bytes\n",(int)(GLOBALNRAD*NSEC*sizeof(real)));
  printf ("gasvrad[i].dat : %d bytes\n",(int)(GLOBALNRAD*NSEC*sizeof(real)));
  printf ("gasvtheta[i].dat : %d bytes\n",(int)(GLOBALNRAD*NSEC*sizeof(real)));
  if (EnergyEquation == YES)
    printf ("gasTemperature[i].dat : %d bytes\n",(int)(GLOBALNRAD*NSEC*sizeof(real)));
  if (AdvecteLabel == YES)
    printf ("gaslabel[i].dat : %d bytes\n",(int)(GLOBALNRAD*NSEC*sizeof(real)));
  printf ("There will be in total %d outputs\n", NTOT/NINTERM);
  printf ("(which correspond to an elapsed time = %.3f or to %.2f orbits)\n", NTOT*DT, TellNbOrbits(NTOT*DT));
  nbfileoutput = 3.0;
  if (EnergyEquation == YES)
    nbfileoutput += 1.0;
  if (AdvecteLabel == YES)
    nbfileoutput += 1.0;
  temp =nbfileoutput*GLOBALNRAD*NSEC*sizeof(real);
  temp *= (real)NTOT/(real)NINTERM;
  temp /= 1024.0*1024.0;
  printf ("So the code will produce ~%.2f Mbytes of data\n", temp);
  printf ("Check (eg by issuing a 'df' command) that you have enough disk space,\n");
  printf ("otherwise you will get a system full and the code will stop.\n");
  fflush (stdout);
}

void GiveTimeInfo (number)
     int number;
{
  struct tms buffer;
  real total, last, mean, totalu;
  Current = times (&buffer);
  CurrentUser = buffer.tms_utime;
  if (FirstStep == YES) {
    First = Current;
    FirstUser = CurrentUser;
    fprintf (stderr, "Time counters initialized\n");
    FirstStep = NO;
    Ticks = sysconf (_SC_CLK_TCK);
  } else {
    total = (real)(Current - First)/Ticks;
    totalu= (real)(CurrentUser-FirstUser)/Ticks;
    last  = (real)(CurrentUser - PreceedingUser)/Ticks;
    number -= begin_i/NINTERM;
    mean  = totalu / number;
    fprintf (stderr, "Total Real Time elapsed    : %.3f s\n", total);
    fprintf (stderr, "Total CPU Time of process  : %.3f s (%.1f %%)\n", totalu, 100.*totalu/total);
    fprintf (stderr, "CPU Time since last time step : %.3f s\n", last);
    fprintf (stderr, "Mean CPU Time between time steps : %.3f s\n", mean);
    fprintf (stderr, "CPU Load on last time step : %.1f %% \n", (real)(CurrentUser-PreceedingUser)/(real)(Current-Preceeding)*100.);
    
  }  
  PreceedingUser = CurrentUser;
  Preceeding = Current;
}

void InitSpecificTime (profiling, process_name, title)
     boolean profiling;
     TimeProcess *process_name;
     char *title;
{
  struct tms buffer;
  if (profiling == NO) return;
  Ticks = sysconf (_SC_CLK_TCK);
  times (&buffer);
  process_name->clicks = buffer.tms_utime;
  strcpy (process_name->name, title);
}

void GiveSpecificTime (profiling, process_name)
     boolean profiling;
     TimeProcess process_name;
{
  struct tms buffer;
  long ticks;
  real t;
  if (profiling == NO) return;
  Ticks = sysconf (_SC_CLK_TCK);
  times (&buffer);
  ticks = buffer.tms_utime - process_name.clicks;
  t = (real)ticks / (real)Ticks;
  fprintf (stderr, "Time spent in %s : %.3f s\n", process_name.name, t);
}

