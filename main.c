/** \file main.c

Main file of the distribution. Manages the call to initialization
functions, then the main loop.

*/

#include "mp.h"

boolean         Restart = NO, OpenInner = NO;
boolean         OneDRun = NO;
int             begin_i = 0, NbRestart = 0, verbose = NO;
int             dimfxy = 2;
static int      StillWriteOneOutput;
extern real     LostMass;
extern boolean  Corotating, ReadPlanetFileAtRestart, Write_OneD_Fields;
extern boolean  CorotateWithOuterPlanet, DiscEvaporation, FargoPlanete;
extern boolean  SelfGravity, SGZeroMode, EnergyEquation, SoftWriting;
real            ScalingFactor = 1.0;
real            Runtime = 0.0;
extern boolean  Write_Sigdot;
real            Mdisc0 = 0.0;


int
main(argc, argv)
     int argc;
     char *argv[];
{
  PolarGrid   *gas_density;
  PolarGrid   *gas_v_rad; 
  PolarGrid   *gas_v_theta; 
  PolarGrid   *gas_energy; 
  PolarGrid   *gas_label;
  int          i;
  real         foostep = 0.;
  real         r, v;
  boolean      disable = NO, TimeInfo = NO, Profiling = NO;
  boolean      Stockholm = NO, updatevelocities = NO;
  TimeProcess  t_Hydro;
  char         ParameterFile[256];
  PlanetarySystem *sys;
  Force *force;
  MPI_Init (&argc, &argv);
  MPI_Comm_rank (MPI_COMM_WORLD, &CPU_Rank);
  MPI_Comm_size (MPI_COMM_WORLD, &CPU_Number);
  CPU_Master = (CPU_Rank == 0 ? 1 : 0);
  setfpe ();  /* Control behavior for floating point
     exceptions trapping (default is not to do anything) */
  if (argc == 1) PrintUsage (argv[0]);
  strcpy (ParameterFile, "");
  for (i = 1; i < argc; i++) {
    if (*(argv[i]) == '-') {
      if (strspn (argv[i], "-xywsecndovtpfamzib0123456789")
       != strlen (argv[i]))
      PrintUsage (argv[0]);
      if (strchr (argv[i], 'n'))
        disable = YES;
      if (strchr (argv[i], 'e'))
        Stockholm = YES;
      if (strchr (argv[i], 'v'))
        verbose = YES;
      if (strchr (argv[i], 't'))
        TimeInfo = YES;
      if (strchr (argv[i], 'c'))
        SloppyCFL = YES;
      if (strchr (argv[i], 'p'))
        Profiling = YES;
      if (strchr (argv[i], 'd'))
        debug = YES;
      if (strchr (argv[i], 'b'))
        CentrifugalBalance = YES;
      if (strchr (argv[i], 'm'))
        Merge = YES;
      if (strchr (argv[i], 'a'))
        MonitorIntegral = YES;
      if (strchr (argv[i], 'z'))
        FakeSequential = YES;
      if (strchr (argv[i], 'i')) {
        StoreSigma = YES;
        if (EnergyEquation)
          StoreEnergy = YES;
      }
      if (strchr (argv[i], '0'))
        OnlyInit = YES;
      if ((argv[i][1] >= '1') && (argv[i][1] <= '9')) {
        GotoNextOutput = YES;
        StillWriteOneOutput = (int)(argv[i][1]-'0');
      }
      if (strchr (argv[i], 's')) {
        Restart = YES;
        i++;
        NbRestart = atoi(argv[i]);
        if ((NbRestart < 0)) {
          masterprint ("Incorrect restart number\n");
          PrintUsage (argv[0]);
        }
      }
      if (strchr (argv[i], 'x')) {
        /* Option where Fargo is run over a specified running time; applies 
        to Planet/Fargo codes coupling. We then set Write_1D_fields to yes! */
        FargoPlanete = YES;
      }
      if (strchr (argv[i], 'w')) {
        /* Option where Fargo is run over a specified running time; applies 
        to Planet/Fargo codes coupling. We then set Write_1D_fields to yes! */
        OneDRun = YES;
      }
      if (strchr (argv[i], 'y')) {
        /* Option where Fargo is run over a specified running time; applies 
        to Planet/Fargo codes coupling. We then set Write_1D_fields to yes! */
        Write_OneD_Fields = YES;
        i++;
        Runtime = atof(argv[i]);
        if ((Runtime < 0.0)) {
          masterprint ("Incorrect NTOT\n");
          PrintUsage (argv[0]);
        }
      }
      if (strchr (argv[i], 'o')) {
        OverridesOutputdir = YES;
        i++;
        sprintf (NewOutputdir, "%s", argv[i]);
      } else {
        if (strchr (argv[i], 'f')) {
          i++;
          ScalingFactor = atof(argv[i]);
          masterprint ("Scaling factor = %g\n", ScalingFactor);
          if ((ScalingFactor <= 0)) {
            masterprint ("Incorrect scaling factor\n");
            PrintUsage (argv[0]);
          }
        }
      }
    } else {
     strcpy (ParameterFile, argv[i]);
   }
  }
  if ( (StoreSigma || StoreEnergy) && !(Restart)) {
    mastererr ("You cannot use tabulated surface density\n");
    mastererr ("or surface internal energy in a non-restart run.\n");
    mastererr ("Aborted\n");
    prs_exit (0);
  }
  if (ParameterFile[0] == 0)
   PrintUsage (argv[0]);
  ReadVariables (ParameterFile);
  SplitDomain ();
  if (verbose == YES) 
    TellEverything ();
  if (disable == YES)
    prs_exit (0);
  DumpSources (argc, argv, ParameterFile, PLANETCONFIG);
  ComputeCodeUnits ();
  masterprint ("Allocating arrays...");
  fflush (stdout);
  gas_density        = CreatePolarGrid(NRAD, NSEC, "dens");
  gas_v_rad          = CreatePolarGrid(NRAD, NSEC, "vrad");
  gas_v_theta        = CreatePolarGrid(NRAD, NSEC, "vtheta");
  gas_energy         = CreatePolarGrid(NRAD, NSEC, "energy");
  gas_label          = CreatePolarGrid(NRAD, NSEC, "label");
  SoundSpeed         = CreatePolarGrid(NRAD, NSEC, "SoundSpeed");
  masterprint ("done.\n");
  FillPolar1DArrays ();
  force = AllocateForce (dimfxy);
  
  /* Gas density initialization */
  InitComputeAccel ();
  InitGasDensity (gas_density);
  floordens = 1e-8; 
  sigcrit = floordens;
  Rhole = GlobalRmed[0];
  ihole = 0;
  /*This has a numerical reason: If the density gets 
  * smaller than 1e-8, some considerable oscillations happens in the 
  * radial velocity that we do not know the reason*/
  
  /* If energy equation is taken into account, we initialize the gas
     thermal energy */
  if ( EnergyEquation )
    InitGasEnergy (gas_energy);
  
  /* Here planets are initialized feeling star potential but they do
     not feel disk potential */
  sys = InitPlanetarySystem (PLANETCONFIG,NbRestart, gas_density, force);   
  //The axi-symmetric part of the SG acceleration on the disc Kley1996
  real *SGAarray;
  SGAarray = (real *)malloc(sizeof(real) * GLOBALNRAD*(GLOBALNRAD+1)); //rows: GLOBALNRAD (rmed locations), cols: GLOBALNRAD+1 (interfaces)
  if (CorrectVgasSG)
    CalculateAxiSGDiskPotentialTools(SGAarray);
  if ( SelfGravity ) {
    /* If SelfGravity = YES or Z, planets are initialized feeling disk
       potential. Only the surface density is required to calculate
       the radial self-gravity acceleration. The disk radial and
       azimutal velocities are not updated */
    /* If this is a restart simulation, we first need to read the 
       density at restart! */
    if (Restart == YES)
      ReadfromFile (gas_density, "gasdens", NbRestart);
    compute_selfgravity (gas_density, gas_v_rad, gas_v_theta, foostep, updatevelocities);
    init_planetarysys_withSG (sys);
  }
  if (Restart == NO)
    ListPlanets (sys);
  OmegaFrame = OMEGAFRAME;
  if (Corotating == YES) {
    if (!SelfGravity) {
      OmegaFrame = GetPsysInfo (sys, FREQUENCY);
    } else {
      if (!CorotateWithOuterPlanet) {
        r = sqrt( sys->x[0]*sys->x[0] + sys->y[0]*sys->y[0] );
        v = sqrt( sys->vx[0]*sys->vx[0] + sys->vy[0]*sys->vy[0] );
      } else {
        r = sqrt( sys->x[1]*sys->x[1] + sys->y[1]*sys->y[1] );
        v = sqrt( sys->vx[1]*sys->vx[1] + sys->vy[1]*sys->vy[1] );
      }
      /* If self-gravity is included, OmegaFrame is first calculated
      as the angular frequency the planet would have if it was on a
      fixed circular orbit at x=sys->x[0], with initial velocity
      sys->vy[0], which includes the radial initial
      self-gravitating acceleration */
      OmegaFrame = v / r;
    }
  }
  /* Only gas velocities remain to be initialized */
  Initialization (gas_density, gas_v_rad, gas_v_theta, gas_energy, gas_label, sys, SGAarray);
  /* Initial gas_density is used to compute the circumplanetary mass
     with initial density field */
  mdcp0 = CircumPlanetaryMass (gas_density, sys);
  
  Mdisc0 = GasTotalMass(gas_density);
  if ((Restart == YES) && (NbRestart != 0)) {
    if (!FargoPlanete)
      begin_i         = NbRestart * NINTERM;
    if (ReadPlanetFileAtRestart) {
      RestartPlanetarySystem (NbRestart, sys);
      ListPlanets (sys);
    }
    Mdisc0 = GetfromTrackMassFile (0, 3);
    LostMass = GetfromPlanetFile (NbRestart, 7, 0); /* 0 refers to planet #0 */
    PhysicalTime  = GetfromPlanetFile (NbRestart, 8, 0);
    if (sys->mass[0]>1e-21)  //Otherwise, it gets nan value and the code crashes
      OmegaFrame  = GetfromPlanetFile (NbRestart, 9, 0);
    DiskMass = GetfromTrackMassFile (NbRestart, 3);
    EvapMass = GetfromTrackMassFile (NbRestart, 4);
    AccMassPls = GetfromTrackMassFile (NbRestart, 5);
    Rhole = GetfromTrackMassFile (NbRestart, 6);
    ihole = ReturnIndex(Rhole);
  } else {      /* We initialize 'planet[i].dat' file */
    EmptyPlanetSystemFile (sys);
  }
  if (MonitorIntegral == YES)
    CheckMomentumConservation (gas_density, gas_v_theta, sys);
  PhysicalTimeInitial = PhysicalTime;
  if (Runtime > 0.0) {
    FinalTime = PhysicalTimeInitial + Runtime;
    masterprint ("Running time is imposed by user, equals %lg\n",Runtime);
  } else {
    FinalTime = NTOT * DT;
  }
  MultiplyPolarGridbyConstant (gas_density, ScalingFactor);
  /* OUTPUTS FIRST TIME */
  if (!FargoPlanete){
    TimeStep = NbRestart+1;
  } else {
    if (NbRestart == 0) {
      DiskMass = GasTotalMass(gas_density);
      SendOutput (TimeStep, gas_density, gas_v_rad, gas_v_theta, gas_energy, gas_label, sys);
      WritePlanetSystemFile (sys, TimeStep);
      WriteMassTrack (TimeStep, DiskMass, EvapMass, AccMassPls);
      WriteBigPlanetSystemFile (sys, TimeStep);
      SolveOrbits (sys,1);
      UpdateLog (force, sys, gas_density, gas_energy, TimeStep, PhysicalTime, dimfxy);
      TimeStep++;
    } else {    
      TimeStep = NbRestart+1;
    }
  }
  /* --------- */
  /* Time loop */
  /* --------- */
  if (FargoPlanete){
    while (PhysicalTime <= FinalTime) {
      /* Algorithm loop begins here */
      /***********************/
      /* Hydrodynamical Part */
      /***********************/
      InitSpecificTime (Profiling, &t_Hydro, "Eulerian Hydro algorithms");
      AlgoGas (force, gas_density, gas_v_rad, gas_v_theta, gas_energy, gas_label, sys, SGAarray);
      GiveSpecificTime (Profiling, t_Hydro);
      
      if (MonitorIntegral == YES) {
        CheckMomentumConservation (gas_density, gas_v_theta, sys);
        masterprint ("Gas Momentum   : %.18g\n", GasMomentum (gas_density, gas_v_theta));
        masterprint ("Gas total Mass : %.18g\n", GasTotalMass (gas_density));
        masterprint ("Gas total Energy : %.18g\n", GasTotalEnergy (gas_density, gas_v_rad, gas_v_theta, gas_energy));
      }
    }
    /* OUTPUTS AT LAST TIMESTEP */
    DiskMass = GasTotalMass(gas_density);
    SendOutput (TimeStep, gas_density, gas_v_rad, gas_v_theta, gas_energy, gas_label,sys);
    WritePlanetSystemFile (sys, TimeStep);
    WriteMassTrack (TimeStep, DiskMass, EvapMass, AccMassPls);
    WriteBigPlanetSystemFile (sys, TimeStep);
    SolveOrbits (sys,1);
    UpdateLog (force, sys, gas_density, gas_energy, TimeStep, PhysicalTime, dimfxy);
  } else {
    for (i = begin_i; i <= NTOT; i++) {
      if (NINTERM * (TimeStep = (i / NINTERM)) == i){
        DiskMass = GasTotalMass(gas_density);
        SendOutput (TimeStep, gas_density, gas_v_rad, gas_v_theta, gas_energy, gas_label,sys);
        WritePlanetSystemFile (sys, TimeStep);
        WriteMassTrack (sys, TimeStep, DiskMass, EvapMass, AccMassPls);
        if (Write_Sigdot)
          WriteSigmaDotFile(TimeStep);
      }
      /* Algorithm loop begins here */
      /***********************/
      /* Hydrodynamical Part */
      /***********************/
      InitSpecificTime (Profiling, &t_Hydro, "Eulerian Hydro algorithms");
      AlgoGas (force, gas_density, gas_v_rad, gas_v_theta, gas_energy, gas_label, sys, SGAarray);
      GiveSpecificTime (Profiling, t_Hydro);
      if (MonitorIntegral == YES) {
        CheckMomentumConservation (gas_density, gas_v_theta, sys);
        masterprint ("Gas Momentum   : %.18g\n", GasMomentum (gas_density, gas_v_theta));
        masterprint ("Gas total Mass : %.18g\n", GasTotalMass (gas_density));
        masterprint ("Gas total Energy : %.18g\n", GasTotalEnergy (gas_density, gas_v_rad, gas_v_theta, gas_energy));
      }
   }
 }
  /* --------------------- */
  FreePlanetary (sys);
  FreeForce (force);
  if ( SelfGravity && !SGZeroMode ) {
    rfftwnd_mpi_destroy_plan(SGP_fftplan_forward);
    rfftwnd_mpi_destroy_plan(SGP_fftplan_backward);
  }
  MPI_Finalize ();
  return 0;
}
