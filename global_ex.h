/********************************/
/*                              */
/* This file is created         */
/* automatically during         */
/* compilation. Do not edit.    */
/* See perl script              */
/* "varparser.pl" for details   */
/*                              */
/********************************/
extern int CPU_Rank;
extern int CPU_Number;
extern boolean CPU_Master;
extern int CPU_Next, CPU_Prev, CPU_Highest;
extern real unit_mass, unit_length, unit_temperature, unit_time, mmw, sigma_SB;
extern int CPU_Friend, CPU_NoFriend;
extern real *dens_friend;
extern real *SGP_buffft_Accr_friend, *SGP_buffft_Acct_friend;
extern real *ffttohydro_transfer, *ffttohydro_transfer_friend;
extern int local_Nx, local_i_start, local_i_start_friend, total_local_size_friend, local_Nx_friend;
extern int local_Ny_after_transpose, local_j_start_after_transpose;
extern int total_local_size, ifront, Zero_or_active_friend;
extern int transfer_size, transfer_size_friend;
extern int hydro_totalsize, active_hydro_totalsize, active_hydro_totalsize_friend;  
extern int IMIN;
extern int IMAX;
extern int Zero_or_active;
extern int Max_or_active;
extern int One_or_active;
extern int MaxMO_or_active;		/* MO: Minus One */
extern int GLOBALNRAD;
extern int ievaporation;
extern real Rinf[MAX1D], Rsup[MAX1D], Rmed[MAX1D], Surf[MAX1D];
extern real InvRmed[MAX1D], InvSurf[MAX1D], InvDiffRmed[MAX1D];
extern real InvDiffRsup[MAX1D], InvRinf[MAX1D], Radii[MAX1D], GlobalRmed[MAX1D], azimuth[MAX1D];
extern real SigmaMed[MAX1D], SigmaInf[MAX1D];
extern real EnergyMed[MAX1D], PrescTimeMed[MAX1D];
extern real FinalPlanetMass[MAX1D], PlanetMassAtRestart[MAX1D];
extern real VMed[MAX1D];
extern real GLOBAL_ImposedAlpha[MAX1D];
extern real OmegaFrame, PhysicalTime, PhysicalTimeInitial, FinalTime;
extern int TimeStep;
extern real HillRadius, mdcp, mdcp0, exces_mdcp;
extern real axics[MAX1D];  // Azimuthally avaraged sound speed
extern real axidens[MAX1D];  // Azimuthally avaraged surface density
extern real axitemp[MAX1D];  // Azimuthally avaraged temperature
extern real opaaxi[MAX1D];  // Azimuthally avaraged opacity
extern boolean Merge, AdvecteLabel, FakeSequential, MonitorIntegral, debug, OnlyInit;
extern boolean	GotoNextOutput, StoreSigma, StoreEnergy, ViscosityAlpha, RocheSmoothing;
extern boolean CentrifugalBalance, ExcludeHill, SloppyCFL;
extern MPI_Status fargostat;
extern PolarGrid *CellAbscissa, *CellOrdinate;
extern PolarGrid *RhoStar, *RhoInt, *Potential, *TurbPotential, *Pressure, *SoundSpeed, *Temperature;
extern PolarGrid *DivergenceVelocity, *TAURR, *TAUPP, *TAURP, *ViscHeat, *ThermHeat, *ThermCool, *Opacity, *ArtViscHeat, *pdvEnergy;
extern PolarGrid *RadDiffusion, *StarIrradiation;
extern PolarGrid *gr, *gtheta;
extern PolarGrid *Test;
extern boolean LogGrid;
extern boolean OverridesOutputdir;
extern char NewOutputdir[1024];
extern real *potturb;
extern real VthetaMed[MAX1D];
extern real VradMed[MAX1D];
extern real SigmaDotW[MAX1D];
extern real SigmaDotV[MAX1D];
extern real floordens;
extern real EvapMass, DiskMass, AccMassPls;
extern real Rhole, sigcrit;
extern int ihole;
extern real Menvelope[MAX1D], MdotEnvelope[MAX1D], MenvRemoved[MAX1D], MenvRemained[MAX1D], MenvAccreted[MAX1D];
extern int MenvCount[MAX1D];
extern boolean IrradStar;
extern real qirrad[MAX1D], tstar, rstar;
extern boolean TeffCal;
extern boolean CorrectVgasSG;
