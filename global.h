int CPU_Rank;
int CPU_Number;
boolean CPU_Master;
int CPU_Next, CPU_Prev, CPU_Highest;
real unit_mass, unit_length, unit_temperature, unit_time, mmw, sigma_SB;
/* ------------------------------------- */
/* Variables specific to fftw mesh split */
/* ------------------------------------- */
int CPU_Friend, CPU_NoFriend;
real *dens_friend;
real *SGP_buffft_Accr_friend, *SGP_buffft_Acct_friend;
real *ffttohydro_transfer, *ffttohydro_transfer_friend;
int local_Nx, local_i_start, local_i_start_friend, total_local_size_friend, local_Nx_friend;
int local_Ny_after_transpose, local_j_start_after_transpose;
int total_local_size, ifront, Zero_or_active_friend;
int transfer_size, transfer_size_friend;
int hydro_totalsize, active_hydro_totalsize, active_hydro_totalsize_friend;  
/* ------------------------------------- */
int IMIN;
int IMAX;
int Zero_or_active;
int Max_or_active;
int One_or_active;
int MaxMO_or_active;		/* MO: Minus One */
int GLOBALNRAD;
int ievaporation = 0;
real Rinf[MAX1D], Rsup[MAX1D], Rmed[MAX1D], Surf[MAX1D];
real InvRmed[MAX1D], InvSurf[MAX1D], InvDiffRmed[MAX1D];
real InvDiffRsup[MAX1D], InvRinf[MAX1D], Radii[MAX1D], GlobalRmed[MAX1D], azimuth[MAX1D];
real SigmaMed[MAX1D], SigmaInf[MAX1D];
real EnergyMed[MAX1D], PrescTimeMed[MAX1D];
real FinalPlanetMass[MAX1D], PlanetMassAtRestart[MAX1D];
real VMed[MAX1D];
real GLOBAL_ImposedAlpha[MAX1D];
real OmegaFrame, PhysicalTime=0.0, PhysicalTimeInitial, FinalTime;
int TimeStep=0;
real HillRadius, mdcp, mdcp0, exces_mdcp;
real axics[MAX1D];  // Azimuthally avaraged sound speed
real axidens[MAX1D];  // Azimuthally avaraged surface density
real axitemp[MAX1D];  // Azimuthally avaraged temperature
boolean Merge, AdvecteLabel, FakeSequential, MonitorIntegral, debug, OnlyInit;
boolean	GotoNextOutput, StoreSigma, StoreEnergy, ViscosityAlpha, RocheSmoothing;
boolean CentrifugalBalance, ExcludeHill, SloppyCFL;
MPI_Status fargostat;
PolarGrid *CellAbscissa, *CellOrdinate;
PolarGrid *RhoStar, *RhoInt, *Potential, *TurbPotential, *Pressure, *SoundSpeed, *Temperature;
PolarGrid *DivergenceVelocity, *TAURR, *TAUPP, *TAURP, *ViscHeat, *ThermHeat, *ThermCool, *Opacity, *ArtViscHeat, *pdvEnergy;
PolarGrid *RadDiffusion, *StarIrradiation;
PolarGrid *gr, *gtheta;
PolarGrid *Test;
boolean LogGrid;
boolean OverridesOutputdir;
char NewOutputdir[1024];
real *potturb;
real VthetaMed[MAX1D];
real VradMed[MAX1D];
real SigmaDotW[MAX1D];
real SigmaDotV[MAX1D];
real floordens;
real EvapMass=0.0, DiskMass=0.0, AccMassPls=0.0;
real Rhole, sigcrit;
int ihole;
real Menvelope[MAX1D], MdotEnvelope[MAX1D], MenvRemoved[MAX1D], MenvRemained[MAX1D], MenvAccreted[MAX1D];
int MenvCount[MAX1D];
boolean IrradStar;
real qirrad[MAX1D], tstar, rstar;
boolean TeffCal;
boolean CorrectVgasSG;
