/** \file proto.h

Declaration of all the functions of the FARGO code

*/

void masterprint (const char *template, ...);
void mastererr (const char *template, ...);
real GetGlobalIFrac ();
void prs_exit ();
void *prs_malloc ();
void erreur ();
void message ();
PolarGrid    *CreatePolarGrid();
void MultiplyPolarGridbyConstant ();
void DumpSources ();
Force *AllocateForce ();
void ComputeForce ();
void FreeForce ();
void UpdateLog ();
void ReadfromFile ();
void InitLabel ();
void Initialization ();
void var ();
void ReadVariables ();
void PrintUsage ();
real TellNbOrbits ();
real TellNbOutputs ();
void TellEverything ();
void GiveTimeInfo ();
void InitSpecificTime ();
void GiveSpecificTime ();
void EmptyPlanetSystemFile ();
void WritePlanetFile ();
void WritePlanetSystemFile ();
void WriteBigPlanetFile ();
void WriteBigPlanetSystemFile ();
real GetfromPlanetFile ();
void RestartPlanetarySystem ();
void WriteDiskPolar();
void WriteDim ();
void SendOutput ();
void ComputeIndirectTerm ();
void FillForcesArrays ();
void AdvanceSystemFromDisk ();
void AdvanceSystemRK5 ();
void SolveOrbits ();
real ConstructSequence ();
void InitGasVelocities ();
void AccreteOntoPlanets ();
void FindOrbitalElements ();
int FindNumberOfPlanets ();
PlanetarySystem *AllocPlanetSystem ();
void FreePlanetary ();
PlanetarySystem *InitPlanetarySystem ();
void ListPlanets ();
real GetPsysInfo ();
void RotatePsys ();
void DerivMotionRK5 ();
void TranslatePlanetRK5 ();
void RungeKunta ();
real GasTotalMass();
real GasMomentum ();
real GasTotalEnergy ();
void DivisePolarGrid ();
void InitComputeAccel ();
Pair ComputeAccel ();
void OpenBoundary ();
void NonReflectingBoundary ();
void ApplyOuterSourceMass ();
void ApplyBoundaryCondition ();
void ApplySubKeplerianBoundary ();
void CorrectVtheta ();
boolean DetectCrash ();
void FillPolar1DArrays ();
void InitEuler ();
real min2 ();
real max2 ();
void ActualiseGas ();
void AlgoGas ();
void SubStep1 ();
void SubStep2 ();
void SubStep3 ();
int ConditionCFL ();
real Sigma ();
void FillSigma ();
void RefillSigma ();
real Energy ();
void FillEnergy ();
void RefillEnergy ();
real PrescTime();
void FillPrescTime();
void Transport ();
void OneWindRad ();
void ComputeThetaElongations ();
void ComputeAverageThetaVelocities ();
void ComputeResiduals ();
void ComputeConstantResidual ();
void AdvectSHIFT ();
void OneWindTheta ();
void QuantitiesAdvection ();
void ComputeExtQty ();
void ComputeSpeQty ();
void InitTransport () ;
void ComputeStarRad ();
void ComputeStarTheta ();
void ComputeLRMomenta ();
void ComputeVelocities ();
real VanLeerRadial ();
void VanLeerTheta ();
void InitViscosity ();
void ComputeViscousTerms ();
void UpdateVelocitiesWithViscosity ();
void AllocateComm ();
void CommunicateBoundaries ();
void handfpe();
void setfpe ();
void merge ();
void ReadPrevDim ();
void CheckRebin ();
void SplitDomain ();
void InitVariables();
real FViscosity ();
real AspectRatio ();
Force ComputeForceStockholm ();
void UpdateLogStockholm ();
Pair MassInOut ();
void CheckMomentumConservation ();
void InitGasDensity();
void InitGasEnergy();
void mpi_make1Dprofile ();
void ComputeSoundSpeed ();
void ComputePressureField ();
void ComputeTemperatureField ();
real CircumPlanetaryMass ();
real compute_smoothing ();
void EvanescentBoundary ();
void ComputeViscousHeating ();
void ComputeThermalDiffusion ();
void ComputeThermalCooling ();
void MakeDir ();
FILE *fopenp ();
void ApplyLSAOnPotential ();
void DampDensity ();
void InitImposedDensity ();
void InitImposedAlpha ();
void ComputeOpacities ();
void ComputeCodeUnits ();
void SubStep0();
real ComputeBetaCooling();
void Evaporation();
void Write1DFields();
real PhotoEvaporation();
void AccretingBoundary();
Pair AccelFromFormula();
int ReturnIndex();
real Ffunc();
real Gfunc();
real Kfunc();
void WriteSigmaDotFile();
real CalculateFracMass();
void WriteMassTrack();
real GetfromTrackMassFile();
void SetRhoFloor();
void SetEnergyFloor();
void mpi_Find1Dprofile();
void LinearInterpole();
real opLBL94();
void EqInitialize();
real CalculateTNew();
real HighMetal();
real LowMetal();
real BitschTemperature();
void ComputeStarIrrad();
void Write1DViscosity();
void CalculateAxiSGDiskPotentialTools();
void CalculateAs();
real rf();
real ellf();
real elle();
real rd();
void nrerror();
void make_azi_average_profile();
double lin();
real OpacityOriginal();
