/** \file var.c

Contains the function that connects the string of the parameter file
to global variables.  The var() function is found in Interpret.c

*/

#define __LOCAL
#include "mp.h"
#undef __LOCAL

void
InitVariables()
{
  var("DT", &DT, REAL, YES, "1.");
  var("SIGMA0", &SIGMA0, REAL, YES, "173.");
  var("NINTERM", &NINTERM, INT, YES, "10.");
  var("NTOT", &NTOT, INT, YES, "1501.");
  var("OUTPUTDIR", OUTPUTDIR, STRING, YES, "~masset");
  var("INNERBOUNDARY", OPENINNERBOUNDARY, STRING, NO, "WALL");
  var("LABELADVECTION", ADVLABEL, STRING, NO, "NO");
  var("TRANSPORT", TRANSPORT, STRING, NO, "FAST");
  var("PLANETCONFIG", PLANETCONFIG, STRING, NO, "Systems/SolarSystem.cfg");
  var("RADIALSPACING", GRIDSPACING, STRING, NO, "ARITHMETIC");
  var("NRAD", &NRAD, INT, YES, "64.0");
  var("NSEC", &NSEC, INT, YES, "64.0");
  var("RMIN", &RMIN, REAL, YES, "1.0");
  var("RMAX", &RMAX, REAL, YES, "1.0");
  var("THICKNESSSMOOTHING", &THICKNESSSMOOTHING, REAL, NO, "0.0");
  var("SGTHICKNESSSMOOTHING", &SGTHICKNESSSMOOTHING, REAL, NO, "0.0");
  var("ROCHESMOOTHING", &ROCHESMOOTHING, REAL, NO, "0.0");
  var("ASPECTRATIO", &ASPECTRATIO, REAL, YES, "0.05");
  var("VISCOSITY", &VISCOSITY, REAL, NO, "0.0");
  var("RELEASEVISCOSITY", &RELEASEVISCOSITY, REAL, NO, "0.0");
  var("RELEASEDATEVISCOSITY", &RELEASEDATEVISCOSITY, REAL, NO, "0.0");
  var("ALPHAVISCOSITY", &ALPHAVISCOSITY, REAL, NO, "0.0");  
  var("SIGMASLOPE", &SIGMASLOPE, REAL, YES, "0.0");  
  var("RELEASERADIUS", &RELEASERADIUS, REAL, NO, "0.0");  
  var("RELEASEDATE", &RELEASEDATE, REAL, NO, "0.0");  
  var("OMEGAFRAME", &OMEGAFRAME, REAL, NO, "0.0");
  var("DISK", DISK, STRING, NO, "YES");
  var("FRAME", FRAME, STRING, NO, "FIXED");
  var("OUTERSOURCEMASS", OUTERSOURCEMASS, STRING, NO, "NO");
  var("WRITEDENSITY", WRITEDENSITY, STRING, NO, "YES");
  var("WRITEVELOCITY", WRITEVELOCITY, STRING, NO, "YES");
  var("WRITEENERGY", WRITEENERGY, STRING, NO, "NO");
  var("WRITETEMPERATURE", WRITETEMPERATURE, STRING, NO, "NO");
  var("WRITEDIVV", WRITEDIVV, STRING, NO, "NO");
  var("WRITEPOTENTIAL", WRITEPOTENTIAL, STRING, NO, "NO");
  var("WRITEVISCHEAT", WRITEVISCHEAT, STRING, NO, "NO");
  var("WRITETHERHEAT", WRITETHERHEAT, STRING, NO, "NO");
  var("WRITETHERCOOL", WRITETHERCOOL, STRING, NO, "NO");
  var("WRITERADDIFF", WRITERADDIFF, STRING, NO, "NO");
  var("WRITESTARIRRAD", WRITESTARIRRAD, STRING, NO, "NO");
  var("WRITEPDV", WRITEPDV, STRING, NO, "NO");
  var("WRITEARTVISC", WRITEARTVISC, STRING, NO, "NO");
  var("WRITEOPACITY", WRITEOPACITY, STRING, NO, "NO");
  var("WRITETEST", WRITETEST, STRING, NO, "NO");
  var("INDIRECTTERM", INDIRECTTERM, STRING, NO, "YES");
  var("EXCLUDEHILL", EXCLUDEHILL, STRING, NO, "NO");
  var("EXCLUDEHILLFACTOR", &EXCLUDEHILLFACTOR, REAL, NO, "1.0");
  var("IMPOSEDDISKDRIFT", &IMPOSEDDISKDRIFT, REAL, NO, "0.0");
  var("FLARINGINDEX", &FLARINGINDEX, REAL, NO, "0.0");
  var("CAVITYRADIUS", &CAVITYRADIUS, REAL, NO, "0.0");
  var("CAVITYRATIO", &CAVITYRATIO, REAL, NO, "1.0");
  var("CAVITYWIDTH", &CAVITYWIDTH, REAL, NO, "1.0");
  var("TRANSITIONRADIUS", &TRANSITIONRADIUS, REAL, NO, "0.0");
  var("TRANSITIONRATIO", &TRANSITIONRATIO, REAL, NO, "1.0");
  var("TRANSITIONWIDTH", &TRANSITIONWIDTH, REAL, NO, "1.0");
  var("LAMBDADOUBLING", &LAMBDADOUBLING, REAL, NO, "0.0");
  var("SELFGRAVITY", SELFGRAVITY, STRING, NO, "NO");
  var("CICPLANET", CICPLANET, STRING, NO, "NO");
  var("FORCEDCIRCULAR", FORCEDCIRCULAR, STRING, NO, "NO");
  var("FORCEDINNERCIRCULAR", FORCEDINNERCIRCULAR, STRING, NO, "NO");
  var("ZMPLUS", ZMPLUS, STRING, NO, "NO");
  var("ENERGYEQUATION", ENERGYEQUATION, STRING, NO, "NO");
  var("MODIFIEDSOUNDSPEED", MODIFIEDSOUNDSPEED, STRING, NO, "NO");  
  var("ADIABATICINDEX", &ADIABATICINDEX, REAL, YES, "1.4");
  var("PLANETASPECTRATIO", &PLANETASPECTRATIO, REAL, NO, "0.5");
  var("THERMALDIFFUSION", THERMALDIFFUSION, STRING, NO, "NO");
  var("THERMALCOOLING", THERMALCOOLING, STRING, NO, "NO");
  var("VISCOUSHEATING", VISCOUSHEATING, STRING, NO, "YES");
  var("DIFFUSIVITY", &DIFFUSIVITY, REAL, NO, "0.000001");
  var("TEMPPRESC", TEMPPRESC, STRING, NO, "NO");
  var("PRESCTIME0", &PRESCTIME0, REAL, NO, "6.28");
  var("MHD", MHD, STRING, NO, "NO");
  var("SOFTWRITING", SOFTWRITING, STRING, NO, "NO");
  var("RETROGRADEPLANET", RETROGRADEPLANET, STRING, NO, "NO");
  var("GAMMATURB", &GAMMATURB, REAL, NO, "0.00001");
  var("LSAMODESPEEDUP", &LSAMODESPEEDUP, REAL, NO, "0.1");
  var("NBTURBMODES", &NBTURBMODES, INT, NO, "50");
  var("HIGHMCUTOFF", HIGHMCUTOFF, STRING, NO, "NO");
  var("PMIN", &PMIN, REAL, NO, "0.0");
  var("PMAX", &PMAX, REAL, NO, "6.2831853071795864");
  var("ADDMASS", ADDMASS, STRING, NO, "NO");
  var("BINARYSEPARATION", &BINARYSEPARATION, REAL, NO, "0.3");
  var("BINARYECCENTRICITY", &BINARYECCENTRICITY, REAL, NO, "0.0");
  var("RETROGRADEBINARY", RETROGRADEBINARY, STRING, NO, "NO");
  var("IMPOSEDDENSITY", IMPOSEDDENSITY, STRING, NO, "NO");
  var("IMPOSEDALPHA", IMPOSEDALPHA, STRING, NO, "NO");
  var("FACTORUNITMASS", &FACTORUNITMASS, REAL, NO, "1.0");
  var("FACTORUNITLENGTH", &FACTORUNITLENGTH, REAL, NO, "1.0");
  var("FACTORMMW", &FACTORMMW, REAL, NO, "1.0");
  var("BETACOOLING", BETACOOLING, STRING, NO, "NO");
  var("BETACOOLINGTIME", &BETACOOLINGTIME, REAL, NO, "10.0");
  var("BETACOOLINGSLOPE", &BETACOOLINGSLOPE, REAL, NO, "0.0");
  var("WRITEGR", WRITEGR, STRING, NO, "NO");
  var("WRITEGTHETA", WRITEGTHETA, STRING, NO, "NO");
  var("ADDNOISE", ADDNOISE, STRING, NO, "NO");
  var("WKZRMIN", &WKZRMIN, REAL, NO, "0.0");
  var("WKZRMAX", &WKZRMAX, REAL, NO, "0.0");
  var("WKZTIN", &WKZTIN, REAL, NO, "0.0");
  var("WKZTOUT", &WKZTOUT, REAL, NO, "0.0");
  var("READPLANETFILEATRESTART", READPLANETFILEATRESTART, STRING, NO, "YES");
  var("DONTAPPLYSUBKEPLERIAN", DONTAPPLYSUBKEPLERIAN, STRING, NO, "NO");
  var("CUSTTQEXC", CUSTTQEXC, STRING, NO, "NO");
  var("CUSTTQRAD", &CUSTTQRAD, REAL, NO, "0.0");
  var("DENSDAMPRAD", &DENSDAMPRAD, REAL, NO, "0.0");
  var("DAMPTOINI", DAMPTOINI, STRING, NO, "YES");
  var("DAMPTOAXI", DAMPTOAXI, STRING, NO, "NO");
  var("COROTATEWITHOUTERPLANET", COROTATEWITHOUTERPLANET, STRING, NO, "NO");
  var("CUTFORCES", CUTFORCES, STRING, NO, "NO");
  var("DISCEVAPORATION", DISCEVAPORATION, STRING, NO, "NO");
  var("TEVAP", &TEVAP, REAL, NO, "0.0");
  var("CUSTIT", CUSTIT, STRING, NO, "NO");
  var("TURBRMIN", &TURBRMIN, REAL, NO, "0.0");
  var("TURBRMAX", &TURBRMAX, REAL, NO, "0.0");
  var("WRITEONEDFIELDS", WRITEONEDFIELDS, STRING, NO, "NO");
  var("FARGOPLANETE", FARGOPLANETE, STRING, NO, "NO");
  var("MCRIFACTOR", &MCRIFACTOR, REAL, NO, "1e-20");  
  var("PHOTOEVAPORATION", PHOTOEVAPORATION, STRING, NO, "NO");
  var("LX", &LX, REAL, NO, "0.0");  
  var("MDOTINIT", &MDOTINIT, REAL, NO, "0.0");
  var("MDOTFINAL", &MDOTFINAL, REAL, NO, "0.0");
  var("MDOTTIME", &MDOTTIME, REAL, NO, "1.0");
  var("WRITESIGDOT", WRITESIGDOT, STRING, NO, "NO");
  var("MDOTHARTMANN", MDOTHARTMANN, STRING, NO, "NO");
  var("THARTMANN", &THARTMANN, REAL, NO, "0.0"); //The time that adjust the value of disc accretion rate to our desired one.
  var("ZEROINNER", ZEROINNER, STRING, NO, "NO"); //It can be D for DecInner or O for OpInner
  var("EXPONENTIALDECAY", EXPONENTIALDECAY, STRING, NO, "NO");
  var("TSTAR", &TSTAR, REAL, NO, "4863.5"); //In Kelvin
  var("RSTAR", &RSTAR, REAL, NO, "1.0"); //In solar radius
  var("STARIRRAD", STARIRRAD, STRING, NO, "NO");
  var("STELLARIRRADIATION", STELLARIRRADIATION, STRING, NO, "NO");
  var("EQUILIBRIUMINIT", EQUILIBRIUMINIT, STRING, NO, "NO");
  var("ZMETAL", &ZMETAL, REAL, NO, "1.0");
  var("ALPHASIGMA", ALPHASIGMA, STRING, NO, "NO");
  var("ALPHADEAD", &ALPHADEAD, REAL, NO, "0.0000001");
  var("ALPHAACTIVE", &ALPHAACTIVE, REAL, NO, "0.01");
  var("SIGMAACTIVE", &SIGMAACTIVE, REAL, NO, "100.0");
  var("WRITEONEDVISC", WRITEONEDVISC, STRING, NO, "NO");
  var("GAMMAVALUE", &GAMMAVALUE, REAL, NO, "0.0");
  var("IMPOSEDGAMMA", IMPOSEDGAMMA, STRING, NO, "NO");
  var("FORCEDCIRCULARTEMPORARY", FORCEDCIRCULARTEMPORARY, STRING, NO, "NO");
  var("SMOOTHINGAT", SMOOTHINGAT, STRING, NO, "Disc"); //If it is p (means planet), the smoothing length is calculated at the planet's location
  var("CORRECTVPLANET", CORRECTVPLANET, STRING, NO, "NO");
  var("WRITEFORCES", WRITEFORCES, STRING, NO, "NO"); //If yes, the force will be written for each planet in a file named forces.
  var("DISKINNEREDGE", &DISKINNEREDGE, REAL, NO, "0.0"); //In radii smaller than this, the surface density is progressively decreased
  var("MULTIMASSTAPER", MULTIMASSTAPER, STRING, NO, "NO");
  var("OPACITYTYPE", OPACITYTYPE, STRING, NO, "NO"); //The type of opacity function: u for unsmoothed Lin&Bell, p for smoothed opacity as alex has in pluto, and else for smoothed opacity as Bert has
  var("PABLOFORCE", PABLOFORCE, STRING, NO, "NO"); // if yes, the torque is calculated using the perturbed surface density
  /* Inned edge model of Willy Kley that has a smoother transtion of sigma and viscsity 
     The Cavity model in Fargo has an abrupt transtion when connected to the background profiles. 
     It might helps producing voctices */
  var("EDGETRANSITION", EDGETRANSITION, STRING, NO, "NO"); // if yes, a transtion zone similar to a cavity is created but this one does not the discontinuity of the Masset's cavit model.
  var("EDGESIGMADROP", &EDGESIGMADROP, REAL, NO, "1.0"); //This factor determines the drop inside the inner disc edge, must be used with EDGETRANSITION YES
  var("EDGERMID", &EDGERMID, REAL, NO, "1.0"); //Middle of the transition for the  disc edge model, must be used with EDGETRANSITION YES
  var("EDGEDELTA", &EDGEDELTA, REAL, NO, "0.0"); //Width of the transition for the  disc edge model, must be used with EDGETRANSITION YES
}
