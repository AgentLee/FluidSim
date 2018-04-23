#include "constants.h"

const int theMillisecondsPerFrame = 10;

// #ifdef _DEBUG
// const int theDim[3] = {4, 4, 1};
// #else
// const int theDim[3] = {32, 32, 1};
// // const int theDim[3] = {4,4,4};
// #endif

// #define SMOKE_SIM
#ifdef SMOKE_SIM
const int theDim[3] = {32, 32, 1};
#else
const int theDim[3] = {5, 5, 5};
#endif

const double theCellSize = 0.5;

const double theAirDensity = 1.0;

const double theBuoyancyAlpha = 0.08; // Gravity's effect on the smoke particles.
const double theBuoyancyBeta = 0.37; // Buoyancy's effect due to temperature difference.	
const double theBuoyancyAmbientTemperature = 0.0; // Ambient temperature.

const double theVorticityEpsilon = 0.10;

const double theGravity = -9.8;