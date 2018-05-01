#ifndef CONSTANTS_H
#define CONSTANTS_H

#define LERP(a,b,t) (1-t)*a + t*b

// Don't modify the values of these here.
// Modify the values of these in Constants.cpp instead.
extern const int theMillisecondsPerFrame;
extern const int theDim[3];
extern const int theContainer[3];
extern const double theCellSize;
extern const double theAirDensity;
extern const double theBuoyancyAlpha;
extern const double theBuoyancyBeta;	
extern const double theBuoyancyAmbientTemperature;
extern const double theVorticityEpsilon;
extern const double theGravity;

#define FOR_EACH_CELL \
   for(int k = 0; k < theDim[MACGrid::Z]; ++k)  \
      for(int j = 0; j < theDim[MACGrid::Y]; ++j) \
         for(int i = 0; i < theDim[MACGrid::X]; ++i) 

#define FOR_EACH_CELL_REVERSE \
   for(int k = theDim[MACGrid::Z] - 1; k >= 0; --k)  \
      for(int j = theDim[MACGrid::Y] - 1; j >= 0; --j) \
         for(int i = theDim[MACGrid::X] - 1; i >= 0; --i) 

#define FOR_EACH_FACE \
   for(int k = 0; k < theDim[MACGrid::Z]+1; ++k) \
      for(int j = 0; j < theDim[MACGrid::Y]+1; ++j) \
         for(int i = 0; i < theDim[MACGrid::X]+1; ++i) 

#define FOR_EACH_FACEX \
   for(int k = 0; k < theDim[MACGrid::Z]; ++k) \
      for(int j = 0; j < theDim[MACGrid::Y]; ++j) \
         for(int i = 0; i < theDim[MACGrid::X] + 1; i) 

#define FOR_EACH_FACEY \
   for(int k = 0; k < theDim[MACGrid::Z]; ++k) \
      for(int j = 0; j < theDim[MACGrid::Y] + 1; ++j) \
         for(int i = 0; i < theDim[MACGrid::X]; ++i) 

#define FOR_EACH_FACEZ \
   for(int k = 0; k < theDim[MACGrid::Z] + 1; ++k) \
      for(int j = 0; j < theDim[MACGrid::Y]; ++j) \
         for(int i = 0; i < theDim[MACGrid::X]; ++i) 

#define FOR_EACH_PARTICLE \
    for(int i = 0; i < particles.size(); i++) \

#endif