#ifndef smokeSim_H_
#define smokeSim_H_

#include "simulator.h"

class Camera;
class SmokeSim : public Simulator
{
public:
    SmokeSim();
    virtual ~SmokeSim() {}

    void reset();
    void step();

    void updateSources(MACGrid &mGrid);
    void advectVelocity(MACGrid &mGrid, double dt);
    void addExternalForces(MACGrid &mGrid, double dt);
    void project(MACGrid &mGrid, double dt);    
    void advectTemperature(MACGrid &mGrid, double dt);    
    void advectDensity(MACGrid &mGrid, double dt);    
    void advectRenderingParticles(MACGrid &mGrid, double dt);    
};

#endif