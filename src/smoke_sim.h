#ifndef smokeSim_H_
#define smokeSim_H_

#include "simulator.h"

class Camera;
class SmokeSim : public Simulator
{
public:
    SmokeSim();
    virtual ~SmokeSim() {}

    int scene;

    void reset();
    void step();

    void updateSources(MACGrid &mGrid, int scene);
    void advectVelocity(MACGrid &mGrid, double dt);
    // External forces
    void addExternalForces(MACGrid &mGrid, double dt);
    void computeBouyancy(MACGrid &mGrid, double dt);
    void computeVorticityConfinement(MACGrid &mGrid, double dt);
    void applyVorticityConfinement(MACGrid &mGrid, vec3 &fConf, int &i, int &j, int &k);
    // Pressure projection
    void project(MACGrid &mGrid, double dt);    
    void calculateAMatrix(MACGrid &mGrid);
    void computeDivergence(MACGrid &mGrid, GridData &d);
    bool preconditionedConjugateGradient(const GridDataMatrix & A, GridData & p, const GridData & d, int maxIterations, double tolerance);
    void apply(MACGrid &mGrid, const GridDataMatrix & matrix, const GridData & vector, GridData & result);
    void calculatePreconditioner(MACGrid &mGrid, const GridDataMatrix & A);
    void applyPreconditioner(MACGrid &mGrid, const GridData & r, const GridDataMatrix & A, GridData & z);
    // Advection
    void advectTemperature(MACGrid &mGrid, double dt);    
    void advectDensity(MACGrid &mGrid, double dt);    
    void advectRenderingParticles(MACGrid &mGrid, double dt);    
};

#endif