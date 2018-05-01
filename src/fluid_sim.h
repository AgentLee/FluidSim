#ifndef fluidSim_H_
#define fluidSim_H_

#include "simulator.h"
#include <Eigen/Core>

class Camera;
class FluidSim : public Simulator
{
public:
    FluidSim();
    virtual ~FluidSim() {}

    void reset();
    void step();

    void initParticles(MACGrid &mGrid);
    void initMarkerGrid(MACGrid &mGrid);
    // Particle to grid
    void particleToGrid(MACGrid &mGrid, double dt);
    void averageVelocities(GridData &velNeighbor, GridData &weight, const vec3 &pos, const double &vel, const vec3 &index);
    double kernelHelper(const double &r);
    vec3 getParticleIndex(const vec3 &pos);
    // FLIP save
    void saveGridVelFLIP(MACGrid &mGrid);
    // External forces
    void addExternalForces(MACGrid &mGrid, double dt);
    void computeGravity(MACGrid &mGrid, double dt);
    // Bounds
    void setBoundsToZero(MACGrid &mGrid);
    // Pressure projection
    void projectPressure(MACGrid &mGrid, double dt);
    void computeDivergence(MACGrid &mGrid, Eigen::VectorXd &rhs);
    void createAMatrix(MACGrid &mGrid, std::vector<Eigen::Triplet<double>> &coefficients, long n);
    void fillPressure(MACGrid &mGrid, Eigen::VectorXd x);
    void subtractPressure(MACGrid &mGrid); 
    // FLIP update
    void updateVelFLIP(MACGrid &mGrid);
    // PIC update
    void gridToParticle(MACGrid &mGrid);
    // LERP PIC/FLIP
    void advectParticle(MACGrid &mGrid, double dt);
    // Clear grid
    void resetGrid(MACGrid &mGrid);
};

#endif