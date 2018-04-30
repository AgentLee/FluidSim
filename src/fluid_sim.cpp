#include "fluid_sim.h"

FluidSim::FluidSim()
{
    std::cout << "PIC/FLIP" << std::endl;
	mFrameNum = 0;
	mTotalFrameNum = 0;
	mRecordEnabled = false;

	reset();
}

void FluidSim::reset()
{
	mGrid.reset();
    mGrid.initParticles();
	mTotalFrameNum = 0;
}

void FluidSim::step()
{
	double dt = 0.04;//0.1;

	// PIC/FLIP
	// Bridson's Sand algorithm
	// mGrid.particleToGrid(dt);
	// mGrid.saveGridVelFLIP(dt);
	// mGrid.addExternalForces(dt);
	// mGrid.project(dt);
	// mGrid.updateVelFLIP(dt);
	// mGrid.gridToParticle(dt);
	// mGrid.advectParticle(dt);

	mTotalFrameNum++;
}