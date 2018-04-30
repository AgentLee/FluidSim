#include "smoke_sim.h"

SmokeSim::SmokeSim()
{
	std::cout << "SMOKE" << std::endl;
	
	mFrameNum = 0;
	mTotalFrameNum = 0;
	mRecordEnabled = false;

	reset();
}

void SmokeSim::reset()
{
	mGrid.reset();
	mTotalFrameNum = 0;
}

void SmokeSim::step()
{
	double dt = 0.04;//0.1;

	// Step0: Gather user forces
	mGrid.updateSources();

	// Step1: Calculate new velocities
	mGrid.advectVelocity(dt);
	mGrid.addExternalForces(dt);
	mGrid.project(dt);

	// Step2: Calculate new temperature
	mGrid.advectTemperature(dt);

	// Step3: Calculate new density 
	mGrid.advectDensity(dt);

	// Step4: Advect rendering particles
	mGrid.advectRenderingParticles(dt);

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