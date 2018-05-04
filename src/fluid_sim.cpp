#include "fluid_sim.h"

double rand(double LO, double HI) {
	return LO + (rand()) / ((RAND_MAX / (HI - LO)));
}

FluidSim::FluidSim()
{
	mGrid.theRenderMode = MACGrid::PARTICLES;	
	mFrameNum = 0;
	mTotalFrameNum = 0;
	mRecordEnabled = false;

	reset();
}

void FluidSim::reset()
{
	mGrid.reset();
	mGrid.particles.clear();
	mGrid.particlesCopy.clear();
	mGrid.particlesCopyPIC.clear();
    initParticles(mGrid);
	mTotalFrameNum = 0;
}

void FluidSim::initParticles(MACGrid &mGrid)
{
    // Bridson recommends 8 particles per cell
    int seed = 10;   
    for(int k = 2; k < theContainer[2]; k++) {
        for(int j = theContainer[1] - 5; j < theContainer[1]; j++) {
            for(int i = 2; i < theContainer[0]; i++) {
                for(int iter = 0; iter < seed; iter++) {
                    vec3 pos = vec3(rand(i, i + 1), rand(j, j + 1), rand(k, k+ 1)) * theCellSize;
                    vec3 vel = vec3(0,0,0);

                    Particle p = Particle(pos, vel);
                    mGrid.particles.push_back(p);
                }
            }
        }
    }
}

void FluidSim::initMarkerGrid(MACGrid &mGrid)
{
	// std::cout << "INIT MARKERS" << std::endl;
    int boundX = theDim[MACGrid::X];
    int boundY = theDim[MACGrid::Y];
    int boundZ = theDim[MACGrid::Z];

	mGrid.markerGrid.initialize(AIR);
    
    FOR_EACH_CELL
    {
        if( i == 0 || i + 1 == boundX || 
            j == 0 || j + 1 == boundY ||
            k == 0 || k + 1 == boundZ) 
        {
            mGrid.markerGrid(i, j, k) = SOLID;
        }
		else {
			mGrid.markerGrid(i, j, k) = AIR;
		}
    }

	for(int i = 0; i < mGrid.particles.size(); i++) {
		vec3 index = getParticleIndex(mGrid.particles.at(i).position);
		mGrid.markerGrid(index[0], index[1], index[2]) = FLUID;
	}
}

double FluidSim::kernelHelper(const double &r)
{
    if(abs(r) <= 1) {
        return 1 - abs(r);
    }
    else {
        std::cout << 0 << std::endl;
        return 0;
    }
}

void FluidSim::averageVelocities(GridData &velNeighbor, GridData &weight, const vec3 &pos, const double &vel, const vec3 &index)
{
	int i = index[0];
    int j = index[1];
    int k = index[2];

    double x = pos[0];
    double y = pos[1];
    double z = pos[2];

    double kernel, rx, ry, rz;

    // Current cell
    rx = x - ((i) * theCellSize);
    ry = y - ((j) * theCellSize);
    rz = z - ((k) * theCellSize);
    kernel = kernelHelper(rx) * kernelHelper(ry) * kernelHelper(rz);
    velNeighbor(i, j, k) += vel * kernel;
    weight(i, j, k) += kernel;   

	// Right cell
    rx = x - ((i + 1) * theCellSize);
    ry = y - ((j) * theCellSize);
    rz = z - ((k) * theCellSize);
    kernel = kernelHelper(rx) * kernelHelper(ry) * kernelHelper(rz);
    velNeighbor(i + 1, j, k) += vel * kernel;
    weight(i + 1, j, k) += kernel;       

    // Top right cell
    rx = x - ((i + 1) * theCellSize);
    ry = y - ((j + 1) * theCellSize);
    rz = z - ((k) * theCellSize);
    kernel = kernelHelper(rx) * kernelHelper(ry) * kernelHelper(rz);
    velNeighbor(i + 1, j + 1, k) += vel * kernel;
    weight(i + 1, j + 1, k) += kernel;    

    // Top cell
    rx = x - ((i) * theCellSize);
    ry = y - ((j + 1) * theCellSize);
    rz = z - ((k) * theCellSize);
    kernel = kernelHelper(rx) * kernelHelper(ry) * kernelHelper(rz);
    velNeighbor(i, j + 1, k) += vel * kernel;
    weight(i, j + 1, k) += kernel;     

    // Top front cell
    rx = x - ((i) * theCellSize);
    ry = y - ((j + 1) * theCellSize);
    rz = z - ((k + 1) * theCellSize);
    kernel = kernelHelper(rx) * kernelHelper(ry) * kernelHelper(rz);
    velNeighbor(i, j + 1, k + 1) += vel * kernel;
    weight(i, j + 1, k + 1) += kernel;     

    // Front cell
    rx = x - ((i) * theCellSize);
    ry = y - ((j) * theCellSize);
    rz = z - ((k + 1) * theCellSize);
    kernel = kernelHelper(rx) * kernelHelper(ry) * kernelHelper(rz);
    velNeighbor(i, j, k + 1) += vel * kernel;
    weight(i, j, k + 1) += kernel;   

    // Front right cell
    rx = x - ((i + 1) * theCellSize);
    ry = y - ((j) * theCellSize);
    rz = z - ((k + 1) * theCellSize);
    kernel = kernelHelper(rx) * kernelHelper(ry) * kernelHelper(rz);
    velNeighbor(i + 1, j, k + 1) += vel * kernel;
    weight(i + 1, j, k + 1) += kernel;  

    // Top front right cell
    rx = x - ((i + 1) * theCellSize);
    ry = y - ((j + 1) * theCellSize);
    rz = z - ((k + 1) * theCellSize);
    kernel = kernelHelper(rx) * kernelHelper(ry) * kernelHelper(rz);
    velNeighbor(i + 1, j + 1, k + 1) += vel * kernel;
    weight(i + 1, j + 1, k + 1) += kernel;  

	// Left cell
    rx = x - ((i - 1) * theCellSize);
    ry = y - ((j) * theCellSize);
    rz = z - ((k) * theCellSize);
    kernel = kernelHelper(rx) * kernelHelper(ry) * kernelHelper(rz);
    velNeighbor(i - 1, j, k) += vel * kernel;
    weight(i + 1, j, k) += kernel;       

    // Bottom left cell
    rx = x - ((i - 1) * theCellSize);
    ry = y - ((j - 1) * theCellSize);
    rz = z - ((k) * theCellSize);
    kernel = kernelHelper(rx) * kernelHelper(ry) * kernelHelper(rz);
    velNeighbor(i - 1, j - 1, k) += vel * kernel;
    weight(i + 1, j + 1, k) += kernel;    

    // Bottom cell
    rx = x - ((i) * theCellSize);
    ry = y - ((j - 1) * theCellSize);
    rz = z - ((k) * theCellSize);
    kernel = kernelHelper(rx) * kernelHelper(ry) * kernelHelper(rz);
    velNeighbor(i, j - 1, k) += vel * kernel;
    weight(i, j - 1, k) += kernel;     

    // Bottom back cell
    rx = x - ((i) * theCellSize);
    ry = y - ((j - 1) * theCellSize);
    rz = z - ((k - 1) * theCellSize);
    kernel = kernelHelper(rx) * kernelHelper(ry) * kernelHelper(rz);
    velNeighbor(i, j - 1, k - 1) += vel * kernel;
    weight(i, j - 1, k - 1) += kernel;     

    // Back cell
    rx = x - ((i) * theCellSize);
    ry = y - ((j) * theCellSize);
    rz = z - ((k - 1) * theCellSize);
    kernel = kernelHelper(rx) * kernelHelper(ry) * kernelHelper(rz);
    velNeighbor(i, j, k - 1) += vel * kernel;
    weight(i, j, k - 1) += kernel;   

    // Back left cell
    rx = x - ((i - 1) * theCellSize);
    ry = y - ((j) * theCellSize);
    rz = z - ((k - 1) * theCellSize);
    kernel = kernelHelper(rx) * kernelHelper(ry) * kernelHelper(rz);
    velNeighbor(i - 1, j, k - 1) += vel * kernel;
    weight(i - 1, j, k - 1) += kernel;  

    // Bottom back left cell
    rx = x - ((i - 1) * theCellSize);
    ry = y - ((j - 1) * theCellSize);
    rz = z - ((k - 1) * theCellSize);
    kernel = kernelHelper(rx) * kernelHelper(ry) * kernelHelper(rz);
    velNeighbor(i - 1, j - 1, k - 1) += vel * kernel;
    weight(i - 1, j - 1, k - 1) += kernel;  
}

vec3 FluidSim::getParticleIndex(const vec3 &pos)
{
    int x = floor(pos[0] / theCellSize);
    int y = floor(pos[1] / theCellSize);
    int z = floor(pos[2] / theCellSize);

    return vec3(x, y, z);
}

void FluidSim::particleToGrid(MACGrid &mGrid, double dt)
{
	GridDataX uNeighbors;
    GridDataY vNeighbors;
    GridDataZ wNeighbors;
    GridDataX uWeight;
    GridDataY vWeight;
    GridDataZ wWeight;

    uNeighbors.initialize(0);
    vNeighbors.initialize(0);
    wNeighbors.initialize(0);
    uWeight.initialize(0);
    vWeight.initialize(0);
    wWeight.initialize(0);

    double rx, ry, rz;
    for(int i = 0; i < mGrid.particles.size(); i++)  
    {
        // Get the particle's position in the velocity components
        // Calculate kernel weight

        Particle p = mGrid.particles.at(i);
        vec3 pPos = p.position;
        vec3 pVel = p.velocity;
        vec3 pIdx = getParticleIndex(pPos);

        // Update marker grid just in case initMarkerGrid() missed it
       	mGrid.markerGrid(pIdx[0], pIdx[1], pIdx[2]) = FLUID;

        // Average x components
        // Get position and index with respect to the mU face
        vec3 pPosX = mGrid.mU.worldToSelf(pPos);
        vec3 pIdxX = getParticleIndex(pPosX);
        double velX = pVel[0];
        averageVelocities(uNeighbors, uWeight, pPosX, velX, pIdxX); 

        // Average y components
        // Get position and index with respect to the mV face
        vec3 pPosY = mGrid.mV.worldToSelf(pPos);
        vec3 pIdxY = getParticleIndex(pPosY);
        double velY = pVel[1];
        averageVelocities(vNeighbors, vWeight, pPosY, velY, pIdxY);

        // Average z components
        // Get position and index with respect to the mW face
        vec3 pPosZ = mGrid.mW.worldToSelf(pPos);
        vec3 pIdxZ = getParticleIndex(pPosZ);
        double velZ = pVel[2];
        averageVelocities(wNeighbors, wWeight, pPosZ, velZ, pIdxZ);
    }

    FOR_EACH_FACE
    {
        if(mGrid.isValidFace(MACGrid::X, i, j, k)) {
            if(uWeight(i, j, k)) {
                mGrid.mU(i, j, k) = uNeighbors(i, j, k) / uWeight(i, j, k);
            }
        }

        if(mGrid.isValidFace(MACGrid::Y, i, j, k)) {
            if(vWeight(i, j, k)) {
                mGrid.mV(i, j, k) = vNeighbors(i, j, k) / vWeight(i, j, k);
            }
        }
        
        if(mGrid.isValidFace(MACGrid::Z, i, j, k)) {
            if(wWeight(i, j, k)) {
                mGrid.mW(i, j, k) = wNeighbors(i, j, k) / wWeight(i, j, k);
            }
        }
    }
}

void FluidSim::saveGridVelFLIP(MACGrid &mGrid)
{
	mGrid.mUcopy = mGrid.mU;
    mGrid.mVcopy = mGrid.mV;
    mGrid.mWcopy = mGrid.mW;
}

void FluidSim::computeGravity(MACGrid &mGrid, double dt)
{
    vec3 gravity = vec3(0, -60, 0) * dt;
    FOR_EACH_FACE
    {
        if(mGrid.isValidFace(MACGrid::Y, i, j, k)) {
            mGrid.mV(i, j, k) += gravity[1];
        }
    }
}

void FluidSim::setBoundsToZero(MACGrid &mGrid) 
{
	FOR_EACH_CELL
    {
        if(i == 0) {
            mGrid.mU(0, j, k) = 0;
            mGrid.mU(1, j, k) = 0;
        }
        else if(i + 1 == theDim[MACGrid::X]) {
            mGrid.mU(theDim[MACGrid::X], j, k) = 0;
            mGrid.mU(theDim[MACGrid::X] - 1, j, k) = 0;
        }

        if(j == 0) {
            mGrid.mV(i, 0, k) = 0;
            mGrid.mV(i, 1, k) = 0;
        }
        else if(j + 1 == theDim[MACGrid::Y]) {
            mGrid.mV(i, theDim[MACGrid::Y], k) = 0;
            mGrid.mV(i, theDim[MACGrid::Y] - 1, k) = 0;
        }

        if(k == 0) {
            mGrid.mW(i, j, 0) = 0;
            mGrid.mW(k, j, 1) = 0;
        }
        else if(k + 1 == theDim[MACGrid::Z]) {
            mGrid.mW(i, j, theDim[MACGrid::Z]) = 0;
            mGrid.mW(i, j, theDim[MACGrid::Z] - 1) = 0;
        }
    }
}

void FluidSim::addExternalForces(MACGrid &mGrid, double dt)
{
	computeGravity(mGrid, dt);
}

void FluidSim::computeDivergence(MACGrid &mGrid, Eigen::VectorXd &rhs)
{
	double scale = 1 / theCellSize;
    double divergence = 0;

    FOR_EACH_CELL
    {
        if(mGrid.markerGrid(i, j, k) == FLUID) {
            double uPlus = mGrid.mU(i + 1, j, k);
            double uMinus = mGrid.mU(i, j, k);
            double vPlus = mGrid.mV(i, j + 1, k);
            double vMinus = mGrid.mV(i, j, k);
            double wPlus = mGrid.mW(i, j, k + 1);
            double wMinus = mGrid.mW(i, j, k);

            // Divergence estimates the rate that fluid coming in and out of a cell.
            // The velocity components at the edges should be 0.

            if(mGrid.markerGrid(i - 1, j, k) == SOLID) {
                uMinus = 0;
            }
            if(mGrid.markerGrid(i + 1, j, k) == SOLID) {
                uPlus = 0;
            }

            if(mGrid.markerGrid(i, j - 1, k) == SOLID) {
                vMinus = 0;
            }
            if(mGrid.markerGrid(i, j + 1, k) == SOLID) {
                vPlus = 0;
            }

            if(mGrid.markerGrid(i, j, k - 1) == SOLID) {
                wMinus = 0;
            }
            if(mGrid.markerGrid(i, j, k + 1) == SOLID) {
                wPlus = 0;
            }

            // RHS of 4.22
            // In the Fluid Simulation BOOK, Bridson says that we're only interested in the negative divergence.
            int id = mGrid.getID(i, j, k);
            rhs[id] = -scale * ((uPlus - uMinus) + (vPlus - vMinus) + (wPlus - wMinus));
        }
    }
}

void FluidSim::createAMatrix(MACGrid &mGrid, std::vector<Eigen::Triplet<double>> &coefficients, long n)
{
	FOR_EACH_CELL
    {
        int id = mGrid.getID(i, j, k);
        double scale = 1 / (theCellSize * theCellSize);
        double ADiag = 0;

        if(mGrid.markerGrid(i, j, k) == FLUID) {
            if(i > 0 && mGrid.markerGrid(i - 1, j, k) == FLUID) {
                ADiag += scale;
                coefficients.push_back(Eigen::Triplet<double>(id, mGrid.getID(i - 1, j, k), -scale));
            }            
            if(i < n && mGrid.markerGrid(i + 1, j, k) == FLUID) {
                ADiag += scale;
                coefficients.push_back(Eigen::Triplet<double>(id, mGrid.getID(i + 1, j, k), -scale));
            }
            if(i > 0 && mGrid.markerGrid(i - 1, j, k) == AIR) {
                ADiag += scale;
            }
            if(i < n && mGrid.markerGrid(i + 1, j, k) == AIR) {
                ADiag += scale;
            }

            if(j > 0 && mGrid.markerGrid(i, j - 1, k) == FLUID) {
                ADiag += scale;
                coefficients.push_back(Eigen::Triplet<double>(id, mGrid.getID(i, j - 1, k), -scale));
            }            
            if(j < n && mGrid.markerGrid(i, j + 1, k) == FLUID) {
                ADiag += scale;
                coefficients.push_back(Eigen::Triplet<double>(id, mGrid.getID(i, j + 1, k), -scale));
            }
            if(j > 0 && mGrid.markerGrid(i, j - 1, k) == AIR) {
                ADiag += scale;
            }
            if(j < n && mGrid.markerGrid(i, j + 1, k) == AIR) {
                ADiag += scale;
            }

            if(k > 0 && mGrid.markerGrid(i, j, k - 1) == FLUID) {
                ADiag += scale;
                coefficients.push_back(Eigen::Triplet<double>(id, mGrid.getID(i, j, k - 1), -scale));
            }            
            if(k < n && mGrid.markerGrid(i, j, k + 1) == FLUID) {
                ADiag += scale;
                coefficients.push_back(Eigen::Triplet<double>(id, mGrid.getID(i, j, k + 1), -scale));
            }
            if(k > 0 && mGrid.markerGrid(i, j, k - 1) == AIR) {
                ADiag += scale;
            }
            if(k < n && mGrid.markerGrid(i, j, k + 1) == AIR) {
                ADiag += scale;
            }
        }

        coefficients.push_back(Eigen::Triplet<double>(id, mGrid.getID(i, j, k), ADiag));
    }	
}

void FluidSim::fillPressure(MACGrid &mGrid, Eigen::VectorXd x)
{
    int id = 0;
    FOR_EACH_CELL
    {
        if(mGrid.markerGrid(i, j, k) == FLUID) {
            id = mGrid.getID(i, j, k);

            if(std::isnan(x[id])) {
                std::cout << "PCG didn't converge" << std::endl;
            }

            mGrid.mP(i, j, k) = (double)x[id];
        }
    }    
}

void FluidSim::subtractPressure(MACGrid &mGrid)
{
    double scale = 1 / theCellSize;

    int x = theDim[0];
    int y = theDim[1];
    int z = theDim[2];

    FOR_EACH_CELL
    {
        if(mGrid.markerGrid(i - 1, j, k) == FLUID || mGrid.markerGrid(i, j, k) == FLUID) {
            if(mGrid.markerGrid(i - 1, j, k) == SOLID || mGrid.markerGrid(i, j, k) == SOLID) {
                mGrid.mU(i, j, k) = 0;
            }
            else {
                mGrid.mU(i, j, k) += -(scale * (mGrid.mP(i, j, k) - mGrid.mP(i - 1, j, k)));
            }
        }

        if(mGrid.markerGrid(i, j - 1, k) == FLUID || mGrid.markerGrid(i, j, k) == FLUID) {
            if(mGrid.markerGrid(i, j - 1, k) == SOLID || mGrid.markerGrid(i, j, k) == SOLID) {
                mGrid.mV(i, j, k) = 0;
            }
            else {
                mGrid.mV(i, j, k) += -(scale * (mGrid.mP(i, j, k) - mGrid.mP(i, j - 1, k)));
            }
        }

        if(mGrid.markerGrid(i, j, k - 1) == FLUID || mGrid.markerGrid(i, j, k) == FLUID) {
            if(mGrid.markerGrid(i, j, k - 1) == SOLID || mGrid.markerGrid(i, j, k) == SOLID) {
                mGrid.mW(i, j, k) = 0;
            }
            else {
                mGrid.mW(i, j, k) += -(scale * (mGrid.mP(i, j, k) - mGrid.mP(i, j, k - 1)));
            }
        }
    }
}


void FluidSim::projectPressure(MACGrid &mGrid, double dt)
{
	long m = theDim[0] * theDim[1] * theDim[2];

	// Compute divergence
	Eigen::VectorXd rhs(m);
	rhs.setZero();
	computeDivergence(mGrid, rhs);

	// Create A Matrix
	std::vector<Eigen::Triplet<double>> coeffs;
	createAMatrix(mGrid, coeffs, m);
	Eigen::SparseMatrix<double> A(m, m);
	A.setZero();
	A.setFromTriplets(coeffs.begin(), coeffs.end());

	// PCG
	Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper, Eigen::IncompleteCholesky<double>> pcg(A);

	// Ap = d
	Eigen::VectorXd p(m);
	p.setZero();
	p = pcg.solve(rhs);

	fillPressure(mGrid,p);
	subtractPressure(mGrid);	
}

void FluidSim::updateVelFLIP(MACGrid &mGrid)
{
	FOR_EACH_FACE
	{
		if(mGrid.isValidFace(MACGrid::X, i, j, k)) {
			mGrid.mUcopy(i, j, k) = mGrid.mU(i, j, k) - mGrid.mUcopy(i, j, k);			
		}
		if(mGrid.isValidFace(MACGrid::Y, i, j, k)) {
			mGrid.mVcopy(i, j, k) = mGrid.mV(i, j, k) - mGrid.mVcopy(i, j, k);
		}
		if(mGrid.isValidFace(MACGrid::Z, i, j, k)) {
			mGrid.mWcopy(i, j, k) = mGrid.mW(i, j, k) - mGrid.mWcopy(i, j, k);
		}
	}

    mGrid.particlesCopy = mGrid.particles;
    for(int i = 0; i < mGrid.particlesCopy.size(); i++) {
        vec3 pos = mGrid.particlesCopy.at(i).position;
        mGrid.particlesCopy.at(i).velocity[0] += mGrid.mUcopy.interpolate(pos);
        mGrid.particlesCopy.at(i).velocity[1] += mGrid.mVcopy.interpolate(pos);
        mGrid.particlesCopy.at(i).velocity[2] += mGrid.mWcopy.interpolate(pos);
    }
}

void FluidSim::gridToParticle(MACGrid &mGrid)
{
	mGrid.particlesCopyPIC = mGrid.particles;
    for(int i = 0; i < mGrid.particlesCopyPIC.size(); i++) {
        vec3 pos = mGrid.particlesCopyPIC.at(i).position;
        mGrid.particlesCopyPIC.at(i).velocity[0] += mGrid.getVelocityX(pos);
        mGrid.particlesCopyPIC.at(i).velocity[1] += mGrid.getVelocityY(pos);
        mGrid.particlesCopyPIC.at(i).velocity[2] += mGrid.getVelocityZ(pos);
    }
}

void FluidSim::advectParticle(MACGrid &mGrid, double dt)
{
    double flipPercentage = .95;
    for(int i = 0; i < mGrid.particles.size(); i++)
    {
        mGrid.particles.at(i).velocity = (1 - flipPercentage) * mGrid.particlesCopyPIC.at(i).velocity + (flipPercentage * mGrid.particlesCopy.at(i).velocity);
        mGrid.particles.at(i).position = mGrid.integrate(mGrid.particles.at(i).position, mGrid.particles.at(i).velocity, dt);
    }
}

void FluidSim::resetGrid(MACGrid &mGrid)
{
	FOR_EACH_FACE
	{
		if(mGrid.isValidFace(MACGrid::X, i, j, k)) {
			mGrid.mU(i, j, k) = 0;
			mGrid.mUcopy(i, j, k) = 0;
		}
		if(mGrid.isValidFace(MACGrid::Y, i, j, k)) {
			mGrid.mV(i, j, k) = 0;
			mGrid.mVcopy(i, j, k) = 0;
		}
		if(mGrid.isValidFace(MACGrid::Z, i, j, k)) {
			mGrid.mW(i, j, k) = 0;
			mGrid.mWcopy(i, j, k) = 0;
		}
	}

	FOR_EACH_CELL
	{
		mGrid.mP(i, j, k) = 0;
		mGrid.markerGrid(i, j, k) = 0;
	}

	mGrid.particlesCopy.clear();
	mGrid.particlesCopyPIC.clear();
}

void FluidSim::step()
{
	double dt = 0.04;//0.1;

	// PIC/FLIP
	// Bridson's Sand algorithm
	initMarkerGrid(mGrid);
	particleToGrid(mGrid, dt);
	saveGridVelFLIP(mGrid);
	addExternalForces(mGrid, dt);
	setBoundsToZero(mGrid);
	projectPressure(mGrid, dt);
	setBoundsToZero(mGrid);
	updateVelFLIP(mGrid);
	gridToParticle(mGrid);
	advectParticle(mGrid, dt);
	resetGrid(mGrid);

	mTotalFrameNum++;
}