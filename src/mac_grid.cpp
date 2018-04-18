#include "mac_grid.h"

#undef max
#undef min 

#define SMOKE_SIM true

// Globals
MACGrid target;
enum cellType 
{
    SOLID   = 0,
    FLUID   = 1,
    AIR     = 2,
};

// NOTE: x -> cols, z -> rows, y -> stacks
MACGrid::RenderMode MACGrid::theRenderMode = SHEETS;
bool MACGrid::theDisplayVel = false;//true

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
    for(int i = 0; i < particles.size(); ++i) 
    
MACGrid::MACGrid()
{
    initialize();
}

MACGrid::MACGrid(const MACGrid& orig)
{
   mU = orig.mU;
   mV = orig.mV;
   mW = orig.mW;
   mP = orig.mP;
   mD = orig.mD;
   mT = orig.mT;
}

MACGrid& MACGrid::operator=(const MACGrid& orig)
{
   if (&orig == this)
   {
      return *this;
   }
   mU = orig.mU;
   mV = orig.mV;
   mW = orig.mW;
   mP = orig.mP;
   mD = orig.mD;
   mT = orig.mT;   

   return *this;
}

MACGrid::~MACGrid()
{
}

void MACGrid::reset()
{
    initMarkerGrid();

    mU.initialize();
    mV.initialize();
    mW.initialize();
    mP.initialize();
    mD.initialize();
    mT.initialize(0.0);

    calculateAMatrix();
    calculatePreconditioner(AMatrix);
}

void MACGrid::initialize()
{
   reset();
}

void MACGrid::updateSources()
{
    // Set initial values for density, temperature, velocity
    for(int i=6; i<12;i++){
        for(int j=0; j<5; j++){
            mV(i,j+1,0) = 2.0;
            mD(i,j,0) = 1.0;
            mT(i,j,0) = 1.0;

            mV(i,j+2,0) = 2.0;
            mD(i,j,0) = 1.0;
            mT(i,j,0) = 1.0;
        }
    }

	// Refresh particles in source.
	for(int i=6; i<12; i++) {
		for (int j = 0; j < 5; j++) {
			for (int k = 0; k <= 0; k++) {
				vec3 cell_center(theCellSize*(i+0.5), theCellSize*(j+0.5), theCellSize*(k+0.5));
				for(int p=0; p<10; p++) {
                    double a = ((float) rand() / RAND_MAX - 0.5) * theCellSize;
                    double b = ((float) rand() / RAND_MAX - 0.5) * theCellSize;
                    double c = ((float) rand() / RAND_MAX - 0.5) * theCellSize;
                    vec3 shift(a, b, c);
                    vec3 xp = cell_center + shift;
                    rendering_particles.push_back(xp);
                }
			}
		}
	}
}

void MACGrid::initParticles()
{
    // Bridson recommends 8 particles per cell
    int seed = 8;   
}

void MACGrid::initMarkerGrid()
{   
    markerGrid.initialize(AIR);
    
    int boundX = theDim[MACGrid::X];
    int boundY = theDim[MACGrid::Y];
    int boundZ = theDim[MACGrid::Z];
    FOR_EACH_CELL
    {
        // 2D
        if( i == 0 || i - 1 == boundX || 
            j == 0 || j - 1 == boundY) 
        {
            markerGrid(i, j, k) = SOLID;
        }

        // 3D 
        if(boundZ > 1) {
            if(k == 0 || k - 1 == boundZ) {
                markerGrid(i, j, k) = SOLID;
            }
        }
    }

    // #define _DEBUG
    #ifdef _DEBUG
    {
        FOR_EACH_CELL
        {
            if(markerGrid(i, j, k) == AIR) {
                std::cout << "AIR: " << i << " " << j << " " << k << std::endl;
            }

            if(markerGrid(i, j, k) == FLUID) {
                std::cout << "FLUID: " << i << " " << j << " " << k << std::endl;
            }

            if(markerGrid(i, j, k) == SOLID) {
                std::cout << "SOLID: " << i << " " << j << " " << k << std::endl;
            }
        }
    }
    #endif
}

double MACGrid::kernelHelper(const double &r)
{
    if(r >= 0 && r <= 1) {
        return 1 - r;
    }
    else if(r >= -1 && r >= 0) {
        return 1 + r;
    }
    else {
        return 0;
    }
}

double MACGrid::kernel(const int &_x, const int &_y, const int &_z) 
{
    double x = _x / theCellSize;
    double y = _y / theCellSize;
    double z = _z / theCellSize;

    double weight = kernelHelper(x) * kernelHelper(y) * kernelHelper(z);

    return weight;
}

/*
 * Calculate a weighted average of nearby particles in a cell
 * using a trilinear interpolation kernel function.
 */
void MACGrid::particleToGrid(double dt)
{
    FOR_EACH_PARTICLE 
    {
        // Get the particle's position in the velocity components
        // Calculate kernel weight
    }

    FOR_EACH_FACEX
    {
        // Divide by kernel weight
    }

    FOR_EACH_FACEY
    {
        // Divide by kernel weight
    }

    FOR_EACH_FACEZ
    {
        // Divide by kernel weight
    }
}

void MACGrid::saveGridVelFLIP(double dt)
{
    // TODO
}

void MACGrid::updateVelFLIP(double dt)
{
    // TODO
}

void MACGrid::gridToParticle(double dt) 
{
    // TODO
}

void MACGrid::advectParticle(double dt)
{
    // TODO
}

/*
 * Want to figure out new value of q at some grid point.
 * To do so we need the old value of q at the point it ends up at.
 * We know where the particle ends up. It moves through the velocity field.
 * To find where it started, we need to back trace from the grid point.
 * 
 * All we need is the value of q at this new grid point and 
 * the value of q at the old grid point
 * If the starting point isn't on the grid then we can interpolate
 * the value from the existing ones on the grid.
 * 
 * xp - where the particle started from
 * xG - where the particle ends up
 * dx/dt = u - how the particle moves
 * Forward Euler - xp = xG - dt(uG)
*/
void MACGrid::advectVelocity(double dt)
{
    // TODO: Calculate new velocities and store in target
    // TODO: Your code is here. It builds target.mU, target.mV and target.mW for all faces
    
    // Velocities are stored at the faces of the MAC grid
    FOR_EACH_FACE
    {
        if(isValidFace(MACGrid::X, i, j, k)) {
            vec3 currPos = getFacePosition(MACGrid::X, i, j, k);
            vec3 startPos = getRewoundPosition(currPos, dt);

            // q_n+1 = interpolate(q_n, xG - dt(uG));
            target.mU(i, j, k) = getVelocityX(startPos);
        }
        
        if(isValidFace(MACGrid::Y, i, j, k)) {
            vec3 currPos = getFacePosition(MACGrid::Y, i, j, k);
            vec3 startPos = getRewoundPosition(currPos, dt);

            // q_n+1 = interpolate(q_n, xG - dt(uG));
            target.mV(i, j, k) = getVelocityY(startPos);
        }

        if(isValidFace(MACGrid::Z, i, j, k)) {
            vec3 currPos = getFacePosition(MACGrid::Z, i, j, k);
            vec3 startPos = getRewoundPosition(currPos, dt);

            // q_n+1 = interpolate(q_n, xG - dt(uG));
            target.mW(i, j, k) = getVelocityZ(startPos);
        }  
    }

    // Then save the result to our object
    mU = target.mU;
    mV = target.mV;
    mW = target.mW;
}

void MACGrid::advectTemperature(double dt)
{
    // TODO: Calculate new temp and store in target
    // TODO: Your code is here. It builds target.mT for all cells.

    // Temperature is contained at the centers
    FOR_EACH_CELL
    {
        vec3 currPos = getCenter(i, j, k);
        vec3 startPos = getRewoundPosition(currPos, dt);

        target.mT(i, j, k) = getTemperature(startPos);
    }

    // Then save the result to our object
    mT = target.mT;
}


void MACGrid::advectRenderingParticles(double dt) {
	rendering_particles_vel.resize(rendering_particles.size());
	for (size_t p = 0; p < rendering_particles.size(); p++) {
		vec3 currentPosition = rendering_particles[p];
        vec3 currentVelocity = getVelocity(currentPosition);
        vec3 nextPosition = currentPosition + currentVelocity * dt;
        vec3 clippedNextPosition = clipToGrid(nextPosition, currentPosition);
        // Keep going...
        vec3 nextVelocity = getVelocity(clippedNextPosition);
        vec3 averageVelocity = (currentVelocity + nextVelocity) / 2.0;
        vec3 betterNextPosition = currentPosition + averageVelocity * dt;
        vec3 clippedBetterNextPosition = clipToGrid(betterNextPosition, currentPosition);
        rendering_particles[p] = clippedBetterNextPosition;
		rendering_particles_vel[p] = averageVelocity;
	}
}

void MACGrid::advectDensity(double dt)
{
    // TODO: Calculate new densitities and store in target
    // TODO: Your code is here. It builds target.mD for all cells.

    // Density is contained at the centers
    FOR_EACH_CELL
    {
        vec3 currPos = getCenter(i, j, k);
        vec3 startPos = getRewoundPosition(currPos, dt);

        target.mD(i, j, k) = getDensity(startPos);
    }

    // Then save the result to our object
    mD = target.mD;
}

/*
 * pg45
 */
void MACGrid::computeBouyancy(double dt)
{
	// TODO: Calculate bouyancy and store in target
    // TODO: Your code is here. It modifies target.mV for all y face velocities.

	FOR_EACH_FACE 
    {
		if (isValidFace(MACGrid::Y, i, j, k)) {
			vec3 currPos = getFacePosition(MACGrid::Y, i, j, k);

            double temp = getTemperature(currPos);
    	    double ambientTemp = 300;
            double s = getDensity(currPos);

            // Equation 5.1
            vec3 fBuoy(0, -theBuoyancyAlpha * s + theBuoyancyBeta * (temp - ambientTemp), 0);

            target.mV(i, j, k) = mV(i, j, k) + fBuoy[1];
		}
	}

    // and then save the result to our object
    mV = target.mV;
}

/*
 * This function prevents the vortices from disappearing too quickly.
 * w = cross(gradient, u) [curl]
 * A vortex is a peak in the vorticity field
 * pg46
 */ 
void MACGrid::computeVorticityConfinement(double dt)
{
    // TODO: Calculate vorticity confinement forces
    // Apply the forces to the current velocity and store the result in target
    // STARTED.
    // TODO: Your code is here. It modifies target.mU,mV,mW for all faces.
    
    // Needed to keep track of the vorticity at the centers
    GridData vorticityX;
    GridData vorticityY;
    GridData vorticityZ;
    GridData vorticityMag;
    vorticityX.initialize();
    vorticityY.initialize();
    vorticityZ.initialize();
    vorticityMag.initialize();

    // Average the velocities to the cell centers 
    FOR_EACH_CELL
    {
        // Equation 5.6
        // Use central differences to approximate vorticity 
        double x = ((mW(i, j + 1, k) - mW(i, j - 1, k)) / (2 * theCellSize)) - 
                   ((mV(i, j, k + 1) - mV(i, j, k - 1)) / (2 * theCellSize));
        double y = ((mU(i, j, k + 1) - mU(i, j, k - 1)) / (2 * theCellSize)) - 
                   ((mW(i + 1, j, k) - mW(i - 1, j, k)) / (2 * theCellSize));
        double z = ((mV(i + 1, j, k) - mV(i - 1, j, k)) / (2 * theCellSize)) - 
                   ((mU(i, j + 1, k) - mU(i, j - 1, k)) / (2 * theCellSize));
        vec3 vorticity(x, y, z);

        vorticityX(i, j, k) = x;
        vorticityY(i, j, k) = y;
        vorticityZ(i, j, k) = z;

        vorticityMag(i, j, k) = vorticity.Length();
    }

    // Need to find the force at the center of the cells
    FOR_EACH_CELL
    {
        vec3 vorticity(vorticityX(i, j, k), vorticityY(i, j, k), vorticityZ(i, j, k));

        // Equation 5.7
        // Calculate gradient |vorticity|
        double x = ((vorticityMag(i + 1, j, k) - vorticityMag(i - 1, j, k)) / (2 * theCellSize));
        double y = ((vorticityMag(i, j + 1, k) - vorticityMag(i, j - 1, k)) / (2 * theCellSize));
        double z = ((vorticityMag(i, j, k + 1) - vorticityMag(i, j, k - 1)) / (2 * theCellSize));
        vec3 gradVorticity(x, y, z);

        // Equation 5.8
        // Normalize and prevent divide by 0
        vec3 N = gradVorticity / (gradVorticity.Length() + pow(10, -20));

        // Equation 5.5
        // Take the cross product to get fConf at the grid centers
        vec3 fConf = theVorticityEpsilon * theCellSize * N.Cross(vorticity);

        applyVorticityConfinement(fConf, i, j, k);
    }

    // Then save the result to our object
    mU = target.mU;
    mV = target.mV;
    mW = target.mW;
}

void MACGrid::applyVorticityConfinement(vec3 &fConf, int &i, int &j, int &k)
{
    // pg47
    // "Take the appropriate averages to apply this 
    // to the different components of velocity on the MAC grid"
    if (isValidFace(MACGrid::X, i, j, k)) {
        target.mU(i, j, k) += fConf[0] / 2;
    } 
    if (isValidFace(MACGrid::X, i + 1, j, k)) {
        target.mU(i + 1, j, k) += fConf[0] / 2;
    }
    if (isValidFace(MACGrid::Y, i, j, k)) {
        target.mV(i, j, k) += fConf[1] / 2;
    }
    if (isValidFace(MACGrid::Y, i, j + 1, k)) {
        target.mV(i, j + 1, k) += fConf[1] / 2;
    }
    if (isValidFace(MACGrid::Z, i, j, k)) {
        target.mW(i, j, k) += fConf[2] / 2;
    }
    if (isValidFace(MACGrid::Z, i, j, k + 1)) {
        target.mW(i, j, k + 1) += fConf[2] / 2;
    }
}

void MACGrid::computeGravity(double dt)
{
    // TODO    
}

void MACGrid::addExternalForces(double dt)
{
    #if SMOKE_SIM
    computeBouyancy(dt);
    computeVorticityConfinement(dt);
    #else
    computeGravity(dt);
    #endif
}

void MACGrid::computeDivergence(GridData &d)
{
    FOR_EACH_FACE
    {
        // Should only do this for fluid cells
        if(isValidCell(i, j, k)) {
            // This isn't so much important for smoke sim. 
            if(markerGrid(i, j, k) == FLUID || SMOKE_SIM) {
                // Use finite differences to approximate the divergence
                double uPlus    = mU(i + 1, j, k);
                double uMinus   = mU(i, j, k);
                double vPlus    = mV(i, j + 1, k);
                double vMinus   = mV(i, j, k);
                double wPlus    = mW(i, j, k + 1);
                double wMinus   = mW(i, j, k);

                // Divergence estimates the rate that fluid coming in and out of a cell.
                // The velocity components at the edges should be 0.
                if(i == 0) {
                    uMinus = 0;
                }
                if(i + 1 == theDim[MACGrid::X]) {
                    uPlus = 0;
                }

                if(j == 0) {
                    vMinus = 0;
                }
                if(j + 1 == theDim[MACGrid::Y]) {
                    vPlus = 0;
                }

                if(k == 0) {
                    wMinus = 0;
                }
                if(k + 1 == theDim[MACGrid::Z]) {
                    wPlus = 0;
                }
                
                // RHS of 4.22
                // In the Fluid Simulation BOOK, Bridson says that we're only interested in the negative divergence.
                double scale = 1 / theCellSize;
                d(i, j, k) = -scale * ((uPlus - uMinus) + (vPlus - vMinus) + (wPlus - wMinus));
            }
        }
    }
}

/*
 * This function makes fluid incompressible and enforces boundary conditions.
 * Subtract off the pressure gradient.
 *      u_n+1 = u - dt(1/dens)(pressure grad)
 * Satisfies incompressibility
 *      dot(grad, u_n+1) = 0
 * Satisfies boundary conditions
 *      dot(u_n+1, n_hat) = dot(u, n_hat)
 */
void MACGrid::project(double dt)
{
    // TODO: Solve Ax = b for pressure
    // 1. Contruct b
    // 2. Construct A 
    // 3. Solve for p
    // 4. Subtract pressure from our velocity and save in target

    // TODO: Your code is here. It solves for a pressure field and modifies target.mU,mV,mW for all faces.
    GridData d;
    d.initialize();
    GridData p;
    p.initialize();

    // Assume that the pressure is 0 outside the fluid
    computeDivergence(d);

    // Solve for p
    if(!preconditionedConjugateGradient(AMatrix, p, d, 150, 0.01)) {
        std::cout << "PCG didn't converge!" << std::endl;
    }

    double pressureConstant = dt / (theAirDensity * theCellSize * theCellSize);
    FOR_EACH_CELL
    {
        // We just solved Ap = d for p.
        // Look at equation 4.22:
        // LHS is (dt / (density * cellsize * cellsize) * p
        // RHS is inv(A) * d
        // Need to divide by the constant value to truly solve for pressure
        p(i, j, k) /= pressureConstant;
        target.mP(i, j, k) = p(i, j, k);
    }

    // BOOK pg71 
    double scale = dt / (theAirDensity * theCellSize);
    FOR_EACH_FACE
    {
        if(isValidFace(MACGrid::X, i, j, k)) {
            if(i - 1 < 0 || i >= theDim[MACGrid::X]) {
                target.mU(i, j, k) = 0;
            }
            else {
                target.mU(i, j, k) -= scale * (p(i, j, k) - p(i - 1, j, k));
            }
        }

        if(isValidFace(MACGrid::Y, i, j, k)) {
            if(j - 1 < 0 || j >= theDim[MACGrid::Y]) {
                target.mV(i, j, k) = 0;
            }
            else {
                target.mV(i, j, k) -= scale * (p(i, j, k) - p(i, j - 1, k));
            }
        }

        if(isValidFace(MACGrid::Z, i, j, k)) {
            if(k - 1 < 0 || k >= theDim[MACGrid::Z]) {
                target.mW(i, j, k) = 0;
            }
            else {
                target.mW(i, j, k) -= scale * (p(i, j, k) - p(i, j, k - 1));
            }
        }
    }

    // #define _DEBUG
    #ifdef _DEBUG
    {
        // Check border velocities:
        FOR_EACH_FACE {
            if (isValidFace(MACGrid::X, i, j, k)) {

                if (i == 0) {
                    if (abs(target.mU(i,j,k)) > 0.0000001) {
                        PRINT_LINE( "LOW X:  " << target.mU(i,j,k) );
                        //target.mU(i,j,k) = 0;
                    }
                }

                if (i == theDim[MACGrid::X]) {
                    if (abs(target.mU(i,j,k)) > 0.0000001) {
                        PRINT_LINE( "HIGH X: " << target.mU(i,j,k) );
                        //target.mU(i,j,k) = 0;
                    }
                }

            }
            if (isValidFace(MACGrid::Y, i, j, k)) {
                

                if (j == 0) {
                    if (abs(target.mV(i,j,k)) > 0.0000001) {
                        PRINT_LINE( "LOW Y:  " << target.mV(i,j,k) );
                        //target.mV(i,j,k) = 0;
                    }
                }

                if (j == theDim[MACGrid::Y]) {
                    if (abs(target.mV(i,j,k)) > 0.0000001) {
                        PRINT_LINE( "HIGH Y: " << target.mV(i,j,k) );
                        //target.mV(i,j,k) = 0;
                    }
                }

            }
            if (isValidFace(MACGrid::Z, i, j, k)) {
                
                if (k == 0) {
                    if (abs(target.mW(i,j,k)) > 0.0000001) {
                        PRINT_LINE( "LOW Z:  " << target.mW(i,j,k) );
                        //target.mW(i,j,k) = 0;
                    }
                }

                if (k == theDim[MACGrid::Z]) {
                    if (abs(target.mW(i,j,k)) > 0.0000001) {
                        PRINT_LINE( "HIGH Z: " << target.mW(i,j,k) );
                        //target.mW(i,j,k) = 0;
                    }
                }
            }
        }
    }
    #endif

    // Then save the result to our object
    {
        mP = target.mP; 
        mU = target.mU;
        mV = target.mV;
        mW = target.mW;
    }

    #ifdef _DEBUG
    {
        // IMPLEMENT THIS AS A SANITY CHECK: assert (checkDivergence());
        // TODO: Fix duplicate code:
        FOR_EACH_CELL {
        // Construct the vector of divergences d:
            double velLowX = mU(i,j,k);
            double velHighX = mU(i+1,j,k);
            double velLowY = mV(i,j,k);
            double velHighY = mV(i,j+1,k);
            double velLowZ = mW(i,j,k);
            double velHighZ = mW(i,j,k+1);
            double divergence = ((velHighX - velLowX) + (velHighY - velLowY) + (velHighZ - velLowZ)) / theCellSize;
            if (abs(divergence) > 0.02 ) {
                PRINT_LINE("WARNING: Divergent! ");
                PRINT_LINE("Divergence: " << divergence);
                PRINT_LINE("Cell: " << i << ", " << j << ", " << k);
            }
        }
    }
    #endif
}

vec3 MACGrid::getVelocity(const vec3& pt)
{
   vec3 vel;
   vel[0] = getVelocityX(pt); 
   vel[1] = getVelocityY(pt); 
   vel[2] = getVelocityZ(pt); 
   return vel;
}

double MACGrid::getVelocityX(const vec3& pt)
{
   return mU.interpolate(pt);
}

double MACGrid::getVelocityY(const vec3& pt)
{
   return mV.interpolate(pt);
}

double MACGrid::getVelocityZ(const vec3& pt)
{
   return mW.interpolate(pt);
}

double MACGrid::getTemperature(const vec3& pt)
{
   return mT.interpolate(pt);
}

double MACGrid::getDensity(const vec3& pt)
{
   return mD.interpolate(pt);
}

vec3 MACGrid::getCenter(int i, int j, int k)
{
   double xstart = theCellSize/2.0;
   double ystart = theCellSize/2.0;
   double zstart = theCellSize/2.0;

   double x = xstart + i*theCellSize;
   double y = ystart + j*theCellSize;
   double z = zstart + k*theCellSize;
   return vec3(x, y, z);
}

vec3 MACGrid::getRewoundPosition(const vec3 & currentPosition, const double dt) 
{
	/*
	// EULER (RK1):
	vec3 currentVelocity = getVelocity(currentPosition);
	vec3 rewoundPosition = currentPosition - currentVelocity * dt;
	vec3 clippedRewoundPosition = clipToGrid(rewoundPosition, currentPosition);
	return clippedRewoundPosition;
	*/

	// HEUN / MODIFIED EULER (RK2):
	vec3 currentVelocity = getVelocity(currentPosition);
	vec3 rewoundPosition = currentPosition - currentVelocity * dt;
	vec3 clippedRewoundPosition = clipToGrid(rewoundPosition, currentPosition);
	// Keep going...
	vec3 rewoundVelocity = getVelocity(clippedRewoundPosition);
	vec3 averageVelocity = (currentVelocity + rewoundVelocity) / 2.0;
	vec3 betterRewoundPosition = currentPosition - averageVelocity * dt;
	vec3 clippedBetterRewoundPosition = clipToGrid(betterRewoundPosition, currentPosition);
	return clippedBetterRewoundPosition;
}

vec3 MACGrid::clipToGrid(const vec3& outsidePoint, const vec3& insidePoint) {
	/*
	// OLD:
	vec3 rewindPosition = outsidePoint;
	if (rewindPosition[0] < 0) rewindPosition[0] = 0; // TEMP!
	if (rewindPosition[1] < 0) rewindPosition[1] = 0; // TEMP!
	if (rewindPosition[2] < 0) rewindPosition[2] = 0; // TEMP!
	if (rewindPosition[0] > theDim[MACGrid::X]) rewindPosition[0] = theDim[MACGrid::X]; // TEMP!
	if (rewindPosition[1] > theDim[MACGrid::Y]) rewindPosition[1] = theDim[MACGrid::Y]; // TEMP!
	if (rewindPosition[2] > theDim[MACGrid::Z]) rewindPosition[2] = theDim[MACGrid::Z]; // TEMP!
	return rewindPosition;
	*/

	vec3 clippedPoint = outsidePoint;

	for (int i = 0; i < 3; i++) {
		if (clippedPoint[i] < 0) {
			vec3 distance = clippedPoint - insidePoint;
			double newDistanceI = 0 - insidePoint[i];
			double ratio = newDistanceI / distance[i];
			clippedPoint = insidePoint + distance * ratio;
		}
		if (clippedPoint[i] > getSize(i)) {
			vec3 distance = clippedPoint - insidePoint;
			double newDistanceI = getSize(i) - insidePoint[i];
			double ratio = newDistanceI / distance[i];
			clippedPoint = insidePoint + distance * ratio;
		}
	}

#ifdef _DEBUG
	// Make sure the point is now in the grid:
	if (clippedPoint[0] < 0 || clippedPoint[1] < 0 || clippedPoint[2] < 0 || clippedPoint[0] > getSize(0) || clippedPoint[1] > getSize(1) || clippedPoint[2] > getSize(2)) {
		PRINT_LINE("WARNING: Clipped point is outside grid!");
	}
#endif

	return clippedPoint;

}


double MACGrid::getSize(int dimension) {
	return theDim[dimension] * theCellSize;
}


int MACGrid::getCellIndex(int i, int j, int k)
{
	return i + j * theDim[MACGrid::X] + k * theDim[MACGrid::Y] * theDim[MACGrid::X];
}


int MACGrid::getNumberOfCells()
{
	return theDim[MACGrid::X] * theDim[MACGrid::Y] * theDim[MACGrid::Z];
}


bool MACGrid::isValidCell(int i, int j, int k)
{
	if (i >= theDim[MACGrid::X] || j >= theDim[MACGrid::Y] || k >= theDim[MACGrid::Z]) {
		return false;
	}

	if (i < 0 || j < 0 || k < 0) {
		return false;
	}

	return true;
}


bool MACGrid::isValidFace(int dimension, int i, int j, int k)
{
	if (dimension == 0) {
		if (i > theDim[MACGrid::X] || j >= theDim[MACGrid::Y] || k >= theDim[MACGrid::Z]) {
			return false;
		}
	} else if (dimension == 1) {
		if (i >= theDim[MACGrid::X] || j > theDim[MACGrid::Y] || k >= theDim[MACGrid::Z]) {
			return false;
		}
	} else if (dimension == 2) {
		if (i >= theDim[MACGrid::X] || j >= theDim[MACGrid::Y] || k > theDim[MACGrid::Z]) {
			return false;
		}
	}

	if (i < 0 || j < 0 || k < 0) {
		return false;
	}

	return true;
}


vec3 MACGrid::getFacePosition(int dimension, int i, int j, int k)
{
	if (dimension == 0) {
		return vec3(i * theCellSize, (j + 0.5) * theCellSize, (k + 0.5) * theCellSize);
	} else if (dimension == 1) {
		return vec3((i + 0.5) * theCellSize, j * theCellSize, (k + 0.5) * theCellSize);
	} else if (dimension == 2) {
		return vec3((i + 0.5) * theCellSize, (j + 0.5) * theCellSize, k * theCellSize);
	}

	return vec3(0,0,0); //???

}

void MACGrid::calculateAMatrix() 
{
	FOR_EACH_CELL 
    {
		int numFluidNeighbors = 0;
		if (i-1 >= 0) {
			AMatrix.plusI(i-1,j,k) = -1;
			numFluidNeighbors++;
		}
		if (i+1 < theDim[MACGrid::X]) {
			AMatrix.plusI(i,j,k) = -1;
			numFluidNeighbors++;
		}
		if (j-1 >= 0) {
			AMatrix.plusJ(i,j-1,k) = -1;
			numFluidNeighbors++;
		}
		if (j+1 < theDim[MACGrid::Y]) {
			AMatrix.plusJ(i,j,k) = -1;
			numFluidNeighbors++;
		}
		if (k-1 >= 0) {
			AMatrix.plusK(i,j,k-1) = -1;
			numFluidNeighbors++;
		}
		if (k+1 < theDim[MACGrid::Z]) {
			AMatrix.plusK(i,j,k) = -1;
			numFluidNeighbors++;
		}
		// Set the diagonal:
		AMatrix.diag(i,j,k) = numFluidNeighbors;
	}
}

bool MACGrid::preconditionedConjugateGradient(const GridDataMatrix & A, GridData & p, const GridData & d, int maxIterations, double tolerance) {
	// Solves Ap = d for p.

	FOR_EACH_CELL {
		p(i,j,k) = 0.0; // Initial guess p = 0.	
	}

	GridData r = d; // Residual vector.

	GridData z; z.initialize();
	applyPreconditioner(r, A, z); // Auxillary vector.

	GridData s = z; // Search vector;

	double sigma = dotProduct(z, r);

	for (int iteration = 0; iteration < maxIterations; iteration++) {
		double rho = sigma; // According to TA. Here???

		apply(A, s, z); // z = applyA(s);

		double alpha = rho/dotProduct(z, s);

		GridData alphaTimesS; alphaTimesS.initialize();
		multiply(alpha, s, alphaTimesS);
		add(p, alphaTimesS, p);
		//p += alpha * s;

		GridData alphaTimesZ; alphaTimesZ.initialize();
		multiply(alpha, z, alphaTimesZ);
		subtract(r, alphaTimesZ, r);
		//r -= alpha * z;

		if (maxMagnitude(r) <= tolerance) {
			//PRINT_LINE("PCG converged in " << (iteration + 1) << " iterations.");
			return true; //return p;
		}

		applyPreconditioner(r, A, z); // z = applyPreconditioner(r);

		double sigmaNew = dotProduct(z, r);
		double beta = sigmaNew / rho;

		GridData betaTimesS; betaTimesS.initialize();
		multiply(beta, s, betaTimesS);
		add(z, betaTimesS, s);
		//s = z + beta * s;

		sigma = sigmaNew;
	}

	PRINT_LINE( "PCG didn't converge!" );
	return false;
}

/*
 * Modified Incomplete Cholesky: pg36 course notes
 */ 
void MACGrid::calculatePreconditioner(const GridDataMatrix & A) {

	precon.initialize();

    // TODO: Build the modified incomplete Cholesky preconditioner following Fig 4.2 on page 36 of Bridson's 2007 SIGGRAPH fluid course notes.
    //       This corresponds to filling in precon(i,j,k) for all cells

    // Tuning constant
    double t = 0.97;

    FOR_EACH_CELL
    {
        // if cell is a fluid 
        {
            double e = A.diag(i, j, k)  - ((A.plusI(i - 1, j, k) * precon(i - 1, j, k)) * (A.plusI(i - 1, j, k) * precon(i - 1, j, k)))
                                        - ((A.plusJ(i, j - 1, k) * precon(i, j - 1, k)) * (A.plusJ(i, j - 1, k) * precon(i, j - 1, k)))
                                        - ((A.plusK(i, j, k - 1) * precon(i, j, k - 1)) * (A.plusK(i, j, k - 1) * precon(i, j, k - 1)))
                                        - t * ( A.plusI(i - 1, j, k) * (A.plusJ(i - 1, j, k) + A.plusK(i - 1, j, k)) * (precon(i - 1, j, k) * precon(i - 1, j, k)) +
                                                A.plusJ(i, j - 1, k) * (A.plusI(i, j - 1, k) + A.plusK(i, j - 1, k)) * (precon(i, j - 1, k) * precon(i, j - 1, k)) +
                                                A.plusK(i, j, k - 1) * (A.plusI(i, j, k - 1) + A.plusJ(i, j, k - 1)) * (precon(i, j, k - 1) * precon(i, j, k - 1)));

            // Set to a small number to prevent divide by 0
            precon(i, j, k) = 1 / sqrt(e + pow(10, -30));
        }
    }
}

/*
 * pg37 in course notes
 */
void MACGrid::applyPreconditioner(const GridData & r, const GridDataMatrix & A, GridData & z) {

    // TODO: change if(0) to if(1) after you implement calculatePreconditoner function.

    if(1) {

        // APPLY THE PRECONDITIONER:
        // Solve Lq = r for q:
        GridData q;
        q.initialize();
        FOR_EACH_CELL {
                    //if (A.diag(i,j,k) != 0.0) { // If cell is a fluid.
                    double t = r(i, j, k) - A.plusI(i - 1, j, k) * precon(i - 1, j, k) * q(i - 1, j, k)
                               - A.plusJ(i, j - 1, k) * precon(i, j - 1, k) * q(i, j - 1, k)
                               - A.plusK(i, j, k - 1) * precon(i, j, k - 1) * q(i, j, k - 1);
                    q(i, j, k) = t * precon(i, j, k);
                    //}
                }
        // Solve L^Tz = q for z:
        FOR_EACH_CELL_REVERSE {
                    //if (A.diag(i,j,k) != 0.0) { // If cell is a fluid.
                    double t = q(i, j, k) - A.plusI(i, j, k) * precon(i, j, k) * z(i + 1, j, k)
                               - A.plusJ(i, j, k) * precon(i, j, k) * z(i, j + 1, k)
                               - A.plusK(i, j, k) * precon(i, j, k) * z(i, j, k + 1);
                    z(i, j, k) = t * precon(i, j, k);
                    //}
                }
    }
    else{
        // Unpreconditioned CG: Bypass preconditioner:
        z = r;
        return;
    }

}

double MACGrid::dotProduct(const GridData & vector1, const GridData & vector2) {
	
	double result = 0.0;

	FOR_EACH_CELL {
		result += vector1(i,j,k) * vector2(i,j,k);
	}

	return result;
}


void MACGrid::add(const GridData & vector1, const GridData & vector2, GridData & result) {
	
	FOR_EACH_CELL {
		result(i,j,k) = vector1(i,j,k) + vector2(i,j,k);
	}

}


void MACGrid::subtract(const GridData & vector1, const GridData & vector2, GridData & result) {
	
	FOR_EACH_CELL {
		result(i,j,k) = vector1(i,j,k) - vector2(i,j,k);
	}

}


void MACGrid::multiply(const double scalar, const GridData & vector, GridData & result) {
	
	FOR_EACH_CELL {
		result(i,j,k) = scalar * vector(i,j,k);
	}

}


double MACGrid::maxMagnitude(const GridData & vector) {
	
	double result = 0.0;

	FOR_EACH_CELL {
		if (abs(vector(i,j,k)) > result) result = abs(vector(i,j,k));
	}

	return result;
}


void MACGrid::apply(const GridDataMatrix & matrix, const GridData & vector, GridData & result) 
{
	FOR_EACH_CELL { // For each row of the matrix.

		double diag = 0;
		double plusI = 0;
		double plusJ = 0;
		double plusK = 0;
		double minusI = 0;
		double minusJ = 0;
		double minusK = 0;

		diag = matrix.diag(i,j,k) * vector(i,j,k);
		if (isValidCell(i+1,j,k)) plusI = matrix.plusI(i,j,k) * vector(i+1,j,k);
		if (isValidCell(i,j+1,k)) plusJ = matrix.plusJ(i,j,k) * vector(i,j+1,k);
		if (isValidCell(i,j,k+1)) plusK = matrix.plusK(i,j,k) * vector(i,j,k+1);
		if (isValidCell(i-1,j,k)) minusI = matrix.plusI(i-1,j,k) * vector(i-1,j,k);
		if (isValidCell(i,j-1,k)) minusJ = matrix.plusJ(i,j-1,k) * vector(i,j-1,k);
		if (isValidCell(i,j,k-1)) minusK = matrix.plusK(i,j,k-1) * vector(i,j,k-1);

		result(i,j,k) = diag + plusI + plusJ + plusK + minusI + minusJ + minusK;
	}

}

void MACGrid::saveSmoke(const char* fileName) 
{
	std::ofstream fileOut(fileName);
	if (fileOut.is_open()) {
		FOR_EACH_CELL {
			fileOut << mD(i,j,k) << std::endl;
		}
		fileOut.close();
	}
}

void MACGrid::saveParticle(std::string filename)
{
	Partio::ParticlesDataMutable *parts = Partio::create();
	Partio::ParticleAttribute posH, vH;
	posH = parts->addAttribute("position", Partio::VECTOR, 3);
	vH = parts->addAttribute("v", Partio::VECTOR, 3);
	for (unsigned int i = 0; i < rendering_particles.size(); i++)
	{
		int idx = parts->addParticle();
		float *p = parts->dataWrite<float>(posH, idx);
		float *v = parts->dataWrite<float>(vH, idx);
		for (int k = 0; k < 3; k++)
		{
			p[k] = rendering_particles[i][k];
			v[k] = rendering_particles_vel[i][k];
		}
	}
	
	Partio::write(filename.c_str(), *parts);
	parts->release();
}

void MACGrid::saveDensity(std::string filename)
{
	Partio::ParticlesDataMutable *density_field = Partio::create();
	Partio::ParticleAttribute posH, rhoH;
	posH = density_field->addAttribute("position", Partio::VECTOR, 3);
	rhoH = density_field->addAttribute("density", Partio::VECTOR, 1);
	FOR_EACH_CELL{
		int idx = density_field->addParticle();
		float *p = density_field->dataWrite<float>(posH, idx);
		float *rho = density_field->dataWrite<float>(rhoH, idx);
		vec3 cellCenter = getCenter(i, j, k);
		for (int l = 0; l < 3; l++)
		{
			p[l] = cellCenter[l];
		}
		rho[0] = getDensity(cellCenter);
	}
	Partio::write(filename.c_str(), *density_field);
	density_field->release();
}

void MACGrid::draw(const Camera& c)
{   
   drawWireGrid();
   if (theDisplayVel) drawVelocities();   
   if (theRenderMode == CUBES) drawSmokeCubes(c);
   else drawSmoke(c);
}

void MACGrid::drawVelocities()
{
   // draw line at each center
   //glColor4f(0.0, 1.0, 0.0, 1.0);
   glBegin(GL_LINES);
      FOR_EACH_CELL
      {
         vec3 pos = getCenter(i,j,k);
         vec3 vel = getVelocity(pos);
         if (vel.Length() > 0.0001)
         {
           //vel.Normalize(); 
           vel *= theCellSize/2.0;
           vel += pos;
		   glColor4f(1.0, 1.0, 0.0, 1.0);
           glVertex3dv(pos.n);
		   glColor4f(0.0, 1.0, 0.0, 1.0);
           glVertex3dv(vel.n);
         }
      }
   glEnd();
}

vec4 MACGrid::getRenderColor(int i, int j, int k)
{
	
	double value = mD(i, j, k); 
	vec4 coldColor(0.5, 0.5, 1.0, value);
	vec4 hotColor(1.0, 0.5, 0.5, value);
    return LERP(coldColor, hotColor, mT(i, j, k));
	

	/*
	// OLD:
    double value = mD(i, j, k); 
    return vec4(1.0, 0.9, 1.0, value);
	*/
}

vec4 MACGrid::getRenderColor(const vec3& pt)
{
	double value = getDensity(pt);
	vec4 coldColor(0.5, 0.5, 1.0, value);
	vec4 hotColor(1.0, 0.5, 0.5, value);
    return LERP(coldColor, hotColor, getTemperature(pt));

	/*
	// OLD:
    double value = getDensity(pt); 
    return vec4(1.0, 1.0, 1.0, value);
	*/
}

void MACGrid::drawZSheets(bool backToFront)
{
   // Draw K Sheets from back to front
   double back =  (theDim[2])*theCellSize;
   double top  =  (theDim[1])*theCellSize;
   double right = (theDim[0])*theCellSize;
  
   double stepsize = theCellSize*0.25;

   double startk = back - stepsize;
   double endk = 0;
   double stepk = -theCellSize;

   if (!backToFront)
   {
      startk = 0;
      endk = back;   
      stepk = theCellSize;
   }

   for (double k = startk; backToFront? k > endk : k < endk; k += stepk)
   {
     for (double j = 0.0; j < top; )
      {
         glBegin(GL_QUAD_STRIP);
         for (double i = 0.0; i <= right; i += stepsize)
         {
            vec3 pos1 = vec3(i,j,k); 
            vec3 pos2 = vec3(i, j+stepsize, k); 

            vec4 color1 = getRenderColor(pos1);
            vec4 color2 = getRenderColor(pos2);

            glColor4dv(color1.n);
            glVertex3dv(pos1.n);

            glColor4dv(color2.n);
            glVertex3dv(pos2.n);
         } 
         glEnd();
         j+=stepsize;

         glBegin(GL_QUAD_STRIP);
         for (double i = right; i >= 0.0; i -= stepsize)
         {
            vec3 pos1 = vec3(i,j,k); 
            vec3 pos2 = vec3(i, j+stepsize, k); 

            vec4 color1 = getRenderColor(pos1);
            vec4 color2 = getRenderColor(pos2);

            glColor4dv(color1.n);
            glVertex3dv(pos1.n);

            glColor4dv(color2.n);
            glVertex3dv(pos2.n);
         } 
         glEnd();
         j+=stepsize;
      }
   }
}

void MACGrid::drawXSheets(bool backToFront)
{
   // Draw K Sheets from back to front
   double back =  (theDim[2])*theCellSize;
   double top  =  (theDim[1])*theCellSize;
   double right = (theDim[0])*theCellSize;
  
   double stepsize = theCellSize*0.25;

   double starti = right - stepsize;
   double endi = 0;
   double stepi = -theCellSize;

   if (!backToFront)
   {
      starti = 0;
      endi = right;   
      stepi = theCellSize;
   }

   for (double i = starti; backToFront? i > endi : i < endi; i += stepi)
   {
     for (double j = 0.0; j < top; )
      {
         glBegin(GL_QUAD_STRIP);
         for (double k = 0.0; k <= back; k += stepsize)
         {
            vec3 pos1 = vec3(i,j,k); 
            vec3 pos2 = vec3(i, j+stepsize, k); 

            vec4 color1 = getRenderColor(pos1);
            vec4 color2 = getRenderColor(pos2);

            glColor4dv(color1.n);
            glVertex3dv(pos1.n);

            glColor4dv(color2.n);
            glVertex3dv(pos2.n);
         } 
         glEnd();
         j+=stepsize;

         glBegin(GL_QUAD_STRIP);
         for (double k = back; k >= 0.0; k -= stepsize)
         {
            vec3 pos1 = vec3(i,j,k); 
            vec3 pos2 = vec3(i, j+stepsize, k); 

            vec4 color1 = getRenderColor(pos1);
            vec4 color2 = getRenderColor(pos2);

            glColor4dv(color1.n);
            glVertex3dv(pos1.n);

            glColor4dv(color2.n);
            glVertex3dv(pos2.n);
         } 
         glEnd();
         j+=stepsize;
      }
   }
}


void MACGrid::drawSmoke(const Camera& c)
{
   vec3 eyeDir = c.getBackward();
   double zresult = fabs(Dot(eyeDir, vec3(1,0,0)));
   double xresult = fabs(Dot(eyeDir, vec3(0,0,1)));
   //double yresult = fabs(Dot(eyeDir, vec3(0,1,0)));

   if (zresult < xresult)
   {      
      drawZSheets(c.getPosition()[2] < 0);
   }
   else 
   {
      drawXSheets(c.getPosition()[0] < 0);
   }
}

void MACGrid::drawSmokeCubes(const Camera& c)
{
   std::multimap<double, MACGrid::Cube, std::greater<double> > cubes;
   FOR_EACH_CELL
   {
      MACGrid::Cube cube;
      cube.color = getRenderColor(i,j,k);
      cube.pos = getCenter(i,j,k);
      cube.dist = DistanceSqr(cube.pos, c.getPosition());
      cubes.insert(make_pair(cube.dist, cube));
   } 

   // Draw cubes from back to front
   std::multimap<double, MACGrid::Cube, std::greater<double> >::const_iterator it;
   for (it = cubes.begin(); it != cubes.end(); ++it)
   {
      drawCube(it->second);
   }
}

void MACGrid::drawWireGrid()
{
   // Display grid in light grey, draw top & bottom

   double xstart = 0.0;
   double ystart = 0.0;
   double zstart = 0.0;
   double xend = theDim[0]*theCellSize;
   double yend = theDim[1]*theCellSize;
   double zend = theDim[2]*theCellSize;

   glPushAttrib(GL_LIGHTING_BIT | GL_LINE_BIT);
      glDisable(GL_LIGHTING);
      glColor3f(0.25, 0.25, 0.25);

      glBegin(GL_LINES);
      for (int i = 0; i <= theDim[0]; i++)
      {
         double x = xstart + i*theCellSize;
         glVertex3d(x, ystart, zstart);
         glVertex3d(x, ystart, zend);

         glVertex3d(x, yend, zstart);
         glVertex3d(x, yend, zend);
      }

      for (int i = 0; i <= theDim[2]; i++)
      {
         double z = zstart + i*theCellSize;
         glVertex3d(xstart, ystart, z);
         glVertex3d(xend, ystart, z);

         glVertex3d(xstart, yend, z);
         glVertex3d(xend, yend, z);
      }

      glVertex3d(xstart, ystart, zstart);
      glVertex3d(xstart, yend, zstart);

      glVertex3d(xend, ystart, zstart);
      glVertex3d(xend, yend, zstart);

      glVertex3d(xstart, ystart, zend);
      glVertex3d(xstart, yend, zend);

      glVertex3d(xend, ystart, zend);
      glVertex3d(xend, yend, zend);
      glEnd();
   glPopAttrib();

   glEnd();
}

#define LEN 0.5
void MACGrid::drawFace(const MACGrid::Cube& cube)
{
   glColor4dv(cube.color.n);
   glPushMatrix();
      glTranslated(cube.pos[0], cube.pos[1], cube.pos[2]);      
      glScaled(theCellSize, theCellSize, theCellSize);
      glBegin(GL_QUADS);
         glNormal3d( 0.0,  0.0, 1.0);
         glVertex3d(-LEN, -LEN, LEN);
         glVertex3d(-LEN,  LEN, LEN);
         glVertex3d( LEN,  LEN, LEN);
         glVertex3d( LEN, -LEN, LEN);
      glEnd();
   glPopMatrix();
}

void MACGrid::drawCube(const MACGrid::Cube& cube)
{
   glColor4dv(cube.color.n);
   glPushMatrix();
      glTranslated(cube.pos[0], cube.pos[1], cube.pos[2]);      
      glScaled(theCellSize, theCellSize, theCellSize);
      glBegin(GL_QUADS);
         glNormal3d( 0.0, -1.0,  0.0);
         glVertex3d(-LEN, -LEN, -LEN);
         glVertex3d(-LEN, -LEN,  LEN);
         glVertex3d( LEN, -LEN,  LEN);
         glVertex3d( LEN, -LEN, -LEN);         

         glNormal3d( 0.0,  0.0, -0.0);
         glVertex3d(-LEN, -LEN, -LEN);
         glVertex3d(-LEN,  LEN, -LEN);
         glVertex3d( LEN,  LEN, -LEN);
         glVertex3d( LEN, -LEN, -LEN);

         glNormal3d(-1.0,  0.0,  0.0);
         glVertex3d(-LEN, -LEN, -LEN);
         glVertex3d(-LEN, -LEN,  LEN);
         glVertex3d(-LEN,  LEN,  LEN);
         glVertex3d(-LEN,  LEN, -LEN);

         glNormal3d( 0.0, 1.0,  0.0);
         glVertex3d(-LEN, LEN, -LEN);
         glVertex3d(-LEN, LEN,  LEN);
         glVertex3d( LEN, LEN,  LEN);
         glVertex3d( LEN, LEN, -LEN);

         glNormal3d( 0.0,  0.0, 1.0);
         glVertex3d(-LEN, -LEN, LEN);
         glVertex3d(-LEN,  LEN, LEN);
         glVertex3d( LEN,  LEN, LEN);
         glVertex3d( LEN, -LEN, LEN);

         glNormal3d(1.0,  0.0,  0.0);
         glVertex3d(LEN, -LEN, -LEN);
         glVertex3d(LEN, -LEN,  LEN);
         glVertex3d(LEN,  LEN,  LEN);
         glVertex3d(LEN,  LEN, -LEN);
      glEnd();
   glPopMatrix();
}