#include "smoke_sim.h"

// MACGrid target;
MACGrid target;

SmokeSim::SmokeSim()
{
	mGrid.theRenderMode = MACGrid::SHEETS;
	mFrameNum = 0;
	mTotalFrameNum = 0;
	mRecordEnabled = false;

	reset();
}

void SmokeSim::reset()
{
	mGrid.reset();
	calculateAMatrix(mGrid);
	calculatePreconditioner(mGrid, mGrid.AMatrix);
	mTotalFrameNum = 0;
}

void SmokeSim::updateSources(MACGrid &mGrid)
{
	// Set initial values for density, temperature, velocity
    for(int i=6; i<12;i++){
        for(int j=0; j<5; j++){
            mGrid.mV(i,j+1,0) = 2.0;
            mGrid.mD(i,j,0) = 1.0;
            mGrid.mT(i,j,0) = 1.0;

            mGrid.mV(i,j+2,0) = 2.0;
            mGrid.mD(i,j,0) = 1.0;
            mGrid.mT(i,j,0) = 1.0;
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
                    mGrid.rendering_particles.push_back(xp);
                }
			}
		}
	}
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
void SmokeSim::advectVelocity(MACGrid &mGrid, double dt)
{
	// TODO: Calculate new velocities and store in target
    // TODO: Your code is here. It builds target.mU, target.mV and target.mW for all faces
    
    // Velocities are stored at the faces of the MAC grid
    FOR_EACH_FACE
    {
        if(mGrid.isValidFace(MACGrid::X, i, j, k)) {
            vec3 currPos = mGrid.getFacePosition(MACGrid::X, i, j, k);
            vec3 startPos = mGrid.getRewoundPosition(currPos, dt);

            // q_n+1 = interpolate(q_n, xG - dt(uG));
            target.mU(i, j, k) = mGrid.getVelocityX(startPos);
        }
        
        if(mGrid.isValidFace(MACGrid::Y, i, j, k)) {
            vec3 currPos = mGrid.getFacePosition(MACGrid::Y, i, j, k);
            vec3 startPos = mGrid.getRewoundPosition(currPos, dt);

            // q_n+1 = interpolate(q_n, xG - dt(uG));
            target.mV(i, j, k) = mGrid.getVelocityY(startPos);
        }

        if(mGrid.isValidFace(MACGrid::Z, i, j, k)) {
            vec3 currPos = mGrid.getFacePosition(MACGrid::Z, i, j, k);
            vec3 startPos = mGrid.getRewoundPosition(currPos, dt);

            // q_n+1 = interpolate(q_n, xG - dt(uG));
            target.mW(i, j, k) = mGrid.getVelocityZ(startPos);
        }  
    }

    // Then save the result to our object
    mGrid.mU = target.mU;
    mGrid.mV = target.mV;
    mGrid.mW = target.mW;
}

/*
 * pg45
 */
void SmokeSim::computeBouyancy(MACGrid &mGrid, double dt)
{
	// TODO: Calculate bouyancy and store in target
    // TODO: Your code is here. It modifies target.mV for all y face velocities.

	FOR_EACH_FACE 
    {
		if (mGrid.isValidFace(MACGrid::Y, i, j, k)) {
			vec3 currPos = mGrid.getFacePosition(MACGrid::Y, i, j, k);

            double temp = mGrid.getTemperature(currPos);
    	    double ambientTemp = 300;
            double s = mGrid.getDensity(currPos);

            // Equation 5.1
            vec3 fBuoy(0, -theBuoyancyAlpha * s + theBuoyancyBeta * (temp - ambientTemp), 0);

            target.mV(i, j, k) = mGrid.mV(i, j, k) + fBuoy[1];
		}
	}

    // and then save the result to our object
    mGrid.mV = target.mV;
}

/*
 * This function prevents the vortices from disappearing too quickly.
 * w = cross(gradient, u) [curl]
 * A vortex is a peak in the vorticity field
 * pg46
 */ 
void SmokeSim::computeVorticityConfinement(MACGrid &mGrid, double dt)
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
        double x = ((mGrid.mW(i, j + 1, k) - mGrid.mW(i, j - 1, k)) / (2 * theCellSize)) - 
                   ((mGrid.mV(i, j, k + 1) - mGrid.mV(i, j, k - 1)) / (2 * theCellSize));
        double y = ((mGrid.mU(i, j, k + 1) - mGrid.mU(i, j, k - 1)) / (2 * theCellSize)) - 
                   ((mGrid.mW(i + 1, j, k) - mGrid.mW(i - 1, j, k)) / (2 * theCellSize));
        double z = ((mGrid.mV(i + 1, j, k) - mGrid.mV(i - 1, j, k)) / (2 * theCellSize)) - 
                   ((mGrid.mU(i, j + 1, k) - mGrid.mU(i, j - 1, k)) / (2 * theCellSize));
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

        applyVorticityConfinement(mGrid, fConf, i, j, k);
    }

    // Then save the result to our object
    mGrid.mU = target.mU;
    mGrid.mV = target.mV;
    mGrid.mW = target.mW;
}

void SmokeSim::applyVorticityConfinement(MACGrid &mGrid, vec3 &fConf, int &i, int &j, int &k)
{
    // pg47
    // "Take the appropriate averages to apply this 
    // to the different components of velocity on the MAC grid"
    if (mGrid.isValidFace(MACGrid::X, i, j, k)) {
        target.mU(i, j, k) += fConf[0] / 2;
    } 
    if (mGrid.isValidFace(MACGrid::X, i + 1, j, k)) {
        target.mU(i + 1, j, k) += fConf[0] / 2;
    }
    if (mGrid.isValidFace(MACGrid::Y, i, j, k)) {
        target.mV(i, j, k) += fConf[1] / 2;
    }
    if (mGrid.isValidFace(MACGrid::Y, i, j + 1, k)) {
        target.mV(i, j + 1, k) += fConf[1] / 2;
    }
    if (mGrid.isValidFace(MACGrid::Z, i, j, k)) {
        target.mW(i, j, k) += fConf[2] / 2;
    }
    if (mGrid.isValidFace(MACGrid::Z, i, j, k + 1)) {
        target.mW(i, j, k + 1) += fConf[2] / 2;
    }
}

void SmokeSim::addExternalForces(MACGrid &mGrid, double dt)
{
   computeBouyancy(mGrid, dt);
   computeVorticityConfinement(mGrid, dt);
}

void SmokeSim::computeDivergence(MACGrid &mGrid, GridData &d)
{
    FOR_EACH_FACE
    {
        // Should only do this for fluid cells
        {
            // Use finite differences to approximate the divergence
            double uPlus    = mGrid.mU(i + 1, j, k);
            double uMinus   = mGrid.mU(i, j, k);
            double vPlus    = mGrid.mV(i, j + 1, k);
            double vMinus   = mGrid.mV(i, j, k);
            double wPlus    = mGrid.mW(i, j, k + 1);
            double wMinus   = mGrid.mW(i, j, k);

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

void SmokeSim::apply(MACGrid &mGrid, const GridDataMatrix & matrix, const GridData & vector, GridData & result) 
{
	FOR_EACH_CELL 
	{ // For each row of the matrix.
		double diag = 0;
		double plusI = 0;
		double plusJ = 0;
		double plusK = 0;
		double minusI = 0;
		double minusJ = 0;
		double minusK = 0;

		diag = matrix.diag(i,j,k) * vector(i,j,k);
		if (mGrid.isValidCell(i+1,j,k)) plusI = matrix.plusI(i,j,k) * vector(i+1,j,k);
		if (mGrid.isValidCell(i,j+1,k)) plusJ = matrix.plusJ(i,j,k) * vector(i,j+1,k);
		if (mGrid.isValidCell(i,j,k+1)) plusK = matrix.plusK(i,j,k) * vector(i,j,k+1);
		if (mGrid.isValidCell(i-1,j,k)) minusI = matrix.plusI(i-1,j,k) * vector(i-1,j,k);
		if (mGrid.isValidCell(i,j-1,k)) minusJ = matrix.plusJ(i,j-1,k) * vector(i,j-1,k);
		if (mGrid.isValidCell(i,j,k-1)) minusK = matrix.plusK(i,j,k-1) * vector(i,j,k-1);

		result(i,j,k) = diag + plusI + plusJ + plusK + minusI + minusJ + minusK;
	}
}

bool SmokeSim::preconditionedConjugateGradient(const GridDataMatrix & A, GridData & p, const GridData & d, int maxIterations, double tolerance) 
{
	// Solves Ap = d for p.

	FOR_EACH_CELL {
		p(i,j,k) = 0.0; // Initial guess p = 0.	
	}

	GridData r = d; // Residual vector.

	GridData z; z.initialize();
	applyPreconditioner(mGrid, r, A, z); // Auxillary vector.

	GridData s = z; // Search vector;

	double sigma = mGrid.dotProduct(z, r);

	for (int iteration = 0; iteration < maxIterations; iteration++) {
		double rho = sigma; // According to TA. Here???

		apply(mGrid, A, s, z); // z = applyA(s);

		double alpha = rho/mGrid.dotProduct(z, s);

		GridData alphaTimesS; alphaTimesS.initialize();
		mGrid.multiply(alpha, s, alphaTimesS);
		mGrid.add(p, alphaTimesS, p);
		//p += alpha * s;

		GridData alphaTimesZ; alphaTimesZ.initialize();
		mGrid.multiply(alpha, z, alphaTimesZ);
		mGrid.subtract(r, alphaTimesZ, r);
		//r -= alpha * z;

		if (mGrid.maxMagnitude(r) <= tolerance) {
			//PRINT_LINE("PCG converged in " << (iteration + 1) << " iterations.");
			return true; //return p;
		}

		applyPreconditioner(mGrid, r, A, z); // z = applyPreconditioner(r);

		double sigmaNew = mGrid.dotProduct(z, r);
		double beta = sigmaNew / rho;

		GridData betaTimesS; betaTimesS.initialize();
		mGrid.multiply(beta, s, betaTimesS);
		mGrid.add(z, betaTimesS, s);
		//s = z + beta * s;

		sigma = sigmaNew;
	}

	PRINT_LINE( "PCG didn't converge!" );
	return false;
}

/*
 * Modified Incomplete Cholesky: pg36 course notes
 */ 
void SmokeSim::calculatePreconditioner(MACGrid &mGrid, const GridDataMatrix & A) {

	mGrid.precon.initialize();

    // TODO: Build the modified incomplete Cholesky preconditioner following Fig 4.2 on page 36 of Bridson's 2007 SIGGRAPH fluid course notes.
    //       This corresponds to filling in precon(i,j,k) for all cells

    // Tuning constant
    double t = 0.97;

    FOR_EACH_CELL
    {
        // if cell is a fluid 
        {
            double e = A.diag(i, j, k)  - ((A.plusI(i - 1, j, k) * mGrid.precon(i - 1, j, k)) * (A.plusI(i - 1, j, k) * mGrid.precon(i - 1, j, k)))
                                        - ((A.plusJ(i, j - 1, k) * mGrid.precon(i, j - 1, k)) * (A.plusJ(i, j - 1, k) * mGrid.precon(i, j - 1, k)))
                                        - ((A.plusK(i, j, k - 1) * mGrid.precon(i, j, k - 1)) * (A.plusK(i, j, k - 1) * mGrid.precon(i, j, k - 1)))
                                        - t * ( A.plusI(i - 1, j, k) * (A.plusJ(i - 1, j, k) + A.plusK(i - 1, j, k)) * (mGrid.precon(i - 1, j, k) * mGrid.precon(i - 1, j, k)) +
                                                A.plusJ(i, j - 1, k) * (A.plusI(i, j - 1, k) + A.plusK(i, j - 1, k)) * (mGrid.precon(i, j - 1, k) * mGrid.precon(i, j - 1, k)) +
                                                A.plusK(i, j, k - 1) * (A.plusI(i, j, k - 1) + A.plusJ(i, j, k - 1)) * (mGrid.precon(i, j, k - 1) * mGrid.precon(i, j, k - 1)));

            // Set to a small number to prevent divide by 0
            mGrid.precon(i, j, k) = 1 / sqrt(e + pow(10, -30));
        }
    }
}

/*
 * pg37 in course notes
 */
void SmokeSim::applyPreconditioner(MACGrid &mGrid, const GridData & r, const GridDataMatrix & A, GridData & z)
{
    // TODO: change if(0) to if(1) after you implement calculatePreconditoner function.
    if(1) {
        // APPLY THE PRECONDITIONER:
        // Solve Lq = r for q:
        GridData q;
        q.initialize();
        FOR_EACH_CELL {
                    //if (A.diag(i,j,k) != 0.0) { // If cell is a fluid.
                    double t = r(i, j, k) - A.plusI(i - 1, j, k) * mGrid.precon(i - 1, j, k) * q(i - 1, j, k)
                               - A.plusJ(i, j - 1, k) * mGrid.precon(i, j - 1, k) * q(i, j - 1, k)
                               - A.plusK(i, j, k - 1) * mGrid.precon(i, j, k - 1) * q(i, j, k - 1);
                    q(i, j, k) = t * mGrid.precon(i, j, k);
                    //}
                }
        // Solve L^Tz = q for z:
        FOR_EACH_CELL_REVERSE {
                    //if (A.diag(i,j,k) != 0.0) { // If cell is a fluid.
                    double t = q(i, j, k) - A.plusI(i, j, k) * mGrid.precon(i, j, k) * z(i + 1, j, k)
                               - A.plusJ(i, j, k) * mGrid.precon(i, j, k) * z(i, j + 1, k)
                               - A.plusK(i, j, k) * mGrid.precon(i, j, k) * z(i, j, k + 1);
                    z(i, j, k) = t * mGrid.precon(i, j, k);
                    //}
                }
    }
    else{
        // Unpreconditioned CG: Bypass preconditioner:
        z = r;
        return;
    }

}

void SmokeSim::calculateAMatrix(MACGrid &mGrid) 
{
	FOR_EACH_CELL 
    {
		int numFluidNeighbors = 0;
		if (i-1 >= 0) {
			mGrid.AMatrix.plusI(i-1,j,k) = -1;
			numFluidNeighbors++;
		}
		if (i+1 < theDim[MACGrid::X]) {
			mGrid.AMatrix.plusI(i,j,k) = -1;
			numFluidNeighbors++;
		}
		if (j-1 >= 0) {
			mGrid.AMatrix.plusJ(i,j-1,k) = -1;
			numFluidNeighbors++;
		}
		if (j+1 < theDim[MACGrid::Y]) {
			mGrid.AMatrix.plusJ(i,j,k) = -1;
			numFluidNeighbors++;
		}
		if (k-1 >= 0) {
			mGrid.AMatrix.plusK(i,j,k-1) = -1;
			numFluidNeighbors++;
		}
		if (k+1 < theDim[MACGrid::Z]) {
			mGrid.AMatrix.plusK(i,j,k) = -1;
			numFluidNeighbors++;
		}
		// Set the diagonal:
		mGrid.AMatrix.diag(i,j,k) = numFluidNeighbors;
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
void SmokeSim::project(MACGrid &mGrid, double dt)
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
    computeDivergence(mGrid, d);

    // Solve for p
    if(!preconditionedConjugateGradient(mGrid.AMatrix, p, d, 150, 0.01)) {
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
        if(mGrid.isValidFace(MACGrid::X, i, j, k)) {
            if(i - 1 < 0 || i >= theDim[MACGrid::X]) {
                target.mU(i, j, k) = 0;
            }
            else {
                target.mU(i, j, k) -= scale * (p(i, j, k) - p(i - 1, j, k));
            }
        }

        if(mGrid.isValidFace(MACGrid::Y, i, j, k)) {
            if(j - 1 < 0 || j >= theDim[MACGrid::Y]) {
                target.mV(i, j, k) = 0;
            }
            else {
                target.mV(i, j, k) -= scale * (p(i, j, k) - p(i, j - 1, k));
            }
        }

        if(mGrid.isValidFace(MACGrid::Z, i, j, k)) {
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
        mGrid.mP = target.mP; 
        mGrid.mU = target.mU;
        mGrid.mV = target.mV;
        mGrid.mW = target.mW;
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

void SmokeSim::advectTemperature(MACGrid &mGrid, double dt)
{
    // TODO: Calculate new temp and store in target
    // TODO: Your code is here. It builds target.mT for all cells.

    // Temperature is contained at the centers
    FOR_EACH_CELL
    {
        vec3 currPos = mGrid.getCenter(i, j, k);
        vec3 startPos = mGrid.getRewoundPosition(currPos, dt);

        target.mT(i, j, k) = mGrid.getTemperature(startPos);
    }

    // Then save the result to our object
    mGrid.mT = target.mT;
}

void SmokeSim::advectDensity(MACGrid &mGrid, double dt)
{
    // TODO: Calculate new densitities and store in target
    // TODO: Your code is here. It builds target.mD for all cells.

    // Density is contained at the centers
    FOR_EACH_CELL
    {
        vec3 currPos = mGrid.getCenter(i, j, k);
        vec3 startPos = mGrid.getRewoundPosition(currPos, dt);

        target.mD(i, j, k) = mGrid.getDensity(startPos);
    }

    // Then save the result to our object
    mGrid.mD = target.mD;
}

void SmokeSim::advectRenderingParticles(MACGrid &mGrid, double dt) 
{
	mGrid.rendering_particles_vel.resize(mGrid.rendering_particles.size());
	for (size_t p = 0; p < mGrid.rendering_particles.size(); p++) {
		vec3 currentPosition = mGrid.rendering_particles[p];
        vec3 currentVelocity = mGrid.getVelocity(currentPosition);
        vec3 nextPosition = currentPosition + currentVelocity * dt;
        vec3 clippedNextPosition = mGrid.clipToGrid(nextPosition, currentPosition);
        // Keep going...
        vec3 nextVelocity = mGrid.getVelocity(clippedNextPosition);
        vec3 averageVelocity = (currentVelocity + nextVelocity) / 2.0;
        vec3 betterNextPosition = currentPosition + averageVelocity * dt;
        vec3 clippedBetterNextPosition = mGrid.clipToGrid(betterNextPosition, currentPosition);
        mGrid.rendering_particles[p] = clippedBetterNextPosition;
		mGrid.rendering_particles_vel[p] = averageVelocity;
	}
}

void SmokeSim::step()
{
	double dt = 0.04;//0.1;

	// Step0: Gather user forces
	updateSources(mGrid);

	// Step1: Calculate new velocities
	advectVelocity(mGrid, dt);
	addExternalForces(mGrid, dt);
	project(mGrid, dt);

	// Step2: Calculate new temperature
	advectTemperature(mGrid, dt);

	// Step3: Calculate new density 
	advectDensity(mGrid, dt);
	
	// Step4: Advect rendering particles
	advectRenderingParticles(mGrid, dt);

	mTotalFrameNum++;
}