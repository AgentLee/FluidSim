# Development 

## Projection
The idea for projection is to make the fluid incompressible and also enforce boundary conditions. 

### Pressure Update
- Apply a pressure update to every velocity component that borders a fluid cell
- Equations 4.4 - 4.8 might affect grid cells outside the fluid region
    - We need to enforce some kind of boundary constraint
- **Dirichlet Boundary Condition**
    - In a water simulation, any cell that doesn't have fluid has **0 pressure** (Air cells)
    - *p<sub>i, j, k</sub>* gets set to 0 in equations 4.4 - 4.8
- **Neumann Boundary Condition**
    - Pressure at solid walls
- the pressure difference between on the boundary faces is really important
    - we don't have to store pressure for a solid cell
    - we can just solve equation 4.10

## Discrete Divergence
- Discrete Divergence is essentially estimating the rate of fluid coming in and out of the cell
- Approximate the incompressibility condition with finite differences
    - Each cell must be 0 

## Pressure Equations
- We know how to update velocities using pressure (Eq. 4.4 - 4.8)
- We know how to estimate divergence (Eq. 4.13)

- The goal is to have the final velocity (u<sup>n + 1</sup>) to be divergence free if inside a fluid cell.
    - We need to find the pressure that can do this.
        - Substitute pressure updates in Eq. 4.4 - 4.8 to divergence formula Eq. 
        - u<sup>n + 1</sup><sub>i + i + &frac12;, j, k</sub> (4.4) for u<sub>i + &frac12;, j, k</sub> (4.13)

- If cell (i, j + 1) is an air cell, we set the pressure to 0
- If cell (i + 1, j) is a solid cell, we set the pressure to Eq. 4.10
- If cells (i - 1, j) and (i, j - 1) are fluid cells
    - Eq. 4.23 and 4.24
    - 



- u<sub>i + &frac12;, j, k</sub><sup>n + 1</sup> = u<sub>i + &frac12;, j, k</sub>



## Animating Sand as Fluid
### 4.2: PIC Methods
- PIC
    - Fluid variables at a grid point were initialized as a weighted average of nearby particles and updated on the grid
    - New particle values get interpolated from the grid
    - The particle gets moved with the velocity field of the grid
    - Bad because of excessive numerical diffusion from the constant averaging and interpolation
- FLIP
    - Particles act as the fundamental fluid representation
    - The grid is used to move the particles
- PIC/FLIP
    ```
    Init particles
    For each time step  {
        Compute weighted average of nearby particle velocities at each staggered MAC grid cell

        FLIP: save grid velocities

        Do non-advection steps on the grid

        FLIP: subtract new grid velocities from saved velocities and add interpolated difference to each particle

        PIC: interpolate new grid velocity to particles

        Move particles with ODE solver

        Write particle positions to target
    }
    ```
    - We don't have to implement grid advection, vorticity confinement to counter numerical dissipation

#### 4.2.2: Transferring to the Grid
-