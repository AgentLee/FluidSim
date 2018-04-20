# PIC/FLIP 

## General Algorithm
```
init particle positions and velocities

for each time step {
    // Particle to Grid
    for each particle {
        compute weighted average of nearby particle velocities
    }

    // FLIP save
    for each face {
        save grid velocities for FLIP
    }

    // Do non advection steps
    account for gravity
    enforce boundary conditions and incompressibility

    // Update FLIP velocities
    subtract new velociites from saved FLIP velocities
    // Update particles
    for each particle {
        particle.vel -= FLIP velocity
    }

    // Grid to Particle
    for each face {
        grid.velocity = interpolate from particle.vel
    }

    // Advect Particles
    for each particle {
        RK2
    }

    // Lerp PIC/FLIP
}
```