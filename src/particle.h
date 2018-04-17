#ifndef PARTICLE_H
#define MACGrid_H_

#pragma warning(disable: 4244 4267 4996)

#include "util/vec.h"

class Particle 
{
public:
    Particle() : position(0, 0, 0), velocity(0, 0, 0) {}
    Particle(vec3 pos, vec3 vel) : position(pos), velocity(vel) {}

    vec3 position;
    vec3 velocity;
};

#endif