#ifndef PARTICLE_H
#define MACGrid_H_

#pragma warning(disable: 4244 4267 4996)

#include "util/vec.h"

class Particle 
{
public:
    Particle() : position(0, 0, 0), velocity(0, 0, 0) {}
    Particle(vec3 pos, vec3 vel, vec3 idx) : position(pos), velocity(vel), index(idx) {}

    vec3 index;
    vec3 position;
    vec3 velocity;
};

#endif