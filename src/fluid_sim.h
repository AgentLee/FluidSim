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
};

#endif