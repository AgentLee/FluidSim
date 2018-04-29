#ifndef smokeSim_H_
#define smokeSim_H_

#include "simulator.h"

class Camera;
class SmokeSim : public Simulator
{
public:
   SmokeSim();
   virtual ~SmokeSim() {}

   void reset();
   void step();
};

#endif