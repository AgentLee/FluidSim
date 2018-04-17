#ifndef smokeSim_H_
#define smokeSim_H_

#include <fstream>
#include <Partio.h>
#include "util/constants.h" 
#include "util/open_gl_headers.h" 
#include "util/stb_image_write.h" 
#include "util/custom_output.h" 
#include "util/basic_math.h"
#include "mac_grid.h"

class Camera;
class SmokeSim
{
public:
   SmokeSim();
   virtual ~SmokeSim();

   virtual void reset();
   virtual void step();
   virtual void draw(const Camera& c);
   //virtual void setGridDimensions(int x, int y, int z); 
   virtual void setRecording(bool on, int width, int height);
   virtual bool isRecording();
	
	
	int getTotalFrames();

protected:
   virtual void drawAxes();
   virtual void grabScreen();

protected:
	MACGrid mGrid;
	bool mRecordEnabled;
	int mFrameNum;
	int mTotalFrameNum; 
	
	
	int recordWidth;
	int recordHeight;
};

#endif