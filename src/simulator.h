#ifndef simulator_H_
#define simulator_H_

#include <fstream>
#include <Partio.h>
#include "util/constants.h" 
#include "util/open_gl_headers.h" 
#include "util/stb_image_write.h" 
#include "util/custom_output.h" 
#include "util/basic_math.h"
#include "mac_grid.h"

class Camera;
class Simulator
{
public:
    Simulator() {}
    
    virtual void reset() = 0;
    virtual void step() = 0;
    
    void draw(const Camera &c);
    void setRecording(bool on, int width, int height);
    bool isRecording();
    int getTotalFrames();
    void drawAxes();
    void grabScreen();

protected:
    MACGrid mGrid;
    bool mRecordEnabled;
    int mFrameNum;
    int mTotalFrameNum; 

    int recordWidth;
    int recordHeight;
};

#endif