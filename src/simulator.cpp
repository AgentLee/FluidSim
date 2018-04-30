#define STB_IMAGE_WRITE_IMPLEMENTATION

#include "simulator.h"

void Simulator::draw(const Camera& c)
{
    drawAxes();
    mGrid.draw(c);
    if(mRecordEnabled) grabScreen();
}

void Simulator::setRecording(bool on, int width, int height)
{
    if (on && ! mRecordEnabled)  // reset counter
    {
        mFrameNum = 0;
    }
    mRecordEnabled = on;
    
    recordWidth = width;
    recordHeight = height;
}

bool Simulator::isRecording() 
{ 
    return mRecordEnabled; 
}

int Simulator::getTotalFrames() 
{ 
    return mTotalFrameNum; 
}

void Simulator::drawAxes()
{
    glPushAttrib(GL_LIGHTING_BIT | GL_LINE_BIT);
        glDisable(GL_LIGHTING);

        glLineWidth(2.0); 
        glBegin(GL_LINES);
            glColor3f(1.0, 0.0, 0.0);
            glVertex3f(0.0, 0.0, 0.0);
            glVertex3f(1.0, 0.0, 0.0);

            glColor3f(0.0, 1.0, 0.0);
            glVertex3f(0.0, 0.0, 0.0);
            glVertex3f(0.0, 1.0, 0.0);

            glColor3f(0.0, 0.0, 1.0);
            glVertex3f(0.0, 0.0, 0.0);
            glVertex3f(0.0, 0.0, 1.0);
        glEnd();
    glPopAttrib();
}

void Simulator::grabScreen()
{
    if (mFrameNum > 9999) exit(0);

    // Save density field to a .bgeo file
    std::string densityFile = "../records/DensityFrame" + std::to_string(mFrameNum) + ".bgeo";
    mGrid.saveDensity(densityFile);

    // Save an image:
    unsigned char* bitmapData = new unsigned char[3 * recordWidth * recordHeight];
    for (int i=0; i<recordHeight; i++) 
    {
        glReadPixels(0,i,recordWidth,1,GL_RGB, GL_UNSIGNED_BYTE, 
            bitmapData + (recordWidth * 3 * ((recordHeight-1)-i)));
    }
    char anim_filename[2048];
    snprintf(anim_filename, 2048, "../records/smoke_%04d.png", mFrameNum); 
    stbi_write_png(anim_filename, recordWidth, recordHeight, 3, bitmapData, recordWidth * 3);
    delete [] bitmapData;

    // Dump out rendering particle data in .bgeo file
    std::string particleFile = "../records/frame" + std::to_string(mFrameNum) + ".bgeo";
    mGrid.saveParticle(particleFile);

    mFrameNum++;
}