#ifndef GEOMETRY
#define GEOMETRY

#include <vector>

using namespace std;

std::vector<float> getStartPos(float horisontal, float vertical, float angle);//return start position of track in the box (x,y,z,angle) [mm,mm,mm,rad]
std::vector<float> solidAngleFactor(const std::vector<float> &startPos,const std::vector<float> &pmtPos);//return solid angle correction factor and path length in the box (factor,length) [,mm]
float solidAngleABH(float A,float B, float H);//return solid angle of pyramid with rectangular base (a,b) and hight of H
int getZone(float x,float y,float z,float pmtX, float pmtY);
float getSolidAngle(float x,float y,float z,float pmtX, float pmtY);

#endif