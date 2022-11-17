//root
#include <TMath.h>

//C, C++
#include <vector>
#include <math.h>

//specific
#include "geometry.h"

using namespace std;

/******** VARIABLES ********/
//box dimensions
float boxL = 500;
float boxH = 250;

//wom dimensions
float womL = 200;//length
float womDout = 70;//outer diameter [mm]
float womDin = 42;//inner diameter [mm]
float womR = 30;//(inner radius = 25)//wom radius [mm]
float womRin = 26;


/******** FUNCTIONS ********/
std::vector<float> getStartPos(float horisontal, float vertical, float angle){
	//return start position of track in the box (x,y,z,angle) [mm,mm,mm,rad]
    std::vector<float> pos(4,0);
    pos[3]=angle*TMath::Pi()/180;
    pos[1]=(boxL/2-vertical*10);
    pos[0]=boxL/2+horisontal*10*TMath::Cos(pos[3])-(boxH/2-horisontal*10*TMath::Sin(pos[3]))*TMath::Tan(pos[3]);
    if(pos[0]<0){
      pos[2]=fabs(pos[0])/TMath::Tan(fabs(pos[3]));
      pos[0]=0;
    }
    else if(pos[0]>boxL){
      pos[2]=(pos[0]-boxL)/TMath::Tan(fabs(pos[3]));
      pos[0]=boxL;
    }
    else {
      pos[2]=0;
    }
    return pos;
}

std::vector<float> solidAngleFactor(const std::vector<float> &startPos,const std::vector<float> &pmtPos){
  std::vector<float> result(2,0);
  float x0 = startPos[0];
  float y0 = startPos[1];
  float z0 = startPos[2];
  float angle = startPos[3];
  
  
  float length = 0;
  float dLength = 1; //mm
  float integratedSolidAngleFactor = 0;
  float x = x0;
  float y = y0;
  float z = z0;  
  while(z<=boxH && x<=boxL && z>=0 && x>=0){
    x=x+TMath::Sin(angle);
    z=z+TMath::Cos(angle);
    length=length+dLength;
    integratedSolidAngleFactor += dLength*getSolidAngle(x,y,z,pmtPos[0],pmtPos[1])/(4*TMath::Pi());
    //printf("solidAngleFactor: %4.2f %4.2f %4.2f %4.2f %d\n",x,y,z,length,getZone(x,y,z,pmt2Pos[0],pmt2Pos[1]));
    //printf("%d ",getZone(x,y,z,pmtPos[0],pmtPos[1]));
  } 
  //printf("integratedSolidAngleFactor: %4.2f\n",integratedSolidAngleFactor);
  result[0]=integratedSolidAngleFactor;
  result[1]=length;
  return result;
}

float solidAngleABH(float A,float B, float H){//return solid angle of pyramid with rectangular base (a,b) and hight of H
  return 4*TMath::ASin((A*B)/sqrt((A*A+H*H)*(B*B+H*H)));
}

int getZone(float x,float y,float z,float pmtX, float pmtY){
  float r2 = (x-pmtX)*(x-pmtX)+(y-pmtY)*(y-pmtY);
  if(z>=(boxH-womL)){
    if(r2>(womDout*womDout/4))return 1;
    else if(r2>(womDout*womDout/4))return 2;
    else return 3;
  }
  else{
    if(r2>=(womR*womR))return 4;
    else return 5;
  }
}

float getSolidAngle(float x,float y,float z,float pmtX, float pmtY){//return solid angle per one step
  float solidAngle = 0;
  float x1=-999,y1=-999,z1=-999;
  float x2=-999,y2=-999,z2=-999;
  float dr = (4-TMath::Pi())*womR/4;
  float r = TMath::Pi()*womR/4;
  int zone = getZone(x,y,z,pmtX,pmtY);
  if(zone==1){
    float h = sqrt((x-pmtX)*(x-pmtX)+(y-pmtY)*(y-pmtY)) - r;
    float l1 = z+womL-boxH;
    float l2 = boxH-z;
    float a1 = sqrt(h*h+l1*l1);
    float a2 = sqrt(h*h+l2*l2);
    float A = a1*sqrt(2*(1-(h*h-l1*l2)/(a1*a2)));
    float B = 2*womR;
    float H = sqrt(a1*a1-A*A/4);
    return solidAngleABH(A,B,H);
  }
  else if(zone==4){
    float h1 = sqrt((x-pmtX)*(x-pmtX)+(y-pmtY)*(y-pmtY)) + r;
    float h2 = sqrt((x-pmtX)*(x-pmtX)+(y-pmtY)*(y-pmtY)) - r;
    float l1 = boxH-womL-z;
    float l2 = boxH-z;
    float a1 = sqrt(h1*h1+l1*l1);
    float a2 = sqrt(h2*h2+l2*l2);
    float cosA = (a1*a1+a2*a2-4*r*r-womL*womL)/(2*a1*a2);
    float A = a1*sqrt(2*(1-cosA));
    float B = 2*womR;
    float H = sqrt(a1*a1-A*A/4);
    return solidAngleABH(A,B,H);
  }
  else if(zone==2){
    float l1 =womRin - sqrt((x-pmtX)*(x-pmtX)+(y-pmtY)*(y-pmtY));
    float l2 =womRin + sqrt((x-pmtX)*(x-pmtX)+(y-pmtY)*(y-pmtY));
    float h = boxH - z;
    float H = womL - h;
    float a1 = sqrt(l1*l1+h*h);
    float a2 = sqrt(l2*l2+h*h);
    float A1 = sqrt(l1*l1+H*H);
    float A2 = sqrt(l2*l2+H*H);
    return 2*TMath::Pi()*(sqrt((1+(h*h-l1*l2)/(a1*a2))/2)+sqrt((1+(H*H-l1*l2)/(A1*A2))/2));
  }
  else if(zone==5){
    float l1 =womR - sqrt((x-pmtX)*(x-pmtX)+(y-pmtY)*(y-pmtY));
    float l2 =womR + sqrt((x-pmtX)*(x-pmtX)+(y-pmtY)*(y-pmtY));
    float h = boxH - z -womL;
    float H = boxH - z;
    float a1 = sqrt(l1*l1+h*h);
    float a2 = sqrt(l2*l2+h*h);
    float A1 = sqrt(l1*l1+H*H);
    float A2 = sqrt(l2*l2+H*H);
    return 2*TMath::Pi()*(sqrt((1+(H*H-l1*l2)/(A1*A2))/2)-sqrt((1+(h*h-l1*l2)/(a1*a2))/2));
  }
  else return 0;
}