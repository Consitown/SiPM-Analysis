//C, C++
#include <stdio.h>
#include <assert.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <sstream>
#include <numeric>
#include <tuple>
#include <map>

#include "meanAngleFuncs.h"

using namespace std;

int main(int argc, char *argv[]) {

  float angles[] = { 0,45,90,135,180,225,270,315 };
  float cartX[] = {10.2361, 8.06099, 7.7071e-16, -11.4438, -10.0545, -8.38805, -2.43603e-15, 8.99222};
  float cartY[] = {0, 8.06099, 12.5867, 11.4438, 1.23133e-15, -8.38805, -13.2611, -8.99222};
  float newAngles[8];

  cartesianToPolar(8, cartX, cartY, newAngles);
  printArray(newAngles, 8);

  int n = 8;
  printArray(angles, n);

  translateAngle(angles, n, 0);

  printArray(angles, n);






}