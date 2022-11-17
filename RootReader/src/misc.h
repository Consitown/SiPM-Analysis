#ifndef CALIB
#define CALIB

#include <TString.h>
#include <TPolyMarker.h>
#include <tuple>

using namespace std;

vector<float> readVector(string workingDir,string _runName,double initValue);
string vectorToString(vector<float> vec);
bool parseBoolean(string in);
pair<vector<float>,vector<float>> readPair(string path, string _runName, double initValueFirst,double initValueSecond);
#endif