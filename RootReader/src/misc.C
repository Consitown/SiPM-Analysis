//root
#include <TH1F.h>
#include <TF1.h>
#include <TMath.h>
#include <TString.h>
#include <TSpectrum.h>
#include <TPolyMarker.h>
#include <TGraphErrors.h>

//C, C++
#include <math.h>
#include <vector>
#include <assert.h>
#include <stdio.h>
#include <sstream>
#include <tuple>
//specific
#include "misc.h"
#include <fstream>

using namespace std;
string extractValues(string str)
{
  unsigned first = str.find("{");
  unsigned last = str.find("}");
  string strNew = str.substr(first + 1, last - first - 1);
  return strNew;
}


bool parseBoolean(string in){
  if(in == "0")return false;
  else if(in=="false") return false;
  else return true;
}

double stringToDouble(string text)
{
  return stod(text);
}
float stringToFloat(string text)
{
  if(text.length()==0)return 0.0;
  float test = 0;
  try {
    test = stof(text);
     
    // cout << "stof test " << test << endl;
  } catch (const std::exception &e)
  {
    std::cerr <<"Error at stof:" <<e.what() << '\n';
    // cout << "\n";
    // cout << text << "\n" << endl;
    throw 20;
    
  }
  return test;
}
string vectorToString(vector<float> vec)
{
  std::ostringstream vts;

  if (!vec.empty())
  {
    // Convert all but the last element to avoid a trailing ","
    std::copy(vec.begin(), vec.end() - 1,
              std::ostream_iterator<float>(vts, ", "));

    // Now add the last element with no delimiter
    vts << vec.back();
  }
  return vts.str();
}

vector<float> readVector(string calib_path, string _runName, double initValue)
{
  std::ifstream infile(calib_path.c_str());
  vector<float> calib_amp;

  std::string line;

  while (std::getline(infile, line))
  {
   // if ((line.find(_runName) != std::string::npos) || (line.find("all") != std::string::npos))
    if ((line.find(_runName) != std::string::npos))
    {

      calib_amp.clear();
      string s = extractValues(line);
      string delimiter = ",";

      size_t pos = 0;
      std::string token;
      string calibValue;
      while ((pos = s.find(delimiter)) != std::string::npos)
      {
        token = s.substr(0, pos);

        calibValue = string(token);
        calib_amp.push_back(stringToDouble(calibValue));

        s.erase(0, pos + delimiter.length());
      }
      calibValue = string(s);

      calib_amp.push_back(stringToDouble(calibValue));
    }
  }

  infile.close();

  for (std::size_t i = calib_amp.size(); i < 32; i++)
  {
    calib_amp.push_back(initValue);
  }

  return calib_amp;
}


pair<vector<float>,vector<float>> readPair(string path, string _runName, double initValueFirst,double initValueSecond)
{
  std::ifstream infile(path.c_str());
  vector<float> peakSignals;
  vector<float> allSignals;

  std::string line;

  while (std::getline(infile, line))
  {
  //  if ((line.find(_runName) != std::string::npos) || (line.find("all") != std::string::npos)) Dangerous, if all is included in names
    if ((line.find(_runName) != std::string::npos) )
    {

      peakSignals.clear();
      allSignals.clear();

      string s = extractValues(line);
      string columnDelimiter = ",";
      string valueDelimiter = "/";

      size_t pos = 0;
      std::string token;
      string calibValue;
      // cout << " string _runname " << _runName << "\n" << endl;
      // cout << " string s " << s << "\n" << endl;
      while ((pos = s.find(columnDelimiter)) != std::string::npos)
      {
        token = s.substr(0, pos); //LIKE: 33/178

        int posOfDelimiter=token.find("/");

        string peak=token.substr(0,posOfDelimiter);
        string all=token.substr(posOfDelimiter+1, token.length()-1);
        // cout << "peak " << peak << " \n" << endl;
        peakSignals.push_back(stringToFloat(peak));
        // cout << " all " << all << " \n " << endl;
        allSignals.push_back(stringToFloat(all));

        s.erase(0, pos + columnDelimiter.length());
      }
      int posOfDelimiter=s.find("/");
        string peak=s.substr(0,posOfDelimiter);
        string all=s.substr(posOfDelimiter+1, s.length()-1);

        // cout << " after while peak " << peak << " \n" << endl;
        peakSignals.push_back(stringToFloat(peak));
        // cout << " after while all " << all << " \n " << endl;
        allSignals.push_back(stringToFloat(all));
    }
  }

  infile.close();

  for (std::size_t i = peakSignals.size(); i < 32; i++)
  {
    peakSignals.push_back(initValueFirst);
    allSignals.push_back(initValueSecond);

  }     
  


  return pair<vector<float>,vector<float>>(peakSignals,allSignals);
}
