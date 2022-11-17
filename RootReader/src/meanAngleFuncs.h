//Including root functionalities:
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TLine.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TString.h>
#include <TCanvas.h>
#include <TPaveLabel.h>
#include <TPad.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TMath.h>
#include <TLegend.h>
#include <TCut.h>
#include <THStack.h>
#include <TGaxis.h>
#include <TF1.h>
#include <TError.h> // root verbosity level
#include <TApplication.h>
#include <TNtuple.h>
#include <TImage.h>
#include <TAttImage.h>
#include <TSpectrum.h>

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

// root
#include <TMath.h>

using namespace std;

float translateAngle(float inputAngle, float zeroAngle); 
void translateAngle(float *inputAngles ,int sizeAngles, float zeroAngle);
void cartesianToPolar(int numberEntries, float* inputX, float* inputY, float* outputAngles);
float cartesianToPolar(float inputX, float inputY);
void arctan2Arr(int numberEntries, float* inputX, float* inputY, float* outputAngles);
float arctan2(float inputX, float inputY);
void printArray(float *array, int arraySize);
Double_t fitf(Double_t *x, Double_t *par);


/**
 * @brief Fitfunction for three fold gaussian with same mean, std dev, height
 * 
 * @param x 
 * @param par - p0 height, p1 mean, p2 std dev, p3 base
 * @return Double_t 
 */
Double_t fitf(Double_t *x, Double_t *par) {
    Double_t p0 = par[0];
    Double_t p1 = par[1];
    Double_t p2 = par[2];
    Double_t p3 = par[3];
    Double_t val = 0;

    val = p0 * (TMath::Exp(-0.5 * (x[0] - p1 - 360)/p2 * (x[0] - p1 - 360)/p2) + TMath::Exp(-0.5 * (x[0] - p1 + 360)/p2 * (x[0] - p1 + 360)/p2) + TMath::Exp(-0.5 * (x[0] - p1)/p2 * (x[0] - p1)/p2)) + p3;

    return val;

}


/**
 * @brief print array, one line no commas
 * 
 * @param array 
 * @param arraySize 
 */
void printArray(float *array, int arraySize) {
    for (int i=0; i<arraySize; i++) {
        cout << array[i] << " ";
    }
    cout << endl;
}

void printArray(Double_t *array, int arraySize) {
    for (int i=0; i<arraySize; i++) {
        cout << array[i] << " ";
    }
    cout << endl;
}

/**
 * @brief Convert cartesian coordinates to polar angle [0, 360)
 * 
 * @param inputX 
 * @param inputY 
 * @return float 
 */
float cartesianToPolar(float inputX, float inputY) {
    float pi = TMath::Pi();
    float outputAngle;
    if (inputX > 0 && inputY >= 0) {
        outputAngle = (TMath::ATan(inputY / inputX) * 180.0 / pi);
    } 
    else if (inputX < 0) {
        outputAngle = (TMath::ATan(inputY / inputX)  + pi) * 180.0 / pi;
    } 
    else if (inputX > 0 && inputY < 0) {
        outputAngle = (TMath::ATan(inputY / inputX)  + 2*pi) * 180.0 / pi;
    } 
    else if (inputX == 0 && inputY > 0) {
        outputAngle = 0.5 * 180.0;
    }
    else if (inputX == 0 && inputY < 0) {
        outputAngle = -1.5* 180.0;
    }
    return outputAngle;
}

/**
 * @brief Convert list of cartesian coordinates to corresponding polar angle [0, 360)
 * 
 * @param numberEntries 
 * @param inputX 
 * @param inputY 
 * @param outputAngles 
 */
void cartesianToPolar(int numberEntries, float* inputX, float* inputY, float* outputAngles) {
    float pi = TMath::Pi();
    for (int i=0; i<numberEntries; i++) {
        // cout << "x = " << inputX[i] << " y = " << inputY[i] << endl;
        if (inputX[i] > 0 && inputY[i] >= 0) {
          outputAngles[i] = (TMath::ATan(inputY[i] / inputX[i]) * 180.0 / pi);
        //   cout << "case 1" << endl;
        } 
        else if (inputX[i] < 0) {
          outputAngles[i] = (TMath::ATan(inputY[i] / inputX[i])  + pi) * 180.0 / pi;
        //   cout << "case 2" << endl;
        } 
        else if (inputX[i] > 0 && inputY[i] < 0) {
          outputAngles[i] = (TMath::ATan(inputY[i] / inputX[i])  + 2*pi) * 180.0 / pi;
        //   cout << "case 3" << endl;
        } 
        else if (inputX[i] == 0 && inputY[i] > 0) {
          outputAngles[i] = 0.5 * 180.0;
        //   cout << "case 4" << endl;
        }
        else if (inputX[i] == 0 && inputY[i] < 0) {
          outputAngles[i] = -1.5* 180.0;
        //   cout << "case 5" << endl;
        } else cout << "no case" << endl;
        // printArray(outputAngles, numberEntries);
    }
}

/**
 * @brief Arctan func in polar angle [0, 360)
 * 
 * @param inputX 
 * @param inputY 
 * @return float 
 */
float arctan2(float inputX, float inputY) {
    float pi = TMath::Pi();
    float outputAngle;
    if (inputX > 0 && inputY >= 0) {
        outputAngle = (TMath::ATan(inputY / inputX) * 180.0 / pi);
    } 
    else if (inputX < 0) {
        outputAngle = (TMath::ATan(inputY / inputX)  + pi) * 180.0 / pi;
    } 
    else if (inputX > 0 && inputY < 0) {
        outputAngle = (TMath::ATan(inputY / inputX)  + 2*pi) * 180.0 / pi;
    } 
    else if (inputX == 0 && inputY > 0) {
        outputAngle = 0.5 * 180.0;
    }
    else if (inputX == 0 && inputY < 0) {
        outputAngle = -1.5* 180.0;
    }
    return outputAngle;
}

/**
 * @brief Convert list of Arctan func in polar angle [0, 360)
 * 
 * @param numberEntries 
 * @param inputX 
 * @param inputY 
 * @param outputAngles 
 */
void arctan2Arr(int numberEntries, float* inputX, float* inputY, float* outputAngles) {
    float pi = TMath::Pi();
    for (int i=0; i<numberEntries; i++) {
        // cout << "x = " << inputX[i] << " y = " << inputY[i] << endl;
        if (inputX[i] > 0 && inputY[i] >= 0) {
          outputAngles[i] = (TMath::ATan(inputY[i] / inputX[i]) * 180.0 / pi);
        //   cout << "case 1" << endl;
        } 
        else if (inputX[i] < 0) {
          outputAngles[i] = (TMath::ATan(inputY[i] / inputX[i])  + pi) * 180.0 / pi;
        //   cout << "case 2" << endl;
        } 
        else if (inputX[i] > 0 && inputY[i] < 0) {
          outputAngles[i] = (TMath::ATan(inputY[i] / inputX[i])  + 2*pi) * 180.0 / pi;
        //   cout << "case 3" << endl;
        } 
        else if (inputX[i] == 0 && inputY[i] > 0) {
          outputAngles[i] = 0.5 * 180.0;
        //   cout << "case 4" << endl;
        }
        else if (inputX[i] == 0 && inputY[i] < 0) {
          outputAngles[i] = -1.5* 180.0;
        //   cout << "case 5" << endl;
        } else cout << "no case" << endl;
        // printArray(outputAngles, numberEntries);
    }
}




/**
 * @brief projects angle into (-180, 180] range, where the zeroAngle is the new 0
 * 
 * @param inputAngle in degrees
 * @param zeroAngle new zero
 * @return float    translated angle
 */
float translateAngle(float inputAngle, float zeroAngle) {
    float newAngle = inputAngle - zeroAngle;
    if (newAngle <= -180.0) newAngle += 360.0;
    else if (newAngle > 180.0) newAngle -= 360;
    return newAngle;
}

/**
 * @brief projects angles into (-180, 180] range, where the zeroAngle is the new 0
 * 
 * @param inputAngles in degrees
 * @param sizeAngles size of angle array
 * @param zeroAngle new zero
 */
void translateAngle(float *inputAngles, int sizeAngles, float zeroAngle) {
    for (int i=0; i<sizeAngles; i++) {
        inputAngles[i] -= zeroAngle;
        if (inputAngles[i] <= -180.0) inputAngles[i] += 360.0;
        else if (inputAngles[i] > 180.0) inputAngles[i] -= 360.0;
    }
}