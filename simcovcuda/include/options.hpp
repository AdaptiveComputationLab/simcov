#pragma once
#include <stdlib.h>
#include <string>
#include <unordered_map>
#include <utility>
#include <fstream>

/**
Options: struct containing parameter variables
**/
struct Options {
    int blocks = 8;
    int threadsPerBlock = 256;
    int dimX = 15000;
    int dimY = 15000;
    int dimZ = 1;
    int tileX = -1;
    int tileY = -1;
    int tileZ = -1;
    int rankDimX = 1;
    int rankDimY = 1;
    int rankDimZ = 1;
    int numInfections = 1;
    long long wholeLungX = 48000;
    long long wholeLungY = 40000;
    long long wholeLungZ = 20000;
    int numPoints = 15000*15000*1;
    int timeSteps = 33120;
    int initialInfections = 1000;
    int incubationPeriod = 480;
    int apoptosisPeriod = 180;
    int expressingPeriod = 900;
    float infectivity = 0.001;
    float virionProduction = 1.1;
    float virionClearance = 0.004;
    float virionDiffusion = 0.15;
    float inflammationProduction = 1.0;
    float inflammationDecay = 0.01;
    float inflammationDiffusion = 1.0;
    float minInflammation = 1e-6;
    float minVirions = 10e-10;
    float antibodyFactor = 1;
    int antibodyPeriod = 5760;
    int tcellGenerationRate = 105000;
    int tcellInitialDelay = 10080;
    int tcellVascularPeriod = 5760;
    int tcellTissuePeriod = 1440;
    int tcellBindingPeriod = 10;
    float maxBindingProb = 1;
    bool tcellsFollowGradient = false;
    bool commTest = false;
    int seed = -1;
    int samplePeriod = 0;
    int sampleResolution = 1;
    int reducePeriod = 0;
    int checkActivesRate = 1;
    std::string output{"simcov.log"};
    std::string gpuconfig{""};
    std::string config{""};
};

Options parseArgs(int argc, char** argv);