#pragma once
#include <string>
#include <cuda.h>
#include <cuda_runtime.h>
#include <curand.h>
#include <curand_kernel.h>
#include "options.hpp"

/**
 * Simulation struct data
*/

enum class EpiType { HEALTHY=0, INCUBATING=1, EXPRESSING=2, APOPTOTIC=3, DEAD=4};
const std::string EpiTypeString[] = {"HEALTHY", "INCUBATING", "EXPRESSING", "APOPTOTIC", "DEAD"};

enum class VoxelType { NORMAL=0, GHOST=1, BOUNDARY=2};
const std::string VoxelTypeString[] = {"NORMAL", "GHOST", "BOUNDARY"};

struct Dimensions {
    int x,y,z,numPoints;
};

struct Voxel {
    //structure data
    int neighborhoodSize;
    int neighborhood[27];
    int x = 0, y = 0, z = 0;
    int myIdx = 0;
    VoxelType voxelType = VoxelType::NORMAL;

    //Concentration data
    float virions, nbVirions, inflammation, nbInflammation;

    //EpiCell data
    EpiType cellType = EpiType::HEALTHY;
    int incubationTimeSteps = -1;
    int expressingTimeSteps = -1;
    int apoptoticTimeSteps = -1;

    //TCell data
    bool hasTCell;
    int bindingPeriod = -1;
    int tissueTimeSteps = -1;
    float bindProb = 0.0f;

    //For Random Walking
    bool flip;
    int fx, fy, fz;
    int tieBreakValue; //tracking for the T Cell pointing to the target

    int winnerID; //tracking at the destination
    int winnerValue;

};

struct Comms {

    int rank;
    int totalRanks;

    int rankX, rankY, rankZ;
    int rankDimX, rankDimY, rankDimZ;

    int xPredecessor = -1;
    int yPredecessor = -1;
    int zPredecessor = -1;
    int xSuccessor = -1;
    int ySuccessor = -1;
    int zSuccessor = -1;

    int xSliceSize = 0;
    int ySliceSize = 0;
    int zSliceSize = 0;

};

/**
 * For communication testing
*/
struct testData {
    int x, y, z, rank;
};

// for device side reductions
struct Globals {
    unsigned long long int totalHealthy;
    unsigned long long int totalIncubating;
    unsigned long long int totalExpressing;
    unsigned long long int totalApoptotic;
    unsigned long long int totalDead;
    unsigned long long int numCirculatingTCells;
    unsigned long long int totalTCells;
    float totalVirions;
    float totalInflammation;
};

void initVoxels(Options opt, Comms comms, Dimensions dims, Dimensions tile, Dimensions tileDims,
                curandState_t* deviceRNG, Voxel* voxels);
                
void initializeInfection(float count, int x, int y, int z, Dimensions dims, Dimensions tile, Dimensions tileDims, Voxel* voxels, Comms comms);
void checkActiveTiles(Options opt, int timeStep,
                    Dimensions dims, Dimensions tile, Dimensions tileDims,
                    Voxel* voxels,
                    int* tileMask, int* newTileMask, int* numActiveTiles);
void updateConcentrations(Options opt, Dimensions dims, Dimensions tile,
                            Voxel* voxels,
                            int* tileMask, int* numActiveTiles);
void updateEpithelialCells(Options opt, Dimensions dims, Dimensions tile,
                                Voxel* voxels, curandState_t* deviceRNG,
                                int* tileMask, int* numActiveTiles);
void reduceGPU(Options opt,
            Dimensions dims, Voxel* voxels,
            Globals* blockGlobals, Globals* globals);

void packSlice(Options& opt,
                        Voxel* data, Voxel* packedBuffer,
                        Dimensions dims, Dimensions tile, Dimensions tileDims, Comms comms,
                        int sliceDim, int sliceIdx);
void unpackSlice(Options& opt,
                        Voxel* data, Voxel* packedBuffer,
                        Dimensions dims, Dimensions tile, Dimensions tileDims, Comms comms,
                        int sliceDim, int sliceIdx);

void ageTCells(Options opt, Dimensions dims, Dimensions tile,
                            Voxel* voxels,
                            int* tileMask, int* numActiveTiles);

void setupBind(Options opt, Dimensions dims, Dimensions tile, Dimensions tileDims,
                            Voxel* voxels, curandState_t* deviceRNG,
                            int* tileMask, int* numActiveTiles);

void executeBind(Options opt, Dimensions dims, Dimensions tile, Dimensions tileDims,
                            Voxel* voxels,
                            int* tileMask, int* numActiveTiles);

void setupMove(Options opt, Dimensions dims, Dimensions tile,
                            Voxel* voxels, curandState_t* deviceRNG,
                            int* tileMask, int* numActiveTiles);

void declareWinners(Options opt, Dimensions dims, Dimensions tile,
                                Voxel* voxels,
                                int* tileMask, int* numActiveTiles);

void flipTCells(Options opt, Dimensions dims, Dimensions tile,
                            Voxel* voxels,
                            int* tileMask, int* numActiveTiles);

void spawnTCells(int numXTravasing,
    Options opt, Comms comms, Dimensions dims, Dimensions tile, Dimensions tileDims,
    Voxel* voxels, Globals* globals, curandState_t* deviceRNG);

void generateTCells(Options opt, Globals* globals);

//Host rng
float drawValue(curandGenerator_t& gen);
uint drawUInt(curandGenerator_t& gen);
int getBetween(curandGenerator_t& gen, int min, int max);
unsigned int drawPoisson(curandGenerator_t& gen, float lambda);
bool trialSuccess(curandGenerator_t& gen, float p);
void setupHostRNG(curandGenerator_t& gen, Options& opt);

//Communication testing
void setCommTestData(Options& opt, testData* data, Dimensions dims, Dimensions tile, Dimensions tileDims, Comms comms);
void packSliceTestData(Options& opt,
                        testData* data, testData* packedBuffer,
                        Dimensions dims, Dimensions tile, Dimensions tileDims, Comms comms,
                        int sliceDim, int sliceIdx);
void unpackSliceTestData(Options& opt,
                        testData* data, testData* packedBuffer,
                        Dimensions dims, Dimensions tile, Dimensions tileDims, Comms comms,
                        int sliceDim, int sliceIdx);