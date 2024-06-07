#pragma once
#include <upcxx/upcxx.hpp>
#include <vector>
#include "simcovcuda_driver.hpp"

const upcxx::memory_kind devMemKind = upcxx::memory_kind::cuda_device;

//short names for distributed objects
typedef upcxx::dist_object<upcxx::global_ptr<Voxel, devMemKind>> distVoxels;
typedef upcxx::dist_object<upcxx::global_ptr<testData, devMemKind>> distTestData;
typedef upcxx::device_allocator<upcxx::cuda_device> devAllocator;


struct GPUMemory {
    //dist objects
    distVoxels dvoxels;
    devAllocator dalloc;

    //global pointers
    upcxx::global_ptr<Voxel, devMemKind> voxels;
    upcxx::global_ptr<curandState_t, devMemKind> deviceRNG;

    //global tile mask pointers
    upcxx::global_ptr<int, devMemKind> tileMask;
    upcxx::global_ptr<int, devMemKind> newTileMask;
    upcxx::global_ptr<int, devMemKind> numActiveTiles;

    //global reduction pointers
    upcxx::global_ptr<Globals, devMemKind> globals;
    upcxx::global_ptr<Globals, devMemKind> blockGlobals;

    //host pointers for extracting GPU memory
    upcxx::global_ptr<Voxel> h_voxels;
    upcxx::global_ptr<Globals> h_globals;

    //for communication
    //Create a map that will hold all of the global pointers to the send/recv buffers
    std::unordered_map<std::string, upcxx::global_ptr<Voxel, devMemKind>> packedBuffers;
    std::unordered_map<std::string, distVoxels> distPackedBuffers;

    upcxx::global_ptr<Voxel, devMemKind> xpredRecv;
    upcxx::global_ptr<Voxel, devMemKind> xpredSend;
    upcxx::global_ptr<Voxel, devMemKind> xsuccRecv;
    upcxx::global_ptr<Voxel, devMemKind> xsuccSend;

    upcxx::global_ptr<Voxel, devMemKind> ypredRecv;
    upcxx::global_ptr<Voxel, devMemKind> ypredSend;
    upcxx::global_ptr<Voxel, devMemKind> ysuccRecv;
    upcxx::global_ptr<Voxel, devMemKind> ysuccSend;

    upcxx::global_ptr<Voxel, devMemKind> zpredRecv;
    upcxx::global_ptr<Voxel, devMemKind> zpredSend;
    upcxx::global_ptr<Voxel, devMemKind> zsuccRecv;
    upcxx::global_ptr<Voxel, devMemKind> zsuccSend;

    size_t bufferSize;
    upcxx::cuda_device gpuDevice;
    devAllocator gpuAlloc;
};

void initComms(const Options& opt, Comms& comms,
                const Dimensions& dims, const Dimensions& simDims, const Dimensions& tile, const Dimensions& tileDims);
void communicateVoxels(Options& opt, GPUMemory& gpumem, Comms& comms, Voxel* data, Dimensions& dims, Dimensions& tile, Dimensions& tileDims);
void setupGPUmemory(GPUMemory& gpumem, Options opt, Dimensions dims, Dimensions tileDims, Comms comms);
void getGPUmemory(Dimensions dims, GPUMemory& gpumem);
void setGPUmemory(Dimensions dims, GPUMemory& gpumem);
void getGPUReduction(GPUMemory& gpumem);
void cleanGPUmemory(GPUMemory& gpumem);
void printActives(GPUMemory& gpumem,  Dimensions tileDims);