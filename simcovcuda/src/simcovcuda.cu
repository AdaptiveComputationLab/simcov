#include "simcovcuda_driver.hpp"
#include "simcovcuda.hpp"

#define RNG_t curandState_t

/**
 * Cuda error checking
*/
inline void lastError(const char* kernel){
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        printf("#!!!CUDA Error in %s: %s\n",kernel, cudaGetErrorString(err));
    }
}

/**
 * Device Functions
*/
__device__ bool inline d_isActive(Voxel* voxels, int i){
    if(voxels[i].virions > 0.0) return true;
    if(voxels[i].inflammation > 0.0) return true;
    if(voxels[i].cellType == EpiType::INCUBATING) return true;
    if(voxels[i].cellType == EpiType::EXPRESSING) return true;
    if(voxels[i].cellType == EpiType::APOPTOTIC) return true;
    if(voxels[i].voxelType == VoxelType::GHOST) return true;
    if(voxels[i].hasTCell) return true;
    return false;
}

__device__ float d_drawValue(RNG_t* state){
    float r = curand_uniform(state);
    return r;
}

__device__ uint d_drawUInt(RNG_t* state){
    unsigned int r = curand(state);
    return r;
}

__device__ int d_getBetween(RNG_t* state, int min, int max){
    if(min == max){
        return min;
    }
    uint r = d_drawUInt(state);
    return (int)((r)%((uint)max - (uint)min) + (uint)min);
}

__device__ unsigned int d_drawPoisson(RNG_t* state, float lambda){
    unsigned int r = curand_poisson(state, lambda);
    return r;
}

__device__ bool d_trialSuccess(RNG_t* state, float p){
    if(p > 1) return true;
    if(p < 0) return false;
    float roll = d_drawValue(state);
    if(roll < p) return true;
    return false;
}

__device__ void d_shuffle(RNG_t* state, int* a, int* b, int n){
    for(int i = 0; i < n; i++)
        b[i] = a[i];
    for(int i = n-1; i >= 1; i--){
        int j = d_getBetween(state, 0, i+1);
        int temp = b[j];
        b[j] = b[i];
        b[i] = temp;
    }
}

// Base functions for computing 3D and 1D coordinate transformations
__device__ inline void d_to3D(int i, Dimensions dims,
                int&dx, int&dy, int& dz){
    dx = i/(dims.z*dims.y);
    dy = (i % (dims.z*dims.y))/dims.z;
    dz = (i % (dims.z*dims.y))%dims.z;
}

__device__ inline int d_to1D(int x, int y, int z, Dimensions dims){
    return z + y*dims.z + x*(dims.y * dims.z);
}

__device__ inline bool d_inBoundaries(int x, int y, int z, Dimensions dims){
    if(x < 0 || x >= dims.x) return false;
    if(y < 0 || y >= dims.y) return false;
    if(z < 0 || z >= dims.z) return false;
    return true;
}

//functions for performing 3D to 1D coordinate transformations with
//tiles present
__device__ inline void d_to3D_tiled(int i,
                        Dimensions dims,
                        Dimensions tile,
                        Dimensions tileDims,
                        int& dx, int& dy, int& dz){
    
    //Determine the tile index
    int i_tile = i/tile.numPoints;

    //Compute the coordinates of the tile in 3D tile space
    int tx, ty, tz;
    d_to3D(i_tile, tileDims, tx, ty, tz);

    //get offset coords in voxel space
    //(i.e. the coords of the first voxel in this tile)
    int ox, oy, oz;
    ox = tx*tile.x;
    oy = ty*tile.y;
    oz = tz*tile.z;

    //Calculate the within tile coords
    int wi = i % tile.numPoints;
    int wx, wy, wz;
    d_to3D(wi, tile, wx, wy, wz);

    //result is offset + within tile coords
    dx = wx + ox;
    dy = wy + oy;
    dz = wz + oz;
}

__device__ inline int d_to1D_tiled(int x, int y, int z,
                    Dimensions dims,
                    Dimensions tile,
                    Dimensions tileDims){
    //get the tile coords in tile space
    int tx, ty, tz;
    tx = x/tile.x;
    ty = y/tile.y;
    tz = z/tile.z;
    //get the tile id
    int ti = d_to1D(tx, ty, tz, tileDims);

    //get the offset within the tile
    int ox, oy, oz;
    ox = x - tx*tile.x;
    oy = y - ty*tile.y;
    oz = z - tz*tile.z;
    int oi = d_to1D(ox, oy, oz, tile);

    //result is the total number of voxels up until this tile + the offest
    return ti*tile.numPoints + oi;
}

/**
 * Kernel definitions
**/
__global__ void k_initVoxels(Options opt, Comms comms, Dimensions dims, Dimensions tile, Dimensions tileDims,
                            curandState_t* deviceRNG,
                            Voxel* voxels){

    int start = blockIdx.x*blockDim.x + threadIdx.x;
    int maximum = dims.numPoints;
    int step = blockDim.x*gridDim.x;

    //initialize rng on the device
    curand_init(opt.seed + (unsigned int)(2048*comms.rank), start, 0, &deviceRNG[start]);
    RNG_t localState = deviceRNG[start];

    for(int i = start; i < maximum; i += step){

        //voxels
        d_to3D_tiled(i, dims, tile, tileDims, voxels[i].x, voxels[i].y, voxels[i].z);
        Voxel* v = &voxels[i];

        if(v->x == 0 && comms.xPredecessor != -1) v->voxelType = VoxelType::GHOST;
        if(v->x == dims.x - 1 && comms.xSuccessor != -1) v->voxelType = VoxelType::GHOST;
        if(v->y == 0 && comms.yPredecessor != -1) v->voxelType = VoxelType::GHOST;
        if(v->y == dims.y - 1 && comms.ySuccessor != -1) v->voxelType = VoxelType::GHOST;
        if(v->z == 0 && comms.zPredecessor != -1) v->voxelType = VoxelType::GHOST;
        if(v->z == dims.z - 1 && comms.zSuccessor != -1) v->voxelType = VoxelType::GHOST;

        v->neighborhoodSize = 0;
        for(int ii = -1; ii <= 1; ii++){
            for(int jj = -1; jj <= 1; jj++){
                for(int kk = -1; kk <= 1; kk++){
                    int id = d_to1D_tiled(v->x + ii, v->y + jj, v->z + kk, dims, tile, tileDims);
                    if(d_inBoundaries(v->x + ii, v->y + jj, v->z + kk, dims)){
                        if(id == i){
                            v->myIdx = v->neighborhoodSize;
                        }
                        v->neighborhood[v->neighborhoodSize] = id;
                        v->neighborhoodSize++;
                    }
                }
            }
        }

        //epicells
        v->cellType = EpiType::HEALTHY;
        v->incubationTimeSteps = d_drawPoisson(&localState, opt.incubationPeriod);
        v->expressingTimeSteps = d_drawPoisson(&localState, opt.expressingPeriod);
        v->apoptoticTimeSteps = d_drawPoisson(&localState, opt.apoptosisPeriod);

        //concentration data
        v->virions = 0.0f;
        v->nbVirions = 0.0f;
        v->inflammation = 0.0f;
        v->nbInflammation = 0.0f;

        //set tcell data
        v->hasTCell = false;
        v->flip = false;
        v->fx = 0;
        v->fy = 0;
        v->fz = 0;
        v->tissueTimeSteps = -1;
        v->bindingPeriod = -1;
        v->bindProb = 2.0f;
        v->tieBreakValue = -1;
        v->winnerID = -1;
        v->winnerValue = -1;
    }
    deviceRNG[start] = localState;
}

__global__ void k_resetActives(Dimensions tileDims, int* tileMask, int* newTileMask){
    int start = blockIdx.x*blockDim.x + threadIdx.x;
    int maximum = tileDims.numPoints;
    int step = blockDim.x*gridDim.x;

    for(int i = start; i < maximum; i+=step){
        newTileMask[i] = -1;
        tileMask[i] = -1;
    }
}

__global__ void k_checkActives(Dimensions dims, Dimensions tile, Dimensions tileDims,
                            Voxel* voxels,
                            int* tileMask, int* newTileMask, int* numActiveTiles){
    int start = blockIdx.x*blockDim.x + threadIdx.x;
    int maximum = dims.numPoints;
    int step = blockDim.x*gridDim.x; 

    for(int i = start; i < maximum; i+=step){
        if(d_isActive(voxels, i)){
            //activate this tile
            int tileID = i/tile.numPoints;
            //must be atomicCAS if threaded
            atomicCAS(&(newTileMask[tileID]), -1, tileID);
            // activate neighbors
            // get x,y,z
            int x,y,z;
            d_to3D(tileID, tileDims, x, y, z);

            for(int dx = -1; dx <= 1; dx++){
                for(int dy = -1; dy <= 1; dy++){
                    for(int dz = -1; dz <= 1; dz++){
                        int nx = x + dx;
                        int ny = y + dy;
                        int nz = z + dz;
                        if(d_inBoundaries(nx, ny, nz, tileDims)){
                            int ni = d_to1D(nx, ny, nz, tileDims);
                            if(ni != i){
                                //must be atomicCAS if threaded
                                atomicCAS(&(newTileMask[ni]), -1, ni);
                            }
                        }
                    }
                }
            }
        }
    }
}

__global__ void k_sortActives(Dimensions tileDims, int* tileMask, int* newTileMask, int* numActiveTiles){
    int start = blockIdx.x*blockDim.x + threadIdx.x;
    int maximum = tileDims.numPoints;
    int step = blockDim.x*gridDim.x; 

    for(int i = start; i < maximum; i+=step){

        //check if I'm active
        if(newTileMask[i] != -1){
            int prev = atomicAdd(numActiveTiles, 1);
            tileMask[prev] = i;
        }
    }
}

__global__ void k_setSingleValue(int* data, int i, int value){
    data[i] = value;
}

__global__ void k_accumulate(Options opt, Dimensions dims, Dimensions tile,
                            Voxel* voxels,
                            int* tileMask, int* numActiveTiles){

    int start = blockIdx.x*blockDim.x + threadIdx.x;
    int maximum = (*numActiveTiles)*tile.numPoints;
    int step = blockDim.x*gridDim.x;

    for(int idx = start; idx < maximum; idx += step){
        //stay within sim bounds even if tiles extend beyond
        if(idx >= dims.numPoints || idx < 0) continue;

        //figure out index into data array
        int tileID = tileMask[idx / (tile.numPoints)];
        if(tileID == -1) continue;
        int i = tileID*tile.numPoints + idx%tile.numPoints;

        // if(voxels[i].voxelType == VoxelType::GHOST) continue;

        Voxel* v = &voxels[i];
        v->nbVirions = 0.0f;
        v->nbInflammation = 0.0f;
        for(int j = 0; j < v->neighborhoodSize; j++){
            int nbID = v->neighborhood[j];
            if(nbID != i){
                v->nbVirions += voxels[nbID].virions;
                v->nbInflammation += voxels[nbID].inflammation;
            }
        }
    }
}

__global__ void  k_spread(Options opt, Dimensions dims, Dimensions tile,
                            Voxel* voxels,
                            int* tileMask, int* numActiveTiles){

    int start = blockIdx.x*blockDim.x + threadIdx.x;
    int maximum = (*numActiveTiles)*tile.numPoints;
    int step = blockDim.x*gridDim.x;

    for(int idx = start; idx < maximum; idx += step){
        //stay within sim bounds even if tiles extend beyond
        if(idx >= dims.numPoints || idx < 0) continue;

        //figure out index into data array
        int tileID = tileMask[idx / (tile.numPoints)];
        if(tileID == -1) continue;
        int i = tileID*tile.numPoints + idx%tile.numPoints;

        // if(voxels[i].voxelType == VoxelType::GHOST) continue;

        Voxel* v = &voxels[i];

        //virion
        float virionsDiffused = v->virions*opt.virionDiffusion;
        float virionsLeft = v->virions - virionsDiffused;
        float avgNBVirions = (virionsDiffused + v->nbVirions*opt.virionDiffusion)/(v->neighborhoodSize);

        v->virions = virionsLeft + avgNBVirions;

        v->virions = (1.0 - opt.virionClearance)*v->virions;
        if(v->virions < opt.minVirions) v->virions = 0.0;
        v->nbVirions = 0.0f;

        //inflammation
        float inflammationDiffused = v->inflammation*opt.inflammationDiffusion;
        float inflammationLeft = v->inflammation - inflammationDiffused;
        float avgNBInflammation = (inflammationDiffused + v->nbInflammation*opt.inflammationDiffusion)/(v->neighborhoodSize);

        v->inflammation = inflammationLeft + avgNBInflammation;
        v->inflammation = (1.0 - opt.inflammationDecay)*v->inflammation;
        if(v->inflammation < opt.minInflammation) v->inflammation = 0.0f;
        v->nbInflammation = 0.0f;

    }
}

__global__ void k_updateEpiCells(Options opt, Dimensions dims, Dimensions tile,
                                Voxel* voxels, curandState_t* deviceRNG,
                                int* tileMask, int* numActiveTiles){

    int start = blockIdx.x*blockDim.x + threadIdx.x;
    int maximum = (*numActiveTiles)*tile.numPoints;
    int step = blockDim.x*gridDim.x;
    RNG_t localState;

    for(int idx = start; idx < maximum; idx += step){
        //stay within sim bounds even if tiles extend beyond
        if(idx >= dims.numPoints || idx < 0) continue;

        //figure out index into data array
        int tileID = tileMask[idx / (tile.numPoints)];
        if(tileID == -1) continue;
        int i = tileID*tile.numPoints + idx%tile.numPoints;

        localState = deviceRNG[i % (opt.threadsPerBlock*opt.blocks)];

        if(voxels[i].voxelType == VoxelType::GHOST) continue;

        bool produceVirions = false;
        switch(voxels[i].cellType) {
            case EpiType::HEALTHY:
                if(voxels[i].virions > 0) {
                    if(d_trialSuccess(&localState, voxels[i].virions*opt.infectivity)) {
                        voxels[i].cellType = EpiType::INCUBATING;
                    }
                }
                break;
            case EpiType::INCUBATING:
                voxels[i].incubationTimeSteps--;
                if(voxels[i].incubationTimeSteps <= 0) {
                    voxels[i].cellType = EpiType::EXPRESSING;
                }
                break;
            case EpiType::EXPRESSING:
                voxels[i].expressingTimeSteps--;
                if(voxels[i].expressingTimeSteps <= 0) {
                    voxels[i].cellType = EpiType::DEAD;
                } else {
                    produceVirions = true;
                }
                break;
            case EpiType::APOPTOTIC:
                voxels[i].apoptoticTimeSteps--;
                if(voxels[i].apoptoticTimeSteps <= 0) {
                    voxels[i].cellType = EpiType::DEAD;
                } else if (voxels[i].incubationTimeSteps<=0) {
                    produceVirions = true;
                }
                break;
            default: break;
        }
        if(produceVirions) {
            voxels[i].virions = opt.virionProduction + voxels[i].virions;
            voxels[i].inflammation = opt.inflammationProduction + voxels[i].inflammation;
            if(voxels[i].inflammation > 1.0) {
                voxels[i].inflammation = 1.0;
            }
        }
        deviceRNG[i % (opt.threadsPerBlock*opt.blocks)] = localState;
    }
    
}

__global__ void k_reduce(Dimensions dims,
                        Voxel* voxels,
                        Globals* blockGlobals,
                        int offset, int threads){
    extern __shared__ float s[];
    float *sharedVirions = s;
    float *sharedInflammation = (float*)&sharedVirions[threads];
    int *sharedTCells = (int*)&sharedInflammation[threads];
    int *sharedHealthy = (int*)&sharedTCells[threads];
    int *sharedIncubating = (int*)&sharedHealthy[threads];
    int *sharedExpressing = (int*)&sharedIncubating[threads];
    int *sharedApoptotic = (int*)&sharedExpressing[threads];
    int *sharedDead = (int*)&sharedApoptotic[threads];

    unsigned int tid = threadIdx.x;
    unsigned int i = blockIdx.x*blockDim.x + threadIdx.x + offset;

    sharedVirions[tid] = 0;
    sharedInflammation[tid] = 0;
    sharedTCells[tid] = 0;
    sharedHealthy[tid] = 0;
    sharedIncubating[tid] = 0;
    sharedExpressing[tid] = 0;
    sharedApoptotic[tid] = 0;
    sharedDead[tid] = 0;

    if(i < dims.numPoints && voxels[i].voxelType != VoxelType::GHOST){
        sharedVirions[tid] = voxels[i].virions;
        sharedInflammation[tid] = voxels[i].inflammation;
        if(voxels[i].hasTCell) sharedTCells[tid] = 1;
        switch(voxels[i].cellType){
            case EpiType::HEALTHY:
                sharedHealthy[tid] = 1;
                break;
            case EpiType::INCUBATING:
                sharedIncubating[tid] = 1;
                break;
            case EpiType::EXPRESSING:
                sharedExpressing[tid] = 1;
                break;
            case EpiType::APOPTOTIC:
                sharedApoptotic[tid] = 1;
                break;
            case EpiType::DEAD:
                sharedDead[tid] = 1;
                break;
            default:
                break;
        }
    }
    __syncthreads();

    for(unsigned int s=blockDim.x/2; s > 0; s >>=1){
        int index = tid;
        if(index < s){
            sharedVirions[index] += sharedVirions[index + s];       
            sharedInflammation[index] += sharedInflammation[index + s];
            sharedTCells[index] += sharedTCells[index + s];
            sharedHealthy[index] += sharedHealthy[index + s];       
            sharedIncubating[index] += sharedIncubating[index + s];       
            sharedExpressing[index] += sharedExpressing[index + s];       
            sharedApoptotic[index] += sharedApoptotic[index + s];       
            sharedDead[index] += sharedDead[index + s];       
        }
        __syncthreads();
    }

    //write back to global memory
    if (tid == 0) {
        blockGlobals[blockIdx.x].totalVirions = sharedVirions[0];
        blockGlobals[blockIdx.x].totalInflammation = sharedInflammation[0];
        blockGlobals[blockIdx.x].totalTCells = sharedTCells[0];
        blockGlobals[blockIdx.x].totalHealthy = sharedHealthy[0];
        blockGlobals[blockIdx.x].totalIncubating = sharedIncubating[0];
        blockGlobals[blockIdx.x].totalExpressing = sharedExpressing[0];
        blockGlobals[blockIdx.x].totalApoptotic = sharedApoptotic[0];
        blockGlobals[blockIdx.x].totalDead = sharedDead[0];
    }
}

__global__ void k_reduceBlocks(Globals* blockGlobals, Globals* globals){
    int i = threadIdx.x;
    atomicAdd(&(globals->totalVirions), blockGlobals[i].totalVirions);
    atomicAdd(&(globals->totalInflammation), blockGlobals[i].totalInflammation);
    atomicAdd(&(globals->totalHealthy), blockGlobals[i].totalHealthy);
    atomicAdd(&(globals->totalIncubating), blockGlobals[i].totalIncubating);
    atomicAdd(&(globals->totalExpressing), blockGlobals[i].totalExpressing);
    atomicAdd(&(globals->totalApoptotic), blockGlobals[i].totalApoptotic);
    atomicAdd(&(globals->totalDead), blockGlobals[i].totalDead);
    atomicAdd(&(globals->totalTCells), blockGlobals[i].totalTCells);
}

__global__ void k_resetGlobals(Globals* globals){
    globals->totalVirions = 0.0;
    globals->totalInflammation = 0.0;
    globals->totalHealthy = 0.0;
    globals->totalIncubating = 0.0;
    globals->totalExpressing = 0.0;
    globals->totalApoptotic = 0.0;
    globals->totalDead = 0.0;
    globals->totalTCells = 0.0;
}

__global__ void k_initializeInfection(float count, int x, int y, int z, Dimensions dims, Dimensions tile, Dimensions tileDims, Voxel* voxels, Comms comms){
    //check if I own this tile
    int minX, maxX;
    int minY, maxY;
    int minZ, maxZ;
    minX = comms.rankX*comms.rankDimX;
    maxX = minX + comms.rankDimX;
    minY = comms.rankY*comms.rankDimY;
    maxY = minY + comms.rankDimY;
    minZ = comms.rankZ*comms.rankDimZ;
    maxZ = minZ + comms.rankDimZ;

    if(x < minX || x >= maxX) return;
    if(y < minY || y >= maxY) return;
    if(z < minZ || z >= maxZ) return;

    int ox = x - minX;
    int oy = y - minY;
    int oz = z - minZ;
    if(comms.xPredecessor != -1) ox++;
    if(comms.yPredecessor != -1) oy++;
    if(comms.zPredecessor != -1) oz++;

    int i = d_to1D_tiled(ox, oy, oz, dims, tile, tileDims);
    voxels[i].virions = count;
    printf("### Rank %d initialized infection %f at %d = (%d,%d,%d) locally at (%d,%d,%d)\n", comms.rank, count, i, x,y,z, ox,oy,oz);
    voxels[i].cellType = EpiType::INCUBATING;
}

__global__ void k_ageTCells(Dimensions dims, Dimensions tile,
                            Voxel* voxels,
                            int* tileMask, int* numActiveTiles){

    int start = blockIdx.x*blockDim.x + threadIdx.x;
    int maximum = (*numActiveTiles)*tile.numPoints;
    int step = blockDim.x*gridDim.x;

    for(int idx = start; idx < maximum; idx += step){
        //stay within sim bounds even if tiles extend beyond
        if(idx >= dims.numPoints || idx < 0) continue;

        //figure out index into data array
        int tileid = tileMask[idx / (tile.numPoints)];
        if(tileid == -1) continue;
        int i = tileid*tile.numPoints + idx%tile.numPoints;

        if(voxels[i].voxelType == VoxelType::GHOST) continue;

        voxels[i].tissueTimeSteps--;
        if(voxels[i].tissueTimeSteps <= 0) {
            //kill the tcell
            voxels[i].tissueTimeSteps = -1;
            voxels[i].bindingPeriod = -1;
            voxels[i].hasTCell = false;
            voxels[i].flip = false;
        }
    }
}

__global__ void k_setupBind(Options opt, Dimensions dims, Dimensions tile, Dimensions tileDims,
                            Voxel* voxels, curandState_t* deviceRNG,
                            int* tileMask, int* numActiveTiles){
    int start = blockIdx.x*blockDim.x + threadIdx.x;
    int maximum = (*numActiveTiles)*tile.numPoints;
    int step = blockDim.x*gridDim.x;
    int stateIDX = start%(opt.threadsPerBlock*opt.blocks);
    RNG_t localState = deviceRNG[stateIDX];

    for(int idx = start; idx < maximum; idx += step){
        //stay within sim bounds even if tiles extend beyond
        if(idx >= dims.numPoints || idx < 0) continue;

        //figure out index into data array
        int tileID = tileMask[idx / (tile.numPoints)];
        if(tileID == -1) continue;
        int i = tileID*tile.numPoints + idx%tile.numPoints;

        if(voxels[i].voxelType == VoxelType::GHOST) continue;

        Voxel* v = &voxels[i];
        v->bindProb = 2.0f; //reset

        if(v->bindingPeriod != -1){
            v->bindingPeriod--;
            if(v->bindingPeriod <= 0){
                v->bindingPeriod = -1;
                v->fx = 0;
                v->fy = 0;
                v->fz = 0;
                v->bindProb = 2.0f; //ensures won't bind this turn
            }
            continue;
        }
        if(!v->hasTCell) continue;

        v->fx = 0;
        v->fy = 0;
        v->fz = 0;
        bool foundTarget = false;
        int dx = 0, dy = 0, dz = 0;

        //shuffle neighborhood
        int* shuffledNeighbors = new int[v->neighborhoodSize];
        d_shuffle(&localState, v->neighborhood, shuffledNeighbors, v->neighborhoodSize);

        for(int ii = 0; ii < v->neighborhoodSize; ii++){
            int ni = shuffledNeighbors[ii];
            Voxel vn = voxels[ni];
            v->bindProb = d_drawValue(&localState);
            dx = vn.x - v->x;
            dy = vn.y - v->y;
            dz = vn.z - v->z;

            if(voxels[ni].cellType == EpiType::HEALTHY || voxels[ni].cellType == EpiType::DEAD){
                continue;
            }

            double scaling = 1.0 - (double)voxels[ni].incubationTimeSteps/opt.incubationPeriod;

            if(voxels[ni].cellType == EpiType::EXPRESSING || voxels[ni].cellType == EpiType::APOPTOTIC){
                scaling = opt.maxBindingProb;
            }

            if(scaling < 0) scaling = 0;
            double prob = opt.maxBindingProb*scaling;
            if(prob >= opt.maxBindingProb) prob = opt.maxBindingProb;

            if(v->bindProb < prob){
                foundTarget = true;
                v->bindProb = -1.0; //i will bind
                break;
            }
        }
        delete [] shuffledNeighbors;

        if(foundTarget){
            v->fx = dx;
            v->fy = dy;
            v->fz = dz;
        }

        if(dims.x == 1) v->fx = 0;
        if(dims.y == 1) v->fy = 0;
        if(dims.z == 1) v->fz = 0;
    }
    deviceRNG[stateIDX] = localState;
}

__device__ EpiType atomicLoad(const EpiType* addr){
    const volatile EpiType *vaddr = addr;
    __threadfence();
    const EpiType value = *vaddr;
    __threadfence();
    return value;
}

__device__ void atomicStore(EpiType* addr, EpiType value){
    volatile EpiType *vaddr = addr;
    __threadfence();
    *vaddr = value;
}

__global__ void k_executeBind(Options opt, Dimensions dims, Dimensions tile, Dimensions tileDims,
                            Voxel* voxels,
                            int* tileMask, int* numActiveTiles){
    int start = blockIdx.x*blockDim.x + threadIdx.x;
    int maximum = (*numActiveTiles)*tile.numPoints;
    int step = blockDim.x*gridDim.x;

    for(int idx = start; idx < maximum; idx += step){
        //stay within sim bounds even if tiles extend beyond
        if(idx >= dims.numPoints || idx < 0) continue;

        //figure out index into data array
        int tileID = tileMask[idx / (tile.numPoints)];
        if(tileID == -1) continue;
        int i = tileID*tile.numPoints + idx%tile.numPoints;
        Voxel* v = &voxels[i];
        int lx, ly, lz;
        d_to3D_tiled(i, dims, tile, tileDims, lx, ly, lz);
        if(v->bindingPeriod == -1 && v->hasTCell){
            if(d_inBoundaries(lx + v->fx, ly + v->fy, lz + v->fz, dims)){
                int id = d_to1D_tiled(lx + v->fx, ly + v->fy, lz + v->fz, dims, tile, tileDims);
                if(v->bindProb < 0.0){
                    v->bindingPeriod = opt.tcellBindingPeriod;
                    atomicStore(&voxels[id].cellType, EpiType::APOPTOTIC);
                }
            }
            v->bindProb = 0.0;

            //dont move if bound
            if(v->bindingPeriod != -1){
                v->fx = 0;
                v->fy = 0;
                v->fz = 0;
            }

        }
    }
}

__global__ void k_setupMove(Dimensions dims, Dimensions tile,
                            Voxel* voxels, curandState_t* deviceRNG,
                            int* tileMask, int* numActiveTiles){
    int start = blockIdx.x*blockDim.x + threadIdx.x;
    int maximum = (*numActiveTiles)*tile.numPoints;
    int step = blockDim.x*gridDim.x;

    RNG_t localState = deviceRNG[start];

    for(int idx = start; idx < maximum; idx += step){
        //stay within sim bounds even if tiles extend beyond
        if(idx >= dims.numPoints || idx < 0) continue;

        //figure out index into data array
        int tileID = tileMask[idx / (tile.numPoints)];
        if(tileID == -1) continue;
        int i = tileID*tile.numPoints + idx%tile.numPoints;

        if(voxels[i].voxelType == VoxelType::GHOST) continue;

        Voxel* v = &voxels[i];

        if(v->hasTCell && v->bindingPeriod == -1){
            //choose a random direction to face
            int nbIDX = d_getBetween(&localState, 0, v->neighborhoodSize-1);
            if(nbIDX >= v->myIdx) nbIDX += 1;
            int nbID = v->neighborhood[nbIDX];
            v->fx = voxels[nbID].x - v->x;
            v->fy = voxels[nbID].y - v->y;
            v->fz = voxels[nbID].z - v->z;

            v->tieBreakValue = d_getBetween(&localState, 0, 1000000);
            if(v->tieBreakValue > voxels[nbID].winnerValue){
                atomicCAS(&voxels[nbID].winnerValue,
                        voxels[nbID].winnerValue,
                        v->tieBreakValue);
            } else {
                v->fx = 0;
                v->fy = 0;
                v->fz = 0;
                v->tieBreakValue = -1;
            }
        }
    }
    deviceRNG[start] = localState;
}

__global__ void k_declareWinners(Dimensions dims, Dimensions tile,
                                Voxel* voxels,
                                int* tileMask, int* numActiveTiles){
    int start = blockIdx.x*blockDim.x + threadIdx.x;
    int maximum = (*numActiveTiles)*tile.numPoints;
    int step = blockDim.x*gridDim.x;

    for(int idx = start; idx < maximum; idx += step){
        //stay within sim bounds even if tiles extend beyond
        if(idx >= dims.numPoints || idx < 0) continue;

        //figure out index into data array
        int tileID = tileMask[idx / (tile.numPoints)];
        if(tileID == -1) continue;
        int i = tileID*tile.numPoints + idx%tile.numPoints;

        Voxel* v = &voxels[i];
        if(v->tieBreakValue != -1){
            for(int ii = 0; ii < v->neighborhoodSize; ii++){
                int nbID = v->neighborhood[ii];
                if(voxels[nbID].winnerValue == v->tieBreakValue){
                    if(!voxels[nbID].hasTCell){
                        voxels[nbID].winnerID = i;
                        v->flip = true;
                        voxels[nbID].flip = true;
                        voxels[nbID].bindingPeriod = v->bindingPeriod;
                        voxels[nbID].tissueTimeSteps = v->tissueTimeSteps;
                        voxels[nbID].bindProb = v->bindProb;
                    }
                    break;
                }
            }
        }

    }
}

__global__ void k_flipTCells(Dimensions dims, Dimensions tile,
                            Voxel* voxels,
                            int* tileMask, int* numActiveTiles){
    int start = blockIdx.x*blockDim.x + threadIdx.x;
    int maximum = (*numActiveTiles)*tile.numPoints;
    int step = blockDim.x*gridDim.x;

    for(int idx = start; idx < maximum; idx += step){
        //stay within sim bounds even if tiles extend beyond
        if(idx >= dims.numPoints || idx < 0) continue;

        //figure out index into data array
        int tileID = tileMask[idx / (tile.numPoints)];
        if(tileID == -1) continue;
        int i = tileID*tile.numPoints + idx%tile.numPoints;

        Voxel* v = &voxels[i];

        if(v->flip){
            if(v->hasTCell){
                //unset tcell values
                v->bindingPeriod = -1;
                v->tissueTimeSteps = -1;
                v->bindProb = 2.0f;
                v->hasTCell = false;
            } else {
                v->hasTCell = true;
            }
        }
        //reset random walk values
        v->tieBreakValue = -1;
        v->winnerValue = -1;
        v->winnerID = -1;
        v->flip = false;
    }
}

__global__ void k_spawnTCells(int numXTravasing,
                            Options opt, Comms comms, Dimensions dims, Dimensions tile, Dimensions tileDims,
                            Voxel* voxels, Globals* globals, curandState_t* deviceRNG){
    int start = blockIdx.x*blockDim.x + threadIdx.x;
    int maximum = numXTravasing;
    int step = blockDim.x*gridDim.x;

    int stateIDX = start%(opt.threadsPerBlock*opt.blocks);
    RNG_t localState = deviceRNG[stateIDX];

    for(int i = start; i < maximum; i += step) {

        int minX = 0, maxX = dims.x;
        int minY = 0, maxY = dims.y;
        int minZ = 0, maxZ = dims.z;

        //deal with ghost voxels
        if(comms.xPredecessor != -1) minX += 1;
        if(comms.xSuccessor != -1) maxX -= 1;
        if(comms.yPredecessor != -1) minY += 1;
        if(comms.ySuccessor != -1) maxY -= 1;
        if(comms.zPredecessor != -1) minZ += 1;
        if(comms.zSuccessor != -1) maxZ -= 1;

        int lx = d_getBetween(&localState, minX, maxX);
        int ly = d_getBetween(&localState, minY, maxY);
        int lz = d_getBetween(&localState, minZ, maxZ);

        int id = d_to1D_tiled(lx,ly,lz,dims, tile, tileDims);
        int lifeTime = d_drawPoisson(&localState, opt.tcellTissuePeriod);

        if(voxels[id].inflammation < opt.minInflammation) continue;
        int prev_id = atomicCAS(&voxels[id].tissueTimeSteps, -1, lifeTime);
        if(prev_id != -1) continue;
        voxels[id].hasTCell = true;
        atomicAdd(&globals->numCirculatingTCells, (unsigned long long int)(-1));
    }
    deviceRNG[stateIDX] = localState;
}

__global__ void k_packSlice(Options& opt,
                        Voxel* data, Voxel* packedBuffer,
                        Dimensions dims, Dimensions tile, Dimensions tileDims, Comms comms,
                        int sliceDim, int sliceIdx){

    int start = blockIdx.x*blockDim.x + threadIdx.x;
    int maximum = 0;
    int step = blockDim.x*gridDim.x;

    //X
    if(sliceDim == 0){
        maximum = comms.xSliceSize;
        for(int i = start; i < maximum; i += step){
            int x, y, z;
            x = sliceIdx;
            y = i / dims.z;
            z = i % dims.z;
            int idx = d_to1D_tiled(x, y, z, dims, tile, tileDims);

            packedBuffer[i] = data[idx];
        }
    }

    //Y
    if(sliceDim == 1){
        maximum = comms.ySliceSize;
        for(int i = start; i < maximum; i += step){
            int x, y, z;
            x = i % dims.x;
            y = sliceIdx;
            z = i / dims.x;
            int idx = d_to1D_tiled(x, y, z, dims, tile, tileDims);
            packedBuffer[i] = data[idx];
        }
    }

    //Z
    if(sliceDim == 2){
        maximum = comms.zSliceSize;
        for(int i = start; i < maximum; i += step){
            int x, y, z;
            x = i / dims.y;
            y = i % dims.y;
            z = sliceIdx;
            int idx = d_to1D_tiled(x, y, z, dims, tile, tileDims);
            packedBuffer[i] = data[idx];
        }
    }
    
}

__global__ void k_unpackSlice(Options& opt,
                        Voxel* data, Voxel* packedBuffer,
                        Dimensions dims, Dimensions tile, Dimensions tileDims, Comms comms,
                        int sliceDim, int sliceIdx){

    int start = blockIdx.x*blockDim.x + threadIdx.x;
    int maximum = 0;
    int step = blockDim.x*gridDim.x;

    //X
    if(sliceDim == 0){
        maximum = comms.xSliceSize;
        for(int i = start; i < maximum; i += step){
            int x, y, z;
            x = sliceIdx;
            y = i / dims.z;
            z = i % dims.z;
            int idx = d_to1D_tiled(x, y, z, dims, tile, tileDims);
            //preserve structure data
            int neighborhoodSize = data[idx].neighborhoodSize;
            int neighborhood[27] = {0};
            for(int ii = 0; ii < 27; ii++){
                neighborhood[ii] = data[idx].neighborhood[ii];
            }
            int myIdx = data[idx].myIdx;
            VoxelType voxelType = data[idx].voxelType;
            int oldx = data[idx].x;
            int oldy = data[idx].y;
            int oldz = data[idx].z;

            //copy
            data[idx] = packedBuffer[i];
            
            //reset structure data
            data[idx].neighborhoodSize = neighborhoodSize;
            for(int ii = 0; ii < 27; ii++){
                data[idx].neighborhood[ii] = neighborhood[ii];
            }
            data[idx].myIdx = myIdx;
            data[idx].voxelType = voxelType;
            data[idx].x = oldx;
            data[idx].y = oldy;
            data[idx].z = oldz;
        }
    }

    //Y
    if(sliceDim == 1){
        maximum = comms.ySliceSize;
        for(int i = start; i < maximum; i += step){
            int x, y, z;
            x = i % dims.x;
            y = sliceIdx;
            z = i / dims.x;
            int idx = d_to1D_tiled(x, y, z, dims, tile, tileDims);
            //preserve structure data
            int neighborhoodSize = data[idx].neighborhoodSize;
            int neighborhood[27] = {0};
            for(int ii = 0; ii < 27; ii++){
                neighborhood[ii] = data[idx].neighborhood[ii];
            }
            int myIdx = data[idx].myIdx;
            VoxelType voxelType = data[idx].voxelType;
            int oldx = data[idx].x;
            int oldy = data[idx].y;
            int oldz = data[idx].z;

            //copy
            data[idx] = packedBuffer[i];
            
            //reset structure data
            data[idx].neighborhoodSize = neighborhoodSize;
            for(int ii = 0; ii < 27; ii++){
                data[idx].neighborhood[ii] = neighborhood[ii];
            }
            data[idx].myIdx = myIdx;
            data[idx].voxelType = voxelType;
            data[idx].x = oldx;
            data[idx].y = oldy;
            data[idx].z = oldz;
        }
    }

    //Z
    if(sliceDim == 2){
        maximum = comms.zSliceSize;
        for(int i = start; i < maximum; i += step){
            int x, y, z;
            x = i / dims.y;
            y = i % dims.y;
            z = sliceIdx;
            int idx = d_to1D_tiled(x, y, z, dims, tile, tileDims);
           //preserve structure data
            int neighborhoodSize = data[idx].neighborhoodSize;
            int neighborhood[27] = {0};
            for(int ii = 0; ii < 27; ii++){
                neighborhood[ii] = data[idx].neighborhood[ii];
            }
            int myIdx = data[idx].myIdx;
            VoxelType voxelType = data[idx].voxelType;
            int oldx = data[idx].x;
            int oldy = data[idx].y;
            int oldz = data[idx].z;

            //copy
            data[idx] = packedBuffer[i];
            
            //reset structure data
            data[idx].neighborhoodSize = neighborhoodSize;
            for(int ii = 0; ii < 27; ii++){
                data[idx].neighborhood[ii] = neighborhood[ii];
            }
            data[idx].myIdx = myIdx;
            data[idx].voxelType = voxelType;
            data[idx].x = oldx;
            data[idx].y = oldy;
            data[idx].z = oldz;
        }
    }
    
}


__global__ void k_generateTCells(Options opt, Globals* globals){
    globals->numCirculatingTCells += (unsigned long long int)(opt.tcellGenerationRate);
}

/**
 * Kernel Wrappers
*/
void initVoxels(Options opt, Comms comms, Dimensions dims, Dimensions tile, Dimensions tileDims,
                curandState_t* deviceRNG, Voxel* voxels){
    k_initVoxels<<<opt.blocks, opt.threadsPerBlock>>>(opt, comms, dims, tile, tileDims,
                                                     deviceRNG, voxels);
    lastError("initVoxels");
    cudaDeviceSynchronize();
}

void checkActiveTiles(Options opt, int timeStep,
                    Dimensions dims, Dimensions tile, Dimensions tileDims,
                    Voxel* voxels,
                    int* tileMask, int* newTileMask, int* numActiveTiles){
    if(timeStep % opt.checkActivesRate == 0){
        k_resetActives<<<opt.blocks, opt.threadsPerBlock>>>(tileDims, tileMask, newTileMask);
        lastError("resetActives");
        
        k_setSingleValue<<<1,1>>>(numActiveTiles, 0, 0); //Might be a smarter way to do this but it isn't very expensive. Faster than copying anyway
        lastError("setSingleValue");

        k_checkActives<<<opt.blocks, opt.threadsPerBlock>>>(dims, tile, tileDims,
                                                            voxels,
                                                            tileMask, newTileMask, numActiveTiles);
        lastError("checkActives");

        k_sortActives<<<opt.blocks, opt.threadsPerBlock>>>(tileDims, tileMask, newTileMask, numActiveTiles);
        lastError("sortActives");
    }
}

void updateConcentrations(Options opt, Dimensions dims, Dimensions tile,
                            Voxel* voxels,
                            int* tileMask, int* numActiveTiles){
    k_accumulate<<<opt.blocks, opt.threadsPerBlock>>>(opt, dims, tile,
                                                    voxels,
                                                    tileMask, numActiveTiles);
    lastError("update concentrations");

    k_spread<<<opt.blocks, opt.threadsPerBlock>>>(opt, dims, tile,
                                                voxels,
                                                tileMask, numActiveTiles);
    lastError("spread");
}
void updateEpithelialCells(Options opt, Dimensions dims, Dimensions tile,
                                Voxel* voxels, curandState_t* deviceRNG,
                                int* tileMask, int* numActiveTiles){
    k_updateEpiCells<<<opt.blocks, opt.threadsPerBlock>>>(opt, dims, tile,
                                voxels, deviceRNG,
                                tileMask, numActiveTiles);
    lastError("updateEpithelialCells");
}

void reduceGPU(Options opt,
            Dimensions dims, Voxel* voxels,
            Globals* blockGlobals, Globals* globals){
    int threads = opt.threadsPerBlock;
    int blocks = opt.blocks;
    int shmem_size = threads*(2*sizeof(float) + 6*sizeof(int));
    k_resetGlobals<<<1,1>>>(globals);
    for(int i = 0; i < dims.numPoints; i+=(threads*blocks))
    {
        k_reduce<<<blocks, threads, shmem_size>>>(dims, voxels, blockGlobals, i, threads);
        lastError("reduce");
        k_reduceBlocks<<<1, blocks>>>(blockGlobals, globals);
        lastError("reduceBlocks");
    }
    cudaDeviceSynchronize();
}

void initializeInfection(float count, int x, int y, int z, Dimensions dims, Dimensions tile, Dimensions tileDims, Voxel* voxels, Comms comms){
    k_initializeInfection<<<1,1>>>(count, x, y, z, dims, tile, tileDims, voxels, comms);
    lastError("initializeInfection");
}

void ageTCells(Options opt, Dimensions dims, Dimensions tile,
                            Voxel* voxels,
                            int* tileMask, int* numActiveTiles){
    k_ageTCells<<<opt.blocks,opt.threadsPerBlock>>>(dims, tile, voxels, tileMask, numActiveTiles);
    lastError("ageTCells");
}

void setupBind(Options opt, Dimensions dims, Dimensions tile, Dimensions tileDims,
                            Voxel* voxels, curandState_t* deviceRNG,
                            int* tileMask, int* numActiveTiles){
    k_setupBind<<<opt.blocks, opt.threadsPerBlock>>>(opt, dims, tile, tileDims,
                                                    voxels, deviceRNG,
                                                    tileMask, numActiveTiles);
}

void executeBind(Options opt, Dimensions dims, Dimensions tile, Dimensions tileDims,
                            Voxel* voxels,
                            int* tileMask, int* numActiveTiles){
    k_executeBind<<<opt.blocks, opt.threadsPerBlock>>>(opt, dims, tile, tileDims,
                                                    voxels,
                                                    tileMask, numActiveTiles);
    lastError("executeBind");
}

void setupMove(Options opt, Dimensions dims, Dimensions tile,
                            Voxel* voxels, curandState_t* deviceRNG,
                            int* tileMask, int* numActiveTiles){
    k_setupMove<<<opt.blocks, opt.threadsPerBlock>>>(dims, tile, voxels, deviceRNG, tileMask, numActiveTiles);
    lastError("setupMove");
}

void declareWinners(Options opt, Dimensions dims, Dimensions tile,
                                Voxel* voxels,
                                int* tileMask, int* numActiveTiles){
    k_declareWinners<<<opt.blocks, opt.threadsPerBlock>>>(dims, tile, voxels, tileMask, numActiveTiles);
    lastError("declare winners");
}

void flipTCells(Options opt, Dimensions dims, Dimensions tile, Voxel* voxels,
                int* tileMask, int* numActiveTiles){
    k_flipTCells<<<opt.blocks, opt.threadsPerBlock>>>(dims, tile, voxels, tileMask, numActiveTiles);
    lastError("flipTCells");
}

void spawnTCells(int numXTravasing,
    Options opt, Comms comms, Dimensions dims, Dimensions tile, Dimensions tileDims,
    Voxel* voxels, Globals* globals, curandState_t* deviceRNG){
    k_spawnTCells<<<opt.blocks, opt.threadsPerBlock>>>(numXTravasing,
                                                        opt, comms, dims, tile, tileDims,
                                                        voxels, globals, deviceRNG);
}

void generateTCells(Options opt, Globals* globals){
    k_generateTCells<<<1, 1>>>(opt, globals);
    lastError("generateTCells");
}

void packSlice(Options& opt,
                        Voxel* data, Voxel* packedBuffer,
                        Dimensions dims, Dimensions tile, Dimensions tileDims, Comms comms,
                        int sliceDim, int sliceIdx){
    k_packSlice<<<opt.blocks, opt.threadsPerBlock>>>(opt, data, packedBuffer, dims, tile, tileDims, comms, sliceDim, sliceIdx);
    cudaDeviceSynchronize();
}

void unpackSlice(Options& opt,
                        Voxel* data, Voxel* packedBuffer,
                        Dimensions dims, Dimensions tile, Dimensions tileDims, Comms comms,
                        int sliceDim, int sliceIdx){
    k_unpackSlice<<<opt.blocks, opt.threadsPerBlock>>>(opt, data, packedBuffer, dims, tile, tileDims, comms, sliceDim, sliceIdx);
    cudaDeviceSynchronize();
}

//host rng functions
float drawValue(curandGenerator_t& gen){
    float r;
    curandGenerateUniform(gen, &r, 1);
    return r;
}

uint drawUInt(curandGenerator_t& gen){
    unsigned int r;
    curandGenerate(gen, &r, 1);
    return r;
}

int getBetween(curandGenerator_t& gen, int min, int max){
    if(min == max){
        return min;
    }
    uint r = drawUInt(gen);
    return (int)((r)%((uint)max - (uint)min) + (uint)min);
}

unsigned int drawPoisson(curandGenerator_t& gen, float lambda){
    unsigned int r;
    curandGeneratePoisson(gen, &r, 1, lambda);
    return r;
}

bool trialSuccess(curandGenerator_t& gen, float p){
    if(p > 1) return true;
    if(p < 0) return false;
    float roll = drawValue(gen);
    if(roll < p) return true;
    return false;
}

void setupHostRNG(curandGenerator_t& gen, Options& opt){
    curandCreateGeneratorHost(&gen, CURAND_RNG_PSEUDO_DEFAULT);
    curandSetPseudoRandomGeneratorSeed(gen, opt.seed);
}

__global__ void k_setCommTestData(Options opt, testData* data,
                                Dimensions dims, Dimensions tile, Dimensions tileDims, Comms comms){
    int start = blockIdx.x*blockDim.x + threadIdx.x;
    int maximum = dims.numPoints;
    int step = blockDim.x*gridDim.x;

    for(int i = start; i < maximum; i += step){
        data[i].rank = comms.rank;
        d_to3D_tiled(i, dims, tile, tileDims, data[i].x, data[i].y, data[i].z);
    }
}

__global__ void k_packSliceTestData(Options& opt,
                        testData* data, testData* packedBuffer,
                        Dimensions dims, Dimensions tile, Dimensions tileDims, Comms comms,
                        int sliceDim, int sliceIdx){

    int start = blockIdx.x*blockDim.x + threadIdx.x;
    int maximum = 0;
    int step = blockDim.x*gridDim.x;

    //X
    if(sliceDim == 0){
        maximum = comms.xSliceSize;
        for(int i = start; i < maximum; i += step){
            int x, y, z;
            x = sliceIdx;
            y = i / dims.z;
            z = i % dims.z;
            int idx = d_to1D_tiled(x, y, z, dims, tile, tileDims);
            packedBuffer[i].rank = data[idx].rank;
        }
    }

    //Y
    if(sliceDim == 1){
        maximum = comms.ySliceSize;
        for(int i = start; i < maximum; i += step){
            int x, y, z;
            x = i % dims.x;
            y = sliceIdx;
            z = i / dims.x;
            int idx = d_to1D_tiled(x, y, z, dims, tile, tileDims);
            packedBuffer[i].rank = data[idx].rank;
        }
    }

    //Z
    if(sliceDim == 2){
        maximum = comms.zSliceSize;
        for(int i = start; i < maximum; i += step){
            int x, y, z;
            x = i / dims.y;
            y = i % dims.y;
            z = sliceIdx;
            int idx = d_to1D_tiled(x, y, z, dims, tile, tileDims);
            packedBuffer[i].rank = data[idx].rank;
        }
    }
    
}

__global__ void k_unpackSliceTestData(Options& opt,
                        testData* data, testData* packedBuffer,
                        Dimensions dims, Dimensions tile, Dimensions tileDims, Comms comms,
                        int sliceDim, int sliceIdx){

    int start = blockIdx.x*blockDim.x + threadIdx.x;
    int maximum = 0;
    int step = blockDim.x*gridDim.x;

    //X
    if(sliceDim == 0){
        maximum = comms.xSliceSize;
        for(int i = start; i < maximum; i += step){
            int x, y, z;
            x = sliceIdx;
            y = i / dims.z;
            z = i % dims.z;
            int idx = d_to1D_tiled(x, y, z, dims, tile, tileDims);
            data[idx].rank = packedBuffer[i].rank;
        }
    }

    //Y
    if(sliceDim == 1){
        maximum = comms.ySliceSize;
        for(int i = start; i < maximum; i += step){
            int x, y, z;
            x = i % dims.x;
            y = sliceIdx;
            z = i / dims.x;
            int idx = d_to1D_tiled(x, y, z, dims, tile, tileDims);
            data[idx].rank = packedBuffer[i].rank;
        }
    }

    //Z
    if(sliceDim == 2){
        maximum = comms.zSliceSize;
        for(int i = start; i < maximum; i += step){
            int x, y, z;
            x = i / dims.y;
            y = i % dims.y;
            z = sliceIdx;
            int idx = d_to1D_tiled(x, y, z, dims, tile, tileDims);
            data[idx].rank = packedBuffer[i].rank;
        }
    }
    
}

//Communication testing
void setCommTestData(Options& opt, testData* data, Dimensions dims, Dimensions tile, Dimensions tileDims, Comms comms){
    k_setCommTestData<<<opt.blocks, opt.threadsPerBlock>>>(opt, data, dims, tile, tileDims, comms);
}

void packSliceTestData(Options& opt,
                        testData* data, testData* packedBuffer,
                        Dimensions dims, Dimensions tile, Dimensions tileDims, Comms comms,
                        int sliceDim, int sliceIdx){
    k_packSliceTestData<<<opt.blocks, opt.threadsPerBlock>>>(opt, data, packedBuffer, dims, tile, tileDims, comms, sliceDim, sliceIdx);
    cudaDeviceSynchronize();
}

void unpackSliceTestData(Options& opt,
                        testData* data, testData* packedBuffer,
                        Dimensions dims, Dimensions tile, Dimensions tileDims, Comms comms,
                        int sliceDim, int sliceIdx){
    k_unpackSliceTestData<<<opt.blocks, opt.threadsPerBlock>>>(opt, data, packedBuffer, dims, tile, tileDims, comms, sliceDim, sliceIdx);
    cudaDeviceSynchronize();
}