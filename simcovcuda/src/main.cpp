#include <upcxx/upcxx.hpp>
#include "simcovcuda_driver.hpp"
#include "simcovcudaMem.hpp"
#include "options.hpp"

#include <fstream>
#include <string>

#ifdef TIMING
#include <chrono>
#endif

// Base functions for computing 3D and 1D coordinate transformations
inline void to3D(int i, Dimensions dims,
                int&dx, int&dy, int& dz){
    dx = i/(dims.z*dims.y);
    dy = (i % (dims.z*dims.y))/dims.z;
    dz = (i % (dims.z*dims.y))%dims.z;
}
inline int to1D(int x, int y, int z, Dimensions dims){
    return z + y*dims.z + x*(dims.y * dims.z);
}
inline int to1D_tiled(int x, int y, int z,
                    Dimensions dims,
                    Dimensions tile,
                    Dimensions tileDims){
    //get the tile coords in tile space
    int tx, ty, tz;
    tx = x/tile.x;
    ty = y/tile.y;
    tz = z/tile.z;
    //get the tile id
    int ti = to1D(tx, ty, tz, tileDims);

    //get the offset within the tile
    int ox, oy, oz;
    ox = x - tx*tile.x;
    oy = y - ty*tile.y;
    oz = z - tz*tile.z;
    int oi = to1D(ox, oy, oz, tile);

    //result is the total number of voxels up until this tile + the offest
    return ti*tile.numPoints + oi;
}

void setupTiles(Options& opt,
                Dimensions& simDims,
                Dimensions& dims,
                Dimensions& tile,
                Dimensions& tileDims,
                Comms& comms){
    // dims.x = comms.myNumPoints;
    // TODO: Update for new communication paradigm (KL: 2023-09-08)
    simDims.x = opt.dimX;
    simDims.y = opt.dimY;
    simDims.z = opt.dimZ;
    simDims.numPoints = simDims.x*simDims.y*simDims.z;
    opt.numPoints = simDims.numPoints;
    tile.x = opt.tileX;
    tile.y = opt.tileY;
    tile.z = opt.tileZ;
    tile.numPoints = tile.x*tile.y*tile.z;

    //compute dims according to our ghost voxel configuration
    //TODO: Allow for ranks to nonevenly divide the simulation space
    dims.x = opt.dimX / opt.rankDimX;
    dims.y = opt.dimY / opt.rankDimY;
    dims.z = opt.dimZ / opt.rankDimZ;
    // add ghost voxels
    if(comms.xPredecessor != -1){
        dims.x++;
    }
    if(comms.xSuccessor != -1){
        dims.x++;
    } else {
        //add remainder to last rank in each dim
        dims.x += opt.dimX % opt.rankDimX;
    }

    if(comms.yPredecessor != -1){
        dims.y++;
    }
    if(comms.ySuccessor != -1){
        dims.y++;
    } else {
        dims.y += opt.dimY % opt.rankDimY;
    }

    if(comms.zPredecessor != -1){
        dims.z++;
    }
    if(comms.zSuccessor != -1){
        dims.z++;
    } else {
        //add remainder to last rank in each dim
        dims.z += opt.dimZ % opt.rankDimZ;
    }
    dims.numPoints = dims.x*dims.y*dims.z;

    #ifdef COMMDEBUG
    printf("#[COMM DEBUG]: Rank %d has dims (%d, %d, %d)\n", comms.rank, dims.x, dims.y, dims.z);
    #endif

    //set tile size to default if none provided
    if(tile.x == -1){
        for(int i = 3; i <= 20; i++){
            if(dims.x % i == 0){
                tile.x = i;
                break;
            }
        }
        if(tile.x == -1) tile.x = 1; //set to 1 if couldn't find alternative
    }
    if(tile.y == -1){
        for(int i = 3; i <= 20; i++){
            if(dims.y % i == 0){
                tile.y = i;
                break;
            }
        }
        if(tile.y == -1) tile.y = 1; //set to 1 if couldn't find alternative
    }
    if(tile.z == -1){
        for(int i = 3; i <= 20; i++){
            if(dims.z % i == 0){
                tile.z = i;
                break;
            }
        }
        if(tile.z == -1) tile.z = 1; //set to 1 if couldn't find alternative
    }
    tile.numPoints = tile.x*tile.y*tile.z;

    if(opt.dimZ == 1)
        opt.checkActivesRate = std::min(tile.x, tile.y);
    else
        opt.checkActivesRate = std::min(tile.x, std::min(tile.y, tile.z));    

    tileDims.x = dims.x/tile.x;
    tileDims.y = dims.y/tile.y;
    tileDims.z = dims.z/tile.z;
    tileDims.numPoints = tileDims.x*tileDims.y*tileDims.z;

    if(dims.x % tile.x != 0) {
        printf("### GPU ERROR! x Tile size does not evenly divide dimensions! %d->%d\n", tile.x, dims.x);
        exit(1);
    }
    if(dims.y % tile.y != 0) {
        printf("### GPU ERROR! y Tile size does not evenly divide dimensions! %d->%d\n", tile.y, dims.y);
        exit(1);
    }
    if(dims.z % tile.z != 0) {
        printf("### GPU ERROR! z Tile size does not evenly divide dimensions! %d->%d\n", tile.z, dims.z);
        exit(1);
    }
}

int updateCirculatingTCells(curandGenerator_t& hostRNG, Options opt,
                            Dimensions dims, upcxx::global_ptr<Globals, devMemKind> globals){
    upcxx::global_ptr<Globals> h_globals = upcxx::new_array<Globals>(1);
    upcxx::copy(globals, h_globals, 1).wait();

    float portionDying = (float)(h_globals.local()->numCirculatingTCells)/(float)opt.tcellVascularPeriod;
    int numDying = floor(portionDying);
    if(trialSuccess(hostRNG, (float)(portionDying - numDying) )) numDying++;
    h_globals.local()->numCirculatingTCells -= numDying;

    long long wholeLungVolume = (long long)opt.wholeLungX*
                                    (long long)opt.wholeLungY*
                                    (long long)opt.wholeLungZ;
    long double extravasateFraction = ((long double)dims.numPoints)/(long double)wholeLungVolume;
    long double portionXTravasing = extravasateFraction * h_globals.local()->numCirculatingTCells;
    int numXTravasing = floor(portionXTravasing);
    if(trialSuccess(hostRNG, (float)((float)portionXTravasing - (float)numXTravasing) )) numXTravasing++;
    upcxx::copy(h_globals, globals, 1).wait();
    upcxx:delete_array(h_globals);
    return numXTravasing;
}

void executeCommTest(Options& opt);
std::vector<std::array<int, 3>> getUniformInfections(int num, int dim_x, int dim_y, int dim_z){
    std::vector<std::array<int, 3>> infections;
    int x_splits = 1, y_splits = 1, z_splits = 1;
    while (x_splits * y_splits * z_splits < num) {
        double x_ratio = (double)dim_x / x_splits;
        double y_ratio = (double)dim_y / y_splits;
        double z_ratio = (double)dim_z / z_splits;
        double ratios[] = {x_ratio, y_ratio, z_ratio};
        std::sort(ratios, ratios + 3);
        if (ratios[2] == z_ratio) {
            z_splits++;
        } else if (ratios[2] == y_ratio) {
            y_splits++;
        } else {
            x_splits++;
        }
    }
    int x_spacing = (double)dim_x / (x_splits + 1);
    int y_spacing = (double)dim_y / (y_splits + 1);
    int z_spacing = (double)dim_z / (z_splits + 1);
    if (dim_z == 1) {
      for (int i = x_spacing; i < dim_x - 1; i += x_spacing) {
        for (int j = y_spacing; j < dim_y - 1; j += y_spacing) {
          if (infections.size() == num) return infections;
          infections.push_back({i, j, 0});
        }
      }
    } else {
      for (int i = x_spacing; i < dim_x - 1; i += x_spacing) {
        for (int j = y_spacing; j < dim_y - 1; j += y_spacing) {
          for (int k = z_spacing; k < dim_z - 1; k += z_spacing) {
            if (infections.size() == num) return infections;
            infections.push_back({i, j, k});
          }
        }
      }
    }
    return infections;
}

void globalSync(){
    cudaDeviceSynchronize();
    upcxx::barrier();
}

void localSync(){
    cudaDeviceSynchronize();
}

int main(int argc, char** argv){

    upcxx::init();

    #ifdef TIMING
    auto startTime = std::chrono::steady_clock::now();
    #endif

    Options opt = parseArgs(argc, argv);

    //COMM TEST
    if(opt.commTest){
        if(upcxx::rank_me() == 0)
            printf("#[COMM TEST]: Performing comm test...\n");

        executeCommTest(opt);

        upcxx::barrier();
        upcxx::finalize();

        exit(0);
    }

    Dimensions simDims, dims, tile, tileDims;
    Comms comms;
    simDims.x = opt.dimX;
    simDims.y = opt.dimY;
    simDims.z = opt.dimZ;
    simDims.numPoints = simDims.x*simDims.y*simDims.z;
    initComms(opt, comms, dims, simDims, tile, tileDims);

    curandGenerator_t hostRNG;
    setupHostRNG(hostRNG, opt);


    setupTiles(opt, simDims, dims, tile, tileDims, comms);
    printf("#Rank %d: tiles: {%d,%d,%d}\n", upcxx::rank_me(), tile.x, tile.y, tile.z);
    upcxx::barrier();
    GPUMemory gpumem;
    setupGPUmemory(gpumem, opt, dims, tileDims, comms);

    //downcast
    Voxel* voxels = gpumem.gpuAlloc.local(gpumem.voxels);
    curandState_t* deviceRNG = gpumem.gpuAlloc.local(gpumem.deviceRNG);

    Globals* globals = gpumem.gpuAlloc.local(gpumem.globals);
    Globals* blockGlobals = gpumem.gpuAlloc.local(gpumem.blockGlobals);
    int* tileMask = gpumem.gpuAlloc.local(gpumem.tileMask);
    int* newTileMask = gpumem.gpuAlloc.local(gpumem.newTileMask);
    int* numActiveTiles = gpumem.gpuAlloc.local(gpumem.numActiveTiles);

    //initialize globals
    Globals* l_globals = gpumem.h_globals.local();
    l_globals->numCirculatingTCells = 0;
    upcxx::copy(gpumem.h_globals, gpumem.globals, 1).wait();

    initVoxels(opt, comms, dims, tile, tileDims,
                deviceRNG,
                voxels);

    upcxx::barrier();

    std::vector<std::array<int, 3>> initialInfections = getUniformInfections(opt.numInfections, simDims.x, simDims.y, simDims.z);
    for(int i = 0; i < opt.numInfections; i++){
        initializeInfection(opt.initialInfections,
                            initialInfections[i][0], initialInfections[i][1], initialInfections[i][2],
                            dims, tile, tileDims, voxels, comms);
    }
    upcxx::barrier();

    //debug state logging
    if(upcxx::rank_me() == 0){
        printf("ts,virs,inf,heal,incb,expr,apop,dead,tcells,tvas\n");
    }

    //initial communication
    communicateVoxels(opt, gpumem, comms, voxels, dims, tile, tileDims);

    for(int i = 0; i < opt.timeSteps; i++){

        globalSync();

        if(i >= opt.tcellInitialDelay) {

            #ifdef PROFILE
            auto beginTime = std::chrono::steady_clock::now();
            #endif

            generateTCells(opt, globals);

            #ifdef PROFILE
            upcxx::barrier();
            auto finishTime = std::chrono::steady_clock::now();
            std::cout << "#[PROFILE]" << "generateTCells," << \
                    std::chrono::duration_cast<std::chrono::microseconds>(finishTime - beginTime).count() << \
                    std::endl;
            #endif

        }

        #ifdef PROFILE
        auto beginTime = std::chrono::steady_clock::now();
        #endif

        checkActiveTiles(opt, i,
                    dims, tile, tileDims,
                    voxels, tileMask, newTileMask, numActiveTiles);

        #ifdef PROFILE
        upcxx::barrier();
        auto finishTime = std::chrono::steady_clock::now();
        std::cout << "#[PROFILE]" << "checkActiveTiles," << \
                std::chrono::duration_cast<std::chrono::microseconds>(finishTime - beginTime).count() << \
                std::endl;
        #endif


        #ifdef PROFILE
        beginTime = std::chrono::steady_clock::now();
        #endif

        updateConcentrations(opt, dims, tile, voxels, tileMask, numActiveTiles);

        #ifdef PROFILE
        upcxx::barrier();
        finishTime = std::chrono::steady_clock::now();
        std::cout << "#[PROFILE]" << "updateConcentrations," << \
                std::chrono::duration_cast<std::chrono::microseconds>(finishTime - beginTime).count() << \
                std::endl;
        #endif

        #ifdef PROFILE
        beginTime = std::chrono::steady_clock::now();
        #endif

        updateEpithelialCells(opt, dims, tile,
                                voxels,
                                deviceRNG,
                                tileMask, numActiveTiles);

        #ifdef PROFILE
        upcxx::barrier();
        finishTime = std::chrono::steady_clock::now();
        std::cout << "#[PROFILE]" << "updateEpithelialCells," << \
                std::chrono::duration_cast<std::chrono::microseconds>(finishTime - beginTime).count() << \
                std::endl;
        #endif

        localSync();

        int numXTravasing = updateCirculatingTCells(hostRNG, opt, simDims, gpumem.globals);
        int rem = numXTravasing % upcxx::rank_n();
        numXTravasing = numXTravasing / upcxx::rank_n();
        if(upcxx::rank_me() == upcxx::rank_n() - 1) numXTravasing += rem;

        #ifdef PROFILE
        beginTime = std::chrono::steady_clock::now();
        #endif

        spawnTCells(numXTravasing, opt, comms, dims, tile, tileDims, voxels, globals, deviceRNG);
        
        #ifdef PROFILE
        upcxx::barrier();
        finishTime = std::chrono::steady_clock::now();
        std::cout << "#[PROFILE]" << "spawnTCells," << \
                std::chrono::duration_cast<std::chrono::microseconds>(finishTime - beginTime).count() << \
                std::endl;
        #endif
        
        #ifdef PROFILE
        beginTime = std::chrono::steady_clock::now();
        #endif

        ageTCells(opt, dims, tile,
                voxels,
                tileMask, numActiveTiles);

        #ifdef PROFILE
        upcxx::barrier();
        finishTime = std::chrono::steady_clock::now();
        std::cout << "#[PROFILE]" << "ageTCells," << \
                std::chrono::duration_cast<std::chrono::microseconds>(finishTime - beginTime).count() << \
                std::endl;
        #endif

        #ifdef PROFILE
        beginTime = std::chrono::steady_clock::now();
        #endif

        setupBind(opt, dims, tile, tileDims,
                voxels,
                deviceRNG,
                tileMask, numActiveTiles);
        
        #ifdef PROFILE
        upcxx::barrier();
        finishTime = std::chrono::steady_clock::now();
        std::cout << "#[PROFILE]" << "setupBind," << \
                std::chrono::duration_cast<std::chrono::microseconds>(finishTime - beginTime).count() << \
                std::endl;
        #endif

        #ifdef PROFILE
        beginTime = std::chrono::steady_clock::now();
        #endif

        globalSync();
        communicateVoxels(opt, gpumem, comms, voxels, dims, tile, tileDims);
        globalSync();

        #ifdef PROFILE
        upcxx::barrier();
        finishTime = std::chrono::steady_clock::now();
        std::cout << "#[PROFILE]" << "communicateVoxels0," << \
                std::chrono::duration_cast<std::chrono::microseconds>(finishTime - beginTime).count() << \
                std::endl;
        #endif

        #ifdef PROFILE
        beginTime = std::chrono::steady_clock::now();
        #endif
        executeBind(opt, dims, tile, tileDims,
                    voxels,
                    tileMask, numActiveTiles);
        #ifdef PROFILE
        upcxx::barrier();
        finishTime = std::chrono::steady_clock::now();
        std::cout << "#[PROFILE]" << "executeBind," << \
                std::chrono::duration_cast<std::chrono::microseconds>(finishTime - beginTime).count() << \
                std::endl;
        #endif
        
        #ifdef PROFILE
        beginTime = std::chrono::steady_clock::now();
        #endif
        setupMove(opt, dims, tile,
                    voxels, deviceRNG,
                    tileMask, numActiveTiles);
        #ifdef PROFILE
        upcxx::barrier();
        finishTime = std::chrono::steady_clock::now();
        std::cout << "#[PROFILE]" << "setupMove," << \
                std::chrono::duration_cast<std::chrono::microseconds>(finishTime - beginTime).count() << \
                std::endl;
        #endif


        #ifdef PROFILE
        beginTime = std::chrono::steady_clock::now();
        #endif

        globalSync();
        communicateVoxels(opt, gpumem, comms, voxels, dims, tile, tileDims);
        globalSync();

        #ifdef PROFILE
        upcxx::barrier();
        finishTime = std::chrono::steady_clock::now();
        std::cout << "#[PROFILE]" << "communicateVoxels1," << \
                std::chrono::duration_cast<std::chrono::microseconds>(finishTime - beginTime).count() << \
                std::endl;
        #endif


        #ifdef PROFILE
        beginTime = std::chrono::steady_clock::now();
        #endif
        declareWinners(opt, dims, tile,
                        voxels,
                        tileMask, numActiveTiles);
        #ifdef PROFILE
        upcxx::barrier();
        finishTime = std::chrono::steady_clock::now();
        std::cout << "#[PROFILE]" << "declareWinners," << \
                std::chrono::duration_cast<std::chrono::microseconds>(finishTime - beginTime).count() << \
                std::endl;
        #endif

        #ifdef PROFILE
        beginTime = std::chrono::steady_clock::now();
        #endif
        flipTCells(opt, dims, tile,
                    voxels,
                    tileMask, numActiveTiles);
        #ifdef PROFILE
        upcxx::barrier();
        finishTime = std::chrono::steady_clock::now();
        std::cout << "#[PROFILE]" << "flipTCells," << \
                std::chrono::duration_cast<std::chrono::microseconds>(finishTime - beginTime).count() << \
                std::endl;
        #endif

        if(opt.reducePeriod > 0)
        if(i % opt.reducePeriod == 0){

            #ifdef PROFILE
            beginTime = std::chrono::steady_clock::now();
            #endif

            reduceGPU(opt, dims,
                    voxels,
                    blockGlobals, globals);
            cudaDeviceSynchronize();
            upcxx::barrier();
            
            getGPUReduction(gpumem);

            cudaDeviceSynchronize();
            upcxx::barrier();
            Globals* l_globals = gpumem.h_globals.local();
            Globals stats;
            upcxx::barrier();
            stats.totalVirions = upcxx::reduce_all(l_globals->totalVirions, upcxx::op_fast_add).wait();
            stats.totalInflammation = upcxx::reduce_all(l_globals->totalInflammation, upcxx::op_fast_add).wait();
            stats.totalHealthy = upcxx::reduce_all(l_globals->totalHealthy, upcxx::op_fast_add).wait();
            stats.totalIncubating = upcxx::reduce_all(l_globals->totalIncubating, upcxx::op_fast_add).wait();
            stats.totalExpressing = upcxx::reduce_all(l_globals->totalExpressing, upcxx::op_fast_add).wait();
            stats.totalApoptotic = upcxx::reduce_all(l_globals->totalApoptotic, upcxx::op_fast_add).wait();
            stats.totalDead = upcxx::reduce_all(l_globals->totalDead, upcxx::op_fast_add).wait();
            stats.totalTCells = upcxx::reduce_all(l_globals->totalTCells, upcxx::op_fast_add).wait();
            stats.numCirculatingTCells = l_globals->numCirculatingTCells;
            if(upcxx::rank_me() == 0){
                printf("%d,%f,%f,%llu,%llu,%llu,%llu,%llu,%llu,%llu\n", i, stats.totalVirions, stats.totalInflammation,
                                                        stats.totalHealthy,
                                                        stats.totalIncubating,
                                                        stats.totalExpressing,
                                                        stats.totalApoptotic,
                                                        stats.totalDead,
                                                        stats.totalTCells,
                                                        stats.numCirculatingTCells);
            }
            #ifdef PROFILE
            upcxx::barrier();
            finishTime = std::chrono::steady_clock::now();
            std::cout << "#[PROFILE]" << "reduction," << \
                    std::chrono::duration_cast<std::chrono::microseconds>(finishTime - beginTime).count() << \
                    std::endl;
            #endif
        }
        
        if(opt.samplePeriod > 0)
        if(i % opt.samplePeriod == 0){

            std::ofstream myStateLog;
            std::string ts = std::to_string(i);
            ts.insert(ts.begin(), 9 - ts.size(), '0');
            std::string fname = "./output/stateFiles/state_" +
                                std::to_string(upcxx::rank_me()) + "_" +
                                ts + ".csv";
            myStateLog.open(fname);
            myStateLog << "ts,x,y,z,virions,inflammation,tcells,epitype\n";

            getGPUmemory(dims, gpumem);
            Voxel* hvoxels = gpumem.h_voxels.local();
            for(int ii = 0; ii < dims.numPoints; ii++){
                int tx, ty, tz;
                to3D(ii, dims, tx, ty, tz);
                int idx = to1D_tiled(tx, ty, tz, dims, tile, tileDims);
                myStateLog << i << "," <<
                        hvoxels[idx].x << "," <<
                        hvoxels[idx].y << "," <<
                        hvoxels[idx].z << "," <<
                        hvoxels[idx].virions << "," <<
                        hvoxels[idx].inflammation << "," <<
                        (int)hvoxels[idx].hasTCell << "," <<
                        (int)hvoxels[idx].cellType << "\n";
            }
            myStateLog.close();
        }
    }

    cleanGPUmemory(gpumem);

    #ifdef TIMING
    if(upcxx::rank_me() == 0){
        auto endTime = std::chrono::steady_clock::now();
        std::cout << "###!TOTALTIME:" << \
                std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count() << \
                " microseconds" << std::endl;
    }
    #endif

    upcxx::finalize();
}


void executeCommTest(Options& opt){

    int rank = upcxx::rank_me();
    int numDevices;
    int myDevice;
    cudaGetDeviceCount(&numDevices);
    myDevice = upcxx::rank_me() % numDevices;

    printf("#[COMM TEST] Hello from rank %d!\n", rank);

    //get number of devices
    printf("#[COMM TEST] Rank %d, Cuda detects %d devices on my node\n", rank, numDevices);

    cudaDeviceProp devProps;
    cudaGetDeviceProperties(&devProps, myDevice);

    std::cout << "#[COMM TEST] Rank " << rank  << " is using device: ";
    std::vector<std::tuple<int, int>> r = {{0,4}, {4,6}, {6,8}, {8,10}, {10,16}};
    for (auto t : r){
	std::cout << "-";
	for (int i = std::get<0>(t); i < std::get<1>(t); i++){
		std::cout << std::hex << std::setfill('0') << std::setw(2) << (unsigned)(unsigned char)devProps.uuid.bytes[i];
	}
    }
    std::cout << "\n";

    Dimensions dims, simDims, tile, tileDims;
    Comms comms;
    simDims.x = opt.dimX;
    simDims.y = opt.dimY;
    simDims.z = opt.dimZ;
    simDims.numPoints = simDims.x*simDims.y*simDims.z;
    initComms(opt, comms, dims, simDims, tile, tileDims);

    printf("#[COMM TEST] Rank %d, rank position (%d,%d,%d)\n",
            rank, comms.rankX, comms.rankY, comms.rankZ);

    printf("#[COMM TEST] Rank %d predecessors = (%d,%d,%d)\n", rank, comms.xPredecessor, comms.yPredecessor, comms.zPredecessor);
    printf("#[COMM TEST] Rank %d successors = (%d,%d,%d)\n", rank, comms.xSuccessor, comms.ySuccessor, comms.zSuccessor);
    printf("#[COMM TEST] Rank %d slice sizes = (%d,%d,%d)\n", rank, comms.xSliceSize, comms.ySliceSize, comms.zSliceSize);

    setupTiles(opt, simDims, dims, tile, tileDims, comms);

    printf("#[COMM TEST] Rank %d dims: (%d,%d,%d) = %d voxels\n", rank, dims.x, dims.y, dims.z, dims.numPoints);

    //Allocate test memory
    const size_t mb = 1024*1024;
    size_t bufferSize = sizeof(testData)*dims.numPoints;

    //We need 2 send and receive buffers along each dimension (12), reduced number if we
    //are missing a neighbor rank
    if(comms.xPredecessor != -1) bufferSize += 2*sizeof(testData)*comms.xSliceSize;
    if(comms.xSuccessor != -1) bufferSize += 2*sizeof(testData)*comms.xSliceSize;
    if(comms.yPredecessor != -1) bufferSize += 2*sizeof(testData)*comms.ySliceSize;
    if(comms.ySuccessor != -1) bufferSize += 2*sizeof(testData)*comms.ySliceSize;
    if(comms.zPredecessor != -1) bufferSize += 2*sizeof(testData)*comms.zSliceSize;
    if(comms.zSuccessor != -1) bufferSize += 2*sizeof(testData)*comms.zSliceSize;

    bufferSize = bufferSize + (mb - (bufferSize%mb));

    upcxx::cuda_device gpuDevice = upcxx::cuda_device(myDevice);
    cudaSetDevice(myDevice);
    printf("#[COMM TEST] Rank %d setting device to %d\n", rank, myDevice);
    upcxx::device_allocator<upcxx::cuda_device> gpuAlloc = upcxx::device_allocator<upcxx::cuda_device>(gpuDevice, bufferSize);
    upcxx::global_ptr<testData, devMemKind> testDataBuff = gpuAlloc.allocate<testData>(dims.numPoints);
    if(!testDataBuff){
        printf("#[COMM TEST] Rank %d failed to allocate testData memory!\n", rank);
        exit(1);
    }

    //Create a map that will hold all of the global pointers to the send/recv buffers
    std::unordered_map<std::string, upcxx::global_ptr<testData, devMemKind>> packedBuffers;
    std::unordered_map<std::string, distTestData> distPackedBuffers;
    if(comms.xPredecessor != -1) {
        packedBuffers.insert( std::make_pair("xpredSend", gpuAlloc.allocate<testData>(comms.xSliceSize) ) );
        distPackedBuffers.insert( std::make_pair( "xpredSend", upcxx::dist_object<upcxx::global_ptr<testData, devMemKind>>( packedBuffers.at("xpredSend") ) ) );
        if(!packedBuffers.at("xpredSend")){
            printf("#[COMM TEST] Rank %d failed to allocate xpredSend packedBuffer\n", rank);
            exit(1);
        }
        packedBuffers.insert( std::make_pair("xpredRecv", gpuAlloc.allocate<testData>(comms.xSliceSize) ) );
        distPackedBuffers.insert( std::make_pair( "xpredRecv", upcxx::dist_object<upcxx::global_ptr<testData, devMemKind>>( packedBuffers.at("xpredRecv") ) ) );
        if(!packedBuffers.at("xpredRecv")){
            printf("#[COMM TEST] Rank %d failed to allocate xpredRecv packedBuffer\n", rank);
            exit(1);
        }
    } else {
        upcxx::global_ptr<testData, devMemKind> empty;
        distPackedBuffers.insert( std::make_pair( "xpredSend", upcxx::dist_object<upcxx::global_ptr<testData, devMemKind>>( empty ) ) );
        distPackedBuffers.insert( std::make_pair( "xpredRecv", upcxx::dist_object<upcxx::global_ptr<testData, devMemKind>>( empty ) ) );
    }
    if(comms.xSuccessor != -1) {
        packedBuffers.insert( std::make_pair("xsuccSend", gpuAlloc.allocate<testData>(comms.xSliceSize) ) );
        distPackedBuffers.insert( std::make_pair( "xsuccSend", upcxx::dist_object<upcxx::global_ptr<testData, devMemKind>>( packedBuffers.at("xsuccSend") ) ) );
        if(!packedBuffers.at("xsuccSend")){
            printf("#[COMM TEST] Rank %d failed to allocate xsuccSend packedBuffer\n", rank);
            exit(1);
        }
        packedBuffers.insert( std::make_pair("xsuccRecv", gpuAlloc.allocate<testData>(comms.xSliceSize) ) );
        distPackedBuffers.insert( std::make_pair( "xsuccRecv", upcxx::dist_object<upcxx::global_ptr<testData, devMemKind>>( packedBuffers.at("xsuccRecv") ) ) );
        if(!packedBuffers.at("xsuccRecv")){
            printf("#[COMM TEST] Rank %d failed to allocate xsuccRecv packedBuffer\n", rank);
            exit(1);
        }
    } else {
        upcxx::global_ptr<testData, devMemKind> empty;
        distPackedBuffers.insert( std::make_pair( "xsuccSend", upcxx::dist_object<upcxx::global_ptr<testData, devMemKind>>( empty ) ) );
        distPackedBuffers.insert( std::make_pair( "xsuccRecv", upcxx::dist_object<upcxx::global_ptr<testData, devMemKind>>( empty ) ) );
    }

    if(comms.yPredecessor != -1) {
        packedBuffers.insert( std::make_pair("ypredSend", gpuAlloc.allocate<testData>(comms.ySliceSize) ) );
        distPackedBuffers.insert( std::make_pair( "ypredSend", upcxx::dist_object<upcxx::global_ptr<testData, devMemKind>>( packedBuffers.at("ypredSend") ) ) );
        if(!packedBuffers.at("ypredSend")){
            printf("#[COMM TEST] Rank %d failed to allocate ypredSend packedBuffer\n", rank);
            exit(1);
        }
        packedBuffers.insert( std::make_pair("ypredRecv", gpuAlloc.allocate<testData>(comms.ySliceSize) ) );
        distPackedBuffers.insert( std::make_pair( "ypredRecv", upcxx::dist_object<upcxx::global_ptr<testData, devMemKind>>( packedBuffers.at("ypredRecv") ) ) );
        if(!packedBuffers.at("ypredRecv")){
            printf("#[COMM TEST] Rank %d failed to allocate ypredRecv packedBuffer\n", rank);
            exit(1);
        }
    } else {
        upcxx::global_ptr<testData, devMemKind> empty;
        distPackedBuffers.insert( std::make_pair( "ypredSend", upcxx::dist_object<upcxx::global_ptr<testData, devMemKind>>( empty ) ) );
        distPackedBuffers.insert( std::make_pair( "ypredRecv", upcxx::dist_object<upcxx::global_ptr<testData, devMemKind>>( empty ) ) );
    }

    if(comms.ySuccessor != -1) {
        packedBuffers.insert( std::make_pair("ysuccSend", gpuAlloc.allocate<testData>(comms.ySliceSize) ) );
        distPackedBuffers.insert( std::make_pair( "ysuccSend", upcxx::dist_object<upcxx::global_ptr<testData, devMemKind>>( packedBuffers.at("ysuccSend") ) ) );
        if(!packedBuffers.at("ysuccSend")){
            printf("#[COMM TEST] Rank %d failed to allocate ysuccSend packedBuffer\n", rank);
            exit(1);
        }
        packedBuffers.insert( std::make_pair("ysuccRecv", gpuAlloc.allocate<testData>(comms.ySliceSize) ) );
        distPackedBuffers.insert( std::make_pair( "ysuccRecv", upcxx::dist_object<upcxx::global_ptr<testData, devMemKind>>( packedBuffers.at("ysuccRecv") ) ) );
        if(!packedBuffers.at("ysuccRecv")){
            printf("#[COMM TEST] Rank %d failed to allocate ysuccRecv packedBuffer\n", rank);
            exit(1);
        }
    } else {
        upcxx::global_ptr<testData, devMemKind> empty;
        distPackedBuffers.insert( std::make_pair( "ysuccSend", upcxx::dist_object<upcxx::global_ptr<testData, devMemKind>>( empty ) ) );
        distPackedBuffers.insert( std::make_pair( "ysuccRecv", upcxx::dist_object<upcxx::global_ptr<testData, devMemKind>>( empty ) ) );
    }


    if(comms.zPredecessor != -1) {
        packedBuffers.insert( std::make_pair("zpredSend", gpuAlloc.allocate<testData>(comms.zSliceSize) ) );
        distPackedBuffers.insert( std::make_pair( "zpredSend", upcxx::dist_object<upcxx::global_ptr<testData, devMemKind>>( packedBuffers.at("zpredSend") ) ) );
        if(!packedBuffers.at("zpredSend")){
            printf("#[COMM TEST] Rank %d failed to allocate zpredSend packedBuffer\n", rank);
            exit(1);
        }
        packedBuffers.insert( std::make_pair("zpredRecv", gpuAlloc.allocate<testData>(comms.zSliceSize) ) );
        distPackedBuffers.insert( std::make_pair( "zpredRecv", upcxx::dist_object<upcxx::global_ptr<testData, devMemKind>>( packedBuffers.at("zpredRecv") ) ) );
        if(!packedBuffers.at("zpredRecv")){
            printf("#[COMM TEST] Rank %d failed to allocate zpredRecv packedBuffer\n", rank);
            exit(1);
        }
    } else {
        upcxx::global_ptr<testData, devMemKind> empty;
        distPackedBuffers.insert( std::make_pair( "zpredSend", upcxx::dist_object<upcxx::global_ptr<testData, devMemKind>>( empty ) ) );
        distPackedBuffers.insert( std::make_pair( "zpredRecv", upcxx::dist_object<upcxx::global_ptr<testData, devMemKind>>( empty ) ) );
    }

    if(comms.zSuccessor != -1) {
        packedBuffers.insert( std::make_pair("zsuccSend", gpuAlloc.allocate<testData>(comms.zSliceSize) ) );
        distPackedBuffers.insert( std::make_pair( "zsuccSend", upcxx::dist_object<upcxx::global_ptr<testData, devMemKind>>( packedBuffers.at("zsuccSend") ) ) );
        if(!packedBuffers.at("zsuccSend")){
            printf("#[COMM TEST] Rank %d failed to allocate zsuccSend packedBuffer\n", rank);
            exit(1);
        }
        packedBuffers.insert( std::make_pair("zsuccRecv", gpuAlloc.allocate<testData>(comms.zSliceSize) ) );
        distPackedBuffers.insert( std::make_pair( "zsuccRecv", upcxx::dist_object<upcxx::global_ptr<testData, devMemKind>>( packedBuffers.at("zsuccRecv") ) ) );
        if(!packedBuffers.at("zsuccRecv")){
            printf("#[COMM TEST] Rank %d failed to allocate zsuccRecv packedBuffer\n", rank);
            exit(1);
        }
    } else {
        upcxx::global_ptr<testData, devMemKind> empty;
        distPackedBuffers.insert( std::make_pair( "zsuccSend", upcxx::dist_object<upcxx::global_ptr<testData, devMemKind>>( empty ) ) );
        distPackedBuffers.insert( std::make_pair( "zsuccRecv", upcxx::dist_object<upcxx::global_ptr<testData, devMemKind>>( empty ) ) );
    }

    printf("#[COMM TEST] Rank %d Successfully Allocated memory of %zu!\n", rank, bufferSize);

    testData* data = gpuAlloc.local(testDataBuff);

    //Set test data
    setCommTestData(opt, data, dims, tile, tileDims, comms);
    cudaDeviceSynchronize();
    upcxx::barrier();

    //wave communication
    //x
    //Pack
    if(comms.xPredecessor != -1){
        packSliceTestData(opt, data, gpuAlloc.local(packedBuffers.at("xpredSend")), dims, tile, tileDims, comms, 0, 1);
        upcxx::global_ptr<testData, devMemKind> predData = distPackedBuffers.at("xsuccRecv").fetch(comms.xPredecessor).wait();
        upcxx::copy(packedBuffers.at("xpredSend"), predData, comms.xSliceSize).wait();
    }
    if(comms.xSuccessor != -1){
        packSliceTestData(opt, data, gpuAlloc.local(packedBuffers.at("xsuccSend")), dims, tile, tileDims, comms, 0, dims.x - 2);
        upcxx::global_ptr<testData, devMemKind> succData = distPackedBuffers.at("xpredRecv").fetch(comms.xSuccessor).wait();
        upcxx::copy(packedBuffers.at("xsuccSend"), succData, comms.xSliceSize).wait();
    }
    upcxx::barrier();

    //unpack x
    if(comms.xPredecessor != -1){
        unpackSliceTestData(opt, data, gpuAlloc.local(packedBuffers.at("xpredRecv")), dims, tile, tileDims, comms, 0, 0);
    }
    if(comms.xSuccessor != -1){
        unpackSliceTestData(opt, data, gpuAlloc.local(packedBuffers.at("xsuccRecv")), dims, tile, tileDims, comms, 0, dims.x - 1);
    }
    upcxx::barrier();

    //pack Y
    if(comms.yPredecessor != -1){
        packSliceTestData(opt, data, gpuAlloc.local(packedBuffers.at("ypredSend")), dims, tile, tileDims, comms, 1, 1);
        upcxx::global_ptr<testData, devMemKind> predData = distPackedBuffers.at("ysuccRecv").fetch(comms.yPredecessor).wait();
        upcxx::copy(packedBuffers.at("ypredSend"), predData, comms.ySliceSize).wait();
    }
    if(comms.ySuccessor != -1){
        packSliceTestData(opt, data, gpuAlloc.local(packedBuffers.at("ysuccSend")), dims, tile, tileDims, comms, 1, dims.y - 2);
        upcxx::global_ptr<testData, devMemKind> succData = distPackedBuffers.at("ypredRecv").fetch(comms.ySuccessor).wait();
        upcxx::copy(packedBuffers.at("ysuccSend"), succData, comms.ySliceSize).wait();
    }
    upcxx::barrier();
    
    //unpack y
    if(comms.yPredecessor != -1){
        unpackSliceTestData(opt, data, gpuAlloc.local(packedBuffers.at("ypredRecv")), dims, tile, tileDims, comms, 1, 0);
    }
    if(comms.ySuccessor != -1){
        unpackSliceTestData(opt, data, gpuAlloc.local(packedBuffers.at("ysuccRecv")), dims, tile, tileDims, comms, 1, dims.y - 1);
    }
    upcxx::barrier();

    //pack Z
    if(comms.zPredecessor != -1){
        packSliceTestData(opt, data, gpuAlloc.local(packedBuffers.at("zpredSend")), dims, tile, tileDims, comms, 2, 1);
        upcxx::global_ptr<testData, devMemKind> predData = distPackedBuffers.at("zsuccRecv").fetch(comms.zPredecessor).wait();
        upcxx::copy(packedBuffers.at("zpredSend"), predData, comms.zSliceSize).wait();
    }
    if(comms.zSuccessor != -1){
        packSliceTestData(opt, data, gpuAlloc.local(packedBuffers.at("zsuccSend")), dims, tile, tileDims, comms, 2, dims.z - 2);
        upcxx::global_ptr<testData, devMemKind> succData = distPackedBuffers.at("zpredRecv").fetch(comms.zSuccessor).wait();
        upcxx::copy(packedBuffers.at("zsuccSend"), succData, comms.zSliceSize).wait();
    }

    //communicate
    //some upcxx copy -> predRecvTestDataZ, etc
    upcxx::barrier();

    //unpack z
    if(comms.zPredecessor != -1){
        unpackSliceTestData(opt, data, gpuAlloc.local(packedBuffers.at("zpredRecv")), dims, tile, tileDims, comms, 2, 0);
    }
    if(comms.zSuccessor != -1){
        unpackSliceTestData(opt, data, gpuAlloc.local(packedBuffers.at("zsuccRecv")), dims, tile, tileDims, comms, 2, dims.z - 1);
    }
    upcxx::barrier();

    //barrier at end of communication
    upcxx::barrier();

    //verify data
    std::string fname = "./notebooks/domain/" +
                                std::to_string(upcxx::rank_me()) + ".csv";
    std::ofstream commTestFile;
    commTestFile.open(fname);

    upcxx::global_ptr<testData> h_data = upcxx::new_array<testData>(dims.numPoints);
    upcxx::copy(testDataBuff, h_data, dims.numPoints).wait();
    testData* l_h_data = h_data.local();
    commTestFile << "x,y,z,rank\n";
    for(int i = 0; i < dims.numPoints; i++){
        commTestFile << l_h_data[i].x << "," << l_h_data[i].y << "," << l_h_data[i].z << "," << l_h_data[i].rank << "\n";
    }
    upcxx::delete_array(h_data);

    commTestFile.close();

    //report communication info


    for(auto& it: packedBuffers){
        gpuAlloc.deallocate(it.second);
    }
    gpuAlloc.deallocate(testDataBuff);
    gpuAlloc.destroy();

}
