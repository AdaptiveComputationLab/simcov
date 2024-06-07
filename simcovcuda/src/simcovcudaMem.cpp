#include "simcovcudaMem.hpp"

template <class T>
upcxx::global_ptr<T, devMemKind> allocateDevBuffer(devAllocator &gpuAlloc, size_t count, std::string name){
    upcxx::global_ptr<T, devMemKind> buff = gpuAlloc.allocate<T>(count);
    if(!buff){
        printf("[GPU ERROR!] Rank %d failed to allocate buffer %s of count %ld\n", upcxx::rank_me(), name.c_str(), count);
        exit(1);
    }
    return buff;
}

void initComms(const Options& opt, Comms& comms,
                const Dimensions& dims, const Dimensions& simDims, const Dimensions& tile, const Dimensions& tileDims){
    comms.rank = upcxx::rank_me();
    comms.totalRanks = upcxx::rank_n();
    comms.rankX = comms.rank/(opt.rankDimZ*opt.rankDimY);
    comms.rankY = (comms.rank % (opt.rankDimZ*opt.rankDimY))/opt.rankDimZ;
    comms.rankZ = (comms.rank % (opt.rankDimZ*opt.rankDimY))%opt.rankDimZ;
    comms.rankDimX = opt.dimX/opt.rankDimX;
    comms.rankDimY = opt.dimY/opt.rankDimY;
    comms.rankDimZ = opt.dimZ/opt.rankDimZ;
    //TODO: Handle non even divides into sim space

    //compute predecessors and successors:
    int numXNeighbors = 0, numYNeighbors = 0, numZNeighbors = 0;
    //pred x
    if(comms.rankX == 0){
        comms.xPredecessor = -1;
    } else {
        int px = comms.rankX - 1;
        int py = comms.rankY;
        int pz = comms.rankZ;
        comms.xPredecessor = pz + py*opt.rankDimZ + px*(opt.rankDimY * opt.rankDimZ);
        numXNeighbors++;
    }

    //succ x
    if(comms.rankX == opt.rankDimX - 1){
        comms.xSuccessor = -1;
    } else {
        int px = comms.rankX + 1;
        int py = comms.rankY;
        int pz = comms.rankZ;
        comms.xSuccessor = pz + py*opt.rankDimZ + px*(opt.rankDimY * opt.rankDimZ);
        numXNeighbors++;
    }

    //pred y
    if(comms.rankY == 0){
        comms.yPredecessor = -1;
    } else {
        int px = comms.rankX;
        int py = comms.rankY - 1;
        int pz = comms.rankZ;
        comms.yPredecessor = pz + py*opt.rankDimZ + px*(opt.rankDimY * opt.rankDimZ);
        numYNeighbors++;
    }

    //succ y
    if(comms.rankY == opt.rankDimY - 1){
        comms.ySuccessor = -1;
    } else {
        int px = comms.rankX;
        int py = comms.rankY + 1;
        int pz = comms.rankZ;
        comms.ySuccessor = pz + py*opt.rankDimZ + px*(opt.rankDimY * opt.rankDimZ);
        numYNeighbors++;
    }

    //pred z
    if(comms.rankZ == 0){
        comms.zPredecessor = -1;
    } else {
        int px = comms.rankX;
        int py = comms.rankY;
        int pz = comms.rankZ - 1;
        comms.zPredecessor = pz + py*opt.rankDimZ + px*(opt.rankDimY * opt.rankDimZ);
        numZNeighbors++;
    }

    //succ z
    if(comms.rankZ == opt.rankDimZ - 1){
        comms.zSuccessor = -1;
    } else {
        int px = comms.rankX;
        int py = comms.rankY;
        int pz = comms.rankZ + 1;
        comms.zSuccessor = pz + py*opt.rankDimZ + px*(opt.rankDimY * opt.rankDimZ);
        numZNeighbors++;
    }

    //compute remainders for uneven divides
    int xrem = 0, yrem = 0, zrem = 0;
    if(comms.xSuccessor == -1) xrem = simDims.x % opt.rankDimX;
    if(comms.ySuccessor == -1) yrem = simDims.y % opt.rankDimY;
    if(comms.zSuccessor == -1) zrem = simDims.z % opt.rankDimZ;

    //Compute slice sizes in each dimension
    comms.xSliceSize = (simDims.z/opt.rankDimZ + numZNeighbors + zrem)*(simDims.y/opt.rankDimY + numYNeighbors + yrem);
    comms.ySliceSize = (simDims.x/opt.rankDimX  + numXNeighbors + xrem)*(simDims.z/opt.rankDimZ + numZNeighbors + zrem);
    comms.zSliceSize = (simDims.x/opt.rankDimX + numXNeighbors + xrem)*(simDims.y/opt.rankDimY + numYNeighbors + yrem);

    #ifdef COMMDEBUG
    printf("#[COMMDEBUG]: rank %d has slice sizes: (%d, %d, %d)\n", comms.rank, comms.xSliceSize, comms.ySliceSize, comms.zSliceSize);
    #endif
    
}


void setupGPUmemory(GPUMemory& gpumem, Options opt, Dimensions dims, Dimensions tileDims, Comms comms) {
    printf("### Allocating simcov device memory\n");

    //host pointers for extracting GPU memory
    gpumem.h_voxels = upcxx::new_array<Voxel>(dims.numPoints);
    gpumem.h_globals = upcxx::new_array<Globals>(1);

    //Allocate test memory
    const size_t mb = 1024*1024;
    gpumem.bufferSize = sizeof(Voxel)*dims.numPoints;
    gpumem.bufferSize += sizeof(curandState_t)*opt.blocks*opt.threadsPerBlock;
    gpumem.bufferSize += sizeof(int)*tileDims.numPoints*2;
    gpumem.bufferSize += sizeof(int);
    gpumem.bufferSize += sizeof(Globals)*(1 + opt.blocks);

    //We need 2 send and receive buffers along each dimension (12), reduced number if we
    //are missing a neighbor rank
    if(comms.xPredecessor != -1) gpumem.bufferSize += 2*sizeof(Voxel)*comms.xSliceSize;
    if(comms.xSuccessor != -1) gpumem.bufferSize += 2*sizeof(Voxel)*comms.xSliceSize;
    if(comms.yPredecessor != -1) gpumem.bufferSize += 2*sizeof(Voxel)*comms.ySliceSize;
    if(comms.ySuccessor != -1) gpumem.bufferSize += 2*sizeof(Voxel)*comms.ySliceSize;
    if(comms.zPredecessor != -1) gpumem.bufferSize += 2*sizeof(Voxel)*comms.zSliceSize;
    if(comms.zSuccessor != -1) gpumem.bufferSize += 2*sizeof(Voxel)*comms.zSliceSize;

    gpumem.bufferSize = gpumem.bufferSize + (mb - (gpumem.bufferSize%mb));

    int rank = upcxx::rank_me();
    int numDevices;
    cudaGetDeviceCount(&numDevices);
    int myDevice = upcxx::rank_me() % numDevices;
    gpumem.gpuDevice = upcxx::cuda_device(myDevice);
    // std::cout << "###" << upcxx::cuda_device::kind_info();
    cudaSetDevice(myDevice);
    printf("#[GPU Memory] Rank %d setting device to %d\n", rank, myDevice);
    gpumem.gpuAlloc = devAllocator(gpumem.gpuDevice, gpumem.bufferSize);

    gpumem.voxels = allocateDevBuffer<Voxel>(gpumem.gpuAlloc, dims.numPoints, "voxels");
    gpumem.deviceRNG = allocateDevBuffer<curandState_t>(gpumem.gpuAlloc, opt.threadsPerBlock*opt.blocks, "deviceRNG");
    gpumem.tileMask = allocateDevBuffer<int>(gpumem.gpuAlloc, tileDims.numPoints, "tileMask");
    gpumem.newTileMask = allocateDevBuffer<int>(gpumem.gpuAlloc, tileDims.numPoints, "newTileMask");
    gpumem.numActiveTiles = allocateDevBuffer<int>(gpumem.gpuAlloc, 1, "numActiveTiles");
    gpumem.globals = allocateDevBuffer<Globals>(gpumem.gpuAlloc, 1, "globals");
    gpumem.blockGlobals = allocateDevBuffer<Globals>(gpumem.gpuAlloc, opt.blocks, "blockGlobals");

    //Create a map that will hold all of the global pointers to the send/recv buffers
    //X
    if(comms.xPredecessor != -1) {
        gpumem.packedBuffers.insert( std::make_pair("xpredSend", allocateDevBuffer<Voxel>(gpumem.gpuAlloc, comms.xSliceSize, "xpredSend") ) );
        gpumem.distPackedBuffers.insert( std::make_pair( "xpredSend", upcxx::dist_object<upcxx::global_ptr<Voxel, devMemKind>>( gpumem.packedBuffers.at("xpredSend") ) ) );
        if(!gpumem.packedBuffers.at("xpredSend")){
            printf("#[GPU Memory] Rank %d failed to allocate xpredSend packedBuffer\n", rank);
            exit(1);
        }

        gpumem.packedBuffers.insert( std::make_pair("xpredRecv", allocateDevBuffer<Voxel>(gpumem.gpuAlloc, comms.xSliceSize, "xpredRecv") ) );
        gpumem.distPackedBuffers.insert( std::make_pair( "xpredRecv", upcxx::dist_object<upcxx::global_ptr<Voxel, devMemKind>>( gpumem.packedBuffers.at("xpredRecv") ) ) );
        if(!gpumem.packedBuffers.at("xpredRecv")){
            printf("#[GPU Memory] Rank %d failed to allocate xpredRecv packedBuffer\n", rank);
            exit(1);
        }


    } else {
        upcxx::global_ptr<Voxel, devMemKind> empty;
        gpumem.distPackedBuffers.insert( std::make_pair( "xpredSend", upcxx::dist_object<upcxx::global_ptr<Voxel, devMemKind>>( empty ) ) );
        gpumem.distPackedBuffers.insert( std::make_pair( "xpredRecv", upcxx::dist_object<upcxx::global_ptr<Voxel, devMemKind>>( empty ) ) );
    }
    if(comms.xSuccessor != -1) {
        gpumem.packedBuffers.insert( std::make_pair("xsuccSend", gpumem.gpuAlloc.allocate<Voxel>(comms.xSliceSize) ) );
        gpumem.distPackedBuffers.insert( std::make_pair( "xsuccSend", upcxx::dist_object<upcxx::global_ptr<Voxel, devMemKind>>( gpumem.packedBuffers.at("xsuccSend") ) ) );
        if(!gpumem.packedBuffers.at("xsuccSend")){
            printf("#[COMM TEST] Rank %d failed to allocate xsuccSend packedBuffer\n", rank);
            exit(1);
        }

        gpumem.packedBuffers.insert( std::make_pair("xsuccRecv", gpumem.gpuAlloc.allocate<Voxel>(comms.xSliceSize) ) );
        gpumem.distPackedBuffers.insert( std::make_pair( "xsuccRecv", upcxx::dist_object<upcxx::global_ptr<Voxel, devMemKind>>( gpumem.packedBuffers.at("xsuccRecv") ) ) );
        if(!gpumem.packedBuffers.at("xsuccRecv")){
            printf("#[GPU Memory] Rank %d failed to allocate xsuccRecv packedBuffer\n", rank);
            exit(1);
        }


    } else {
        upcxx::global_ptr<Voxel, devMemKind> empty;
        gpumem.distPackedBuffers.insert( std::make_pair( "xsuccSend", upcxx::dist_object<upcxx::global_ptr<Voxel, devMemKind>>( empty ) ) );
        gpumem.distPackedBuffers.insert( std::make_pair( "xsuccRecv", upcxx::dist_object<upcxx::global_ptr<Voxel, devMemKind>>( empty ) ) );
    }

    //Y
    if(comms.yPredecessor != -1) {
        gpumem.packedBuffers.insert( std::make_pair("ypredSend", gpumem.gpuAlloc.allocate<Voxel>(comms.ySliceSize) ) );
        gpumem.distPackedBuffers.insert( std::make_pair( "ypredSend", upcxx::dist_object<upcxx::global_ptr<Voxel, devMemKind>>( gpumem.packedBuffers.at("ypredSend") ) ) );
        if(!gpumem.packedBuffers.at("ypredSend")){
            printf("#[GPU Memory] Rank %d failed to allocate ypredSend packedBuffer\n", rank);
            exit(1);
        }

        gpumem.packedBuffers.insert( std::make_pair("ypredRecv", gpumem.gpuAlloc.allocate<Voxel>(comms.ySliceSize) ) );
        gpumem.distPackedBuffers.insert( std::make_pair( "ypredRecv", upcxx::dist_object<upcxx::global_ptr<Voxel, devMemKind>>( gpumem.packedBuffers.at("ypredRecv") ) ) );
        if(!gpumem.packedBuffers.at("ypredRecv")){
            printf("#[GPU Memory] Rank %d failed to allocate ypredRecv packedBuffer\n", rank);
            exit(1);
        }

    } else {
        upcxx::global_ptr<Voxel, devMemKind> empty;
        gpumem.distPackedBuffers.insert( std::make_pair( "ypredSend", upcxx::dist_object<upcxx::global_ptr<Voxel, devMemKind>>( empty ) ) );
        gpumem.distPackedBuffers.insert( std::make_pair( "ypredRecv", upcxx::dist_object<upcxx::global_ptr<Voxel, devMemKind>>( empty ) ) );
    }
    if(comms.ySuccessor != -1) {
        gpumem.packedBuffers.insert( std::make_pair("ysuccSend", gpumem.gpuAlloc.allocate<Voxel>(comms.ySliceSize) ) );
        gpumem.distPackedBuffers.insert( std::make_pair( "ysuccSend", upcxx::dist_object<upcxx::global_ptr<Voxel, devMemKind>>( gpumem.packedBuffers.at("ysuccSend") ) ) );
        if(!gpumem.packedBuffers.at("ysuccSend")){
            printf("#[GPU Memory] Rank %d failed to allocate ysuccSend packedBuffer\n", rank);
            exit(1);
        }


        gpumem.packedBuffers.insert( std::make_pair("ysuccRecv", gpumem.gpuAlloc.allocate<Voxel>(comms.ySliceSize) ) );
        gpumem.distPackedBuffers.insert( std::make_pair( "ysuccRecv", upcxx::dist_object<upcxx::global_ptr<Voxel, devMemKind>>( gpumem.packedBuffers.at("ysuccRecv") ) ) );
        if(!gpumem.packedBuffers.at("ysuccRecv")){
            printf("#[GPU Memory] Rank %d failed to allocate ysuccRecv packedBuffer\n", rank);
            exit(1);
        }


    } else {
        upcxx::global_ptr<Voxel, devMemKind> empty;
        gpumem.distPackedBuffers.insert( std::make_pair( "ysuccSend", upcxx::dist_object<upcxx::global_ptr<Voxel, devMemKind>>( empty ) ) );
        gpumem.distPackedBuffers.insert( std::make_pair( "ysuccRecv", upcxx::dist_object<upcxx::global_ptr<Voxel, devMemKind>>( empty ) ) );
    }
    
    //Z
    if(comms.zPredecessor != -1) {
        gpumem.packedBuffers.insert( std::make_pair("zpredSend", gpumem.gpuAlloc.allocate<Voxel>(comms.zSliceSize) ) );
        gpumem.distPackedBuffers.insert( std::make_pair( "zpredSend", upcxx::dist_object<upcxx::global_ptr<Voxel, devMemKind>>( gpumem.packedBuffers.at("zpredSend") ) ) );
        if(!gpumem.packedBuffers.at("zpredSend")){
            printf("#[GPU Memory] Rank %d failed to allocate zpredSend packedBuffer\n", rank);
            exit(1);
        }

        gpumem.packedBuffers.insert( std::make_pair("zpredRecv", gpumem.gpuAlloc.allocate<Voxel>(comms.zSliceSize) ) );
        gpumem.distPackedBuffers.insert( std::make_pair( "zpredRecv", upcxx::dist_object<upcxx::global_ptr<Voxel, devMemKind>>( gpumem.packedBuffers.at("zpredRecv") ) ) );
        if(!gpumem.packedBuffers.at("zpredRecv")){
            printf("#[GPU Memory] Rank %d failed to allocate zpredRecv packedBuffer\n", rank);
            exit(1);
        }

    } else {
        upcxx::global_ptr<Voxel, devMemKind> empty;
        gpumem.distPackedBuffers.insert( std::make_pair( "zpredSend", upcxx::dist_object<upcxx::global_ptr<Voxel, devMemKind>>( empty ) ) );
        gpumem.distPackedBuffers.insert( std::make_pair( "zpredRecv", upcxx::dist_object<upcxx::global_ptr<Voxel, devMemKind>>( empty ) ) );
    }
    if(comms.zSuccessor != -1) {
        gpumem.packedBuffers.insert( std::make_pair("zsuccSend", gpumem.gpuAlloc.allocate<Voxel>(comms.zSliceSize) ) );
        gpumem.distPackedBuffers.insert( std::make_pair( "zsuccSend", upcxx::dist_object<upcxx::global_ptr<Voxel, devMemKind>>( gpumem.packedBuffers.at("zsuccSend") ) ) );
        if(!gpumem.packedBuffers.at("zsuccSend")){
            printf("#[GPU Memory] Rank %d failed to allocate zsuccSend packedBuffer\n", rank);
            exit(1);
        }

        gpumem.packedBuffers.insert( std::make_pair("zsuccRecv", gpumem.gpuAlloc.allocate<Voxel>(comms.zSliceSize) ) );
        gpumem.distPackedBuffers.insert( std::make_pair( "zsuccRecv", upcxx::dist_object<upcxx::global_ptr<Voxel, devMemKind>>( gpumem.packedBuffers.at("zsuccRecv") ) ) );
        if(!gpumem.packedBuffers.at("zsuccRecv")){
            printf("#[GPU Memory] Rank %d failed to allocate zsuccRecv packedBuffer\n", rank);
            exit(1);
        }


    } else {
        upcxx::global_ptr<Voxel, devMemKind> empty;
        gpumem.distPackedBuffers.insert( std::make_pair( "zsuccSend", upcxx::dist_object<upcxx::global_ptr<Voxel, devMemKind>>( empty ) ) );
        gpumem.distPackedBuffers.insert( std::make_pair( "zsuccRecv", upcxx::dist_object<upcxx::global_ptr<Voxel, devMemKind>>( empty ) ) );
    }


    printf("#[GPU Memory] Rank %d Successfully Allocated memory of %zu!\n", rank, gpumem.bufferSize);

    upcxx::barrier();
    if(comms.xPredecessor != -1){
        gpumem.xpredSend = gpumem.distPackedBuffers.at("xsuccSend").fetch(comms.xPredecessor).wait();
        gpumem.xpredRecv = gpumem.distPackedBuffers.at("xsuccRecv").fetch(comms.xPredecessor).wait();
        std::cout << "# Rank " << comms.rank << " has xpredSend = " << gpumem.xpredSend << "\n";
        std::cout << "# Rank " << comms.rank << " has xpredRecv = " << gpumem.xpredRecv << "\n";
    }
    if(comms.xSuccessor != -1){
        gpumem.xsuccSend = gpumem.distPackedBuffers.at("xpredSend").fetch(comms.xSuccessor).wait();
        gpumem.xsuccRecv = gpumem.distPackedBuffers.at("xpredRecv").fetch(comms.xSuccessor).wait();
        std::cout << "# Rank " << comms.rank << " has xsuccSend = " << gpumem.xsuccSend << "\n";
        std::cout << "# Rank " << comms.rank << " has xsuccRecv = " << gpumem.xsuccRecv << "\n";
    }

    if(comms.yPredecessor != -1){
        gpumem.ypredSend = gpumem.distPackedBuffers.at("ysuccSend").fetch(comms.yPredecessor).wait();
        gpumem.ypredRecv = gpumem.distPackedBuffers.at("ysuccRecv").fetch(comms.yPredecessor).wait();
    }
    if(comms.ySuccessor != -1){
        gpumem.ysuccSend = gpumem.distPackedBuffers.at("ypredSend").fetch(comms.ySuccessor).wait();
        gpumem.ysuccRecv = gpumem.distPackedBuffers.at("ypredRecv").fetch(comms.ySuccessor).wait();
    }

    if(comms.zPredecessor != -1){
        gpumem.zpredSend = gpumem.distPackedBuffers.at("zsuccSend").fetch(comms.zPredecessor).wait();
        gpumem.zpredRecv = gpumem.distPackedBuffers.at("zsuccRecv").fetch(comms.zPredecessor).wait();
    }
    if(comms.zSuccessor != -1){
        gpumem.zsuccSend = gpumem.distPackedBuffers.at("zpredSend").fetch(comms.zSuccessor).wait();
        gpumem.zsuccRecv = gpumem.distPackedBuffers.at("zpredRecv").fetch(comms.zSuccessor).wait();
    }

}

void getGPUmemory(Dimensions dims, GPUMemory& gpumem){
    upcxx::copy(gpumem.voxels, gpumem.h_voxels, dims.numPoints).wait();
}

void setGPUmemory(Dimensions dims, GPUMemory& gpumem){
    upcxx::copy(gpumem.h_voxels, gpumem.voxels, dims.numPoints).wait();
}

void cleanGPUmemory(GPUMemory& gpumem){
    gpumem.gpuAlloc.deallocate(gpumem.voxels);

    gpumem.gpuAlloc.deallocate(gpumem.deviceRNG);
    gpumem.gpuAlloc.deallocate(gpumem.tileMask);
    gpumem.gpuAlloc.deallocate(gpumem.newTileMask);
    gpumem.gpuAlloc.deallocate(gpumem.numActiveTiles);
    gpumem.gpuAlloc.deallocate(gpumem.globals);
    gpumem.gpuAlloc.deallocate(gpumem.blockGlobals);

    upcxx::delete_array(gpumem.h_voxels);
    upcxx::delete_array(gpumem.h_globals);

    for(auto& it: gpumem.packedBuffers){
        gpumem.gpuAlloc.deallocate(it.second);
    }

    gpumem.gpuDevice.destroy();
}

void getGPUReduction(GPUMemory& gpumem){
    upcxx::copy(gpumem.globals, gpumem.h_globals, 1).wait();
}

void printActives(GPUMemory& gpumem,  Dimensions tileDims){
    upcxx::global_ptr<int> h_tileMask = upcxx::new_array<int>(tileDims.numPoints);
    upcxx::copy(gpumem.tileMask, h_tileMask, tileDims.numPoints).wait();

    for(int i = 0; i < tileDims.numPoints; i++){
        printf("%d,%d\n", i, h_tileMask.local()[i]);
    }

    upcxx::delete_array(h_tileMask);
}

void communicateVoxels(Options& opt, GPUMemory& gpumem, Comms& comms, Voxel* data, Dimensions& dims, Dimensions& tile, Dimensions& tileDims){
    //wave communication
    //x
    //Pack
    upcxx::future<> futureX = upcxx::make_future();
    if(comms.xPredecessor != -1){

        #ifdef COMMDEBUG
        printf("#[COMMDEBUG]: Rank %d Copying to x predecessor %d.\n", comms.rank, comms.xPredecessor);
        #endif

        packSlice(opt, data, gpumem.gpuAlloc.local(gpumem.packedBuffers.at("xpredSend")), dims, tile, tileDims, comms, 0, 1);
        upcxx::future<> f = upcxx::copy(gpumem.packedBuffers.at("xpredSend"), gpumem.xpredRecv, comms.xSliceSize);
        futureX = upcxx::when_all(futureX, f);
    }
    if(comms.xSuccessor != -1){

        #ifdef COMMDEBUG
        printf("#[COMMDEBUG]: Rank %d Copying to x successor %d.\n", comms.rank, comms.xSuccessor);
        #endif

        packSlice(opt, data, gpumem.gpuAlloc.local(gpumem.packedBuffers.at("xsuccSend")), dims, tile, tileDims, comms, 0, dims.x - 2);
        upcxx::future<> f = upcxx::copy(gpumem.packedBuffers.at("xsuccSend"), gpumem.xsuccRecv, comms.xSliceSize);
        futureX = upcxx::when_all(futureX, f);
    }
    futureX.wait();
    upcxx::barrier();

    //unpack x
    if(comms.xPredecessor != -1){
        unpackSlice(opt, data, gpumem.gpuAlloc.local(gpumem.packedBuffers.at("xpredRecv")), dims, tile, tileDims, comms, 0, 0);
    }
    if(comms.xSuccessor != -1){
        unpackSlice(opt, data, gpumem.gpuAlloc.local(gpumem.packedBuffers.at("xsuccRecv")), dims, tile, tileDims, comms, 0, dims.x - 1);
    }
    upcxx::barrier();

    //pack Y
    upcxx::future<> futureY = upcxx::make_future();
    if(comms.yPredecessor != -1){

        #ifdef COMMDEBUG
        printf("#[COMMDEBUG]: Rank %d Copying to y predecessor %d.\n", comms.rank, comms.yPredecessor);
        #endif

        packSlice(opt, data, gpumem.gpuAlloc.local(gpumem.packedBuffers.at("ypredSend")), dims, tile, tileDims, comms, 1, 1);
        upcxx::future<> f = upcxx::copy(gpumem.packedBuffers.at("ypredSend"), gpumem.ypredRecv, comms.ySliceSize);
        futureY = upcxx::when_all(futureY, f);
    }
    if(comms.ySuccessor != -1){

        #ifdef COMMDEBUG
        printf("#[COMMDEBUG]: Rank %d Copying to y successor %d.\n", comms.rank, comms.ySuccessor);
        #endif

        packSlice(opt, data, gpumem.gpuAlloc.local(gpumem.packedBuffers.at("ysuccSend")), dims, tile, tileDims, comms, 1, dims.y - 2);
        upcxx::future<> f = upcxx::copy(gpumem.packedBuffers.at("ysuccSend"), gpumem.ysuccRecv, comms.ySliceSize);
        futureY = upcxx::when_all(futureY, f);
    }
    futureY.wait();
    upcxx::barrier();
    
    //unpack y
    if(comms.yPredecessor != -1){
        unpackSlice(opt, data, gpumem.gpuAlloc.local(gpumem.packedBuffers.at("ypredRecv")), dims, tile, tileDims, comms, 1, 0);
    }
    if(comms.ySuccessor != -1){
        unpackSlice(opt, data, gpumem.gpuAlloc.local(gpumem.packedBuffers.at("ysuccRecv")), dims, tile, tileDims, comms, 1, dims.y - 1);
    }
    upcxx::barrier();

    //pack Z
    upcxx::future<> futureZ = upcxx::make_future();
    if(comms.zPredecessor != -1){

        #ifdef COMMDEBUG
        printf("#[COMMDEBUG]: Rank %d Copying to z predecessor %d.\n", comms.rank, comms.zPredecessor);
        #endif

        packSlice(opt, data, gpumem.gpuAlloc.local(gpumem.packedBuffers.at("zpredSend")), dims, tile, tileDims, comms, 2, 1);
        upcxx::future<> f = upcxx::copy(gpumem.packedBuffers.at("zpredSend"), gpumem.zpredRecv, comms.zSliceSize);
        futureZ = upcxx::when_all(futureZ, f);
    }
    if(comms.zSuccessor != -1){

        #ifdef COMMDEBUG
        printf("#[COMMDEBUG]: Rank %d Copying to z successor %d.\n", comms.rank, comms.zSuccessor);
        #endif

        packSlice(opt, data, gpumem.gpuAlloc.local(gpumem.packedBuffers.at("zsuccSend")), dims, tile, tileDims, comms, 2, dims.z - 2);
        upcxx::future<> f = upcxx::copy(gpumem.packedBuffers.at("zsuccSend"), gpumem.zsuccRecv, comms.zSliceSize);
        futureZ = upcxx::when_all(futureZ, f);
    }
    futureZ.wait();

    //communicate
    //some upcxx copy -> predRecvTestDataZ, etc
    upcxx::barrier();

    //unpack z
    if(comms.zPredecessor != -1){
        unpackSlice(opt, data, gpumem.gpuAlloc.local(gpumem.packedBuffers.at("zpredRecv")), dims, tile, tileDims, comms, 2, 0);
    }
    if(comms.zSuccessor != -1){
        unpackSlice(opt, data, gpumem.gpuAlloc.local(gpumem.packedBuffers.at("zsuccRecv")), dims, tile, tileDims, comms, 2, dims.z - 1);
    }
    upcxx::barrier();
}