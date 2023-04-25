#pragma once

#ifdef __cplusplus
extern "C" {
#endif

extern const char *SIMCOV_VERSION;
extern const char *SIMCOV_VERSION_DATE;
extern const char *SIMCOV_BUILD_DATE;
extern const char *SIMCOV_BRANCH;

#define MAX_BUILD_KMER_STR "Maximum kmer len=MAX_BUILD_KMER"

#ifdef USE_CUDA
#define WITH_CUDA "CUDA + CPU SW"
#else
#define WITH_CUDA "CPU SW Only"
#endif



#ifdef __cplusplus
}
#endif


