#include <cuda_runtime_api.h>
typedef enum cudaMemHostAlloc_option_enum {
//  CUDA_MEMHOSTALLOC_OPTION_DEFAULT        = cudaHostAllocDefault,
    CUDA_MEMHOSTALLOC_OPTION_DEVICE_MAPPED  = cudaHostAllocMapped,
    CUDA_MEMHOSTALLOC_OPTION_PORTABLE       = cudaHostAllocPortable,
    CUDA_MEMHOSTALLOC_OPTION_WRITE_COMBINED = cudaHostAllocWriteCombined
} cudaMemHostAlloc_option;



