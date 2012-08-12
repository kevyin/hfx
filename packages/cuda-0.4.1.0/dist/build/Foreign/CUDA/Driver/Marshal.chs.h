#include <cuda.h>
#include "cbits/stubs.h"
typedef enum CUmemhostalloc_option_enum {
    CU_MEMHOSTALLOC_OPTION_PORTABLE       = CU_MEMHOSTALLOC_PORTABLE,
    CU_MEMHOSTALLOC_OPTION_DEVICE_MAPPED  = CU_MEMHOSTALLOC_DEVICEMAP,
    CU_MEMHOSTALLOC_OPTION_WRITE_COMBINED = CU_MEMHOSTALLOC_WRITECOMBINED
} CUmemhostalloc_option;



#if CUDA_VERSION < 4000
struct C2HS_COND_SENTRY_0;
#else
#endif
#if CUDA_VERSION < 4000
struct C2HS_COND_SENTRY_1;
#else
#endif
#if CUDA_VERSION < 4000
struct C2HS_COND_SENTRY_2;
#else
#endif
#if CUDA_VERSION < 4000
struct C2HS_COND_SENTRY_3;
#else
#endif
#if CUDA_VERSION < 3020
struct C2HS_COND_SENTRY_4;
#else
#endif
