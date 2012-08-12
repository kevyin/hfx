#include <cuda_runtime_api.h>
typedef enum cudaEvent_option_enum {
//  CUDA_EVENT_OPTION_DEFAULT       = cuduEventDefault,
    CUDA_EVENT_OPTION_BLOCKING_SYNC = cudaEventBlockingSync
} cudaEvent_option;


#if CUDART_VERSION < 3020
struct C2HS_COND_SENTRY_0;
#else
#endif
