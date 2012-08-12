#include <cuda_runtime_api.h>
typedef struct cudaDeviceProp   cudaDeviceProp;

typedef enum
{
    cudaDeviceFlagScheduleAuto    = cudaDeviceScheduleAuto,
    cudaDeviceFlagScheduleSpin    = cudaDeviceScheduleSpin,
    cudaDeviceFlagScheduleYield   = cudaDeviceScheduleYield,
    cudaDeviceFlagBlockingSync    = cudaDeviceBlockingSync,
    cudaDeviceFlagMapHost         = cudaDeviceMapHost,
#if CUDART_VERSION >= 3000
    cudaDeviceFlagLMemResizeToMax = cudaDeviceLmemResizeToMax
#endif
} cudaDeviceFlags;



#if CUDART_VERSION >= 3000
struct C2HS_COND_SENTRY_0;
#endif
#if CUDART_VERSION >= 3010
struct C2HS_COND_SENTRY_1;
#endif
#if CUDART_VERSION >= 4000
struct C2HS_COND_SENTRY_2;
#endif
#if CUDART_VERSION >= 3000
struct C2HS_COND_SENTRY_3;
#endif
#if CUDART_VERSION >= 3000
struct C2HS_COND_SENTRY_4;
#endif
#if CUDART_VERSION >= 3010
struct C2HS_COND_SENTRY_5;
#endif
#if CUDART_VERSION >= 3000 && CUDART_VERSION < 3010
struct C2HS_COND_SENTRY_6;
#endif
#if CUDART_VERSION >= 4000
struct C2HS_COND_SENTRY_7;
#endif
#if CUDART_VERSION < 4000
struct C2HS_COND_SENTRY_8;
#else
#endif
#if CUDART_VERSION >= 4000
struct C2HS_COND_SENTRY_9;
#else
#endif
#if CUDART_VERSION < 4000
struct C2HS_COND_SENTRY_10;
#else
#endif
#if CUDART_VERSION < 4000
struct C2HS_COND_SENTRY_11;
#else
#endif
#if CUDART_VERSION < 4000
struct C2HS_COND_SENTRY_12;
#else
#endif
#if CUDART_VERSION < 3010
struct C2HS_COND_SENTRY_13;
#else
#endif
#if   CUDART_VERSION < 3010
struct C2HS_COND_SENTRY_15;
#elif CUDART_VERSION < 4000
struct C2HS_COND_SENTRY_14;
#else
#endif
#if   CUDART_VERSION < 3010
struct C2HS_COND_SENTRY_17;
#elif CUDART_VERSION < 4000
struct C2HS_COND_SENTRY_16;
#else
#endif
