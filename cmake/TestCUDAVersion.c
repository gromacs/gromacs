#include <cuda_runtime_api.h>
int main()
{
/* This macro checks the generated GMX_CUDA_VERSION against the CUDA runtime
 * API version in cuda_runtime_api.h.
 */
#if CUDART_VERSION != GMX_CUDA_VERSION
#error CUDA version mismatch: CUDART_VERSION != GMX_CUDA_VERSION
#endif
    return 0;
}
