#include <cuda_runtime_api.h>
int main()
{
/* This macro checks the generated GMX_CUDA_VERSION against the CUDA runtime
 * API version in cuda_runtime_api.h.
 *
 * From CUDA v7.5 it is expected that nvcc will define its own version; a check
 * of that version should be implemented here later.
 */
#if CUDART_VERSION != GMX_CUDA_VERSION
#error CUDA version mismatch: CUDART_VERSION != GMX_CUDA_VERSION
#endif
    return 0;
}
