#include <cuda.h>
int main()
{
/* This macro checks the generated GMX_CUDA_VERSION against the value of the
 * CUDA_VERSION macro defined in cuda.h.
 */
#if CUDA_VERSION != GMX_CUDA_VERSION
#error CUDA version mismatch: CUDA_VERSION != GMX_CUDA_VERSION
#endif
    return 0;
}
