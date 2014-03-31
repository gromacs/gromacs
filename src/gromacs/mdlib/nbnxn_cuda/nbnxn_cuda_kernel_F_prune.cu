#include "types/simple.h"
#include "types/ishift.h"
#include "../nbnxn_consts.h"

#include "nbnxn_cuda_types.h"
#include "../../gmxlib/cuda_tools/cudautils.cuh"
#include "nbnxn_cuda.h"

#include "nbnxn_cuda_kernel_utils.cuh"

/* Top-level kernel generation: will generate through multiple
 * inclusion the following flavors for all kernel:
 * force-only output with pair list pruning;
 */
#define PRUNE_NBL
#include "nbnxn_cuda_kernels.cuh"
#undef PRUNE_NBL
