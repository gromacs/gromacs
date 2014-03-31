#include "types/simple.h"
#include "types/ishift.h"
#include "../nbnxn_consts.h"

#include "nbnxn_cuda_types.h"
#include "../../gmxlib/cuda_tools/cudautils.cuh"
#include "nbnxn_cuda.h"

#include "nbnxn_cuda_kernel_utils.cuh"

/* Top-level kernel generation: will generate through multiple
 * inclusion the following flavors for all kernel:
 * force and energy output without pair list pruning;
 */
#define PRUNE_NBL
#define CALC_ENERGIES
#include "nbnxn_cuda_kernels.cuh"
#undef CALC_ENERGIES
#undef PRUNE_NBL
