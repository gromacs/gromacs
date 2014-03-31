#include "gromacs/gmxlib/cuda_tools/cudautils.cuh"
#include "gromacs/legacyheaders/types/simple.h"
#include "gromacs/mdlib/nbnxn_consts.h"
#include "gromacs/pbcutil/ishift.h"

#include "nbnxn_cuda_types.h"
#include "nbnxn_cuda.h"

#include "nbnxn_cuda_kernel_utils.cuh"

/* Top-level kernel generation: will generate through multiple
 * inclusion the following flavors for all kernel:
 * force-only output without pair list pruning;
 */
#include "nbnxn_cuda_kernels.cuh"
