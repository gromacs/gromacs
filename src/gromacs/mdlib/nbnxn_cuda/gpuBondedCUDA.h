#include "gromacs/mdlib/nbnxn_gpu_types.h"
#include "gromacs/ewald/pme-internal.h"
#include "gromacs/ewald/pme.h"
#include "gromacs/ewald/pme-grid.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/enerdata.h"

//#include "gmxpre.h"
#include "gromacs/mdlib/nbnxn_atomdata.h"
//#include "gromacs/mdlib/nbnxn_consts.h"
//#include "gromacs/mdlib/nbnxn_internal.h"
//#include "gromacs/mdlib/nbnxn_search.h"
//#include "gromacs/mdlib/nbnxn_util.h"

#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/pbcutil/mshift.h"

//CUDA threads per block
#define TPB_BONDED 256

//Macro for checking cuda errors following a cuda launch or api call
#define cudaCheckError() {                                          \
    cudaError_t e=cudaGetLastError();                                 \
    if(e!=cudaSuccess) {                                              \
      printf("Cuda failure %s:%d: '%s'\n",__FILE__,__LINE__,cudaGetErrorString(e));           \
      exit(0); \
    }                                                                 \
  }
  
#if defined(__CUDACC__)
#include "nbnxn_cuda.h"
#include "nbnxn_cuda_types.h"
#include "cuda.h"
#include "cuda_runtime.h"
#endif

// new stuff for bonded
void update_gpu_bonded( const t_idef *idef,  const t_forcerec *fr, matrix box,
                        int size,  const t_mdatoms *md, const real *lambda,
                        gmx_grppairener_t *grppener );                       
void do_bonded_gpu(t_forcerec *fr, const t_inputrec *ir, const t_idef *idef, 
                   int flags, const t_graph *graph , int natoms, rvec x[], 
                   real *lambda, const t_mdatoms *md, 
                   rvec *input_force, t_lambda *fepvals, gmx_enerdata_t *enerd);

void do_bonded_gpu_finalize(t_forcerec *fr, const t_inputrec *ir, const t_idef *idef,
                   int flags, const t_graph *graph , int natoms, rvec x[],
                   real *lambda, const t_mdatoms *md,
                   rvec *input_force, t_lambda *fepvals, gmx_enerdata_t *enerd);

void reset_gpu_bonded(const int size, const int nener);

#define STANDALONE

//#define DEBUG_JV




