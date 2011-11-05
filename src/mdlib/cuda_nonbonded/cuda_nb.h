#ifndef GPU_NB_H
#define GPU_NB_H

#include "cutypedefs_ext.h"

#ifdef __cplusplus
extern "C" {
#endif

/*! Launch asynchronously the nonbonded force calculations. 
 *  This consists of the following (async) steps launched: 
 *  - upload x and q;
 *  - upload shift vector;
 *  - launch kernel;
 *  The local and non-local interation calculations are launched in two 
 *  separate streams.
 */
void cu_nb_launch_kernel(cu_nonbonded_t /*cu_nb*/, 
                         const gmx_nb_atomdata_t * /*nbdata*/,                  
                         int /*flags*/,
                         int /*iloc*/);

/*! Launch asynchronously the dowload of nonbonded forces from the GPU 
 *  (also energies/shift forces if required). 
 */
void cu_nb_launch_cpyback(cu_nonbonded_t /*cu_nb*/,
                          const gmx_nb_atomdata_t * /*nbatom*/,
                          int /*flags*/,
                          int /*aloc*/);

/*! Wait for the asynchrounously launched nonbonded calculations and data 
 *  transfers to finish. 
 */
void cu_nb_wait_gpu(cu_nonbonded_t /*cu_nb*/,
                    const gmx_nb_atomdata_t * /*nbatom*/,
                    int /*flags*/, int /*aloc*/,
                    float * /*e_lj*/, float * /*e_el*/, rvec * /*fshift*/);

#ifdef __cplusplus
}
#endif

#endif
