#ifndef CUTYPEDEFS_EXT_H
#define CUTYPEDEFS_EXT_H

#define GPU_NS_CELL_SIZE    8

#define GPU_MIN_CI_BALANCED_FACTOR 40

#ifdef __cplusplus
extern "C" {
#endif
/* Abstract types for CUDA nonbonded and device info data structures */
typedef struct cu_nonbonded * cu_nonbonded_t;

typedef struct cu_timings cu_timings_t;
typedef struct nb_kernel_time nb_kernel_time_t;

struct nb_kernel_time
{
    double  t;
    int     c;
}; 

struct cu_timings 
{
    nb_kernel_time_t k_time[2][2]; /* table containing the timings of the four 
                                      version of the nonbonded kernels: force-only, 
                                      force+energy, force+pruning, and force+energy+pruning */
    double  nb_h2d_time;    /* host to device transfer time of data */
    double  nb_d2h_time;    /* device to host transfer time of data */
    int     nb_count;       /* total call count of the nonbonded gpu operations */

    double  nbl_h2d_time;   /* total time of the data trasnfer after a neighbor search step */
    int     nbl_h2d_count;  /* cll count   - || - */
};

#ifdef __cplusplus
}
#endif

#endif
