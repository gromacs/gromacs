#ifndef CUTYPEDEFS_EXT_H
#define CUTYPEDEFS_EXT_H

#define GPU_NS_CELL_SIZE    8

#ifdef __cplusplus
extern "C" {
#endif

typedef struct cudata * t_cudata;
typedef struct gpu_times gpu_times_t;

struct nb_kernel_time
{
    double  t;
    int     c;
}; 

typedef struct nb_kernel_time nb_kernel_time_t;


struct gpu_times
{
    /*float   nb_total_time;  *//* total execution time of the nonbonded gpu operations:
                               - trasfer to/from GPU: x, q, shifts, f
                               - kernel exection */
    nb_kernel_time_t k_time[2][2]; /* TODO */
    double  nb_h2d_time;    /* host to device transfer time of data */
    double  nb_d2h_time;    /* device to host transfer time of data */
    int     nb_count;       /* total call count of the nonbonded gpu operations */

    double  atomdt_h2d_total_time;  /* total time of the data trasnfer after a neighbor search step */
    int     atomdt_count;           /* cll count   - || - */
    
};

#ifdef __cplusplus
}
#endif

#endif
