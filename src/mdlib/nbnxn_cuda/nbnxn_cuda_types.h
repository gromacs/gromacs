#ifndef CUTYPEDEFS_H
#define CUTYPEDEFS_H

#include "types/nbnxn_pairlist.h"
#include "types/nbnxn_cuda_types_ext.h"
#include "cudautils.cuh"

#ifdef __cplusplus
extern "C" {
#endif

/*! Types of electrostatics available in the CUDA nonbonded force kernels. */
enum {
    cu_eelEWALD, cu_eelRF, cu_eelCUT
};

typedef struct cu_nblist    cu_nblist_t;
typedef struct cu_atomdata  cu_atomdata_t;
typedef struct cu_nb_params cu_nb_params_t;
typedef struct cu_timers    cu_timers_t;
typedef struct gpu_tmp_data gpu_tmp_data_t;

/*! Staging area for temporary data. The energies get downloaded here first, 
 *  before getting added to the CPU-side aggregate values.
 */
struct nb_tmp_data
{
    float   *e_lj;      /* LJ energy */
    float   *e_el;      /* electrostatic energy */
    float4  *f_shift;   /* shift forces */
};

/*! Nonbonded atom data */
struct cu_atomdata
{
    int     natoms;             /* number of atoms                      */
    int     natoms_local;       /* number of local atoms                */
    int     nalloc;             /* allocation size for the atom data (xq, f) */
    
    float4  *xq;                /* atom coordinates + charges, size natoms  */
    float4  *f;                 /* force output array, size natoms          */
    /* TODO: try float2 for the energies */
    float   *e_lj,              /* LJ energy output, size 1                 */
            *e_el;              /* Electrostatics energy intput, size 1     */

    float4  *f_shift;           /* shift forces */

    int     ntypes;             /* number of atom types             */
    int     *atom_types;        /* atom type indices, size natoms   */
 
    float3  *shift_vec;         /* shifts */
    gmx_bool shift_vec_copied;  /* has the shift vector already been transfered? */
};

/*! Nonbonded paramters */
struct cu_nb_params
{
    int     eeltype;        /* type of electrostatics */
    
    float   epsfac;         /* charge multiplication factor */
    float   c_rf, two_k_rf; /* Reaction-Field constants */
    float   ewald_beta;     /* Ewald/PME parameter */
    float   rvdw_sq;        /* VdW cut-off */
    float   rcoulomb_sq;    /* Coulomb cut-off */
    float   rlist_sq;       /* neighborlist cut-off */
    float   lj_shift;       /* LJ potential correction term */

    float   *nbfp;          /* nonbonded parameter table with C6/C12 pairs  */

    /* Ewald Coulomb force table */
    int     coulomb_tab_size;
    float   coulomb_tab_scale;
    float   *coulomb_tab;
};

/*! Neighbor list data */
struct cu_nblist 
{
    int             na_c;       /* number of atoms per cluster               */
    
    int             nsci;       /* size of sci, # of i clusters in the list  */
    int             sci_nalloc; /* allocation size for sci                   */
    nbnxn_sci_t     *sci;       /* list of i-cluster ("superclusters")       */

    int             ncj4;       /* total # of i-j cell subcell pairs         */
    int             cj4_nalloc; /* allocation size for sj                    */
    nbnxn_cj4_t     *cj4;       /* j cluster list, contains j cluster number 
                                   and index into the i cluster list         */
    nbnxn_excl_t    *excl;      /* Exclusions                                */
    int             nexcl;      /* The count for excl                        */
    int             excl_nalloc; /* The allocation size for excl      */

    gmx_bool        prune_nbl;  /* true if neighbor list pruning needs to be 
                                   done during the  current step                */
};

/* CUDA events for timing GPU kernels and H2D/D2H transfers and for 
 * synchronization purposes (not neded for this purpose ATM). 
 * The two-sized arrays hold the local and non-local values (respectively)
 * and should always be indexed with eintLocal or eintNonlocal.
 */
struct cu_timers
{
    cudaEvent_t start_atdat, stop_atdat;            /* events for atom data transfer (every NS step)    */

    cudaEvent_t start_nb_h2d[2], stop_nb_h2d[2];    /* events for H2D transfer (every step)             */
    cudaEvent_t start_nb_d2h[2], stop_nb_d2h[2];    /* events for D2H transfer (every step)             */
    cudaEvent_t start_nbl_h2d[2], stop_nbl_h2d[2];  /* events for pair-list H2D strnasfer (every nstlist step) */

    cudaEvent_t  start_nb_k[2], stop_nb_k[2];       /* events for non-bonded kernels                        */
};

/* Main data structure for CUDA nonbonded force calculations. */
struct cu_nonbonded 
{
    cu_dev_info_t   *dev_info;
    gmx_bool        use_stream_sync; /* if true use memory polling-based waiting instread 
                                        of cudaStreamSynchronize */
    cu_atomdata_t   *atomdata;
    cu_nb_params_t  *nb_params; 
    cu_nblist_t     *nblist[2]; /* nbl data structures, local and non-local (only with DD) */
    nb_tmp_data     tmpdata;
    
    cudaStream_t    stream[2];  /* local and non-local GPU streams */

    gmx_bool        do_time;    /* true if CUDA event-based timing is enabled, off with DD */
    cu_timers_t     *timers;
    cu_timings_t    *timings;
};

#ifdef __cplusplus
}
#endif

#endif	/* _CUTYPEDEFS_H_ */
