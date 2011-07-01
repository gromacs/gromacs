#ifndef CUTYPEDEFS_H
#define CUTYPEDEFS_H

#include "types/nblist_box.h"
#include "cutypedefs_ext.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Types of electrostatics available in the CUDA nonbonded force kernels. */
enum {
    cu_eelEWALD, cu_eelRF, cu_eelCUT
};

typedef struct cu_nblist    cu_nblist_t;
typedef struct cu_atomdata  cu_atomdata_t;
typedef struct cu_nb_params cu_nb_params_t;
typedef struct cu_timers    cu_timers_t;
typedef struct gpu_tmp_data gpu_tmp_data_t;

/* Staging area for temporary data. The energies get downloaded here first, 
 * before getting added to the CPU-side aggregate values.
 * */
struct nb_tmp_data
{
    float   *e_lj;      /* LJ energy */
    float   *e_el;      /* electrostatic energy */
    float4  *f_shift;   /* shift forces */
};

struct cu_atomdata
{
    int     natoms;         /* number of atoms                      */
    int     natoms_local;   /* number of local atoms                */
    int     nalloc;     /* allocation size for the atom data (xq, f) */ 
    
    float4  *xq;        /* atom coordinates + charges, size natoms  */
    float4  *f;         /* force output array, size natoms          */
    /* TODO: try float2 for the energies */
    float   *e_lj,      /* LJ energy output, size 1                 */
            *e_el;      /* Electrostatics energy intput, size 1     */

    float4  *f_shift;   /* shift forces */

    int     ntypes;     /* number of atom types             */
    int     *atom_types;/* atom type indices, size natoms   */
 
    float3  *shift_vec;  /* shifts */    
    gmx_bool shift_vec_copied;   /* indicates whether shift vector has already been transfered */
};

/* nonbonded paramters */
struct cu_nb_params
{
    int  eeltype;       /* type of electrostatics */ 
    
    float   eps_r; /* TODO rename this to epsfac  and get rid of the FACEL constant */
    float   c_rf;       
    float   two_k_rf;   
    float   ewald_beta; 
    float   cutoff_sq;  /* cut-off */
    float   rlist_sq;   /* neighborlist cut-off */
    float   lj_shift;   /* LJ potential correction term */

    float   *nbfp;      /* nonbonded parameters C12, C6 */

    /* Ewald Coulomb force table */
    int     coulomb_tab_size;
    float   coulomb_tab_scale;
    float   *coulomb_tab;
};

/* neighbor list data */
struct cu_nblist 
{
    int             naps;       /* number of atoms per subcell                  */
    
    int             nci;        /* size of ci, # of i cells in the list         */
    int             ci_nalloc;  /* allocation size for ci                       */
    gmx_nbl_ci_t     *ci;       /* list of i-cells ("supercells")               */

    int             nsj4;       /* total # of i-j cell subcell pairs           */
    int             sj4_nalloc; /* allocation size for sj                      */
    gmx_nbl_sj4_t   *sj4;       /* j subcell list, contains j subcell number 
                                    and index into the i subcell list           */
    gmx_nbl_excl_t  *excl;      /* Exclusions                               */
    int             nexcl;      /* The count for excl                       */
    int             excl_nalloc;/* The allocation size for excl             */

    gmx_bool        prune_nbl;  /* true if neighbor list pruning needs to be 
                                   done during the  current step                */
};

/* timers for the GPU kernels and H2D/D2H trasfers */
struct cu_timers
{
    cudaEvent_t start_nb, stop_nb;          /* events for timing nonbonded calculation + related 
                                                   data transfers                                           */
    cudaEvent_t start_nb_nl, stop_nb_nl;    /* events for timing nonbonded calculation on non-local data
                                                   + related data transfers                                 */   
    gmx_bool    time_transfers;             /* enable/disable separate host-device data trasnfer timing     */
    cudaEvent_t start_nb_h2d, stop_nb_h2d;  /* events for timing host to device transfer (every step)       */
    cudaEvent_t start_nb_h2d_nl, stop_nb_h2d_nl;  /* events for timing host to device transfer (every step) */
    cudaEvent_t start_nb_d2h, stop_nb_d2h;  /* events for timing device to host transfer (every step)       */
    cudaEvent_t start_nb_d2h_nl, stop_nb_d2h_nl; /* events for timing device to host transfer (every step) */
    cudaEvent_t start_atdat, stop_atdat;    /* events for timing atom data transfer (every NS step)         */
    cudaEvent_t start_atdat_nl, stop_atdat_nl;    /* events for timing atom data transfer (every NS step)   */

    cudaStream_t    nbstream, nbstream_nl;      /* local and non-local calculation streams */
};

/* main data structure for CUDA nonbonded force evaluation */
struct cu_nonbonded 
{
    gmx_bool        streamGPU;  /* XXX are we streaming or not (debugging)              */
    
    cu_atomdata_t   *atomdata;
    cu_nb_params_t  *nb_params; 
    cu_nblist_t     *nblist; 
    cu_nblist_t     *nblist_nl; 
    cu_timers_t     *timers;
    cu_timings_t    *timings;
    nb_tmp_data     tmpdata;
};

#ifdef __cplusplus
}
#endif

#endif	/* _CUTYPEDEFS_H_ */
