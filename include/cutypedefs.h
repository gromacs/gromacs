#ifndef CUTYPEDEFS_H
#define CUTYPEDEFS_H

#include "types/nblist_box.h"
#include "cutypedefs_ext.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Types of electrostatics available in the CUDA GPU imlementation. */
enum {
    cu_eelEWALD, cu_eelRF, cu_eelCUT
};

typedef struct gpu_tmp_data gpu_tmp_data_t;

/* Staging area and temporary data for GPU ops */
struct gpu_tmp_data
{
    float *e_lj;
    float *e_el; 
};

struct cudata 
{
    int     natoms;     /* number of atoms for all 8 neighbouring domains 
                           !!! ATM only one value, with MPI it'll be 8                      */
    int     nalloc;     /* allocation size for the atom data (xq, f), 
                           when needed it's reallocated to natoms * 20% + 100 buffer zone   */ 
    
    float4  *xq;        /* atom coordinates + charges, size natoms  */
    float4  *f;         /* forces, size natoms                      */
    /* TODO: try float2 for the energies */
    float   *e_lj,      /* LJ energy                                */
            *e_el;      /* Electrostatics energy                    */

    int     ntypes;     /* number of atom types             */
    int     *atom_types;/* atom type indices, size natoms   */
    
    /* nonbonded paramters 
       TODO -> constant so some of them should be moved to constant memory */
    float   eps_r; /* TODO rename this to epsfac  and get rid of the FACEL constant */
    float   c_rf;
    float   two_k_rf;
    float   ewald_beta;
    float   cutoff_sq;
    float   rlist_sq;   /* neighborlist cut-off */
    float   *nbfp;      /* nonbonded parameters C12, C6 */

    int  eeltype;       /* type of electrostatics */ 

    /* Ewald Coulomb tabulated force */
    int     coulomb_tab_size;
    float   coulomb_tab_scale;
    float   *coulomb_tab;

    /* async execution stuff */
    gmx_bool        streamGPU;                  /* are we streaming of not (debugging)              */
    cudaStream_t    nb_stream;                  /* XXX nonbonded calculation stream - not in use    */
    cudaEvent_t     start_nb, stop_nb;          /* events for timing nonbonded calculation + related 
                                                   data transfers                                   */
    gmx_bool        time_transfers;             /* enable/disable separate host-device data trasnfer timing */
    cudaEvent_t     start_nb_h2d, stop_nb_h2d;  /* events for timing host to device transfer (every step) */
    cudaEvent_t     start_nb_d2h, stop_nb_d2h;  /* events for timing device to host transfer (every step) */
    cudaEvent_t     start_atdat, stop_atdat;    /* events for timing atom data transfer (every NS step) */

    /* neighbor list data */
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

    float3          *shift_vec;     /* shifts */    
    gmx_bool    shift_vec_copied;   /* indicates whether shift vector has already been transfered */

    gpu_tmp_data_t  tmpdata;    
    gpu_times_t     timings;
};


#ifdef __cplusplus
}
#endif


#endif	/* _CUTYPEDEFS_H_ */
