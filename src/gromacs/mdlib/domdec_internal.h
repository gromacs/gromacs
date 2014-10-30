/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2003 David van der Spoel, Erik Lindahl, University of Groningen.
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
/*! \file
 *  \brief
 * Declares implementation types for the domain decomposition module.
 */
#ifndef GMX_DOMDEC_DOMDEC_IMPL_H
#define GMX_DOMDEC_DOMDEC_IMPL_H


#include "gromacs/legacyheaders/types/commrec.h"
#include "gromacs/legacyheaders/domdec.h"


#define DDRANK(dd, rank)    (rank)
#define DDMASTERRANK(dd)   (dd->masterrank)

typedef struct gmx_domdec_master
{
    /* The cell boundaries */
    real **cell_x;
    /* The global charge group division */
    int   *ncg;    /* Number of home charge groups for each node */
    int   *index;  /* Index of nnodes+1 into cg */
    int   *cg;     /* Global charge group index */
    int   *nat;    /* Number of home atoms for each node. */
    int   *ibuf;   /* Buffer for communication */
    rvec  *vbuf;   /* Buffer for state scattering and gathering */
} gmx_domdec_master_t;

typedef struct
{
    /* The numbers of charge groups to send and receive for each cell
     * that requires communication, the last entry contains the total
     * number of atoms that needs to be communicated.
     */
    int  nsend[DD_MAXIZONE+2];
    int  nrecv[DD_MAXIZONE+2];
    /* The charge groups to send */
    int *index;
    int  nalloc;
    /* The atom range for non-in-place communication */
    int  cell2at0[DD_MAXIZONE];
    int  cell2at1[DD_MAXIZONE];
} gmx_domdec_ind_t;

typedef struct
{
    int               np;       /* Number of grid pulses in this dimension */
    int               np_dlb;   /* For dlb, for use with edlbAUTO          */
    gmx_domdec_ind_t *ind;      /* The indices to communicate, size np     */
    int               np_nalloc;
    gmx_bool          bInPlace; /* Can we communicate in place?            */
} gmx_domdec_comm_dim_t;

typedef struct
{
    gmx_bool *bCellMin;    /* Temp. var.: is this cell size at the limit     */
    real     *cell_f;      /* State var.: cell boundaries, box relative      */
    real     *old_cell_f;  /* Temp. var.: old cell size                      */
    real     *cell_f_max0; /* State var.: max lower boundary, incl neighbors */
    real     *cell_f_min1; /* State var.: min upper boundary, incl neighbors */
    real     *bound_min;   /* Temp. var.: lower limit for cell boundary      */
    real     *bound_max;   /* Temp. var.: upper limit for cell boundary      */
    gmx_bool  bLimited;    /* State var.: is DLB limited in this dim and row */
    real     *buf_ncd;     /* Temp. var.                                     */
} gmx_domdec_root_t;

#define DD_NLOAD_MAX 9

/* Here floats are accurate enough, since these variables
 * only influence the load balancing, not the actual MD results.
 */
typedef struct
{
    int    nload;
    float *load;
    float  sum;
    float  max;
    float  sum_m;
    float  cvol_min;
    float  mdf;
    float  pme;
    int    flags;
} gmx_domdec_load_t;

typedef struct
{
    int  nsc;
    int  ind_gl;
    int  ind;
} gmx_cgsort_t;

typedef struct
{
    gmx_cgsort_t *sort;
    gmx_cgsort_t *sort2;
    int           sort_nalloc;
    gmx_cgsort_t *sort_new;
    int           sort_new_nalloc;
    int          *ibuf;
    int           ibuf_nalloc;
} gmx_domdec_sort_t;

typedef struct
{
    rvec *v;
    int   nalloc;
} vec_rvec_t;

/* This enum determines the order of the coordinates.
 * ddnatHOME and ddnatZONE should be first and second,
 * the others can be ordered as wanted.
 */
enum {
    ddnatHOME, ddnatZONE, ddnatVSITE, ddnatCON, ddnatNR
};

typedef struct
{
    int      dim;       /* The dimension                                          */
    gmx_bool dim_match; /* Tells if DD and PME dims match                         */
    int      nslab;     /* The number of PME slabs in this dimension              */
    real    *slb_dim_f; /* Cell sizes for determining the PME comm. with SLB    */
    int     *pp_min;    /* The minimum pp node location, size nslab               */
    int     *pp_max;    /* The maximum pp node location,size nslab                */
    int      maxshift;  /* The maximum shift for coordinate redistribution in PME */
} gmx_ddpme_t;

typedef struct
{
    real min0;    /* The minimum bottom of this zone                        */
    real max1;    /* The maximum top of this zone                           */
    real min1;    /* The minimum top of this zone                           */
    real mch0;    /* The maximum bottom communicaton height for this zone   */
    real mch1;    /* The maximum top communicaton height for this zone      */
    real p1_0;    /* The bottom value of the first cell in this zone        */
    real p1_1;    /* The top value of the first cell in this zone           */
} gmx_ddzone_t;

/* In this struct we put all non-pointer elements at the start,
 * so we can communicate them with one MPI call.
 */
typedef struct {
    int      rank;              /* The rank we communicate with */
    gmx_bool bPBC;              /* Do we communicate over PBC? */
    gmx_bool bScrew;            /* Do we need screw PBC? */
    int      shift_ind;         /* Shift vector index for PBC */
    int      ncx, ncy;
    int      natoms;
    ivec     rCheckBonded;      /* The receiver should check bonded distances */
    rvec     corner0, corner1;
    real     size_x, size_y;
    int     *cxy_natoms;
    int      ncolumn;           /* Only for sending */
    int     *column_atom_range; /* Only for sending */
    int     *index_gl;          /* Global atom indices */
    rvec    *at_buf;            /* Atom vector send+recv buffer */
    int      cxy_nalloc;        /* Allocation size of cxy_natoms, column_atom_range */
    int      at_nalloc;         /* Allocation size of index_gl, at_buf */
} zone_comm_t;

typedef struct {
#ifdef GMX_MPI
    /* We store the requests of all zones in one array to pass to MPI_Waitall */
    MPI_Request req[DD_MAXZONE];
#endif
    int         nreq;
} zones_mpi_dir_t;

typedef struct {
    zones_mpi_dir_t recv;
    zones_mpi_dir_t send;
} zones_mpi_t;

typedef struct
{
    gmx_domdec_ind_t ind;
    int             *ibuf;
    int              ibuf_nalloc;
    vec_rvec_t       vbuf;
    int              nsend;
    int              nat;
    int              nsend_zone;
} dd_comm_setup_work_t;

typedef struct gmx_domdec_comm
{
    /* All arrays are indexed with 0 to dd->ndim (not Cartesian indexing),
     * unless stated otherwise.
     */

    /* The number of decomposition dimensions for PME, 0: no PME */
    int         npmedecompdim;
    /* The number of nodes doing PME (PP/PME or only PME) */
    int         npmenodes;
    int         npmenodes_x;
    int         npmenodes_y;
    /* The communication setup including the PME only nodes */
    gmx_bool    bCartesianPP_PME;
    ivec        ntot;
    int         cartpmedim;
    int        *pmenodes;          /* size npmenodes                         */
    int        *ddindex2simnodeid; /* size npmenodes, only with bCartesianPP
                                    * but with bCartesianPP_PME              */
    gmx_ddpme_t ddpme[2];

    /* The DD particle-particle nodes only */
    gmx_bool bCartesianPP;
    int     *ddindex2ddnodeid; /* size npmenode, only with bCartesianPP_PME */

    /* The global charge groups */
    t_block cgs_gl;

    /* Should we sort the cgs */
    int                nstSortCG;
    gmx_domdec_sort_t *sort;

    /* Are there charge groups? */
    gmx_bool bCGs;

    /* Are there bonded and multi-body interactions between charge groups? */
    gmx_bool bInterCGBondeds;
    gmx_bool bInterCGMultiBody;

    /* Data for the optional bonded interaction atom communication range */
    gmx_bool  bBondComm;
    t_blocka *cglink;
    /* TODO: When the group scheme is removed, this should go */
    char     *bLocalCG;

    /* The DLB option */
    int      eDLB;
    /* Is eDLB=edlbAUTO locked such that we currently can't turn it on? */
    gmx_bool bDLB_locked;
    /* Are we actually using DLB? */
    gmx_bool bDynLoadBal;

    /* Cell sizes for static load balancing, first index cartesian */
    real **slb_frac;

    /* The width of the communicated boundaries */
    real     cutoff_mbody;
    real     cutoff;
    /* The minimum cell size (including triclinic correction) */
    rvec     cellsize_min;
    /* For dlb, for use with edlbAUTO */
    rvec     cellsize_min_dlb;
    /* The lower limit for the DD cell size with DLB */
    real     cellsize_limit;
    /* Effectively no NB cut-off limit with DLB for systems without PBC? */
    gmx_bool bVacDLBNoLimit;

    /* With PME load balancing we set limits on DLB */
    gmx_bool bPMELoadBalDLBLimits;
    /* DLB needs to take into account that we want to allow this maximum
     * cut-off (for PME load balancing), this could limit cell boundaries.
     */
    real PMELoadBal_max_cutoff;

    /* tric_dir is only stored here because dd_get_ns_ranges needs it */
    ivec tric_dir;
    /* box0 and box_size are required with dim's without pbc and -gcom */
    rvec box0;
    rvec box_size;

    /* The cell boundaries */
    rvec cell_x0;
    rvec cell_x1;

    /* The old location of the cell boundaries, to check cg displacements */
    rvec old_cell_x0;
    rvec old_cell_x1;

    /* The communication setup and charge group boundaries for the zones */
    gmx_domdec_zones_t zones;

    /* The eighth shell x/f communication, backward & forward along the grid */
    zone_comm_t  zone_bw[DD_MAXZONE];
    zone_comm_t  zone_fw[DD_MAXZONE];
    zones_mpi_t  zones_mpi_x;
    zones_mpi_t  zones_mpi_f;

    /* The zone limits for DD dimensions 1 and 2 (not 0), determined from
     * cell boundaries of neighboring cells for staggered grids when using
     * dynamic load balancing.
     */
    gmx_ddzone_t zone_d1[2];
    gmx_ddzone_t zone_d2[2][2];

    /* The coordinate/force communication setup and indices */
    gmx_domdec_comm_dim_t cd[DIM];
    /* The maximum number of cells to communicate with in one dimension */
    int                   maxpulse;

    /* Should we check distances when assigning bonded interactions?
     * This is checked by the ranks we communicate with and reduced.
     */
    ivec rCheckBonded;

    /* Which cg distribution is stored on the master node */
    int master_cg_ddp_count;

    /* The number of cg's received from the direct neighbors */
    int  zone_ncg1[DD_MAXZONE];

    /* The atom counts, the range for each type t is nat[t-1] <= at < nat[t] */
    int  nat[ddnatNR];

    /* Array for signalling if atoms have moved to another domain */
    int  *moved;
    int   moved_nalloc;

    /* Communication buffer for general use */
    int  *buf_int;
    int   nalloc_int;

    /* Communication buffer for general use */
    vec_rvec_t vbuf;

    /* Temporary storage for thread parallel communication setup */
    int                   nth;
    dd_comm_setup_work_t *dth;

    /* Communication buffers only used with multiple grid pulses */
    int       *buf_int2;
    int        nalloc_int2;
    vec_rvec_t vbuf2;

    /* Communication buffers for local redistribution */
    int  **cggl_flag;
    int    cggl_flag_nalloc[DIM*2];
    rvec **cgcm_state;
    int    cgcm_state_nalloc[DIM*2];

    /* Cell sizes for dynamic load balancing */
    gmx_domdec_root_t **root;
    real               *cell_f_row;
    real                cell_f0[DIM];
    real                cell_f1[DIM];
    real                cell_f_max0[DIM];
    real                cell_f_min1[DIM];

    /* Stuff for load communication */
    gmx_bool           bRecordLoad;
    gmx_domdec_load_t *load;
    int                nrank_gpu_shared;
#ifdef GMX_MPI
    MPI_Comm          *mpi_comm_load;
    MPI_Comm           mpi_comm_gpu_shared;
#endif

    /* Maximum DLB scaling per load balancing step in percent */
    int dlb_scale_lim;

    /* Cycle counters */
    float  cycl[ddCyclNr];
    int    cycl_n[ddCyclNr];
    float  cycl_max[ddCyclNr];
    /* Flop counter (0=no,1=yes,2=with (eFlop-1)*5% noise */
    int    eFlop;
    double flop;
    int    flop_n;
    /* How many times did we have load measurements */
    int    n_load_have;
    /* How many times have we collected the load measurements */
    int    n_load_collect;

    /* Statistics */
    double sum_nat[ddnatNR-ddnatZONE];
    int    ndecomp;
    int    nload;
    double load_step;
    double load_sum;
    double load_max;
    ivec   load_lim;
    double load_mdf;
    double load_pme;

    /* The last partition step */
    gmx_int64_t partition_step;

    /* Debugging */
    int  nstDDDump;
    int  nstDDDumpGrid;
    int  DD_debug;
} gmx_domdec_comm_t;

/* The size per charge group of the cggl_flag buffer in gmx_domdec_comm_t */
#define DD_CGIBS 2

/* The flags for the cggl_flag buffer in gmx_domdec_comm_t */
#define DD_FLAG_NRCG  65535
#define DD_FLAG_FW(d) (1<<(16+(d)*2))
#define DD_FLAG_BW(d) (1<<(16+(d)*2+1))

/* Zone permutation required to obtain consecutive charge groups
 * for neighbor searching.
 */
static const int zone_perm[3][4] = { {0, 0, 0, 0}, {1, 0, 0, 0}, {3, 0, 1, 2} };

/* dd_zo and dd_zp3/dd_zp2 are set up such that i zones with non-zero
 * components see only j zones with that component 0.
 */

/* The DD zone order */
static const ivec dd_zo[DD_MAXZONE] =
{{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 1, 1}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1}};

/* The 3D setup */
#define dd_z3n  8
#define dd_zp3n 4
static const ivec dd_zp3[dd_zp3n] = {{0, 0, 8}, {1, 3, 6}, {2, 5, 6}, {3, 5, 7}};

/* The 2D setup */
#define dd_z2n  4
#define dd_zp2n 2
static const ivec dd_zp2[dd_zp2n] = {{0, 0, 4}, {1, 3, 4}};

/* The 1D setup */
#define dd_z1n  2
#define dd_zp1n 1
static const ivec dd_zp1[dd_zp1n] = {{0, 0, 2}};

/* Factors used to avoid problems due to rounding issues */
#define DD_CELL_MARGIN       1.0001
#define DD_CELL_MARGIN2      1.00005
/* Factor to account for pressure scaling during nstlist steps */
#define DD_PRES_SCALE_MARGIN 1.02

/* Turn on DLB when the load imbalance causes this amount of total loss.
 * There is a bit of overhead with DLB and it's difficult to achieve
 * a load imbalance of less than 2% with DLB.
 */
#define DD_PERF_LOSS_DLB_ON  0.02

/* Warn about imbalance due to PP or PP/PME load imbalance at this loss */
#define DD_PERF_LOSS_WARN    0.05

#define DD_CELL_F_SIZE(dd, di) ((dd)->nc[(dd)->dim[(di)]]+1+(di)*2+1+(di))

/* Use separate MPI send and receive commands
 * when nnodes <= GMX_DD_NNODES_SENDRECV.
 * This saves memory (and some copying for small nnodes).
 * For high parallelization scatter and gather calls are used.
 */
#define GMX_DD_NNODES_SENDRECV 4

#endif
