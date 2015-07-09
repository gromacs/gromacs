/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015, by the GROMACS development team, led by
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
#ifndef _commrec_h
#define _commrec_h

#include <stddef.h>

#include "gromacs/legacyheaders/types/commrec_fwd.h" // IWYU pragma: export
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/gmxmpi.h"
#include "gromacs/utility/real.h"

#ifdef __cplusplus
extern "C" {
#endif


#define DD_MAXZONE  8
#define DD_MAXIZONE 4

typedef struct gmx_domdec_master *gmx_domdec_master_p_t;

typedef struct {
    int  j0;     /* j-zone start               */
    int  j1;     /* j-zone end                 */
    int  cg1;    /* i-charge-group end         */
    int  jcg0;   /* j-charge-group start       */
    int  jcg1;   /* j-charge-group end         */
    ivec shift0; /* Minimum shifts to consider */
    ivec shift1; /* Maximum shifts to consider */
} gmx_domdec_ns_ranges_t;

typedef struct {
    rvec x0;     /* Zone lower corner in triclinic coordinates         */
    rvec x1;     /* Zone upper corner in triclinic coordinates         */
    rvec bb_x0;  /* Zone bounding box lower corner in Cartesian coords */
    rvec bb_x1;  /* Zone bounding box upper corner in Cartesian coords */
} gmx_domdec_zone_size_t;

struct gmx_domdec_zones_t {
    /* The number of zones including the home zone */
    int                    n;
    /* The shift of the zones with respect to the home zone */
    ivec                   shift[DD_MAXZONE];
    /* The charge group boundaries for the zones */
    int                    cg_range[DD_MAXZONE+1];
    /* The number of neighbor search zones with i-particles */
    int                    nizone;
    /* The neighbor search charge group ranges for each i-zone */
    gmx_domdec_ns_ranges_t izone[DD_MAXIZONE];
    /* Boundaries of the zones */
    gmx_domdec_zone_size_t size[DD_MAXZONE];
    /* The cg density of the home zone */
    real                   dens_zone0;
};

typedef struct gmx_ga2la *gmx_ga2la_t;

typedef struct gmx_hash *gmx_hash_t;

typedef struct gmx_reverse_top *gmx_reverse_top_p_t;

typedef struct gmx_domdec_constraints *gmx_domdec_constraints_p_t;

typedef struct gmx_domdec_specat_comm *gmx_domdec_specat_comm_p_t;

typedef struct gmx_domdec_comm *gmx_domdec_comm_p_t;

typedef struct gmx_pme_comm_n_box *gmx_pme_comm_n_box_p_t;

struct gmx_ddbox_t {
    int  npbcdim;
    int  nboundeddim;
    rvec box0;
    rvec box_size;
    /* Tells if the box is skewed for each of the three cartesian directions */
    ivec tric_dir;
    rvec skew_fac;
    /* Orthogonal vectors for triclinic cells, Cartesian index */
    rvec v[DIM][DIM];
    /* Normal vectors for the cells walls */
    rvec normal[DIM];
};


typedef struct {
    /* these buffers are used as destination buffers if MPI_IN_PLACE isn't
       supported.*/
    int             *ibuf; /* for ints */
    int              ibuf_alloc;

    gmx_int64_t     *libuf;
    int              libuf_alloc;

    float           *fbuf; /* for floats */
    int              fbuf_alloc;

    double          *dbuf; /* for doubles */
    int              dbuf_alloc;
} mpi_in_place_buf_t;


struct gmx_domdec_t {
    /* The DD particle-particle nodes only */
    /* The communication setup within the communicator all
     * defined in dd->comm in domdec.c
     */
    int                    nnodes;
    MPI_Comm               mpi_comm_all;
    /* Use MPI_Sendrecv communication instead of non-blocking calls */
    gmx_bool               bSendRecv2;
    /* The local DD cell index and rank */
    ivec                   ci;
    int                    rank;
    ivec                   master_ci;
    int                    masterrank;
    /* Communication with the PME only nodes */
    int                    pme_nodeid;
    gmx_bool               pme_receive_vir_ener;
    gmx_pme_comm_n_box_p_t cnb;
    int                    nreq_pme;
    MPI_Request            req_pme[8];


    /* The communication setup, identical for each cell, cartesian index */
    ivec     nc;
    int      ndim;
    ivec     dim; /* indexed by 0 to ndim */

    /* PBC from dim 0 to npbcdim */
    int npbcdim;

    /* Screw PBC? */
    gmx_bool bScrewPBC;

    /* Forward and backward neighboring cells, indexed by 0 to ndim */
    int  neighbor[DIM][2];

    /* Only available on the master node */
    gmx_domdec_master_p_t ma;

    /* Are there inter charge group constraints */
    gmx_bool bInterCGcons;
    gmx_bool bInterCGsettles;

    /* Global atom number to interaction list */
    gmx_reverse_top_p_t reverse_top;
    int                 nbonded_global;
    int                 nbonded_local;

    /* The number of inter charge-group exclusions */
    int  n_intercg_excl;

    /* Vsite stuff */
    gmx_hash_t                 ga2la_vsite;
    gmx_domdec_specat_comm_p_t vsite_comm;

    /* Constraint stuff */
    gmx_domdec_constraints_p_t constraints;
    gmx_domdec_specat_comm_p_t constraint_comm;

    /* The local to gobal charge group index and local cg to local atom index */
    int   ncg_home;
    int   ncg_tot;
    int  *index_gl;
    int  *cgindex;
    int   cg_nalloc;
    /* Local atom to local cg index, only for special cases */
    int  *la2lc;
    int   la2lc_nalloc;

    /* The number of home atoms */
    int   nat_home;
    /* The total number of atoms: home and received zones */
    int   nat_tot;
    /* Index from the local atoms to the global atoms */
    int  *gatindex;
    int   gatindex_nalloc;

    /* Global atom number to local atom number list */
    gmx_ga2la_t ga2la;

    /* Communication stuff */
    gmx_domdec_comm_p_t comm;

    /* The partioning count, to keep track of the state */
    gmx_int64_t ddp_count;


    /* gmx_pme_recv_f buffer */
    int   pme_recv_f_alloc;
    rvec *pme_recv_f_buf;

};

struct gmx_multisim_t {
    int       nsim;
    int       sim;
    MPI_Group mpi_group_masters;
    MPI_Comm  mpi_comm_masters;
    /* these buffers are used as destination buffers if MPI_IN_PLACE isn't
       supported.*/
    mpi_in_place_buf_t *mpb;
};

#define DUTY_PP  (1<<0)
#define DUTY_PME (1<<1)

typedef struct {
    int      bUse;
    MPI_Comm comm_intra;
    int      rank_intra;
    MPI_Comm comm_inter;

} gmx_nodecomm_t;

struct t_commrec {
    /* The nodeids in one sim are numbered sequentially from 0.
     * All communication within some simulation should happen
     * in mpi_comm_mysim, or its subset mpi_comm_mygroup.
     */
    int sim_nodeid, nnodes, npmenodes;

    /* thread numbers: */
    /* Not used yet: int threadid, nthreads; */
    /* The nodeid in the PP/PME, PP or PME group */
    int      nodeid;
    MPI_Comm mpi_comm_mysim;
    MPI_Comm mpi_comm_mygroup;

    /* MPI ranks within a physical node for hardware access */
    int            nrank_intranode;    /* nr of ranks on this physical node */
    int            rank_intranode;     /* our rank on this physical node */
    int            nrank_pp_intranode; /* as nrank_intranode, for particle-particle only */
    int            rank_pp_intranode;  /* as rank_intranode, for particle-particle only */

    gmx_nodecomm_t nc;

    /* For domain decomposition */
    struct gmx_domdec_t *dd;

    /* The duties of this node, see the defines above */
    int                    duty;

    struct gmx_multisim_t *ms;

    /* these buffers are used as destination buffers if MPI_IN_PLACE isn't
       supported.*/
    mpi_in_place_buf_t *mpb;
};

#define MASTERNODE(cr)     (((cr)->nodeid == 0) || !PAR(cr))
/* #define MASTERTHREAD(cr)   ((cr)->threadid == 0) */
/* #define MASTER(cr)         (MASTERNODE(cr) && MASTERTHREAD(cr)) */
#define MASTER(cr)         MASTERNODE(cr)
#define SIMMASTER(cr)      ((MASTER(cr) && ((cr)->duty & DUTY_PP)) || !PAR(cr))
#define NODEPAR(cr)        ((cr)->nnodes > 1)
/* #define THREADPAR(cr)      ((cr)->nthreads > 1) */
/* #define PAR(cr)            (NODEPAR(cr) || THREADPAR(cr)) */
#define PAR(cr)            NODEPAR(cr)
#define RANK(cr, nodeid)    (nodeid)
#define MASTERRANK(cr)     (0)

/* Note that even with particle decomposition removed, the use of
 * non-DD parallelization in TPI, NM and multi-simulations means that
 * PAR(cr) and DOMAINDECOMP(cr) are not universally synonymous. In
 * particular, DOMAINDECOMP(cr) == true indicates that there is more
 * than one domain, not just that the dd algorithm is active. */
#define DOMAINDECOMP(cr)   (((cr)->dd != NULL) && PAR(cr))
#define DDMASTER(dd)       ((dd)->rank == (dd)->masterrank)

#define MULTISIM(cr)       ((cr)->ms)
#define MSRANK(ms, nodeid)  (nodeid)
#define MASTERSIM(ms)      ((ms)->sim == 0)

/* The master of all (the node that prints the remaining run time etc.) */
#define MULTIMASTER(cr)    (SIMMASTER(cr) && (!MULTISIM(cr) || MASTERSIM((cr)->ms)))

#ifdef __cplusplus
}
#endif
#endif
