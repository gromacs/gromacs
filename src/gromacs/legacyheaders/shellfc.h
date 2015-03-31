/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team.
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

#include "../timing/wallcycle.h"

#include "../pbcutil/mshift.h"
#include "../pbcutil/ishift.h"
#include "../pbcutil/pbc.h"

#include "typedefs.h"
#include "vsite.h"
#include "update.h"

#ifndef _shellfc_h
#define _shellfc_h

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int         nnucl;
    atom_id     shell;               /* The shell id                */
    atom_id     nucl1, nucl2, nucl3; /* The nuclei connected to the shell   */
    /* gmx_bool    bInterCG; */            /* Coupled to nuclei outside cg?        */
    real        k;                   /* force constant              */
    real        k_1;                 /* 1 over force constant       */
    rvec        xold;
    rvec        fold;
    rvec        step;
    real        k11, k22, k33;       /* anisotropic polarization force constants */
} t_shell;

typedef struct {
    t_ilist ilist[F_NRE];     /* shell ilists for this thread            */
    rvec    fshift[SHIFTS];   /* fshift accumulation buffer              */
    matrix  dxdf;             /* virial dx*df accumulation buffer        */
} gmx_shell_thread_t;

typedef struct gmx_shellfc {
    int                 nshell_gl;            /* The number of shells in the system       */
    t_shell            *shell_gl;             /* All the shells (for DD only)             */
    int                *shell_index_gl;       /* Global shell index (for DD only)         */
    int                 nshell;               /* The number of local shells               */
    t_shell            *shell;                /* The local shells                         */
    int                 shell_nalloc;         /* The allocation size of shell             */
    gmx_bool            bPredict;             /* Predict shell positions                  */
    gmx_bool            bRequireInit;         /* Require initialization of shell positions  */
    int                 nflexcon;             /* The number of flexible constraints       */
    rvec               *x[2];                 /* Array for iterative minimization         */
    rvec               *f[2];                 /* Array for iterative minimization         */
    int                 x_nalloc;             /* The allocation size of x and f           */
    rvec               *acc_dir;              /* Acceleration direction for flexcon       */
    rvec               *x_old;                /* Old coordinates for flexcon              */
    int                 flex_nalloc;          /* The allocation size of acc_dir and x_old */
    rvec               *adir_xnold;           /* Work space for init_adir                 */
    rvec               *adir_xnew;            /* Work space for init_adir                 */
    int                 adir_nalloc;          /* Work space for init_adir                 */
    gmx_bool            bInterCG;             /* Are there inter charge-group shells?     */
    int                 n_intercg_shells;     /* inter-charge group shells                */
    int                 nshell_pbc_molt;      /* The array size of shell_pbc_molt         */
    int              ***shell_pbc_molt;       /* The pbc atoms for intercg shells         */
    int               **shell_pbc_loc;        /* The local pbc atoms                      */
    int                *shell_pbc_loc_nalloc; /* Sizes of shell_pbc_loc                   */
    int                 nthreads;             /* Number of threads used for shells        */
    gmx_shell_thread_t *tdata;                /* Thread local shells and work structs     */
    int                *th_ind;               /* Work array                               */
    int                 th_ind_nalloc;        /* Size of th_ind                           */
} t_gmx_shellfc;

/* Abstract type for shells */
typedef struct gmx_shellfc *gmx_shellfc_t;

struct t_graph;
struct t_pbc;

/* Initialization function, also predicts the initial shell postions.
 * If x!=NULL, the shells are predict for the global coordinates x.
 */
/* TODO: restructuring */
void init_shell_flexcon(FILE *fplog, gmx_shellfc_t shfc, t_inputrec *ir,
                        gmx_mtop_t *mtop, int nflexcon,
                        rvec *x);

/* Get the local shell with domain decomposition */
void make_local_shells(t_commrec *cr, t_mdatoms *md,
                       gmx_shellfc_t shfc);

/* Optimize shell positions */
void relax_shell_flexcon(FILE *log, t_commrec *cr, gmx_bool bVerbose,
                        gmx_int64_t mdstep, t_inputrec *inputrec,
                        gmx_bool bDoNS, int force_flags,
                        gmx_localtop_t *top,
                        gmx_constr_t constr,
                        gmx_enerdata_t *enerd, t_fcdata *fcd,
                        t_state *state, rvec f[],
                        tensor force_vir,
                        t_mdatoms *md,
                        t_nrnb *nrnb, gmx_wallcycle_t wcycle,
                        struct t_graph *graph,
                        gmx_groups_t *groups,
                        gmx_shellfc_t shfc,
                        t_forcerec *fr,
                        gmx_bool bBornRadii,
                        double t, rvec mu_tot,
                        gmx_bool *bConverged,
                        gmx_vsite_t *vsite,
                        FILE *fp_field,
                        int *count);

/* functions for DD */
void add_quartic_restraint_force(t_inputrec *ir, gmx_shellfc_t shfc, rvec x[], rvec f[]);

void apply_drude_hardwall(t_commrec *cr, t_idef *idef, t_inputrec *ir, t_mdatoms *md,       
                          t_state *state, tensor force_vir);

static void spread_shell(t_iatom ia[],
                         rvec x[], rvec f[], rvec fshift[],
                         t_pbc *pbc, t_graph *g);

static void spread_shell_f_thread(gmx_shellfc_t shell,
                                  rvec x[], rvec f[], rvec *fshift,
                                  gmx_bool VirCorr, matrix dxdf,
                                  t_ilist ilist[],
                                  t_graph *g, t_pbc *pbc_null);

void spread_shell_f(gmx_shellfc_t shell,
                    rvec x[], rvec f[], rvec *fshift,
                    gmx_bool VirCorr, matrix vir,
                    t_nrnb *nrnb, t_idef *idef,
                    int ePBC, gmx_bool bMolPBC, t_graph *g, matrix box,
                    t_commrec *cr);

gmx_shellfc_t init_shell(gmx_mtop_t *mtop, t_commrec *cr,
                         gmx_bool bSerial_NoPBC);

void split_shells_over_threads(const t_ilist   *ilist,
                               const t_mdatoms *mdatoms,
                               gmx_bool         bLimitRange,
                               gmx_shellfc_t    shfc);

void set_shell_top(gmx_shellfc_t shfc, gmx_localtop_t *top, t_mdatoms *md,
                   t_commrec *cr);

#ifdef __cplusplus
}
#endif

#endif
