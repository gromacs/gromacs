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

#ifndef GMX_LEGACYHEADERS_TYPES_FORCEREC_H
#define GMX_LEGACYHEADERS_TYPES_FORCEREC_H

#include "gromacs/legacyheaders/types/enums.h"
#include "gromacs/legacyheaders/types/genborn.h"
#include "gromacs/legacyheaders/types/hw_info.h"
#include "gromacs/legacyheaders/types/interaction_const.h"
#include "gromacs/legacyheaders/types/nblist.h"
#include "gromacs/legacyheaders/types/ns.h"
#include "gromacs/legacyheaders/types/qmmmrec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/topology/idef.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

#ifdef __cplusplus
extern "C" {
#endif
#if 0
} /* fixes auto-indentation problems */
#endif

/* Abstract type for PME that is defined only in the routine that use them. */
struct gmx_pme_t;
struct nonbonded_verlet_t;

/* Structure describing the data in a single table */
typedef struct
{
    enum gmx_table_interaction  interaction; /* Types of interactions stored in this table */
    enum gmx_table_format       format;      /* Interpolation type and data format */

    real                        r;           /* range of the table */
    int                         n;           /* n+1 is the number of table points */
    real                        scale;       /* distance (nm) between two table points */
    real                        scale_exp;   /* distance for exponential part of VdW table, not always used */
    real *                      data;        /* the actual table data */

    /* Some information about the table layout. This can also be derived from the interpolation
     * type and the table interactions, but it is convenient to have here for sanity checks, and it makes it
     * much easier to access the tables in the nonbonded kernels when we can set the data from variables.
     * It is always true that stride = formatsize*ninteractions
     */
    int                         formatsize;    /* Number of fp variables for each table point (1 for F, 2 for VF, 4 for YFGH, etc.) */
    int                         ninteractions; /* Number of interactions in table, 1 for coul-only, 3 for coul+rep+disp. */
    int                         stride;        /* Distance to next table point (number of fp variables per table point in total) */
} t_forcetable;

typedef struct
{
    t_forcetable   table_elec;
    t_forcetable   table_vdw;
    t_forcetable   table_elec_vdw;

    /* The actual neighbor lists, short and long range, see enum above
     * for definition of neighborlist indices.
     */
    t_nblist nlist_sr[eNL_NR];
    t_nblist nlist_lr[eNL_NR];
} t_nblists;

/* macros for the cginfo data in forcerec
 *
 * Since the tpx format support max 256 energy groups, we do the same here.
 * Note that we thus have bits 8-14 still unused.
 *
 * The maximum cg size in cginfo is 63
 * because we only have space for 6 bits in cginfo,
 * this cg size entry is actually only read with domain decomposition.
 * But there is a smaller limit due to the t_excl data structure
 * which is defined in nblist.h.
 */
#define SET_CGINFO_GID(cgi, gid)     (cgi) = (((cgi)  &  ~255) | (gid))
#define GET_CGINFO_GID(cgi)        ( (cgi)            &   255)
#define SET_CGINFO_FEP(cgi)          (cgi) =  ((cgi)  |  (1<<15))
#define GET_CGINFO_FEP(cgi)        ( (cgi)            &  (1<<15))
#define SET_CGINFO_EXCL_INTRA(cgi)   (cgi) =  ((cgi)  |  (1<<16))
#define GET_CGINFO_EXCL_INTRA(cgi) ( (cgi)            &  (1<<16))
#define SET_CGINFO_EXCL_INTER(cgi)   (cgi) =  ((cgi)  |  (1<<17))
#define GET_CGINFO_EXCL_INTER(cgi) ( (cgi)            &  (1<<17))
#define SET_CGINFO_SOLOPT(cgi, opt)  (cgi) = (((cgi)  & ~(3<<18)) | ((opt)<<18))
#define GET_CGINFO_SOLOPT(cgi)     (((cgi)>>18)       &   3)
#define SET_CGINFO_CONSTR(cgi)       (cgi) =  ((cgi)  |  (1<<20))
#define GET_CGINFO_CONSTR(cgi)     ( (cgi)            &  (1<<20))
#define SET_CGINFO_SETTLE(cgi)       (cgi) =  ((cgi)  |  (1<<21))
#define GET_CGINFO_SETTLE(cgi)     ( (cgi)            &  (1<<21))
/* This bit is only used with bBondComm in the domain decomposition */
#define SET_CGINFO_BOND_INTER(cgi)   (cgi) =  ((cgi)  |  (1<<22))
#define GET_CGINFO_BOND_INTER(cgi) ( (cgi)            &  (1<<22))
#define SET_CGINFO_HAS_VDW(cgi)      (cgi) =  ((cgi)  |  (1<<23))
#define GET_CGINFO_HAS_VDW(cgi)    ( (cgi)            &  (1<<23))
#define SET_CGINFO_HAS_Q(cgi)        (cgi) =  ((cgi)  |  (1<<24))
#define GET_CGINFO_HAS_Q(cgi)      ( (cgi)            &  (1<<24))
#define SET_CGINFO_NATOMS(cgi, opt)  (cgi) = (((cgi)  & ~(63<<25)) | ((opt)<<25))
#define GET_CGINFO_NATOMS(cgi)     (((cgi)>>25)       &   63)


/* Value to be used in mdrun for an infinite cut-off.
 * Since we need to compare with the cut-off squared,
 * this value should be slighlty smaller than sqrt(GMX_FLOAT_MAX).
 */
#define GMX_CUTOFF_INF 1E+18

/* enums for the neighborlist type */
enum {
    enbvdwNONE, enbvdwLJ, enbvdwBHAM, enbvdwTAB, enbvdwNR
};
/* OOR is "one over r" -- standard coul */
enum {
    enbcoulNONE, enbcoulOOR, enbcoulRF, enbcoulTAB, enbcoulGB, enbcoulFEWALD, enbcoulNR
};

enum {
    egCOULSR, egLJSR, egBHAMSR, egCOULLR, egLJLR, egBHAMLR,
    egCOUL14, egLJ14, egGB, egNR
};

typedef struct {
    int   nener;      /* The number of energy group pairs     */
    real *ener[egNR]; /* Energy terms for each pair of groups */
} gmx_grppairener_t;

typedef struct {
    real              term[F_NRE];         /* The energies for all different interaction types */
    gmx_grppairener_t grpp;
    double            dvdl_lin[efptNR];    /* Contributions to dvdl with linear lam-dependence */
    double            dvdl_nonlin[efptNR]; /* Idem, but non-linear dependence                  */
    int               n_lambda;
    int               fep_state;           /*current fep state -- just for printing */
    double           *enerpart_lambda;     /* Partial energy for lambda and flambda[] */
    real              foreign_term[F_NRE]; /* alternate array for storing foreign lambda energies */
    gmx_grppairener_t foreign_grpp;        /* alternate array for storing foreign lambda energies */
} gmx_enerdata_t;
/* The idea is that dvdl terms with linear lambda dependence will be added
 * automatically to enerpart_lambda. Terms with non-linear lambda dependence
 * should explicitly determine the energies at foreign lambda points
 * when n_lambda > 0.
 */

typedef struct {
    int  cg_start;
    int  cg_end;
    int  cg_mod;
    int *cginfo;
} cginfo_mb_t;


/* Forward declaration of type for managing Ewald tables */
struct gmx_ewald_tab_t;

typedef struct f_thread_t f_thread_t;

typedef struct {
    interaction_const_t *ic;

    /* Domain Decomposition */
    gmx_bool bDomDec;

    /* PBC stuff */
    int                  ePBC;
    gmx_bool             bMolPBC;
    int                  rc_scaling;
    rvec                 posres_com;
    rvec                 posres_comB;

    const gmx_hw_info_t *hwinfo;
    const gmx_gpu_opt_t *gpu_opt;
    gmx_bool             use_simd_kernels;

    /* Interaction for calculated in kernels. In many cases this is similar to
     * the electrostatics settings in the inputrecord, but the difference is that
     * these variables always specify the actual interaction in the kernel - if
     * we are tabulating reaction-field the inputrec will say reaction-field, but
     * the kernel interaction will say cubic-spline-table. To be safe we also
     * have a kernel-specific setting for the modifiers - if the interaction is
     * tabulated we already included the inputrec modification there, so the kernel
     * modification setting will say 'none' in that case.
     */
    int nbkernel_elec_interaction;
    int nbkernel_vdw_interaction;
    int nbkernel_elec_modifier;
    int nbkernel_vdw_modifier;

    /* Use special N*N kernels? */
    gmx_bool bAllvsAll;
    /* Private work data */
    void    *AllvsAll_work;
    void    *AllvsAll_workgb;

    /* Cut-Off stuff.
     * Infinite cut-off's will be GMX_CUTOFF_INF (unlike in t_inputrec: 0).
     */
    real rlist, rlistlong;

    /* Dielectric constant resp. multiplication factor for charges */
    real zsquare, temp;
    real epsilon_r, epsilon_rf, epsfac;

    /* Constants for reaction fields */
    real kappa, k_rf, c_rf;

    /* Charge sum and dipole for topology A/B ([0]/[1]) for Ewald corrections */
    double qsum[2];
    double q2sum[2];
    double c6sum[2];
    rvec   mu_tot[2];

    /* Dispersion correction stuff */
    int  eDispCorr;

    /* The shift of the shift or user potentials */
    real enershiftsix;
    real enershifttwelve;
    /* Integrated differces for energy and virial with cut-off functions */
    real enerdiffsix;
    real enerdifftwelve;
    real virdiffsix;
    real virdifftwelve;
    /* Constant for long range dispersion correction (average dispersion)
     * for topology A/B ([0]/[1]) */
    real avcsix[2];
    /* Constant for long range repulsion term. Relative difference of about
     * 0.1 percent with 0.8 nm cutoffs. But hey, it's cheap anyway...
     */
    real avctwelve[2];

    /* Fudge factors */
    real fudgeQQ;

    /* Table stuff */
    gmx_bool     bcoultab;
    gmx_bool     bvdwtab;
    /* The normal tables are in the nblists struct(s) below */
    t_forcetable tab14; /* for 1-4 interactions only */

    /* PPPM & Shifting stuff */
    int   coulomb_modifier;
    real  rcoulomb_switch, rcoulomb;
    real *phi;

    /* VdW stuff */
    int    vdw_modifier;
    double reppow;
    real   rvdw_switch, rvdw;
    real   bham_b_max;

    /* Free energy */
    int      efep;
    real     sc_alphavdw;
    real     sc_alphacoul;
    int      sc_power;
    real     sc_r_power;
    real     sc_sigma6_def;
    real     sc_sigma6_min;

    /* NS Stuff */
    int  eeltype;
    int  vdwtype;
    int  cg0, hcg;
    /* solvent_opt contains the enum for the most common solvent
     * in the system, which will be optimized.
     * It can be set to esolNO to disable all water optimization */
    int          solvent_opt;
    int          nWatMol;
    gmx_bool     bGrid;
    gmx_bool     bExcl_IntraCGAll_InterCGNone;
    cginfo_mb_t *cginfo_mb;
    int         *cginfo;
    rvec        *cg_cm;
    int          cg_nalloc;
    rvec        *shift_vec;

    /* The neighborlists including tables */
    int                        nnblists;
    int                       *gid2nblists;
    t_nblists                 *nblists;

    int                        cutoff_scheme; /* group- or Verlet-style cutoff */
    gmx_bool                   bNonbonded;    /* true if nonbonded calculations are *not* turned off */
    struct nonbonded_verlet_t *nbv;

    /* The wall tables (if used) */
    int            nwall;
    t_forcetable **wall_tab;

    /* The number of charge groups participating in do_force_lowlevel */
    int ncg_force;
    /* The number of atoms participating in do_force_lowlevel */
    int natoms_force;
    /* The number of atoms participating in force and constraints */
    int natoms_force_constr;
    /* The allocation size of vectors of size natoms_force */
    int nalloc_force;

    /* Twin Range stuff, f_twin has size natoms_force */
    gmx_bool bTwinRange;
    int      nlr;
    rvec    *f_twin;
    /* Constraint virial correction for multiple time stepping */
    tensor   vir_twin_constr;

    /* Forces that should not enter into the virial summation:
     * PPPM/PME/Ewald/posres
     */
    gmx_bool bF_NoVirSum;
    int      f_novirsum_n;
    int      f_novirsum_nalloc;
    rvec    *f_novirsum_alloc;
    /* Pointer that points to f_novirsum_alloc when pressure is calcaluted,
     * points to the normal force vectors wen pressure is not requested.
     */
    rvec *f_novirsum;

    /* Long-range forces and virial for PPPM/PME/Ewald */
    struct gmx_pme_t *pmedata;
    int               ljpme_combination_rule;
    tensor            vir_el_recip;
    tensor            vir_lj_recip;

    /* PME/Ewald stuff */
    gmx_bool                bEwald;
    real                    ewaldcoeff_q;
    real                    ewaldcoeff_lj;
    struct gmx_ewald_tab_t *ewald_table;

    /* Virial Stuff */
    rvec *fshift;
    rvec  vir_diag_posres;
    dvec  vir_wall_z;

    /* Non bonded Parameter lists */
    int      ntype; /* Number of atom types */
    gmx_bool bBHAM;
    real    *nbfp;
    real    *ljpme_c6grid; /* C6-values used on grid in LJPME */

    /* Energy group pair flags */
    int *egp_flags;

    /* Shell molecular dynamics flexible constraints */
    real fc_stepsize;

    /* Generalized born implicit solvent */
    gmx_bool       bGB;
    /* Generalized born stuff */
    real           gb_epsilon_solvent;
    /* Table data for GB */
    t_forcetable   gbtab;
    /* VdW radius for each atomtype (dim is thus ntype) */
    real          *atype_radius;
    /* Effective radius (derived from effective volume) for each type */
    real          *atype_vol;
    /* Implicit solvent - surface tension for each atomtype */
    real          *atype_surftens;
    /* Implicit solvent - radius for GB calculation */
    real          *atype_gb_radius;
    /* Implicit solvent - overlap for HCT model */
    real          *atype_S_hct;
    /* Generalized born interaction data */
    gmx_genborn_t *born;

    /* Table scale for GB */
    real gbtabscale;
    /* Table range for GB */
    real gbtabr;
    /* GB neighborlists (the sr list will contain for each atom all other atoms
     * (for use in the SA calculation) and the lr list will contain
     * for each atom all atoms 1-4 or greater (for use in the GB calculation)
     */
    t_nblist gblist_sr;
    t_nblist gblist_lr;
    t_nblist gblist;

    /* Inverse square root of the Born radii for implicit solvent */
    real *invsqrta;
    /* Derivatives of the potential with respect to the Born radii */
    real *dvda;
    /* Derivatives of the Born radii with respect to coordinates */
    real *dadx;
    real *dadx_rawptr;
    int   nalloc_dadx; /* Allocated size of dadx */

    /* If > 0 signals Test Particle Insertion,
     * the value is the number of atoms of the molecule to insert
     * Only the energy difference due to the addition of the last molecule
     * should be calculated.
     */
    gmx_bool n_tpi;

    /* Neighbor searching stuff */
    gmx_ns_t ns;

    /* QMMM stuff */
    gmx_bool         bQMMM;
    t_QMMMrec       *qr;

    /* QM-MM neighborlists */
    t_nblist QMMMlist;

    /* Limit for printing large forces, negative is don't print */
    real print_force;

    /* coarse load balancing time measurement */
    double t_fnbf;
    double t_wait;
    int    timesteps;

    /* parameter needed for AdResS simulation */
    int             adress_type;
    gmx_bool        badress_tf_full_box;
    real            adress_const_wf;
    real            adress_ex_width;
    real            adress_hy_width;
    int             adress_icor;
    int             adress_site;
    rvec            adress_refs;
    int             n_adress_tf_grps;
    int           * adress_tf_table_index;
    int            *adress_group_explicit;
    t_forcetable *  atf_tabs;
    real            adress_ex_forcecap;
    gmx_bool        adress_do_hybridpairs;

    /* User determined parameters, copied from the inputrec */
    int  userint1;
    int  userint2;
    int  userint3;
    int  userint4;
    real userreal1;
    real userreal2;
    real userreal3;
    real userreal4;

    /* Thread local force and energy data */
    /* FIXME move to bonded_thread_data_t */
    int         nthreads;
    int         red_ashift;
    int         red_nblock;
    f_thread_t *f_t;

    /* Maximum thread count for uniform distribution of bondeds over threads */
    int   bonded_max_nthread_uniform;

    /* Exclusion load distribution over the threads */
    int  *excl_load;
} t_forcerec;

/* Important: Starting with Gromacs-4.6, the values of c6 and c12 in the nbfp array have
 * been scaled by 6.0 or 12.0 to save flops in the kernels. We have corrected this everywhere
 * in the code, but beware if you are using these macros externally.
 */
#define C6(nbfp, ntp, ai, aj)     (nbfp)[2*((ntp)*(ai)+(aj))]
#define C12(nbfp, ntp, ai, aj)    (nbfp)[2*((ntp)*(ai)+(aj))+1]
#define BHAMC(nbfp, ntp, ai, aj)  (nbfp)[3*((ntp)*(ai)+(aj))]
#define BHAMA(nbfp, ntp, ai, aj)  (nbfp)[3*((ntp)*(ai)+(aj))+1]
#define BHAMB(nbfp, ntp, ai, aj)  (nbfp)[3*((ntp)*(ai)+(aj))+2]

#ifdef __cplusplus
}
#endif
#endif
