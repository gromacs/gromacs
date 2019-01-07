/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017,2018,2019, by the GROMACS development team, led by
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
#ifndef GMX_MDTYPES_TYPES_FORCEREC_H
#define GMX_MDTYPES_TYPES_FORCEREC_H

#include <array>
#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/interaction_const.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

struct ForceProviders;

/* Abstract type for PME that is defined only in the routine that use them. */
struct gmx_ns_t;
struct gmx_pme_t;
struct nonbonded_verlet_t;
struct bonded_threading_t;
struct t_forcetable;
struct t_nblist;
struct t_nblists;
struct t_QMMMrec;

namespace gmx
{
class GpuBonded;
}

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

struct cginfo_mb_t
{
    int  cg_start;
    int  cg_end;
    int  cg_mod;
    int *cginfo;
};


/* Forward declaration of type for managing Ewald tables */
struct gmx_ewald_tab_t;

struct ewald_corr_thread_t;

struct t_forcerec { // NOLINT (clang-analyzer-optin.performance.Padding)
    struct interaction_const_t *ic;

    /* PBC stuff */
    int                         ePBC;
    //! Whether PBC must be considered for bonded interactions.
    gmx_bool                    bMolPBC;
    int                         rc_scaling;
    rvec                        posres_com;
    rvec                        posres_comB;

    gmx_bool                    use_simd_kernels;

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

    /* Cut-Off stuff.
     * Infinite cut-off's will be GMX_CUTOFF_INF (unlike in t_inputrec: 0).
     */
    real rlist;

    /* Parameters for generalized reaction field */
    real zsquare, temp;

    /* Charge sum and dipole for topology A/B ([0]/[1]) for Ewald corrections */
    double qsum[2];
    double q2sum[2];
    double c6sum[2];
    rvec   mu_tot[2];

    /* Dispersion correction stuff */
    int                  eDispCorr;
    int                  numAtomsForDispersionCorrection;
    struct t_forcetable *dispersionCorrectionTable;

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
    gmx_bool             bcoultab;
    gmx_bool             bvdwtab;
    /* The normal tables are in the nblists struct(s) below */

    struct t_forcetable *pairsTable; /* for 1-4 interactions, [pairs] and [pairs_nb] */

    /* Free energy */
    int      efep;
    real     sc_alphavdw;
    real     sc_alphacoul;
    int      sc_power;
    real     sc_r_power;
    real     sc_sigma6_def;
    real     sc_sigma6_min;

    /* NS Stuff */
    int  cg0, hcg;
    /* solvent_opt contains the enum for the most common solvent
     * in the system, which will be optimized.
     * It can be set to esolNO to disable all water optimization */
    int                 solvent_opt;
    int                 nWatMol;
    gmx_bool            bGrid;
    gmx_bool            bExcl_IntraCGAll_InterCGNone;
    struct cginfo_mb_t *cginfo_mb;
    int                *cginfo;
    rvec               *cg_cm;
    int                 cg_nalloc;
    rvec               *shift_vec;

    /* The neighborlists including tables */
    int                        nnblists;
    int                       *gid2nblists;
    struct t_nblists          *nblists;

    int                        cutoff_scheme; /* group- or Verlet-style cutoff */
    gmx_bool                   bNonbonded;    /* true if nonbonded calculations are *not* turned off */
    struct nonbonded_verlet_t *nbv;

    /* The wall tables (if used) */
    int                    nwall;
    struct t_forcetable ***wall_tab;

    /* The number of charge groups participating in do_force_lowlevel */
    int ncg_force;
    /* The number of atoms participating in do_force_lowlevel */
    int natoms_force;
    /* The number of atoms participating in force and constraints */
    int natoms_force_constr;
    /* The allocation size of vectors of size natoms_force */
    int nalloc_force;

    /* Forces that should not enter into the coord x force virial summation:
     * PPPM/PME/Ewald/posres/ForceProviders
     */
    /* True when we have contributions that are directly added to the virial */
    gmx_bool                 haveDirectVirialContributions;
    /* TODO: Replace the pointer by an object once we got rid of C */
    std::vector<gmx::RVec>  *forceBufferForDirectVirialContributions;

    /* Data for PPPM/PME/Ewald */
    struct gmx_pme_t *pmedata;
    int               ljpme_combination_rule;

    /* PME/Ewald stuff */
    struct gmx_ewald_tab_t *ewald_table;

    /* Shift force array for computing the virial */
    rvec *fshift;

    /* Non bonded Parameter lists */
    int      ntype; /* Number of atom types */
    gmx_bool bBHAM;
    real    *nbfp;
    real    *ljpme_c6grid; /* C6-values used on grid in LJPME */

    /* Energy group pair flags */
    int *egp_flags;

    /* Shell molecular dynamics flexible constraints */
    real fc_stepsize;

    /* If > 0 signals Test Particle Insertion,
     * the value is the number of atoms of the molecule to insert
     * Only the energy difference due to the addition of the last molecule
     * should be calculated.
     */
    int n_tpi;

    /* Neighbor searching stuff */
    struct gmx_ns_t *ns;

    /* QMMM stuff */
    gmx_bool          bQMMM;
    struct t_QMMMrec *qr;

    /* QM-MM neighborlists */
    struct t_nblist        *QMMMlist;

    /* Limit for printing large forces, negative is don't print */
    real print_force;

    /* coarse load balancing time measurement */
    double t_fnbf;
    double t_wait;
    int    timesteps;

    /* User determined parameters, copied from the inputrec */
    int  userint1;
    int  userint2;
    int  userint3;
    int  userint4;
    real userreal1;
    real userreal2;
    real userreal3;
    real userreal4;

    /* Pointer to struct for managing threading of bonded force calculation */
    struct bonded_threading_t *bondedThreading;

    /* TODO: Replace the pointer by an object once we got rid of C */
    gmx::GpuBonded *gpuBonded;

    /* Ewald correction thread local virial and energy data */
    int                         nthread_ewc;
    struct ewald_corr_thread_t *ewc_t;

    struct ForceProviders      *forceProviders;
};

/* Important: Starting with Gromacs-4.6, the values of c6 and c12 in the nbfp array have
 * been scaled by 6.0 or 12.0 to save flops in the kernels. We have corrected this everywhere
 * in the code, but beware if you are using these macros externally.
 */
#define C6(nbfp, ntp, ai, aj)     (nbfp)[2*((ntp)*(ai)+(aj))]
#define C12(nbfp, ntp, ai, aj)    (nbfp)[2*((ntp)*(ai)+(aj))+1]
#define BHAMC(nbfp, ntp, ai, aj)  (nbfp)[3*((ntp)*(ai)+(aj))]
#define BHAMA(nbfp, ntp, ai, aj)  (nbfp)[3*((ntp)*(ai)+(aj))+1]
#define BHAMB(nbfp, ntp, ai, aj)  (nbfp)[3*((ntp)*(ai)+(aj))+2]

#endif
