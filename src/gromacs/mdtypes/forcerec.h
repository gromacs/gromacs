/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017 by the GROMACS development team.
 * Copyright (c) 2018,2019,2020, by the GROMACS development team, led by
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
#include <memory>
#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

/* Abstract type for PME that is defined only in the routine that use them. */
struct gmx_pme_t;
struct nonbonded_verlet_t;
struct bonded_threading_t;
class DeviceContext;
class DispersionCorrection;
struct t_forcetable;
struct t_QMMMrec;

namespace gmx
{
class DeviceStreamManager;
class GpuBonded;
class ForceProviders;
class StatePropagatorDataGpu;
class PmePpCommGpu;
class WholeMoleculeTransform;
} // namespace gmx

/* macros for the cginfo data in forcerec
 *
 * Since the tpx format support max 256 energy groups, we do the same here.
 * Note that we thus have bits 8-14 still unused.
 *
 * The maximum cg size in cginfo is 63
 * because we only have space for 6 bits in cginfo,
 * this cg size entry is actually only read with domain decomposition.
 */
#define SET_CGINFO_GID(cgi, gid) (cgi) = (((cgi) & ~255) | (gid))
#define GET_CGINFO_GID(cgi) ((cgi)&255)
#define SET_CGINFO_FEP(cgi) (cgi) = ((cgi) | (1 << 15))
#define GET_CGINFO_FEP(cgi) ((cgi) & (1 << 15))
#define SET_CGINFO_EXCL_INTER(cgi) (cgi) = ((cgi) | (1 << 17))
#define GET_CGINFO_EXCL_INTER(cgi) ((cgi) & (1 << 17))
#define SET_CGINFO_CONSTR(cgi) (cgi) = ((cgi) | (1 << 20))
#define GET_CGINFO_CONSTR(cgi) ((cgi) & (1 << 20))
#define SET_CGINFO_SETTLE(cgi) (cgi) = ((cgi) | (1 << 21))
#define GET_CGINFO_SETTLE(cgi) ((cgi) & (1 << 21))
/* This bit is only used with bBondComm in the domain decomposition */
#define SET_CGINFO_BOND_INTER(cgi) (cgi) = ((cgi) | (1 << 22))
#define GET_CGINFO_BOND_INTER(cgi) ((cgi) & (1 << 22))
#define SET_CGINFO_HAS_VDW(cgi) (cgi) = ((cgi) | (1 << 23))
#define GET_CGINFO_HAS_VDW(cgi) ((cgi) & (1 << 23))
#define SET_CGINFO_HAS_Q(cgi) (cgi) = ((cgi) | (1 << 24))
#define GET_CGINFO_HAS_Q(cgi) ((cgi) & (1 << 24))


/* Value to be used in mdrun for an infinite cut-off.
 * Since we need to compare with the cut-off squared,
 * this value should be slighlty smaller than sqrt(GMX_FLOAT_MAX).
 */
#define GMX_CUTOFF_INF 1E+18

/* enums for the neighborlist type */
enum
{
    enbvdwNONE,
    enbvdwLJ,
    enbvdwBHAM,
    enbvdwTAB,
    enbvdwNR
};

struct cginfo_mb_t
{
    int              cg_start = 0;
    int              cg_end   = 0;
    int              cg_mod   = 0;
    std::vector<int> cginfo;
};


/* Forward declaration of type for managing Ewald tables */
struct gmx_ewald_tab_t;

struct ewald_corr_thread_t;

struct t_forcerec
{ // NOLINT (clang-analyzer-optin.performance.Padding)
    // Declare an explicit constructor and destructor, so they can be
    // implemented in a single source file, so that not every source
    // file that includes this one needs to understand how to find the
    // destructors of the objects pointed to by unique_ptr members.
    t_forcerec();
    ~t_forcerec();

    struct interaction_const_t* ic = nullptr;

    /* PBC stuff */
    PbcType pbcType = PbcType::Xyz;
    //! Tells whether atoms inside a molecule can be in different periodic images,
    //  i.e. whether we need to take into account PBC when computing distances inside molecules.
    //  This determines whether PBC must be considered for e.g. bonded interactions.
    gmx_bool bMolPBC     = FALSE;
    int      rc_scaling  = 0;
    rvec     posres_com  = { 0 };
    rvec     posres_comB = { 0 };

    gmx_bool use_simd_kernels = FALSE;

    /* Interaction for calculated in kernels. In many cases this is similar to
     * the electrostatics settings in the inputrecord, but the difference is that
     * these variables always specify the actual interaction in the kernel - if
     * we are tabulating reaction-field the inputrec will say reaction-field, but
     * the kernel interaction will say cubic-spline-table. To be safe we also
     * have a kernel-specific setting for the modifiers - if the interaction is
     * tabulated we already included the inputrec modification there, so the kernel
     * modification setting will say 'none' in that case.
     */
    int nbkernel_elec_interaction = 0;
    int nbkernel_vdw_interaction  = 0;
    int nbkernel_elec_modifier    = 0;
    int nbkernel_vdw_modifier     = 0;

    /* Cut-Off stuff.
     * Infinite cut-off's will be GMX_CUTOFF_INF (unlike in t_inputrec: 0).
     */
    real rlist = 0;

    /* Charge sum for topology A/B ([0]/[1]) for Ewald corrections */
    double qsum[2]  = { 0 };
    double q2sum[2] = { 0 };
    double c6sum[2] = { 0 };

    /* Dispersion correction stuff */
    std::unique_ptr<DispersionCorrection> dispersionCorrection;

    /* Fudge factors */
    real fudgeQQ = 0;

    /* Table stuff */
    gmx_bool bcoultab = FALSE;
    gmx_bool bvdwtab  = FALSE;

    t_forcetable* pairsTable = nullptr; /* for 1-4 interactions, [pairs] and [pairs_nb] */

    /* Free energy */
    int  efep          = 0;
    real sc_alphavdw   = 0;
    real sc_alphacoul  = 0;
    int  sc_power      = 0;
    real sc_r_power    = 0;
    real sc_sigma6_def = 0;
    real sc_sigma6_min = 0;

    /* Information about atom properties for the molecule blocks in the system */
    std::vector<cginfo_mb_t> cginfo_mb;
    /* Information about atom properties for local and non-local atoms */
    std::vector<int> cginfo;

    rvec* shift_vec = nullptr;

    std::unique_ptr<gmx::WholeMoleculeTransform> wholeMoleculeTransform;

    int      cutoff_scheme = 0;     /* group- or Verlet-style cutoff */
    gmx_bool bNonbonded    = FALSE; /* true if nonbonded calculations are *not* turned off */

    /* The Nbnxm Verlet non-bonded machinery */
    std::unique_ptr<nonbonded_verlet_t> nbv;

    /* The wall tables (if used) */
    int             nwall    = 0;
    t_forcetable*** wall_tab = nullptr;

    /* The number of atoms participating in do_force_lowlevel */
    int natoms_force = 0;
    /* The number of atoms participating in force calculation and constraints */
    int natoms_force_constr = 0;
    /* Forces that should not enter into the coord x force virial summation:
     * PPPM/PME/Ewald/posres/ForceProviders
     */
    /* True when we have contributions that are directly added to the virial */
    bool haveDirectVirialContributions = false;
    /* Force buffer for force computation with direct virial contributions */
    std::vector<gmx::RVec> forceBufferForDirectVirialContributions;

    /* Data for PPPM/PME/Ewald */
    struct gmx_pme_t* pmedata                = nullptr;
    int               ljpme_combination_rule = 0;

    /* PME/Ewald stuff */
    struct gmx_ewald_tab_t* ewald_table = nullptr;

    /* Shift force array for computing the virial, size SHIFTS */
    std::vector<gmx::RVec> shiftForces;

    /* Non bonded Parameter lists */
    int               ntype = 0; /* Number of atom types */
    gmx_bool          bBHAM = FALSE;
    std::vector<real> nbfp;
    real*             ljpme_c6grid = nullptr; /* C6-values used on grid in LJPME */

    /* Energy group pair flags */
    int* egp_flags = nullptr;

    /* Shell molecular dynamics flexible constraints */
    real fc_stepsize = 0;

    /* If > 0 signals Test Particle Insertion,
     * the value is the number of atoms of the molecule to insert
     * Only the energy difference due to the addition of the last molecule
     * should be calculated.
     */
    int n_tpi = 0;

    /* Limit for printing large forces, negative is don't print */
    real print_force = 0;

    /* User determined parameters, copied from the inputrec */
    int  userint1  = 0;
    int  userint2  = 0;
    int  userint3  = 0;
    int  userint4  = 0;
    real userreal1 = 0;
    real userreal2 = 0;
    real userreal3 = 0;
    real userreal4 = 0;

    /* Pointer to struct for managing threading of bonded force calculation */
    struct bonded_threading_t* bondedThreading = nullptr;

    /* TODO: Replace the pointer by an object once we got rid of C */
    gmx::GpuBonded* gpuBonded = nullptr;

    /* Ewald correction thread local virial and energy data */
    int                         nthread_ewc = 0;
    struct ewald_corr_thread_t* ewc_t       = nullptr;

    gmx::ForceProviders* forceProviders = nullptr;

    // The stateGpu object is created in runner, forcerec just keeps the copy of the pointer.
    // TODO: This is not supposed to be here. StatePropagatorDataGpu should be a part of
    //       general StatePropagatorData object that is passed around
    gmx::StatePropagatorDataGpu* stateGpu = nullptr;
    // TODO: Should not be here. This is here only to pass the pointer around.
    gmx::DeviceStreamManager* deviceStreamManager = nullptr;

    //! GPU device context
    DeviceContext* deviceContext = nullptr;

    /* For PME-PP GPU communication */
    std::unique_ptr<gmx::PmePpCommGpu> pmePpCommGpu;
};

/* Important: Starting with Gromacs-4.6, the values of c6 and c12 in the nbfp array have
 * been scaled by 6.0 or 12.0 to save flops in the kernels. We have corrected this everywhere
 * in the code, but beware if you are using these macros externally.
 */
#define C6(nbfp, ntp, ai, aj) (nbfp)[2 * ((ntp) * (ai) + (aj))]
#define C12(nbfp, ntp, ai, aj) (nbfp)[2 * ((ntp) * (ai) + (aj)) + 1]
#define BHAMC(nbfp, ntp, ai, aj) (nbfp)[3 * ((ntp) * (ai) + (aj))]
#define BHAMA(nbfp, ntp, ai, aj) (nbfp)[3 * ((ntp) * (ai) + (aj)) + 1]
#define BHAMB(nbfp, ntp, ai, aj) (nbfp)[3 * ((ntp) * (ai) + (aj)) + 2]

#endif
