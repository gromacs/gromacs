/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2012- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */

/* We could use nbat->xstride and nbat->fstride, but macros might be faster */
#define X_STRIDE 3
#define F_STRIDE 3
/* Local i-atom buffer strides */
#define XI_STRIDE 3
#define FI_STRIDE 3


/* All functionality defines are set here, except for:
 * CALC_ENERGIES, ENERGY_GROUPS which are defined before.
 * CHECK_EXCLS, which is set just before including the inner loop contents.
 */

/* We always calculate shift forces, because it's cheap anyhow */
#define CALC_SHIFTFORCES

#ifdef CALC_COUL_RF
#    define NBK_FUNC_NAME2(ui, uj, ljt, feg) nbnxn_kernel_##ui##x##uj##_ElecRF##ljt##feg##_ref
#endif
#ifdef CALC_COUL_TAB
#    ifndef VDW_CUTOFF_CHECK
#        define NBK_FUNC_NAME2(ui, uj, ljt, feg) \
            nbnxn_kernel_##ui##x##uj##_ElecQSTab##ljt##feg##_ref
#    else
#        define NBK_FUNC_NAME2(ui, uj, ljt, feg) \
            nbnxn_kernel_##ui##x##uj##_ElecQSTabTwinCut##ljt##feg##_ref
#    endif
#endif

#ifdef CALC_COUL_NONE
#    define NBK_FUNC_NAME2(ui, uj, ljt, feg) nbnxn_kernel_##ui##x##uj##_ElecNone##ljt##feg##_ref
#    define HAVE_ELECTROSTATICS 0
#else
#    define HAVE_ELECTROSTATICS 1
#endif

#if defined LJ_CUT && !defined LJ_EWALD
#    define NBK_FUNC_NAME(ui, uj, feg) NBK_FUNC_NAME2(ui, uj, _VdwLJ, feg)
#elif defined LJ_FORCE_SWITCH
#    define NBK_FUNC_NAME(ui, uj, feg) NBK_FUNC_NAME2(ui, uj, _VdwLJFsw, feg)
#elif defined LJ_POT_SWITCH
#    define NBK_FUNC_NAME(ui, uj, feg) NBK_FUNC_NAME2(ui, uj, _VdwLJPsw, feg)
#elif defined LJ_EWALD
#    ifdef LJ_EWALD_COMB_GEOM
#        define NBK_FUNC_NAME(ui, uj, feg) NBK_FUNC_NAME2(ui, uj, _VdwLJEwCombGeom, feg)
#    else
#        define NBK_FUNC_NAME(ui, uj, feg) NBK_FUNC_NAME2(ui, uj, _VdwLJEwCombLB, feg)
#    endif
#else
#    error "No VdW type defined"
#endif

void
#ifndef CALC_ENERGIES
        NBK_FUNC_NAME(UNROLLI, UNROLLJ, _F) // NOLINT(misc-definitions-in-headers)
#else
#    ifndef ENERGY_GROUPS
        NBK_FUNC_NAME(UNROLLI, UNROLLJ, _VF) // NOLINT(misc-definitions-in-headers)
#    else
        NBK_FUNC_NAME(UNROLLI, UNROLLJ, _VgrpF) // NOLINT(misc-definitions-in-headers)
#    endif
#endif
#undef NBK_FUNC_NAME
#undef NBK_FUNC_NAME2
        (const NbnxnPairlistCpu&    pairlist,
         const nbnxn_atomdata_t&    nbat,
         const interaction_const_t& ic,
         const rvec*                shift_vec,
         nbnxn_atomdata_output_t*   out)
{
    GMX_RELEASE_ASSERT(UNROLLI == pairlist.na_ci && UNROLLJ == pairlist.na_cj,
                       "Kernel and list cluster sizes should match");

    static_assert(GMX_USE_EXT_FMM || HAVE_ELECTROSTATICS,
                  "Reference kernels that do not compute Coulomb interactions are supported only "
                  "with an FMM build configuration");

    /* Unpack pointers for output */
    real* f = out->f.data();
#ifdef CALC_SHIFTFORCES
    real* fshift = out->fshift.data();
#endif
#ifdef CALC_ENERGIES
    real* Vvdw = out->Vvdw.data();
#    if HAVE_ELECTROSTATICS
    real* Vc = out->Vc.data();
#    endif
#endif

    real xi[UNROLLI * XI_STRIDE];
    real fi[UNROLLI * FI_STRIDE];
#if HAVE_ELECTROSTATICS
    real qi[UNROLLI];
#endif

#ifdef COUNT_PAIRS
    int npair = 0;
#endif

#ifdef LJ_POT_SWITCH
    const real swV3 = ic.vdw.switchConstants.c3;
    const real swV4 = ic.vdw.switchConstants.c4;
    const real swV5 = ic.vdw.switchConstants.c5;
    const real swF2 = 3 * ic.vdw.switchConstants.c3;
    const real swF3 = 4 * ic.vdw.switchConstants.c4;
    const real swF4 = 5 * ic.vdw.switchConstants.c5;
#endif

    const nbnxn_atomdata_t::Params& nbatParams = nbat.params();

#ifdef LJ_EWALD
    const real lje_coeff2   = gmx::square(ic.vdw.ewaldCoeff);
    const real lje_coeff6_6 = lje_coeff2 * lje_coeff2 * lje_coeff2 / 6.0;
#    ifdef CALC_ENERGIES
    const real lje_vc = ic.vdw.ewaldShift;
#    endif

    const real* ljc = nbatParams.nbfp_comb.data();
#endif

#ifdef CALC_COUL_RF
    const real k_rf2 = 2 * ic.coulomb.reactionFieldCoefficient;
#    ifdef CALC_ENERGIES
    const real reactionFieldCoefficient = ic.coulomb.reactionFieldCoefficient;
    const real reactionFieldShift       = ic.coulomb.reactionFieldShift;
#    endif
#endif
#ifdef CALC_COUL_TAB
    const real tab_coul_scale = ic.coulombEwaldTables->scale;
#    ifdef CALC_ENERGIES
    const real halfsp = 0.5 / tab_coul_scale;
#    endif

#    if !GMX_DOUBLE
    const real* tab_coul_FDV0 = ic.coulombEwaldTables->tableFDV0.data();
#    else
    const real* tab_coul_F = ic.coulombEwaldTables->tableF.data();
#        ifdef CALC_ENERGIES
    const real* tab_coul_V = ic.coulombEwaldTables->tableV.data();
#        endif
#    endif
#endif

    const real rcut2 = gmx::square(ic.coulomb.cutoff);
#ifdef VDW_CUTOFF_CHECK
    const real rvdw2 = gmx::square(ic.vdw.cutoff);
#endif

    const int   ntype2 = nbatParams.numTypes * 2;
    const real* nbfp   = nbatParams.nbfp.data();
#if HAVE_ELECTROSTATICS
    const real* q     = nbatParams.q.data();
    const real  facel = ic.coulomb.epsfac;
#endif
    const int* type = nbatParams.type.data();

    const real* shiftvec = shift_vec[0];
    const real* x        = nbat.x().data();

    const nbnxn_cj_t* l_cj = pairlist.cj.list_.data();

    for (const nbnxn_ci_t& ciEntry : pairlist.ci)
    {
        const int ish = (ciEntry.shift & NBNXN_CI_SHIFT);
        /* x, f and fshift are assumed to be stored with stride 3 */
        const int ishf   = ish * DIM;
        const int cjind0 = ciEntry.cj_ind_start;
        const int cjind1 = ciEntry.cj_ind_end;
        /* Currently only works super-cells equal to sub-cells */
        const int ci = ciEntry.ci;
#if HAVE_ELECTROSTATICS || defined(LJ_EWALD)
        const int ci_sh = (ish == gmx::c_centralShiftIndex ? ci : -1);
#endif
/* We have 5 LJ/C combinations, but use only three inner loops,
 * as the other combinations are unlikely and/or not much faster:
 * inner half-LJ + C for half-LJ + C / no-LJ + C
 * inner LJ + C      for full-LJ + C
 * inner LJ          for full-LJ + no-C / half-LJ + no-C
 */
#if HAVE_ELECTROSTATICS
        const bool do_LJ   = ((ciEntry.shift & NBNXN_CI_DO_LJ(0)) != 0);
        const bool do_coul = ((ciEntry.shift & NBNXN_CI_DO_COUL(0)) != 0);
        const bool half_LJ = (((ciEntry.shift & NBNXN_CI_HALF_LJ(0)) != 0) || !do_LJ) && do_coul;
#endif
#ifdef CALC_ENERGIES

#    ifdef LJ_EWALD
        const bool do_self = true;
#    elif HAVE_ELECTROSTATICS
        const bool do_self = do_coul;
#    endif


#    ifndef ENERGY_GROUPS
        real Vvdw_ci = 0;
#        if HAVE_ELECTROSTATICS
        real Vc_ci = 0;
#        endif
#    else
        int egp_sh_i[UNROLLI];
        for (int i = 0; i < UNROLLI; i++)
        {
            egp_sh_i[i] = nbatParams.energyGroupsPerCluster->getEnergyGroup(ci, i) * nbatParams.numEnergyGroups;
        }
#    endif
#endif

        for (int i = 0; i < UNROLLI; i++)
        {
            for (int d = 0; d < DIM; d++)
            {
                xi[i * XI_STRIDE + d] = x[(ci * UNROLLI + i) * X_STRIDE + d] + shiftvec[ishf + d];
                fi[i * FI_STRIDE + d] = 0;
            }

#if HAVE_ELECTROSTATICS
            qi[i] = facel * q[ci * UNROLLI + i];
#endif
        }

#if defined CALC_ENERGIES && (HAVE_ELECTROSTATICS || defined LJ_EWALD)
        if (do_self)
        {
#    ifdef CALC_COUL_RF
            const real Vc_sub_self = 0.5 * reactionFieldShift;
#    endif
#    ifdef CALC_COUL_TAB
#        if GMX_DOUBLE
            const real Vc_sub_self = 0.5 * tab_coul_V[0];
#        else
            const real Vc_sub_self = 0.5 * tab_coul_FDV0[2];
#        endif
#    endif

            if (l_cj[ciEntry.cj_ind_start].cj == ci_sh)
            {
                for (int i = 0; i < UNROLLI; i++)
                {
#    ifdef ENERGY_GROUPS
                    const int egp_ind =
                            egp_sh_i[i] + nbatParams.energyGroupsPerCluster->getEnergyGroup(ci, i);
#    else
                    const int egp_ind = 0;
#    endif

#    if HAVE_ELECTROSTATICS
                    /* Coulomb self interaction */
                    Vc[egp_ind] -= qi[i] * q[ci * UNROLLI + i] * Vc_sub_self;
#    endif

#    ifdef LJ_EWALD
                    /* LJ Ewald self interaction */
                    Vvdw[egp_ind] +=
                            0.5
                            * nbatParams.nbfp[nbatParams.type[ci * UNROLLI + i] * (nbatParams.numTypes + 1) * 2]
                            / 6 * lje_coeff6_6;
#    endif
                }
            }
        }
#endif /* CALC_ENERGIES */

        // Without loop vectorization we first loop over all pairs with exclusions and then loop
        // over the remaining pairs without exclusion without checking for exclusions
        int cjind = cjind0;
#if VECTORIZE_JLOOP && defined __clang__
#    pragma clang loop vectorize(assume_safety)
#endif
        for (; cjind < cjind1 && (VECTORIZE_JLOOP || pairlist.cj.excl(cjind) != 0xffff); cjind++)
        {
#define CHECK_EXCLS
#if HAVE_ELECTROSTATICS
            if (half_LJ)
            {
#    define CALC_COULOMB
#    define HALF_LJ
#    include "kernel_ref_inner.h"
#    undef HALF_LJ
#    undef CALC_COULOMB
            }
            else if (do_coul)
            {
#    define CALC_COULOMB
#    include "kernel_ref_inner.h"
#    undef CALC_COULOMB
            }
            else
            {
#    include "kernel_ref_inner.h"
            }
#else
#    include "kernel_ref_inner.h"
#endif
#undef CHECK_EXCLS
        }

#if !VECTORIZE_JLOOP
        // Second part of the j-loop, does not check for exclusions
        for (; (cjind < cjind1); cjind++)
        {
#    if HAVE_ELECTROSTATICS
            if (half_LJ)
            {
#        define CALC_COULOMB
#        define HALF_LJ
#        include "kernel_ref_inner.h"
#        undef HALF_LJ
#        undef CALC_COULOMB
            }
            else if (do_coul)
            {
#        define CALC_COULOMB
#        include "kernel_ref_inner.h"
#        undef CALC_COULOMB
            }
            else
            {
#        include "kernel_ref_inner.h"
            }
#    else
#        include "kernel_ref_inner.h"
#    endif
        }
#endif // !VECTORIZE_JLOOP

        /* Add accumulated i-forces to the force array */
        for (int i = 0; i < UNROLLI; i++)
        {
            for (int d = 0; d < DIM; d++)
            {
                f[(ci * UNROLLI + i) * F_STRIDE + d] += fi[i * FI_STRIDE + d];
            }
        }
#ifdef CALC_SHIFTFORCES
        if (fshift != nullptr)
        {
            /* Add i forces to shifted force list */
            for (int i = 0; i < UNROLLI; i++)
            {
                for (int d = 0; d < DIM; d++)
                {
                    fshift[ishf + d] += fi[i * FI_STRIDE + d];
                }
            }
        }
#endif

#ifdef CALC_ENERGIES
#    ifndef ENERGY_GROUPS
        *Vvdw += Vvdw_ci;

#        if HAVE_ELECTROSTATICS
        *Vc += Vc_ci;
#        endif
#    endif
#endif
    }

#ifdef COUNT_PAIRS
    printf("atom pairs %d\n", npair);
#endif
}

#undef CALC_SHIFTFORCES

#undef X_STRIDE
#undef F_STRIDE
#undef XI_STRIDE
#undef FI_STRIDE
#undef HAVE_ELECTROSTATICS
