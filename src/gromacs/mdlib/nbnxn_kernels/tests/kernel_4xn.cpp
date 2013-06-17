/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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
/*! \internal \file
 * \brief
 * Tests for 4xn kernels
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_nbnxn
 */

#include "../simd_4xn/nbnxn_kernel_simd_4xn.h"

#include "gromacs/legacyheaders/types/nbnxn_pairlist.h"
#include "gromacs/legacyheaders/types/nb_verlet.h"
#include "gromacs/simd/macros.h"

// TODO relax this to test other simd widths
#if defined GMX_NBNXN_SIMD_4XN && (4 == GMX_SIMD_WIDTH_HERE)

#include <gmock/gmock.h>
#include "testutils/testexceptions.h"
#include "testutils/refdata.h"

#include "typedefs.h"
#include "smalloc.h"
#include "gromacs/utility/gmxassert.h"

#include "gromacs/legacyheaders/types/simple.h"
#include "gromacs/legacyheaders/types/ishift.h"
#include "gromacs/legacyheaders/coulomb.h"
#include "gromacs/legacyheaders/string2.h"
#include "gromacs/legacyheaders/names.h"
#include "gromacs/legacyheaders/force.h"
#include "gromacs/legacyheaders/pbc.h"
#include "gromacs/legacyheaders/gmx_omp_nthreads.h"
#include "gromacs/mdlib/nbnxn_atomdata.h"
#include "gromacs/mdlib/nbnxn_search.h"
#include "gromacs/simd/tests/utils.h"

namespace KernelTests
{

class Test4xnKernelObject
{
    public:
        Test4xnKernelObject(t_forcerec *fr) :
            nbl_list(build_nbl_list()),
            nbat(build_nbat(fr->ntype, fr->nbfp)),
            ic(build_ic(fr)),
            ewald_excl(ewaldexclTable),
            shift_vec(NULL),
            force_flags(0),
            clearF(false),
            fshift(NULL),
            Vc(NULL),
            Vvdw(NULL),
            forceArraySize(0)
        {
            // TODO fix this?
            //            box({{0,0,0}, {0,0,0}, {0,0,0}}),
        };

        ~Test4xnKernelObject()
        {
            sfree(nbl_list->nbl[0]->ci);
            sfree(nbl_list->nbl[0]->cj);
            sfree(nbl_list->nbl[0]);
            sfree(nbl_list->nbl);
            sfree(nbl_list);
            sfree(nbat->x);
            sfree(nbat->q);
            sfree(nbat->type);
            sfree(nbat->out[0].Vc);
            sfree(nbat->out[0].Vvdw);
            sfree(nbat->out[0].f);
            sfree(nbat->out[0].fshift);
            sfree(nbat->out);
            sfree(nbat->simd_4xn_diagonal_j_minus_i);
            sfree(nbat->simd_2xnn_diagonal_j_minus_i);
            sfree(nbat->simd_exclusion_filter1);
            sfree(nbat->simd_exclusion_filter2);
#ifdef GMX_CPU_ACCELERATION_IBM_QPX
            sfree(nbat->simd_interaction_array);
#endif
            sfree(nbat->shift_vec);
            sfree(nbat->nbfp); // TODO why this as well as fr->nbfp?
            sfree(nbat->nbfp_comb);
            sfree(nbat->lj_comb);
            sfree(nbat);
            sfree(ic->tabq_coul_FDV0);
            sfree(ic->tabq_coul_F);
            sfree(ic->tabq_coul_V);
            sfree(ic);
            sfree(shift_vec);
            sfree(fshift);
        };

        interaction_const_t *build_ic(t_forcerec const *fr)
        {
            interaction_const_t *ic;
            real                 rtab = 1.5;
            init_interaction_const(NULL, &ic, fr, rtab);
            return ic;
        }

        nbnxn_atomdata_t *build_nbat(int ntype, const real *nbfp)
        {
            nbnxn_atomdata_t *nbat;
            snew(nbat, 1);
            int               nb_kernel_type = nbnxnk4xN_SIMD_4xN;
            int               n_energygroups = 0;
            int               nout           = 1;
            nbnxn_alloc_t    *alloc          = nbnxn_alloc_aligned;
            nbnxn_free_t     *free           = nbnxn_free_aligned;
            nbnxn_atomdata_init(NULL, nbat, nb_kernel_type, ntype,
                                nbfp, n_energygroups, nout, alloc, free);
            return nbat;
        }

        /*

           nbnxn_pairlist_t *build_nbnxn_pairlist()
           {
            nbnxn_pairlist_t *nbl;
            snew(nbl, 1);
            gmx_bool bSimple = true;
            nbnxn_init_pairlist(nbl, bSimple, NULL, NULL);
            return nbl;
           }
         */

        nbnxn_pairlist_set_t *build_nbl_list()
        {
            nbnxn_pairlist_set_t *nbl_list;
            snew(nbl_list, 1);
            gmx_bool              bCombined = true;
            gmx_bool              bSimple   = true;
            nbnxn_init_pairlist_set(nbl_list, bCombined, bSimple, nbnxn_alloc_aligned, nbnxn_free_aligned);
            /*
               nbl_list->nnbl = 1;
               snew(nbl_list->nbl, nbl_list->nnbl);
               nbl_list->nbl[0] = build_nbnxn_pairlist();
               nbl_list->bCombined = true;
               nbl_list->bSimple = true;
               nbl_list->natpair_ljq = 0;
               nbl_list->natpair_lj = 0;
               nbl_list->natpair_q = 1;
             */
            return nbl_list;
        }

        void set_nbat()
        {
            int i, j;
            GMX_ASSERT(nbatX4 == nbat->XFormat, "XFormat must be nbatX4");
            /* Thus XXXXYYYYZZZZ layout */

            // four atoms for each of one i and j cluster
            nbat->natoms = 8;
            int xsize = nbat->natoms * nbat->xstride;
            nbnxn_alloc_aligned((void **)&nbat->x, xsize * sizeof(nbat->x));
            for (i = 0, j = 0; j != 4 && i < xsize; ++i, ++j)
            {
                nbat->x[i] = j;
            }
            for (; i < xsize/2; ++i, ++j)
            {
                nbat->x[i] = j - 0.2;
            }
            for (j = 0; j != 4 && i < xsize; ++i, ++j)
            {
                nbat->x[i] = real(j) + 0.1;
            }
            for (; i < xsize; ++i, ++j)
            {
                nbat->x[i] = j - 0.15;
            }

            // one charged atom in i and j
            int qsize = nbat->natoms;
            nbnxn_alloc_aligned((void**)&nbat->q, qsize * sizeof(nbat->q));
            for (i = 0, j = 0; j != 1 && i < qsize; ++i, ++j)
            {
                nbat->q[i] = -1;
            }
            for (; i < qsize/2; ++i, ++j)
            {
                nbat->q[i] = i * -0.08;
            }
            for (j = 0; j != 1 && i < qsize; ++i, ++j)
            {
                nbat->q[i] = 0.5;
            }
            for (; i < qsize; ++i, ++j)
            {
                nbat->q[i] = i * 0.13;
            }

            snew(nbat->type, nbat->natoms);
            for (int i = 0; i < nbat->natoms; ++i)
            {
                nbat->type[i] = 0;
            }

            nbat->alloc((void**) &nbat->lj_comb, nbat->natoms*2*sizeof(*nbat->lj_comb));
            copy_lj_to_nbat_lj_comb_x4(nbat->nbfp_comb, nbat->type, nbat->natoms, nbat->lj_comb);

            // nbat->nout = 1 was set earlier

            // The checking code needs to know this size later
            forceArraySize = nbat->natoms * nbat->fstride;

            nbat->alloc((void**) &nbat->out[0].f,
                        forceArraySize * sizeof(*nbat->out[0].f));
            //            nbat->out[0].fshift = NULL;
            nbat->out[0].nV = 1;
            snew(nbat->out[0].Vc, nbat->out[0].nV);
            snew(nbat->out[0].Vvdw, nbat->out[0].nV);
            nbat->out[0].nVS = 0;
        }

        void build_nbnxn_ci_t(nbnxn_ci_t &ci)
        {
            ci.ci           = 0;
            ci.shift        = XYZ2IS(0, 0, 0) | NBNXN_CI_DO_LJ(0) | NBNXN_CI_DO_COUL(0);
            ci.cj_ind_start = 0;
            ci.cj_ind_end   = 1;
        }

        void build_nbnxn_cj_t(nbnxn_cj_t &cj)
        {
            cj.cj   = 1;
            cj.excl = ~0; // exclude no atoms
#ifdef GMX_CPU_ACCELERATION_IBM_QPX
            cj.interaction_mask_indices[0] = 15;
            cj.interaction_mask_indices[1] = 15;
            cj.interaction_mask_indices[2] = 15;
            cj.interaction_mask_indices[3] = 15;
#endif
        }

        void set_nbl_list()
        {
            nbl_list->nnbl = 1;
            snew(nbl_list->nbl, nbl_list->nnbl);
            snew(nbl_list->nbl[0], 1);
            set_nblist(nbl_list->nbl[0][0]);
        }

        void set_nblist(nbnxn_pairlist_t &nbl)
        {
            nbl.na_ci = 4;
            nbl.na_cj = 4;
            nbl.na_sc = 4;
            //            nbl.rlist = 1;
            nbl.nci = 1;
            snew(nbl.ci, nbl.nci);
            build_nbnxn_ci_t(nbl.ci[0]);
            nbl.ci_nalloc = 1;
            /*
               nbl.nsci = 0;
               nbl.sci = NULL;
               nbl.sci_nalloc = 0;
             */
            nbl.ncj = 1;
            snew(nbl.cj, nbl.ncj);
            build_nbnxn_cj_t(nbl.cj[0]);
            nbl.cj_nalloc = 1;
            /*
               nbl.ncj4 = 0;
               nbl.cj4 = NULL;
               nbl.cj4_nalloc = 0;
               nbl.nexcl = 0;
               nbl.excl = build_excl();
               nbl.excl_nalloc = 1;
             */
            nbl.nci_tot = 1;
        }

        void set_box_and_shifts()
        {
            // TODO use nbat->shift_vec?
            snew(shift_vec, SHIFTS);
            box[XX][XX] = 100;
            box[XX][YY] = 0;
            box[XX][ZZ] = 0;
            box[YY][XX] = 0;
            box[YY][YY] = 100;
            box[YY][ZZ] = 0;
            box[ZZ][XX] = 0;
            box[ZZ][YY] = 0;
            box[ZZ][ZZ] = 100;
            calc_shifts(box, shift_vec);
        }

        void callFunction()
        {
            nbnxn_kernel_simd_4xn(nbl_list, nbat, ic, ewald_excl,
                                  shift_vec,
                                  force_flags, clearF, fshift,
                                  Vc, Vvdw);
        }

        nbnxn_pairlist_set_t       *nbl_list;
        nbnxn_atomdata_t           *nbat;
        interaction_const_t        *ic;
        int                         ewald_excl;
        rvec                       *shift_vec;
        int                         force_flags;
        int                         clearF;
        real                       *fshift;
        real                       *Vc;
        real                       *Vvdw;
        matrix                      box;
        int                         forceArraySize;
};

t_forcerec *init_forcerec_for_test(int eel)
{
    t_forcerec *fr;
    snew(fr, 1);

    /* init modifiers */
    fr->vdw_modifier     = eintmodNONE;
    fr->coulomb_modifier = eintmodNONE;

    /* init cutoffs */
    fr->rlist      = 1;
    fr->rlistlong  = 1;
    fr->rvdw       = 1;
    fr->rcoulomb   = 1;

    /* init electrostatics */
    fr->eeltype          = eel;
    fr->epsilon_r        = 1;
    fr->epsfac           = 2.0; // keep numbers nice if possible

    /* init Ewald */
    real ewald_rtol = 1e-5;
    fr->ewaldcoeff = calc_ewaldcoeff(fr->rcoulomb, ewald_rtol);

    /* init reaction-field */
    fr->epsilon_rf = 80;
    fr->temp       = 298; // Don't care about Generalized RF
    fr->zsquare    = 0;   // Don't care about Generalized RF
    matrix dummy_box;     // Don't care about Generalized RF
    calc_rffac(NULL, fr->eeltype, fr->epsilon_r, fr->epsilon_rf,
               fr->rcoulomb, fr->temp, fr->zsquare, dummy_box,
               &fr->kappa, &fr->k_rf, &fr->c_rf);

    /* init Vdw */
    fr->vdwtype = evdwCUT;
    fr->ntype   = 1;
    snew(fr->nbfp, 2 * fr->ntype * fr->ntype);
    fr->nbfp[0] = 2.0e-6;  // C6
    fr->nbfp[1] = 4.2e-11; // C12

    fr->nbv = NULL;

    return fr;
}

int
lookUpEnumValueFromName(const char **defs, const char *name)
{
    int i;
    for (i = 0; (defs[i] != NULL); i++)
    {
        if (0 == gmx_strcasecmp_min(defs[i], name))
        {
            return i;
        }
    }
    GMX_THROW(::gmx::test::TestException("Invalid eel name"));
    return i;
}

typedef gmx::test::TestReferenceChecker Checker;

typedef ::testing::TestWithParam<const char *> Test4xnKernelDoublePrecision;
typedef ::testing::TestWithParam<const char *> Test4xnKernelSinglePrecision;

/* TODO Needing to name the test parameterization fixtures differently
   is nasty. The need arises from the fact that the two tests need
   different reference data. The reference data comes from a file
   named by the test fixture and test case. Those strings are used in
   macros that I understand are not expanded correctly for
   parameterized tests. This could be avoided if there was a way to
   embed both kinds of precision in the same reference data file and
   compare only with the right one. */
#ifdef GMX_DOUBLE
TEST_P(Test4xnKernelDoublePrecision, SingleClusterPairInteractionCalculatesEnergies)
#else
TEST_P(Test4xnKernelSinglePrecision, SingleClusterPairInteractionCalculatesEnergies)
#endif
{
    const char         *eelName = GetParam();
    int                 eel     = lookUpEnumValueFromName(eel_names, eelName);
    t_forcerec         *fr      = init_forcerec_for_test(eel);
    Test4xnKernelObject tester(fr);

    tester.set_nbat();
    tester.set_nbl_list();
    tester.force_flags = GMX_FORCE_ENERGY; // Forces are always calculated
    tester.clearF      = true;
    tester.set_box_and_shifts();
    gmx_omp_nthreads_set(emntNonbonded, 1);
    tester.callFunction();

    gmx::test::TestReferenceData results;
    Checker checker(results.rootChecker());
    checker.checkReal(tester.nbat->out[0].Vc[0], "Vc");
    checker.checkReal(tester.nbat->out[0].Vvdw[0], "Vvdw");

    sfree(fr->nbfp);
    sfree(fr);
}

//TODO test multiple i and j clusters
//TODO test the real Vc and Vvdw after reduction

#ifdef GMX_DOUBLE
TEST_P(Test4xnKernelDoublePrecision, SingleClusterPairInteractionCalculatesForces)
#else
TEST_P(Test4xnKernelSinglePrecision, SingleClusterPairInteractionCalculatesForces)
#endif
{
    const char         *eelName = GetParam();
    int                 eel     = lookUpEnumValueFromName(eel_names, eelName);
    t_forcerec         *fr      = init_forcerec_for_test(eel);
    Test4xnKernelObject tester(fr);

    tester.set_nbat();
    tester.set_nbl_list();
    tester.force_flags = GMX_FORCE_ENERGY; // Forces are always calculated
    tester.clearF      = true;
    tester.set_box_and_shifts();
    gmx_omp_nthreads_set(emntNonbonded, 1);
    tester.callFunction();

    gmx::test::TestReferenceData results;
    Checker checker(results.rootChecker());
    checker.checkSequenceArray(tester.forceArraySize, tester.nbat->out[0].f, "Forces");

    /*
       #ifdef GMX_DOUBLE
       real maxUlps = 5e9;
       #else
       real maxUlps = 1e6;
       #endif
     */

    sfree(fr->nbfp);
    sfree(fr);
}

// TODO Use strings in ::testing::ValuesIn(eel_names). Eventually, use
// all of them for which the kernels work.

#ifdef GMX_DOUBLE
INSTANTIATE_TEST_CASE_P(eelTypes,
                        Test4xnKernelDoublePrecision,
                            ::testing::Values("Cutoff",
                                              "PME",
                                              "Reaction-Field"
                                              ));
#else
INSTANTIATE_TEST_CASE_P(eelTypes,
                        Test4xnKernelSinglePrecision,
                            ::testing::Values("Cutoff",
                                              "PME",
                                              "Reaction-Field"
                                              ));
#endif

} // namespace

#endif
