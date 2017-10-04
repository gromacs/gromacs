/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017, by the GROMACS development team, led by
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
/* This file is completely threadsafe - keep it that way! */
#include "gmxpre.h"

#include "state.h"

#include <cstring>

#include <algorithm>

#include "gromacs/math/paddedvector.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/veccompare.h"
#include "gromacs/mdtypes/df_history.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/swaphistory.h"
#include "gromacs/pbcutil/boxutilities.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/utility/compare.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"

/* The source code in this file should be thread-safe.
      Please keep it that way. */

history_t::history_t() : disre_initf(0),
                         ndisrepairs(0),
                         disre_rm3tav(nullptr),
                         orire_initf(0),
                         norire_Dtav(0),
                         orire_Dtav(nullptr)
{
};

ekinstate_t::ekinstate_t() : bUpToDate(FALSE),
                             ekin_n(0),
                             ekinh(nullptr),
                             ekinf(nullptr),
                             ekinh_old(nullptr),
                             ekin_total(),
                             ekinscalef_nhc(),
                             ekinscaleh_nhc(),
                             vscale_nhc(),
                             dekindl(0),
                             mvcos(0)
{
    clear_mat(ekin_total);
};

void init_gtc_state(t_state *state, int ngtc, int nnhpres, int nhchainlength)
{
    state->ngtc          = ngtc;
    state->nnhpres       = nnhpres;
    state->nhchainlength = nhchainlength;
    state->nosehoover_xi.resize(state->nhchainlength*state->ngtc, 0);
    state->nosehoover_vxi.resize(state->nhchainlength*state->ngtc, 0);
    state->therm_integral.resize(state->ngtc, 0);
    state->baros_integral = 0.0;
    state->nhpres_xi.resize(state->nhchainlength*nnhpres, 0);
    state->nhpres_vxi.resize(state->nhchainlength*nnhpres, 0);
}


/* Checkpoint code relies on this function having no effect if
   state->natoms is > 0 and passed as natoms. */
void state_change_natoms(t_state *state, int natoms)
{
    state->natoms = natoms;

    /* We need padding, since we might use SIMD access */
    const size_t paddedSize = gmx::paddedRVecVectorSize(state->natoms);

    if (state->flags & (1 << estX))
    {
        state->x.resize(paddedSize);
    }
    if (state->flags & (1 << estV))
    {
        state->v.resize(paddedSize);
    }
    if (state->flags & (1 << estCGP))
    {
        state->cg_p.resize(paddedSize);
    }
}

void init_dfhist_state(t_state *state, int dfhistNumLambda)
{
    if (dfhistNumLambda > 0)
    {
        snew(state->dfhist, 1);
        init_df_history(state->dfhist, dfhistNumLambda);
    }
    else
    {
        state->dfhist = nullptr;
    }
}

void comp_state(const t_state *st1, const t_state *st2,
                gmx_bool bRMSD, real ftol, real abstol)
{
    int i, j, nc;

    fprintf(stdout, "comparing flags\n");
    cmp_int(stdout, "flags", -1, st1->flags, st2->flags);
    fprintf(stdout, "comparing box\n");
    cmp_rvecs(stdout, "box", DIM, st1->box, st2->box, FALSE, ftol, abstol);
    fprintf(stdout, "comparing box_rel\n");
    cmp_rvecs(stdout, "box_rel", DIM, st1->box_rel, st2->box_rel, FALSE, ftol, abstol);
    fprintf(stdout, "comparing boxv\n");
    cmp_rvecs(stdout, "boxv", DIM, st1->boxv, st2->boxv, FALSE, ftol, abstol);
    if (st1->flags & (1<<estSVIR_PREV))
    {
        fprintf(stdout, "comparing shake vir_prev\n");
        cmp_rvecs(stdout, "svir_prev", DIM, st1->svir_prev, st2->svir_prev, FALSE, ftol, abstol);
    }
    if (st1->flags & (1<<estFVIR_PREV))
    {
        fprintf(stdout, "comparing force vir_prev\n");
        cmp_rvecs(stdout, "fvir_prev", DIM, st1->fvir_prev, st2->fvir_prev, FALSE, ftol, abstol);
    }
    if (st1->flags & (1<<estPRES_PREV))
    {
        fprintf(stdout, "comparing prev_pres\n");
        cmp_rvecs(stdout, "pres_prev", DIM, st1->pres_prev, st2->pres_prev, FALSE, ftol, abstol);
    }
    cmp_int(stdout, "ngtc", -1, st1->ngtc, st2->ngtc);
    cmp_int(stdout, "nhchainlength", -1, st1->nhchainlength, st2->nhchainlength);
    if (st1->ngtc == st2->ngtc && st1->nhchainlength == st2->nhchainlength)
    {
        for (i = 0; i < st1->ngtc; i++)
        {
            nc = i*st1->nhchainlength;
            for (j = 0; j < nc; j++)
            {
                cmp_real(stdout, "nosehoover_xi",
                         i, st1->nosehoover_xi[nc+j], st2->nosehoover_xi[nc+j], ftol, abstol);
            }
        }
    }
    cmp_int(stdout, "nnhpres", -1, st1->nnhpres, st2->nnhpres);
    if (st1->nnhpres == st2->nnhpres && st1->nhchainlength == st2->nhchainlength)
    {
        for (i = 0; i < st1->nnhpres; i++)
        {
            nc = i*st1->nhchainlength;
            for (j = 0; j < nc; j++)
            {
                cmp_real(stdout, "nosehoover_xi",
                         i, st1->nhpres_xi[nc+j], st2->nhpres_xi[nc+j], ftol, abstol);
            }
        }
    }

    cmp_int(stdout, "natoms", -1, st1->natoms, st2->natoms);
    if (st1->natoms == st2->natoms)
    {
        if ((st1->flags & (1<<estX)) && (st2->flags & (1<<estX)))
        {
            fprintf(stdout, "comparing x\n");
            cmp_rvecs(stdout, "x", st1->natoms, as_rvec_array(st1->x.data()), as_rvec_array(st2->x.data()), bRMSD, ftol, abstol);
        }
        if ((st1->flags & (1<<estV)) && (st2->flags & (1<<estV)))
        {
            fprintf(stdout, "comparing v\n");
            cmp_rvecs(stdout, "v", st1->natoms, as_rvec_array(st1->v.data()), as_rvec_array(st2->v.data()), bRMSD, ftol, abstol);
        }
    }
}

rvec *getRvecArrayFromPaddedRVecVector(const PaddedRVecVector *v,
                                       unsigned int            n)
{
    GMX_ASSERT(v->size() >= n, "We can't copy more elements than the vector size");

    rvec *dest;

    snew(dest, n);

    const rvec *vPtr = as_rvec_array(v->data());
    for (unsigned int i = 0; i < n; i++)
    {
        copy_rvec(vPtr[i], dest[i]);
    }

    return dest;
}

t_state::t_state() : natoms(0),
                     ngtc(0),
                     nnhpres(0),
                     nhchainlength(0),
                     flags(0),
                     fep_state(0),
                     lambda(),
                     nosehoover_xi(),
                     nosehoover_vxi(),
                     nhpres_xi(),
                     nhpres_vxi(),
                     therm_integral(),
                     baros_integral(0),
                     veta(0),
                     vol0(0),
                     x(),
                     v(),
                     cg_p(),
                     ekinstate(),
                     hist(),
                     dfhist(nullptr),
                     ddp_count(0),
                     ddp_count_cg_gl(0),
                     cg_gl()
{
    // It would be nicer to initialize these with {} or {{0}} in the
    // above initialization list, but uncrustify doesn't understand
    // that.
    // TODO Fix this if we switch to clang-format some time.
    // cppcheck-suppress useInitializationList
    lambda = {{ 0 }};
    clear_mat(box);
    clear_mat(box_rel);
    clear_mat(boxv);
    clear_mat(pres_prev);
    clear_mat(svir_prev);
    clear_mat(fvir_prev);
}

void set_box_rel(const t_inputrec *ir, t_state *state)
{
    /* Make sure the box obeys the restrictions before we fix the ratios */
    correct_box(nullptr, 0, state->box, nullptr);

    clear_mat(state->box_rel);

    if (inputrecPreserveShape(ir))
    {
        const int ndim = ir->epct == epctSEMIISOTROPIC ? 2 : 3;
        do_box_rel(ndim, ir->deform, state->box_rel, state->box, true);
    }
}

void preserve_box_shape(const t_inputrec *ir, matrix box_rel, matrix box)
{
    if (inputrecPreserveShape(ir))
    {
        const int ndim = ir->epct == epctSEMIISOTROPIC ? 2 : 3;
        do_box_rel(ndim, ir->deform, box_rel, box, false);
    }
}
