/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017,2018, by the GROMACS development team, led by
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
#include <gromacs/utility/fatalerror.h>

#include "gromacs/math/paddedvector.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/veccompare.h"
#include "gromacs/mdtypes/awh-history.h"
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

rvec *makeRvecArray(gmx::ArrayRef<const gmx::RVec> v,
                    gmx::index                     n)
{
    GMX_ASSERT(v.size() >= n, "We can't copy more elements than the vector size");

    rvec *dest;

    snew(dest, n);

    const rvec *vPtr = as_rvec_array(v.data());
    for (gmx::index i = 0; i < n; i++)
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

                     baros_integral(0),
                     veta(0),
                     vol0(0),

                     ekinstate(),
                     hist(),
                     dfhist(nullptr),
                     awhHistory(nullptr),
                     ddp_count(0),
                     ddp_count_cg_gl(0)

{
    // It would be nicer to initialize these with {} or {{0}} in the
    // above initialization list, but uncrustify doesn't understand
    // that.
    // TODO Fix this if we switch to clang-format some time.
    lambda = {{ 0 }};
    clear_mat(box);
    clear_mat(box_rel);
    clear_mat(boxv);
    clear_mat(pres_prev);
    clear_mat(svir_prev);
    clear_mat(fvir_prev);
    
    intVector_.resize(getNumInt());
    int64Vector_.resize(getNumInt64());
    realVector_.resize(getNumReal());
    doubleVector_.resize(getNumDouble());

    intView_ = intVector_;
    int64View_ = int64Vector_;
    realView_ = realVector_;
    doubleView_ = doubleVector_;
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

const int t_state::checkpointVersion = 1;

CheckpointKeyword t_state::getKeyword()
{
    return CheckpointKeyword::state;
}
int t_state::getVersion()
{
    return checkpointVersion;
}

size_t t_state::getNumInt()
{
    return 8 + cg_gl.size();
}
size_t t_state::getNumInt64()
{
    return 0;
}
size_t t_state::getNumReal()
{
    return (lambda.size() +
            6*DIM*DIM + 2 +  // box (x3), pressure, virial matrices(x2), veta, vol0
            natoms*DIM*2);   // x, v (not including padding)
}
size_t t_state::getNumDouble()
{
    return (nosehoover_xi.size() + nosehoover_vxi.size() +
            nhpres_xi.size() + nhpres_vxi.size() +
            therm_integral.size() + 1);
}



ArrayRef<int> t_state::getIntView()
{
    return intView_;
}
ArrayRef<int64_t> t_state::getInt64View()
{
    return int64View_;
}
ArrayRef<real> t_state::getRealView()
{
    return realView_;
}
ArrayRef<double> t_state::getDoubleView()
{
    return doubleView_;
}

void t_state::notifyRead()
{
    // Obviously, we could just use getters and setters and actually use these vectors to store
    // our data instead of copying it forth and back, but that would require a lot of changes all
    // over the code. To demonstrate the approach, this should work.
    if (natoms != intView_[0])
    {
        gmx_fatal(FARGS, "Checkpoint file is for a system of %d atoms, while the current system "
                         "consists of %d atoms", intView_[0], natoms);
    }
    if (ngtc != intView_[1])
    {
        gmx_fatal(FARGS, "Checkpoint file is for a system of %d T-coupling groups, while the "
                         "current system consists of %d T-coupling groups", intView_[1], ngtc);
    }
    if (nnhpres != intView_[2])
    {
        gmx_fatal(FARGS, "Checkpoint file is for a system of %d NH-pressure-coupling variables, "
                         "while the current system consists of %d NH-pressure-coupling variables",
                  intView_[2], nnhpres);
    }

    /* write over whatever was read; we use the number of Nose-Hoover chains from the checkpoint */
    nhchainlength = intView_[3];

    /* need to keep this here to keep the tpr format working */
    init_gtc_state(this, ngtc, nnhpres, nhchainlength);

    if (flags != intView_[4])
    {
        gmx_fatal(FARGS, "Cannot change a simulation algorithm during a checkpoint restart. "
                         "Perhaps you should make a new .tpr with grompp -f new.mdp -t checkpoint.cpt");
    }

    // TODO: Is this really needed? We checked that natoms didn't change, no?
    // If reading, state->natoms was probably just read, so
    // allocations need to be managed. If writing, this won't change
    // anything that matters.
    state_change_natoms(this, natoms);

    fep_state = intView_[5];
    ddp_count = intView_[6];
    ddp_count_cg_gl = intView_[7];
    for (unsigned int i = 0; i < cg_gl.size(); ++i)
    {
        cg_gl[i] = intView_[7+i];
    }

    size_t currentReal;
    for (currentReal = 0; currentReal < efptNR; ++currentReal)
    {
        lambda[currentReal] = realView_[currentReal];
    }

    copyMatrix(realView_.subArray(currentReal, currentReal+DIM*DIM), box);
    currentReal += DIM*DIM;
    copyMatrix(realView_.subArray(currentReal, currentReal+DIM*DIM), box_rel);
    currentReal += DIM*DIM;
    copyMatrix(realView_.subArray(currentReal, currentReal+DIM*DIM), boxv);
    currentReal += DIM*DIM;
    copyMatrix(realView_.subArray(currentReal, currentReal+DIM*DIM), pres_prev);
    currentReal += DIM*DIM;
    copyMatrix(realView_.subArray(currentReal, currentReal+DIM*DIM), svir_prev);
    currentReal += DIM*DIM;
    copyMatrix(realView_.subArray(currentReal, currentReal+DIM*DIM), fvir_prev);
    currentReal += DIM*DIM;

    veta = realView_[currentReal++];
    vol0 = realView_[currentReal++];

    copyRvec(x, realView_.subArray(currentReal, DIM*natoms));
    currentReal += DIM*natoms;
    copyRvec(v, realView_.subArray(currentReal, DIM*natoms));

    size_t currentDouble = 0;
    gmx::copyVec(doubleView_.subArray(currentDouble, nosehoover_xi.size()), nosehoover_xi);
    currentDouble += nosehoover_xi.size();
    gmx::copyVec(doubleView_.subArray(currentDouble, nosehoover_vxi.size()), nosehoover_vxi);
    currentDouble += nosehoover_vxi.size();
    gmx::copyVec(doubleView_.subArray(currentDouble, nhpres_xi.size()), nhpres_xi);
    currentDouble += nhpres_xi.size();
    gmx::copyVec(doubleView_.subArray(currentDouble, nhpres_vxi.size()), nhpres_vxi);
    currentDouble += nhpres_vxi.size();
    gmx::copyVec(doubleView_.subArray(currentDouble, therm_integral.size()), therm_integral);
    currentDouble += therm_integral.size();
    baros_integral = doubleView_[currentDouble];

    // TODO: What about the distance / orientation restraint history data?
}

void t_state::notifyWrite()
{
    /* TODO: We will register the global state, which will be notified that checkpoint writing
     *       is imminent. Would it be our job to collect the state? For now, this is done in
     *       mdoutf_write_to_trajectory_files()...
     */

    // Obviously, we could just use getters and setters and actually use these vectors to store
    // our data instead of copying it forth and back, but that would require a lot of changes all
    // over the code. To demonstrate the approach, this should work.

    intView_[0] = natoms;
    intView_[1] = ngtc;
    intView_[2] = nnhpres;
    intView_[3] = nhchainlength;
    intView_[4] = flags;
    intView_[5] = fep_state;
    intView_[6] = ddp_count;
    intView_[7] = ddp_count_cg_gl;
    for (unsigned int i = 0; i < cg_gl.size(); ++i)
    {
        intView_[7+i] = cg_gl[i];
    }

    size_t currentReal;
    for (currentReal = 0; currentReal < efptNR; ++currentReal)
    {
        realView_[currentReal] = lambda[currentReal];
    }

    ArrayRef<real> subArrayReal;
    subArrayReal = realView_.subArray(currentReal, currentReal+DIM*DIM);
    copyMatrix(box, subArrayReal);
    currentReal += DIM*DIM;
    subArrayReal = realView_.subArray(currentReal, currentReal+DIM*DIM);
    copyMatrix(box_rel, subArrayReal);
    currentReal += DIM*DIM;
    subArrayReal = realView_.subArray(currentReal, currentReal+DIM*DIM);
    copyMatrix(boxv, subArrayReal);
    currentReal += DIM*DIM;
    subArrayReal = realView_.subArray(currentReal, currentReal+DIM*DIM);
    copyMatrix(pres_prev, subArrayReal);
    currentReal += DIM*DIM;
    subArrayReal = realView_.subArray(currentReal, currentReal+DIM*DIM);
    copyMatrix(svir_prev, subArrayReal);
    currentReal += DIM*DIM;
    subArrayReal = realView_.subArray(currentReal, currentReal+DIM*DIM);
    copyMatrix(fvir_prev, subArrayReal);
    currentReal += DIM*DIM;

    realView_[currentReal++] = veta;
    realView_[currentReal++] = vol0;

    copyRvec(realView_.subArray(currentReal, DIM*natoms), x);
    currentReal += DIM*natoms;
    copyRvec(realView_.subArray(currentReal, DIM*natoms), v);

    size_t currentDouble = 0;
    ArrayRef<double> subArrayDouble;
    subArrayDouble = doubleView_.subArray(currentDouble, nosehoover_xi.size());
    gmx::copyVec(nosehoover_xi, subArrayDouble);
    currentDouble += nosehoover_xi.size();
    subArrayDouble = doubleView_.subArray(currentDouble, nosehoover_vxi.size());
    gmx::copyVec(nosehoover_vxi, subArrayDouble);
    currentDouble += nosehoover_vxi.size();
    subArrayDouble = doubleView_.subArray(currentDouble, nhpres_xi.size());
    gmx::copyVec(nhpres_xi, subArrayDouble);
    currentDouble += nhpres_xi.size();
    subArrayDouble = doubleView_.subArray(currentDouble, nhpres_vxi.size());
    gmx::copyVec(nhpres_vxi, subArrayDouble);
    currentDouble += nhpres_vxi.size();
    subArrayDouble = doubleView_.subArray(currentDouble, therm_integral.size());
    gmx::copyVec(therm_integral, subArrayDouble);
    currentDouble += therm_integral.size();
    doubleView_[currentDouble] = baros_integral;

    // TODO: What about the distance / orientation restraint history data?
}
