/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2016- The GROMACS Authors
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
#include "gmxpre.h"

#include "gromacs/trajectory/trajectoryframe.h"

#include <cinttypes>
#include <cstdint>
#include <cstdio>

#include <algorithm>
#include <array>
#include <string>

#include "gromacs/math/veccompare.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/compare.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"

enum class PbcType : int;

void comp_frame(FILE* fp, t_trxframe* fr1, t_trxframe* fr2, gmx_bool bRMSD, real ftol, real abstol)
{
    fprintf(fp, "\n");
    cmp_int(fp, "not_ok", -1, fr1->not_ok, fr2->not_ok);
    cmp_int(fp, "natoms", -1, fr1->natoms, fr2->natoms);
    if (cmp_bool(fp, "bStep", -1, fr1->bStep, fr2->bStep))
    {
        cmp_int(fp, "step", -1, fr1->step, fr2->step);
    }
    cmp_int(fp, "step", -1, fr1->step, fr2->step);
    if (cmp_bool(fp, "bTime", -1, fr1->bTime, fr2->bTime))
    {
        cmp_real(fp, "time", -1, fr1->time, fr2->time, ftol, abstol);
    }
    if (cmp_bool(fp, "bLambda", -1, fr1->bLambda, fr2->bLambda))
    {
        cmp_real(fp, "lambda", -1, fr1->lambda, fr2->lambda, ftol, abstol);
    }
    if (cmp_bool(fp, "bAtoms", -1, fr1->bAtoms, fr2->bAtoms))
    {
        compareAtoms(fp, fr1->atoms, fr2->atoms, ftol, abstol);
    }
    if (cmp_bool(fp, "bPrec", -1, fr1->bPrec, fr2->bPrec))
    {
        cmp_real(fp, "prec", -1, fr1->prec, fr2->prec, ftol, abstol);
    }
    if (cmp_bool(fp, "bX", -1, fr1->bX, fr2->bX))
    {
        cmp_rvecs(fp, "x", std::min(fr1->natoms, fr2->natoms), fr1->x, fr2->x, bRMSD, ftol, abstol);
    }
    if (cmp_bool(fp, "bV", -1, fr1->bV, fr2->bV))
    {
        cmp_rvecs(fp, "v", std::min(fr1->natoms, fr2->natoms), fr1->v, fr2->v, bRMSD, ftol, abstol);
    }
    if (cmp_bool(fp, "bF", -1, fr1->bF, fr2->bF))
    {
        cmp_rvecs(fp, "f", std::min(fr1->natoms, fr2->natoms), fr1->f, fr2->f, bRMSD, ftol, abstol);
    }
    if (cmp_bool(fp, "bBox", -1, fr1->bBox, fr2->bBox))
    {
        cmp_rvecs(fp, "box", 3, fr1->box, fr2->box, FALSE, ftol, abstol);
    }
}

void done_frame(t_trxframe* frame)
{
    if (frame->atoms)
    {
        done_atom(frame->atoms);
        sfree(frame->atoms);
    }
    sfree(frame->x);
    sfree(frame->v);
    sfree(frame->f);
}

namespace gmx
{

TrajectoryFrame::TrajectoryFrame(const t_trxframe& frame) : frame_(frame), box_{ { { { 0 } } } }
{
    if (!frame.bStep)
    {
        GMX_THROW(APIError("Cannot handle trajectory frame that lacks a step number"));
    }
    if (!frame.bTime)
    {
        GMX_THROW(APIError("Cannot handle trajectory frame that lacks a time"));
    }
    if (frame.bBox)
    {
        for (int d = 0; d < DIM; ++d)
        {
            for (int dd = 0; dd < DIM; ++dd)
            {
                box_[d][dd] = frame.box[d][dd];
            }
        }
    }
}

std::string TrajectoryFrame::frameName() const
{
    return formatString("Time %f Step %" PRId64, frame_.time, frame_.step);
}

std::int64_t TrajectoryFrame::step() const
{
    return frame_.step;
}

double TrajectoryFrame::time() const
{
    return frame_.time;
}

PbcType TrajectoryFrame::pbc() const
{
    return frame_.pbcType;
}

ArrayRef<const RVec> TrajectoryFrame::x() const
{
    if (frame_.bX)
    {
        return arrayRefFromArray(reinterpret_cast<RVec*>(frame_.x), frame_.natoms);
    }
    else
    {
        return {};
    }
}

ArrayRef<const RVec> TrajectoryFrame::v() const
{
    if (frame_.bV)
    {
        return arrayRefFromArray(reinterpret_cast<RVec*>(frame_.v), frame_.natoms);
    }
    else
    {
        return {};
    }
}

ArrayRef<const RVec> TrajectoryFrame::f() const
{
    if (frame_.bF)
    {
        return arrayRefFromArray(reinterpret_cast<RVec*>(frame_.f), frame_.natoms);
    }
    else
    {
        return {};
    }
}

bool TrajectoryFrame::hasBox() const
{
    return frame_.bBox;
}

const BoxMatrix& TrajectoryFrame::box() const
{
    return box_;
}

} // namespace gmx
