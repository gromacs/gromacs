/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
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
/*! \internal \file
 * \brief
 * Implements routines in pbc.h.
 *
 * Utility functions for handling periodic boundary conditions.
 * Mainly used in analysis tools.
 */
#include "gmxpre.h"

#include "gromacs/pbcutil/pbc.h"

#include <cinttypes>
#include <cmath>
#include <cstdio>

#include <algorithm>
#include <filesystem>
#include <string>
#include <vector>

#include "gromacs/math/functions.h"
#include "gromacs/math/units.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vecdump.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/pbcutil/mshift.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"

const gmx::EnumerationArray<PbcType, std::string> c_pbcTypeNames = {
    { "xyz", "no", "xy", "screw", "unset" }
};

/* Skip 0 so we have more chance of detecting if we forgot to call set_pbc. */
enum
{
    epbcdxRECTANGULAR = 1,
    epbcdxTRICLINIC,
    epbcdx2D_RECT,
    epbcdx2D_TRIC,
    epbcdx1D_RECT,
    epbcdx1D_TRIC,
    epbcdxSCREW_RECT,
    epbcdxSCREW_TRIC,
    epbcdxNOPBC,
    epbcdxUNSUPPORTED
};

//! Margin for correction when the box is too skewed
static constexpr real sc_skewnessMargin = 1.001;
/*! \brief Margin factor for warning/error message
 *
 * The term 0.004 is just sufficient to not warn for a box that
 * is deformed with a shear rate of 0.001 ps^-1 over a pairlist
 * lifetime of 0.2 ps.
 *
 * Note that there can be errors in periodic images for atom pairs
 * close to half the box size. In nearly all cases the cut-off
 * distance is more then 0.5% shorter than half the minimum periodic
 * vector length, so there is no issue.
 */
static constexpr real sc_boxSkewnessMarginForWarning = sc_skewnessMargin + 0.004;

int numPbcDimensions(PbcType pbcType)
{
    int npbcdim = 0;

    switch (pbcType)
    {
        case PbcType::Unset:
            GMX_RELEASE_ASSERT(false,
                               "Number of PBC dimensions was requested before the PBC type set.");
            break;
        case PbcType::Xyz: npbcdim = 3; break;
        case PbcType::XY: npbcdim = 2; break;
        case PbcType::Screw: npbcdim = 3; break;
        case PbcType::No: npbcdim = 0; break;
        default: GMX_RELEASE_ASSERT(false, "Invalid pbcType in numPbcDimensions");
    }

    return npbcdim;
}

void dump_pbc(FILE* fp, t_pbc* pbc)
{
    rvec sum_box;

    fprintf(fp, "pbcTypeDX = %d\n", pbc->pbcTypeDX);
    pr_rvecs(fp, 0, "box", pbc->box, DIM);
    pr_rvecs(fp, 0, "fbox_diag", &pbc->fbox_diag, 1);
    pr_rvecs(fp, 0, "hbox_diag", &pbc->hbox_diag, 1);
    pr_rvecs(fp, 0, "mhbox_diag", &pbc->mhbox_diag, 1);
    rvec_add(pbc->hbox_diag, pbc->mhbox_diag, sum_box);
    pr_rvecs(fp, 0, "sum of the above two", &sum_box, 1);
    fprintf(fp, "max_cutoff2 = %g\n", pbc->max_cutoff2);
    fprintf(fp, "ntric_vec = %d\n", pbc->ntric_vec);
    if (pbc->ntric_vec > 0)
    {
        pr_ivecs(fp, 0, "tric_shift", pbc->tric_shift, pbc->ntric_vec, FALSE);
        pr_rvecs(fp, 0, "tric_vec", pbc->tric_vec, pbc->ntric_vec);
    }
}

const char* check_box(PbcType pbcType, const matrix box)
{
    const char* ptr;

    if (pbcType == PbcType::Unset)
    {
        pbcType = guessPbcType(box);
    }

    if (pbcType == PbcType::No)
    {
        return nullptr;
    }

    GMX_ASSERT(box != nullptr, "check_box requires a valid box unless pbcType is No");

    if (pbcType == PbcType::Xyz && box[XX][XX] == 0 && box[YY][YY] == 0 && box[ZZ][ZZ] == 0)
    {
        ptr = "Empty diagonal for a 3-dimensional periodic box";
    }
    else if (pbcType == PbcType::XY && box[XX][XX] == 0 && box[YY][YY] == 0)
    {
        ptr = "Empty diagonal for a 2-dimensional periodic box";
    }
    else if ((box[XX][YY] != 0) || (box[XX][ZZ] != 0) || (box[YY][ZZ] != 0))
    {
        ptr = "Only triclinic boxes with the first vector parallel to the x-axis and the second "
              "vector in the xy-plane are supported.";
    }
    else if (pbcType == PbcType::Screw && (box[YY][XX] != 0 || box[ZZ][XX] != 0))
    {
        ptr = "The unit cell can not have off-diagonal x-components with screw pbc";
    }
    else if (std::fabs(box[YY][XX]) > sc_boxSkewnessMarginForWarning * 0.5_real * box[XX][XX]
             || (pbcType != PbcType::XY
                 && (std::fabs(box[ZZ][XX]) > sc_boxSkewnessMarginForWarning * 0.5_real * box[XX][XX]
                     || std::fabs(box[ZZ][YY]) > sc_boxSkewnessMarginForWarning * 0.5_real * box[YY][YY])))
    {
        ptr = "Triclinic box is too skewed.";
    }
    else
    {
        ptr = nullptr;
    }

    return ptr;
}

void matrix_convert(matrix box, const rvec vec, const rvec angleInDegrees)
{
    rvec angle;
    svmul(gmx::c_deg2Rad, angleInDegrees, angle);
    box[XX][XX] = vec[XX];
    box[YY][XX] = vec[YY] * std::cos(angle[ZZ]);
    box[YY][YY] = vec[YY] * std::sin(angle[ZZ]);
    box[ZZ][XX] = vec[ZZ] * std::cos(angle[YY]);
    box[ZZ][YY] = vec[ZZ] * (std::cos(angle[XX]) - std::cos(angle[YY]) * std::cos(angle[ZZ]))
                  / std::sin(angle[ZZ]);
    box[ZZ][ZZ] =
            std::sqrt(gmx::square(vec[ZZ]) - box[ZZ][XX] * box[ZZ][XX] - box[ZZ][YY] * box[ZZ][YY]);
}

real max_cutoff2(PbcType pbcType, const matrix box)
{
    real       min_hv2, min_ss;
    const real oneFourth = 0.25;

    /* Physical limitation of the cut-off
     * by half the length of the shortest box vector.
     */
    min_hv2 = oneFourth * std::min(norm2(box[XX]), norm2(box[YY]));
    if (pbcType != PbcType::XY)
    {
        min_hv2 = std::min(min_hv2, oneFourth * norm2(box[ZZ]));
    }

    /* Limitation to the smallest diagonal element due to optimizations:
     * checking only linear combinations of single box-vectors (2 in x)
     * in the grid search and pbc_dx is a lot faster
     * than checking all possible combinations.
     */
    if (pbcType == PbcType::XY)
    {
        min_ss = std::min(box[XX][XX], box[YY][YY]);
    }
    else
    {
        min_ss = std::min(box[XX][XX], std::min(box[YY][YY] - std::fabs(box[ZZ][YY]), box[ZZ][ZZ]));
    }

    return std::min(min_hv2, min_ss * min_ss);
}

//! Set to true if warning has been printed
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
static gmx_bool bWarnedGuess = FALSE;

PbcType guessPbcType(const matrix box)
{
    PbcType pbcType;
    GMX_RELEASE_ASSERT(box != nullptr, "guessPbcType requires a valid box");

    if (box[XX][XX] > 0 && box[YY][YY] > 0 && box[ZZ][ZZ] > 0)
    {
        pbcType = PbcType::Xyz;
    }
    else if (box[XX][XX] > 0 && box[YY][YY] > 0 && box[ZZ][ZZ] == 0)
    {
        pbcType = PbcType::XY;
    }
    else if (box[XX][XX] == 0 && box[YY][YY] == 0 && box[ZZ][ZZ] == 0)
    {
        pbcType = PbcType::No;
    }
    else
    {
        if (!bWarnedGuess)
        {
            fprintf(stderr,
                    "WARNING: Unsupported box diagonal %f %f %f, "
                    "will not use periodic boundary conditions\n\n",
                    box[XX][XX],
                    box[YY][YY],
                    box[ZZ][ZZ]);
            bWarnedGuess = TRUE;
        }
        pbcType = PbcType::No;
    }

    if (debug)
    {
        fprintf(debug, "Guessed pbc = %s from the box matrix\n", c_pbcTypeNames[pbcType].c_str());
    }

    return pbcType;
}

//! Check if the box still obeys the restrictions, if not, correct it
static int correct_box_elem(FILE* fplog, const int64_t step, matrix box, const int v, const int d)
{
    int shift, maxshift = 10;

    shift = 0;

    /* correct elem d of vector v with vector d */
    while (box[v][d] > sc_skewnessMargin * 0.5_real * box[d][d])
    {
        if (fplog)
        {
            fprintf(fplog, "Step %" PRId64 ": correcting invalid box:\n", step);
            pr_rvecs(fplog, 0, "old box", box, DIM);
        }
        rvec_dec(box[v], box[d]);
        shift--;
        if (fplog)
        {
            pr_rvecs(fplog, 0, "new box", box, DIM);
        }
        if (shift <= -maxshift)
        {
            gmx_fatal(FARGS, "Box was shifted at least %d times. Please see log-file.", maxshift);
        }
    }
    while (box[v][d] < -sc_skewnessMargin * 0.5_real * box[d][d])
    {
        if (fplog)
        {
            fprintf(fplog, "Step %" PRId64 ": correcting invalid box:\n", step);
            pr_rvecs(fplog, 0, "old box", box, DIM);
        }
        rvec_inc(box[v], box[d]);
        shift++;
        if (fplog)
        {
            pr_rvecs(fplog, 0, "new box", box, DIM);
        }
        if (shift >= maxshift)
        {
            gmx_fatal(FARGS, "Box was shifted at least %d times. Please see log-file.", maxshift);
        }
    }

    return shift;
}

gmx_bool correct_box(FILE* fplog, const int64_t step, matrix box)
{
    int      zy, zx, yx;
    gmx_bool bCorrected;

    zy = correct_box_elem(fplog, step, box, ZZ, YY);
    zx = correct_box_elem(fplog, step, box, ZZ, XX);
    yx = correct_box_elem(fplog, step, box, YY, XX);

    bCorrected = ((zy != 0) || (zx != 0) || (yx != 0));

    return bCorrected;
}

//! Do the real arithmetic for filling the pbc struct
static void low_set_pbc(t_pbc* pbc, PbcType pbcType, const ivec dd_pbc, const matrix box)
{
    int         order[3] = { 0, -1, 1 };
    ivec        bPBC;
    const char* ptr;

    pbc->pbcType   = pbcType;
    pbc->ndim_ePBC = numPbcDimensions(pbcType);

    if (pbc->pbcType == PbcType::No)
    {
        pbc->pbcTypeDX = epbcdxNOPBC;

        return;
    }

    copy_mat(box, pbc->box);
    pbc->max_cutoff2 = 0;
    pbc->dim         = -1;
    pbc->ntric_vec   = 0;

    for (int i = 0; (i < DIM); i++)
    {
        pbc->fbox_diag[i]  = box[i][i];
        pbc->hbox_diag[i]  = pbc->fbox_diag[i] * 0.5_real;
        pbc->mhbox_diag[i] = -pbc->hbox_diag[i];
    }

    ptr = check_box(pbcType, box);
    if (ptr)
    {
        fprintf(stderr, "Warning: %s\n", ptr);
        pr_rvecs(stderr, 0, "         Box", box, DIM);
        fprintf(stderr, "         Can not fix pbc.\n\n");
        pbc->pbcTypeDX = epbcdxUNSUPPORTED;
    }
    else
    {
        if (pbcType == PbcType::Screw && nullptr != dd_pbc)
        {
            /* This combinated should never appear here */
            gmx_incons("low_set_pbc called with screw pbc and dd_nc != NULL");
        }

        int npbcdim = 0;
        for (int i = 0; i < DIM; i++)
        {
            if ((dd_pbc && dd_pbc[i] == 0) || (pbcType == PbcType::XY && i == ZZ))
            {
                bPBC[i] = 0;
            }
            else
            {
                bPBC[i] = 1;
                npbcdim++;
            }
        }
        switch (npbcdim)
        {
            case 1:
                /* 1D pbc is not an mdp option and it is therefore only used
                 * with single shifts.
                 */
                pbc->pbcTypeDX = epbcdx1D_RECT;
                for (int i = 0; i < DIM; i++)
                {
                    if (bPBC[i])
                    {
                        pbc->dim = i;
                    }
                }
                GMX_ASSERT(pbc->dim < DIM, "Dimension for PBC incorrect");
                for (int i = 0; i < pbc->dim; i++)
                {
                    if (pbc->box[pbc->dim][i] != 0)
                    {
                        pbc->pbcTypeDX = epbcdx1D_TRIC;
                    }
                }
                break;
            case 2:
                pbc->pbcTypeDX = epbcdx2D_RECT;
                for (int i = 0; i < DIM; i++)
                {
                    if (!bPBC[i])
                    {
                        pbc->dim = i;
                    }
                }
                for (int i = 0; i < DIM; i++)
                {
                    if (bPBC[i])
                    {
                        for (int j = 0; j < i; j++)
                        {
                            if (pbc->box[i][j] != 0)
                            {
                                pbc->pbcTypeDX = epbcdx2D_TRIC;
                            }
                        }
                    }
                }
                break;
            case 3:
                if (pbcType != PbcType::Screw)
                {
                    if (TRICLINIC(box))
                    {
                        pbc->pbcTypeDX = epbcdxTRICLINIC;
                    }
                    else
                    {
                        pbc->pbcTypeDX = epbcdxRECTANGULAR;
                    }
                }
                else
                {
                    pbc->pbcTypeDX = (box[ZZ][YY] == 0 ? epbcdxSCREW_RECT : epbcdxSCREW_TRIC);
                    if (pbc->pbcTypeDX == epbcdxSCREW_TRIC)
                    {
                        fprintf(stderr,
                                "Screw pbc is not yet implemented for triclinic boxes.\n"
                                "Can not fix pbc.\n");
                        pbc->pbcTypeDX = epbcdxUNSUPPORTED;
                    }
                }
                break;
            default: gmx_fatal(FARGS, "Incorrect number of pbc dimensions with DD: %d", npbcdim);
        }
        pbc->max_cutoff2 = max_cutoff2(pbcType, box);

        if (pbc->pbcTypeDX == epbcdxTRICLINIC || pbc->pbcTypeDX == epbcdx2D_TRIC
            || pbc->pbcTypeDX == epbcdxSCREW_TRIC)
        {
            if (debug)
            {
                pr_rvecs(debug, 0, "Box", box, DIM);
                fprintf(debug, "max cutoff %.3f\n", std::sqrt(pbc->max_cutoff2));
            }
            /* We will only need single shifts here */
            for (int kk = 0; kk < 3; kk++)
            {
                int k = order[kk];
                if (!bPBC[ZZ] && k != 0)
                {
                    continue;
                }
                for (int jj = 0; jj < 3; jj++)
                {
                    int j = order[jj];
                    if (!bPBC[YY] && j != 0)
                    {
                        continue;
                    }
                    for (int ii = 0; ii < 3; ii++)
                    {
                        int i = order[ii];
                        if (!bPBC[XX] && i != 0)
                        {
                            continue;
                        }
                        /* A shift is only useful when it is trilinic */
                        if (j != 0 || k != 0)
                        {
                            rvec trial;
                            rvec pos;
                            real d2old = 0;
                            real d2new = 0;

                            for (int d = 0; d < DIM; d++)
                            {
                                trial[d] = i * box[XX][d] + j * box[YY][d] + k * box[ZZ][d];
                                /* Choose the vector within the brick around 0,0,0 that
                                 * will become the shortest due to shift try.
                                 */
                                if (d == pbc->dim)
                                {
                                    trial[d] = 0;
                                    pos[d]   = 0;
                                }
                                else
                                {
                                    if (trial[d] < 0)
                                    {
                                        pos[d] = std::min(pbc->hbox_diag[d], -trial[d]);
                                    }
                                    else
                                    {
                                        pos[d] = std::max(-pbc->hbox_diag[d], -trial[d]);
                                    }
                                }
                                d2old += gmx::square(pos[d]);
                                d2new += gmx::square(pos[d] + trial[d]);
                            }
                            if (sc_skewnessMargin * d2new < d2old)
                            {
                                /* Check if shifts with one box vector less do better */
                                gmx_bool bUse = TRUE;
                                for (int dd = 0; dd < DIM; dd++)
                                {
                                    int shift = (dd == 0 ? i : (dd == 1 ? j : k));
                                    if (shift)
                                    {
                                        real d2new_c = 0;
                                        for (int d = 0; d < DIM; d++)
                                        {
                                            d2new_c += gmx::square(pos[d] + trial[d] - shift * box[dd][d]);
                                        }
                                        if (d2new_c <= sc_skewnessMargin * d2new)
                                        {
                                            bUse = FALSE;
                                        }
                                    }
                                }
                                if (bUse)
                                {
                                    /* Accept this shift vector. */
                                    if (pbc->ntric_vec >= MAX_NTRICVEC)
                                    {
                                        fprintf(stderr,
                                                "\nWARNING: Found more than %d triclinic "
                                                "correction vectors, ignoring some.\n"
                                                "  There is probably something wrong with your "
                                                "box.\n",
                                                MAX_NTRICVEC);
                                        pr_rvecs(stderr, 0, "         Box", box, DIM);
                                    }
                                    else
                                    {
                                        copy_rvec(trial, pbc->tric_vec[pbc->ntric_vec]);
                                        pbc->tric_shift[pbc->ntric_vec][XX] = i;
                                        pbc->tric_shift[pbc->ntric_vec][YY] = j;
                                        pbc->tric_shift[pbc->ntric_vec][ZZ] = k;
                                        pbc->ntric_vec++;

                                        if (debug)
                                        {
                                            fprintf(debug,
                                                    "  tricvec %2d = %2d %2d %2d  %5.2f %5.2f  "
                                                    "%5.2f %5.2f %5.2f  %5.2f %5.2f %5.2f\n",
                                                    pbc->ntric_vec,
                                                    i,
                                                    j,
                                                    k,
                                                    std::sqrt(d2old),
                                                    std::sqrt(d2new),
                                                    trial[XX],
                                                    trial[YY],
                                                    trial[ZZ],
                                                    pos[XX],
                                                    pos[YY],
                                                    pos[ZZ]);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void set_pbc(t_pbc* pbc, PbcType pbcType, const matrix box)
{
    if (pbcType == PbcType::Unset)
    {
        pbcType = guessPbcType(box);
    }

    low_set_pbc(pbc, pbcType, nullptr, box);
}

t_pbc* set_pbc_dd(t_pbc* pbc, PbcType pbcType, const gmx::IVec* domdecCells, gmx_bool bSingleDir, const matrix box)
{
    if (pbcType == PbcType::No)
    {
        pbc->pbcType = pbcType;

        return nullptr;
    }

    if (nullptr == domdecCells)
    {
        low_set_pbc(pbc, pbcType, nullptr, box);
    }
    else
    {
        if (pbcType == PbcType::Screw && (*domdecCells)[XX] > 1)
        {
            /* The rotation has been taken care of during coordinate communication */
            pbcType = PbcType::Xyz;
        }

        ivec usePBC;
        int  npbcdim = 0;
        for (int i = 0; i < DIM; i++)
        {
            usePBC[i] = 0;
            if ((*domdecCells)[i] <= (bSingleDir ? 1 : 2) && !(pbcType == PbcType::XY && i == ZZ))
            {
                usePBC[i] = 1;
                npbcdim++;
            }
        }

        if (npbcdim > 0)
        {
            low_set_pbc(pbc, pbcType, usePBC, box);
        }
        else
        {
            pbc->pbcType = PbcType::No;
        }
    }

    return (pbc->pbcType != PbcType::No ? pbc : nullptr);
}

void pbc_dx(const t_pbc* pbc, const rvec x1, const rvec x2, rvec dx)
{
    int      i, j;
    rvec     dx_start, trial;
    real     d2min, d2trial;
    gmx_bool bRot;

    rvec_sub(x1, x2, dx);

    switch (pbc->pbcTypeDX)
    {
        case epbcdxRECTANGULAR:
            for (i = 0; i < DIM; i++)
            {
                while (dx[i] > pbc->hbox_diag[i])
                {
                    dx[i] -= pbc->fbox_diag[i];
                }
                while (dx[i] <= pbc->mhbox_diag[i])
                {
                    dx[i] += pbc->fbox_diag[i];
                }
            }
            break;
        case epbcdxTRICLINIC:
            for (i = DIM - 1; i >= 0; i--)
            {
                while (dx[i] > pbc->hbox_diag[i])
                {
                    for (j = i; j >= 0; j--)
                    {
                        dx[j] -= pbc->box[i][j];
                    }
                }
                while (dx[i] <= pbc->mhbox_diag[i])
                {
                    for (j = i; j >= 0; j--)
                    {
                        dx[j] += pbc->box[i][j];
                    }
                }
            }
            /* dx is the distance in a rectangular box */
            d2min = norm2(dx);
            if (d2min > pbc->max_cutoff2)
            {
                copy_rvec(dx, dx_start);
                d2min = norm2(dx);
                /* Now try all possible shifts, when the distance is within max_cutoff
                 * it must be the shortest possible distance.
                 */
                i = 0;
                while ((d2min > pbc->max_cutoff2) && (i < pbc->ntric_vec))
                {
                    rvec_add(dx_start, pbc->tric_vec[i], trial);
                    d2trial = norm2(trial);
                    if (d2trial < d2min)
                    {
                        copy_rvec(trial, dx);
                        d2min = d2trial;
                    }
                    i++;
                }
            }
            break;
        case epbcdx2D_RECT:
            for (i = 0; i < DIM; i++)
            {
                if (i != pbc->dim)
                {
                    while (dx[i] > pbc->hbox_diag[i])
                    {
                        dx[i] -= pbc->fbox_diag[i];
                    }
                    while (dx[i] <= pbc->mhbox_diag[i])
                    {
                        dx[i] += pbc->fbox_diag[i];
                    }
                }
            }
            break;
        case epbcdx2D_TRIC:
            d2min = 0;
            for (i = DIM - 1; i >= 0; i--)
            {
                if (i != pbc->dim)
                {
                    while (dx[i] > pbc->hbox_diag[i])
                    {
                        for (j = i; j >= 0; j--)
                        {
                            dx[j] -= pbc->box[i][j];
                        }
                    }
                    while (dx[i] <= pbc->mhbox_diag[i])
                    {
                        for (j = i; j >= 0; j--)
                        {
                            dx[j] += pbc->box[i][j];
                        }
                    }
                    d2min += dx[i] * dx[i];
                }
            }
            if (d2min > pbc->max_cutoff2)
            {
                copy_rvec(dx, dx_start);
                d2min = norm2(dx);
                /* Now try all possible shifts, when the distance is within max_cutoff
                 * it must be the shortest possible distance.
                 */
                i = 0;
                while ((d2min > pbc->max_cutoff2) && (i < pbc->ntric_vec))
                {
                    rvec_add(dx_start, pbc->tric_vec[i], trial);
                    d2trial = 0;
                    for (j = 0; j < DIM; j++)
                    {
                        if (j != pbc->dim)
                        {
                            d2trial += trial[j] * trial[j];
                        }
                    }
                    if (d2trial < d2min)
                    {
                        copy_rvec(trial, dx);
                        d2min = d2trial;
                    }
                    i++;
                }
            }
            break;
        case epbcdxSCREW_RECT:
            /* The shift definition requires x first */
            bRot = FALSE;
            while (dx[XX] > pbc->hbox_diag[XX])
            {
                dx[XX] -= pbc->fbox_diag[XX];
                bRot = !bRot;
            }
            while (dx[XX] <= pbc->mhbox_diag[XX])
            {
                dx[XX] += pbc->fbox_diag[YY];
                bRot = !bRot;
            }
            if (bRot)
            {
                /* Rotate around the x-axis in the middle of the box */
                dx[YY] = pbc->box[YY][YY] - x1[YY] - x2[YY];
                dx[ZZ] = pbc->box[ZZ][ZZ] - x1[ZZ] - x2[ZZ];
            }
            /* Normal pbc for y and z */
            for (i = YY; i <= ZZ; i++)
            {
                while (dx[i] > pbc->hbox_diag[i])
                {
                    dx[i] -= pbc->fbox_diag[i];
                }
                while (dx[i] <= pbc->mhbox_diag[i])
                {
                    dx[i] += pbc->fbox_diag[i];
                }
            }
            break;
        case epbcdxNOPBC:
        case epbcdxUNSUPPORTED: break;
        default: gmx_fatal(FARGS, "Internal error in pbc_dx, set_pbc has not been called");
    }
}

int pbc_dx_aiuc(const t_pbc* pbc, const rvec x1, const rvec x2, rvec dx)
{
    int  i, j, is;
    rvec dx_start, trial;
    real d2min, d2trial;
    ivec ishift, ishift_start;

    rvec_sub(x1, x2, dx);
    clear_ivec(ishift);

    switch (pbc->pbcTypeDX)
    {
        case epbcdxRECTANGULAR:
            for (i = 0; i < DIM; i++)
            {
                if (dx[i] > pbc->hbox_diag[i])
                {
                    dx[i] -= pbc->fbox_diag[i];
                    ishift[i]--;
                }
                else if (dx[i] <= pbc->mhbox_diag[i])
                {
                    dx[i] += pbc->fbox_diag[i];
                    ishift[i]++;
                }
            }
            break;
        case epbcdxTRICLINIC:
            /* For triclinic boxes the performance difference between
             * if/else and two while loops is negligible.
             * However, the while version can cause extreme delays
             * before a simulation crashes due to large forces which
             * can cause unlimited displacements.
             * Also allowing multiple shifts would index fshift beyond bounds.
             */
            for (i = DIM - 1; i >= 1; i--)
            {
                if (dx[i] > pbc->hbox_diag[i])
                {
                    for (j = i; j >= 0; j--)
                    {
                        dx[j] -= pbc->box[i][j];
                    }
                    ishift[i]--;
                }
                else if (dx[i] <= pbc->mhbox_diag[i])
                {
                    for (j = i; j >= 0; j--)
                    {
                        dx[j] += pbc->box[i][j];
                    }
                    ishift[i]++;
                }
            }
            /* Allow 2 shifts in x */
            if (dx[XX] > pbc->hbox_diag[XX])
            {
                dx[XX] -= pbc->fbox_diag[XX];
                ishift[XX]--;
                if (dx[XX] > pbc->hbox_diag[XX])
                {
                    dx[XX] -= pbc->fbox_diag[XX];
                    ishift[XX]--;
                }
            }
            else if (dx[XX] <= pbc->mhbox_diag[XX])
            {
                dx[XX] += pbc->fbox_diag[XX];
                ishift[XX]++;
                if (dx[XX] <= pbc->mhbox_diag[XX])
                {
                    dx[XX] += pbc->fbox_diag[XX];
                    ishift[XX]++;
                }
            }
            /* dx is the distance in a rectangular box */
            d2min = norm2(dx);
            if (d2min > pbc->max_cutoff2)
            {
                copy_rvec(dx, dx_start);
                copy_ivec(ishift, ishift_start);
                d2min = norm2(dx);
                /* Now try all possible shifts, when the distance is within max_cutoff
                 * it must be the shortest possible distance.
                 */
                i = 0;
                while ((d2min > pbc->max_cutoff2) && (i < pbc->ntric_vec))
                {
                    rvec_add(dx_start, pbc->tric_vec[i], trial);
                    d2trial = norm2(trial);
                    if (d2trial < d2min)
                    {
                        copy_rvec(trial, dx);
                        ivec_add(ishift_start, pbc->tric_shift[i], ishift);
                        d2min = d2trial;
                    }
                    i++;
                }
            }
            break;
        case epbcdx2D_RECT:
            for (i = 0; i < DIM; i++)
            {
                if (i != pbc->dim)
                {
                    if (dx[i] > pbc->hbox_diag[i])
                    {
                        dx[i] -= pbc->fbox_diag[i];
                        ishift[i]--;
                    }
                    else if (dx[i] <= pbc->mhbox_diag[i])
                    {
                        dx[i] += pbc->fbox_diag[i];
                        ishift[i]++;
                    }
                }
            }
            break;
        case epbcdx2D_TRIC:
            d2min = 0;
            for (i = DIM - 1; i >= 1; i--)
            {
                if (i != pbc->dim)
                {
                    if (dx[i] > pbc->hbox_diag[i])
                    {
                        for (j = i; j >= 0; j--)
                        {
                            dx[j] -= pbc->box[i][j];
                        }
                        ishift[i]--;
                    }
                    else if (dx[i] <= pbc->mhbox_diag[i])
                    {
                        for (j = i; j >= 0; j--)
                        {
                            dx[j] += pbc->box[i][j];
                        }
                        ishift[i]++;
                    }
                    d2min += dx[i] * dx[i];
                }
            }
            if (pbc->dim != XX)
            {
                /* Allow 2 shifts in x */
                if (dx[XX] > pbc->hbox_diag[XX])
                {
                    dx[XX] -= pbc->fbox_diag[XX];
                    ishift[XX]--;
                    if (dx[XX] > pbc->hbox_diag[XX])
                    {
                        dx[XX] -= pbc->fbox_diag[XX];
                        ishift[XX]--;
                    }
                }
                else if (dx[XX] <= pbc->mhbox_diag[XX])
                {
                    dx[XX] += pbc->fbox_diag[XX];
                    ishift[XX]++;
                    if (dx[XX] <= pbc->mhbox_diag[XX])
                    {
                        dx[XX] += pbc->fbox_diag[XX];
                        ishift[XX]++;
                    }
                }
                d2min += dx[XX] * dx[XX];
            }
            if (d2min > pbc->max_cutoff2)
            {
                copy_rvec(dx, dx_start);
                copy_ivec(ishift, ishift_start);
                /* Now try all possible shifts, when the distance is within max_cutoff
                 * it must be the shortest possible distance.
                 */
                i = 0;
                while ((d2min > pbc->max_cutoff2) && (i < pbc->ntric_vec))
                {
                    rvec_add(dx_start, pbc->tric_vec[i], trial);
                    d2trial = 0;
                    for (j = 0; j < DIM; j++)
                    {
                        if (j != pbc->dim)
                        {
                            d2trial += trial[j] * trial[j];
                        }
                    }
                    if (d2trial < d2min)
                    {
                        copy_rvec(trial, dx);
                        ivec_add(ishift_start, pbc->tric_shift[i], ishift);
                        d2min = d2trial;
                    }
                    i++;
                }
            }
            break;
        case epbcdx1D_RECT:
            i = pbc->dim;
            if (dx[i] > pbc->hbox_diag[i])
            {
                dx[i] -= pbc->fbox_diag[i];
                ishift[i]--;
            }
            else if (dx[i] <= pbc->mhbox_diag[i])
            {
                dx[i] += pbc->fbox_diag[i];
                ishift[i]++;
            }
            break;
        case epbcdx1D_TRIC:
            i = pbc->dim;
            if (dx[i] > pbc->hbox_diag[i])
            {
                rvec_dec(dx, pbc->box[i]);
                ishift[i]--;
            }
            else if (dx[i] <= pbc->mhbox_diag[i])
            {
                rvec_inc(dx, pbc->box[i]);
                ishift[i]++;
            }
            break;
        case epbcdxSCREW_RECT:
            /* The shift definition requires x first */
            if (dx[XX] > pbc->hbox_diag[XX])
            {
                dx[XX] -= pbc->fbox_diag[XX];
                ishift[XX]--;
            }
            else if (dx[XX] <= pbc->mhbox_diag[XX])
            {
                dx[XX] += pbc->fbox_diag[XX];
                ishift[XX]++;
            }
            if (ishift[XX] == 1 || ishift[XX] == -1)
            {
                /* Rotate around the x-axis in the middle of the box */
                dx[YY] = pbc->box[YY][YY] - x1[YY] - x2[YY];
                dx[ZZ] = pbc->box[ZZ][ZZ] - x1[ZZ] - x2[ZZ];
            }
            /* Normal pbc for y and z */
            for (i = YY; i <= ZZ; i++)
            {
                if (dx[i] > pbc->hbox_diag[i])
                {
                    dx[i] -= pbc->fbox_diag[i];
                    ishift[i]--;
                }
                else if (dx[i] <= pbc->mhbox_diag[i])
                {
                    dx[i] += pbc->fbox_diag[i];
                    ishift[i]++;
                }
            }
            break;
        case epbcdxNOPBC:
        case epbcdxUNSUPPORTED: break;
        default:
            gmx_fatal(FARGS,
                      "Internal error in pbc_dx_aiuc, set_pbc_dd or set_pbc has not been called");
    }

    is = gmx::ivecToShiftIndex(ishift);
    if (debug)
    {
        range_check_mesg(is, 0, gmx::c_numShiftVectors, "PBC shift vector index range check.");
    }

    return is;
}

//! Compute distance vector in double precision
void pbc_dx_d(const t_pbc* pbc, const dvec x1, const dvec x2, dvec dx)
{
    int      i, j;
    dvec     dx_start, trial;
    double   d2min, d2trial;
    gmx_bool bRot;

    dvec_sub(x1, x2, dx);

    switch (pbc->pbcTypeDX)
    {
        case epbcdxRECTANGULAR:
        case epbcdx2D_RECT:
            for (i = 0; i < DIM; i++)
            {
                if (i != pbc->dim)
                {
                    while (dx[i] > pbc->hbox_diag[i])
                    {
                        dx[i] -= pbc->fbox_diag[i];
                    }
                    while (dx[i] <= pbc->mhbox_diag[i])
                    {
                        dx[i] += pbc->fbox_diag[i];
                    }
                }
            }
            break;
        case epbcdxTRICLINIC:
        case epbcdx2D_TRIC:
            d2min = 0;
            for (i = DIM - 1; i >= 0; i--)
            {
                if (i != pbc->dim)
                {
                    while (dx[i] > pbc->hbox_diag[i])
                    {
                        for (j = i; j >= 0; j--)
                        {
                            dx[j] -= pbc->box[i][j];
                        }
                    }
                    while (dx[i] <= pbc->mhbox_diag[i])
                    {
                        for (j = i; j >= 0; j--)
                        {
                            dx[j] += pbc->box[i][j];
                        }
                    }
                    d2min += dx[i] * dx[i];
                }
            }
            if (d2min > pbc->max_cutoff2)
            {
                copy_dvec(dx, dx_start);
                /* Now try all possible shifts, when the distance is within max_cutoff
                 * it must be the shortest possible distance.
                 */
                i = 0;
                while ((d2min > pbc->max_cutoff2) && (i < pbc->ntric_vec))
                {
                    for (j = 0; j < DIM; j++)
                    {
                        trial[j] = dx_start[j] + pbc->tric_vec[i][j];
                    }
                    d2trial = 0;
                    for (j = 0; j < DIM; j++)
                    {
                        if (j != pbc->dim)
                        {
                            d2trial += trial[j] * trial[j];
                        }
                    }
                    if (d2trial < d2min)
                    {
                        copy_dvec(trial, dx);
                        d2min = d2trial;
                    }
                    i++;
                }
            }
            break;
        case epbcdxSCREW_RECT:
            /* The shift definition requires x first */
            bRot = FALSE;
            while (dx[XX] > pbc->hbox_diag[XX])
            {
                dx[XX] -= pbc->fbox_diag[XX];
                bRot = !bRot;
            }
            while (dx[XX] <= pbc->mhbox_diag[XX])
            {
                dx[XX] += pbc->fbox_diag[YY];
                bRot = !bRot;
            }
            if (bRot)
            {
                /* Rotate around the x-axis in the middle of the box */
                dx[YY] = pbc->box[YY][YY] - x1[YY] - x2[YY];
                dx[ZZ] = pbc->box[ZZ][ZZ] - x1[ZZ] - x2[ZZ];
            }
            /* Normal pbc for y and z */
            for (i = YY; i <= ZZ; i++)
            {
                while (dx[i] > pbc->hbox_diag[i])
                {
                    dx[i] -= pbc->fbox_diag[i];
                }
                while (dx[i] <= pbc->mhbox_diag[i])
                {
                    dx[i] += pbc->fbox_diag[i];
                }
            }
            break;
        case epbcdxNOPBC:
        case epbcdxUNSUPPORTED: break;
        default: gmx_fatal(FARGS, "Internal error in pbc_dx, set_pbc has not been called");
    }
}

void calc_shifts(const matrix box, gmx::ArrayRef<gmx::RVec> shift_vec)
{
    for (int n = 0, m = -gmx::c_dBoxZ; m <= gmx::c_dBoxZ; m++)
    {
        for (int l = -gmx::c_dBoxY; l <= gmx::c_dBoxY; l++)
        {
            for (int k = -gmx::c_dBoxX; k <= gmx::c_dBoxX; k++, n++)
            {
                for (int d = 0; d < DIM; d++)
                {
                    shift_vec[n][d] = k * box[XX][d] + l * box[YY][d] + m * box[ZZ][d];
                }
            }
        }
    }
}

void calc_box_center(int ecenter, const matrix box, rvec box_center)
{
    int d, m;

    clear_rvec(box_center);
    switch (ecenter)
    {
        case ecenterTRIC:
            for (m = 0; (m < DIM); m++)
            {
                for (d = 0; d < DIM; d++)
                {
                    box_center[d] += 0.5_real * box[m][d];
                }
            }
            break;
        case ecenterRECT:
            for (d = 0; d < DIM; d++)
            {
                box_center[d] = 0.5_real * box[d][d];
            }
            break;
        case ecenterZERO: break;
        default: gmx_fatal(FARGS, "Unsupported value %d for ecenter", ecenter);
    }
}

void calc_triclinic_images(const matrix box, rvec img[])
{
    int i;

    /* Calculate 3 adjacent images in the xy-plane */
    copy_rvec(box[0], img[0]);
    copy_rvec(box[1], img[1]);
    if (img[1][XX] < 0)
    {
        svmul(-1, img[1], img[1]);
    }
    rvec_sub(img[1], img[0], img[2]);

    /* Get the next 3 in the xy-plane as mirror images */
    for (i = 0; i < 3; i++)
    {
        svmul(-1, img[i], img[3 + i]);
    }

    /* Calculate the first 4 out of xy-plane images */
    copy_rvec(box[2], img[6]);
    if (img[6][XX] < 0)
    {
        svmul(-1, img[6], img[6]);
    }
    for (i = 0; i < 3; i++)
    {
        rvec_add(img[6], img[i + 1], img[7 + i]);
    }

    /* Mirror the last 4 from the previous in opposite rotation */
    for (i = 0; i < 4; i++)
    {
        svmul(-1, img[6 + (2 + i) % 4], img[10 + i]);
    }
}

void calc_compact_unitcell_vertices(int ecenter, const matrix box, rvec vert[])
{
    rvec       img[NTRICIMG], box_center;
    int        n, i, j, tmp[4], d;
    const real oneFourth = 0.25;

    calc_triclinic_images(box, img);

    n = 0;
    for (i = 2; i <= 5; i += 3)
    {
        tmp[0] = i - 1;
        if (i == 2)
        {
            tmp[1] = 8;
        }
        else
        {
            tmp[1] = 6;
        }
        tmp[2] = (i + 1) % 6;
        tmp[3] = tmp[1] + 4;
        for (j = 0; j < 4; j++)
        {
            for (d = 0; d < DIM; d++)
            {
                vert[n][d] = img[i][d] + img[tmp[j]][d] + img[tmp[(j + 1) % 4]][d];
            }
            n++;
        }
    }
    for (i = 7; i <= 13; i += 6)
    {
        tmp[0] = (i - 7) / 2;
        tmp[1] = tmp[0] + 1;
        if (i == 7)
        {
            tmp[2] = 8;
        }
        else
        {
            tmp[2] = 10;
        }
        tmp[3] = i - 1;
        for (j = 0; j < 4; j++)
        {
            for (d = 0; d < DIM; d++)
            {
                vert[n][d] = img[i][d] + img[tmp[j]][d] + img[tmp[(j + 1) % 4]][d];
            }
            n++;
        }
    }
    for (i = 9; i <= 11; i += 2)
    {
        if (i == 9)
        {
            tmp[0] = 3;
        }
        else
        {
            tmp[0] = 0;
        }
        tmp[1] = tmp[0] + 1;
        if (i == 9)
        {
            tmp[2] = 6;
        }
        else
        {
            tmp[2] = 12;
        }
        tmp[3] = i - 1;
        for (j = 0; j < 4; j++)
        {
            for (d = 0; d < DIM; d++)
            {
                vert[n][d] = img[i][d] + img[tmp[j]][d] + img[tmp[(j + 1) % 4]][d];
            }
            n++;
        }
    }

    calc_box_center(ecenter, box, box_center);
    for (i = 0; i < NCUCVERT; i++)
    {
        for (d = 0; d < DIM; d++)
        {
            vert[i][d] = vert[i][d] * oneFourth + box_center[d];
        }
    }
}

int* compact_unitcell_edges()
{
    /* this is an index in vert[] (see calc_box_vertices) */
    /*static int edge[NCUCEDGE*2];*/
    int*             edge;
    static const int hexcon[24] = { 0, 9,  1, 19, 2, 15, 3,  21, 4,  17, 5,  11,
                                    6, 23, 7, 13, 8, 20, 10, 18, 12, 16, 14, 22 };
    int              e, i, j;

    snew(edge, NCUCEDGE * 2);

    e = 0;
    for (i = 0; i < 6; i++)
    {
        for (j = 0; j < 4; j++)
        {
            edge[e++] = 4 * i + j;
            edge[e++] = 4 * i + (j + 1) % 4;
        }
    }
    for (i = 0; i < 12 * 2; i++)
    {
        edge[e++] = hexcon[i];
    }

    return edge;
}

template<bool haveBoxDeformation>
static void putAtomsInBoxTemplated(PbcType                  pbcType,
                                   const matrix             box,
                                   const matrix             boxDeformation,
                                   gmx::ArrayRef<gmx::RVec> x,
                                   gmx::ArrayRef<gmx::RVec> v)
{
    if constexpr (haveBoxDeformation)
    {
        GMX_ASSERT(v.size() == x.size(), "Need velocities for box deformation");
    }

    // NOLINTNEXTLINE(readability-misleading-indentation)
    int npbcdim;

    if (pbcType == PbcType::Screw)
    {
        gmx_fatal(FARGS, "Sorry, %s pbc is not yet supported", c_pbcTypeNames[pbcType].c_str());
    }

    if (pbcType == PbcType::XY)
    {
        npbcdim = 2;
    }
    else
    {
        npbcdim = 3;
    }

    gmx::RVec invBox;
    for (int m = 0; m < npbcdim; ++m)
    {
        invBox[m] = 1 / box[m][m];
    }

    if (TRICLINIC(box))
    {
        for (gmx::Index i = 0; i < x.ssize(); i++)
        {
            for (int m = npbcdim - 1; m >= 0; m--)
            {
                const auto boxVectorShift = std::floor(x[i][m] * invBox[m]);
                for (int d = 0; d <= m; d++)
                {
                    x[i][d] -= boxVectorShift * box[m][d];
                    if constexpr (haveBoxDeformation)
                    {
                        v[i][d] -= boxVectorShift * boxDeformation[m][d];
                    }
                }
            }
        }
    }
    else
    {
        for (gmx::Index i = 0; i < x.ssize(); i++)
        {
            for (int d = 0; d < npbcdim; d++)
            {
                const auto boxVectorShift = std::floor(x[i][d] * invBox[d]);
                x[i][d] -= boxVectorShift * box[d][d];
                if constexpr (haveBoxDeformation)
                {
                    for (int d2 = 0; d2 <= d; d2++)
                    {
                        v[i][d2] -= boxVectorShift * boxDeformation[d][d2];
                    }
                }
            }
        }
    }

    if constexpr (!haveBoxDeformation)
    {
        GMX_UNUSED_VALUE(v);
    }
}


void put_atoms_in_box(PbcType pbcType, const matrix box, gmx::ArrayRef<gmx::RVec> x)
{
    putAtomsInBoxTemplated<false>(pbcType, box, nullptr, x, {});
}

void put_atoms_in_box_omp(PbcType                  pbcType,
                          const matrix             box,
                          const bool               haveBoxDeformation,
                          const matrix             boxDeformation,
                          gmx::ArrayRef<gmx::RVec> x,
                          gmx::ArrayRef<gmx::RVec> v,
                          gmx_unused int           nth)
{
#pragma omp parallel for num_threads(nth) schedule(static)
    for (int t = 0; t < nth; t++)
    {
        try
        {
            size_t natoms = x.size();
            size_t offset = (natoms * t) / nth;
            size_t len    = (natoms * (t + 1)) / nth - offset;
            if (haveBoxDeformation)
            {
                putAtomsInBoxTemplated<true>(
                        pbcType, box, boxDeformation, x.subArray(offset, len), v.subArray(offset, len));
            }
            else
            {
                putAtomsInBoxTemplated<false>(pbcType, box, boxDeformation, x.subArray(offset, len), {});
            }
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
    }
}

void put_atoms_in_triclinic_unitcell(int ecenter, const matrix box, gmx::ArrayRef<gmx::RVec> x)
{
    rvec box_center, shift_center;
    real shm01, shm02, shm12, shift;
    int  m, d;

    calc_box_center(ecenter, box, box_center);

    /* The product of matrix shm with a coordinate gives the shift vector
       which is required determine the periodic cell position */
    shm01 = box[1][0] / box[1][1];
    shm02 = (box[1][1] * box[2][0] - box[2][1] * box[1][0]) / (box[1][1] * box[2][2]);
    shm12 = box[2][1] / box[2][2];

    clear_rvec(shift_center);
    for (d = 0; d < DIM; d++)
    {
        rvec_inc(shift_center, box[d]);
    }
    svmul(0.5_real, shift_center, shift_center);
    rvec_sub(box_center, shift_center, shift_center);

    shift_center[0] = shm01 * shift_center[1] + shm02 * shift_center[2];
    shift_center[1] = shm12 * shift_center[2];
    shift_center[2] = 0;

    for (gmx::Index i = 0; (i < x.ssize()); ++i)
    {
        for (m = DIM - 1; m >= 0; m--)
        {
            shift = shift_center[m];
            if (m == 0)
            {
                shift += shm01 * x[i][1] + shm02 * x[i][2];
            }
            else if (m == 1)
            {
                shift += shm12 * x[i][2];
            }
            while (x[i][m] - shift < 0)
            {
                for (d = 0; d <= m; d++)
                {
                    x[i][d] += box[m][d];
                }
            }
            while (x[i][m] - shift >= box[m][m])
            {
                for (d = 0; d <= m; d++)
                {
                    x[i][d] -= box[m][d];
                }
            }
        }
    }
}

void put_atoms_in_compact_unitcell(PbcType pbcType, int ecenter, const matrix box, gmx::ArrayRef<gmx::RVec> x)
{
    t_pbc pbc;
    rvec  box_center, dx;

    set_pbc(&pbc, pbcType, box);

    if (pbc.pbcTypeDX == epbcdxUNSUPPORTED)
    {
        gmx_fatal(FARGS, "Can not put atoms in compact unitcell with unsupported PBC");
    }

    calc_box_center(ecenter, box, box_center);
    for (gmx::Index i = 0; (i < x.ssize()); ++i)
    {
        pbc_dx(&pbc, x[i], box_center, dx);
        rvec_add(box_center, dx, x[i]);
    }
}

/*! \brief Make molecules whole by shifting positions
 *
 * \param[in]     fplog     Log file
 * \param[in]     pbcType   The PBC type
 * \param[in]     correctVelocitiesForBoxDeformation  Whether to correct the velocities for
 *                                                    continuous box deformation
 * \param[in]     boxDeformation  The box deformation velocity
 * \param[in]     box       The simulation box
 * \param[in]     mtop      System topology definition
 * \param[in,out] x         The coordinates of the atoms
 * \param[in]     v         The velocities of the atoms
 * \param[in]     bFirst    Specifier for first-time PBC removal
 */
static void low_do_pbc_mtop(FILE*                    fplog,
                            PbcType                  pbcType,
                            const bool               correctVelocitiesForBoxDeformation,
                            const matrix             boxDeformation,
                            const matrix             box,
                            const gmx_mtop_t*        mtop,
                            gmx::ArrayRef<gmx::RVec> x,
                            gmx::ArrayRef<gmx::RVec> v,
                            gmx_bool                 bFirst)
{
    int as, mol;

    if (bFirst && fplog)
    {
        fprintf(fplog, "Removing pbc first time\n");
    }

    matrix boxDeformationRate;
    if (correctVelocitiesForBoxDeformation)
    {
        GMX_RELEASE_ASSERT(v.size() == x.size(), "Need velocities with box deformation");

        setBoxDeformationRate(boxDeformation, box, boxDeformationRate);
    }

    as = 0;
    for (const gmx_molblock_t& molb : mtop->molblock)
    {
        const gmx_moltype_t& moltype = mtop->moltype[molb.type];
        if (moltype.atoms.nr == 1 || (!bFirst && moltype.atoms.nr == 1))
        {
            /* Just one atom or charge group in the molecule, no PBC required */
            as += molb.nmol * moltype.atoms.nr;
        }
        else
        {
            t_graph graph = mk_graph_moltype(moltype);

            std::vector<gmx::RVec> xOrig(correctVelocitiesForBoxDeformation ? moltype.atoms.nr : 0);

            for (mol = 0; mol < molb.nmol; mol++)
            {
                auto xMol = x.subArray(as, moltype.atoms.nr);

                mk_mshift(fplog, &graph, pbcType, box, as_rvec_array(xMol.data()));

                if (correctVelocitiesForBoxDeformation)
                {
                    // Store a copy of the original coordinates, so we can compute displacements
                    std::copy(xMol.begin(), xMol.end(), xOrig.begin());
                }

                shift_self(graph, box, as_rvec_array(xMol.data()));
                /* The molecule is whole now.
                 * We don't need the second mk_mshift call as in do_pbc_first,
                 * since we no longer need this graph.
                 */

                if (correctVelocitiesForBoxDeformation)
                {
                    for (int i = 0; i < moltype.atoms.nr; i++)
                    {
                        correctVelocityForDisplacement<true>(
                                boxDeformationRate, v[as + i], xMol[i] - xOrig[i]);
                    }
                }

                as += moltype.atoms.nr;
            }
        }
    }
}

void do_pbc_first_mtop(FILE*                    fplog,
                       PbcType                  pbcType,
                       const bool               correctVelocitiesForBoxDeformation,
                       const matrix             boxDeformation,
                       const matrix             box,
                       const gmx_mtop_t*        mtop,
                       gmx::ArrayRef<gmx::RVec> x,
                       gmx::ArrayRef<gmx::RVec> v)
{
    low_do_pbc_mtop(
            fplog, pbcType, correctVelocitiesForBoxDeformation, boxDeformation, box, mtop, x, v, TRUE);
}

void do_pbc_mtop(PbcType pbcType, const matrix box, const gmx_mtop_t* mtop, rvec x[])
{
    low_do_pbc_mtop(nullptr,
                    pbcType,
                    false,
                    nullptr,
                    box,
                    mtop,
                    gmx::arrayRefFromArray<gmx::RVec>(reinterpret_cast<gmx::RVec*>(x), mtop->natoms),
                    {},
                    FALSE);
}

void setBoxDeformationRate(const matrix boxDeformation, const matrix box, matrix boxDeformationRate)
{
    clear_mat(boxDeformationRate);
    for (int d1 = 0; d1 < DIM; d1++)
    {
        for (int d2 = 0; d2 <= d1; d2++)
        {
            if (box[d1][d1] > 0)
            {
                boxDeformationRate[d1][d2] = boxDeformation[d1][d2] / box[d1][d1];
            }
        }
    }
}
