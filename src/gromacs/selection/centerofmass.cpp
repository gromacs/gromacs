/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2009- The GROMACS Authors
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
 * Implements functions in centerofmass.h.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_selection
 */
#include "gmxpre.h"

#include "centerofmass.h"

#include <cmath>

#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/block.h"
#include "gromacs/topology/mtop_lookup.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/real.h"

void gmx_calc_cog(const gmx_mtop_t* /* top */, rvec x[], int nrefat, const int index[], rvec xout)
{
    int m, ai;

    clear_rvec(xout);
    for (m = 0; m < nrefat; ++m)
    {
        ai = index[m];
        rvec_inc(xout, x[ai]);
    }
    svmul(1.0 / nrefat, xout, xout);
}

/*!
 * \param[in]  top    Topology structure with masses.
 * \param[in]  x      Position vectors of all atoms.
 * \param[in]  nrefat Number of atoms in the index.
 * \param[in]  index  Indices of atoms.
 * \param[out] xout   COM position for the indexed atoms.
 *
 * Works exactly as gmx_calc_cog() with the exception that a center of
 * mass are calculated, and hence a topology with masses is required.
 */
void gmx_calc_com(const gmx_mtop_t* top, rvec x[], int nrefat, const int index[], rvec xout)
{
    GMX_RELEASE_ASSERT(gmx_mtop_has_masses(top),
                       "No masses available while mass weighting was requested");
    clear_rvec(xout);
    real mtot = 0;
    int  molb = 0;
    for (int m = 0; m < nrefat; ++m)
    {
        const int  ai   = index[m];
        const real mass = mtopGetAtomMass(*top, ai, &molb);
        for (int j = 0; j < DIM; ++j)
        {
            xout[j] += mass * x[ai][j];
        }
        mtot += mass;
    }
    svmul(1.0 / mtot, xout, xout);
}

/*!
 * \param[in]  top    Topology structure with masses.
 * \param[in]  f      Forces on all atoms.
 * \param[in]  nrefat Number of atoms in the index.
 * \param[in]  index  Indices of atoms.
 * \param[out] fout   Force on the COG position for the indexed atoms.
 */
void gmx_calc_cog_f(const gmx_mtop_t* top, rvec f[], int nrefat, const int index[], rvec fout)
{
    GMX_RELEASE_ASSERT(gmx_mtop_has_masses(top),
                       "No masses available while mass weighting was requested");
    clear_rvec(fout);
    real mtot = 0;
    int  molb = 0;
    for (int m = 0; m < nrefat; ++m)
    {
        const int  ai   = index[m];
        const real mass = mtopGetAtomMass(*top, ai, &molb);
        for (int j = 0; j < DIM; ++j)
        {
            fout[j] += f[ai][j] / mass;
        }
        mtot += mass;
    }
    svmul(mtot / nrefat, fout, fout);
}

void gmx_calc_com_f(const gmx_mtop_t* /* top */, rvec f[], int nrefat, const int index[], rvec fout)
{
    clear_rvec(fout);
    for (int m = 0; m < nrefat; ++m)
    {
        const int ai = index[m];
        rvec_inc(fout, f[ai]);
    }
}

/*!
 * \param[in]  top   Topology structure with masses
 *   (can be NULL if \p bMASS==false).
 * \param[in]  x     Position vectors of all atoms.
 * \param[in]  nrefat Number of atoms in the index.
 * \param[in]  index Indices of atoms.
 * \param[in]  bMass If true, mass weighting is used.
 * \param[out] xout  COM/COG position for the indexed atoms.
 *
 * Calls either gmx_calc_com() or gmx_calc_cog() depending on the value of
 * \p bMass.
 * Other parameters are passed unmodified to these functions.
 */
void gmx_calc_comg(const gmx_mtop_t* top, rvec x[], int nrefat, const int index[], bool bMass, rvec xout)
{
    if (bMass)
    {
        gmx_calc_com(top, x, nrefat, index, xout);
    }
    else
    {
        gmx_calc_cog(top, x, nrefat, index, xout);
    }
}

/*!
 * \param[in]  top   Topology structure with masses
 *   (can be NULL if \p bMASS==true).
 * \param[in]  f     Forces on all atoms.
 * \param[in]  nrefat Number of atoms in the index.
 * \param[in]  index Indices of atoms.
 * \param[in]  bMass If true, force on COM is calculated.
 * \param[out] fout  Force on the COM/COG position for the indexed atoms.
 *
 * Calls either gmx_calc_cog_f() or gmx_calc_com_f() depending on the value of
 * \p bMass.
 * Other parameters are passed unmodified to these functions.
 */
void gmx_calc_comg_f(const gmx_mtop_t* top, rvec f[], int nrefat, const int index[], bool bMass, rvec fout)
{
    if (bMass)
    {
        gmx_calc_com_f(top, f, nrefat, index, fout);
    }
    else
    {
        gmx_calc_cog_f(top, f, nrefat, index, fout);
    }
}


/*!
 * \param[in]  top    Topology structure (unused, can be NULL).
 * \param[in]  x      Position vectors of all atoms.
 * \param[in]  pbc    Periodic boundary conditions structure.
 * \param[in]  nrefat Number of atoms in the index.
 * \param[in]  index  Indices of atoms.
 * \param[out] xout   COG position for the indexed atoms.
 *
 * Works exactly as gmx_calc_com_pbc(), but calculates the center of geometry.
 */
void gmx_calc_cog_pbc(const gmx_mtop_t* top, rvec x[], const t_pbc* pbc, int nrefat, const int index[], rvec xout)
{
    const real tol = 1e-4;
    bool       bChanged;
    int        m, j, ai;
    rvec       dx, xtest;

    /* First simple calculation */
    gmx_calc_cog(top, x, nrefat, index, xout);
    /* Now check if any atom is more than half the box from the COM */
    if (pbc)
    {
        do
        {
            bChanged = false;
            for (m = 0; m < nrefat; ++m)
            {
                ai = index[m];
                pbc_dx(pbc, x[ai], xout, dx);
                rvec_add(xout, dx, xtest);
                for (j = 0; j < DIM; ++j)
                {
                    if (std::fabs(xtest[j] - x[ai][j]) > tol)
                    {
                        /* Here we have used the wrong image for contributing to the COM */
                        xout[j] += (xtest[j] - x[ai][j]) / nrefat;
                        x[ai][j] = xtest[j];
                        bChanged = true;
                    }
                }
            }
        } while (bChanged);
    }
}

/*!
 * \param[in]  top    Topology structure with masses.
 * \param[in]  x      Position vectors of all atoms.
 * \param[in]  pbc    Periodic boundary conditions structure.
 * \param[in]  nrefat Number of atoms in the index.
 * \param[in]  index  Indices of atoms.
 * \param[out] xout   COM position for the indexed atoms.
 *
 * Works as gmx_calc_com(), but takes into account periodic boundary
 * conditions: If any atom is more than half the box from the COM,
 * it is wrapped around and a new COM is calculated. This is repeated
 * until no atoms violate the condition.
 *
 * Modified from src/tools/gmx_sorient.c in Gromacs distribution.
 */
void gmx_calc_com_pbc(const gmx_mtop_t* top, rvec x[], const t_pbc* pbc, int nrefat, const int index[], rvec xout)
{
    GMX_RELEASE_ASSERT(gmx_mtop_has_masses(top),
                       "No masses available while mass weighting was requested");
    /* First simple calculation */
    clear_rvec(xout);
    real mtot = 0;
    int  molb = 0;
    for (int m = 0; m < nrefat; ++m)
    {
        const int  ai   = index[m];
        const real mass = mtopGetAtomMass(*top, ai, &molb);
        for (int j = 0; j < DIM; ++j)
        {
            xout[j] += mass * x[ai][j];
        }
        mtot += mass;
    }
    svmul(1.0 / mtot, xout, xout);
    /* Now check if any atom is more than half the box from the COM */
    if (pbc)
    {
        const real tol = 1e-4;
        bool       bChanged;
        do
        {
            bChanged = false;
            molb     = 0;
            for (int m = 0; m < nrefat; ++m)
            {
                rvec       dx, xtest;
                const int  ai   = index[m];
                const real mass = mtopGetAtomMass(*top, ai, &molb) / mtot;
                pbc_dx(pbc, x[ai], xout, dx);
                rvec_add(xout, dx, xtest);
                for (int j = 0; j < DIM; ++j)
                {
                    if (std::fabs(xtest[j] - x[ai][j]) > tol)
                    {
                        /* Here we have used the wrong image for contributing to the COM */
                        xout[j] += mass * (xtest[j] - x[ai][j]);
                        x[ai][j] = xtest[j];
                        bChanged = true;
                    }
                }
            }
        } while (bChanged);
    }
}

/*!
 * \param[in]  top   Topology structure with masses
 *   (can be NULL if \p bMASS==false).
 * \param[in]  x     Position vectors of all atoms.
 * \param[in]  pbc    Periodic boundary conditions structure.
 * \param[in]  nrefat Number of atoms in the index.
 * \param[in]  index Indices of atoms.
 * \param[in]  bMass If true, mass weighting is used.
 * \param[out] xout  COM/COG position for the indexed atoms.
 *
 * Calls either gmx_calc_com() or gmx_calc_cog() depending on the value of
 * \p bMass.
 * Other parameters are passed unmodified to these functions.
 */
void gmx_calc_comg_pbc(const gmx_mtop_t* top, rvec x[], const t_pbc* pbc, int nrefat, const int index[], bool bMass, rvec xout)
{
    if (bMass)
    {
        gmx_calc_com_pbc(top, x, pbc, nrefat, index, xout);
    }
    else
    {
        gmx_calc_cog_pbc(top, x, pbc, nrefat, index, xout);
    }
}


void gmx_calc_cog_block(const gmx_mtop_t* /* top */, rvec x[], const t_block* block, const int index[], rvec xout[])
{
    int  b, i, ai;
    rvec xb;

    for (b = 0; b < block->nr; ++b)
    {
        clear_rvec(xb);
        for (i = block->index[b]; i < block->index[b + 1]; ++i)
        {
            ai = index[i];
            rvec_inc(xb, x[ai]);
        }
        svmul(1.0 / (block->index[b + 1] - block->index[b]), xb, xout[b]);
    }
}

/*!
 * \param[in]  top   Topology structure with masses.
 * \param[in]  x     Position vectors of all atoms.
 * \param[in]  block t_block structure that divides \p index into blocks.
 * \param[in]  index Indices of atoms.
 * \param[out] xout  \p block->nr COM positions.
 *
 * Works exactly as gmx_calc_cog_block() with the exception that centers of
 * mass are calculated, and hence a topology with masses is required.
 */
void gmx_calc_com_block(const gmx_mtop_t* top, rvec x[], const t_block* block, const int index[], rvec xout[])
{
    GMX_RELEASE_ASSERT(gmx_mtop_has_masses(top),
                       "No masses available while mass weighting was requested");
    int molb = 0;
    for (int b = 0; b < block->nr; ++b)
    {
        rvec xb;
        clear_rvec(xb);
        real mtot = 0;
        for (int i = block->index[b]; i < block->index[b + 1]; ++i)
        {
            const int  ai   = index[i];
            const real mass = mtopGetAtomMass(*top, ai, &molb);
            for (int d = 0; d < DIM; ++d)
            {
                xb[d] += mass * x[ai][d];
            }
            mtot += mass;
        }
        svmul(1.0 / mtot, xb, xout[b]);
    }
}

/*!
 * \param[in]  top   Topology structure with masses.
 * \param[in]  f     Forces on all atoms.
 * \param[in]  block t_block structure that divides \p index into blocks.
 * \param[in]  index Indices of atoms.
 * \param[out] fout  \p block->nr Forces on COG positions.
 */
void gmx_calc_cog_f_block(const gmx_mtop_t* top, rvec f[], const t_block* block, const int index[], rvec fout[])
{
    GMX_RELEASE_ASSERT(gmx_mtop_has_masses(top),
                       "No masses available while mass weighting was requested");
    int molb = 0;
    for (int b = 0; b < block->nr; ++b)
    {
        rvec fb;
        clear_rvec(fb);
        real mtot = 0;
        for (int i = block->index[b]; i < block->index[b + 1]; ++i)
        {
            const int  ai   = index[i];
            const real mass = mtopGetAtomMass(*top, ai, &molb);
            for (int d = 0; d < DIM; ++d)
            {
                fb[d] += f[ai][d] / mass;
            }
            mtot += mass;
        }
        svmul(mtot / (block->index[b + 1] - block->index[b]), fb, fout[b]);
    }
}

void gmx_calc_com_f_block(const gmx_mtop_t* /* top */,
                          rvec           f[],
                          const t_block* block,
                          const int      index[],
                          rvec           fout[])
{
    for (int b = 0; b < block->nr; ++b)
    {
        rvec fb;
        clear_rvec(fb);
        for (int i = block->index[b]; i < block->index[b + 1]; ++i)
        {
            const int ai = index[i];
            rvec_inc(fb, f[ai]);
        }
        copy_rvec(fb, fout[b]);
    }
}

/*!
 * \param[in]  top   Topology structure with masses
 *   (can be NULL if \p bMASS==false).
 * \param[in]  x     Position vectors of all atoms.
 * \param[in]  block t_block structure that divides \p index into blocks.
 * \param[in]  index Indices of atoms.
 * \param[in]  bMass If true, mass weighting is used.
 * \param[out] xout  \p block->nr COM/COG positions.
 *
 * Calls either gmx_calc_com_block() or gmx_calc_cog_block() depending on the
 * value of \p bMass.
 * Other parameters are passed unmodified to these functions.
 */
void gmx_calc_comg_block(const gmx_mtop_t* top,
                         rvec              x[],
                         const t_block*    block,
                         const int         index[],
                         bool              bMass,
                         rvec              xout[])
{
    if (bMass)
    {
        gmx_calc_com_block(top, x, block, index, xout);
    }
    else
    {
        gmx_calc_cog_block(top, x, block, index, xout);
    }
}

/*!
 * \param[in]  top   Topology structure with masses
 *   (can be NULL if \p bMASS==false).
 * \param[in]  f     Forces on all atoms.
 * \param[in]  block t_block structure that divides \p index into blocks.
 * \param[in]  index Indices of atoms.
 * \param[in]  bMass If true, force on COM is calculated.
 * \param[out] fout  \p block->nr forces on the COM/COG positions.
 *
 * Calls either gmx_calc_com_f_block() or gmx_calc_cog_f_block() depending on
 * the value of \p bMass.
 * Other parameters are passed unmodified to these functions.
 */
void gmx_calc_comg_f_block(const gmx_mtop_t* top,
                           rvec              f[],
                           const t_block*    block,
                           const int         index[],
                           bool              bMass,
                           rvec              fout[])
{
    if (bMass)
    {
        gmx_calc_com_f_block(top, f, block, index, fout);
    }
    else
    {
        gmx_calc_cog_f_block(top, f, block, index, fout);
    }
}

/*!
 * \param[in]  top   Topology structure with masses
 *   (can be NULL if \p bMASS==false).
 * \param[in]  x     Position vectors of all atoms.
 * \param[in]  block Blocks for calculation.
 * \param[in]  bMass If true, mass weighting is used.
 * \param[out] xout  \p block->nr COM/COG positions.
 *
 * Calls gmx_calc_comg_block(), converting the \p t_blocka structure into
 * a \p t_block and an index. Other parameters are passed unmodified.
 *
 * \attention
 * This function assumes that a pointer to \c t_blocka can be safely typecast
 * into \c t_block such that the index fields can still be referenced.
 * With the present Gromacs defitions of these types, this is the case,
 * but if the layout of these structures is changed, this may lead to strange
 * crashes.
 */
void gmx_calc_comg_blocka(const gmx_mtop_t* top, rvec x[], const t_blocka* block, bool bMass, rvec xout[])
{
    /* TODO: It would probably be better to do this without the type cast */
    gmx_calc_comg_block(top, x, reinterpret_cast<const t_block*>(block), block->a, bMass, xout);
}

/*!
 * \param[in]  top   Topology structure with masses
 *   (can be NULL if \p bMASS==true).
 * \param[in]  f     Forces on all atoms.
 * \param[in]  block Blocks for calculation.
 * \param[in]  bMass If true, force on COM is calculated.
 * \param[out] fout  \p block->nr forces on the COM/COG positions.
 *
 * Calls gmx_calc_comg_f_block(), converting the \p t_blocka structure into
 * a \p t_block and an index. Other parameters are passed unmodified.
 *
 * \attention
 * This function assumes that a pointer to \c t_blocka can be safely typecast
 * into \c t_block such that the index fields can still be referenced.
 * With the present Gromacs defitions of these types, this is the case,
 * but if the layout of these structures is changed, this may lead to strange
 * crashes.
 */
void gmx_calc_comg_f_blocka(const gmx_mtop_t* top, rvec f[], const t_blocka* block, bool bMass, rvec fout[])
{
    /* TODO: It would probably be better to do this without the type cast */
    gmx_calc_comg_f_block(top, f, reinterpret_cast<const t_block*>(block), block->a, bMass, fout);
}
