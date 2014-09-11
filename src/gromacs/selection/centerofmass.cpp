/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2009,2010,2011,2012,2013,2014, by the GROMACS development team, led by
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
/*! \internal \file
 * \brief
 * Implements functions in centerofmass.h.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_selection
 */
#include "gmxpre.h"

#include "centerofmass.h"

#include <errno.h>

#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/block.h"
#include "gromacs/topology/topology.h"

int
gmx_calc_cog(t_topology * /* top */, rvec x[], int nrefat, atom_id index[], rvec xout)
{
    int                 m, ai;

    clear_rvec(xout);
    for (m = 0; m < nrefat; ++m)
    {
        ai = index[m];
        rvec_inc(xout, x[ai]);
    }
    svmul(1.0/nrefat, xout, xout);
    return 0;
}

/*!
 * \param[in]  top    Topology structure with masses.
 * \param[in]  x      Position vectors of all atoms.
 * \param[in]  nrefat Number of atoms in the index.
 * \param[in]  index  Indices of atoms.
 * \param[out] xout   COM position for the indexed atoms.
 * \returns    0 on success, EINVAL if \p top is NULL.
 *
 * Works exactly as gmx_calc_cog() with the exception that a center of
 * mass are calculated, and hence a topology with masses is required.
 */
int
gmx_calc_com(t_topology *top, rvec x[], int nrefat, atom_id index[], rvec xout)
{
    int                 m, j, ai;
    real                mass, mtot;

    if (!top)
    {
        gmx_incons("no masses available while mass weighting was requested");
        return EINVAL;
    }
    clear_rvec(xout);
    mtot = 0;
    for (m = 0; m < nrefat; ++m)
    {
        ai   = index[m];
        mass = top->atoms.atom[ai].m;
        for (j = 0; j < DIM; ++j)
        {
            xout[j] += mass * x[ai][j];
        }
        mtot += mass;
    }
    svmul(1.0/mtot, xout, xout);
    return 0;
}

/*!
 * \param[in]  top    Topology structure with masses.
 * \param[in]  f      Forces on all atoms.
 * \param[in]  nrefat Number of atoms in the index.
 * \param[in]  index  Indices of atoms.
 * \param[out] fout   Force on the COG position for the indexed atoms.
 * \returns    0 on success, EINVAL if \p top is NULL.
 */
int
gmx_calc_cog_f(t_topology *top, rvec f[], int nrefat, atom_id index[], rvec fout)
{
    int                 m, j, ai;
    real                mass, mtot;

    if (!top)
    {
        gmx_incons("no masses available while mass weighting was needed");
        return EINVAL;
    }
    clear_rvec(fout);
    mtot = 0;
    for (m = 0; m < nrefat; ++m)
    {
        ai   = index[m];
        mass = top->atoms.atom[ai].m;
        for (j = 0; j < DIM; ++j)
        {
            fout[j] += f[ai][j] / mass;
        }
        mtot += mass;
    }
    svmul(mtot / nrefat, fout, fout);
    return 0;
}

int
gmx_calc_com_f(t_topology * /* top */, rvec f[], int nrefat, atom_id index[], rvec fout)
{
    clear_rvec(fout);
    for (int m = 0; m < nrefat; ++m)
    {
        const int ai = index[m];
        rvec_inc(fout, f[ai]);
    }
    return 0;
}

/*!
 * \param[in]  top   Topology structure with masses
 *   (can be NULL if \p bMASS==false).
 * \param[in]  x     Position vectors of all atoms.
 * \param[in]  nrefat Number of atoms in the index.
 * \param[in]  index Indices of atoms.
 * \param[in]  bMass If true, mass weighting is used.
 * \param[out] xout  COM/COG position for the indexed atoms.
 * \returns    0 on success, EINVAL if \p top is NULL and \p bMass is true.
 *
 * Calls either gmx_calc_com() or gmx_calc_cog() depending on the value of
 * \p bMass.
 * Other parameters are passed unmodified to these functions.
 */
int
gmx_calc_comg(t_topology *top, rvec x[], int nrefat, atom_id index[],
              bool bMass, rvec xout)
{
    if (bMass)
    {
        return gmx_calc_com(top, x, nrefat, index, xout);
    }
    else
    {
        return gmx_calc_cog(top, x, nrefat, index, xout);
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
 * \returns    0 on success, EINVAL if \p top is NULL and \p bMass is false.
 *
 * Calls either gmx_calc_cog_f() or gmx_calc_com_f() depending on the value of
 * \p bMass.
 * Other parameters are passed unmodified to these functions.
 */
int
gmx_calc_comg_f(t_topology *top, rvec f[], int nrefat, atom_id index[],
                bool bMass, rvec fout)
{
    if (bMass)
    {
        return gmx_calc_com_f(top, f, nrefat, index, fout);
    }
    else
    {
        return gmx_calc_cog_f(top, f, nrefat, index, fout);
    }
}


/*!
 * \param[in]  top    Topology structure (unused, can be NULL).
 * \param[in]  x      Position vectors of all atoms.
 * \param[in]  pbc    Periodic boundary conditions structure.
 * \param[in]  nrefat Number of atoms in the index.
 * \param[in]  index  Indices of atoms.
 * \param[out] xout   COG position for the indexed atoms.
 * \returns    0 on success.
 *
 * Works exactly as gmx_calc_com_pbc(), but calculates the center of geometry.
 */
int
gmx_calc_cog_pbc(t_topology *top, rvec x[], t_pbc *pbc,
                 int nrefat, atom_id index[], rvec xout)
{
    const real          tol = 1e-4;
    bool                bChanged;
    int                 m, j, ai, iter;
    rvec                dx, xtest;

    /* First simple calculation */
    gmx_calc_cog(top, x, nrefat, index, xout);
    /* Now check if any atom is more than half the box from the COM */
    if (pbc)
    {
        iter = 0;
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
                    if (fabs(xtest[j] - x[ai][j]) > tol)
                    {
                        /* Here we have used the wrong image for contributing to the COM */
                        xout[j] += (xtest[j] - x[ai][j]) / nrefat;
                        x[ai][j] = xtest[j];
                        bChanged = true;
                    }
                }
            }
            iter++;
        }
        while (bChanged);
    }
    return 0;
}

/*!
 * \param[in]  top    Topology structure with masses.
 * \param[in]  x      Position vectors of all atoms.
 * \param[in]  pbc    Periodic boundary conditions structure.
 * \param[in]  nrefat Number of atoms in the index.
 * \param[in]  index  Indices of atoms.
 * \param[out] xout   COM position for the indexed atoms.
 * \returns    0 on success, EINVAL if \p top is NULL.
 *
 * Works as gmx_calc_com(), but takes into account periodic boundary
 * conditions: If any atom is more than half the box from the COM,
 * it is wrapped around and a new COM is calculated. This is repeated
 * until no atoms violate the condition.
 *
 * Modified from src/tools/gmx_sorient.c in Gromacs distribution.
 */
int
gmx_calc_com_pbc(t_topology *top, rvec x[], t_pbc *pbc,
                 int nrefat, atom_id index[], rvec xout)
{
    const real          tol = 1e-4;
    bool                bChanged;
    int                 m, j, ai, iter;
    real                mass, mtot;
    rvec                dx, xtest;

    if (!top)
    {
        gmx_incons("no masses available while mass weighting was requested");
        return EINVAL;
    }
    /* First simple calculation */
    clear_rvec(xout);
    mtot = 0;
    for (m = 0; m < nrefat; ++m)
    {
        ai   = index[m];
        mass = top->atoms.atom[ai].m;
        for (j = 0; j < DIM; ++j)
        {
            xout[j] += mass * x[ai][j];
        }
        mtot += mass;
    }
    svmul(1.0/mtot, xout, xout);
    /* Now check if any atom is more than half the box from the COM */
    if (pbc)
    {
        iter = 0;
        do
        {
            bChanged = false;
            for (m = 0; m < nrefat; ++m)
            {
                ai   = index[m];
                mass = top->atoms.atom[ai].m / mtot;
                pbc_dx(pbc, x[ai], xout, dx);
                rvec_add(xout, dx, xtest);
                for (j = 0; j < DIM; ++j)
                {
                    if (fabs(xtest[j] - x[ai][j]) > tol)
                    {
                        /* Here we have used the wrong image for contributing to the COM */
                        xout[j] += mass * (xtest[j] - x[ai][j]);
                        x[ai][j] = xtest[j];
                        bChanged = true;
                    }
                }
            }
            iter++;
        }
        while (bChanged);
    }
    return 0;
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
 * \returns    0 on success, EINVAL if \p top is NULL and \p bMass is true.
 *
 * Calls either gmx_calc_com() or gmx_calc_cog() depending on the value of
 * \p bMass.
 * Other parameters are passed unmodified to these functions.
 */
int
gmx_calc_comg_pbc(t_topology *top, rvec x[], t_pbc *pbc,
                  int nrefat, atom_id index[], bool bMass, rvec xout)
{
    if (bMass)
    {
        return gmx_calc_com_pbc(top, x, pbc, nrefat, index, xout);
    }
    else
    {
        return gmx_calc_cog_pbc(top, x, pbc, nrefat, index, xout);
    }
}


int
gmx_calc_cog_block(t_topology * /* top */, rvec x[], t_block *block, atom_id index[],
                   rvec xout[])
{
    int                 b, i, ai;
    rvec                xb;

    for (b = 0; b < block->nr; ++b)
    {
        clear_rvec(xb);
        for (i = block->index[b]; i < block->index[b+1]; ++i)
        {
            ai = index[i];
            rvec_inc(xb, x[ai]);
        }
        svmul(1.0/(block->index[b+1] - block->index[b]), xb, xout[b]);
    }
    return 0;
}

/*!
 * \param[in]  top   Topology structure with masses.
 * \param[in]  x     Position vectors of all atoms.
 * \param[in]  block t_block structure that divides \p index into blocks.
 * \param[in]  index Indices of atoms.
 * \param[out] xout  \p block->nr COM positions.
 * \returns    0 on success, EINVAL if \p top is NULL.
 *
 * Works exactly as gmx_calc_cog_block() with the exception that centers of
 * mass are calculated, and hence a topology with masses is required.
 */
int
gmx_calc_com_block(t_topology *top, rvec x[], t_block *block, atom_id index[],
                   rvec xout[])
{
    int                 b, i, ai, d;
    rvec                xb;
    real                mass, mtot;

    if (!top)
    {
        gmx_incons("no masses available while mass weighting was requested");
        return EINVAL;
    }
    for (b = 0; b < block->nr; ++b)
    {
        clear_rvec(xb);
        mtot = 0;
        for (i = block->index[b]; i < block->index[b+1]; ++i)
        {
            ai   = index[i];
            mass = top->atoms.atom[ai].m;
            for (d = 0; d < DIM; ++d)
            {
                xb[d] += mass * x[ai][d];
            }
            mtot += mass;
        }
        svmul(1.0/mtot, xb, xout[b]);
    }
    return 0;
}

/*!
 * \param[in]  top   Topology structure with masses.
 * \param[in]  f     Forces on all atoms.
 * \param[in]  block t_block structure that divides \p index into blocks.
 * \param[in]  index Indices of atoms.
 * \param[out] fout  \p block->nr Forces on COG positions.
 * \returns    0 on success, EINVAL if \p top is NULL.
 */
int
gmx_calc_cog_f_block(t_topology *top, rvec f[], t_block *block, atom_id index[],
                     rvec fout[])
{
    int                 b, i, ai, d;
    rvec                fb;
    real                mass, mtot;

    if (!top)
    {
        gmx_incons("no masses available while mass weighting was needed");
        return EINVAL;
    }
    for (b = 0; b < block->nr; ++b)
    {
        clear_rvec(fb);
        mtot = 0;
        for (i = block->index[b]; i < block->index[b+1]; ++i)
        {
            ai   = index[i];
            mass = top->atoms.atom[ai].m;
            for (d = 0; d < DIM; ++d)
            {
                fb[d] += f[ai][d] / mass;
            }
            mtot += mass;
        }
        svmul(mtot / (block->index[b+1] - block->index[b]), fb, fout[b]);
    }
    return 0;
}

int
gmx_calc_com_f_block(t_topology * /* top */, rvec f[], t_block *block, atom_id index[],
                     rvec fout[])
{
    for (int b = 0; b < block->nr; ++b)
    {
        rvec fb;
        clear_rvec(fb);
        for (int i = block->index[b]; i < block->index[b+1]; ++i)
        {
            const int ai = index[i];
            rvec_inc(fb, f[ai]);
        }
        copy_rvec(fb, fout[b]);
    }
    return 0;
}

/*!
 * \param[in]  top   Topology structure with masses
 *   (can be NULL if \p bMASS==false).
 * \param[in]  x     Position vectors of all atoms.
 * \param[in]  block t_block structure that divides \p index into blocks.
 * \param[in]  index Indices of atoms.
 * \param[in]  bMass If true, mass weighting is used.
 * \param[out] xout  \p block->nr COM/COG positions.
 * \returns    0 on success, EINVAL if \p top is NULL and \p bMass is true.
 *
 * Calls either gmx_calc_com_block() or gmx_calc_cog_block() depending on the
 * value of \p bMass.
 * Other parameters are passed unmodified to these functions.
 */
int
gmx_calc_comg_block(t_topology *top, rvec x[], t_block *block, atom_id index[],
                    bool bMass, rvec xout[])
{
    if (bMass)
    {
        return gmx_calc_com_block(top, x, block, index, xout);
    }
    else
    {
        return gmx_calc_cog_block(top, x, block, index, xout);
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
 * \returns    0 on success, EINVAL if \p top is NULL and \p bMass is true.
 *
 * Calls either gmx_calc_com_f_block() or gmx_calc_cog_f_block() depending on
 * the value of \p bMass.
 * Other parameters are passed unmodified to these functions.
 */
int
gmx_calc_comg_f_block(t_topology *top, rvec f[], t_block *block, atom_id index[],
                      bool bMass, rvec fout[])
{
    if (bMass)
    {
        return gmx_calc_com_f_block(top, f, block, index, fout);
    }
    else
    {
        return gmx_calc_cog_f_block(top, f, block, index, fout);
    }
}

/*!
 * \param[in]  top   Topology structure with masses
 *   (can be NULL if \p bMASS==false).
 * \param[in]  x     Position vectors of all atoms.
 * \param[in]  block Blocks for calculation.
 * \param[in]  bMass If true, mass weighting is used.
 * \param[out] xout  \p block->nr COM/COG positions.
 * \returns    0 on success, EINVAL if \p top is NULL and \p bMass is true.
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
int
gmx_calc_comg_blocka(t_topology *top, rvec x[], t_blocka *block,
                     bool bMass, rvec xout[])
{
    /* TODO: It would probably be better to do this without the type cast */
    return gmx_calc_comg_block(top, x, (t_block *)block, block->a, bMass, xout);
}

/*!
 * \param[in]  top   Topology structure with masses
 *   (can be NULL if \p bMASS==true).
 * \param[in]  f     Forces on all atoms.
 * \param[in]  block Blocks for calculation.
 * \param[in]  bMass If true, force on COM is calculated.
 * \param[out] fout  \p block->nr forces on the COM/COG positions.
 * \returns    0 on success, EINVAL if \p top is NULL and \p bMass is false.
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
int
gmx_calc_comg_f_blocka(t_topology *top, rvec f[], t_blocka *block,
                       bool bMass, rvec fout[])
{
    /* TODO: It would probably be better to do this without the type cast */
    return gmx_calc_comg_f_block(top, f, (t_block *)block, block->a, bMass, fout);
}
