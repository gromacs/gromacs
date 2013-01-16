/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
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
/*! \page nbsearch Neighborhood search routines
 *
 * Functions to find particles within a neighborhood of a set of particles
 * are defined in nbsearch.h.
 * The usage is simple: a data structure is allocated with
 * gmx_ana_nbsearch_create(), and the box shape and reference positions for a
 * frame are set using gmx_ana_nbsearch_init() or gmx_ana_nbsearch_pos_init().
 * Searches can then be performed with gmx_ana_nbsearch_is_within() and
 * gmx_ana_nbsearch_mindist(), or with versions that take the \c gmx_ana_pos_t
 * data structure.
 * When the data structure is no longer required, it can be freed with
 * gmx_ana_nbsearch_free().
 *
 * \internal
 *
 * \todo
 * The grid implementation could still be optimized in several different ways:
 *   - Triclinic grid cells are not the most efficient shape, but make PBC
 *     handling easier.
 *   - Precalculating the required PBC shift for a pair of cells outside the
 *     inner loop. After this is done, it should be quite straightforward to
 *     move to rectangular cells.
 *   - Pruning grid cells from the search list if they are completely outside
 *     the sphere that is being considered.
 *   - A better heuristic could be added for falling back to simple loops for a
 *     small number of reference particles.
 *   - A better heuristic for selecting the grid size.
 *   - A multi-level grid implementation could be used to be able to use small
 *     grids for short cutoffs with very inhomogeneous particle distributions
 *     without a memory cost.
 */
/*! \internal \file
 * \brief Implementation of functions in nbsearch.h.
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>

#include <smalloc.h>
#include <typedefs.h>
#include <pbc.h>
#include <vec.h>

#include <nbsearch.h>
#include <position.h>

/*! \internal \brief
 * Data structure for neighborhood searches.
 */
struct gmx_ana_nbsearch_t
{
    /** The cutoff. */
    real           cutoff;
    /** The cutoff squared. */
    real           cutoff2;
    /** Maximum number of reference points. */
    int            maxnref;

    /** Number of reference points for the current frame. */
    int            nref;
    /** Reference point positions. */
    rvec          *xref;
    /** Reference position ids (NULL if not available). */
    int           *refid;
    /** PBC data. */
    t_pbc         *pbc;

    /** Number of excluded reference positions for current test particle. */
    int            nexcl;
    /** Exclusions for current test particle. */
    int           *excl;

    /** Whether to try grid searching. */
    gmx_bool           bTryGrid;
    /** Whether grid searching is actually used for the current positions. */
    gmx_bool           bGrid;
    /** Array allocated for storing in-unit-cell reference positions. */
    rvec              *xref_alloc;
    /** FALSE if the box is rectangular. */
    gmx_bool           bTric;
    /** Box vectors of a single grid cell. */
    matrix             cellbox;
    /** The reciprocal cell vectors as columns; the inverse of \p cellbox. */
    matrix             recipcell;
    /** Number of cells along each dimension. */
    ivec               ncelldim;
    /** Total number of cells. */
    int                ncells;
    /** Number of reference positions in each cell. */
    int               *ncatoms;
    /** List of reference positions in each cell. */
    atom_id          **catom;
    /** Allocation counts for each \p catom[i]. */
    int               *catom_nalloc;
    /** Allocation count for the per-cell arrays. */
    int                cells_nalloc;
    /** Number of neighboring cells to consider. */
    int                ngridnb;
    /** Offsets of the neighboring cells to consider. */
    ivec              *gnboffs;
    /** Allocation count for \p gnboffs. */
    int                gnboffs_nalloc;

    /** Stores test position during a pair loop. */
    rvec           xtest;
    /** Stores the previous returned position during a pair loop. */
    int            previ;
    /** Stores the current exclusion index during loops. */
    int            exclind;
    /** Stores the test particle cell index during loops. */
    ivec           testcell;
    /** Stores the current cell neighbor index during pair loops. */
    int            prevnbi;
    /** Stores the index within the current cell during pair loops. */
    int            prevcai;
};

/*!
 * \param[out] data   Neighborhood search data structure pointer to initialize.
 * \param[in]  cutoff Cutoff distance for the search
 *   (<=0 stands for no cutoff).
 * \param[in]  maxn   Maximum number of reference particles.
 * \returns    0 on success.
 */
int
gmx_ana_nbsearch_create(gmx_ana_nbsearch_t **data, real cutoff, int maxn)
{
    gmx_ana_nbsearch_t *d;

    snew(d, 1);
    d->bTryGrid = TRUE;
    if (cutoff <= 0)
    {
        cutoff      = GMX_REAL_MAX;
        d->bTryGrid = FALSE;
    }
    d->cutoff  = cutoff;
    d->cutoff2 = sqr(cutoff);
    d->maxnref = maxn;

    d->xref    = NULL;
    d->nexcl   = 0;
    d->exclind = 0;

    d->xref_alloc   = NULL;
    d->ncells       = 0;
    d->ncatoms      = NULL;
    d->catom        = NULL;
    d->catom_nalloc = 0;
    d->cells_nalloc = 0;

    d->ngridnb        = 0;
    d->gnboffs        = NULL;
    d->gnboffs_nalloc = 0;

    *data = d;
    return 0;
}

/*!
 * \param     d Data structure to free.
 *
 * After the call, the pointer \p d is no longer valid.
 */
void
gmx_ana_nbsearch_free(gmx_ana_nbsearch_t *d)
{
    sfree(d->xref_alloc);
    sfree(d->ncatoms);
    if (d->catom)
    {
        int ci;

        for (ci = 0; ci < d->ncells; ++ci)
        {
            sfree(d->catom[ci]);
        }
        sfree(d->catom);
    }
    sfree(d->catom_nalloc);
    sfree(d->gnboffs);
    sfree(d);
}

/*! \brief
 * Calculates offsets to neighboring grid cells that should be considered.
 *
 * \param[in,out] d    Grid information.
 * \param[in]     pbc  Information about the box.
 */
static void
grid_init_cell_nblist(gmx_ana_nbsearch_t *d, t_pbc *pbc)
{
    int   maxx, maxy, maxz;
    int   x, y, z, i;
    real  rvnorm;

    /* Find the extent of the sphere in triclinic coordinates */
    maxz   = (int)(d->cutoff * d->recipcell[ZZ][ZZ]) + 1;
    rvnorm = sqrt(sqr(d->recipcell[YY][YY]) + sqr(d->recipcell[ZZ][YY]));
    maxy   = (int)(d->cutoff * rvnorm) + 1;
    rvnorm = sqrt(sqr(d->recipcell[XX][XX]) + sqr(d->recipcell[YY][XX])
                  + sqr(d->recipcell[ZZ][XX]));
    maxx = (int)(d->cutoff * rvnorm) + 1;

    /* Calculate the number of cells and reallocate if necessary */
    d->ngridnb = (2 * maxx + 1) * (2 * maxy + 1) * (2 * maxz + 1);
    if (d->gnboffs_nalloc < d->ngridnb)
    {
        d->gnboffs_nalloc = d->ngridnb;
        srenew(d->gnboffs, d->gnboffs_nalloc);
    }

    /* Store the whole cube */
    /* TODO: Prune off corners that are not needed */
    i = 0;
    for (x = -maxx; x <= maxx; ++x)
    {
        for (y = -maxy; y <= maxy; ++y)
        {
            for (z = -maxz; z <= maxz; ++z)
            {
                d->gnboffs[i][XX] = x;
                d->gnboffs[i][YY] = y;
                d->gnboffs[i][ZZ] = z;
                ++i;
            }
        }
    }
}

/*! \brief
 * Determines a suitable grid size.
 *
 * \param[in,out] d    Grid information.
 * \param[in]     pbc  Information about the box.
 * \returns  FALSE if grid search is not suitable.
 */
static gmx_bool
grid_setup_cells(gmx_ana_nbsearch_t *d, t_pbc *pbc)
{
    real targetsize;
    int  dd;

#ifdef HAVE_CBRT
    targetsize = cbrt(pbc->box[XX][XX] * pbc->box[YY][YY] * pbc->box[ZZ][ZZ]
                      * 10 / d->nref);
#else
    targetsize = pow(pbc->box[XX][XX] * pbc->box[YY][YY] * pbc->box[ZZ][ZZ]
                     * 10 / d->nref, 1./3.);
#endif

    d->ncells = 1;
    for (dd = 0; dd < DIM; ++dd)
    {
        d->ncelldim[dd] = (int)(pbc->box[dd][dd] / targetsize);
        d->ncells      *= d->ncelldim[dd];
        if (d->ncelldim[dd] < 3)
        {
            return FALSE;
        }
    }
    /* Reallocate if necessary */
    if (d->cells_nalloc < d->ncells)
    {
        int  i;

        srenew(d->ncatoms, d->ncells);
        srenew(d->catom, d->ncells);
        srenew(d->catom_nalloc, d->ncells);
        for (i = d->cells_nalloc; i < d->ncells; ++i)
        {
            d->catom[i]        = NULL;
            d->catom_nalloc[i] = 0;
        }
        d->cells_nalloc = d->ncells;
    }
    return TRUE;
}

/*! \brief
 * Sets ua a search grid for a given box.
 *
 * \param[in,out] d    Grid information.
 * \param[in]     pbc  Information about the box.
 * \returns  FALSE if grid search is not suitable.
 */
static gmx_bool
grid_set_box(gmx_ana_nbsearch_t *d, t_pbc *pbc)
{
    int dd;

    /* TODO: This check could be improved. */
    if (0.5*pbc->max_cutoff2 < d->cutoff2)
    {
        return FALSE;
    }

    if (!grid_setup_cells(d, pbc))
    {
        return FALSE;
    }

    d->bTric = TRICLINIC(pbc->box);
    if (d->bTric)
    {
        for (dd = 0; dd < DIM; ++dd)
        {
            svmul(1.0 / d->ncelldim[dd], pbc->box[dd], d->cellbox[dd]);
        }
        m_inv_ur0(d->cellbox, d->recipcell);
    }
    else
    {
        for (dd = 0; dd < DIM; ++dd)
        {
            d->cellbox[dd][dd]   = pbc->box[dd][dd] / d->ncelldim[dd];
            d->recipcell[dd][dd] = 1 / d->cellbox[dd][dd];
        }
    }
    grid_init_cell_nblist(d, pbc);
    return TRUE;
}

/*! \brief
 * Maps a point into a grid cell.
 *
 * \param[in]  d    Grid information.
 * \param[in]  x    Point to map.
 * \param[out] cell Indices of the grid cell in which \p x lies.
 *
 * \p x should be in the triclinic unit cell.
 */
static void
grid_map_onto(gmx_ana_nbsearch_t *d, const rvec x, ivec cell)
{
    int dd;

    if (d->bTric)
    {
        rvec tx;

        tmvmul_ur0(d->recipcell, x, tx);
        for (dd = 0; dd < DIM; ++dd)
        {
            cell[dd] = (int)tx[dd];
        }
    }
    else
    {
        for (dd = 0; dd < DIM; ++dd)
        {
            cell[dd] = (int)(x[dd] * d->recipcell[dd][dd]);
        }
    }
}

/*! \brief
 * Calculates linear index of a grid cell.
 *
 * \param[in]  d    Grid information.
 * \param[in]  cell Cell indices.
 * \returns    Linear index of \p cell.
 */
static int
grid_index(gmx_ana_nbsearch_t *d, const ivec cell)
{
    return cell[XX] + cell[YY] * d->ncelldim[XX]
           + cell[ZZ] * d->ncelldim[XX] * d->ncelldim[YY];
}

/*! \brief
 * Clears all grid cells.
 *
 * \param[in,out] d    Grid information.
 */
static void
grid_clear_cells(gmx_ana_nbsearch_t *d)
{
    int  ci;

    for (ci = 0; ci < d->ncells; ++ci)
    {
        d->ncatoms[ci] = 0;
    }
}

/*! \brief
 * Adds an index into a grid cell.
 *
 * \param[in,out] d    Grid information.
 * \param[in]     cell Cell into which \p i should be added.
 * \param[in]     i    Index to add.
 */
static void
grid_add_to_cell(gmx_ana_nbsearch_t *d, const ivec cell, int i)
{
    int ci = grid_index(d, cell);

    if (d->ncatoms[ci] == d->catom_nalloc[ci])
    {
        d->catom_nalloc[ci] += 10;
        srenew(d->catom[ci], d->catom_nalloc[ci]);
    }
    d->catom[ci][d->ncatoms[ci]++] = i;
}

/*!
 * \param[in,out] d   Neighborhood search data structure.
 * \param[in]     pbc PBC information for the frame.
 * \param[in]     n   Number of reference positions for the frame.
 * \param[in]     x   \p n reference positions for the frame.
 * \returns       0 on success.
 *
 * Initializes the data structure \p d such that it can be used to search
 * for the neighbors of \p x.
 */
int
gmx_ana_nbsearch_init(gmx_ana_nbsearch_t *d, t_pbc *pbc, int n, rvec x[])
{
    d->pbc  = pbc;
    d->nref = n;
    if (!pbc)
    {
        d->bGrid = FALSE;
    }
    else if (d->bTryGrid)
    {
        d->bGrid = grid_set_box(d, pbc);
    }
    if (d->bGrid)
    {
        int  i;

        if (!d->xref_alloc)
        {
            snew(d->xref_alloc, d->maxnref);
        }
        d->xref = d->xref_alloc;
        grid_clear_cells(d);

        for (i = 0; i < n; ++i)
        {
            copy_rvec(x[i], d->xref[i]);
        }
        put_atoms_in_triclinic_unitcell(ecenterTRIC, pbc->box, n, d->xref);
        for (i = 0; i < n; ++i)
        {
            ivec refcell;

            grid_map_onto(d, d->xref[i], refcell);
            grid_add_to_cell(d, refcell, i);
        }
    }
    else
    {
        d->xref = x;
    }
    d->refid = NULL;
    return 0;
}

/*!
 * \param[in,out] d   Neighborhood search data structure.
 * \param[in]     pbc PBC information for the frame.
 * \param[in]     p   Reference positions for the frame.
 * \returns       0 on success.
 *
 * A convenience wrapper for gmx_ana_nbsearch_init().
 */
int
gmx_ana_nbsearch_pos_init(gmx_ana_nbsearch_t *d, t_pbc *pbc, gmx_ana_pos_t *p)
{
    int rc;

    rc       = gmx_ana_nbsearch_init(d, pbc, p->nr, p->x);
    d->refid = (p->nr < d->maxnref ? p->m.refid : NULL);
    return rc;
}

/*!
 * \param[in,out] d     Neighborhood search data structure.
 * \param[in]     nexcl Number of reference positions to exclude from next
 *      search.
 * \param[in]     excl  Indices of reference positions to exclude.
 * \returns       0 on success.
 *
 * The set exclusions remain in effect until the next call of this function.
 */
int
gmx_ana_nbsearch_set_excl(gmx_ana_nbsearch_t *d, int nexcl, int excl[])
{

    d->nexcl = nexcl;
    d->excl  = excl;
    return 0;
}

/*! \brief
 * Helper function to check whether a reference point should be excluded.
 */
static gmx_bool
is_excluded(gmx_ana_nbsearch_t *d, int j)
{
    if (d->exclind < d->nexcl)
    {
        if (d->refid)
        {
            while (d->exclind < d->nexcl &&d->refid[j] > d->excl[d->exclind])
            {
                ++d->exclind;
            }
            if (d->exclind < d->nexcl && d->refid[j] == d->excl[d->exclind])
            {
                ++d->exclind;
                return TRUE;
            }
        }
        else
        {
            while (d->bGrid && d->exclind < d->nexcl && d->excl[d->exclind] < j)
            {
                ++d->exclind;
            }
            if (d->excl[d->exclind] == j)
            {
                ++d->exclind;
                return TRUE;
            }
        }
    }
    return FALSE;
}

/*! \brief
 * Initializes a grid search to find reference positions neighboring \p x.
 */
static void
grid_search_start(gmx_ana_nbsearch_t *d, rvec x)
{
    copy_rvec(x, d->xtest);
    if (d->bGrid)
    {
        put_atoms_in_triclinic_unitcell(ecenterTRIC, d->pbc->box, 1, &d->xtest);
        grid_map_onto(d, d->xtest, d->testcell);
        d->prevnbi = 0;
        d->prevcai = -1;
    }
    else
    {
        d->previ = -1;
    }
    d->exclind = 0;
}

/*! \brief
 * Does a grid search.
 */
static gmx_bool
grid_search(gmx_ana_nbsearch_t *d,
            gmx_bool (*action)(gmx_ana_nbsearch_t *d, int i, real r2))
{
    int  i;
    rvec dx;
    real r2;

    if (d->bGrid)
    {
        int  nbi, ci, cai;

        nbi = d->prevnbi;
        cai = d->prevcai + 1;

        for (; nbi < d->ngridnb; ++nbi)
        {
            ivec cell;

            ivec_add(d->testcell, d->gnboffs[nbi], cell);
            /* TODO: Support for 2D and screw PBC */
            cell[XX] = (cell[XX] + d->ncelldim[XX]) % d->ncelldim[XX];
            cell[YY] = (cell[YY] + d->ncelldim[YY]) % d->ncelldim[YY];
            cell[ZZ] = (cell[ZZ] + d->ncelldim[ZZ]) % d->ncelldim[ZZ];
            ci       = grid_index(d, cell);
            /* TODO: Calculate the required PBC shift outside the inner loop */
            for (; cai < d->ncatoms[ci]; ++cai)
            {
                i = d->catom[ci][cai];
                if (is_excluded(d, i))
                {
                    continue;
                }
                pbc_dx_aiuc(d->pbc, d->xtest, d->xref[i], dx);
                r2 = norm2(dx);
                if (r2 <= d->cutoff2)
                {
                    if (action(d, i, r2))
                    {
                        d->prevnbi = nbi;
                        d->prevcai = cai;
                        d->previ   = i;
                        return TRUE;
                    }
                }
            }
            d->exclind = 0;
            cai        = 0;
        }
    }
    else
    {
        i = d->previ + 1;
        for (; i < d->nref; ++i)
        {
            if (is_excluded(d, i))
            {
                continue;
            }
            if (d->pbc)
            {
                pbc_dx(d->pbc, d->xtest, d->xref[i], dx);
            }
            else
            {
                rvec_sub(d->xtest, d->xref[i], dx);
            }
            r2 = norm2(dx);
            if (r2 <= d->cutoff2)
            {
                if (action(d, i, r2))
                {
                    d->previ = i;
                    return TRUE;
                }
            }
        }
    }
    return FALSE;
}

/*! \brief
 * Helper function to use with grid_search() to find the next neighbor.
 *
 * Simply breaks the loop on the first found neighbor.
 */
static gmx_bool
within_action(gmx_ana_nbsearch_t *d, int i, real r2)
{
    return TRUE;
}

/*! \brief
 * Helper function to use with grid_search() to find the minimum distance.
 */
static gmx_bool
mindist_action(gmx_ana_nbsearch_t *d, int i, real r2)
{
    d->cutoff2 = r2;
    return FALSE;
}

/*!
 * \param[in] d   Neighborhood search data structure.
 * \param[in] x   Test position.
 * \returns   TRUE if \p x is within the cutoff of any reference position,
 *   FALSE otherwise.
 */
gmx_bool
gmx_ana_nbsearch_is_within(gmx_ana_nbsearch_t *d, rvec x)
{
    grid_search_start(d, x);
    return grid_search(d, &within_action);
}

/*!
 * \param[in] d   Neighborhood search data structure.
 * \param[in] p   Test positions.
 * \param[in] i   Use the i'th position in \p p for testing.
 * \returns   TRUE if the test position is within the cutoff of any reference
 *   position, FALSE otherwise.
 */
gmx_bool
gmx_ana_nbsearch_pos_is_within(gmx_ana_nbsearch_t *d, gmx_ana_pos_t *p, int i)
{
    return gmx_ana_nbsearch_is_within(d, p->x[i]);
}

/*!
 * \param[in] d   Neighborhood search data structure.
 * \param[in] x   Test position.
 * \returns   The distance to the nearest reference position, or the cutoff
 *   value if there are no reference positions within the cutoff.
 */
real
gmx_ana_nbsearch_mindist(gmx_ana_nbsearch_t *d, rvec x)
{
    real mind;

    grid_search_start(d, x);
    grid_search(d, &mindist_action);
    mind       = sqrt(d->cutoff2);
    d->cutoff2 = sqr(d->cutoff);
    return mind;
}

/*!
 * \param[in] d   Neighborhood search data structure.
 * \param[in] p   Test positions.
 * \param[in] i   Use the i'th position in \p p for testing.
 * \returns   The distance to the nearest reference position, or the cutoff
 *   value if there are no reference positions within the cutoff.
 */
real
gmx_ana_nbsearch_pos_mindist(gmx_ana_nbsearch_t *d, gmx_ana_pos_t *p, int i)
{
    return gmx_ana_nbsearch_mindist(d, p->x[i]);
}

/*!
 * \param[in]  d   Neighborhood search data structure.
 * \param[in]  n   Number of test positions in \p x.
 * \param[in]  x   Test positions.
 * \param[out] jp  Index of the reference position in the first pair.
 * \returns    TRUE if there are positions within the cutoff.
 */
gmx_bool
gmx_ana_nbsearch_first_within(gmx_ana_nbsearch_t *d, rvec x, int *jp)
{
    grid_search_start(d, x);
    return gmx_ana_nbsearch_next_within(d, jp);
}

/*!
 * \param[in]  d   Neighborhood search data structure.
 * \param[in]  p   Test positions.
 * \param[in]  i   Use the i'th position in \p p.
 * \param[out] jp  Index of the reference position in the first pair.
 * \returns    TRUE if there are positions within the cutoff.
 */
gmx_bool
gmx_ana_nbsearch_pos_first_within(gmx_ana_nbsearch_t *d, gmx_ana_pos_t *p,
                                  int i, int *jp)
{
    return gmx_ana_nbsearch_first_within(d, p->x[i], jp);
}

/*!
 * \param[in]  d   Neighborhood search data structure.
 * \param[out] jp  Index of the test position in the next pair.
 * \returns    TRUE if there are positions within the cutoff.
 */
gmx_bool
gmx_ana_nbsearch_next_within(gmx_ana_nbsearch_t *d, int *jp)
{
    if (grid_search(d, &within_action))
    {
        *jp = d->previ;
        return TRUE;
    }
    *jp = -1;
    return FALSE;
}
