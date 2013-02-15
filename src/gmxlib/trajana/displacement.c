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
/*! \page displacements Displacement calculation routines
 *
 * Functions to calculate particle displacements as a function of time are
 * defined in displacement.h.
 * These functions allow one to easily calculate the displacement during a
 * single pass over the trajectory.
 * A new calculation is initialized with gmx_ana_displ_create().
 * For each frame, gmx_ana_displ_start_frame() needs to be called, followed by
 * a series of gmx_ana_displ_store() or a single call to
 * gmx_ana_displ_store_array(), gmx_ana_displ_store_all() or
 * gmx_ana_displ_store_pos().
 * The displacements for the current frame can then be calculated with
 * gmx_ana_displ_vector(), gmx_ana_displ_vectors() and/or
 * gmx_ana_displ_vectors_all().
 * These calculate the displacements from old positions to the ones on the
 * current frame.
 * The input frames should be evenly spaced, otherwise a fatal error is
 * reported. A time interval can be converted to a number of steps for
 * the calculation functions using gmx_ana_displ_time_to_steps().
 * When the structure is no longer required, gmx_ana_displ_free() can be used
 * to free the memory.
 *
 * The functions support calculation for dynamic selections:
 * when storing the positions, it is possible to specify whether a particle
 * is present or not.
 * This data is then accessible in the \p pout array returned by the
 * calculation functions.
 * Note that if you only call gmx_ana_displ_store() for the particles that are
 * present, use gmx_ana_displ_store_array(), or use positions calculated without
 * \ref POS_MASKONLY for gmx_ana_displ_store_pos(), you should ensure that you
 * do not use the displacements for which \p pout is FALSE (the values cannot
 * be calculated based on the provided data, and are undefined).
 */
/*! \internal \file
 * \brief Implementation of functions in displacement.h.
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>

#include <pbc.h>
#include <smalloc.h>
#include <vec.h>

#include <displacement.h>
#include <position.h>

/*! \internal \brief
 * Internal structure used by the displacement calculation routines in
 * displacement.h.
 */
typedef struct gmx_ana_displpos_t
{
    /** Stored position vector. */
    rvec                     x;
    /** TRUE if there is something stored. */
    gmx_bool                 bPres;
} gmx_ana_displpos_t;

/*! \internal \brief
 * Stores displacement calculation data and parameters.
 *
 * There data can be accessed only through the functions defined in
 * displacement.h.
 */
struct gmx_ana_displ_t
{
    /** Maximum number of particles for which the displacements are calculated. */
    int                  nmax;
    /** Maximum time for which the displacements are needed. */
    real                 tmax;

    /** TRUE if no frames have been read. */
    gmx_bool                 bFirst;
    /** Stores the time of the first frame. */
    real                     t0;
    /** Stores the time interval between frames. */
    real                     dt;
    /** Stores the time of the current frame. */
    real                     t;
    /** Stores the index in the store for the current positions. */
    int                      ci;

    /** Maximum number of positions to store for a particle. */
    int                  max_store;
    /** The total number of positions ever stored (can be larger than \p max_store). */
    int                  nstored;
    /** Two-dimensional array of stored positions, second index is the particle position. */
    gmx_ana_displpos_t **p;
    /** Pointer to the memory allocated for \p p. */
    gmx_ana_displpos_t  *palloc;
};

/*!
 * \param[out] data Displacement calculation data structure pointer to
 *   initialize.
 * \param[in]  nmax Number of particles
 *   for which displacements should be calculated.
 * \param[in]  tmax Maximum interval (in same units as in the trajectory)
 *   for which displacements are required.
 * \returns    0 on success.
 *
 * This function only allocates memory for the first frame;
 * gmx_ana_displ_start_frame() allocates a sufficient amount of memory for
 * storing all the required positions once it known the interval between
 * the trajectory frames.
 */
int
gmx_ana_displ_create(gmx_ana_displ_t **data, int nmax, real tmax)
{
    gmx_ana_displ_t *d;

    snew(d, 1);
    d->nmax        = nmax;
    d->tmax        = tmax;

    d->bFirst      = TRUE;
    d->t0          = 0;
    d->dt          = 0;
    d->t           = 0;
    d->ci          = -1;
    d->nstored     = 0;

    d->max_store   = -1;
    snew(d->palloc, nmax);
    snew(d->p,      1);
    d->p[0]        = d->palloc;

    *data = d;
    return 0;
}

/*!
 * \param[in,out] d   Displacement calculation data.
 * \param[in]     t   Time for the new frame.
 * \returns       0 on success, -1 on failure.
 *
 * This function should be called at the beginning of each frame for which
 * atomic positions should be stored.
 * It is required that the times for which this function is called are
 * evenly spaced; otherwise, it reports a fatal error.
 */
int
gmx_ana_displ_start_frame(gmx_ana_displ_t *d, real t)
{
    int   i;

    /* Initialize times */
    if (d->bFirst)
    {
        d->t0 = t;
    }
    else if (d->dt <= 0)
    {
        d->dt = t - d->t0;
    }
    else
    {
        if (!gmx_within_tol(t - d->t, d->dt, GMX_REAL_EPS))
        {
            gmx_input("Trajectory not evenly spaced");
            return -1;
        }
    }
    d->t = t;

    /* Allocate memory for all the positions once it is possible */
    if (d->max_store == -1 && !d->bFirst)
    {
        d->max_store = (int)(d->tmax/d->dt + 1);
        srenew(d->palloc, d->nmax * d->max_store);
        sfree(d->p);
        snew(d->p, d->max_store);
        for (i = 0; i < d->max_store; ++i)
        {
            d->p[i] = &d->palloc[d->nmax*i];
        }
    }

    /* Increment the index where current positions are stored */
    d->ci++;
    if (d->ci >= d->max_store)
    {
        d->ci = 0;
    }

    for (i = 0; i < d->nmax; ++i)
    {
        d->p[d->ci][i].bPres = FALSE;
    }
    d->nstored++;
    d->bFirst = FALSE;
    return 0;
}

/*!
 * \param[in]     d     Displacement calculation data.
 * \param[in]     time  Time to convert to steps.
 * \param[out]    steps The number of steps between the frames separated by
 *   \p time, 0 if the number of steps cannot be determined.
 * \returns       0 on success (also if \p *steps is 0), -1 on error.
 *
 * \p *steps is set to zero if gmx_ana_displ_start_frame() has not yet been
 * called twice.
 * -1 is returned if \p time is not an integer multiple of the interval
 * between frames.
 */
int
gmx_ana_displ_time_to_steps(gmx_ana_displ_t *d, real time, int *steps)
{
    if (d->dt <= 0)
    {
        *steps = 0;
        return 0;
    }
    if (!gmx_within_tol(fmod(time, d->dt), 0, GMX_REAL_EPS))
    {
        gmx_input("Interval not multiple of trajectory time step");
        return -1;
    }
    *steps = (int)(time / d->dt + 0.5);
    return 0;
}

/*!
 * \param[in,out] d     Displacement calculation data.
 * \param[in]     id    Particle ID.
 * \param[in]     x     Particle position.
 * \param[in]     bPres TRUE if the particle should be marked as present.
 * \returns       0 on success.
 */
int
gmx_ana_displ_store(gmx_ana_displ_t *d, atom_id id, rvec x, gmx_bool bPres)
{
    copy_rvec(x, d->p[d->ci][id].x);
    d->p[d->ci][id].bPres = bPres;
    return 0;
}

/*!
 * \param[in,out] d     Displacement calculation data.
 * \param[in]     n     Number of elements in the \p id and \p x arrays.
 * \param[in]     id    Particle IDs.
 * \param[in]     x     Particle positions.
 * \returns       0 on success.
 *
 * The positions of the given \p n particles are stored, and all of them
 * are marked as present.
 */
int
gmx_ana_displ_store_array(gmx_ana_displ_t *d, int n, atom_id id[], rvec x[])
{
    int i;

    for (i = 0; i < n; ++i)
    {
        gmx_ana_displ_store(d, id[i], x[i], TRUE);
    }
    return 0;
}

/*!
 * \param[in,out] d     Displacement calculation data.
 * \param[in]     id    Particle IDs.
 * \param[in]     x     Particle positions.
 * \returns       0 on success.
 *
 * Both the \p id and \p x arrays should have
 * \p nmax items, where \p nmax was provided for gmx_ana_displ_create(), and
 * the \p id[i] should be either \p i or -1.
 * If \p id[i]==-1, then that particle is marked as not present,
 * but the position is still stored.
 */
int
gmx_ana_displ_store_all(gmx_ana_displ_t *d, atom_id id[], rvec x[])
{
    int i;

    for (i = 0; i < d->nmax; ++i)
    {
        gmx_ana_displ_store(d, i, x[i], id[i] >= 0);
    }
    return 0;
}

/*!
 * \param[in,out] d     Displacement calculation data.
 * \param[in]     p     Particle positions.
 * \returns       0 on success.
 *
 * \p p should have at most \p nmax positions, where \p nmax was provided for
 * gmx_ana_displ_create().
 * If \p has exactly \p nmax positions, gmx_ana_displ_store_all() is called,
 * otherwise gmx_ana_displ_store_array() is used.
 */
int
gmx_ana_displ_store_pos(gmx_ana_displ_t *d, gmx_ana_pos_t *p)
{
    if (p->nr == d->nmax)
    {
        gmx_ana_displ_store_all(d, p->m.refid, p->x);
    }
    else
    {
        gmx_ana_displ_store_array(d, p->nr, p->m.refid, p->x);
    }
    return 0;
}

/*! \brief
 * Finds the index in the displacement position storage array.
 *
 * \param[in] d    Displacement calculation data.
 * \param[in] step Calculate the index of positions \p steps back.
 * \returns   Index into the \ref gmx_ana_displ_t::p "d->p" array
 *   that contains the positions \p steps back, -1 on error.
 *
 * If \p step is too large, -1 is returned.
 */
static int
find_store_index(gmx_ana_displ_t *d, int step)
{
    int si;

    si = d->ci - step;
    if (si < 0)
    {
        si += d->max_store;
    }
    if (si < 0)
    {
        gmx_incons("Displacement requested for an interval too long");
        return -1;
    }

    return si;
}

/*!
 * \param[in]     d     Displacement calculation data.
 * \param[in]     step  Displacement is calculated from the location
 *   \p step steps back.
 * \param[in]     pbc   Periodic boundary structure
 *   (if NULL, PBC are not used in distance calculation).
 * \param[in]     id    Particle ID for which the displacement is calculated.
 * \param[in]     x     Current particle position.
 * \param[out]    xout  Displacement of the particle.
 * \param[out]    pout  TRUE if the old particle position was marked present.
 *   If \p *pout is FALSE and the old particle position was not stored,
 *   the value of \p xout is undefined.
 *   Can be NULL, in which case the value is not stored.
 * \return        0 on success, -1 if less than \p step frames have been
 *   stored or \p step <= 0, or EINVAL if \p step is too large.
 */
int
gmx_ana_displ_vector(gmx_ana_displ_t *d, int step, t_pbc *pbc,
                     atom_id id, rvec x, rvec xout, gmx_bool *pout)
{
    int si;

    if (step >= d->nstored || step <= 0)
    {
        return -1;
    }
    si = find_store_index(d, step);
    if (si == -1)
    {
        return EINVAL;
    }
    if (pout)
    {
        *pout = d->p[si][id].bPres;
    }
    if (pbc)
    {
        pbc_dx(pbc, x, d->p[si][id].x, xout);
    }
    else
    {
        rvec_sub(x, d->p[si][id].x, xout);
    }
    return 0;
}

/*!
 * \param[in]     d     Displacement calculation data.
 * \param[in]     step  Displacement is calculated from the location
 *   \p step steps back.
 * \param[in]     pbc   Periodic boundary structure
 *   (if NULL, PBC are not used in distance calculation).
 * \param[in]     n     Number of particles for which to calculate.
 * \param[in]     id    Particle IDs for which the displacement is calculated.
 * \param[in]     x     Current particle positions.
 * \param[out]    xout  Displacement of the particles.
 * \param[out]    pout  TRUE if the old particle position was marked present.
 *   If \p pout[i] is FALSE and the old particle position was not stored,
 *   the value of \p xout[i] is undefined.
 *   Can be NULL, in which case the value is not stored.
 * \return        0 on success, -1 if less than \p step frames have been
 *   stored or \p step <= 0, or EINVAL if \p step is too large.
 */
int
gmx_ana_displ_vectors(gmx_ana_displ_t *d, int step, t_pbc *pbc,
                      int n, atom_id id[], rvec x[], rvec xout[], gmx_bool *pout)
{
    int si, i;

    if (step >= d->nstored || step <= 0)
    {
        return -1;
    }
    si = find_store_index(d, step);
    if (si == -1)
    {
        return EINVAL;
    }
    for (i = 0; i < n; ++i)
    {
        if (pout)
        {
            pout[i] =  d->p[si][id[i]].bPres;
        }
        if (pbc)
        {
            pbc_dx(pbc, x[i], d->p[si][id[i]].x, xout[i]);
        }
        else
        {
            rvec_sub(x[i], d->p[si][id[i]].x, xout[i]);
        }
    }
    return 0;
}

/*!
 * \param[in]     d     Displacement calculation data.
 * \param[in]     step  Displacement is calculated from the location
 *   \p step steps back.
 * \param[in]     pbc   Periodic boundary structure
 *   (if NULL, PBC are not used in distance calculation).
 * \param[in]     x     Current particle positions.
 * \param[out]    xout  Displacement of the particles.
 * \param[out]    pout  TRUE if the old particle position was marked present.
 *   If \p pout[i] is FALSE and the old particle position was not stored,
 *   the value of \p xout[i] is undefined.
 *   Can be NULL, in which case the value is not stored.
 * \return        0 on success, -1 if less than \p step frames have been
 *   stored or \p step <= 0, or EINVAL if \p step is too large.
 *
 * The \p x array should have \p nmax items, where \p nmax is the one provided
 * to init_displ_calc().
 */
int
gmx_ana_displ_vectors_all(gmx_ana_displ_t *d, int step, t_pbc *pbc,
                          rvec x[], rvec xout[], gmx_bool *pout)
{
    int si, i;

    if (step >= d->nstored || step <= 0)
    {
        return -1;
    }
    si = find_store_index(d, step);
    if (si == -1)
    {
        return EINVAL;
    }
    for (i = 0; i < d->nmax; ++i)
    {
        if (pout)
        {
            pout[i] =  d->p[si][i].bPres;
        }
        if (pbc)
        {
            pbc_dx(pbc, x[i], d->p[si][i].x, xout[i]);
        }
        else
        {
            rvec_sub(x[i], d->p[si][i].x, xout[i]);
        }
    }
    return 0;
}

/*!
 * \param[in,out] d  Displacement calculation data.
 *
 * After the call, the pointer \p d is invalid.
 */
void
gmx_ana_displ_free(gmx_ana_displ_t *d)
{
    sfree(d->p);
    sfree(d->palloc);
    sfree(d);
}
