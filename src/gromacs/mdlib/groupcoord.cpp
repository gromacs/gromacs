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
#include "gmxpre.h"

#include "groupcoord.h"

#include "gromacs/domdec/ga2la.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/vec.h"

#define MIN(a, b) (((a) < (b)) ? (a) : (b))


/* Select the indices of the group's atoms which are local and store them in
 * anrs_loc[0..nr_loc]. The indices are saved in coll_ind[] for later reduction
 * in communicate_group_positions()
 */
void dd_make_local_group_indices(const gmx_ga2la_t* ga2la,
                                 const int nr,    /* IN:  Total number of atoms in the group */
                                 int  anrs[],     /* IN:  Global atom numbers of the groups atoms */
                                 int* nr_loc,     /* OUT: Number of group atoms found locally */
                                 int* anrs_loc[], /* OUT: Local atom numbers of the group  */
                                 int* nalloc_loc, /* IN+OUT: Allocation size of anrs_loc */
                                 int coll_ind[]) /* OUT (opt): Where is this position found in the collective array? */
{
    GMX_ASSERT(ga2la, "We need a valid ga2la object");
    GMX_RELEASE_ASSERT(anrs != *anrs_loc, "Can not update indices in-place");

    /* Loop over all the atom indices of the group to check
     * which ones are on the local node */
    int localnr = 0;
    for (int i = 0; i < nr; i++)
    {
        if (const int* ii = ga2la->findHome(anrs[i]))
        {
            /* The atom with this index is a home atom */
            if (localnr >= *nalloc_loc) /* Check whether memory suffices */
            {
                *nalloc_loc = over_alloc_dd(localnr + 1);
                /* We never need more memory than the number of atoms in the group */
                *nalloc_loc = MIN(*nalloc_loc, nr);
                srenew(*anrs_loc, *nalloc_loc);
            }
            /* Save the atoms index in the local atom numbers array */
            (*anrs_loc)[localnr] = *ii;

            if (coll_ind != nullptr)
            {
                /* Keep track of where this local atom belongs in the collective index array.
                 * This is needed when reducing the local arrays to a collective/global array
                 * in communicate_group_positions */
                coll_ind[localnr] = i;
            }

            /* add one to the local atom count */
            localnr++;
        }
    }

    /* Return the number of local atoms that were found */
    *nr_loc = localnr;
}


static void get_shifts_group(int          npbcdim,
                             const matrix box,
                             rvec*        xcoll, /* IN:  Collective set of positions [0..nr] */
                             int          nr,    /* IN:  Total number of atoms in the group */
                             rvec* xcoll_old, /* IN:  Positions from the last time step [0...nr] */
                             ivec* shifts)    /* OUT: Shifts for xcoll */
{
    int  i, m, d;
    rvec dx;


    /* Get the shifts such that each atom is within closest
     * distance to its position at the last NS time step after shifting.
     * If we start with a whole group, and always keep track of
     * shift changes, the group will stay whole this way */
    for (i = 0; i < nr; i++)
    {
        clear_ivec(shifts[i]);
    }

    for (i = 0; i < nr; i++)
    {
        /* The distance this atom moved since the last time step */
        /* If this is more than just a bit, it has changed its home pbc box */
        rvec_sub(xcoll[i], xcoll_old[i], dx);

        for (m = npbcdim - 1; m >= 0; m--)
        {
            while (dx[m] < -0.5 * box[m][m])
            {
                for (d = 0; d < DIM; d++)
                {
                    dx[d] += box[m][d];
                }
                shifts[i][m]++;
            }
            while (dx[m] >= 0.5 * box[m][m])
            {
                for (d = 0; d < DIM; d++)
                {
                    dx[d] -= box[m][d];
                }
                shifts[i][m]--;
            }
        }
    }
}


static void shift_positions_group(const matrix box,
                                  rvec         x[], /* The positions [0..nr] */
                                  ivec*        is,  /* The shifts [0..nr] */
                                  int          nr)           /* The number of positions and shifts */
{
    int i, tx, ty, tz;


    /* Loop over the group's atoms */
    if (TRICLINIC(box))
    {
        for (i = 0; i < nr; i++)
        {
            tx = is[i][XX];
            ty = is[i][YY];
            tz = is[i][ZZ];

            x[i][XX] = x[i][XX] + tx * box[XX][XX] + ty * box[YY][XX] + tz * box[ZZ][XX];
            x[i][YY] = x[i][YY] + ty * box[YY][YY] + tz * box[ZZ][YY];
            x[i][ZZ] = x[i][ZZ] + tz * box[ZZ][ZZ];
        }
    }
    else
    {
        for (i = 0; i < nr; i++)
        {
            tx = is[i][XX];
            ty = is[i][YY];
            tz = is[i][ZZ];

            x[i][XX] = x[i][XX] + tx * box[XX][XX];
            x[i][YY] = x[i][YY] + ty * box[YY][YY];
            x[i][ZZ] = x[i][ZZ] + tz * box[ZZ][ZZ];
        }
    }
}


/* Assemble the positions of the group such that every node has all of them.
 * The atom indices are retrieved from anrs_loc[0..nr_loc]
 * Note that coll_ind[i] = i is needed in the serial case */
extern void
communicate_group_positions(const gmx::MpiComm& mpiComm, /* Reference to MPI communication data */
                            rvec*               xcoll,   /* Collective array of positions */
                            ivec* shifts, /* Collective array of shifts for xcoll (can be NULL) */
                            ivec* extra_shifts,   /* (optional) Extra shifts since last time step */
                            const gmx_bool bNS,   /* (optional) NS step, the shifts have changed */
                            const rvec*    x_loc, /* Local positions on this node */
                            const int      nr,    /* Total number of atoms in the group */
                            const int      nr_loc,    /* Local number of atoms in the group */
                            const int*     anrs_loc,  /* Local atom numbers */
                            const int*     coll_ind,  /* Collective index */
                            rvec*          xcoll_old, /* (optional) Positions from the last time
                                                         step,  used to make group whole */
                            const matrix box)         /* (optional) The box */
{
    int i;


    /* Zero out the groups' global position array */
    clear_rvecs(nr, xcoll);

    /* Put the local positions that this node has into the right place of
     * the collective array. Note that in the serial case, coll_ind[i] = i */
    for (i = 0; i < nr_loc; i++)
    {
        copy_rvec(x_loc[anrs_loc[i]], xcoll[coll_ind[i]]);
    }

    if (mpiComm.size() > 1)
    {
        /* Add the arrays from all nodes together */
        mpiComm.sumReduce(nr * 3, xcoll[0]);
    }
    /* Now we have all the positions of the group in the xcoll array present on all
     * nodes.
     *
     * The rest of the code is for making the group whole again in case atoms changed
     * their PBC representation / crossed a box boundary. We only do that if the
     * shifts array is allocated. */
    if (nullptr != shifts)
    {
        /* To make the group whole, start with a whole group and each
         * step move the assembled positions at closest distance to the positions
         * from the last step. First shift the positions with the saved shift
         * vectors (these are 0 when this routine is called for the first time!) */
        shift_positions_group(box, xcoll, shifts, nr);

        /* Now check if some shifts changed since the last step.
         * This only needs to be done when the shifts are expected to have changed,
         * i.e. after neighbor searching */
        if (bNS)
        {
            get_shifts_group(3, box, xcoll, nr, xcoll_old, extra_shifts);

            /* Shift with the additional shifts such that we get a whole group now */
            shift_positions_group(box, xcoll, extra_shifts, nr);

            /* Add the shift vectors together for the next time step */
            for (i = 0; i < nr; i++)
            {
                shifts[i][XX] += extra_shifts[i][XX];
                shifts[i][YY] += extra_shifts[i][YY];
                shifts[i][ZZ] += extra_shifts[i][ZZ];
            }

            /* Store current correctly-shifted positions for comparison in the next NS time step */
            for (i = 0; i < nr; i++)
            {
                copy_rvec(xcoll[i], xcoll_old[i]);
            }
        }
    }
}


/* Determine the (weighted) sum vector from positions x */
extern double get_sum_of_positions(const rvec x[], real weight[], const int nat, dvec dsumvec)
{
    int    i;
    rvec   x_weighted;
    double weight_sum = 0.0;


    /* Zero out the center */
    clear_dvec(dsumvec);

    /* Loop over all atoms and add their weighted position vectors */
    if (weight != nullptr)
    {
        for (i = 0; i < nat; i++)
        {
            weight_sum += weight[i];
            svmul(weight[i], x[i], x_weighted);
            dsumvec[XX] += x_weighted[XX];
            dsumvec[YY] += x_weighted[YY];
            dsumvec[ZZ] += x_weighted[ZZ];
        }
    }
    else
    {
        for (i = 0; i < nat; i++)
        {
            dsumvec[XX] += x[i][XX];
            dsumvec[YY] += x[i][YY];
            dsumvec[ZZ] += x[i][ZZ];
        }
    }
    return weight_sum;
}


/* Determine center of structure from collective positions x */
extern void get_center(const rvec x[], real weight[], const int nr, rvec rcenter)
{
    dvec   dcenter;
    double weight_sum, denom;


    weight_sum = get_sum_of_positions(x, weight, nr, dcenter);

    if (weight != nullptr)
    {
        denom = weight_sum; /* Divide by the sum of weight */
    }
    else
    {
        denom = nr; /* Divide by the number of atoms */
    }
    dsvmul(1.0 / denom, dcenter, dcenter);

    rcenter[XX] = dcenter[XX];
    rcenter[YY] = dcenter[YY];
    rcenter[ZZ] = dcenter[ZZ];
}


/* Get the center from local positions that already have the correct
 * PBC representation */
extern void get_center_comm(const gmx::MpiComm& mpiComm,
                            rvec                x_loc[],      /* Local positions */
                            real                weight_loc[], /* Local masses or other weights */
                            int                 nr_loc,       /* Local number of atoms */
                            int                 nr_group, /* Total number of atoms of the group */
                            rvec                center)                  /* Weighted center */
{
    double weight_sum, denom;
    dvec   dsumvec;

    weight_sum = get_sum_of_positions(x_loc, weight_loc, nr_loc, dsumvec);

    /* Add the local contributions from all nodes. Put the sum vector and the
     * weight in a buffer array so that we get along with a single communication
     * call. */
    if (mpiComm.size() > 1)
    {
        std::array<double, 4> buf;

        buf[0] = dsumvec[XX];
        buf[1] = dsumvec[YY];
        buf[2] = dsumvec[ZZ];
        buf[3] = weight_sum;

        /* Communicate buffer */
        mpiComm.sumReduce(buf);

        dsumvec[XX] = buf[0];
        dsumvec[YY] = buf[1];
        dsumvec[ZZ] = buf[2];
        weight_sum  = buf[3];
    }

    if (weight_loc != nullptr)
    {
        denom = 1.0 / weight_sum; /* Divide by the sum of weight to get center of mass e.g. */
    }
    else
    {
        denom = 1.0 / nr_group; /* Divide by the number of atoms to get the geometrical center */
    }
    center[XX] = dsumvec[XX] * denom;
    center[YY] = dsumvec[YY] * denom;
    center[ZZ] = dsumvec[ZZ] * denom;
}


/* Translate x with transvec */
extern void translate_x(rvec x[], const int nr, const rvec transvec)
{
    int i;


    for (i = 0; i < nr; i++)
    {
        rvec_inc(x[i], transvec);
    }
}


extern void rotate_x(rvec x[], const int nr, matrix rmat)
{
    int  i, j, k;
    rvec x_old;


    /* Apply the rotation matrix */
    for (i = 0; i < nr; i++)
    {
        for (j = 0; j < 3; j++)
        {
            x_old[j] = x[i][j];
        }
        for (j = 0; j < 3; j++)
        {
            x[i][j] = 0;
            for (k = 0; k < 3; k++)
            {
                x[i][j] += rmat[k][j] * x_old[k];
            }
        }
    }
}
