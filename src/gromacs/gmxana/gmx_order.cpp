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

#include <cmath>
#include <cstring>

#include <algorithm>
#include <optional>

#include "gromacs/commandline/pargs.h"
#include "gromacs/commandline/viewit.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/gmxana/cmat.h"
#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/gmxana/gstat.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/units.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pbcutil/rmpbc.h"
#include "gromacs/topology/index.h"
#include "gromacs/topology/topology.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"

/****************************************************************************/
/* This program calculates the order parameter per atom for an interface or */
/* bilayer, averaged over time.                                             */
/* S = 1/2 * (3 * cos(i)cos(j) - delta(ij))                                 */
/* It is assumed that the order parameter with respect to a box-axis        */
/* is calculated. In that case i = j = axis, and delta(ij) = 1.             */
/*                                                                          */
/* Peter Tieleman,  April 1995                                              */
/* P.J. van Maaren, November 2005     Added tetrahedral stuff               */
/****************************************************************************/

static void find_nearest_neighbours(PbcType     pbcType,
                                    int         natoms,
                                    matrix      box,
                                    rvec        x[],
                                    int         maxidx,
                                    const int   index[],
                                    real*       sgmean,
                                    real*       skmean,
                                    int         nslice,
                                    int         slice_dim,
                                    real        sgslice[],
                                    real        skslice[],
                                    gmx_rmpbc_t gpbc)
{
    int   ix, jx, nsgbin, *sgbin;
    int   i, ibin, j, k, *nn[4];
    rvec  dx, rj, rk, urk, urj;
    real  cost, cost2, *sgmol, *skmol, rmean, rmean2, r2, box2, *r_nn[4];
    t_pbc pbc;
    int   sl_index;
    real* sl_count;
    real  onethird = 1.0 / 3.0;
    /*  dmat = init_mat(maxidx, FALSE); */
    box2 = box[XX][XX] * box[XX][XX];
    snew(sl_count, nslice);
    for (i = 0; (i < 4); i++)
    {
        snew(r_nn[i], natoms);
        snew(nn[i], natoms);

        for (j = 0; (j < natoms); j++)
        {
            r_nn[i][j] = box2;
        }
    }

    snew(sgmol, maxidx);
    snew(skmol, maxidx);

    /* Must init pbc every step because of pressure coupling */
    set_pbc(&pbc, pbcType, box);

    gmx_rmpbc_apply(gpbc, natoms, box, x);

    nsgbin = 2001; // Calculated as (1 + 1/0.0005)
    snew(sgbin, nsgbin);

    *sgmean = 0.0;
    *skmean = 0.0;
    for (i = 0; (i < maxidx); i++) /* loop over index file */
    {
        ix = index[i];
        for (j = 0; (j < maxidx); j++)
        {
            if (i == j)
            {
                continue;
            }

            jx = index[j];

            pbc_dx(&pbc, x[ix], x[jx], dx);
            r2 = iprod(dx, dx);

            /* set_mat_entry(dmat,i,j,r2); */

            /* determine the nearest neighbours */
            if (r2 < r_nn[0][i])
            {
                r_nn[3][i] = r_nn[2][i];
                nn[3][i]   = nn[2][i];
                r_nn[2][i] = r_nn[1][i];
                nn[2][i]   = nn[1][i];
                r_nn[1][i] = r_nn[0][i];
                nn[1][i]   = nn[0][i];
                r_nn[0][i] = r2;
                nn[0][i]   = j;
            }
            else if (r2 < r_nn[1][i])
            {
                r_nn[3][i] = r_nn[2][i];
                nn[3][i]   = nn[2][i];
                r_nn[2][i] = r_nn[1][i];
                nn[2][i]   = nn[1][i];
                r_nn[1][i] = r2;
                nn[1][i]   = j;
            }
            else if (r2 < r_nn[2][i])
            {
                r_nn[3][i] = r_nn[2][i];
                nn[3][i]   = nn[2][i];
                r_nn[2][i] = r2;
                nn[2][i]   = j;
            }
            else if (r2 < r_nn[3][i])
            {
                r_nn[3][i] = r2;
                nn[3][i]   = j;
            }
        }


        /* calculate mean distance between nearest neighbours */
        rmean = 0;
        for (j = 0; (j < 4); j++)
        {
            r_nn[j][i] = std::sqrt(r_nn[j][i]);
            rmean += r_nn[j][i];
        }
        rmean /= 4;

        sgmol[i] = 0.0;
        skmol[i] = 0.0;

        /* Chau1998a eqn 3 */
        /* angular part tetrahedrality order parameter per atom */
        for (j = 0; (j < 3); j++)
        {
            for (k = j + 1; (k < 4); k++)
            {
                pbc_dx(&pbc, x[ix], x[index[nn[k][i]]], rk);
                pbc_dx(&pbc, x[ix], x[index[nn[j][i]]], rj);

                unitv(rk, urk);
                unitv(rj, urj);

                cost  = iprod(urk, urj) + onethird;
                cost2 = cost * cost;

                /* sgmol[i] += 3*cost2/32;  */
                sgmol[i] += cost2;

                /* determine distribution */
                ibin = static_cast<int>(static_cast<real>(nsgbin) * cost2);
                if (ibin < nsgbin)
                {
                    sgbin[ibin]++;
                }
                /* printf("%d %d %f %d %d\n", j, k, cost * cost, ibin, sgbin[ibin]);*/
            }
        }

        /* normalize sgmol between 0.0 and 1.0 */
        sgmol[i] = 3 * sgmol[i] / 32;
        *sgmean += sgmol[i];

        /* distance part tetrahedrality order parameter per atom */
        rmean2 = 4 * 3 * rmean * rmean;
        for (j = 0; (j < 4); j++)
        {
            skmol[i] += (rmean - r_nn[j][i]) * (rmean - r_nn[j][i]) / rmean2;
            /*      printf("%d %f (%f %f %f %f) \n",
                i, skmol[i], rmean, rmean2, r_nn[j][i], (rmean - r_nn[j][i]) );
             */
        }

        *skmean += skmol[i];

        /* Compute sliced stuff */
        sl_index = static_cast<int>(std::round((1 + x[i][slice_dim] / box[slice_dim][slice_dim])
                                               * static_cast<real>(nslice)))
                   % nslice;
        sgslice[sl_index] += sgmol[i];
        skslice[sl_index] += skmol[i];
        sl_count[sl_index]++;
    } /* loop over entries in index file */

    *sgmean /= static_cast<real>(maxidx);
    *skmean /= static_cast<real>(maxidx);

    for (i = 0; (i < nslice); i++)
    {
        if (sl_count[i] > 0)
        {
            sgslice[i] /= sl_count[i];
            skslice[i] /= sl_count[i];
        }
    }
    sfree(sl_count);
    sfree(sgbin);
    sfree(sgmol);
    sfree(skmol);
    for (i = 0; (i < 4); i++)
    {
        sfree(r_nn[i]);
        sfree(nn[i]);
    }
}


static void calc_tetra_order_parm(const char*             fnNDX,
                                  const char*             fnTPS,
                                  const char*             fnTRX,
                                  const char*             sgfn,
                                  const char*             skfn,
                                  int                     nslice,
                                  int                     slice_dim,
                                  const char*             sgslfn,
                                  const char*             skslfn,
                                  const gmx_output_env_t* oenv)
{
    FILE *       fpsg = nullptr, *fpsk = nullptr;
    t_topology   top;
    PbcType      pbcType;
    t_trxstatus* status;
    int          natoms;
    real         t;
    rvec *       xtop, *x;
    matrix       box;
    real         sg, sk;
    int**        index;
    char**       grpname;
    int          i, *isize, ng, nframes;
    real *       sg_slice, *sg_slice_tot, *sk_slice, *sk_slice_tot;
    gmx_rmpbc_t  gpbc = nullptr;


    read_tps_conf(fnTPS, &top, &pbcType, &xtop, nullptr, box, FALSE);

    snew(sg_slice, nslice);
    snew(sk_slice, nslice);
    snew(sg_slice_tot, nslice);
    snew(sk_slice_tot, nslice);
    ng = 1;
    /* get index groups */
    printf("Select the group that contains the atoms you want to use for the tetrahedrality order "
           "parameter calculation:\n");
    snew(grpname, ng);
    snew(index, ng);
    snew(isize, ng);
    get_index(&top.atoms, fnNDX, ng, isize, index, grpname);

    /* Analyze trajectory */
    natoms = read_first_x(oenv, &status, fnTRX, &t, &x, box);
    if (natoms > top.atoms.nr)
    {
        gmx_fatal(FARGS, "Topology (%d atoms) does not match trajectory (%d atoms)", top.atoms.nr, natoms);
    }
    check_index(nullptr, ng, index[0], nullptr, natoms);

    fpsg = xvgropen(sgfn, "S\\sg\\N Angle Order Parameter", "Time (ps)", "S\\sg\\N", oenv);
    fpsk = xvgropen(skfn, "S\\sk\\N Distance Order Parameter", "Time (ps)", "S\\sk\\N", oenv);

    /* loop over frames */
    gpbc    = gmx_rmpbc_init(&top.idef, pbcType, natoms);
    nframes = 0;
    do
    {
        find_nearest_neighbours(
                pbcType, natoms, box, x, isize[0], index[0], &sg, &sk, nslice, slice_dim, sg_slice, sk_slice, gpbc);
        for (i = 0; (i < nslice); i++)
        {
            sg_slice_tot[i] += sg_slice[i];
            sk_slice_tot[i] += sk_slice[i];
        }
        fprintf(fpsg, "%f %f\n", t, sg);
        fprintf(fpsk, "%f %f\n", t, sk);
        nframes++;
    } while (read_next_x(oenv, status, &t, x, box));
    close_trx(status);
    gmx_rmpbc_done(gpbc);

    sfree(grpname);
    sfree(index);
    sfree(isize);

    xvgrclose(fpsg);
    xvgrclose(fpsk);

    fpsg = xvgropen(sgslfn, "S\\sg\\N Angle Order Parameter / Slab", "(nm)", "S\\sg\\N", oenv);
    fpsk = xvgropen(skslfn, "S\\sk\\N Distance Order Parameter / Slab", "(nm)", "S\\sk\\N", oenv);
    for (i = 0; (i < nslice); i++)
    {
        fprintf(fpsg,
                "%10g  %10g\n",
                (i + 0.5) * box[slice_dim][slice_dim] / nslice,
                sg_slice_tot[i] / static_cast<real>(nframes));
        fprintf(fpsk,
                "%10g  %10g\n",
                (i + 0.5) * box[slice_dim][slice_dim] / nslice,
                sk_slice_tot[i] / static_cast<real>(nframes));
    }
    xvgrclose(fpsg);
    xvgrclose(fpsk);
}


/* Print name of first atom in all groups in index file */
static void print_types(gmx::ArrayRef<const IndexGroup> indexGroups, const t_topology* top)
{
    fprintf(stderr, "Using following groups: \n");
    for (const auto& indexGroup : indexGroups)
    {
        fprintf(stderr,
                "Groupname: %s First atomname: %s First atomnr %d\n",
                indexGroup.name.c_str(),
                *(top->atoms.atomname[indexGroup.particleIndices[0]]),
                1 + indexGroup.particleIndices[0]);
    }
    fprintf(stderr, "\n");
}

static void check_length(real length, int a, int b)
{
    if (length > 0.3)
    {
        fprintf(stderr,
                "WARNING: distance between atoms %d and "
                "%d > 0.3 nm (%f). Index file might be corrupt.\n",
                a,
                b,
                length);
    }
}

static void calc_order(const char*                     fn,
                       gmx::ArrayRef<const IndexGroup> indexGroups,
                       rvec**                          order,
                       real***                         slOrder,
                       real*                           slWidth,
                       int                             nslices,
                       gmx_bool                        bSliced,
                       const t_topology*               top,
                       PbcType                         pbcType,
                       int                             axis,
                       gmx_bool                        permolecule,
                       gmx_bool                        radial,
                       gmx_bool                        distcalc,
                       const char*                     radfn,
                       real***                         distvals,
                       const gmx_output_env_t*         oenv)
{
    /* if permolecule = TRUE, order parameters will be calculed per molecule
     * and stored in slOrder with #slices = # molecules */
    rvec *x0,         /* coordinates with pbc                           */
            *x1;      /* coordinates without pbc                        */
    matrix       box; /* box (3x3)                                      */
    t_trxstatus* status;
    rvec         cossum,       /* sum of vector angles for three axes            */
            Sx, Sy, Sz,        /* the three molecular axes                       */
            tmp1, tmp2,        /* temp. rvecs for calculating dot products       */
            frameorder;        /* order parameters for one frame                 */
    real* slFrameorder;        /* order parameter for one frame, per slice      */
    real  length,              /* total distance between two atoms               */
            t,                 /* time from trajectory                           */
            z_ave, z1, z2;     /* average z, used to det. which slice atom is in */
    int natoms,                /* nr. atoms in trj                               */
            nr_tails,          /* nr tails, to check if index file is correct    */
            size = 0,          /* nr. of atoms in group. same as nr_tails        */
            i, j, m, k, slice; /* current slice number                           */
    real nr_frames = 0;
    int* slCount;                    /* nr. of atoms in one slice                      */
    gmx_bool use_unitvector = FALSE; /* use a specified unit vector instead of axis to specify unit normal*/
    rvec        direction, com;
    int         comsize, distsize;
    int *       comidx = nullptr, *distidx = nullptr;
    char*       grpname = nullptr;
    t_pbc       pbc;
    real        arcdist, tmpdist;
    gmx_rmpbc_t gpbc = nullptr;

    /* PBC added for center-of-mass vector*/
    /* Initiate the pbc structure */
    std::memset(&pbc, 0, sizeof(pbc));

    if ((natoms = read_first_x(oenv, &status, fn, &t, &x0, box)) == 0)
    {
        gmx_fatal(FARGS, "Could not read coordinates from statusfile\n");
    }

    nr_tails = gmx::ssize(indexGroups[0].particleIndices);
    fprintf(stderr, "Number of elements in first group: %d\n", nr_tails);
    /* take first group as standard. Not rocksolid, but might catch error in index*/

    if (permolecule)
    {
        nslices = nr_tails;
        bSliced = FALSE; /*force slices off */
        fprintf(stderr, "Calculating order parameters for each of %d molecules\n", nslices);
    }

    if (radial)
    {
        use_unitvector = TRUE;
        fprintf(stderr, "Select an index group to calculate the radial membrane normal\n");
        get_index(&top->atoms, radfn, 1, &comsize, &comidx, &grpname);
    }
    if (distcalc)
    {
        if (grpname != nullptr)
        {
            sfree(grpname);
        }
        fprintf(stderr, "Select an index group to use as distance reference\n");
        get_index(&top->atoms, radfn, 1, &distsize, &distidx, &grpname);
        bSliced = FALSE; /*force slices off*/
    }

    if (use_unitvector && bSliced)
    {
        fprintf(stderr,
                "Warning:  slicing and specified unit vectors are not currently compatible\n");
    }

    const int ngrps = gmx::ssize(indexGroups);

    snew(slCount, nslices);
    snew(*slOrder, nslices);
    for (i = 0; i < nslices; i++)
    {
        snew((*slOrder)[i], ngrps);
    }
    if (distcalc)
    {
        snew(*distvals, nslices);
        for (i = 0; i < nslices; i++)
        {
            snew((*distvals)[i], ngrps);
        }
    }
    snew(*order, ngrps);
    snew(slFrameorder, nslices);
    snew(x1, natoms);

    if (bSliced)
    {
        *slWidth = box[axis][axis] / static_cast<real>(nslices);
        fprintf(stderr, "Box divided in %d slices. Initial width of slice: %f\n", nslices, *slWidth);
    }


#if 0
    nr_tails = index[1] - index[0];
    fprintf(stderr, "Number of elements in first group: %d\n", nr_tails);
    /* take first group as standard. Not rocksolid, but might catch error
       in index*/
#endif

    gpbc = gmx_rmpbc_init(&top->idef, pbcType, natoms);
    /*********** Start processing trajectory ***********/
    do
    {
        if (bSliced)
        {
            *slWidth = box[axis][axis] / static_cast<real>(nslices);
        }

        set_pbc(&pbc, pbcType, box);
        gmx_rmpbc_copy(gpbc, natoms, box, x0, x1);

        /* Now loop over all groups. There are ngrps groups, the order parameter can
           be calculated for grp 1 to grp ngrps - 1. For each group, loop over all
           atoms in group, which is index[i] to (index[i+1] - 1) See block.h. Of
           course, in this case index[i+1] -index[i] has to be the same for all
           groups, namely the number of tails. i just runs over all atoms in a tail,
           so for DPPC ngrps = 16 and i runs from 1 to 14, including 14
         */


        if (radial)
        {
            /*center-of-mass determination*/
            com[XX] = 0.0;
            com[YY] = 0.0;
            com[ZZ] = 0.0;
            for (j = 0; j < comsize; j++)
            {
                rvec_inc(com, x1[comidx[j]]);
            }
            svmul(1.0 / comsize, com, com);
        }
        rvec displacementFromReference;
        if (distcalc)
        {
            rvec dref = { 0.0, 0.0, 0.0 };
            for (j = 0; j < distsize; j++)
            {
                rvec_inc(dref, x1[distidx[j]]);
            }
            svmul(1.0 / distsize, dref, dref);
            if (radial)
            {
                pbc_dx(&pbc, dref, com, displacementFromReference);
                unitv(displacementFromReference, displacementFromReference);
            }
        }

        for (i = 1; i < ngrps - 1; i++)
        {
            clear_rvec(frameorder);

            size = gmx::ssize(indexGroups[i].particleIndices);
            if (size != nr_tails)
            {
                gmx_fatal(FARGS,
                          "grp %d does not have same number of"
                          " elements as grp 1\n",
                          i);
            }

            for (j = 0; j < size; j++)
            {
                if (radial)
                /*create unit vector*/
                {
                    pbc_dx(&pbc, x1[indexGroups[i].particleIndices[j]], com, direction);
                    unitv(direction, direction);
                    /*DEBUG*/
                    /*if (j==0)
                        fprintf(stderr,"X %f %f %f\tcom %f %f %f\tdirection %f %f %f\n",x1[a[index[i]+j]][0],x1[a[index[i]+j]][1],x1[a[index[i]+j]][2],com[0],com[1],com[2],
                            direction[0],direction[1],direction[2]);*/
                }

                rvec dist;
                /* get vector dist(Cn-1,Cn+1) for tail atoms */
                rvec_sub(x1[indexGroups[i + 1].particleIndices[j]],
                         x1[indexGroups[i - 1].particleIndices[j]],
                         dist);
                length = norm(dist); /* determine distance between two atoms */
                check_length(length,
                             indexGroups[i - 1].particleIndices[j],
                             indexGroups[i + 1].particleIndices[j]);
                svmul(1.0 / length, dist, Sz);
                /* Sz is now the molecular axis Sz, normalized and all that */

                /* now get Sx. Sx is normal to the plane of Cn-1, Cn and Cn+1 so
                   we can use the outer product of Cn-1->Cn and Cn+1->Cn, I hope */
                rvec_sub(x1[indexGroups[i + 1].particleIndices[j]],
                         x1[indexGroups[i].particleIndices[j]],
                         tmp1);
                rvec_sub(x1[indexGroups[i - 1].particleIndices[j]],
                         x1[indexGroups[i].particleIndices[j]],
                         tmp2);
                cprod(tmp1, tmp2, Sx);
                svmul(1.0 / norm(Sx), Sx, Sx);

                /* now we can get Sy from the outer product of Sx and Sz   */
                cprod(Sz, Sx, Sy);
                svmul(1.0 / norm(Sy), Sy, Sy);

                /* the square of cosine of the angle between dist and the axis.
                   Using the innerproduct, but two of the three elements are zero
                   Determine the sum of the orderparameter of all atoms in group
                 */
                if (use_unitvector)
                {
                    cossum[XX] = gmx::square(iprod(Sx, direction)); /* this is allowed, since Sa is normalized */
                    cossum[YY] = gmx::square(iprod(Sy, direction));
                    cossum[ZZ] = gmx::square(iprod(Sz, direction));
                }
                else
                {
                    cossum[XX] = gmx::square(Sx[axis]); /* this is allowed, since Sa is normalized */
                    cossum[YY] = gmx::square(Sy[axis]);
                    cossum[ZZ] = gmx::square(Sz[axis]);
                }

                for (m = 0; m < DIM; m++)
                {
                    frameorder[m] += 0.5 * (3.0 * cossum[m] - 1.0);
                }

                if (bSliced)
                {
                    /* get average coordinate in box length for slicing,
                       determine which slice atom is in, increase count for that
                       slice. slFrameorder and slOrder are reals, not
                       rvecs. Only the component [axis] of the order tensor is
                       kept, until I find it necessary to know the others too
                     */

                    z1    = x1[indexGroups[i - 1].particleIndices[j]][axis];
                    z2    = x1[indexGroups[i + 1].particleIndices[j]][axis];
                    z_ave = 0.5 * (z1 + z2);
                    slice = static_cast<int>((static_cast<real>(nslices) * z_ave) / box[axis][axis]);
                    while (slice < 0)
                    {
                        slice += static_cast<real>(nslices);
                    }
                    slice = slice % nslices;

                    slCount[slice]++; /* determine slice, increase count */

                    slFrameorder[slice] += 0.5 * (3 * cossum[axis] - 1);
                }
                else if (permolecule)
                {
                    /*  store per-molecule order parameter
                     *  To just track single-axis order: (*slOrder)[j][i] += 0.5 * (3 *
                     * iprod(cossum,direction) - 1); following is for Scd order: */
                    (*slOrder)[j][i] += -1
                                        * (1.0 / 3.0 * (3 * cossum[XX] - 1)
                                           + 1.0 / 3.0 * 0.5 * (3.0 * cossum[YY] - 1));
                }
                if (distcalc)
                {
                    if (radial)
                    {
                        /* bin order parameter by arc distance from reference group*/
                        arcdist = gmx_angle(displacementFromReference, direction);
                        (*distvals)[j][i] += arcdist;
                    }
                    else if (i == 1)
                    {
                        /* Want minimum lateral distance to first group calculated */
                        tmpdist = trace(box); /* should be max value */
                        for (k = 0; k < distsize; k++)
                        {
                            rvec displacement;
                            pbc_dx(&pbc, x1[distidx[k]], x1[indexGroups[i].particleIndices[j]], displacement);
                            /* at the moment, just remove displacement[axis] */
                            displacement[axis] = 0;
                            tmpdist            = std::min(tmpdist, norm2(displacement));
                        }
                        // fprintf(stderr, "Min dist %f; trace %f\n", tmpdist, trace(box));
                        (*distvals)[j][i] += std::sqrt(tmpdist);
                    }
                }
            } /* end loop j, over all atoms in group */

            for (m = 0; m < DIM; m++)
            {
                (*order)[i][m] += (frameorder[m] / static_cast<real>(size));
            }

            if (!permolecule)
            { /*Skip following if doing per-molecule*/
                for (k = 0; k < nslices; k++)
                {
                    if (slCount[k]) /* if no elements, nothing has to be added */
                    {
                        (*slOrder)[k][i] += slFrameorder[k] / static_cast<real>(slCount[k]);
                        slFrameorder[k] = 0;
                        slCount[k]      = 0;
                    }
                }
            } /* end loop i, over all groups in indexfile */
        }
        nr_frames++;

    } while (read_next_x(oenv, status, &t, x0, box));
    /*********** done with status file **********/

    fprintf(stderr, "\nRead trajectory. Printing parameters to file\n");
    gmx_rmpbc_done(gpbc);

    /* average over frames */
    for (i = 1; i < ngrps - 1; i++)
    {
        svmul(1.0 / nr_frames, (*order)[i], (*order)[i]);
        fprintf(stderr,
                "Atom %d Tensor: x=%g , y=%g, z=%g\n",
                i,
                (*order)[i][XX],
                (*order)[i][YY],
                (*order)[i][ZZ]);
        if (bSliced || permolecule)
        {
            for (k = 0; k < nslices; k++)
            {
                (*slOrder)[k][i] /= nr_frames;
            }
        }
        if (distcalc)
        {
            for (k = 0; k < nslices; k++)
            {
                (*distvals)[k][i] /= nr_frames;
            }
        }
    }

    sfree(x0); /* free memory used by coordinate arrays */
    sfree(x1);
    if (comidx != nullptr)
    {
        sfree(comidx);
    }
    if (distidx != nullptr)
    {
        sfree(distidx);
    }
    if (grpname != nullptr)
    {
        sfree(grpname);
    }
}


static void order_plot(rvec                    order[],
                       real*                   slOrder[],
                       const char*             afile,
                       const char*             bfile,
                       const char*             cfile,
                       int                     ngrps,
                       int                     nslices,
                       real                    slWidth,
                       gmx_bool                bSzonly,
                       gmx_bool                permolecule,
                       real**                  distvals,
                       const gmx_output_env_t* oenv)
{
    FILE *ord, *slOrd; /* xvgr files with order parameters  */
    int   atom, slice; /* atom corresponding to order para.*/
    char  buf[256];    /* for xvgr title */
    real  S;           /* order parameter averaged over all atoms */

    if (permolecule)
    {
        sprintf(buf, "Scd order parameters");
        ord = xvgropen(afile, buf, "Atom", "S", oenv);
        sprintf(buf, "Orderparameters per atom per slice");
        slOrd = xvgropen(bfile, buf, "Molecule", "S", oenv);
        for (atom = 1; atom < ngrps - 1; atom++)
        {
            fprintf(ord,
                    "%12d   %12g\n",
                    atom,
                    -1.0 * (2.0 / 3.0 * order[atom][XX] + 1.0 / 3.0 * order[atom][YY]));
        }

        for (slice = 0; slice < nslices; slice++)
        {
            fprintf(slOrd, "%12d\t", slice);
            if (distvals)
            {
                fprintf(slOrd, "%12g\t", distvals[slice][1]); /*use distance value at second carbon*/
            }
            for (atom = 1; atom < ngrps - 1; atom++)
            {
                fprintf(slOrd, "%12g\t", slOrder[slice][atom]);
            }
            fprintf(slOrd, "\n");
        }
    }
    else if (bSzonly)
    {
        sprintf(buf, "Orderparameters Sz per atom");
        ord = xvgropen(afile, buf, "Atom", "S", oenv);
        fprintf(stderr, "ngrps = %d, nslices = %d", ngrps, nslices);

        sprintf(buf, "Orderparameters per atom per slice");
        slOrd = xvgropen(bfile, buf, "Slice", "S", oenv);

        for (atom = 1; atom < ngrps - 1; atom++)
        {
            fprintf(ord, "%12d       %12g\n", atom, order[atom][ZZ]);
        }

        for (slice = 0; slice < nslices; slice++)
        {
            S = 0;
            for (atom = 1; atom < ngrps - 1; atom++)
            {
                S += slOrder[slice][atom];
            }
            fprintf(slOrd, "%12g     %12g\n", static_cast<real>(slice) * slWidth, S / static_cast<real>(atom));
        }
    }
    else
    {
        sprintf(buf, "Order tensor diagonal elements");
        ord = xvgropen(afile, buf, "Atom", "S", oenv);
        sprintf(buf, "Deuterium order parameters");
        slOrd = xvgropen(cfile, buf, "Atom", "Scd", oenv);

        for (atom = 1; atom < ngrps - 1; atom++)
        {
            fprintf(ord, "%12d   %12g   %12g   %12g\n", atom, order[atom][XX], order[atom][YY], order[atom][ZZ]);
            fprintf(slOrd,
                    "%12d   %12g\n",
                    atom,
                    -1.0 * (2.0 / 3.0 * order[atom][XX] + 1.0 / 3.0 * order[atom][YY]));
        }
    }
    xvgrclose(ord);
    xvgrclose(slOrd);
}

static void write_bfactors(t_filenm*                       fnm,
                           int                             nfile,
                           gmx::ArrayRef<const IndexGroup> indexGroups,
                           int                             nslices,
                           real**                          order,
                           const t_topology*               top,
                           real**                          distvals,
                           gmx_output_env_t*               oenv)
{
    /*function to write order parameters as B factors in PDB file using
          first frame of trajectory*/
    t_trxstatus* status;
    t_trxframe   fr, frout;
    t_atoms      useatoms;
    int          i, j, ctr, nout;

    /* we don't have an order parameter for the first or last atom in each chain */
    const int ngrps = gmx::ssize(indexGroups) - 2;

    nout = nslices * ngrps;
    read_first_frame(oenv, &status, ftp2fn(efTRX, nfile, fnm), &fr, TRX_NEED_X);

    close_trx(status);
    frout        = fr;
    frout.natoms = nout;
    frout.bF     = FALSE;
    frout.bV     = FALSE;
    frout.x      = nullptr;
    snew(frout.x, nout);

    init_t_atoms(&useatoms, nout, TRUE);
    useatoms.nr = nout;

    /*initialize PDBinfo*/
    for (i = 0; i < useatoms.nr; ++i)
    {
        useatoms.pdbinfo[i].type         = PdbRecordType::Atom;
        useatoms.pdbinfo[i].occup        = 0.0;
        useatoms.pdbinfo[i].bfac         = 0.0;
        useatoms.pdbinfo[i].bAnisotropic = FALSE;
    }

    for (j = 0, ctr = 0; j < nslices; j++)
    {
        for (i = 0; i < ngrps; i++, ctr++)
        {
            /*iterate along each chain*/
            useatoms.pdbinfo[ctr].bfac = order[j][i + 1];
            if (distvals)
            {
                useatoms.pdbinfo[ctr].occup = distvals[j][i + 1];
            }
            const int atomIndex = indexGroups[i + 1].particleIndices[j];
            copy_rvec(fr.x[atomIndex], frout.x[ctr]);
            useatoms.atomname[ctr] = top->atoms.atomname[atomIndex];
            useatoms.atom[ctr]     = top->atoms.atom[atomIndex];
            useatoms.nres          = std::max(useatoms.nres, useatoms.atom[ctr].resind + 1);
            useatoms.resinfo[useatoms.atom[ctr].resind] =
                    top->atoms.resinfo[useatoms.atom[ctr].resind]; /*copy resinfo*/
        }
    }

    write_sto_conf(
            opt2fn("-ob", nfile, fnm), "Order parameters", &useatoms, frout.x, nullptr, frout.pbcType, frout.box);

    sfree(frout.x);
    done_atom(&useatoms);
}

int gmx_order(int argc, char* argv[])
{
    const char* desc[] = {
        "[THISMODULE] computes the order parameter per atom for carbon tails. For atom i the",
        "vector i-1, i+1 is used together with an axis. ",
        "The index file should contain only the groups to be used for calculations,",
        "with each group of equivalent carbons along the relevant acyl chain in its own",
        "group. There should not be any generic groups (like System, Protein) in the index",
        "file to avoid confusing the program (this is not relevant to tetrahedral order",
        "parameters however, which only work for water anyway).[PAR]",
        "[THISMODULE] can also give all",
        "diagonal elements of the order tensor and even calculate the deuterium",
        "order parameter Scd (default). If the option [TT]-szonly[tt] is given, only one",
        "order tensor component (specified by the [TT]-d[tt] option) is given and the",
        "order parameter per slice is calculated as well. If [TT]-szonly[tt] is not",
        "selected, all diagonal elements and the deuterium order parameter is",
        "given.[PAR]",
        "The tetrahedrality order parameters can be determined",
        "around an atom. Both angle an distance order parameters are calculated. See",
        "P.-L. Chau and A.J. Hardwick, Mol. Phys., 93, (1998), 511-518.",
        "for more details."
    };

    const char* bugs[] = {
        "This tool only works for saturated carbons and united atom force fields.",
        "For anything else, it is highly recommended to use a different analysis method!",
        "The option [TT]-unsat[tt] claimed to do analysis for unsaturated carbons",
        "this but hasn't worked ever since it was added and has thus been removed."
    };

    static int         nslices       = 1;     /* nr of slices defined       */
    static gmx_bool    bSzonly       = FALSE; /* True if only Sz is wanted  */
    static bool        bUnsatRemoved = false; /* Removed because it doesn't work. */
    static const char* normal_axis[] = { nullptr, "z", "x", "y", nullptr };
    static gmx_bool    permolecule   = FALSE; /*compute on a per-molecule basis */
    static gmx_bool    radial        = FALSE; /*compute a radial membrane normal */
    static gmx_bool    distcalc      = FALSE; /*calculate distance from a reference group */
    t_pargs            pa[]          = {
        { "-d", FALSE, etENUM, { normal_axis }, "Direction of the normal on the membrane" },
        { "-sl",
          FALSE,
          etINT,
          { &nslices },
          "Calculate order parameter as function of box length, dividing the box"
          " into this number of slices." },
        { "-szonly",
          FALSE,
          etBOOL,
          { &bSzonly },
          "Only give Sz element of order tensor. (axis can be specified with [TT]-d[tt])" },
        { "-unsat",
          FALSE,
          etBOOL,
          { &bUnsatRemoved },
          "HIDDENThis option has been removed as it didn't ever properly work." },
        { "-permolecule",
          FALSE,
          etBOOL,
          { &permolecule },
          "Compute per-molecule Scd order parameters" },
        { "-radial", FALSE, etBOOL, { &radial }, "Compute a radial membrane normal" },
        { "-calcdist", FALSE, etBOOL, { &distcalc }, "Compute distance from a reference" },
    };

    rvec*       order;         /* order par. for each atom   */
    real**      slOrder;       /* same, per slice            */
    real        slWidth = 0.0; /* width of a slice           */
    int         i, axis = 0;   /* normal axis                */
    t_topology* top;           /* topology                   */
    PbcType     pbcType;       /* type of periodic boundary conditions */

    t_filenm fnm[] = {
        /* files for gmx order    */
        { efTRX, "-f", nullptr, ffREAD },     /* trajectory file              */
        { efNDX, "-n", nullptr, ffREAD },     /* index file           */
        { efNDX, "-nr", nullptr, ffOPTRD },   /* index for radial axis calculation */
        { efTPR, nullptr, nullptr, ffREAD },  /* topology file                */
        { efXVG, "-o", "order", ffWRITE },    /* xvgr output file     */
        { efXVG, "-od", "deuter", ffWRITE },  /* xvgr output file           */
        { efPDB, "-ob", nullptr, ffOPTWR },   /* write Scd as B factors to PDB if permolecule   */
        { efXVG, "-os", "sliced", ffWRITE },  /* xvgr output file           */
        { efXVG, "-Sg", "sg-ang", ffOPTWR },  /* xvgr output file           */
        { efXVG, "-Sk", "sk-dist", ffOPTWR }, /* xvgr output file           */
        { efXVG, "-Sgsl", "sg-ang-slice", ffOPTWR },  /* xvgr output file           */
        { efXVG, "-Sksl", "sk-dist-slice", ffOPTWR }, /* xvgr output file           */
    };
    gmx_bool bSliced = FALSE; /* True if box is sliced      */
#define NFILE asize(fnm)
    real**            distvals = nullptr;
    const char *      sgfnm, *skfnm, *ndxfnm, *tpsfnm, *trxfnm;
    gmx_output_env_t* oenv;

    if (!parse_common_args(
                &argc, argv, PCA_CAN_VIEW | PCA_CAN_TIME, NFILE, fnm, asize(pa), pa, asize(desc), desc, asize(bugs), bugs, &oenv))
    {
        return 0;
    }
    if (nslices < 1)
    {
        gmx_fatal(FARGS, "Can not have nslices < 1");
    }
    sgfnm  = opt2fn_null("-Sg", NFILE, fnm);
    skfnm  = opt2fn_null("-Sk", NFILE, fnm);
    ndxfnm = opt2fn_null("-n", NFILE, fnm);
    tpsfnm = ftp2fn(efTPR, NFILE, fnm);
    trxfnm = ftp2fn(efTRX, NFILE, fnm);

    /* Calculate axis */
    GMX_RELEASE_ASSERT(normal_axis[0] != nullptr, "Options inconsistency; normal_axis[0] is NULL");
    if (std::strcmp(normal_axis[0], "x") == 0)
    {
        axis = XX;
    }
    else if (std::strcmp(normal_axis[0], "y") == 0)
    {
        axis = YY;
    }
    else if (std::strcmp(normal_axis[0], "z") == 0)
    {
        axis = ZZ;
    }
    else
    {
        gmx_fatal(FARGS, "Invalid axis, use x, y or z");
    }

    switch (axis)
    {
        case 0: fprintf(stderr, "Taking x axis as normal to the membrane\n"); break;
        case 1: fprintf(stderr, "Taking y axis as normal to the membrane\n"); break;
        case 2: fprintf(stderr, "Taking z axis as normal to the membrane\n"); break;
    }

    /* tetraheder order parameter */
    if (skfnm || sgfnm)
    {
        /* If either of theoptions is set we compute both */
        sgfnm = opt2fn("-Sg", NFILE, fnm);
        skfnm = opt2fn("-Sk", NFILE, fnm);
        calc_tetra_order_parm(ndxfnm,
                              tpsfnm,
                              trxfnm,
                              sgfnm,
                              skfnm,
                              nslices,
                              axis,
                              opt2fn("-Sgsl", NFILE, fnm),
                              opt2fn("-Sksl", NFILE, fnm),
                              oenv);
        /* view xvgr files */
        do_view(oenv, opt2fn("-Sg", NFILE, fnm), nullptr);
        do_view(oenv, opt2fn("-Sk", NFILE, fnm), nullptr);
        if (nslices > 1)
        {
            do_view(oenv, opt2fn("-Sgsl", NFILE, fnm), nullptr);
            do_view(oenv, opt2fn("-Sksl", NFILE, fnm), nullptr);
        }
    }
    else
    {
        /* tail order parameter */

        if (nslices > 1)
        {
            bSliced = TRUE;
            fprintf(stderr, "Dividing box in %d slices.\n\n", nslices);
        }

        if (bSzonly)
        {
            fprintf(stderr, "Only calculating Sz\n");
        }
        if (bUnsatRemoved)
        {
            gmx_fatal(FARGS,
                      "The option to process unsaturated carbons has been removed because it never "
                      "properly worked. Please use a different tool to analyse your data!\n");
        }

        top = read_top(ftp2fn(efTPR, NFILE, fnm), &pbcType); /* read topology file */

        auto indexGroups = init_index(ftp2fn(efNDX, NFILE, fnm));

        if (permolecule)
        {
            nslices = gmx::ssize(indexGroups[0].particleIndices);
            fprintf(stderr, "Calculating Scd order parameters for each of %d molecules\n", nslices);
        }

        if (radial)
        {
            fprintf(stderr, "Calculating radial distances\n");
            if (!permolecule)
            {
                gmx_fatal(FARGS, "Cannot yet output radial distances without permolecule\n");
            }
        }

        /* show atomtypes, to check if index file is correct */
        print_types(indexGroups, top);

        calc_order(ftp2fn(efTRX, NFILE, fnm),
                   indexGroups,
                   &order,
                   &slOrder,
                   &slWidth,
                   nslices,
                   bSliced,
                   top,
                   pbcType,
                   axis,
                   permolecule,
                   radial,
                   distcalc,
                   opt2fn_null("-nr", NFILE, fnm),
                   &distvals,
                   oenv);

        if (radial)
        {
            /* don't print the last group--was used for center-of-mass determination */
            indexGroups.pop_back();
        }
        order_plot(order,
                   slOrder,
                   opt2fn("-o", NFILE, fnm),
                   opt2fn("-os", NFILE, fnm),
                   opt2fn("-od", NFILE, fnm),
                   gmx::ssize(indexGroups),
                   nslices,
                   slWidth,
                   bSzonly,
                   permolecule,
                   distvals,
                   oenv);

        if (opt2bSet("-ob", NFILE, fnm))
        {
            if (!permolecule)
            {
                fprintf(stderr,
                        "Won't write B-factors with averaged order parameters; use -permolecule\n");
            }
            else
            {
                write_bfactors(fnm, NFILE, indexGroups, nslices, slOrder, top, distvals, oenv);
            }
        }


        do_view(oenv, opt2fn("-o", NFILE, fnm), nullptr);  /* view xvgr file */
        do_view(oenv, opt2fn("-os", NFILE, fnm), nullptr); /* view xvgr file */
        do_view(oenv, opt2fn("-od", NFILE, fnm), nullptr); /* view xvgr file */
    }

    if (distvals != nullptr)
    {
        for (i = 0; i < nslices; ++i)
        {
            sfree(distvals[i]);
        }
        sfree(distvals);
    }

    return 0;
}
