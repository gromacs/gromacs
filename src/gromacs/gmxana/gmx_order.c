/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015, by the GROMACS development team, led by
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
#include "gmxpre.h"

#include <math.h>
#include <string.h>

#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/gmxana/cmat.h"
#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/gmxana/gstat.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/legacyheaders/viewit.h"
#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pbcutil/rmpbc.h"
#include "gromacs/topology/index.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
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

static void find_nearest_neighbours(int ePBC,
                                    int natoms, matrix box,
                                    rvec x[], int maxidx, atom_id index[],
                                    real *sgmean, real *skmean,
                                    int nslice, int slice_dim,
                                    real sgslice[], real skslice[],
                                    gmx_rmpbc_t gpbc)
{
    FILE    *fpoutdist;
    char     fnsgdist[32];
    int      ix, jx, nsgbin, *sgbin;
    int      i1, i2, i, ibin, j, k, l, n, *nn[4];
    rvec     dx, dx1, dx2, rj, rk, urk, urj;
    real     cost, cost2, *sgmol, *skmol, rmean, rmean2, r2, box2, *r_nn[4];
    t_pbc    pbc;
    t_mat   *dmat;
    t_dist  *d;
    int      m1, mm, sl_index;
    int    **nnb, *sl_count;
    real     onethird = 1.0/3.0;
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
    set_pbc(&pbc, ePBC, box);

    gmx_rmpbc(gpbc, natoms, box, x);

    nsgbin = 1 + 1/0.0005;
    snew(sgbin, nsgbin);

    *sgmean = 0.0;
    *skmean = 0.0;
    l       = 0;
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
                r_nn[3][i] = r_nn[2][i]; nn[3][i] = nn[2][i];
                r_nn[2][i] = r_nn[1][i]; nn[2][i] = nn[1][i];
                r_nn[1][i] = r_nn[0][i]; nn[1][i] = nn[0][i];
                r_nn[0][i] = r2;         nn[0][i] = j;
            }
            else if (r2 < r_nn[1][i])
            {
                r_nn[3][i] = r_nn[2][i]; nn[3][i] = nn[2][i];
                r_nn[2][i] = r_nn[1][i]; nn[2][i] = nn[1][i];
                r_nn[1][i] = r2;         nn[1][i] = j;
            }
            else if (r2 < r_nn[2][i])
            {
                r_nn[3][i] = r_nn[2][i]; nn[3][i] = nn[2][i];
                r_nn[2][i] = r2;         nn[2][i] = j;
            }
            else if (r2 < r_nn[3][i])
            {
                r_nn[3][i] = r2;         nn[3][i] = j;
            }
        }


        /* calculate mean distance between nearest neighbours */
        rmean = 0;
        for (j = 0; (j < 4); j++)
        {
            r_nn[j][i] = sqrt(r_nn[j][i]);
            rmean     += r_nn[j][i];
        }
        rmean /= 4;

        n        = 0;
        sgmol[i] = 0.0;
        skmol[i] = 0.0;

        /* Chau1998a eqn 3 */
        /* angular part tetrahedrality order parameter per atom */
        for (j = 0; (j < 3); j++)
        {
            for (k = j+1; (k < 4); k++)
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
                ibin = nsgbin * cost2;
                if (ibin < nsgbin)
                {
                    sgbin[ibin]++;
                }
                /* printf("%d %d %f %d %d\n", j, k, cost * cost, ibin, sgbin[ibin]);*/
                l++;
                n++;
            }
        }

        /* normalize sgmol between 0.0 and 1.0 */
        sgmol[i] = 3*sgmol[i]/32;
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
        sl_index           = gmx_nint((1+x[i][slice_dim]/box[slice_dim][slice_dim])*nslice) % nslice;
        sgslice[sl_index] += sgmol[i];
        skslice[sl_index] += skmol[i];
        sl_count[sl_index]++;
    } /* loop over entries in index file */

    *sgmean /= maxidx;
    *skmean /= maxidx;

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


static void calc_tetra_order_parm(const char *fnNDX, const char *fnTPS,
                                  const char *fnTRX, const char *sgfn,
                                  const char *skfn,
                                  int nslice, int slice_dim,
                                  const char *sgslfn, const char *skslfn,
                                  const output_env_t oenv)
{
    FILE        *fpsg = NULL, *fpsk = NULL;
    t_topology   top;
    int          ePBC;
    char         title[STRLEN], fn[STRLEN], subtitle[STRLEN];
    t_trxstatus *status;
    int          natoms;
    real         t;
    rvec        *xtop, *x;
    matrix       box;
    real         sg, sk;
    atom_id    **index;
    char       **grpname;
    int          i, *isize, ng, nframes;
    real        *sg_slice, *sg_slice_tot, *sk_slice, *sk_slice_tot;
    gmx_rmpbc_t  gpbc = NULL;


    read_tps_conf(fnTPS, title, &top, &ePBC, &xtop, NULL, box, FALSE);

    snew(sg_slice, nslice);
    snew(sk_slice, nslice);
    snew(sg_slice_tot, nslice);
    snew(sk_slice_tot, nslice);
    ng = 1;
    /* get index groups */
    printf("Select the group that contains the atoms you want to use for the tetrahedrality order parameter calculation:\n");
    snew(grpname, ng);
    snew(index, ng);
    snew(isize, ng);
    get_index(&top.atoms, fnNDX, ng, isize, index, grpname);

    /* Analyze trajectory */
    natoms = read_first_x(oenv, &status, fnTRX, &t, &x, box);
    if (natoms > top.atoms.nr)
    {
        gmx_fatal(FARGS, "Topology (%d atoms) does not match trajectory (%d atoms)",
                  top.atoms.nr, natoms);
    }
    check_index(NULL, ng, index[0], NULL, natoms);

    fpsg = xvgropen(sgfn, "S\\sg\\N Angle Order Parameter", "Time (ps)", "S\\sg\\N",
                    oenv);
    fpsk = xvgropen(skfn, "S\\sk\\N Distance Order Parameter", "Time (ps)", "S\\sk\\N",
                    oenv);

    /* loop over frames */
    gpbc    = gmx_rmpbc_init(&top.idef, ePBC, natoms);
    nframes = 0;
    do
    {
        find_nearest_neighbours(ePBC, natoms, box, x, isize[0], index[0],
                                &sg, &sk, nslice, slice_dim, sg_slice, sk_slice, gpbc);
        for (i = 0; (i < nslice); i++)
        {
            sg_slice_tot[i] += sg_slice[i];
            sk_slice_tot[i] += sk_slice[i];
        }
        fprintf(fpsg, "%f %f\n", t, sg);
        fprintf(fpsk, "%f %f\n", t, sk);
        nframes++;
    }
    while (read_next_x(oenv, status, &t, x, box));
    close_trj(status);
    gmx_rmpbc_done(gpbc);

    sfree(grpname);
    sfree(index);
    sfree(isize);

    xvgrclose(fpsg);
    xvgrclose(fpsk);

    fpsg = xvgropen(sgslfn,
                    "S\\sg\\N Angle Order Parameter / Slab", "(nm)", "S\\sg\\N",
                    oenv);
    fpsk = xvgropen(skslfn,
                    "S\\sk\\N Distance Order Parameter / Slab", "(nm)", "S\\sk\\N",
                    oenv);
    for (i = 0; (i < nslice); i++)
    {
        fprintf(fpsg, "%10g  %10g\n", (i+0.5)*box[slice_dim][slice_dim]/nslice,
                sg_slice_tot[i]/nframes);
        fprintf(fpsk, "%10g  %10g\n", (i+0.5)*box[slice_dim][slice_dim]/nslice,
                sk_slice_tot[i]/nframes);
    }
    xvgrclose(fpsg);
    xvgrclose(fpsk);
}


/* Print name of first atom in all groups in index file */
static void print_types(atom_id index[], atom_id a[], int ngrps,
                        char *groups[], t_topology *top)
{
    int i;

    fprintf(stderr, "Using following groups: \n");
    for (i = 0; i < ngrps; i++)
    {
        fprintf(stderr, "Groupname: %s First atomname: %s First atomnr %d\n",
                groups[i], *(top->atoms.atomname[a[index[i]]]), a[index[i]]);
    }
    fprintf(stderr, "\n");
}

static void check_length(real length, int a, int b)
{
    if (length > 0.3)
    {
        fprintf(stderr, "WARNING: distance between atoms %d and "
                "%d > 0.3 nm (%f). Index file might be corrupt.\n",
                a, b, length);
    }
}

void calc_order(const char *fn, atom_id *index, atom_id *a, rvec **order,
                real ***slOrder, real *slWidth, int nslices, gmx_bool bSliced,
                gmx_bool bUnsat, t_topology *top, int ePBC, int ngrps, int axis,
                gmx_bool permolecule, gmx_bool radial, gmx_bool distcalc, const char *radfn,
                real ***distvals,
                const output_env_t oenv)
{
    /* if permolecule = TRUE, order parameters will be calculed per molecule
     * and stored in slOrder with #slices = # molecules */
    rvec *x0,                                    /* coordinates with pbc                           */
    *x1,                                         /* coordinates without pbc                        */
          dist;                                  /* vector between two atoms                       */
    matrix       box;                            /* box (3x3)                                      */
    t_trxstatus *status;
    rvec         cossum,                         /* sum of vector angles for three axes            */
                 Sx, Sy, Sz,                     /* the three molecular axes                       */
                 tmp1, tmp2,                     /* temp. rvecs for calculating dot products       */
                 frameorder;                     /* order parameters for one frame                 */
    real *slFrameorder;                          /* order parameter for one frame, per slice      */
    real  length,                                /* total distance between two atoms               */
          t,                                     /* time from trajectory                           */
          z_ave, z1, z2;                         /* average z, used to det. which slice atom is in */
    int natoms,                                  /* nr. atoms in trj                               */
        nr_tails,                                /* nr tails, to check if index file is correct    */
        size = 0,                                /* nr. of atoms in group. same as nr_tails        */
        i, j, m, k, l, teller = 0,
        slice,                                   /* current slice number                           */
        nr_frames = 0;
    int         *slCount;                        /* nr. of atoms in one slice                      */
    real         dbangle                = 0,     /* angle between double bond and  axis            */
                 sdbangle               = 0;     /* sum of these angles                            */
    gmx_bool     use_unitvector         = FALSE; /* use a specified unit vector instead of axis to specify unit normal*/
    rvec         direction, com, dref, dvec;
    int          comsize, distsize;
    atom_id     *comidx  = NULL, *distidx = NULL;
    char        *grpname = NULL;
    t_pbc        pbc;
    real         arcdist, tmpdist;
    gmx_rmpbc_t  gpbc = NULL;

    /* PBC added for center-of-mass vector*/
    /* Initiate the pbc structure */
    memset(&pbc, 0, sizeof(pbc));

    if ((natoms = read_first_x(oenv, &status, fn, &t, &x0, box)) == 0)
    {
        gmx_fatal(FARGS, "Could not read coordinates from statusfile\n");
    }

    nr_tails = index[1] - index[0];
    fprintf(stderr, "Number of elements in first group: %d\n", nr_tails);
    /* take first group as standard. Not rocksolid, but might catch error in index*/

    if (permolecule)
    {
        nslices = nr_tails;
        bSliced = FALSE; /*force slices off */
        fprintf(stderr, "Calculating order parameters for each of %d molecules\n",
                nslices);
    }

    if (radial)
    {
        use_unitvector = TRUE;
        fprintf(stderr, "Select an index group to calculate the radial membrane normal\n");
        get_index(&top->atoms, radfn, 1, &comsize, &comidx, &grpname);
    }
    if (distcalc)
    {
        if (grpname != NULL)
        {
            sfree(grpname);
        }
        fprintf(stderr, "Select an index group to use as distance reference\n");
        get_index(&top->atoms, radfn, 1, &distsize, &distidx, &grpname);
        bSliced = FALSE; /*force slices off*/
    }

    if (use_unitvector && bSliced)
    {
        fprintf(stderr, "Warning:  slicing and specified unit vectors are not currently compatible\n");
    }

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
        *slWidth = box[axis][axis]/nslices;
        fprintf(stderr, "Box divided in %d slices. Initial width of slice: %f\n",
                nslices, *slWidth);
    }


#if 0
    nr_tails = index[1] - index[0];
    fprintf(stderr, "Number of elements in first group: %d\n", nr_tails);
    /* take first group as standard. Not rocksolid, but might catch error
       in index*/
#endif

    teller = 0;

    gpbc = gmx_rmpbc_init(&top->idef, ePBC, natoms);
    /*********** Start processing trajectory ***********/
    do
    {
        if (bSliced)
        {
            *slWidth = box[axis][axis]/nslices;
        }
        teller++;

        set_pbc(&pbc, ePBC, box);
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
            com[XX] = 0.0; com[YY] = 0.0; com[ZZ] = 0.0;
            for (j = 0; j < comsize; j++)
            {
                rvec_inc(com, x1[comidx[j]]);
            }
            svmul(1.0/comsize, com, com);
        }
        if (distcalc)
        {
            dref[XX] = 0.0; dref[YY] = 0.0; dref[ZZ] = 0.0;
            for (j = 0; j < distsize; j++)
            {
                rvec_inc(dist, x1[distidx[j]]);
            }
            svmul(1.0/distsize, dref, dref);
            if (radial)
            {
                pbc_dx(&pbc, dref, com, dvec);
                unitv(dvec, dvec);
            }
        }

        for (i = 1; i < ngrps - 1; i++)
        {
            clear_rvec(frameorder);

            size = index[i+1] - index[i];
            if (size != nr_tails)
            {
                gmx_fatal(FARGS, "grp %d does not have same number of"
                          " elements as grp 1\n", i);
            }

            for (j = 0; j < size; j++)
            {
                if (radial)
                /*create unit vector*/
                {
                    pbc_dx(&pbc, x1[a[index[i]+j]], com, direction);
                    unitv(direction, direction);
                    /*DEBUG*/
                    /*if (j==0)
                        fprintf(stderr,"X %f %f %f\tcom %f %f %f\tdirection %f %f %f\n",x1[a[index[i]+j]][0],x1[a[index[i]+j]][1],x1[a[index[i]+j]][2],com[0],com[1],com[2],
                            direction[0],direction[1],direction[2]);*/
                }

                if (bUnsat)
                {
                    /* Using convention for unsaturated carbons */
                    /* first get Sz, the vector from Cn to Cn+1 */
                    rvec_sub(x1[a[index[i+1]+j]], x1[a[index[i]+j]], dist);
                    length = norm(dist);
                    check_length(length, a[index[i]+j], a[index[i+1]+j]);
                    svmul(1/length, dist, Sz);

                    /* this is actually the cosine of the angle between the double bond
                       and axis, because Sz is normalized and the two other components of
                       the axis on the bilayer are zero */
                    if (use_unitvector)
                    {
                        sdbangle += gmx_angle(direction, Sz); /*this can probably be optimized*/
                    }
                    else
                    {
                        sdbangle += acos(Sz[axis]);
                    }
                }
                else
                {
                    /* get vector dist(Cn-1,Cn+1) for tail atoms */
                    rvec_sub(x1[a[index[i+1]+j]], x1[a[index[i-1]+j]], dist);
                    length = norm(dist); /* determine distance between two atoms */
                    check_length(length, a[index[i-1]+j], a[index[i+1]+j]);

                    svmul(1/length, dist, Sz);
                    /* Sz is now the molecular axis Sz, normalized and all that */
                }

                /* now get Sx. Sx is normal to the plane of Cn-1, Cn and Cn+1 so
                   we can use the outer product of Cn-1->Cn and Cn+1->Cn, I hope */
                rvec_sub(x1[a[index[i+1]+j]], x1[a[index[i]+j]], tmp1);
                rvec_sub(x1[a[index[i-1]+j]], x1[a[index[i]+j]], tmp2);
                cprod(tmp1, tmp2, Sx);
                svmul(1/norm(Sx), Sx, Sx);

                /* now we can get Sy from the outer product of Sx and Sz   */
                cprod(Sz, Sx, Sy);
                svmul(1/norm(Sy), Sy, Sy);

                /* the square of cosine of the angle between dist and the axis.
                   Using the innerproduct, but two of the three elements are zero
                   Determine the sum of the orderparameter of all atoms in group
                 */
                if (use_unitvector)
                {
                    cossum[XX] = sqr(iprod(Sx, direction)); /* this is allowed, since Sa is normalized */
                    cossum[YY] = sqr(iprod(Sy, direction));
                    cossum[ZZ] = sqr(iprod(Sz, direction));
                }
                else
                {
                    cossum[XX] = sqr(Sx[axis]); /* this is allowed, since Sa is normalized */
                    cossum[YY] = sqr(Sy[axis]);
                    cossum[ZZ] = sqr(Sz[axis]);
                }

                for (m = 0; m < DIM; m++)
                {
                    frameorder[m] += 0.5 * (3 * cossum[m] - 1);
                }

                if (bSliced)
                {
                    /* get average coordinate in box length for slicing,
                       determine which slice atom is in, increase count for that
                       slice. slFrameorder and slOrder are reals, not
                       rvecs. Only the component [axis] of the order tensor is
                       kept, until I find it necessary to know the others too
                     */

                    z1    = x1[a[index[i-1]+j]][axis];
                    z2    = x1[a[index[i+1]+j]][axis];
                    z_ave = 0.5 * (z1 + z2);
                    if (z_ave < 0)
                    {
                        z_ave += box[axis][axis];
                    }
                    if (z_ave > box[axis][axis])
                    {
                        z_ave -= box[axis][axis];
                    }

                    slice  = (int)(0.5 + (z_ave / (*slWidth))) - 1;
                    slCount[slice]++;     /* determine slice, increase count */

                    slFrameorder[slice] += 0.5 * (3 * cossum[axis] - 1);
                }
                else if (permolecule)
                {
                    /*  store per-molecule order parameter
                     *  To just track single-axis order: (*slOrder)[j][i] += 0.5 * (3 * iprod(cossum,direction) - 1);
                     *  following is for Scd order: */
                    (*slOrder)[j][i] += -1* (0.3333 * (3 * cossum[XX] - 1) + 0.3333 * 0.5 * (3 * cossum[YY] - 1));
                }
                if (distcalc)
                {
                    if (radial)
                    {
                        /* bin order parameter by arc distance from reference group*/
                        arcdist            = gmx_angle(dvec, direction);
                        (*distvals)[j][i] += arcdist;
                    }
                    else if (i == 1)
                    {
                        /* Want minimum lateral distance to first group calculated */
                        tmpdist = trace(box);  /* should be max value */
                        for (k = 0; k < distsize; k++)
                        {
                            pbc_dx(&pbc, x1[distidx[k]], x1[a[index[i]+j]], dvec);
                            /* at the moment, just remove dvec[axis] */
                            dvec[axis] = 0;
                            tmpdist    = min(tmpdist, norm2(dvec));
                        }
                        //fprintf(stderr, "Min dist %f; trace %f\n", tmpdist, trace(box));
                        (*distvals)[j][i] += sqrt(tmpdist);
                    }
                }
            } /* end loop j, over all atoms in group */

            for (m = 0; m < DIM; m++)
            {
                (*order)[i][m] += (frameorder[m]/size);
            }

            if (!permolecule)
            {   /*Skip following if doing per-molecule*/
                for (k = 0; k < nslices; k++)
                {
                    if (slCount[k]) /* if no elements, nothing has to be added */
                    {
                        (*slOrder)[k][i] += slFrameorder[k]/slCount[k];
                        slFrameorder[k]   = 0; slCount[k] = 0;
                    }
                }
            } /* end loop i, over all groups in indexfile */
        }
        nr_frames++;

    }
    while (read_next_x(oenv, status, &t, x0, box));
    /*********** done with status file **********/

    fprintf(stderr, "\nRead trajectory. Printing parameters to file\n");
    gmx_rmpbc_done(gpbc);

    /* average over frames */
    for (i = 1; i < ngrps - 1; i++)
    {
        svmul(1.0/nr_frames, (*order)[i], (*order)[i]);
        fprintf(stderr, "Atom %d Tensor: x=%g , y=%g, z=%g\n", i, (*order)[i][XX],
                (*order)[i][YY], (*order)[i][ZZ]);
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

    if (bUnsat)
    {
        fprintf(stderr, "Average angle between double bond and normal: %f\n",
                180*sdbangle/(nr_frames * size*M_PI));
    }

    sfree(x0); /* free memory used by coordinate arrays */
    sfree(x1);
    if (comidx != NULL)
    {
        sfree(comidx);
    }
    if (distidx != NULL)
    {
        sfree(distidx);
    }
    if (grpname != NULL)
    {
        sfree(grpname);
    }
}


void order_plot(rvec order[], real *slOrder[], const char *afile, const char *bfile,
                const char *cfile, int ngrps, int nslices, real slWidth, gmx_bool bSzonly,
                gmx_bool permolecule, real **distvals, const output_env_t oenv)
{
    FILE       *ord, *slOrd;      /* xvgr files with order parameters  */
    int         atom, slice;      /* atom corresponding to order para.*/
    char        buf[256];         /* for xvgr title */
    real        S;                /* order parameter averaged over all atoms */

    if (permolecule)
    {
        sprintf(buf, "Scd order parameters");
        ord = xvgropen(afile, buf, "Atom", "S", oenv);
        sprintf(buf, "Orderparameters per atom per slice");
        slOrd = xvgropen(bfile, buf, "Molecule", "S", oenv);
        for (atom = 1; atom < ngrps - 1; atom++)
        {
            fprintf(ord, "%12d   %12g\n", atom, -1 * (0.6667 * order[atom][XX] +
                                                      0.333 * order[atom][YY]));
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
            fprintf(slOrd, "%12g     %12g\n", slice*slWidth, S/atom);
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
            fprintf(ord, "%12d   %12g   %12g   %12g\n", atom, order[atom][XX],
                    order[atom][YY], order[atom][ZZ]);
            fprintf(slOrd, "%12d   %12g\n", atom, -1 * (0.6667 * order[atom][XX] +
                                                        0.333 * order[atom][YY]));
        }

    }
    xvgrclose(ord);
    xvgrclose(slOrd);
}

void write_bfactors(t_filenm  *fnm, int nfile, atom_id *index, atom_id *a, int nslices, int ngrps, real **order, t_topology *top, real **distvals, output_env_t oenv)
{
    /*function to write order parameters as B factors in PDB file using
          first frame of trajectory*/
    t_trxstatus *status;
    int          natoms;
    t_trxframe   fr, frout;
    t_atoms      useatoms;
    int          i, j, ctr, nout;

    ngrps -= 2;  /*we don't have an order parameter for the first or
                       last atom in each chain*/
    nout   = nslices*ngrps;
    natoms = read_first_frame(oenv, &status, ftp2fn(efTRX, nfile, fnm), &fr,
                              TRX_NEED_X);
    close_trj(status);
    frout        = fr;
    frout.natoms = nout;
    frout.bF     = FALSE;
    frout.bV     = FALSE;
    frout.x      = 0;
    snew(frout.x, nout);

    init_t_atoms(&useatoms, nout, TRUE);
    useatoms.nr = nout;

    /*initialize PDBinfo*/
    for (i = 0; i < useatoms.nr; ++i)
    {
        useatoms.pdbinfo[i].type         = 0;
        useatoms.pdbinfo[i].occup        = 0.0;
        useatoms.pdbinfo[i].bfac         = 0.0;
        useatoms.pdbinfo[i].bAnisotropic = FALSE;
    }

    for (j = 0, ctr = 0; j < nslices; j++)
    {
        for (i = 0; i < ngrps; i++, ctr++)
        {
            /*iterate along each chain*/
            useatoms.pdbinfo[ctr].bfac = order[j][i+1];
            if (distvals)
            {
                useatoms.pdbinfo[ctr].occup = distvals[j][i+1];
            }
            copy_rvec(fr.x[a[index[i+1]+j]], frout.x[ctr]);
            useatoms.atomname[ctr] = top->atoms.atomname[a[index[i+1]+j]];
            useatoms.atom[ctr]     = top->atoms.atom[a[index[i+1]+j]];
            useatoms.nres          = max(useatoms.nres, useatoms.atom[ctr].resind+1);
            useatoms.resinfo[useatoms.atom[ctr].resind] = top->atoms.resinfo[useatoms.atom[ctr].resind]; /*copy resinfo*/
        }
    }

    write_sto_conf(opt2fn("-ob", nfile, fnm), "Order parameters", &useatoms, frout.x, NULL, frout.ePBC, frout.box);

    sfree(frout.x);
    free_t_atoms(&useatoms, FALSE);
}

int gmx_order(int argc, char *argv[])
{
    const char        *desc[] = {
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
        "given.[PAR]"
        "The tetrahedrality order parameters can be determined",
        "around an atom. Both angle an distance order parameters are calculated. See",
        "P.-L. Chau and A.J. Hardwick, Mol. Phys., 93, (1998), 511-518.",
        "for more details."
    };

    static int         nslices       = 1;     /* nr of slices defined       */
    static gmx_bool    bSzonly       = FALSE; /* True if only Sz is wanted  */
    static gmx_bool    bUnsat        = FALSE; /* True if carbons are unsat. */
    static const char *normal_axis[] = { NULL, "z", "x", "y", NULL };
    static gmx_bool    permolecule   = FALSE; /*compute on a per-molecule basis */
    static gmx_bool    radial        = FALSE; /*compute a radial membrane normal */
    static gmx_bool    distcalc      = FALSE; /*calculate distance from a reference group */
    t_pargs            pa[]          = {
        { "-d",      FALSE, etENUM, {normal_axis},
          "Direction of the normal on the membrane" },
        { "-sl",     FALSE, etINT, {&nslices},
          "Calculate order parameter as function of box length, dividing the box"
          " into this number of slices." },
        { "-szonly", FALSE, etBOOL, {&bSzonly},
          "Only give Sz element of order tensor. (axis can be specified with [TT]-d[tt])" },
        { "-unsat",  FALSE, etBOOL, {&bUnsat},
          "Calculate order parameters for unsaturated carbons. Note that this can"
          "not be mixed with normal order parameters." },
        { "-permolecule", FALSE, etBOOL, {&permolecule},
          "Compute per-molecule Scd order parameters" },
        { "-radial", FALSE, etBOOL, {&radial},
          "Compute a radial membrane normal" },
        { "-calcdist", FALSE, etBOOL, {&distcalc},
          "Compute distance from a reference" },
    };

    rvec              *order;                         /* order par. for each atom   */
    real             **slOrder;                       /* same, per slice            */
    real               slWidth = 0.0;                 /* width of a slice           */
    char             **grpname;                       /* groupnames                 */
    int                ngrps,                         /* nr. of groups              */
                       i,
                       axis = 0;                      /* normal axis                */
    t_topology   *top;                                /* topology         */
    int           ePBC;
    atom_id      *index,                              /* indices for a              */
    *a;                                               /* atom numbers in each group */
    t_blocka     *block;                              /* data from index file       */
    t_filenm      fnm[] = {                           /* files for g_order    */
        { efTRX, "-f", NULL,  ffREAD },               /* trajectory file              */
        { efNDX, "-n", NULL,  ffREAD },               /* index file           */
        { efNDX, "-nr", NULL,  ffREAD },              /* index for radial axis calculation	  */
        { efTPR, NULL, NULL,  ffREAD },               /* topology file                */
        { efXVG, "-o", "order", ffWRITE },            /* xvgr output file     */
        { efXVG, "-od", "deuter", ffWRITE },          /* xvgr output file           */
        { efPDB, "-ob", NULL, ffWRITE },              /* write Scd as B factors to PDB if permolecule           */
        { efXVG, "-os", "sliced", ffWRITE },          /* xvgr output file           */
        { efXVG, "-Sg", "sg-ang", ffOPTWR },          /* xvgr output file           */
        { efXVG, "-Sk", "sk-dist", ffOPTWR },         /* xvgr output file           */
        { efXVG, "-Sgsl", "sg-ang-slice", ffOPTWR },  /* xvgr output file           */
        { efXVG, "-Sksl", "sk-dist-slice", ffOPTWR }, /* xvgr output file           */
    };
    gmx_bool      bSliced = FALSE;                    /* True if box is sliced      */
#define NFILE asize(fnm)
    real        **distvals = NULL;
    const char   *sgfnm, *skfnm, *ndxfnm, *tpsfnm, *trxfnm;
    output_env_t  oenv;

    if (!parse_common_args(&argc, argv, PCA_CAN_VIEW | PCA_CAN_TIME,
                           NFILE, fnm, asize(pa), pa, asize(desc), desc, 0, NULL, &oenv))
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
    if (strcmp(normal_axis[0], "x") == 0)
    {
        axis = XX;
    }
    else if (strcmp(normal_axis[0], "y") == 0)
    {
        axis = YY;
    }
    else if (strcmp(normal_axis[0], "z") == 0)
    {
        axis = ZZ;
    }
    else
    {
        gmx_fatal(FARGS, "Invalid axis, use x, y or z");
    }

    switch (axis)
    {
        case 0:
            fprintf(stderr, "Taking x axis as normal to the membrane\n");
            break;
        case 1:
            fprintf(stderr, "Taking y axis as normal to the membrane\n");
            break;
        case 2:
            fprintf(stderr, "Taking z axis as normal to the membrane\n");
            break;
    }

    /* tetraheder order parameter */
    if (skfnm || sgfnm)
    {
        /* If either of theoptions is set we compute both */
        sgfnm = opt2fn("-Sg", NFILE, fnm);
        skfnm = opt2fn("-Sk", NFILE, fnm);
        calc_tetra_order_parm(ndxfnm, tpsfnm, trxfnm, sgfnm, skfnm, nslices, axis,
                              opt2fn("-Sgsl", NFILE, fnm), opt2fn("-Sksl", NFILE, fnm),
                              oenv);
        /* view xvgr files */
        do_view(oenv, opt2fn("-Sg", NFILE, fnm), NULL);
        do_view(oenv, opt2fn("-Sk", NFILE, fnm), NULL);
        if (nslices > 1)
        {
            do_view(oenv, opt2fn("-Sgsl", NFILE, fnm), NULL);
            do_view(oenv, opt2fn("-Sksl", NFILE, fnm), NULL);
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
        if (bUnsat)
        {
            fprintf(stderr, "Taking carbons as unsaturated!\n");
        }

        top = read_top(ftp2fn(efTPR, NFILE, fnm), &ePBC); /* read topology file */

        block = init_index(ftp2fn(efNDX, NFILE, fnm), &grpname);
        index = block->index;                   /* get indices from t_block block */
        a     = block->a;                       /* see block.h                    */
        ngrps = block->nr;

        if (permolecule)
        {
            nslices = index[1] - index[0]; /*I think this assumes contiguous lipids in topology*/
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
        print_types(index, a, ngrps, grpname, top);

        calc_order(ftp2fn(efTRX, NFILE, fnm), index, a, &order,
                   &slOrder, &slWidth, nslices, bSliced, bUnsat,
                   top, ePBC, ngrps, axis, permolecule, radial, distcalc, opt2fn_null("-nr", NFILE, fnm), &distvals, oenv);

        if (radial)
        {
            ngrps--; /*don't print the last group--was used for
                               center-of-mass determination*/

        }
        order_plot(order, slOrder, opt2fn("-o", NFILE, fnm), opt2fn("-os", NFILE, fnm),
                   opt2fn("-od", NFILE, fnm), ngrps, nslices, slWidth, bSzonly, permolecule, distvals, oenv);

        if (opt2bSet("-ob", NFILE, fnm))
        {
            if (!permolecule)
            {
                fprintf(stderr,
                        "Won't write B-factors with averaged order parameters; use -permolecule\n");
            }
            else
            {
                write_bfactors(fnm, NFILE, index, a, nslices, ngrps, slOrder, top, distvals, oenv);
            }
        }


        do_view(oenv, opt2fn("-o", NFILE, fnm), NULL);  /* view xvgr file */
        do_view(oenv, opt2fn("-os", NFILE, fnm), NULL); /* view xvgr file */
        do_view(oenv, opt2fn("-od", NFILE, fnm), NULL); /* view xvgr file */
    }

    if (distvals != NULL)
    {
        for (i = 0; i < nslices; ++i)
        {
            sfree(distvals[i]);
        }
        sfree(distvals);
    }

    return 0;
}
