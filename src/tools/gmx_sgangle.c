/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <math.h>

#include "sysstuff.h"
#include <string.h>
#include "typedefs.h"
#include "smalloc.h"
#include "macros.h"
#include "rmpbc.h"
#include "vec.h"
#include "xvgr.h"
#include "pbc.h"
#include "copyrite.h"
#include "futil.h"
#include "statutil.h"
#include "index.h"
#include "gmx_ana.h"


/* this version only works correctly if one of the entries in the index file
   is a plane (three atoms specified) and the other a vector. Distance
   is calculated from the center of the plane to both atoms of the vector */

static void print_types(atom_id index1[], int gnx1, char *group1,
                        atom_id index2[], int gnx2, char *group2,
                        t_topology *top)
{
    int i, j;

    fprintf(stderr, "\n");
    fprintf(stderr, "Group %s contains the following atoms: \n", group1);
    for (i = 0; i < gnx1; i++)
    {
        fprintf(stderr, "Atomname %d: %s\n", i, *(top->atoms.atomname[index1[i]]));
    }
    fprintf(stderr, "\n");

    fprintf(stderr, "Group %s contains the following atoms: \n", group2);
    for (j = 0; j < gnx2; j++)
    {
        fprintf(stderr, "Atomname %d: %s\n", j, *(top->atoms.atomname[index2[j]]));
    }
    fprintf(stderr, "\n");

    fprintf(stderr, "Careful: distance only makes sense in some situations.\n\n");
}

static void calculate_normal(atom_id index[], rvec x[], rvec result, rvec center)
{
    rvec c1, c2;
    int  i;

    /* calculate centroid of triangle spanned by the three points */
    for (i = 0; i < 3; i++)
    {
        center[i] = (x[index[0]][i] + x[index[1]][i] + x[index[2]][i])/3;
    }

    /* use P1P2 x P1P3 to calculate normal, given three points P1-P3 */
    rvec_sub(x[index[1]], x[index[0]], c1); /* find two vectors */
    rvec_sub(x[index[2]], x[index[0]], c2);

    cprod(c1, c2, result);                /* take crossproduct between these */
}

/* calculate the angle and distance between the two groups */
static void calc_angle(int ePBC, matrix box, rvec x[], atom_id index1[],
                       atom_id index2[], int gnx1, int gnx2,
                       real *angle,      real *distance,
                       real *distance1,  real *distance2)

/* distance is distance between centers, distance 1 between center of plane
   and atom one of vector, distance 2 same for atom two
 */

{
    rvec
        normal1, normal2,   /* normals on planes of interest */
        center1, center2,   /* center of triangle of points given to define plane,*/
                            /* or center of vector if a vector is given */
        h1, h2, h3, h4, h5; /* temp. vectors */
    t_pbc pbc;

    set_pbc(&pbc, ePBC, box);

    switch (gnx1)
    {
        case 3:       /* group 1 defines plane */
            calculate_normal(index1, x, normal1, center1);
            break;
        case 2:       /* group 1 defines vector */
            rvec_sub(x[index1[0]], x[index1[1]], normal1);
            rvec_add(x[index1[0]], x[index1[1]], h1);
            svmul(0.5, h1, center1); /* center is geometric mean */
            break;
        default:                     /* group 1 does none of the above */
            gmx_fatal(FARGS, "Something wrong with contents of index file. Groups should contain 2 or 3 atoms.\n");
    }

    switch (gnx2)
    {
        case 3:      /* group 2 defines plane */
            calculate_normal(index2, x, normal2, center2);
            break;
        case 2:      /* group 2 defines vector */
            rvec_sub(x[index2[0]], x[index2[1]], normal2);
            rvec_add(x[index2[0]], x[index2[1]], h2);
            svmul(0.5, h2, center2); /* center is geometric mean */
            break;
        case 0:
            normal2[XX] = 0;
            normal2[YY] = 0;
            normal2[ZZ] = 1;
            center2[XX] = 0;
            center2[YY] = 0;
            center2[ZZ] = 0;
            break;
        default:     /* group 2 does none of the above */
            gmx_fatal(FARGS, "Something wrong with contents of index file.\n");
    }

    *angle = cos_angle(normal1, normal2);

    if (box)
    {
        pbc_dx(&pbc, center1, center2, h3);
    }
    else
    {
        rvec_sub(center1, center2, h3);
    }
    *distance = norm(h3);

    if (gnx1 == 3 && gnx2 == 2)
    {
        if (box)
        {
            pbc_dx(&pbc, center1, x[index2[0]], h4);
            pbc_dx(&pbc, center1, x[index2[1]], h5);
        }
        else
        {
            rvec_sub(center1, x[index2[0]], h4);
            rvec_sub(center1, x[index2[1]], h5);
        }
        *distance1 = norm(h4);
        *distance2 = norm(h5);
    }
    else if (gnx1 == 2 && gnx2 == 3)
    {
        rvec_sub(center1, x[index1[0]], h4);
        rvec_sub(center1, x[index1[1]], h5);
        *distance1 = norm(h4);
        *distance2 = norm(h5);
    }
    else
    {
        *distance1 = 0; *distance2 = 0;
    }
}

void sgangle_plot(const char *fn, const char *afile, const char *dfile,
                  const char *d1file, const char *d2file,
                  atom_id index1[], int gnx1, char *grpn1,
                  atom_id index2[], int gnx2, char *grpn2,
                  t_topology *top, int ePBC, const output_env_t oenv)
{
    FILE
    *sg_angle,            /* xvgr file with angles */
    *sg_distance  = NULL, /* xvgr file with distances */
    *sg_distance1 = NULL, /* xvgr file with distance between plane and atom */
    *sg_distance2 = NULL; /* xvgr file with distance between plane and atom2 */
    real
        t,                /* time */
        angle,            /* cosine of angle between two groups */
        distance,         /* distance between two groups. */
        distance1,        /* distance between plane and one of two atoms */
        distance2;        /* same for second of two atoms */
    t_trxstatus *status;
    int          natoms, teller = 0;
    rvec        *x0;       /* coordinates, and coordinates corrected for pb */
    matrix       box;
    char         buf[256]; /* for xvgr title */
    gmx_rmpbc_t  gpbc = NULL;
    const char*  aleg[2] = { "cos(Angle)", "Angle (degrees)" };     /* legends for sg_angle output file */

    if ((natoms = read_first_x(oenv, &status, fn, &t, &x0, box)) == 0)
    {
        gmx_fatal(FARGS, "Could not read coordinates from statusfile\n");
    }

    sprintf(buf, "Angle between %s and %s", grpn1, grpn2);
    sg_angle = xvgropen(afile, buf, "Time (ps)", "Angle (degrees)", oenv);
    xvgr_legend(sg_angle, 2, aleg, oenv);

    if (dfile)
    {
        sprintf(buf, "Distance between %s and %s", grpn1, grpn2);
        sg_distance = xvgropen(dfile, buf, "Time (ps)", "Distance (nm)", oenv);
    }

    if (d1file)
    {
        sprintf(buf, "Distance between plane and first atom of vector");
        sg_distance1 = xvgropen(d1file, buf, "Time (ps)", "Distance (nm)", oenv);
    }

    if (d2file)
    {
        sprintf(buf, "Distance between plane and second atom of vector");
        sg_distance2 = xvgropen(d2file, buf, "Time (ps", "Distance (nm)", oenv);
    }

    gpbc = gmx_rmpbc_init(&(top->idef), ePBC, natoms, box);

    do
    {
        teller++;

        gmx_rmpbc(gpbc, natoms, box, x0);

        calc_angle(ePBC, box, x0, index1, index2, gnx1, gnx2, &angle,
                   &distance, &distance1, &distance2);

        fprintf(sg_angle, "%12g  %12g  %12g\n", t, angle, acos(angle)*180.0/M_PI);
        if (dfile)
        {
            fprintf(sg_distance, "%12g  %12g\n", t, distance);
        }
        if (d1file)
        {
            fprintf(sg_distance1, "%12g  %12g\n", t, distance1);
        }
        if (d2file)
        {
            fprintf(sg_distance2, "%12g  %12g\n", t, distance1);
        }

    }
    while (read_next_x(oenv, status, &t, natoms, x0, box));

    gmx_rmpbc_done(gpbc);

    fprintf(stderr, "\n");
    close_trj(status);
    ffclose(sg_angle);
    if (dfile)
    {
        ffclose(sg_distance);
    }
    if (d1file)
    {
        ffclose(sg_distance1);
    }
    if (d2file)
    {
        ffclose(sg_distance2);
    }

    sfree(x0);
}

static void calc_angle_single(int     ePBC,
                              matrix  box,
                              rvec    xzero[],
                              rvec    x[],
                              atom_id index1[],
                              atom_id index2[],
                              int     gnx1,
                              int     gnx2,
                              real   *angle,
                              real   *distance,
                              real   *distance1,
                              real   *distance2)
{
    t_pbc pbc;

    /* distance is distance between centers, distance 1 between center of plane
       and atom one of vector, distance 2 same for atom two
     */

    rvec  normal1, normal2; /* normals on planes of interest */
    rvec  center1, center2;
    /* center of triangle of pts to define plane,
     * or center of vector if a vector is given
     */
    rvec  h1, h2, h3, h4, h5; /* temp. vectors */

    if (box)
    {
        set_pbc(&pbc, ePBC, box);
    }

    switch (gnx1)
    {
        case 3:     /* group 1 defines plane */
            calculate_normal(index1, xzero, normal1, center1);
            break;
        case 2:     /* group 1 defines vector */
            rvec_sub(xzero[index1[0]], xzero[index1[1]], normal1);
            rvec_add(xzero[index1[0]], xzero[index1[1]], h1);
            svmul(0.5, h1, center1); /* center is geometric mean */
            break;
        default:                     /* group 1 does none of the above */
            gmx_fatal(FARGS, "Something wrong with contents of index file.\n");
    }

    switch (gnx2)
    {
        case 3:    /* group 2 defines plane */
            calculate_normal(index2, x, normal2, center2);
            break;
        case 2:    /* group 2 defines vector */
            rvec_sub(x[index2[0]], x[index2[1]], normal2);
            rvec_add(x[index2[0]], x[index2[1]], h2);
            svmul(0.5, h2, center2); /* center is geometric mean */
            break;
        default:                     /* group 2 does none of the above */
            gmx_fatal(FARGS, "Something wrong with contents of index file.\n");
    }

    *angle = cos_angle(normal1, normal2);

    if (box)
    {
        pbc_dx(&pbc, center1, center2, h3);
    }
    else
    {
        rvec_sub(center1, center2, h3);
    }
    *distance = norm(h3);

    if (gnx1 == 3 && gnx2 == 2)
    {
        if (box)
        {
            pbc_dx(&pbc, center1, x[index2[0]], h4);
            pbc_dx(&pbc, center1, x[index2[1]], h5);
        }
        else
        {
            rvec_sub(center1, x[index2[0]], h4);
            rvec_sub(center1, x[index2[1]], h5);
        }
        *distance1 = norm(h4);
        *distance2 = norm(h5);
    }
    else if (gnx1 == 2 && gnx2 == 3)
    {
        rvec_sub(center1, xzero[index1[0]], h4);
        rvec_sub(center1, xzero[index1[1]], h5);
        *distance1 = norm(h4);
        *distance2 = norm(h5);
    }
    else
    {
        *distance1 = 0; *distance2 = 0;
    }
}


void sgangle_plot_single(const char *fn, const char *afile, const char *dfile,
                         const char *d1file, const char *d2file,
                         atom_id index1[], int gnx1, char *grpn1,
                         atom_id index2[], int gnx2, char *grpn2,
                         t_topology *top, int ePBC, const output_env_t oenv)
{
    FILE
    *sg_angle,            /* xvgr file with angles */
    *sg_distance  = NULL, /* xvgr file with distances */
    *sg_distance1 = NULL, /* xvgr file with distance between plane and atom */
    *sg_distance2 = NULL; /* xvgr file with distance between plane and atom2 */
    real
        t,                /* time */
        angle,            /* cosine of angle between two groups */
        distance,         /* distance between two groups. */
        distance1,        /* distance between plane and one of two atoms */
        distance2;        /* same for second of two atoms */
    t_trxstatus *status;
    int          natoms, teller = 0;
    int          i;
    rvec        *x0; /* coordinates, and coordinates corrected for pb */
    rvec        *xzero;
    matrix       box;
    char         buf[256]; /* for xvgr title */
    gmx_rmpbc_t  gpbc = NULL;


    if ((natoms = read_first_x(oenv, &status, fn, &t, &x0, box)) == 0)
    {
        gmx_fatal(FARGS, "Could not read coordinates from statusfile\n");
    }

    sprintf(buf, "Angle between %s and %s", grpn1, grpn2);
    sg_angle = xvgropen(afile, buf, "Time (ps)", "Cos(angle) ", oenv);

    if (dfile)
    {
        sprintf(buf, "Distance between %s and %s", grpn1, grpn2);
        sg_distance = xvgropen(dfile, buf, "Time (ps)", "Distance (nm)", oenv);
    }

    if (d1file)
    {
        sprintf(buf, "Distance between plane and first atom of vector");
        sg_distance1 = xvgropen(d1file, buf, "Time (ps)", "Distance (nm)", oenv);
    }

    if (d2file)
    {
        sprintf(buf, "Distance between plane and second atom of vector");
        sg_distance2 = xvgropen(d2file, buf, "Time (ps", "Distance (nm)", oenv);
    }

    snew(xzero, natoms);
    gpbc = gmx_rmpbc_init(&top->idef, ePBC, natoms, box);

    do
    {
        teller++;

        gmx_rmpbc(gpbc, natoms, box, x0);
        if (teller == 1)
        {
            for (i = 0; i < natoms; i++)
            {
                copy_rvec(x0[i], xzero[i]);
            }
        }


        calc_angle_single(ePBC, box,
                          xzero, x0, index1, index2, gnx1, gnx2, &angle,
                          &distance, &distance1, &distance2);

        fprintf(sg_angle, "%12g  %12g  %12g\n", t, angle, acos(angle)*180.0/M_PI);
        if (dfile)
        {
            fprintf(sg_distance, "%12g  %12g\n", t, distance);
        }
        if (d1file)
        {
            fprintf(sg_distance1, "%12g  %12g\n", t, distance1);
        }
        if (d2file)
        {
            fprintf(sg_distance2, "%12g  %12g\n", t, distance1);
        }

    }
    while (read_next_x(oenv, status, &t, natoms, x0, box));
    gmx_rmpbc_done(gpbc);

    fprintf(stderr, "\n");
    close_trj(status);
    ffclose(sg_angle);
    if (dfile)
    {
        ffclose(sg_distance);
    }
    if (d1file)
    {
        ffclose(sg_distance1);
    }
    if (d2file)
    {
        ffclose(sg_distance2);
    }

    sfree(x0);
}



int gmx_sgangle(int argc, char *argv[])
{
    const char     *desc[] = {
        "Compute the angle and distance between two groups. ",
        "The groups are defined by a number of atoms given in an index file and",
        "may be two or three atoms in size.",
        "If [TT]-one[tt] is set, only one group should be specified in the index",
        "file and the angle between this group at time 0 and t will be computed.",
        "The angles calculated depend on the order in which the atoms are ",
        "given. Giving, for instance, 5 6 will rotate the vector 5-6 with ",
        "180 degrees compared to giving 6 5. [PAR]If three atoms are given, ",
        "the normal on the plane spanned by those three atoms will be",
        "calculated, using the formula  P1P2 x P1P3.",
        "The cos of the angle is calculated, using the inproduct of the two",
        "normalized vectors.[PAR]",
        "Here is what some of the file options do:[BR]",
        "[TT]-oa[tt]: Angle between the two groups specified in the index file. If a group contains three atoms the normal to the plane defined by those three atoms will be used. If a group contains two atoms, the vector defined by those two atoms will be used.[BR]",
        "[TT]-od[tt]: Distance between two groups. Distance is taken from the center of one group to the center of the other group.[BR]",
        "[TT]-od1[tt]: If one plane and one vector is given, the distances for each of the atoms from the center of the plane is given separately.[BR]",
        "[TT]-od2[tt]: For two planes this option has no meaning."
    };

    output_env_t    oenv;
    const char     *fna, *fnd, *fnd1, *fnd2;
    char     *      grpname[2];         /* name of the two groups */
    int             gnx[2];             /* size of the two groups */
    t_topology     *top;                /* topology         */
    int             ePBC;
    atom_id        *index[2];
    static gmx_bool bOne = FALSE, bZ = FALSE;
    t_pargs         pa[] = {
        { "-one", FALSE, etBOOL, {&bOne},
          "Only one group compute angle between vector at time zero and time t"},
        { "-z", FALSE, etBOOL, {&bZ},
          "Use the [IT]z[it]-axis as reference" }
    };
#define NPA asize(pa)

    t_filenm  fnm[] = {                         /* files for g_sgangle  */
        { efTRX, "-f", NULL,  ffREAD },         /* trajectory file  */
        { efNDX, NULL, NULL,  ffREAD },         /* index file       */
        { efTPX, NULL, NULL,  ffREAD },         /* topology file    */
        { efXVG, "-oa", "sg_angle", ffWRITE },  /* xvgr output file     */
        { efXVG, "-od", "sg_dist", ffOPTWR },   /* xvgr output file     */
        { efXVG, "-od1", "sg_dist1", ffOPTWR }, /* xvgr output file     */
        { efXVG, "-od2", "sg_dist2", ffOPTWR }  /* xvgr output file     */
    };

#define NFILE asize(fnm)

    CopyRight(stderr, argv[0]);
    parse_common_args(&argc, argv, PCA_CAN_VIEW | PCA_CAN_TIME | PCA_BE_NICE,
                      NFILE, fnm, NPA, pa, asize(desc), desc, 0, NULL, &oenv);


    top = read_top(ftp2fn(efTPX, NFILE, fnm), &ePBC); /* read topology file */

    fna  = opt2fn("-oa", NFILE, fnm);
    fnd  = opt2fn_null("-od", NFILE, fnm);
    fnd1 = opt2fn_null("-od1", NFILE, fnm);
    fnd2 = opt2fn_null("-od2", NFILE, fnm);

    /* read index file. */
    if (bOne)
    {
        rd_index(ftp2fn(efNDX, NFILE, fnm), 1, gnx, index, grpname);
        print_types(index[0], gnx[0], grpname[0],
                    index[0], gnx[0], grpname[0], top);

        sgangle_plot_single(ftp2fn(efTRX, NFILE, fnm), fna, fnd, fnd1, fnd2,
                            index[0], gnx[0], grpname[0],
                            index[0], gnx[0], grpname[0],
                            top, ePBC, oenv);
    }
    else
    {
        rd_index(ftp2fn(efNDX, NFILE, fnm), bZ ? 1 : 2, gnx, index, grpname);
        if (!bZ)
        {
            print_types(index[0], gnx[0], grpname[0],
                        index[1], gnx[1], grpname[1], top);
        }
        else
        {
            gnx[1]     = 0;
            grpname[1] = strdup("Z-axis");
        }
        sgangle_plot(ftp2fn(efTRX, NFILE, fnm), fna, fnd, fnd1, fnd2,
                     index[0], gnx[0], grpname[0],
                     index[1], gnx[1], grpname[1],
                     top, ePBC, oenv);
    }

    do_view(oenv, fna, "-nxy");  /* view xvgr file */
    do_view(oenv, fnd, "-nxy");  /* view xvgr file */
    do_view(oenv, fnd1, "-nxy"); /* view xvgr file */
    do_view(oenv, fnd2, "-nxy"); /* view xvgr file */

    thanx(stderr);
    return 0;
}
