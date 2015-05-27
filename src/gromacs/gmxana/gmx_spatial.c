/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2007,2008,2009,2010,2011,2012,2013,2014,2015, by the GROMACS development team, led by
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
#include <stdlib.h>

#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pbcutil/rmpbc.h"
#include "gromacs/topology/index.h"
#include "gromacs/utility/smalloc.h"

static const double bohr = 0.529177249;  /* conversion factor to compensate for VMD plugin conversion... */

int gmx_spatial(int argc, char *argv[])
{
    const char     *desc[] = {
        "[THISMODULE] calculates the spatial distribution function and",
        "outputs it in a form that can be read by VMD as Gaussian98 cube format.",
        "For a system of 32,000 atoms and a 50 ns trajectory, the SDF can be generated",
        "in about 30 minutes, with most of the time dedicated to the two runs through",
        "[TT]trjconv[tt] that are required to center everything properly.",
        "This also takes a whole bunch of space (3 copies of the trajectory file).",
        "Still, the pictures are pretty and very informative when the fitted selection is properly made.",
        "3-4 atoms in a widely mobile group (like a free amino acid in solution) works",
        "well, or select the protein backbone in a stable folded structure to get the SDF",
        "of solvent and look at the time-averaged solvation shell.",
        "It is also possible using this program to generate the SDF based on some arbitrary",
        "Cartesian coordinate. To do that, simply omit the preliminary [gmx-trjconv] steps.",
        "",
        "Usage:",
        "",
        "1. Use [gmx-make_ndx] to create a group containing the atoms around which you want the SDF",
        "2. [TT]gmx trjconv -s a.tpr -f a.tng -o b.tng -boxcenter tric -ur compact -pbc none[tt]",
        "3. [TT]gmx trjconv -s a.tpr -f b.tng -o c.tng -fit rot+trans[tt]",
        "4. run [THISMODULE] on the [TT]c.tng[tt] output of step #3.",
        "5. Load [TT]grid.cube[tt] into VMD and view as an isosurface.",
        "",
        "[BB]Note[bb] that systems such as micelles will require [TT]gmx trjconv -pbc cluster[tt] between steps 1 and 2.",
        "",
        "Warnings",
        "^^^^^^^^",
        "",
        "The SDF will be generated for a cube that contains all bins that have some non-zero occupancy.",
        "However, the preparatory [TT]-fit rot+trans[tt] option to [gmx-trjconv] implies that your system will be rotating",
        "and translating in space (in order that the selected group does not). Therefore the values that are",
        "returned will only be valid for some region around your central group/coordinate that has full overlap",
        "with system volume throughout the entire translated/rotated system over the course of the trajectory.",
        "It is up to the user to ensure that this is the case.",
        "",
        "Risky options",
        "^^^^^^^^^^^^^",
        "",
        "To reduce the amount of space and time required, you can output only the coords",
        "that are going to be used in the first and subsequent run through [gmx-trjconv].",
        "However, be sure to set the [TT]-nab[tt] option to a sufficiently high value since",
        "memory is allocated for cube bins based on the initial coordinates and the [TT]-nab[tt]",
        "option value."
    };
    const char     *bugs[] = {
        "When the allocated memory is not large enough, a segmentation fault may occur. This is usually detected "
        "and the program is halted prior to the fault while displaying a warning message suggesting the use of the [TT]-nab[tt] (Number of Additional Bins) "
        "option. However, the program does not detect all such events. If you encounter a segmentation fault, run it again "
        "with an increased [TT]-nab[tt] value."
    };

    static gmx_bool bPBC         = FALSE;
    static gmx_bool bSHIFT       = FALSE;
    static int      iIGNOREOUTER = -1;   /*Positive values may help if the surface is spikey */
    static gmx_bool bCUTDOWN     = TRUE;
    static real     rBINWIDTH    = 0.05; /* nm */
    static gmx_bool bCALCDIV     = TRUE;
    static int      iNAB         = 4;

    t_pargs         pa[] = {
        { "-pbc",      FALSE, etBOOL, {&bPBC},
          "Use periodic boundary conditions for computing distances" },
        { "-div",      FALSE, etBOOL, {&bCALCDIV},
          "Calculate and apply the divisor for bin occupancies based on atoms/minimal cube size. Set as TRUE for visualization and as FALSE ([TT]-nodiv[tt]) to get accurate counts per frame" },
        { "-ign",      FALSE, etINT, {&iIGNOREOUTER},
          "Do not display this number of outer cubes (positive values may reduce boundary speckles; -1 ensures outer surface is visible)" },
        /*    { "-cut",      bCUTDOWN, etBOOL, {&bCUTDOWN},*/
        /*      "Display a total cube that is of minimal size" }, */
        { "-bin",      FALSE, etREAL, {&rBINWIDTH},
          "Width of the bins (nm)" },
        { "-nab",      FALSE, etINT, {&iNAB},
          "Number of additional bins to ensure proper memory allocation" }
    };

    double          MINBIN[3];
    double          MAXBIN[3];
    t_topology      top;
    int             ePBC;
    char            title[STRLEN];
    t_trxframe      fr;
    rvec           *xtop, *shx[26];
    matrix          box, box_pbc;
    t_trxstatus    *status;
    int             flags = TRX_READ_X;
    t_pbc           pbc;
    t_atoms        *atoms;
    int             natoms;
    char           *grpnm, *grpnmp;
    atom_id        *index, *indexp;
    int             i, nidx, nidxp;
    int             v;
    int             j, k;
    long         ***bin = (long ***)NULL;
    long            nbin[3];
    FILE           *flp;
    long            x, y, z, minx, miny, minz, maxx, maxy, maxz;
    long            numfr, numcu;
    long            tot, max, min;
    double          norm;
    output_env_t    oenv;
    gmx_rmpbc_t     gpbc = NULL;

    t_filenm        fnm[] = {
        { efTPS,  NULL,  NULL, ffREAD }, /* this is for the topology */
        { efTRX, "-f", NULL, ffREAD },   /* and this for the trajectory */
        { efNDX, NULL, NULL, ffOPTRD }
    };

#define NFILE asize(fnm)

    /* This is the routine responsible for adding default options,
     * calling the X/motif interface, etc. */
    if (!parse_common_args(&argc, argv, PCA_CAN_TIME | PCA_CAN_VIEW,
                           NFILE, fnm, asize(pa), pa, asize(desc), desc, asize(bugs), bugs, &oenv))
    {
        return 0;
    }

    read_tps_conf(ftp2fn(efTPS, NFILE, fnm), title, &top, &ePBC, &xtop, NULL, box, TRUE);
    sfree(xtop);

    atoms = &(top.atoms);
    printf("Select group to generate SDF:\n");
    get_index(atoms, ftp2fn_null(efNDX, NFILE, fnm), 1, &nidx, &index, &grpnm);
    printf("Select group to output coords (e.g. solute):\n");
    get_index(atoms, ftp2fn_null(efNDX, NFILE, fnm), 1, &nidxp, &indexp, &grpnmp);

    /* The first time we read data is a little special */
    natoms = read_first_frame(oenv, &status, ftp2fn(efTRX, NFILE, fnm), &fr, flags);

    /* Memory Allocation */
    MINBIN[XX] = MAXBIN[XX] = fr.x[0][XX];
    MINBIN[YY] = MAXBIN[YY] = fr.x[0][YY];
    MINBIN[ZZ] = MAXBIN[ZZ] = fr.x[0][ZZ];
    for (i = 1; i < top.atoms.nr; ++i)
    {
        if (fr.x[i][XX] < MINBIN[XX])
        {
            MINBIN[XX] = fr.x[i][XX];
        }
        if (fr.x[i][XX] > MAXBIN[XX])
        {
            MAXBIN[XX] = fr.x[i][XX];
        }
        if (fr.x[i][YY] < MINBIN[YY])
        {
            MINBIN[YY] = fr.x[i][YY];
        }
        if (fr.x[i][YY] > MAXBIN[YY])
        {
            MAXBIN[YY] = fr.x[i][YY];
        }
        if (fr.x[i][ZZ] < MINBIN[ZZ])
        {
            MINBIN[ZZ] = fr.x[i][ZZ];
        }
        if (fr.x[i][ZZ] > MAXBIN[ZZ])
        {
            MAXBIN[ZZ] = fr.x[i][ZZ];
        }
    }
    for (i = ZZ; i >= XX; --i)
    {
        MAXBIN[i]  = (ceil((MAXBIN[i]-MINBIN[i])/rBINWIDTH)+(double)iNAB)*rBINWIDTH+MINBIN[i];
        MINBIN[i] -= (double)iNAB*rBINWIDTH;
        nbin[i]    = (long)ceil((MAXBIN[i]-MINBIN[i])/rBINWIDTH);
    }
    snew(bin, nbin[XX]);
    for (i = 0; i < nbin[XX]; ++i)
    {
        snew(bin[i], nbin[YY]);
        for (j = 0; j < nbin[YY]; ++j)
        {
            snew(bin[i][j], nbin[ZZ]);
        }
    }
    copy_mat(box, box_pbc);
    numfr = 0;
    minx  = miny = minz = 999;
    maxx  = maxy = maxz = 0;

    if (bPBC)
    {
        gpbc = gmx_rmpbc_init(&top.idef, ePBC, natoms);
    }
    /* This is the main loop over frames */
    do
    {
        /* Must init pbc every step because of pressure coupling */

        copy_mat(box, box_pbc);
        if (bPBC)
        {
            gmx_rmpbc_trxfr(gpbc, &fr);
            set_pbc(&pbc, ePBC, box_pbc);
        }

        for (i = 0; i < nidx; i++)
        {
            if (fr.x[index[i]][XX] < MINBIN[XX] || fr.x[index[i]][XX] > MAXBIN[XX] ||
                fr.x[index[i]][YY] < MINBIN[YY] || fr.x[index[i]][YY] > MAXBIN[YY] ||
                fr.x[index[i]][ZZ] < MINBIN[ZZ] || fr.x[index[i]][ZZ] > MAXBIN[ZZ])
            {
                printf("There was an item outside of the allocated memory. Increase the value given with the -nab option.\n");
                printf("Memory was allocated for [%f,%f,%f]\tto\t[%f,%f,%f]\n", MINBIN[XX], MINBIN[YY], MINBIN[ZZ], MAXBIN[XX], MAXBIN[YY], MAXBIN[ZZ]);
                printf("Memory was required for [%f,%f,%f]\n", fr.x[index[i]][XX], fr.x[index[i]][YY], fr.x[index[i]][ZZ]);
                exit(1);
            }
            x = (long)ceil((fr.x[index[i]][XX]-MINBIN[XX])/rBINWIDTH);
            y = (long)ceil((fr.x[index[i]][YY]-MINBIN[YY])/rBINWIDTH);
            z = (long)ceil((fr.x[index[i]][ZZ]-MINBIN[ZZ])/rBINWIDTH);
            ++bin[x][y][z];
            if (x < minx)
            {
                minx = x;
            }
            if (x > maxx)
            {
                maxx = x;
            }
            if (y < miny)
            {
                miny = y;
            }
            if (y > maxy)
            {
                maxy = y;
            }
            if (z < minz)
            {
                minz = z;
            }
            if (z > maxz)
            {
                maxz = z;
            }
        }
        numfr++;
        /* printf("%f\t%f\t%f\n",box[XX][XX],box[YY][YY],box[ZZ][ZZ]); */

    }
    while (read_next_frame(oenv, status, &fr));

    if (bPBC)
    {
        gmx_rmpbc_done(gpbc);
    }

    if (!bCUTDOWN)
    {
        minx = miny = minz = 0;
        maxx = nbin[XX];
        maxy = nbin[YY];
        maxz = nbin[ZZ];
    }

    /* OUTPUT */
    flp = gmx_ffopen("grid.cube", "w");
    fprintf(flp, "Spatial Distribution Function\n");
    fprintf(flp, "test\n");
    fprintf(flp, "%5d%12.6f%12.6f%12.6f\n", nidxp, (MINBIN[XX]+(minx+iIGNOREOUTER)*rBINWIDTH)*10./bohr, (MINBIN[YY]+(miny+iIGNOREOUTER)*rBINWIDTH)*10./bohr, (MINBIN[ZZ]+(minz+iIGNOREOUTER)*rBINWIDTH)*10./bohr);
    fprintf(flp, "%5ld%12.6f%12.6f%12.6f\n", maxx-minx+1-(2*iIGNOREOUTER), rBINWIDTH*10./bohr, 0., 0.);
    fprintf(flp, "%5ld%12.6f%12.6f%12.6f\n", maxy-miny+1-(2*iIGNOREOUTER), 0., rBINWIDTH*10./bohr, 0.);
    fprintf(flp, "%5ld%12.6f%12.6f%12.6f\n", maxz-minz+1-(2*iIGNOREOUTER), 0., 0., rBINWIDTH*10./bohr);
    for (i = 0; i < nidxp; i++)
    {
        v = 2;
        if (*(top.atoms.atomname[indexp[i]][0]) == 'C')
        {
            v = 6;
        }
        if (*(top.atoms.atomname[indexp[i]][0]) == 'N')
        {
            v = 7;
        }
        if (*(top.atoms.atomname[indexp[i]][0]) == 'O')
        {
            v = 8;
        }
        if (*(top.atoms.atomname[indexp[i]][0]) == 'H')
        {
            v = 1;
        }
        if (*(top.atoms.atomname[indexp[i]][0]) == 'S')
        {
            v = 16;
        }
        fprintf(flp, "%5d%12.6f%12.6f%12.6f%12.6f\n", v, 0., (double)fr.x[indexp[i]][XX]*10./bohr, (double)fr.x[indexp[i]][YY]*10./bohr, (double)fr.x[indexp[i]][ZZ]*10./bohr);
    }

    tot = 0;
    for (k = 0; k < nbin[XX]; k++)
    {
        if (!(k < minx || k > maxx))
        {
            continue;
        }
        for (j = 0; j < nbin[YY]; j++)
        {
            if (!(j < miny || j > maxy))
            {
                continue;
            }
            for (i = 0; i < nbin[ZZ]; i++)
            {
                if (!(i < minz || i > maxz))
                {
                    continue;
                }
                if (bin[k][j][i] != 0)
                {
                    printf("A bin was not empty when it should have been empty. Programming error.\n");
                    printf("bin[%d][%d][%d] was = %ld\n", k, j, i, bin[k][j][i]);
                    exit(1);
                }
            }
        }
    }

    min = 999;
    max = 0;
    for (k = 0; k < nbin[XX]; k++)
    {
        if (k < minx+iIGNOREOUTER || k > maxx-iIGNOREOUTER)
        {
            continue;
        }
        for (j = 0; j < nbin[YY]; j++)
        {
            if (j < miny+iIGNOREOUTER || j > maxy-iIGNOREOUTER)
            {
                continue;
            }
            for (i = 0; i < nbin[ZZ]; i++)
            {
                if (i < minz+iIGNOREOUTER || i > maxz-iIGNOREOUTER)
                {
                    continue;
                }
                tot += bin[k][j][i];
                if (bin[k][j][i] > max)
                {
                    max = bin[k][j][i];
                }
                if (bin[k][j][i] < min)
                {
                    min = bin[k][j][i];
                }
            }
        }
    }

    numcu = (maxx-minx+1-(2*iIGNOREOUTER))*(maxy-miny+1-(2*iIGNOREOUTER))*(maxz-minz+1-(2*iIGNOREOUTER));
    if (bCALCDIV)
    {
        norm = ((double)numcu*(double)numfr) / (double)tot;
    }
    else
    {
        norm = 1.0;
    }

    for (k = 0; k < nbin[XX]; k++)
    {
        if (k < minx+iIGNOREOUTER || k > maxx-iIGNOREOUTER)
        {
            continue;
        }
        for (j = 0; j < nbin[YY]; j++)
        {
            if (j < miny+iIGNOREOUTER || j > maxy-iIGNOREOUTER)
            {
                continue;
            }
            for (i = 0; i < nbin[ZZ]; i++)
            {
                if (i < minz+iIGNOREOUTER || i > maxz-iIGNOREOUTER)
                {
                    continue;
                }
                fprintf(flp, "%12.6f ", norm*(double)bin[k][j][i]/(double)numfr);
            }
            fprintf(flp, "\n");
        }
        fprintf(flp, "\n");
    }
    gmx_ffclose(flp);

    /* printf("x=%d to %d\n",minx,maxx); */
    /* printf("y=%d to %d\n",miny,maxy); */
    /* printf("z=%d to %d\n",minz,maxz); */

    if (bCALCDIV)
    {
        printf("Counts per frame in all %ld cubes divided by %le\n", numcu, 1.0/norm);
        printf("Normalized data: average %le, min %le, max %le\n", 1.0, norm*(double)min/(double)numfr, norm*(double)max/(double)numfr);
    }
    else
    {
        printf("grid.cube contains counts per frame in all %ld cubes\n", numcu);
        printf("Raw data: average %le, min %le, max %le\n", 1.0/norm, (double)min/(double)numfr, (double)max/(double)numfr);
    }

    return 0;
}
