/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2007- The GROMACS Authors
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
#include <cstdlib>

#include <limits>

#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdspan/mdspan.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pbcutil/rmpbc.h"
#include "gromacs/topology/index.h"
#include "gromacs/topology/topology.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"

static const double bohr =
        0.529177249; /* conversion factor to compensate for VMD plugin conversion... */

int gmx_spatial(int argc, char* argv[])
{
    const char* desc[] = {
        "[THISMODULE] calculates the spatial distribution function and",
        "outputs it in a form that can be read by VMD as Gaussian98 cube format.",
        "For a system of 32,000 atoms and a 50 ns trajectory, the SDF can be generated",
        "in about 30 minutes, with most of the time dedicated to the two runs through",
        "[TT]trjconv[tt] that are required to center everything properly.",
        "This also takes a whole bunch of space (3 copies of the trajectory file).",
        "Still, the pictures are pretty and very informative when the fitted selection is ",
        "properly ",
        "made.",
        "3-4 atoms in a widely mobile group (like a free amino acid in solution) works",
        "well, or select the protein backbone in a stable folded structure to get the SDF",
        "of solvent and look at the time-averaged solvation shell.",
        "It is also possible using this program to generate the SDF based on some arbitrary",
        "Cartesian coordinate. To do that, simply omit the preliminary [gmx-trjconv] steps.",
        "",
        "Usage:",
        "",
        ("  1. Use [gmx-make_ndx] to create a group containing the atoms around which you want the "
         "SDF"),
        "  2. [TT]gmx trjconv -s a.tpr -f a.tng -o b.tng -boxcenter tric -ur compact -pbc none[tt]",
        "  3. [TT]gmx trjconv -s a.tpr -f b.tng -o c.tng -fit rot+trans[tt]",
        "  4. run [THISMODULE] on the [TT]c.tng[tt] output of step #3.",
        "  5. Load [TT]grid.cube[tt] into VMD and view as an isosurface.",
        "",
        "[BB]Note[bb] that systems such as micelles will require [TT]gmx trjconv -pbc cluster[tt] ",
        "between steps 1 and 2.",
        "",
        "Warnings",
        "^^^^^^^^",
        "",
        "The SDF will be generated for a cube that contains all bins that have some non-zero ",
        "occupancy.",
        "However, the preparatory [TT]-fit rot+trans[tt] option to [gmx-trjconv] implies that ",
        "your system will be rotating",
        "and translating in space (in order that the selected group does not). Therefore the ",
        "values that are",
        "returned will only be valid for some region around your central group/coordinate that ",
        "has full overlap",
        "with system volume throughout the entire translated/rotated system over the course of ",
        "the trajectory.",
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
    const char* bugs[] = {
        "When the allocated memory is not large enough, an error may occur "
        "suggesting the use of the [TT]-nab[tt] (Number of Additional Bins) "
        "option or increasing the [TT]-nab[tt] value."
    };

    static gmx_bool bPBC         = FALSE;
    static int      iIGNOREOUTER = -1; /*Positive values may help if the surface is spikey */
    static gmx_bool bCUTDOWN     = TRUE;
    static real     rBINWIDTH    = 0.05; /* nm */
    static gmx_bool bCALCDIV     = TRUE;
    static int      iNAB         = 16;

    t_pargs pa[] = { { "-pbc",
                       FALSE,
                       etBOOL,
                       { &bPBC },
                       "Use periodic boundary conditions for computing distances" },
                     { "-div",
                       FALSE,
                       etBOOL,
                       { &bCALCDIV },
                       "Calculate and apply the divisor for bin occupancies based on atoms/minimal "
                       "cube size. Set as TRUE for visualization and as FALSE ([TT]-nodiv[tt]) to "
                       "get accurate counts per frame" },
                     { "-ign",
                       FALSE,
                       etINT,
                       { &iIGNOREOUTER },
                       "Do not display this number of outer cubes (positive values may reduce "
                       "boundary speckles; -1 ensures outer surface is visible)" },
                     /*    { "-cut",      bCUTDOWN, etBOOL, {&bCUTDOWN},*/
                     /*      "Display a total cube that is of minimal size" }, */
                     { "-bin", FALSE, etREAL, { &rBINWIDTH }, "Width of the bins (nm)" },
                     { "-nab",
                       FALSE,
                       etINT,
                       { &iNAB },
                       "Number of additional bins to ensure proper memory allocation" } };

    double            MINBIN[3];
    double            MAXBIN[3];
    t_topology        top;
    PbcType           pbcType;
    t_trxframe        fr;
    rvec*             xtop;
    matrix            box, box_pbc;
    t_trxstatus*      status;
    int               flags = TRX_READ_X;
    t_pbc             pbc;
    t_atoms*          atoms;
    int               natoms;
    char *            grpnm, *grpnmp;
    int *             index, *indexp;
    int               nidx, nidxp;
    int               v;
    int               nbin[3];
    FILE*             flp;
    int               minx, miny, minz, maxx, maxy, maxz;
    int               numfr, numcu;
    int               maxval, minval;
    int64_t           tot;
    double            norm;
    gmx_output_env_t* oenv;
    gmx_rmpbc_t       gpbc = nullptr;

    t_filenm fnm[] = { { efTPS, nullptr, nullptr, ffREAD }, /* this is for the topology */
                       { efTRX, "-f", nullptr, ffREAD },    /* and this for the trajectory */
                       { efNDX, nullptr, nullptr, ffOPTRD } };

#define NFILE asize(fnm)

    /* This is the routine responsible for adding default options,
     * calling the X/motif interface, etc. */
    if (!parse_common_args(
                &argc, argv, PCA_CAN_TIME | PCA_CAN_VIEW, NFILE, fnm, asize(pa), pa, asize(desc), desc, asize(bugs), bugs, &oenv))
    {
        return 0;
    }

    read_tps_conf(ftp2fn(efTPS, NFILE, fnm), &top, &pbcType, &xtop, nullptr, box, TRUE);
    sfree(xtop);

    atoms = &(top.atoms);
    printf("Select group to generate SDF:\n");
    get_index(atoms, ftp2fn_null(efNDX, NFILE, fnm), 1, &nidx, &index, &grpnm);
    printf("Select group to output coords (e.g. solute):\n");
    get_index(atoms, ftp2fn_null(efNDX, NFILE, fnm), 1, &nidxp, &indexp, &grpnmp);

    /* The first time we read data is a little special */
    read_first_frame(oenv, &status, ftp2fn(efTRX, NFILE, fnm), &fr, flags);
    natoms = fr.natoms;

    /* Memory Allocation */
    MINBIN[XX] = MAXBIN[XX] = fr.x[0][XX];
    MINBIN[YY] = MAXBIN[YY] = fr.x[0][YY];
    MINBIN[ZZ] = MAXBIN[ZZ] = fr.x[0][ZZ];
    for (int i = 1; i < top.atoms.nr; ++i)
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
    for (int i = ZZ; i >= XX; --i)
    {
        MAXBIN[i] = (std::ceil((MAXBIN[i] - MINBIN[i]) / rBINWIDTH) + iNAB) * rBINWIDTH + MINBIN[i];
        MINBIN[i] -= iNAB * rBINWIDTH;
        nbin[i] = static_cast<int>(std::ceil((MAXBIN[i] - MINBIN[i]) / rBINWIDTH));
    }
    std::vector<int> binData(nbin[XX] * nbin[YY] * nbin[ZZ], 0);
    gmx::basic_mdspan<int, gmx::extents<gmx::dynamic_extent, gmx::dynamic_extent, gmx::dynamic_extent>> bin(
            binData.data(), nbin[XX], nbin[YY], nbin[ZZ]);
    copy_mat(box, box_pbc);
    numfr = 0;
    minx = miny = minz = std::numeric_limits<int>::max();
    maxx = maxy = maxz = std::numeric_limits<int>::min();

    if (bPBC)
    {
        gpbc = gmx_rmpbc_init(&top.idef, pbcType, natoms);
    }
    /* This is the main loop over frames */
    do
    {
        /* Must init pbc every step because of pressure coupling */

        copy_mat(box, box_pbc);
        if (bPBC)
        {
            gmx_rmpbc_trxfr(gpbc, &fr);
            set_pbc(&pbc, pbcType, box_pbc);
        }

        for (int i = 0; i < nidx; i++)
        {
            int x = static_cast<int>(std::floor((fr.x[index[i]][XX] - MINBIN[XX]) / rBINWIDTH));
            int y = static_cast<int>(std::floor((fr.x[index[i]][YY] - MINBIN[YY]) / rBINWIDTH));
            int z = static_cast<int>(std::floor((fr.x[index[i]][ZZ] - MINBIN[ZZ]) / rBINWIDTH));
            if (x < 0 || x >= nbin[XX] || y < 0 || y >= nbin[YY] || z < 0 || z >= nbin[ZZ])
            {
                printf("There was an item outside of the allocated memory. Increase the value "
                       "given with the -nab option.\n");
                printf("Memory was allocated for [%f,%f,%f]\tto\t[%f,%f,%f]\n",
                       MINBIN[XX],
                       MINBIN[YY],
                       MINBIN[ZZ],
                       MAXBIN[XX],
                       MAXBIN[YY],
                       MAXBIN[ZZ]);
                printf("Memory was required for [%f,%f,%f]\n",
                       fr.x[index[i]][XX],
                       fr.x[index[i]][YY],
                       fr.x[index[i]][ZZ]);
                exit(1);
            }

            bin[x][y][z]++;
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

    } while (read_next_frame(oenv, status, &fr));

    if (bPBC)
    {
        gmx_rmpbc_done(gpbc);
    }

    if (!bCUTDOWN)
    {
        minx = miny = minz = 0;
        maxx               = nbin[XX] - 1;
        maxy               = nbin[YY] - 1;
        maxz               = nbin[ZZ] - 1;
    }

    iIGNOREOUTER = std::max(iIGNOREOUTER, 0);
    int outputStarts[DIM], outputEnds[DIM];

    outputStarts[XX] = minx + iIGNOREOUTER;
    outputStarts[YY] = miny + iIGNOREOUTER;
    outputStarts[ZZ] = minz + iIGNOREOUTER;
    outputEnds[XX]   = maxx - iIGNOREOUTER;
    outputEnds[YY]   = maxy - iIGNOREOUTER;
    outputEnds[ZZ]   = maxz - iIGNOREOUTER;

    /* OUTPUT */
    flp = gmx_ffopen("grid.cube", "w");
    fprintf(flp, "Spatial Distribution Function\n");
    fprintf(flp, "test\n");
    /*
      Values in .cube file represent the density at the grid point.
     Corresponding coordinates to the binned value is the center of the bin, i.e. + 0.5 to the bin index.
   */
    fprintf(flp,
            "%5d%12.6f%12.6f%12.6f\n",
            nidxp,
            (MINBIN[XX] + (outputStarts[XX] + 0.5) * rBINWIDTH) * 10. / bohr,
            (MINBIN[YY] + (outputStarts[YY] + 0.5) * rBINWIDTH) * 10. / bohr,
            (MINBIN[ZZ] + (outputStarts[ZZ] + 0.5) * rBINWIDTH) * 10. / bohr);
    fprintf(flp, "%5d%12.6f%12.6f%12.6f\n", outputEnds[XX] - outputStarts[XX], rBINWIDTH * 10. / bohr, 0., 0.);
    fprintf(flp, "%5d%12.6f%12.6f%12.6f\n", outputEnds[YY] - outputStarts[YY], 0., rBINWIDTH * 10. / bohr, 0.);
    fprintf(flp, "%5d%12.6f%12.6f%12.6f\n", outputEnds[ZZ] - outputStarts[ZZ], 0., 0., rBINWIDTH * 10. / bohr);

    for (int i = 0; i < nidxp; i++)
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
        fprintf(flp,
                "%5d%12.6f%12.6f%12.6f%12.6f\n",
                v,
                0.,
                fr.x[indexp[i]][XX] * 10.0 / bohr,
                fr.x[indexp[i]][YY] * 10.0 / bohr,
                fr.x[indexp[i]][ZZ] * 10.0 / bohr);
    }

    tot = 0;
    for (int i = 0; i < nbin[XX]; i++)
    {
        if (!(i < minx || i > maxx))
        {
            continue;
        }
        for (int j = 0; j < nbin[YY]; j++)
        {
            if (!(j < miny || j > maxy))
            {
                continue;
            }
            for (int k = 0; k < nbin[ZZ]; k++)
            {
                if (!(k < minz || k > maxz))
                {
                    continue;
                }
                int binValue = bin[i][j][k];
                GMX_RELEASE_ASSERT(
                        binValue == 0,
                        gmx::formatString("A bin was not empty when it should have been empty. "
                                          "Programming error.\n bin[%d][%d][%d] was = %d\n",
                                          i,
                                          j,
                                          k,
                                          binValue)
                                .c_str());
            }
        }
    }

    minval = 999;
    maxval = 0;
    for (int i = outputStarts[XX]; i < outputEnds[XX]; i++)
    {
        for (int j = outputStarts[YY]; j < outputEnds[YY]; j++)
        {
            for (int k = outputStarts[ZZ]; k < outputEnds[ZZ]; k++)
            {
                int binValue = bin[i][j][k];
                tot += binValue;
                if (binValue > maxval)
                {
                    maxval = binValue;
                }
                if (binValue < minval)
                {
                    minval = binValue;
                }
            }
        }
    }

    numcu = (outputEnds[XX] - outputStarts[XX]) * (outputEnds[YY] - outputStarts[YY])
            * (outputEnds[ZZ] - outputStarts[ZZ]);
    if (bCALCDIV)
    {
        norm = double(numcu) * numfr / tot;
        GMX_ASSERT(norm >= 0, "The norm should be non-negative.");
    }
    else
    {
        norm = 1.0;
    }

    for (int i = outputStarts[XX]; i < outputEnds[XX]; i++)
    {
        for (int j = outputStarts[YY]; j < outputEnds[YY]; j++)
        {
            for (int k = outputStarts[ZZ]; k < outputEnds[ZZ]; k++)
            {
                fprintf(flp, "%12.6f ", static_cast<double>(norm * bin[i][j][k]) / numfr);
            }
            fprintf(flp, "\n");
        }
        fprintf(flp, "\n");
    }
    gmx_ffclose(flp);

    if (bCALCDIV)
    {
        printf("Counts per frame in all %d cubes divided by %le\n", numcu, 1.0 / norm);
        printf("Normalized data: average %le, min %le, max %le\n",
               1.0,
               minval * norm / numfr,
               maxval * norm / numfr);
    }
    else
    {
        printf("grid.cube contains counts per frame in all %d cubes\n", numcu);
        printf("Raw data: average %le, min %le, max %le\n",
               static_cast<double>(tot) / numfr / numcu,
               static_cast<double>(minval) / numfr,
               static_cast<double>(maxval) / numfr);
    }

    return 0;
}
