/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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

#include <stdio.h>

#include <cstring>

#include "gromacs/applied-forces/maputil.h"
#include "gromacs/applied-forces/densityfitting/densfit.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/oenv.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/fileio/griddataio.h"
#include "gromacs/fileio/mrcmetadata.h"
#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/index.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/pleasecite.h"
#include "gromacs/utility/smalloc.h"

#define round(a) (int)(a+0.5)

///* Can be used to set density data */
//static void transform(t_mapdata *map)
//{
//    int i;
//
//    for (i=0; i<map->grid[XX]; i++)
//        map->d[0][0][i] = 1.0;
//
//    for (i=0; i<map->grid[YY]; i++)
//        map->d[0][i][0] = 1.0;
//
//    for (i=0; i<map->grid[ZZ]; i++)
//        map->d[i][0][0] = 1.0;
//
//}



///* Put the origin in the center of the protein */
//    clear_rvec(com);
//    sum_of_weights=0.0;
//    for (i=0; i<natoms; i++)
//    {
//        svmul(w[i], x[i], add);
//        rvec_inc(com, add);
//        sum_of_weights += w[i];
//    }
//    /* Divide by the sum of weights ("the total mass") */
//    svmul(1.0/sum_of_weights, com, com);
//
//    /* Subtract the center */
//    for (i=0; i<natoms; i++)
//        rvec_dec(x[i], com);
//
//    fprintf(stdout, "Protein center was at %g %g %g nm.\n", com[XX], com[YY], com[ZZ]);



//! \brief Do we have a reasonable box?
static void verify_box(int ePBC, matrix box)
{
    int i;


    /* Check whether we have box values at all, at least the diagonal elements
     * must have values */
    for (i = 0; i < 3; i++)
    {
        if (box[i][i] < GMX_REAL_EPS)
        {
            gmx_fatal(FARGS, "Invalid box! Diagonal is %g %g %g.\n",
                      box[XX][XX], box[YY][YY], box [ZZ][ZZ]);
        }
    }

    /* Make shure the box makes sense */
    check_box(ePBC, box);
}


/* Returns TRUE if we are dealing with a whole trajectory, and FALSE if
 * the format just supports a single structure */
static gmx_bool is_trajectory(const char *fn)
{
    gmx_bool bTraj;

    int      ftp = fn2ftp(fn);


    switch (ftp)
    {
        case efXTC:
        case efTRX:
        case efTRO:
        case efTRN:
        case efTRR:
//    case efTRJ:
            bTraj = TRUE;
            break;
        default:
            bTraj = FALSE;
            break;
    }

    return bTraj;
}


static void print_matrix(matrix m, FILE *out)
{
    int i;


    for (i = 0; i < DIM; i++)
    {
        fprintf(out, "%5.3f  %5.3f  %5.3f\n", m[XX][i], m[YY][i], m[ZZ][i]);
    }
}

/* Only allow for input file combinations that are tested: */
static void check_input_combinations(int nfile, const t_filenm fnm[])
{
    if (opt2bSet("-f", nfile, fnm) && opt2bSet("-s", nfile, fnm))
    {
        gmx_fatal(FARGS, "Please use either -f or -s");
    }

}


int gmx_map(int argc, char *argv[])
{
    const char *desc[] = {
        "[THISMODULE] is a tool to read in and write out (electron) density maps.",
        "With this tool you can produce a density map from an input ",
        "coordinate file to be embedded in a [REF].tpr[ref] file with [gmx-grompp] for",
        "example.[PAR]",
        "Possible uses are:",
        "",
        "* Provide an input structure with [TT]-f[tt] and output a density map based on",
        "  the grid settings. Note that you can also input a whole trajectory, in that",
        "  case, a series of map files will be output. Use and index file to spread just",
        "  a part of the atoms",
        "* Provide a map with [TT]-mi[tt] and and output some map characteristics",
        "* Provide a [REF].tpr[ref] file with [TT]-s[tt] and output the embedded map with [TT]-mo[tt]",
        ""
    };

    const char *bugs[] =
    {
        "This tool is under construction :)"
    };

    /* Command line options ! */
    gmx_bool bVerbose    = FALSE;
    gmx_bool bRemoveH    = FALSE;
    gmx_bool bListHeader = FALSE;
    gmx_bool bTranslate  = TRUE;  /* Translate PDB to first quadrant          */
    gmx_bool bPositive   = FALSE; /* Translate all voxels to values >= 0      */

    real     sigma         = 0.4; /* Gaussian width (nm) used to tranform the
                                     discrete atomic positions into a density */
    real     grid_spacing  = 0.2; /* Spacing (nm) of the density grid         */
    real     ref_spacing   = 0.0; /* Spacing (nm) determined from map         */
    real     sigma_dist    = 4.0; /* Calculate the Gaussians in this range    */
    real     scale         = 1.0; /* Rescaling factor for x, y, and z dim     */
    int      ndigit        = 0;
    int      skip          = 1;


    t_filenm fnm[] = {
        { efDENSITYMAP, "-mi", "ccp4in", ffOPTRD },    /* CCP4 density map input file */
        { efDENSITYMAP, "-mo", "ccp4out", ffOPTWR },   /* CCP4 density map output file */
        { efSTX, "-f", NULL, ffOPTRD },
        { efPDB, "-o", "processed", ffOPTWR },         /* Optionally processed PDB file */
        { efNDX, "-n", NULL, ffOPTRD },                /* Optionally select just a subset of the atoms to spread */
        { efTPR, NULL, NULL, ffOPTRD }                 /* You can also output the map from the .tpr file */
    };
#define NFILE asize(fnm)

    t_pargs pa[] = {
        { "-v",      FALSE, etBOOL, {&bVerbose},
          "Be loud and noisy" },
        { "-list",   FALSE, etBOOL, {&bListHeader},
          "List density fitting parameters in [TT].tpr[tt] and map files" },
        { "-ignh",   FALSE, etBOOL, {&bRemoveH},
          "Ignore hydrogen atoms that are in the PDB coordinate file" },
        { "-sigma",  FALSE, etREAL, {&sigma},
          "Create a simulated density by replacing the atoms by Gaussian functions of width sigma (nm)" },
        { "-dist",   FALSE, etREAL, {&sigma_dist},
          "Only calculate Gaussians within this distance of an atomic position (sigmas),"
          "i.e. 4 means +/- 4 times sigma."},
        { "-sp",     FALSE, etREAL, {&grid_spacing},
          "Spacing of the density grid (nm)" },
        { "-trans",  FALSE, etBOOL, {&bTranslate},
          "Translate PDB positions to first quadrant" },
        { "-skip",   FALSE, etINT,  {&skip},
          "Only analyse every nr-th frame" },
        { "-nzero",  FALSE, etINT,  {&ndigit},
          "When writing separate output maps, use these many digits "
          "for the file numbers and prepend zeros as needed" },
        { "-scale",  FALSE, etREAL, {&scale},
          "Rescale the input density map resolution in each dimension by this factor" },
        { "-positive", FALSE, etBOOL, {&bPositive},
          "If negative, subtract the smallest voxel value from all voxels. This will result in a map with all voxel values >= 0."}
    };
#define NPA asize(pa)

    char              title[10000];
    char              fn_map[STRLEN];
    int               i, ePBC;
    matrix            box;
    t_inputrec        ir;
    t_state           state;
    gmx_mtop_t        mtop;
    gmx_bool          bMoreFrames;

    gmx_bool          bWroteMap       = FALSE;
    t_trxstatus      *status          = NULL;
    t_atoms          *atoms           = NULL;
    gmx_output_env_t *oenv            = NULL;
    char             *grpname         = NULL;
    const char       *fn_in           = NULL;
    rvec             *pos             = NULL; /* Positions, part of these are
                                                 for spreading                  */
    int              *ind             = NULL; /* Indices of atoms for spreading */
    int               natoms          = 0;    /* Size of the x array            */
    int               iframe          = 0;
    int               iout            = 0;
    real              time            = 0;


    /* Parse the command line arguments */
    if (!parse_common_args(&argc, argv, PCA_CAN_VIEW, NFILE, fnm, NPA, pa,
                           asize(desc), desc, asize(bugs), bugs, &oenv))
    {
        return 0;
    }

    /* Only allow for input file combinations that are tested: */
    check_input_combinations(NFILE, fnm);

    fn_in = ftp2fn(efSTX, NFILE, fnm);

    /**************************************************************************/
    /* 1. CHECK FILES AND READ IN THE DATA ************************************/
    /**************************************************************************/

    /* Read coordinates */
    if (opt2bSet("-f", NFILE, fnm) )
    {
        if (is_trajectory(fn_in))
        {
            natoms = read_first_x(oenv, &status, fn_in, &time, &pos, box);
        }
        else
        {
            bool haveTopology;
            readConfAndTopology(fn_in, &haveTopology, &mtop, &ePBC, &pos, NULL, box);
        }

        /* Read in index file, if supplied */
        if (opt2bSet("-n", NFILE, fnm) )
        {
            if (gmx_fexist(opt2fn("-n", NFILE, fnm)) )
            {
                fprintf(stdout, "Select the atoms for the density map:\n");
                get_index(atoms, opt2fn("-n", NFILE, fnm), 1, &natoms, &ind, &grpname);
            }
            else
            {
                gmx_fatal(FARGS, "Index file %s not found!\n", opt2fn("-n", NFILE, fnm));
            }
        }
    }

    /* Read the .tpr file, if supplied */
    if (opt2bSet("-s", NFILE, fnm) )
    {
        read_tpx_state( opt2fn("-s", NFILE, fnm), &ir, &state, &mtop);

        snew(atoms, 1);
        *atoms = gmx_mtop_global_atoms(&mtop);

        /* For the map embedded in the .tpr file, print header and output map to file */
        if (ir.bDensityFitting)
        {
            if (bListHeader)
            {
                fprintf(stdout, "Density map parameters:\n");
                pr_density_fitting(stdout, 3, *(ir.densfit));
            }
            if (opt2bSet("-mo", NFILE, fnm) )
            {
                fprintf(stderr, "Will output the map that is embedded in the .tpr file to %s\n",
                        opt2fn("-mo", NFILE, fnm) );
            }
        }
    }

    /* Initialize global and local density fitting data structures */
    if (ir.densfit == nullptr)
    {
        ir.densfit    = std::unique_ptr<gmx::Densfit>(new gmx::Densfit(sigma, sigma_dist, 0, natoms, ind, grid_spacing, bVerbose));
    }

    auto  &densfit = *(ir.densfit);

    /* Was an input density map provided? */
    if (opt2bSet("-mi", NFILE, fnm) )
    {
        /* Read input file (map) */
        gmx::MrcMetaData metadata;
        auto             result = gmx::MrcFile().read_with_meta(opt2fn("-mi", NFILE, fnm), &metadata);
        densfit.map_ref = gmx::gridDataToMapData(result);
        fprintf(stdout, "%s\n", metadata.to_string().c_str());

        ref_spacing = gmx::get_map_spacing(densfit.map_ref, stderr);
    }


    /**************************************************************************/
    /* 2. DO WHAT NEEDS TO BE DONE ********************************************/
    /**************************************************************************/

    if (opt2bSet("-f", NFILE, fnm) )
    {
        please_cite(stdout, "Tama2008");

        if (densfit.map_ref.vox.size() == 0)
        {
            /* If no reference map is provided, create a cell and allocate a grid
             * for the spreading routine to use */

            /* Do we have a reasonable box? */
            verify_box(ir.ePBC, box);

            for (i = 0; i < 3; i++)
            {
                densfit.map_ref.cell   [i]   = NM2A*box[i][i];
                /* Set box angles to a default of 90 degrees */
                densfit.map_ref.cell   [i+3] = 90.0;
                densfit.map_ref.map_dim[i]   = densfit.map_ref.cell[i] / (NM2A*grid_spacing);
                /* Choose grid equal to map_dim */
                densfit.map_ref.grid   [i] = densfit.map_ref.map_dim[i];
                /* Set a default axes order */
                densfit.map_ref.axes_order[i] = i+1;
            }
            gmx::allocate_density_grid(densfit.map_ref.map_dim, &densfit.map_ref.vox);
        }
        else
        {
            fprintf(stderr, "Spreading atomic positions on the provided input map.\n");
            fprintf(stderr, "Old box from input structure:\n");
            print_matrix(box, stderr);
//            put_atoms_in_box(ePBC, box, natoms, pos);
            clear_mat(box);
            for (i = 0; i < 3; i++)
            {
                box[i][i] = A2NM * densfit.map_ref.cell[i];
            }
            fprintf(stderr, "Using box from the reference map:\n");
            print_matrix(box, stderr);
            if (opt2bSet("-s", NFILE, fnm) )
            {
                write_sto_conf(opt2fn("-o", NFILE, fnm), title, atoms, pos, NULL, ePBC, box);
            }
            fprintf(stderr, "Spacing: %g\n", ref_spacing);
        }

        /* The simulated map copies the grid settings from the reference map. */
        gmx::new_map_sim_from_ref(&densfit, densfit.map_ref);


        /* We may have a whole trajectory of frames to work on now */
        strncpy(fn_map, opt2fn("-mo", NFILE, fnm), STRLEN);
        bMoreFrames = FALSE;
        iframe      = 0;
        iout        = 0;
        do
        {
            if (0 == iframe % skip)
            {
                /* Construct the simulated density! */
                densfit.assemble_atoms_for_spread(pos);

                /* The reference map grows and shrinks with the box */
                if (opt2bSet("-mi", NFILE, fnm))
                {
                    gmx::get_map_spacing(densfit.map_ref, stdout);
                    densfit.couple_map_spacing_to_box(box);
                    gmx::get_map_spacing(densfit.map_ref, stdout);
                }

                densfit.spread_atoms(box);

                /* Calculate the correlation coefficient to the input map */
                if (opt2bSet("-mi", NFILE, fnm) )
                {
                    densfit.calc_correlation_coeff(stdout);
                }

                /* If needed, attach frame numbers to output files */
                if (is_trajectory(fn_in))
                {
                    gmx::make_filename(opt2fn("-mo", NFILE, fnm), ndigit, iout++, fn_map);
                }

                /* Write out simulated map */
                auto gridData = gmx::mapDataToGridData(densfit.map_sim);
                gmx::MrcFile().write(fn_map, gridData);
                if (bVerbose)
                {
                    gmx::MrcMetaData meta;
                    meta.fromGrid(gridData.getGrid());
                    meta.set_grid_stats(gridData);
                    fprintf(stdout, "%s\n", meta.to_string().c_str());
                }
                bWroteMap = TRUE;
            }
            /* Read the next frame */
            if (is_trajectory(fn_in))
            {
                bMoreFrames = read_next_x(oenv, status, &time, pos, box);
            }
            iframe++;

        }
        while (bMoreFrames);

        if (is_trajectory(fn_in))
        {
            close_trx(status);
        }
    }
    else if (opt2bSet("-mi", NFILE, fnm) && ((scale != 1.0) || bPositive ) )
    {
        if (bPositive)
        {
            fprintf(stdout, "\nWill translate the voxels to positive values.\n");
            gmx::make_positive(&(densfit.map_ref));
        }
        if (scale != 1.0)
        {
            fprintf(stdout, "\nWill rescale the input density map by a factor of %g in each dimension.\n", scale);
            densfit.map_ref = gmx::rescale_map(densfit.map_ref, scale);
        }
    }

    /* Output the density map */
    if (opt2bSet("-mo", NFILE, fnm) && !bWroteMap)
    {
        if (densfit.map_ref.vox.size() != 0)
        {
            /* Write out map that was provided with -mi or -s files */
            auto gridData = gmx::mapDataToGridData(densfit.map_ref);
            gmx::MrcFile().write(opt2fn("-mo", NFILE, fnm), gridData);
            if (bVerbose)
            {
                gmx::MrcMetaData meta;
                meta.fromGrid(gridData.getGrid());
                meta.set_grid_stats(gridData);
                fprintf(stdout, "%s\n", meta.to_string().c_str());
            }
        }
        else
        {
            /* What went wrong? */
            if (opt2bSet("-s", NFILE, fnm) )
            {
                if (ir.bDensityFitting)
                {
                    gmx_fatal(FARGS, "Could not read density fitting data from .tpr file.");
                }
                else
                {
                    fprintf(stderr, "This .tpr file does not contain a density map!\n");
                }
            }
            else
            {
                fprintf(stderr, "WARNING: Don't know what to write into file '%s'\n", opt2fn("-mo", NFILE, fnm));
                fprintf(stderr, "         Please provide an input file containing a map (with -s or -mi) or generate one with -f.\n");
            }
        }
    }

    return 0;
}
