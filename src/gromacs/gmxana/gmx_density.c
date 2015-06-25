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

#include <ctype.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/gmxana/gstat.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/legacyheaders/viewit.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pbcutil/rmpbc.h"
#include "gromacs/topology/index.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"

typedef struct {
    char *atomname;
    int   nr_el;
} t_electron;

/****************************************************************************/
/* This program calculates the partial density across the box.              */
/* Peter Tieleman, Mei 1995                                                 */
/****************************************************************************/

/* used for sorting the list */
int compare(void *a, void *b)
{
    t_electron *tmp1, *tmp2;
    tmp1 = (t_electron *)a; tmp2 = (t_electron *)b;

    return strcmp(tmp1->atomname, tmp2->atomname);
}

int get_electrons(t_electron **eltab, const char *fn)
{
    char  buffer[256];  /* to read in a line   */
    char  tempname[80]; /* buffer to hold name */
    int   tempnr;

    FILE *in;
    int   nr;        /* number of atomstypes to read */
    int   i;

    if (!(in = gmx_ffopen(fn, "r")))
    {
        gmx_fatal(FARGS, "Couldn't open %s. Exiting.\n", fn);
    }

    if (NULL == fgets(buffer, 255, in))
    {
        gmx_fatal(FARGS, "Error reading from file %s", fn);
    }

    if (sscanf(buffer, "%d", &nr) != 1)
    {
        gmx_fatal(FARGS, "Invalid number of atomtypes in datafile\n");
    }

    snew(*eltab, nr);

    for (i = 0; i < nr; i++)
    {
        if (fgets(buffer, 255, in) == NULL)
        {
            gmx_fatal(FARGS, "reading datafile. Check your datafile.\n");
        }
        if (sscanf(buffer, "%s = %d", tempname, &tempnr) != 2)
        {
            gmx_fatal(FARGS, "Invalid line in datafile at line %d\n", i+1);
        }
        (*eltab)[i].nr_el    = tempnr;
        (*eltab)[i].atomname = gmx_strdup(tempname);
    }
    gmx_ffclose(in);

    /* sort the list */
    fprintf(stderr, "Sorting list..\n");
    qsort ((void*)*eltab, nr, sizeof(t_electron),
           (int(*)(const void*, const void*))compare);

    return nr;
}

void center_coords(t_atoms *atoms, atom_id *index_center, int ncenter,
                   matrix box, rvec x0[])
{
    int  i, k, m;
    real tmass, mm;
    rvec com, shift, box_center;

    tmass = 0;
    clear_rvec(com);
    for (k = 0; (k < ncenter); k++)
    {
        i = index_center[k];
        if (i >= atoms->nr)
        {
            gmx_fatal(FARGS, "Index %d refers to atom %d, which is larger than natoms (%d).",
                      k+1, i+1, atoms->nr);
        }
        mm     = atoms->atom[i].m;
        tmass += mm;
        for (m = 0; (m < DIM); m++)
        {
            com[m] += mm*x0[i][m];
        }
    }
    for (m = 0; (m < DIM); m++)
    {
        com[m] /= tmass;
    }
    calc_box_center(ecenterDEF, box, box_center);
    rvec_sub(com, box_center, shift);

    /* Important - while the center was calculated based on a group, we should move all atoms */
    for (i = 0; (i < atoms->nr); i++)
    {
        rvec_dec(x0[i], shift);
    }
}

void calc_electron_density(const char *fn, atom_id **index, int gnx[],
                           double ***slDensity, int *nslices, t_topology *top,
                           int ePBC,
                           int axis, int nr_grps, real *slWidth,
                           t_electron eltab[], int nr, gmx_bool bCenter,
                           atom_id *index_center, int ncenter,
                           gmx_bool bRelative, const output_env_t oenv)
{
    rvec        *x0;            /* coordinates without pbc */
    matrix       box;           /* box (3x3) */
    double       invvol;
    int          natoms;        /* nr. atoms in trj */
    t_trxstatus *status;
    int          i, n,          /* loop indices */
                 nr_frames = 0, /* number of frames */
                 slice;         /* current slice */
    t_electron  *found;         /* found by bsearch */
    t_electron   sought;        /* thingie thought by bsearch */
    real         boxSz, aveBox;
    gmx_rmpbc_t  gpbc = NULL;

    real         t,
                 z;

    if (axis < 0 || axis >= DIM)
    {
        gmx_fatal(FARGS, "Invalid axes. Terminating\n");
    }

    if ((natoms = read_first_x(oenv, &status, fn, &t, &x0, box)) == 0)
    {
        gmx_fatal(FARGS, "Could not read coordinates from statusfile\n");
    }

    aveBox = 0;

    if (!*nslices)
    {
        *nslices = (int)(box[axis][axis] * 10); /* default value */
        fprintf(stderr, "\nDividing the box in %d slices\n", *nslices);
    }

    snew(*slDensity, nr_grps);
    for (i = 0; i < nr_grps; i++)
    {
        snew((*slDensity)[i], *nslices);
    }

    gpbc = gmx_rmpbc_init(&top->idef, ePBC, top->atoms.nr);
    /*********** Start processing trajectory ***********/
    do
    {
        gmx_rmpbc(gpbc, natoms, box, x0);

        /* Translate atoms so the com of the center-group is in the
         * box geometrical center.
         */
        if (bCenter)
        {
            center_coords(&top->atoms, index_center, ncenter, box, x0);
        }

        invvol   = *nslices/(box[XX][XX]*box[YY][YY]*box[ZZ][ZZ]);

        if (bRelative)
        {
            *slWidth = 1.0/(*nslices);
            boxSz    = 1.0;
        }
        else
        {
            *slWidth = box[axis][axis]/(*nslices);
            boxSz    = box[axis][axis];
        }

        aveBox += box[axis][axis];

        for (n = 0; n < nr_grps; n++)
        {
            for (i = 0; i < gnx[n]; i++) /* loop over all atoms in index file */
            {
                z = x0[index[n][i]][axis];
                while (z < 0)
                {
                    z += box[axis][axis];
                }
                while (z > box[axis][axis])
                {
                    z -= box[axis][axis];
                }

                if (bRelative)
                {
                    z = z/box[axis][axis];
                }

                /* determine which slice atom is in */
                if (bCenter)
                {
                    slice = floor( (z-(boxSz/2.0)) / (*slWidth) ) + *nslices/2;
                }
                else
                {
                    slice = (z / (*slWidth));
                }
                sought.nr_el    = 0;
                sought.atomname = gmx_strdup(*(top->atoms.atomname[index[n][i]]));

                /* now find the number of electrons. This is not efficient. */
                found = (t_electron *)
                    bsearch((const void *)&sought,
                            (const void *)eltab, nr, sizeof(t_electron),
                            (int(*)(const void*, const void*))compare);

                if (found == NULL)
                {
                    fprintf(stderr, "Couldn't find %s. Add it to the .dat file\n",
                            *(top->atoms.atomname[index[n][i]]));
                }
                else
                {
                    (*slDensity)[n][slice] += (found->nr_el -
                                               top->atoms.atom[index[n][i]].q)*invvol;
                }
                free(sought.atomname);
            }
        }
        nr_frames++;
    }
    while (read_next_x(oenv, status, &t, x0, box));
    gmx_rmpbc_done(gpbc);

    /*********** done with status file **********/
    close_trj(status);

/* slDensity now contains the total number of electrons per slice, summed
   over all frames. Now divide by nr_frames and volume of slice
 */

    fprintf(stderr, "\nRead %d frames from trajectory. Counting electrons\n",
            nr_frames);

    if (bRelative)
    {
        aveBox  /= nr_frames;
        *slWidth = aveBox/(*nslices);
    }

    for (n = 0; n < nr_grps; n++)
    {
        for (i = 0; i < *nslices; i++)
        {
            (*slDensity)[n][i] /= nr_frames;
        }
    }

    sfree(x0); /* free memory used by coordinate array */
}

void calc_density(const char *fn, atom_id **index, int gnx[],
                  double ***slDensity, int *nslices, t_topology *top, int ePBC,
                  int axis, int nr_grps, real *slWidth, gmx_bool bCenter,
                  atom_id *index_center, int ncenter,
                  gmx_bool bRelative, const output_env_t oenv)
{
    rvec        *x0;            /* coordinates without pbc */
    matrix       box;           /* box (3x3) */
    double       invvol;
    int          natoms;        /* nr. atoms in trj */
    t_trxstatus *status;
    int        **slCount,       /* nr. of atoms in one slice for a group */
                 i, j, n,       /* loop indices */
                 ax1       = 0, ax2 = 0,
                 nr_frames = 0, /* number of frames */
                 slice;         /* current slice */
    real         t,
                 z;
    real         boxSz, aveBox;
    char        *buf;    /* for tmp. keeping atomname */
    gmx_rmpbc_t  gpbc = NULL;

    if (axis < 0 || axis >= DIM)
    {
        gmx_fatal(FARGS, "Invalid axes. Terminating\n");
    }

    if ((natoms = read_first_x(oenv, &status, fn, &t, &x0, box)) == 0)
    {
        gmx_fatal(FARGS, "Could not read coordinates from statusfile\n");
    }

    aveBox = 0;

    if (!*nslices)
    {
        *nslices = (int)(box[axis][axis] * 10); /* default value */
        fprintf(stderr, "\nDividing the box in %d slices\n", *nslices);
    }

    snew(*slDensity, nr_grps);
    for (i = 0; i < nr_grps; i++)
    {
        snew((*slDensity)[i], *nslices);
    }

    gpbc = gmx_rmpbc_init(&top->idef, ePBC, top->atoms.nr);
    /*********** Start processing trajectory ***********/
    do
    {
        gmx_rmpbc(gpbc, natoms, box, x0);

        /* Translate atoms so the com of the center-group is in the
         * box geometrical center.
         */
        if (bCenter)
        {
            center_coords(&top->atoms, index_center, ncenter, box, x0);
        }

        invvol   = *nslices/(box[XX][XX]*box[YY][YY]*box[ZZ][ZZ]);

        if (bRelative)
        {
            *slWidth = 1.0/(*nslices);
            boxSz    = 1.0;
        }
        else
        {
            *slWidth = box[axis][axis]/(*nslices);
            boxSz    = box[axis][axis];
        }

        aveBox += box[axis][axis];

        for (n = 0; n < nr_grps; n++)
        {
            for (i = 0; i < gnx[n]; i++) /* loop over all atoms in index file */
            {
                z = x0[index[n][i]][axis];
                while (z < 0)
                {
                    z += box[axis][axis];
                }
                while (z > box[axis][axis])
                {
                    z -= box[axis][axis];
                }

                if (bRelative)
                {
                    z = z/box[axis][axis];
                }

                /* determine which slice atom is in */
                if (bCenter)
                {
                    slice = floor( (z-(boxSz/2.0)) / (*slWidth) ) + *nslices/2;
                }
                else
                {
                    slice = floor(z / (*slWidth));
                }

                /* Slice should already be 0<=slice<nslices, but we just make
                 * sure we are not hit by IEEE rounding errors since we do
                 * math operations after applying PBC above.
                 */
                if (slice < 0)
                {
                    slice += *nslices;
                }
                else if (slice >= *nslices)
                {
                    slice -= *nslices;
                }

                (*slDensity)[n][slice] += top->atoms.atom[index[n][i]].m*invvol;
            }
        }
        nr_frames++;
    }
    while (read_next_x(oenv, status, &t, x0, box));
    gmx_rmpbc_done(gpbc);

    /*********** done with status file **********/
    close_trj(status);

    /* slDensity now contains the total mass per slice, summed over all
       frames. Now divide by nr_frames and volume of slice
     */

    fprintf(stderr, "\nRead %d frames from trajectory. Calculating density\n",
            nr_frames);

    if (bRelative)
    {
        aveBox  /= nr_frames;
        *slWidth = aveBox/(*nslices);
    }

    for (n = 0; n < nr_grps; n++)
    {
        for (i = 0; i < *nslices; i++)
        {
            (*slDensity)[n][i] /= nr_frames;
        }
    }

    sfree(x0); /* free memory used by coordinate array */
}

void plot_density(double *slDensity[], const char *afile, int nslices,
                  int nr_grps, char *grpname[], real slWidth,
                  const char **dens_opt,
                  gmx_bool bCenter, gmx_bool bRelative, gmx_bool bSymmetrize,
                  const output_env_t oenv)
{
    FILE       *den;
    const char *title  = NULL;
    const char *xlabel = NULL;
    const char *ylabel = NULL;
    int         slice, n;
    real        ddd;
    real        axispos;

    title = bSymmetrize ? "Symmetrized partial density" : "Partial density";

    if (bCenter)
    {
        xlabel = bRelative ?
            "Average relative position from center (nm)" :
            "Relative position from center (nm)";
    }
    else
    {
        xlabel = bRelative ? "Average coordinate (nm)" : "Coordinate (nm)";
    }

    switch (dens_opt[0][0])
    {
        case 'm': ylabel = "Density (kg m\\S-3\\N)"; break;
        case 'n': ylabel = "Number density (nm\\S-3\\N)"; break;
        case 'c': ylabel = "Charge density (e nm\\S-3\\N)"; break;
        case 'e': ylabel = "Electron density (e nm\\S-3\\N)"; break;
    }

    den = xvgropen(afile,
                   title, xlabel, ylabel, oenv);

    xvgr_legend(den, nr_grps, (const char**)grpname, oenv);

    for (slice = 0; (slice < nslices); slice++)
    {
        if (bCenter)
        {
            axispos = (slice - nslices/2.0 + 0.5) * slWidth;
        }
        else
        {
            axispos = (slice + 0.5) * slWidth;
        }
        fprintf(den, "%12g  ", axispos);
        for (n = 0; (n < nr_grps); n++)
        {
            if (bSymmetrize)
            {
                ddd = (slDensity[n][slice]+slDensity[n][nslices-slice-1])*0.5;
            }
            else
            {
                ddd = slDensity[n][slice];
            }
            if (dens_opt[0][0] == 'm')
            {
                fprintf(den, "   %12g", ddd*AMU/(NANO*NANO*NANO));
            }
            else
            {
                fprintf(den, "   %12g", ddd);
            }
        }
        fprintf(den, "\n");
    }

    xvgrclose(den);
}

int gmx_density(int argc, char *argv[])
{
    const char        *desc[] = {
        "[THISMODULE] computes partial densities across the box, using an index file.[PAR]",
        "For the total density of NPT simulations, use [gmx-energy] instead.",
        "[PAR]",

        "Option [TT]-center[tt] performs the histogram binning relative to the center",
        "of an arbitrary group, in absolute box coordinates. If you are calculating",
        "profiles along the Z axis box dimension bZ, output would be from -bZ/2 to",
        "bZ/2 if you center based on the entire system.",
        "Note that this behaviour has changed in GROMACS 5.0; earlier versions",
        "merely performed a static binning in (0,bZ) and shifted the output. Now",
        "we compute the center for each frame and bin in (-bZ/2,bZ/2).[PAR]",

        "Option [TT]-symm[tt] symmetrizes the output around the center. This will",
        "automatically turn on [TT]-center[tt] too.",

        "Option [TT]-relative[tt] performs the binning in relative instead of absolute",
        "box coordinates, and scales the final output with the average box dimension",
        "along the output axis. This can be used in combination with [TT]-center[tt].[PAR]",

        "Densities are in kg/m^3, and number densities or electron densities can also be",
        "calculated. For electron densities, a file describing the number of",
        "electrons for each type of atom should be provided using [TT]-ei[tt].",
        "It should look like::",
        "",
        "   2",
        "   atomname = nrelectrons",
        "   atomname = nrelectrons",
        "",
        "The first line contains the number of lines to read from the file.",
        "There should be one line for each unique atom name in your system.",
        "The number of electrons for each atom is modified by its atomic",
        "partial charge.[PAR]",

        "IMPORTANT CONSIDERATIONS FOR BILAYERS[PAR]",
        "One of the most common usage scenarios is to calculate the density of various",
        "groups across a lipid bilayer, typically with the z axis being the normal",
        "direction. For short simulations, small systems, and fixed box sizes this",
        "will work fine, but for the more general case lipid bilayers can be complicated.",
        "The first problem that while both proteins and lipids have low volume",
        "compressibility, lipids have quite high area compressiblity. This means the",
        "shape of the box (thickness and area/lipid) will fluctuate substantially even",
        "for a fully relaxed system. Since GROMACS places the box between the origin",
        "and positive coordinates, this in turn means that a bilayer centered in the",
        "box will move a bit up/down due to these fluctuations, and smear out your",
        "profile. The easiest way to fix this (if you want pressure coupling) is",
        "to use the [TT]-center[tt] option that calculates the density profile with",
        "respect to the center of the box. Note that you can still center on the",
        "bilayer part even if you have a complex non-symmetric system with a bilayer",
        "and, say, membrane proteins - then our output will simply have more values",
        "on one side of the (center) origin reference.[PAR]",

        "Even the centered calculation will lead to some smearing out the output",
        "profiles, as lipids themselves are compressed and expanded. In most cases",
        "you probably want this (since it corresponds to macroscopic experiments),",
        "but if you want to look at molecular details you can use the [TT]-relative[tt]",
        "option to attempt to remove even more of the effects of volume fluctuations.[PAR]",

        "Finally, large bilayers that are not subject to a surface tension will exhibit",
        "undulatory fluctuations, where there are 'waves' forming in the system.",
        "This is a fundamental property of the biological system, and if you are",
        "comparing against experiments you likely want to include the undulation",
        "smearing effect.",
        "",
    };

    output_env_t       oenv;
    static const char *dens_opt[] =
    { NULL, "mass", "number", "charge", "electron", NULL };
    static int         axis        = 2;  /* normal to memb. default z  */
    static const char *axtitle     = "Z";
    static int         nslices     = 50; /* nr of slices defined       */
    static int         ngrps       = 1;  /* nr. of groups              */
    static gmx_bool    bSymmetrize = FALSE;
    static gmx_bool    bCenter     = FALSE;
    static gmx_bool    bRelative   = FALSE;

    t_pargs            pa[]        = {
        { "-d", FALSE, etSTR, {&axtitle},
          "Take the normal on the membrane in direction X, Y or Z." },
        { "-sl",  FALSE, etINT, {&nslices},
          "Divide the box in this number of slices." },
        { "-dens",    FALSE, etENUM, {dens_opt},
          "Density"},
        { "-ng",       FALSE, etINT, {&ngrps},
          "Number of groups of which to compute densities." },
        { "-center",   FALSE, etBOOL, {&bCenter},
          "Perform the binning relative to the center of the (changing) box. Useful for bilayers." },
        { "-symm",     FALSE, etBOOL, {&bSymmetrize},
          "Symmetrize the density along the axis, with respect to the center. Useful for bilayers." },
        { "-relative", FALSE, etBOOL, {&bRelative},
          "Use relative coordinates for changing boxes and scale output by average dimensions." }
    };

    const char        *bugs[] = {
        "When calculating electron densities, atomnames are used instead of types. This is bad.",
    };

    double           **density;        /* density per slice          */
    real               slWidth;        /* width of one slice         */
    char              *grpname_center; /* centering group name     */
    char             **grpname;        /* groupnames                 */
    int                nr_electrons;   /* nr. electrons              */
    int                ncenter;        /* size of centering group    */
    int               *ngx;            /* sizes of groups            */
    t_electron        *el_tab;         /* tabel with nr. of electrons*/
    t_topology        *top;            /* topology               */
    int                ePBC;
    atom_id           *index_center;   /* index for centering group  */
    atom_id          **index;          /* indices for all groups     */
    int                i;

    t_filenm           fnm[] = { /* files for g_density       */
        { efTRX, "-f", NULL,  ffREAD },
        { efNDX, NULL, NULL,  ffOPTRD },
        { efTPR, NULL, NULL,  ffREAD },
        { efDAT, "-ei", "electrons", ffOPTRD }, /* file with nr. of electrons */
        { efXVG, "-o", "density", ffWRITE },
    };

#define NFILE asize(fnm)

    if (!parse_common_args(&argc, argv, PCA_CAN_VIEW | PCA_CAN_TIME,
                           NFILE, fnm, asize(pa), pa, asize(desc), desc, asize(bugs), bugs,
                           &oenv))
    {
        return 0;
    }

    if (bSymmetrize && !bCenter)
    {
        fprintf(stderr, "Can not symmetrize without centering. Turning on -center\n");
        bCenter = TRUE;
    }
    /* Calculate axis */
    axis = toupper(axtitle[0]) - 'X';

    top = read_top(ftp2fn(efTPR, NFILE, fnm), &ePBC); /* read topology file */
    if (dens_opt[0][0] == 'n')
    {
        for (i = 0; (i < top->atoms.nr); i++)
        {
            top->atoms.atom[i].m = 1;
        }
    }
    else if (dens_opt[0][0] == 'c')
    {
        for (i = 0; (i < top->atoms.nr); i++)
        {
            top->atoms.atom[i].m = top->atoms.atom[i].q;
        }
    }

    snew(grpname, ngrps);
    snew(index, ngrps);
    snew(ngx, ngrps);

    if (bCenter)
    {
        fprintf(stderr,
                "\nNote: that the center of mass is calculated inside the box without applying\n"
                "any special periodicity. If necessary, it is your responsibility to first use\n"
                "trjconv to make sure atoms in this group are placed in the right periodicity.\n\n"
                "Select the group to center density profiles around:\n");
        get_index(&top->atoms, ftp2fn_null(efNDX, NFILE, fnm), 1, &ncenter,
                  &index_center, &grpname_center);
    }
    else
    {
        ncenter = 0;
    }

    fprintf(stderr, "\nSelect %d group%s to calculate density for:\n", ngrps, (ngrps > 1) ? "s" : "");
    get_index(&top->atoms, ftp2fn_null(efNDX, NFILE, fnm), ngrps, ngx, index, grpname);

    if (dens_opt[0][0] == 'e')
    {
        nr_electrons =  get_electrons(&el_tab, ftp2fn(efDAT, NFILE, fnm));
        fprintf(stderr, "Read %d atomtypes from datafile\n", nr_electrons);

        calc_electron_density(ftp2fn(efTRX, NFILE, fnm), index, ngx, &density,
                              &nslices, top, ePBC, axis, ngrps, &slWidth, el_tab,
                              nr_electrons, bCenter, index_center, ncenter,
                              bRelative, oenv);
    }
    else
    {
        calc_density(ftp2fn(efTRX, NFILE, fnm), index, ngx, &density, &nslices, top,
                     ePBC, axis, ngrps, &slWidth, bCenter, index_center, ncenter,
                     bRelative, oenv);
    }

    plot_density(density, opt2fn("-o", NFILE, fnm),
                 nslices, ngrps, grpname, slWidth, dens_opt,
                 bCenter, bRelative, bSymmetrize, oenv);

    do_view(oenv, opt2fn("-o", NFILE, fnm), "-nxy");  /* view xvgr file */
    return 0;
}
