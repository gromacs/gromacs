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
#include <ctype.h>

#include "sysstuff.h"
#include <string.h>
#include "string2.h"
#include "typedefs.h"
#include "smalloc.h"
#include "macros.h"
#include "gstat.h"
#include "vec.h"
#include "xvgr.h"
#include "pbc.h"
#include "copyrite.h"
#include "futil.h"
#include "statutil.h"
#include "index.h"
#include "tpxio.h"
#include "physics.h"
#include "gmx_ana.h"

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

    if (!(in = ffopen(fn, "r")))
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
        (*eltab)[i].atomname = strdup(tempname);
    }
    ffclose(in);

    /* sort the list */
    fprintf(stderr, "Sorting list..\n");
    qsort ((void*)*eltab, nr, sizeof(t_electron),
           (int(*)(const void*, const void*))compare);

    return nr;
}

void center_coords(t_atoms *atoms, matrix box, rvec x0[], int axis)
{
    int  i, m;
    real tmass, mm;
    rvec com, shift, box_center;

    tmass = 0;
    clear_rvec(com);
    for (i = 0; (i < atoms->nr); i++)
    {
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
    rvec_sub(box_center, com, shift);
    shift[axis] -= box_center[axis];

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
                           const output_env_t oenv)
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

    if (!*nslices)
    {
        *nslices = (int)(box[axis][axis] * 10); /* default value */
    }
    fprintf(stderr, "\nDividing the box in %d slices\n", *nslices);

    snew(*slDensity, nr_grps);
    for (i = 0; i < nr_grps; i++)
    {
        snew((*slDensity)[i], *nslices);
    }

    gpbc = gmx_rmpbc_init(&top->idef, ePBC, top->atoms.nr, box);
    /*********** Start processing trajectory ***********/
    do
    {
        gmx_rmpbc(gpbc, natoms, box, x0);

        if (bCenter)
        {
            center_coords(&top->atoms, box, x0, axis);
        }

        *slWidth = box[axis][axis]/(*nslices);
        invvol   = *nslices/(box[XX][XX]*box[YY][YY]*box[ZZ][ZZ]);

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

                /* determine which slice atom is in */
                slice           = (z / (*slWidth));
                sought.nr_el    = 0;
                sought.atomname = strdup(*(top->atoms.atomname[index[n][i]]));

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
    while (read_next_x(oenv, status, &t, natoms, x0, box));
    gmx_rmpbc_done(gpbc);

    /*********** done with status file **********/
    close_trj(status);

/* slDensity now contains the total number of electrons per slice, summed
   over all frames. Now divide by nr_frames and volume of slice
 */

    fprintf(stderr, "\nRead %d frames from trajectory. Counting electrons\n",
            nr_frames);

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
                  const output_env_t oenv)
{
    rvec        *x0;      /* coordinates without pbc */
    matrix       box;     /* box (3x3) */
    double       invvol;
    int          natoms;  /* nr. atoms in trj */
    t_trxstatus *status;
    int        **slCount, /* nr. of atoms in one slice for a group */
                 i, j, n, /* loop indices */
                 teller    = 0,
                 ax1       = 0, ax2 = 0,
                 nr_frames = 0, /* number of frames */
                 slice;         /* current slice */
    real         t,
                 z;
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

    gpbc = gmx_rmpbc_init(&top->idef, ePBC, top->atoms.nr, box);
    /*********** Start processing trajectory ***********/
    do
    {
        gmx_rmpbc(gpbc, natoms, box, x0);

        if (bCenter)
        {
            center_coords(&top->atoms, box, x0, axis);
        }

        *slWidth = box[axis][axis]/(*nslices);
        invvol   = *nslices/(box[XX][XX]*box[YY][YY]*box[ZZ][ZZ]);
        teller++;

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

                /* determine which slice atom is in */
                slice                   = (int)(z / (*slWidth));
                (*slDensity)[n][slice] += top->atoms.atom[index[n][i]].m*invvol;
            }
        }
        nr_frames++;
    }
    while (read_next_x(oenv, status, &t, natoms, x0, box));
    gmx_rmpbc_done(gpbc);

    /*********** done with status file **********/
    close_trj(status);

    /* slDensity now contains the total mass per slice, summed over all
       frames. Now divide by nr_frames and volume of slice
     */

    fprintf(stderr, "\nRead %d frames from trajectory. Calculating density\n",
            nr_frames);

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
                  gmx_bool bSymmetrize, const output_env_t oenv)
{
    FILE       *den;
    const char *ylabel = NULL;
    int         slice, n;
    real        ddd;

    switch (dens_opt[0][0])
    {
        case 'm': ylabel = "Density (kg m\\S-3\\N)"; break;
        case 'n': ylabel = "Number density (nm\\S-3\\N)"; break;
        case 'c': ylabel = "Charge density (e nm\\S-3\\N)"; break;
        case 'e': ylabel = "Electron density (e nm\\S-3\\N)"; break;
    }

    den = xvgropen(afile, "Partial densities", "Box (nm)", ylabel, oenv);

    xvgr_legend(den, nr_grps, (const char**)grpname, oenv);

    for (slice = 0; (slice < nslices); slice++)
    {
        fprintf(den, "%12g  ", slice * slWidth);
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

    ffclose(den);
}

int gmx_density(int argc, char *argv[])
{
    const char        *desc[] = {
        "Compute partial densities across the box, using an index file.[PAR]",
        "For the total density of NPT simulations, use [TT]g_energy[tt] instead.",
        "[PAR]",
        "Densities are in kg/m^3, and number densities or electron densities can also be",
        "calculated. For electron densities, a file describing the number of",
        "electrons for each type of atom should be provided using [TT]-ei[tt].",
        "It should look like:[BR]",
        "   [TT]2[tt][BR]",
        "   [TT]atomname = nrelectrons[tt][BR]",
        "   [TT]atomname = nrelectrons[tt][BR]",
        "The first line contains the number of lines to read from the file.",
        "There should be one line for each unique atom name in your system.",
        "The number of electrons for each atom is modified by its atomic",
        "partial charge."
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
    t_pargs            pa[]        = {
        { "-d", FALSE, etSTR, {&axtitle},
          "Take the normal on the membrane in direction X, Y or Z." },
        { "-sl",  FALSE, etINT, {&nslices},
          "Divide the box in this number of slices." },
        { "-dens",    FALSE, etENUM, {dens_opt},
          "Density"},
        { "-ng",       FALSE, etINT, {&ngrps},
          "Number of groups of which to compute densities." },
        { "-symm",    FALSE, etBOOL, {&bSymmetrize},
          "Symmetrize the density along the axis, with respect to the center. Useful for bilayers." },
        { "-center",  FALSE, etBOOL, {&bCenter},
          "Shift the center of mass along the axis to zero. This means if your axis is Z and your box is bX, bY, bZ, the center of mass will be at bX/2, bY/2, 0."}
    };

    const char        *bugs[] = {
        "When calculating electron densities, atomnames are used instead of types. This is bad.",
    };

    double           **density;      /* density per slice          */
    real               slWidth;      /* width of one slice         */
    char             **grpname;      /* groupnames                 */
    int                nr_electrons; /* nr. electrons              */
    int               *ngx;          /* sizes of groups            */
    t_electron        *el_tab;       /* tabel with nr. of electrons*/
    t_topology        *top;          /* topology               */
    int                ePBC;
    atom_id          **index;        /* indices for all groups     */
    int                i;

    t_filenm           fnm[] = { /* files for g_density       */
        { efTRX, "-f", NULL,  ffREAD },
        { efNDX, NULL, NULL,  ffOPTRD },
        { efTPX, NULL, NULL,  ffREAD },
        { efDAT, "-ei", "electrons", ffOPTRD }, /* file with nr. of electrons */
        { efXVG, "-o", "density", ffWRITE },
    };

#define NFILE asize(fnm)

    CopyRight(stderr, argv[0]);

    parse_common_args(&argc, argv, PCA_CAN_VIEW | PCA_CAN_TIME | PCA_BE_NICE,
                      NFILE, fnm, asize(pa), pa, asize(desc), desc, asize(bugs), bugs,
                      &oenv);

    if (bSymmetrize && !bCenter)
    {
        fprintf(stderr, "Can not symmetrize without centering. Turning on -center\n");
        bCenter = TRUE;
    }
    /* Calculate axis */
    axis = toupper(axtitle[0]) - 'X';

    top = read_top(ftp2fn(efTPX, NFILE, fnm), &ePBC); /* read topology file */
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

    get_index(&top->atoms, ftp2fn_null(efNDX, NFILE, fnm), ngrps, ngx, index, grpname);

    if (dens_opt[0][0] == 'e')
    {
        nr_electrons =  get_electrons(&el_tab, ftp2fn(efDAT, NFILE, fnm));
        fprintf(stderr, "Read %d atomtypes from datafile\n", nr_electrons);

        calc_electron_density(ftp2fn(efTRX, NFILE, fnm), index, ngx, &density,
                              &nslices, top, ePBC, axis, ngrps, &slWidth, el_tab,
                              nr_electrons, bCenter, oenv);
    }
    else
    {
        calc_density(ftp2fn(efTRX, NFILE, fnm), index, ngx, &density, &nslices, top,
                     ePBC, axis, ngrps, &slWidth, bCenter, oenv);
    }

    plot_density(density, opt2fn("-o", NFILE, fnm),
                 nslices, ngrps, grpname, slWidth, dens_opt,
                 bSymmetrize, oenv);

    do_view(oenv, opt2fn("-o", NFILE, fnm), "-nxy");  /* view xvgr file */
    thanx(stderr);
    return 0;
}
