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

#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <filesystem>
#include <string>
#include <vector>

#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/commandline/viewit.h"
#include "gromacs/fileio/filetypes.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pbcutil/rmpbc.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/index.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"

enum class PbcType : int;
struct gmx_output_env_t;

typedef struct
{
    char* atomname;
    int   nr_el;
} t_electron;

/****************************************************************************/
/* This program calculates the partial density across the box.              */
/* Peter Tieleman, Mei 1995                                                 */
/****************************************************************************/

/* used for sorting the list */
static int compare(const void* a, const void* b)
{
    const t_electron *tmp1, *tmp2;
    tmp1 = static_cast<const t_electron*>(a);
    tmp2 = static_cast<const t_electron*>(b);

    return std::strcmp(tmp1->atomname, tmp2->atomname);
}

static int get_electrons(t_electron** eltab, const char* fn)
{
    char buffer[256];  /* to read in a line   */
    char tempname[80]; /* buffer to hold name */
    int  tempnr;

    FILE* in;
    int   nr; /* number of atomstypes to read */
    int   i;

    if ((in = gmx_ffopen(fn, "r")) == nullptr)
    {
        gmx_fatal(FARGS, "Couldn't open %s. Exiting.\n", fn);
    }

    if (nullptr == fgets(buffer, 255, in))
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
        if (fgets(buffer, 255, in) == nullptr)
        {
            gmx_fatal(FARGS, "reading datafile. Check your datafile.\n");
        }
        if (sscanf(buffer, "%s = %d", tempname, &tempnr) != 2)
        {
            gmx_fatal(FARGS, "Invalid line in datafile at line %d\n", i + 1);
        }
        (*eltab)[i].nr_el    = tempnr;
        (*eltab)[i].atomname = gmx_strdup(tempname);
    }
    gmx_ffclose(in);

    /* sort the list */
    fprintf(stderr, "Sorting list..\n");
    qsort(*eltab, nr, sizeof(t_electron), compare);

    return nr;
}

static void center_coords(t_atoms* atoms, const int* index_center, int ncenter, matrix box, rvec x0[])
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
            gmx_fatal(FARGS,
                      "Index %d refers to atom %d, which is larger than natoms (%d).",
                      k + 1,
                      i + 1,
                      atoms->nr);
        }
        mm = atoms->atom[i].m;
        tmass += mm;
        for (m = 0; (m < DIM); m++)
        {
            com[m] += mm * x0[i][m];
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

static void calc_electron_density(const char*             fn,
                                  int**                   index,
                                  const int               gnx[],
                                  double***               slDensity,
                                  int*                    nslices,
                                  t_topology*             top,
                                  PbcType                 pbcType,
                                  int                     axis,
                                  int                     nr_grps,
                                  real*                   slWidth,
                                  t_electron              eltab[],
                                  int                     nr,
                                  gmx_bool                bCenter,
                                  int*                    index_center,
                                  int                     ncenter,
                                  const gmx_output_env_t* oenv)
{
    rvec*        x0;  /* coordinates without pbc */
    matrix       box; /* box (3x3) */
    double       invvol;
    int          natoms; /* nr. atoms in trj */
    t_trxstatus* status;
    int          i, n;
    int          nr_frames = 0;
    t_electron*  found;  /* found by bsearch */
    t_electron   sought; /* thingie thought by bsearch */
    int          sliceIndex;
    real         boxSize;
    real         sliceWidth;
    double       averageBoxSize;
    gmx_rmpbc_t  gpbc = nullptr;

    real t, z;

    if (axis < 0 || axis >= DIM)
    {
        gmx_fatal(FARGS, "Invalid axes. Terminating\n");
    }

    if ((natoms = read_first_x(oenv, &status, fn, &t, &x0, box)) == 0)
    {
        gmx_fatal(FARGS, "Could not read coordinates from statusfile\n");
    }

    averageBoxSize = 0;

    if (!*nslices)
    {
        *nslices = static_cast<int>(box[axis][axis] * 10); /* default value */
        fprintf(stderr, "\nDividing the box in %d slices\n", *nslices);
    }

    snew(*slDensity, nr_grps);
    for (i = 0; i < nr_grps; i++)
    {
        snew((*slDensity)[i], *nslices);
    }

    gpbc = gmx_rmpbc_init(&top->idef, pbcType, top->atoms.nr);
    /*********** Start processing trajectory ***********/
    do
    {
        gmx_rmpbc_apply(gpbc, natoms, box, x0);

        /* Translate atoms so the com of the center-group is in the
         * box geometrical center.
         */
        if (bCenter)
        {
            center_coords(&top->atoms, index_center, ncenter, box, x0);
        }

        invvol = *nslices / (box[XX][XX] * box[YY][YY] * box[ZZ][ZZ]);

        boxSize    = box[axis][axis];
        sliceWidth = boxSize / *nslices;
        averageBoxSize += boxSize;

        for (n = 0; n < nr_grps; n++)
        {
            for (i = 0; i < gnx[n]; i++) /* loop over all atoms in index file */
            {
                z = x0[index[n][i]][axis];
                while (z < 0)
                {
                    z += boxSize;
                }
                while (z > boxSize)
                {
                    z -= boxSize;
                }

                /* determine which slice atom is in */
                if (bCenter)
                {
                    const real positionRelativeToCenter = z - boxSize / 2.0;
                    // Always round down since relative position might be negative.
                    const real sliceIndexOffset = std::floor(positionRelativeToCenter / sliceWidth);
                    // We kept sliceIndexOffset as floating-point in case nslices was odd
                    sliceIndex = static_cast<int>(sliceIndexOffset + *nslices / 2.0);
                }
                else
                {
                    sliceIndex = static_cast<int>(z / sliceWidth);
                }
                // Safeguard to avoid potential rounding errors during truncation
                // Add nslices first (in case sliceIndex was negative), then clamp with modulo operation.
                sliceIndex = (sliceIndex + *nslices) % *nslices;

                sought.nr_el    = 0;
                sought.atomname = gmx_strdup(*(top->atoms.atomname[index[n][i]]));

                /* now find the number of electrons. This is not efficient. */
                found = static_cast<t_electron*>(bsearch(&sought, eltab, nr, sizeof(t_electron), compare));

                if (found == nullptr)
                {
                    fprintf(stderr,
                            "Couldn't find %s. Add it to the .dat file\n",
                            *(top->atoms.atomname[index[n][i]]));
                }
                else
                {
                    (*slDensity)[n][sliceIndex] += (found->nr_el - top->atoms.atom[index[n][i]].q) * invvol;
                }
                free(sought.atomname);
            }
        }
        nr_frames++;
    } while (read_next_x(oenv, status, &t, x0, box));
    gmx_rmpbc_done(gpbc);

    /*********** done with status file **********/
    close_trx(status);

    /* slDensity now contains the total number of electrons per slice, summed
       over all frames. Now divide by nr_frames and volume of slice
     */

    fprintf(stderr, "\nRead %d frames from trajectory. Counting electrons\n", nr_frames);

    averageBoxSize /= nr_frames;
    *slWidth = averageBoxSize / (*nslices);

    for (n = 0; n < nr_grps; n++)
    {
        for (i = 0; i < *nslices; i++)
        {
            (*slDensity)[n][i] /= nr_frames;
        }
    }

    sfree(x0); /* free memory used by coordinate array */
}

static void calc_density(const char*             fn,
                         int**                   index,
                         const int               gnx[],
                         double***               slDensity,
                         int*                    nslices,
                         t_topology*             top,
                         PbcType                 pbcType,
                         int                     axis,
                         int                     nr_grps,
                         real*                   slWidth,
                         gmx_bool                bCenter,
                         int*                    index_center,
                         int                     ncenter,
                         const gmx_output_env_t* oenv,
                         const char**            dens_opt)
{
    rvec*        x0;  /* coordinates without pbc */
    matrix       box; /* box (3x3) */
    double       invvol;
    int          natoms; /* nr. atoms in trj */
    t_trxstatus* status;
    int          i, n;
    int          nr_frames = 0;
    real         t, z;
    real*        den_val; /* values from which the density is calculated */
    int          sliceIndex;
    real         boxSize;
    real         sliceWidth;
    double       averageBoxSize;
    gmx_rmpbc_t  gpbc = nullptr;

    if (axis < 0 || axis >= DIM)
    {
        gmx_fatal(FARGS, "Invalid axes. Terminating\n");
    }

    if ((natoms = read_first_x(oenv, &status, fn, &t, &x0, box)) == 0)
    {
        gmx_fatal(FARGS, "Could not read coordinates from statusfile\n");
    }

    averageBoxSize = 0;

    if (!*nslices)
    {
        *nslices = static_cast<int>(box[axis][axis] * 10); /* default value */
        fprintf(stderr, "\nDividing the box in %d slices\n", *nslices);
    }

    snew(*slDensity, nr_grps);
    for (i = 0; i < nr_grps; i++)
    {
        snew((*slDensity)[i], *nslices);
    }

    gpbc = gmx_rmpbc_init(&top->idef, pbcType, top->atoms.nr);
    /*********** Start processing trajectory ***********/

    snew(den_val, top->atoms.nr);
    if (dens_opt[0][0] == 'n')
    {
        for (i = 0; (i < top->atoms.nr); i++)
        {
            den_val[i] = 1;
        }
    }
    else if (dens_opt[0][0] == 'c')
    {
        for (i = 0; (i < top->atoms.nr); i++)
        {
            den_val[i] = top->atoms.atom[i].q;
        }
    }
    else
    {
        for (i = 0; (i < top->atoms.nr); i++)
        {
            den_val[i] = top->atoms.atom[i].m;
        }
    }

    do
    {
        gmx_rmpbc_apply(gpbc, natoms, box, x0);

        /* Translate atoms so the com of the center-group is in the
         * box geometrical center.
         */
        if (bCenter)
        {
            center_coords(&top->atoms, index_center, ncenter, box, x0);
        }

        invvol = *nslices / (box[XX][XX] * box[YY][YY] * box[ZZ][ZZ]);

        boxSize    = box[axis][axis];
        sliceWidth = boxSize / *nslices;
        averageBoxSize += boxSize;

        for (n = 0; n < nr_grps; n++)
        {
            for (i = 0; i < gnx[n]; i++) /* loop over all atoms in index file */
            {
                z = x0[index[n][i]][axis];
                while (z < 0)
                {
                    z += boxSize;
                }
                while (z > boxSize)
                {
                    z -= boxSize;
                }

                /* determine which slice atom is in */
                if (bCenter)
                {
                    const real positionRelativeToCenter = z - boxSize / 2.0;
                    // Always round down since relative position might be negative.
                    const real sliceIndexOffset = std::floor(positionRelativeToCenter / sliceWidth);
                    // We kept sliceIndexOffset as floating-point in case nslices was odd
                    sliceIndex = static_cast<int>(sliceIndexOffset + *nslices / 2.0);
                }
                else
                {
                    sliceIndex = static_cast<int>(z / sliceWidth);
                }
                // Safeguard to avoid potential rounding errors during truncation
                // Add nslices first (in case sliceIndex was negative), then clamp with modulo operation.
                sliceIndex = (sliceIndex + *nslices) % *nslices;

                (*slDensity)[n][sliceIndex] += den_val[index[n][i]] * invvol;
            }
        }
        nr_frames++;
    } while (read_next_x(oenv, status, &t, x0, box));
    gmx_rmpbc_done(gpbc);

    /*********** done with status file **********/
    close_trx(status);

    /* slDensity now contains the total mass per slice, summed over all
       frames. Now divide by nr_frames and volume of slice
     */

    fprintf(stderr, "\nRead %d frames from trajectory. Calculating density\n", nr_frames);

    averageBoxSize /= nr_frames;
    *slWidth = averageBoxSize / (*nslices);

    for (n = 0; n < nr_grps; n++)
    {
        for (i = 0; i < *nslices; i++)
        {
            (*slDensity)[n][i] /= nr_frames;
        }
    }

    sfree(x0); /* free memory used by coordinate array */
    sfree(den_val);
}

static void plot_density(double*                          slDensity[],
                         const char*                      afile,
                         int                              nslices,
                         gmx::ArrayRef<const std::string> grpname,
                         real                             slWidth,
                         const char**                     dens_opt,
                         gmx_bool                         bCenter,
                         gmx_bool                         bSymmetrize,
                         const gmx_output_env_t*          oenv)
{
    FILE*       den;
    const char* title  = nullptr;
    const char* xlabel = nullptr;
    const char* ylabel = nullptr;
    int         slice;
    real        ddd;
    real        axispos;

    title = bSymmetrize ? "Symmetrized partial density" : "Partial density";

    xlabel = bCenter ? "Average relative position from center (nm)" : "Average coordinate (nm)";

    switch (dens_opt[0][0])
    {
        case 'm': ylabel = "Density (kg m\\S-3\\N)"; break;
        case 'n': ylabel = "Number density (nm\\S-3\\N)"; break;
        case 'c': ylabel = "Charge density (e nm\\S-3\\N)"; break;
        case 'e': ylabel = "Electron density (e nm\\S-3\\N)"; break;
    }

    den = xvgropen(afile, title, xlabel, ylabel, oenv);

    xvgrLegend(den, grpname, oenv);

    for (slice = 0; (slice < nslices); slice++)
    {
        if (bCenter)
        {
            axispos = (slice - nslices / 2.0 + 0.5) * slWidth;
        }
        else
        {
            axispos = (slice + 0.5) * slWidth;
        }
        fprintf(den, "%12g  ", axispos);
        for (int n = 0; (n < gmx::ssize(grpname)); n++)
        {
            if (bSymmetrize)
            {
                ddd = (slDensity[n][slice] + slDensity[n][nslices - slice - 1]) * 0.5;
            }
            else
            {
                ddd = slDensity[n][slice];
            }
            if (dens_opt[0][0] == 'm')
            {
                fprintf(den, "   %12g", ddd * gmx::c_amu / (gmx::c_nano * gmx::c_nano * gmx::c_nano));
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

int gmx_density(int argc, char* argv[])
{
    const char* desc[] = {
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

        "The binning is now always performed in relative coordinates to account",
        "for changing box dimensions with pressure coupling, with the output",
        "scaled to the average box dimension along the output axis.[PAR]",

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

        "Finally, large bilayers that are not subject to a surface tension will exhibit",
        "undulatory fluctuations, where there are 'waves' forming in the system.",
        "This is a fundamental property of the biological system, and if you are",
        "comparing against experiments you likely want to include the undulation",
        "smearing effect.",
        "",
    };

    gmx_output_env_t*  oenv;
    static const char* dens_opt[]  = { nullptr, "mass", "number", "charge", "electron", nullptr };
    static int         axis        = 2; /* normal to memb. default z  */
    static const char* axtitle     = "Z";
    static int         nslices     = 50; /* nr of slices defined       */
    static int         ngrps       = 1;  /* nr. of groups              */
    static gmx_bool    bSymmetrize = FALSE;
    static gmx_bool    bCenter     = FALSE;

    t_pargs pa[] = {
        { "-d",
          FALSE,
          etSTR,
          { &axtitle },
          "Take the normal on the membrane in direction X, Y or Z." },
        { "-sl", FALSE, etINT, { &nslices }, "Divide the box in this number of slices." },
        { "-dens", FALSE, etENUM, { dens_opt }, "Density" },
        { "-ng", FALSE, etINT, { &ngrps }, "Number of groups of which to compute densities." },
        { "-center",
          FALSE,
          etBOOL,
          { &bCenter },
          "Perform the binning relative to the center of the (changing) box. Useful for "
          "bilayers." },
        { "-symm",
          FALSE,
          etBOOL,
          { &bSymmetrize },
          "Symmetrize the density along the axis, with respect to the center. Useful for "
          "bilayers." }
    };

    const char* bugs[] = {
        "When calculating electron densities, atomnames are used instead of types. This is bad.",
    };

    double**    density;        /* density per slice          */
    real        slWidth;        /* width of one slice         */
    char*       grpname_center; /* centering group name     */
    char**      grpname;        /* groupnames                 */
    int         nr_electrons;   /* nr. electrons              */
    int         ncenter;        /* size of centering group    */
    int*        ngx;            /* sizes of groups            */
    t_electron* el_tab;         /* tabel with nr. of electrons*/
    t_topology* top;            /* topology               */
    PbcType     pbcType;
    int*        index_center; /* index for centering group  */
    int**       index;        /* indices for all groups     */

    t_filenm fnm[] = {
        /* files for gmx density       */
        { efTRX, "-f", nullptr, ffREAD },
        { efNDX, nullptr, nullptr, ffOPTRD },
        { efTPR, nullptr, nullptr, ffREAD },
        { efDAT, "-ei", "electrons", ffOPTRD }, /* file with nr. of electrons */
        { efXVG, "-o", "density", ffWRITE },
    };

#define NFILE asize(fnm)

    if (!parse_common_args(
                &argc, argv, PCA_CAN_VIEW | PCA_CAN_TIME, NFILE, fnm, asize(pa), pa, asize(desc), desc, asize(bugs), bugs, &oenv))
    {
        return 0;
    }

    GMX_RELEASE_ASSERT(dens_opt[0] != nullptr, "Option setting inconsistency; dens_opt[0] is NULL");

    if (bSymmetrize && !bCenter)
    {
        fprintf(stderr, "Can not symmetrize without centering. Turning on -center\n");
        bCenter = TRUE;
    }
    /* Calculate axis */
    axis = toupper(axtitle[0]) - 'X';

    top = read_top(ftp2fn(efTPR, NFILE, fnm), &pbcType); /* read topology file */

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
        get_index(&top->atoms, ftp2fn_null(efNDX, NFILE, fnm), 1, &ncenter, &index_center, &grpname_center);
    }
    else
    {
        ncenter      = 0;
        index_center = nullptr;
    }

    fprintf(stderr, "\nSelect %d group%s to calculate density for:\n", ngrps, (ngrps > 1) ? "s" : "");
    get_index(&top->atoms, ftp2fn_null(efNDX, NFILE, fnm), ngrps, ngx, index, grpname);

    if (dens_opt[0][0] == 'e')
    {
        nr_electrons = get_electrons(&el_tab, ftp2fn(efDAT, NFILE, fnm));
        fprintf(stderr, "Read %d atomtypes from datafile\n", nr_electrons);

        calc_electron_density(ftp2fn(efTRX, NFILE, fnm),
                              index,
                              ngx,
                              &density,
                              &nslices,
                              top,
                              pbcType,
                              axis,
                              ngrps,
                              &slWidth,
                              el_tab,
                              nr_electrons,
                              bCenter,
                              index_center,
                              ncenter,
                              oenv);
    }
    else
    {
        calc_density(ftp2fn(efTRX, NFILE, fnm),
                     index,
                     ngx,
                     &density,
                     &nslices,
                     top,
                     pbcType,
                     axis,
                     ngrps,
                     &slWidth,
                     bCenter,
                     index_center,
                     ncenter,
                     oenv,
                     dens_opt);
    }

    std::vector<std::string> names;
    names.resize(ngrps);
    for (int i = 0; i < ngrps; ++i)
    {
        names[i] = grpname[i];
    }
    plot_density(density, opt2fn("-o", NFILE, fnm), nslices, names, slWidth, dens_opt, bCenter, bSymmetrize, oenv);

    do_view(oenv, opt2fn("-o", NFILE, fnm), "-nxy"); /* view xvgr file */
    return 0;
}
