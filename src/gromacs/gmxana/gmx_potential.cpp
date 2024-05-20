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

#include <vector>

#include "gromacs/commandline/pargs.h"
#include "gromacs/commandline/viewit.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/gmxana/princ.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pbcutil/rmpbc.h"
#include "gromacs/topology/index.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

#define EPS0 8.85419E-12
#define ELC 1.60219E-19

/****************************************************************************/
/* This program calculates the electrostatic potential across the box by    */
/* determining the charge density in slices of the box and integrating these*/
/* twice.                                                                   */
/* Peter Tieleman, April 1995                                               */
/* It now also calculates electrostatic potential in spherical micelles,    */
/* using \frac{1}{r}\frac{d^2r\Psi}{r^2} = - \frac{\rho}{\epsilon_0}        */
/* This probably sucks but it seems to work.                                */
/****************************************************************************/

/* this routine integrates the array data and returns the resulting array */
/* routine uses simple trapezoid rule                                     */
static void p_integrate(double* result, const double data[], int ndata, double slWidth, int cb, int ce)
{
    int    i, slice;
    double sum;

    if (ndata <= 2)
    {
        fprintf(stderr,
                "Warning: nr of slices very small. This will result"
                "in nonsense.\n");
    }

    fprintf(stderr, "Integrating from slice %d to slice %d\n", cb, ndata - ce);

    for (slice = cb; slice < (ndata - ce); slice++)
    {
        sum = 0;
        for (i = cb; i < slice; i++)
        {
            sum += slWidth * (data[i] + 0.5 * (data[i + 1] - data[i]));
        }
        result[slice] = sum;
    }
}

static void center_coords(const t_atoms* atoms, const int* index_center, int ncenter, matrix box, rvec x0[])
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

static void calc_potential(const char*             fn,
                           int**                   index,
                           int                     gnx[],
                           double***               slPotential,
                           double***               slCharge,
                           double***               slField,
                           int*                    nslices,
                           const t_topology*       top,
                           PbcType                 pbcType,
                           int                     axis,
                           int                     nr_grps,
                           double*                 slWidth,
                           double                  fudge_z,
                           gmx_bool                bSpherical,
                           gmx_bool                bCenter,
                           const int*              index_center,
                           int                     ncenter,
                           gmx_bool                bCorrect,
                           int                     cb,
                           int                     ce,
                           const gmx_output_env_t* oenv)
{
    rvec*        x0;     /* coordinates without pbc */
    matrix       box;    /* box (3x3) */
    int          natoms; /* nr. atoms in trj */
    t_trxstatus* status;
    int          i, n;
    int          nr_frames = 0;
    int          slice;
    double       qsum, nn;
    real         t;
    double       z;
    rvec         xcm;
    real         boxSize;
    real         sliceWidth;
    double       averageBoxSize;
    gmx_rmpbc_t  gpbc = nullptr;

    if ((natoms = read_first_x(oenv, &status, fn, &t, &x0, box)) == 0)
    {
        gmx_fatal(FARGS, "Could not read coordinates from statusfile\n");
    }

    averageBoxSize = 0;

    if (!*nslices)
    {
        *nslices = static_cast<int>(box[axis][axis] * 10.0); /* default value */
        fprintf(stderr, "\nDividing the box in %d slices\n", *nslices);
    }

    snew(*slField, nr_grps);
    snew(*slCharge, nr_grps);
    snew(*slPotential, nr_grps);

    for (i = 0; i < nr_grps; i++)
    {
        snew((*slField)[i], *nslices);
        snew((*slCharge)[i], *nslices);
        snew((*slPotential)[i], *nslices);
    }


    gpbc = gmx_rmpbc_init(&top->idef, pbcType, natoms);

    /*********** Start processing trajectory ***********/
    do
    {
        gmx_rmpbc_apply(gpbc, natoms, box, x0);

        // Translate atoms so the com of the center-group is in the
        // box geometrical center.
        if (bCenter)
        {
            center_coords(&top->atoms, index_center, ncenter, box, x0);
        }

        /* calculate position of center of mass based on group 1 */
        calc_xcm(x0, gnx[0], index[0], top->atoms.atom, xcm, FALSE);
        svmul(-1, xcm, xcm);

        boxSize    = box[axis][axis];
        sliceWidth = boxSize / *nslices;
        averageBoxSize += boxSize;

        for (n = 0; n < nr_grps; n++)
        {
            /* Check whether we actually have all positions of the requested index
             * group in the trajectory file */
            if (gnx[n] > natoms)
            {
                gmx_fatal(FARGS,
                          "You selected a group with %d atoms, but only %d atoms\n"
                          "were found in the trajectory.\n",
                          gnx[n],
                          natoms);
            }
            for (i = 0; i < gnx[n]; i++) /* loop over all atoms in index file */
            {
                if (bSpherical)
                {
                    rvec_add(x0[index[n][i]], xcm, x0[index[n][i]]);
                    /* only distance from origin counts, not sign */
                    slice = static_cast<int>(norm(x0[index[n][i]]) / sliceWidth);

                    /* this is a nice check for spherical groups but not for
                       all water in a cubic box since a lot will fall outside
                       the sphere
                       if (slice > (*nslices))
                       {
                       fprintf(stderr,"Warning: slice = %d\n",slice);
                       }
                     */
                    (*slCharge)[n][slice] += top->atoms.atom[index[n][i]].q;
                }
                else
                {
                    z = x0[index[n][i]][axis];
                    z = z + fudge_z;
                    if (z < 0)
                    {
                        z += boxSize;
                    }
                    if (z > boxSize)
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
                        slice = static_cast<int>(sliceIndexOffset + *nslices / 2.0);
                    }
                    else
                    {
                        slice = static_cast<int>(z / sliceWidth);
                    }
                    // Safeguard to avoid potential rounding errors during truncation
                    // Add nslices first (in case sliceIndex was negative),
                    // then clamp with modulo operation.
                    slice = (slice + *nslices) % *nslices;

                    (*slCharge)[n][slice] += top->atoms.atom[index[n][i]].q;
                }
            }
        }
        nr_frames++;
    } while (read_next_x(oenv, status, &t, x0, box));

    gmx_rmpbc_done(gpbc);

    /*********** done with status file **********/
    close_trx(status);

    /* slCharge now contains the total charge per slice, summed over all
       frames. Now divide by nr_frames and integrate twice
     */

    averageBoxSize /= nr_frames;
    *slWidth = averageBoxSize / (*nslices);

    if (bSpherical)
    {
        fprintf(stderr,
                "\n\nRead %d frames from trajectory. Calculating potential"
                "in spherical coordinates\n",
                nr_frames);
    }
    else
    {
        fprintf(stderr, "\n\nRead %d frames from trajectory. Calculating potential\n", nr_frames);
    }

    for (n = 0; n < nr_grps; n++)
    {
        for (i = 0; i < *nslices; i++)
        {
            if (bSpherical)
            {
                /* charge per volume is now the summed charge, divided by the nr
                   of frames and by the volume of the slice it's in, 4pi r^2 dr
                 */
                double sliceVolume = 4 * M_PI * gmx::square(i) * gmx::square(*slWidth) * *slWidth;
                if (sliceVolume == 0)
                {
                    (*slCharge)[n][i] = 0;
                }
                else
                {
                    (*slCharge)[n][i] /= (nr_frames * sliceVolume);
                }
            }
            else
            {
                double sliceVolume = (box[XX][XX] * box[YY][YY] * box[ZZ][ZZ]) / (*nslices);
                (*slCharge)[n][i] /= (nr_frames * sliceVolume);
            }
        }
        /* Now we have charge densities */
    }

    if (bCorrect && !bSpherical)
    {
        for (n = 0; n < nr_grps; n++)
        {
            nn   = 0;
            qsum = 0;
            for (i = 0; i < *nslices; i++)
            {
                if (std::abs((*slCharge)[n][i]) >= GMX_DOUBLE_MIN)
                {
                    nn++;
                    qsum += (*slCharge)[n][i];
                }
            }
            qsum /= nn;
            for (i = 0; i < *nslices; i++)
            {
                if (std::abs((*slCharge)[n][i]) >= GMX_DOUBLE_MIN)
                {
                    (*slCharge)[n][i] -= qsum;
                }
            }
        }
    }

    for (n = 0; n < nr_grps; n++)
    {
        /* integrate twice to get field and potential */
        p_integrate((*slField)[n], (*slCharge)[n], *nslices, *slWidth, cb, ce);
    }


    if (bCorrect && !bSpherical)
    {
        for (n = 0; n < nr_grps; n++)
        {
            nn   = 0;
            qsum = 0;
            for (i = 0; i < *nslices; i++)
            {
                if (std::abs((*slCharge)[n][i]) >= GMX_DOUBLE_MIN)
                {
                    nn++;
                    qsum += (*slField)[n][i];
                }
            }
            qsum /= nn;
            for (i = 0; i < *nslices; i++)
            {
                if (std::abs((*slCharge)[n][i]) >= GMX_DOUBLE_MIN)
                {
                    (*slField)[n][i] -= qsum;
                }
            }
        }
    }

    for (n = 0; n < nr_grps; n++)
    {
        p_integrate((*slPotential)[n], (*slField)[n], *nslices, *slWidth, cb, ce);
    }

    /* Now correct for eps0 and in spherical case for r*/
    for (n = 0; n < nr_grps; n++)
    {
        for (i = 0; i < *nslices; i++)
        {
            if (bSpherical)
            {
                (*slPotential)[n][i] = ELC * (*slPotential)[n][i] * -1.0E9 / (EPS0 * i * (*slWidth));
                (*slField)[n][i]     = ELC * (*slField)[n][i] * 1E18 / (EPS0 * i * (*slWidth));
            }
            else
            {
                (*slPotential)[n][i] = ELC * (*slPotential)[n][i] * -1.0E9 / EPS0;
                (*slField)[n][i]     = ELC * (*slField)[n][i] * 1E18 / EPS0;
            }
        }
    }

    sfree(x0); /* free memory used by coordinate array */
}

static void plot_potential(double*                          potential[],
                           double*                          charge[],
                           double*                          field[],
                           const char*                      afile,
                           const char*                      bfile,
                           const char*                      cfile,
                           int                              nslices,
                           int                              nr_grps,
                           gmx::ArrayRef<const std::string> grpname,
                           double                           slWidth,
                           gmx_bool                         bCenter,
                           gmx_bool                         bSymmetrize,
                           int                              cb,
                           int                              ce,
                           const gmx_output_env_t*          oenv)
{
    FILE *pot,    /* xvgr file with potential */
            *cha, /* xvgr file with charges   */
            *fie; /* xvgr files with fields   */
    int         slice, n;
    const char* title  = nullptr;
    const char* xlabel = nullptr;

    xlabel = bCenter ? "Average relative position from center (nm)" : "Average coordinate (nm)";

    title = bSymmetrize ? "Symmetrized electrostatic potential" : "Electrostatic Potential";
    pot   = xvgropen(afile, title, xlabel, "Potential (V)", oenv);
    xvgrLegend(pot, grpname, oenv);

    title = bSymmetrize ? "Symmetrized charge distribution" : "Charge Distribution";
    cha   = xvgropen(bfile, title, xlabel, "Charge density (q/nm\\S3\\N)", oenv);
    xvgrLegend(cha, grpname, oenv);

    title = bSymmetrize ? "Symmetrized electric field" : "Electric Field";
    fie   = xvgropen(cfile, title, xlabel, "Field (V/nm)", oenv);
    xvgrLegend(fie, grpname, oenv);

    for (slice = cb; slice < (nslices - ce); slice++)
    {
        float axisPosition;
        if (bCenter)
        {
            axisPosition = (slice - nslices / 2.0) * slWidth;
        }
        else
        {
            axisPosition = slice * slWidth;
        }

        fprintf(pot, "%20.16g  ", axisPosition);
        fprintf(cha, "%20.16g  ", axisPosition);
        fprintf(fie, "%20.16g  ", axisPosition);
        for (n = 0; n < nr_grps; n++)
        {
            float potentialValue;
            float fieldValue;
            float chargeValue;

            if (bSymmetrize)
            {
                potentialValue = (potential[n][slice] + potential[n][nslices - slice - 1]) / 2.0;
                fieldValue     = (field[n][slice] + field[n][nslices - slice - 1]) / 2.0;
                chargeValue    = (charge[n][slice] + charge[n][nslices - slice - 1]) / 2.0;
            }
            else
            {
                potentialValue = potential[n][slice];
                fieldValue     = field[n][slice];
                chargeValue    = charge[n][slice];
            }
            fprintf(pot, "   %20.16g", potentialValue);
            fprintf(fie, "   %20.16g", fieldValue / 1e9); // convert to V/nm
            fprintf(cha, "   %20.16g", chargeValue);
        }
        fprintf(pot, "\n");
        fprintf(cha, "\n");
        fprintf(fie, "\n");
    }

    xvgrclose(pot);
    xvgrclose(cha);
    xvgrclose(fie);
}

int gmx_potential(int argc, char* argv[])
{
    const char* desc[] = {
        "[THISMODULE] computes the electrostatical potential across the box. The potential is",
        "calculated by first summing the charges per slice and then integrating",
        "twice of this charge distribution. Periodic boundaries are not taken",
        "into account. Reference of potential is taken to be the left side of",
        "the box. It is also possible to calculate the potential in spherical",
        "coordinates as function of r by calculating a charge distribution in",
        "spherical slices and twice integrating them. epsilon_r is taken as 1,",
        "but 2 is more appropriate in many cases.",
        "",
        "Option [TT]-center[tt] performs the histogram binning and potential",
        "calculation relative to the center of an arbitrary group, in absolute box",
        "coordinates. If you are calculating profiles along the Z axis box dimension bZ,",
        "output would be from -bZ/2 to bZ/2 if you center based on the entire system.",

        "Option [TT]-symm[tt] symmetrizes the output around the center. This will",
        "automatically turn on [TT]-center[tt] too.",

    };
    gmx_output_env_t*  oenv;
    static int         axis        = 2; /* normal to memb. default z  */
    static const char* axtitle     = "Z";
    static int         nslices     = 10; /* nr of slices defined       */
    static int         ngrps       = 1;
    static gmx_bool    bSpherical  = FALSE; /* default is bilayer types   */
    static real        fudge_z     = 0;     /* translate coordinates      */
    static gmx_bool    bCorrect    = false;
    int                cb          = 0;
    int                ce          = 0;
    static gmx_bool    bSymmetrize = false;
    static gmx_bool    bCenter     = false;

    t_pargs pa[] = {
        { "-d",
          FALSE,
          etSTR,
          { &axtitle },
          "Take the normal on the membrane in direction X, Y or Z." },
        { "-sl",
          FALSE,
          etINT,
          { &nslices },
          "Calculate potential as function of boxlength, dividing the box"
          " in this number of slices." },
        { "-cb",
          FALSE,
          etINT,
          { &cb },
          "Discard this number of  first slices of box for integration" },
        { "-ce",
          FALSE,
          etINT,
          { &ce },
          "Discard this number of last slices of box for integration" },
        { "-tz",
          FALSE,
          etREAL,
          { &fudge_z },
          "Translate all coordinates by this distance in the direction of the box" },
        { "-spherical", FALSE, etBOOL, { &bSpherical }, "Calculate in spherical coordinates" },
        { "-ng", FALSE, etINT, { &ngrps }, "Number of groups to consider" },
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
          "bilayers." },
        { "-correct",
          FALSE,
          etBOOL,
          { &bCorrect },
          "Assume net zero charge of groups to improve accuracy" }
    };
    const char* bugs[] = { "Discarding slices for integration should not be necessary." };

    double **potential,         /* potential per slice        */
            **charge,           /* total charge per slice     */
            **field,            /* field per slice            */
            slWidth;            /* width of one slice         */
    char**      grpname;        /* groupnames                 */
    char*       grpname_center; /* centering group name     */
    int*        ngx;            /* sizes of groups            */
    t_topology* top;            /* topology        */
    PbcType     pbcType;
    int*        index_center; /* index for centering group  */
    int         ncenter;      /* size of centering group    */
    int**       index;        /* indices for all groups     */
    t_filenm    fnm[] = {
        /* files for gmx order       */
        { efTRX, "-f", nullptr, ffREAD },      /* trajectory file             */
        { efNDX, nullptr, nullptr, ffREAD },   /* index file          */
        { efTPR, nullptr, nullptr, ffREAD },   /* topology file               */
        { efXVG, "-o", "potential", ffWRITE }, /* xvgr output file    */
        { efXVG, "-oc", "charge", ffWRITE },   /* xvgr output file    */
        { efXVG, "-of", "field", ffWRITE },    /* xvgr output file    */
    };

#define NFILE asize(fnm)

    if (!parse_common_args(
                &argc, argv, PCA_CAN_VIEW | PCA_CAN_TIME, NFILE, fnm, asize(pa), pa, asize(desc), desc, asize(bugs), bugs, &oenv))
    {
        return 0;
    }

    if (bSpherical && (bCenter || bSymmetrize))
    {
        fprintf(stderr,
                "Centering/symmetrization not supported for spherical potential. Disabling.\n");
        bCenter     = false;
        bSymmetrize = false;
    }
    /* Calculate axis */
    axis = toupper(axtitle[0]) - 'X';

    top = read_top(ftp2fn(efTPR, NFILE, fnm), &pbcType); /* read topology file */

    snew(grpname, ngrps);
    snew(index, ngrps);
    snew(ngx, ngrps);

    rd_index(ftp2fn(efNDX, NFILE, fnm), ngrps, ngx, index, grpname);

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

    calc_potential(ftp2fn(efTRX, NFILE, fnm),
                   index,
                   ngx,
                   &potential,
                   &charge,
                   &field,
                   &nslices,
                   top,
                   pbcType,
                   axis,
                   ngrps,
                   &slWidth,
                   fudge_z,
                   bSpherical,
                   bCenter,
                   index_center,
                   ncenter,
                   bCorrect,
                   cb,
                   ce,
                   oenv);

    std::vector<std::string> names;
    names.resize(ngrps);
    for (int i = 0; i < ngrps; ++i)
    {
        names[i] = grpname[i];
    }
    plot_potential(potential,
                   charge,
                   field,
                   opt2fn("-o", NFILE, fnm),
                   opt2fn("-oc", NFILE, fnm),
                   opt2fn("-of", NFILE, fnm),
                   nslices,
                   ngrps,
                   names,
                   slWidth,
                   bCenter,
                   bSymmetrize,
                   cb,
                   ce,
                   oenv);

    do_view(oenv, opt2fn("-o", NFILE, fnm), nullptr);  /* view xvgr file */
    do_view(oenv, opt2fn("-oc", NFILE, fnm), nullptr); /* view xvgr file */
    do_view(oenv, opt2fn("-of", NFILE, fnm), nullptr); /* view xvgr file */

    return 0;
}
