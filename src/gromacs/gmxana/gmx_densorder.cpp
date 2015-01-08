/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2010,2011,2012,2013,2014,2015, by the GROMACS development team, led by
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
#include <string.h>

#include "gromacs/commandline/pargs.h"
#include "gromacs/correlationfunctions/autocorr.h"
#include "gromacs/correlationfunctions/expfit.h"
#include "gromacs/fileio/matio.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/gmxana/binsearch.h"
#include "gromacs/gmxana/dens_filter.h"
#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/gmxana/gstat.h"
#include "gromacs/gmxana/powerspect.h"
#include "gromacs/legacyheaders/copyrite.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pbcutil/rmpbc.h"
#include "gromacs/topology/index.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"

#ifdef GMX_DOUBLE
#define FLOOR(x) ((int) floor(x))
#else
#define FLOOR(x) ((int) floorf(x))
#endif

enum {
    methSEL, methBISECT, methFUNCFIT, methNR
};

static void center_coords(t_atoms *atoms, matrix box, rvec x0[], int axis)
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


static void density_in_time (const char *fn, atom_id **index, int gnx[], real bw, real bwz, int nsttblock, real *****Densdevel, int *xslices, int *yslices, int *zslices, int *tblock, t_topology *top, int ePBC, int axis, gmx_bool bCenter, gmx_bool bps1d, const output_env_t oenv)

{
/*
 * *****Densdevel pointer to array of density values in slices and frame-blocks Densdevel[*nsttblock][*xslices][*yslices][*zslices]
 * Densslice[x][y][z]
 * nsttblock - nr of frames in each time-block
 * bw  widths of normal slices
 *
 * axis	 - axis direction (normal to slices)
 * nndx - number ot atoms in **index
 * grpn	 - group number in index
 */
    t_trxstatus *status;
    gmx_rmpbc_t  gpbc = NULL;
    matrix       box;                    /* Box - 3x3 -each step*/
    rvec        *x0;                     /* List of Coord without PBC*/
    int          i, j,                   /* loop indices, checks etc*/
                 ax1     = 0, ax2 = 0,   /* tangent directions */
                 framenr = 0,            /* frame number in trajectory*/
                 slicex, slicey, slicez; /*slice # of x y z position */
    real ***Densslice = NULL;            /* Density-slice in one frame*/
    real    dscale;                      /*physical scaling factor*/
    real    t, x, y, z;                  /* time and coordinates*/
    rvec    bbww;

    *tblock = 0; /* blocknr in block average - initialise to 0*/
    /* Axis: X=0, Y=1,Z=2 */
    switch (axis)
    {
        case 0:
            ax1 = YY; ax2 = ZZ; /*Surface: YZ*/
            break;
        case 1:
            ax1 = ZZ; ax2 = XX; /* Surface : XZ*/
            break;
        case 2:
            ax1 = XX; ax2 = YY; /* Surface XY*/
            break;
        default:
            gmx_fatal(FARGS, "Invalid axes. Terminating\n");
    }

    if (read_first_x(oenv, &status, fn, &t, &x0, box) == 0)
    {
        gmx_fatal(FARGS, "Could not read coordinates from file"); /* Open trajectory for read*/


    }
    *zslices = 1+FLOOR(box[axis][axis]/bwz);
    *yslices = 1+FLOOR(box[ax2][ax2]/bw);
    *xslices = 1+FLOOR(box[ax1][ax1]/bw);
    if (bps1d)
    {
        if (*xslices < *yslices)
        {
            *xslices = 1;
        }
        else
        {
            *yslices = 1;
        }
    }
    fprintf(stderr,
            "\nDividing the box in %5d x %5d x %5d slices with binw %f along axis %d\n", *xslices, *yslices, *zslices, bw, axis );


    /****Start trajectory processing***/

    /*Initialize Densdevel and PBC-remove*/
    gpbc = gmx_rmpbc_init(&top->idef, ePBC, top->atoms.nr);

    *Densdevel = NULL;

    do
    {
        bbww[XX] = box[ax1][ax1]/ *xslices;
        bbww[YY] = box[ax2][ax2]/ *yslices;
        bbww[ZZ] = box[axis][axis]/ *zslices;
        gmx_rmpbc(gpbc, top->atoms.nr, box, x0);
        /*Reset Densslice every nsttblock steps*/
        /* The first conditional is for clang to understand that this branch is
         * always taken the first time. */
        if (Densslice == NULL || framenr % nsttblock == 0)
        {
            snew(Densslice, *xslices);
            for (i = 0; i < *xslices; i++)
            {
                snew(Densslice[i], *yslices);
                for (j = 0; j < *yslices; j++)
                {
                    snew(Densslice[i][j], *zslices);
                }
            }

            /* Allocate Memory to  extra frame in Densdevel -  rather stupid approach:
             * A single frame each time, although only every nsttblock steps.
             */
            srenew(*Densdevel, *tblock+1);
            (*Densdevel)[*tblock] = Densslice;
        }

        dscale = (*xslices)*(*yslices)*(*zslices)*AMU/ (box[ax1][ax1]*box[ax2][ax2]*box[axis][axis]*nsttblock*(NANO*NANO*NANO));

        if (bCenter)
        {
            center_coords(&top->atoms, box, x0, axis);
        }


        for (j = 0; j < gnx[0]; j++)
        {   /*Loop over all atoms in selected index*/
            x = x0[index[0][j]][ax1];
            y = x0[index[0][j]][ax2];
            z = x0[index[0][j]][axis];
            while (x < 0)
            {
                x += box[ax1][ax1];
            }
            while (x > box[ax1][ax1])
            {
                x -= box[ax1][ax1];
            }

            while (y < 0)
            {
                y += box[ax2][ax2];
            }
            while (y > box[ax2][ax2])
            {
                y -= box[ax2][ax2];
            }

            while (z < 0)
            {
                z += box[axis][axis];
            }
            while (z > box[axis][axis])
            {
                z -= box[axis][axis];
            }

            slicex = ((int) (x/bbww[XX])) % *xslices;
            slicey = ((int) (y/bbww[YY])) % *yslices;
            slicez = ((int) (z/bbww[ZZ])) % *zslices;
            Densslice[slicex][slicey][slicez] += (top->atoms.atom[index[0][j]].m*dscale);
        }

        framenr++;

        if (framenr % nsttblock == 0)
        {
            /*Implicit incrementation of Densdevel via renewal of Densslice*/
            /*only every nsttblock steps*/
            (*tblock)++;
        }

    }
    while (read_next_x(oenv, status, &t, x0, box));


    /*Free memory we no longer need and exit.*/
    gmx_rmpbc_done(gpbc);
    close_trj(status);

    if (0)
    {
        FILE *fp;
        fp = fopen("koko.xvg", "w");
        for (j = 0; (j < *zslices); j++)
        {
            fprintf(fp, "%5d", j);
            for (i = 0; (i < *tblock); i++)
            {
                fprintf(fp, "  %10g", (*Densdevel)[i][9][1][j]);
            }
            fprintf(fp, "\n");
        }
        fclose(fp);
    }

}

static void outputfield(const char *fldfn, real ****Densmap,
                        int xslices, int yslices, int zslices, int tdim)
{
/*Debug-filename and filehandle*/
    FILE *fldH;
    int   n, i, j, k;
    int   dim[4];
    real  totdens = 0;

    dim[0] = tdim;
    dim[1] = xslices;
    dim[2] = yslices;
    dim[3] = zslices;

    fldH = gmx_ffopen(fldfn, "w");
    fwrite(dim, sizeof(int), 4, fldH);
    for (n = 0; n < tdim; n++)
    {
        for (i = 0; i < xslices; i++)
        {
            for (j = 0; j < yslices; j++)
            {
                for (k = 0; k < zslices; k++)
                {
                    fwrite(&(Densmap[n][i][j][k]), sizeof(real), 1, fldH);
                    totdens += (Densmap[n][i][j][k]);
                }
            }
        }
    }
    totdens /= (xslices*yslices*zslices*tdim);
    fprintf(stderr, "Total density [kg/m^3]  %8f", totdens);
    gmx_ffclose(fldH);
}

static void filterdensmap(real ****Densmap, int xslices, int yslices, int zslices, int tblocks, int ftsize)
{
    real *kernel;
    real  std, var;
    int   i, j, n, order;
    order = ftsize/2;
    std   = ((real)order/2.0);
    var   = std*std;
    snew(kernel, ftsize);
    gausskernel(kernel, ftsize, var);
    for (n = 0; n < tblocks; n++)
    {
        for (i = 0; i < xslices; i++)
        {
            for (j = 0; j < yslices; j++)
            {
                periodic_convolution(zslices, Densmap[n][i][j], ftsize, kernel);
            }
        }
    }
}




static void interfaces_txy (real ****Densmap, int xslices, int yslices, int zslices,
                            int tblocks, real binwidth, int method,
                            real dens1, real dens2, t_interf ****intf1,
                            t_interf ****intf2, const output_env_t oenv)
{
    /*Returns two pointers to 3D arrays of t_interf structs containing (position,thickness) of the interface(s)*/
    FILE         *xvg;
    real         *zDensavg; /* zDensavg[z]*/
    int           i, j, k, n;
    int           xysize;
    int           ndx1, ndx2, *zperm;
    real          densmid;
    real          splitpoint, startpoint, endpoint;
    real         *sigma1, *sigma2;
    double        beginfit1[4];
    double        beginfit2[4];
    double       *fit1 = NULL, *fit2 = NULL;
    const double *avgfit1;
    const double *avgfit2;
    const real    onehalf = 1.00/2.00;
    t_interf   ***int1    = NULL, ***int2 = NULL; /*Interface matrices [t][x,y] - last index in row-major order*/
    /*Create int1(t,xy) and int2(t,xy) arrays with correct number of interf_t elements*/
    xysize = xslices*yslices;
    snew(int1, tblocks);
    snew(int2, tblocks);
    for (i = 0; i < tblocks; i++)
    {
        snew(int1[i], xysize);
        snew(int2[i], xysize);
        for (j = 0; j < xysize; j++)
        {
            snew(int1[i][j], 1);
            snew(int2[i][j], 1);
            init_interf(int1[i][j]);
            init_interf(int2[i][j]);
        }
    }

    if (method == methBISECT)
    {
        densmid = onehalf*(dens1+dens2);
        snew(zperm, zslices);
        for (n = 0; n < tblocks; n++)
        {
            for (i = 0; i < xslices; i++)
            {
                for (j = 0; j < yslices; j++)
                {
                    rangeArray(zperm, zslices); /*reset permutation array to identity*/
                    /*Binsearch returns slice-nr where the order param is  <= setpoint sgmid*/
                    ndx1 = start_binsearch(Densmap[n][i][j], zperm, 0, zslices/2-1, densmid, 1);
                    ndx2 = start_binsearch(Densmap[n][i][j], zperm, zslices/2, zslices-1, densmid, -1);

                    /* Linear interpolation (for use later if time allows)
                     * rho_1s= Densmap[n][i][j][zperm[ndx1]]
                     * rho_1e =Densmap[n][i][j][zperm[ndx1+1]] - in worst case might be far off
                     * rho_2s =Densmap[n][i][j][zperm[ndx2+1]]
                     * rho_2e =Densmap[n][i][j][zperm[ndx2]]
                     * For 1st interface we have:
                       densl= Densmap[n][i][j][zperm[ndx1]];
                       densr= Densmap[n][i][j][zperm[ndx1+1]];
                       alpha=(densmid-densl)/(densr-densl);
                       deltandx=zperm[ndx1+1]-zperm[ndx1];

                       if(debug){
                       printf("Alpha, Deltandx  %f %i\n", alpha,deltandx);
                       }
                       if(abs(alpha)>1.0 || abs(deltandx)>3){
                       pos=zperm[ndx1];
                       spread=-1;
                       }
                       else {
                       pos=zperm[ndx1]+alpha*deltandx;
                       spread=binwidth*deltandx;
                       }
                     * For the 2nd interface  can use the same formulation, since alpha should become negative ie:
                     * alpha=(densmid-Densmap[n][i][j][zperm[ndx2]])/(Densmap[n][i][j][zperm[nxd2+1]]-Densmap[n][i][j][zperm[ndx2]]);
                     * deltandx=zperm[ndx2+1]-zperm[ndx2];
                     * pos=zperm[ndx2]+alpha*deltandx;   */

                    /*After filtering we use the direct approach	*/
                    int1[n][j+(i*yslices)]->Z = (zperm[ndx1]+onehalf)*binwidth;
                    int1[n][j+(i*yslices)]->t = binwidth;
                    int2[n][j+(i*yslices)]->Z = (zperm[ndx2]+onehalf)*binwidth;
                    int2[n][j+(i*yslices)]->t = binwidth;
                }
            }
        }
    }

    if (method == methFUNCFIT)
    {
        /*Assume a box divided in 2 along midpoint of z for starters*/
        startpoint = 0.0;
        endpoint   = binwidth*zslices;
        splitpoint = (startpoint+endpoint)/2.0;
        /*Initial fit proposals*/
        beginfit1[0] = dens1;
        beginfit1[1] = dens2;
        beginfit1[2] = (splitpoint/2);
        beginfit1[3] = 0.5;

        beginfit2[0] = dens2;
        beginfit2[1] = dens1;
        beginfit2[2] = (3*splitpoint/2);
        beginfit2[3] = 0.5;

        snew(zDensavg, zslices);
        snew(sigma1, zslices);
        snew(sigma2, zslices);

        for (k = 0; k < zslices; k++)
        {
            sigma1[k] = sigma2[k] = 1;
        }
        /*Calculate average density along z - avoid smoothing by using coarse-grained-mesh*/
        for (k = 0; k < zslices; k++)
        {
            for (n = 0; n < tblocks; n++)
            {
                for (i = 0; i < xslices; i++)
                {
                    for (j = 0; j < yslices; j++)
                    {
                        zDensavg[k] += (Densmap[n][i][j][k]/(xslices*yslices*tblocks));
                    }
                }
            }
        }

        if (debug)
        {
            xvg = xvgropen("DensprofileonZ.xvg", "Averaged Densityprofile on Z", "z[nm]", "Density[kg/m^3]", oenv);
            for (k = 0; k < zslices; k++)
            {
                fprintf(xvg, "%4f.3   %8f.4\n", k*binwidth, zDensavg[k]);
            }
            xvgrclose(xvg);
        }

        /*Fit average density in z over whole trajectory to obtain tentative fit-parameters in fit1 and fit2*/

        /*Fit 1st half of box*/
        do_lmfit(zslices, zDensavg, sigma1, binwidth, NULL, startpoint, splitpoint, oenv, FALSE, effnERF, beginfit1, 8, NULL);
        /*Fit 2nd half of box*/
        do_lmfit(zslices, zDensavg, sigma2, binwidth, NULL, splitpoint, endpoint, oenv, FALSE, effnERF, beginfit2, 8, NULL);

        /*Initialise the const arrays for storing the average fit parameters*/
        avgfit1 = beginfit1;
        avgfit2 = beginfit2;



        /*Now do fit over each x  y and t slice to get Zint(x,y,t) - loop is very large, we potentially should average over time directly*/
        for (n = 0; n < tblocks; n++)
        {
            for (i = 0; i < xslices; i++)
            {
                for (j = 0; j < yslices; j++)
                {
                    /*Reinitialise fit for each mesh-point*/
                    srenew(fit1, 4);
                    srenew(fit2, 4);
                    for (k = 0; k < 4; k++)
                    {
                        fit1[k] = avgfit1[k];
                        fit2[k] = avgfit2[k];
                    }
                    /*Now fit and store in structures in row-major order int[n][i][j]*/
                    do_lmfit(zslices, Densmap[n][i][j], sigma1, binwidth, NULL, startpoint, splitpoint, oenv, FALSE, effnERF, fit1, 0, NULL);
                    int1[n][j+(yslices*i)]->Z = fit1[2];
                    int1[n][j+(yslices*i)]->t = fit1[3];
                    do_lmfit(zslices, Densmap[n][i][j], sigma2, binwidth, NULL, splitpoint, endpoint, oenv, FALSE, effnERF, fit2, 0, NULL);
                    int2[n][j+(yslices*i)]->Z = fit2[2];
                    int2[n][j+(yslices*i)]->t = fit2[3];
                }
            }
        }
    }


    *intf1 = int1;
    *intf2 = int2;

}

static void writesurftoxpms(t_interf ***surf1, t_interf ***surf2, int tblocks, int xbins, int ybins, int zbins, real bw, real bwz, char **outfiles, int maplevels )
{
    char   numbuf[13];
    int    n, i, j;
    real **profile1, **profile2;
    real   max1, max2, min1, min2, *xticks, *yticks;
    t_rgb  lo = {0, 0, 0};
    t_rgb  hi = {1, 1, 1};
    FILE  *xpmfile1, *xpmfile2;

/*Prepare xpm structures for output*/

/*Allocate memory to tick's and matrices*/
    snew (xticks, xbins+1);
    snew (yticks, ybins+1);

    profile1 = mk_matrix(xbins, ybins, FALSE);
    profile2 = mk_matrix(xbins, ybins, FALSE);

    for (i = 0; i < xbins+1; i++)
    {
        xticks[i] += bw;
    }
    for (j = 0; j < ybins+1; j++)
    {
        yticks[j] += bw;
    }

    xpmfile1 = gmx_ffopen(outfiles[0], "w");
    xpmfile2 = gmx_ffopen(outfiles[1], "w");

    max1 = max2 = 0.0;
    min1 = min2 = zbins*bwz;

    for (n = 0; n < tblocks; n++)
    {
        sprintf(numbuf, "tblock: %4i", n);
/*Filling matrices for inclusion in xpm-files*/
        for (i = 0; i < xbins; i++)
        {
            for (j = 0; j < ybins; j++)
            {
                profile1[i][j] = (surf1[n][j+ybins*i])->Z;
                profile2[i][j] = (surf2[n][j+ybins*i])->Z;
                /*Finding max and min values*/
                if (profile1[i][j] > max1)
                {
                    max1 = profile1[i][j];
                }
                if (profile1[i][j] < min1)
                {
                    min1 = profile1[i][j];
                }
                if (profile2[i][j] > max2)
                {
                    max2 = profile2[i][j];
                }
                if (profile2[i][j] < min2)
                {
                    min2 = profile2[i][j];
                }
            }
        }

        write_xpm(xpmfile1, 3, numbuf, "Height", "x[nm]", "y[nm]", xbins, ybins, xticks, yticks, profile1, min1, max1, lo, hi, &maplevels);
        write_xpm(xpmfile2, 3, numbuf, "Height", "x[nm]", "y[nm]", xbins, ybins, xticks, yticks, profile2, min2, max2, lo, hi, &maplevels);
    }

    gmx_ffclose(xpmfile1);
    gmx_ffclose(xpmfile2);


    sfree(profile1);
    sfree(profile2);
    sfree(xticks);
    sfree(yticks);
}

static void writeraw(t_interf ***int1, t_interf ***int2, int tblocks,
                     int xbins, int ybins, char **fnms,
                     const output_env_t oenv)
{
    FILE *raw1, *raw2;
    int   i, j, n;

    raw1 = gmx_ffopen(fnms[0], "w");
    raw2 = gmx_ffopen(fnms[1], "w");
    try
    {
        gmx::BinaryInformationSettings settings;
        settings.generatedByHeader(true);
        settings.linePrefix("# ");
        gmx::printBinaryInformation(raw1, output_env_get_program_context(oenv),
                                    settings);
        gmx::printBinaryInformation(raw2, output_env_get_program_context(oenv),
                                    settings);
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
    fprintf(raw1, "# Legend: nt nx ny\n# Xbin Ybin Z t\n");
    fprintf(raw2, "# Legend: nt nx ny\n# Xbin Ybin Z t\n");
    fprintf(raw1, "%i %i %i\n", tblocks, xbins, ybins);
    fprintf(raw2, "%i %i %i\n", tblocks, xbins, ybins);
    for (n = 0; n < tblocks; n++)
    {
        for (i = 0; i < xbins; i++)
        {
            for (j = 0; j < ybins; j++)
            {
                fprintf(raw1, "%i  %i  %8.5f  %6.4f\n", i, j, (int1[n][j+ybins*i])->Z, (int1[n][j+ybins*i])->t);
                fprintf(raw2, "%i  %i  %8.5f  %6.4f\n", i, j, (int2[n][j+ybins*i])->Z, (int2[n][j+ybins*i])->t);
            }
        }
    }

    gmx_ffclose(raw1);
    gmx_ffclose(raw2);
}



int gmx_densorder(int argc, char *argv[])
{
    static const char *desc[] = {
        "[THISMODULE] reduces a two-phase density distribution",
        "along an axis, computed over a MD trajectory,",
        "to 2D surfaces fluctuating in time, by a fit to",
        "a functional profile for interfacial densities.",
        "A time-averaged spatial representation of the",
        "interfaces can be output with the option [TT]-tavg[tt]."
    };

    /* Extra arguments - but note how you always get the begin/end
     * options when running the program, without mentioning them here!
     */

    output_env_t       oenv;
    t_topology        *top;
    char             **grpname;
    int                ePBC, *ngx;
    static real        binw      = 0.2;
    static real        binwz     = 0.05;
    static real        dens1     = 0.00;
    static real        dens2     = 1000.00;
    static int         ftorder   = 0;
    static int         nsttblock = 100;
    static int         axis      = 2;
    static const char *axtitle   = "Z";
    atom_id          **index; /* Index list for single group*/
    int                xslices, yslices, zslices, tblock;
    static gmx_bool    bGraph   = FALSE;
    static gmx_bool    bCenter  = FALSE;
    static gmx_bool    bFourier = FALSE;
    static gmx_bool    bRawOut  = FALSE;
    static gmx_bool    bOut     = FALSE;
    static gmx_bool    b1d      = FALSE;
    static int         nlevels  = 100;
    /*Densitymap - Densmap[t][x][y][z]*/
    real           ****Densmap = NULL;
    /* Surfaces surf[t][surf_x,surf_y]*/
    t_interf        ***surf1, ***surf2;

    static const char *meth[] = {NULL, "bisect", "functional", NULL};
    int                eMeth;

    char             **graphfiles, **rawfiles, **spectra; /* Filenames for xpm-surface maps, rawdata and powerspectra */
    int                nfxpm = -1, nfraw, nfspect;        /* # files for interface maps and spectra = # interfaces */

    t_pargs            pa[] = {
        { "-1d", FALSE, etBOOL, {&b1d},
          "Pseudo-1d interface geometry"},
        { "-bw", FALSE, etREAL, {&binw},
          "Binwidth of density distribution tangential to interface"},
        { "-bwn", FALSE, etREAL, {&binwz},
          "Binwidth of density distribution normal to interface"},
        { "-order", FALSE, etINT, {&ftorder},
          "Order of Gaussian filter, order 0 equates to NO filtering"},
        {"-axis", FALSE, etSTR, {&axtitle},
         "Axis Direction - X, Y or Z"},
        {"-method", FALSE, etENUM, {meth},
         "Interface location method"},
        {"-d1", FALSE, etREAL, {&dens1},
         "Bulk density phase 1 (at small z)"},
        {"-d2", FALSE, etREAL, {&dens2},
         "Bulk density phase 2 (at large z)"},
        { "-tblock", FALSE, etINT, {&nsttblock},
          "Number of frames in one time-block average"},
        { "-nlevel", FALSE, etINT, {&nlevels},
          "Number of Height levels in 2D - XPixMaps"}
    };


    t_filenm fnm[] = {
        { efTPR, "-s",  NULL, ffREAD },               /* this is for the topology */
        { efTRX, "-f", NULL, ffREAD },                /* and this for the trajectory */
        { efNDX, "-n", NULL, ffREAD},                 /* this is to select groups */
        { efDAT, "-o", "Density4D", ffOPTWR},         /* This is for outputting the entire 4D densityfield in binary format */
        { efOUT, "-or", NULL, ffOPTWRMULT},           /* This is for writing out the entire information in the t_interf arrays */
        { efXPM, "-og", "interface", ffOPTWRMULT},    /* This is for writing out the interface meshes - one xpm-file per tblock*/
        { efOUT, "-Spect", "intfspect", ffOPTWRMULT}, /* This is for the trajectory averaged Fourier-spectra*/
    };

#define NFILE asize(fnm)

    /* This is the routine responsible for adding default options,
     * calling the X/motif interface, etc. */
    if (!parse_common_args(&argc, argv, PCA_CAN_TIME | PCA_CAN_VIEW,
                           NFILE, fnm, asize(pa), pa, asize(desc), desc, 0, NULL, &oenv))
    {
        return 0;
    }


    eMeth    = nenum(meth);
    bFourier = opt2bSet("-Spect", NFILE, fnm);
    bRawOut  = opt2bSet("-or", NFILE, fnm);
    bGraph   = opt2bSet("-og", NFILE, fnm);
    bOut     = opt2bSet("-o", NFILE, fnm);
    top      = read_top(ftp2fn(efTPR, NFILE, fnm), &ePBC);
    snew(grpname, 1);
    snew(index, 1);
    snew(ngx, 1);

/* Calculate axis */
    axis = toupper(axtitle[0]) - 'X';

    get_index(&top->atoms, ftp2fn_null(efNDX, NFILE, fnm), 1, ngx, index, grpname);

    density_in_time(ftp2fn(efTRX, NFILE, fnm), index, ngx, binw, binwz, nsttblock, &Densmap, &xslices, &yslices, &zslices, &tblock, top, ePBC, axis, bCenter, b1d, oenv);

    if (ftorder > 0)
    {
        filterdensmap(Densmap, xslices, yslices, zslices, tblock, 2*ftorder+1);
    }

    if (bOut)
    {
        outputfield(opt2fn("-o", NFILE, fnm), Densmap, xslices, yslices, zslices, tblock);
    }

    interfaces_txy(Densmap, xslices, yslices, zslices, tblock, binwz, eMeth, dens1, dens2, &surf1, &surf2, oenv);

    if (bGraph)
    {

        /*Output surface-xpms*/
        nfxpm = opt2fns(&graphfiles, "-og", NFILE, fnm);
        if (nfxpm != 2)
        {
            gmx_fatal(FARGS, "No or not correct number (2) of output-files: %d", nfxpm);
        }
        writesurftoxpms(surf1, surf2, tblock, xslices, yslices, zslices, binw, binwz, graphfiles, zslices);
    }





/*Output raw-data*/
    if (bRawOut)
    {
        nfraw = opt2fns(&rawfiles, "-or", NFILE, fnm);
        if (nfraw != 2)
        {
            gmx_fatal(FARGS, "No or not correct number (2) of output-files: %d", nfxpm);
        }
        writeraw(surf1, surf2, tblock, xslices, yslices, rawfiles, oenv);
    }



    if (bFourier)
    {
        nfspect = opt2fns(&spectra, "-Spect", NFILE, fnm);
        if (nfspect != 2)
        {
            gmx_fatal(FARGS, "No or not correct number (2) of output-file-series: %d",
                      nfspect);
        }
        powerspectavg_intf(surf1, surf2, tblock, xslices, yslices, spectra);
    }

    sfree(Densmap);
    if (bGraph || bFourier || bRawOut)
    {
        sfree(surf1);
        sfree(surf2);
    }

    return 0;
}
