/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2010- The GROMACS Authors
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
#include <cstdio>
#include <cstring>

#include <filesystem>
#include <string>

#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/filetypes.h"
#include "gromacs/fileio/matio.h"
#include "gromacs/fileio/rgb.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/gmxana/binsearch.h"
#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/gmxana/gstat.h"
#include "gromacs/gmxana/powerspect.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pbcutil/rmpbc.h"
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
#include "gromacs/utility/vec.h"
#include "gromacs/utility/vectypes.h"

enum class PbcType : int;
struct gmx_output_env_t;

static void find_tetra_order_grid(t_topology top,
                                  PbcType    pbcType,
                                  int        natoms,
                                  matrix     box,
                                  rvec       x[],
                                  int        maxidx,
                                  const int  index[],
                                  real*      sgmean,
                                  real*      skmean,
                                  int        nslicex,
                                  int        nslicey,
                                  int        nslicez,
                                  real***    sggrid,
                                  real***    skgrid)
{
    int         ix, jx, i, j, k, *nn[4];
    rvec        dx, rj, rk, urk, urj;
    real        cost, cost2, *sgmol, *skmol, rmean, rmean2, r2, box2, *r_nn[4];
    t_pbc       pbc;
    int         slindex_x, slindex_y, slindex_z;
    int***      sl_count;
    real        onethird = 1.0 / 3.0;
    gmx_rmpbc_t gpbc;

    /*  dmat = init_mat(maxidx, FALSE); */

    box2 = box[XX][XX] * box[XX][XX];

    /* Initialize expanded sl_count array */
    snew(sl_count, nslicex);
    for (i = 0; i < nslicex; i++)
    {
        snew(sl_count[i], nslicey);
        for (j = 0; j < nslicey; j++)
        {
            snew(sl_count[i][j], nslicez);
        }
    }


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
    set_pbc(&pbc, pbcType, box);
    gpbc = gmx_rmpbc_init(&top.idef, pbcType, natoms);
    gmx_rmpbc_apply(gpbc, natoms, box, x);

    *sgmean = 0.0;
    *skmean = 0.0;
    for (i = 0; (i < maxidx); i++)
    { /* loop over index file */
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
                r_nn[3][i] = r_nn[2][i];
                nn[3][i]   = nn[2][i];
                r_nn[2][i] = r_nn[1][i];
                nn[2][i]   = nn[1][i];
                r_nn[1][i] = r_nn[0][i];
                nn[1][i]   = nn[0][i];
                r_nn[0][i] = r2;
                nn[0][i]   = j;
            }
            else if (r2 < r_nn[1][i])
            {
                r_nn[3][i] = r_nn[2][i];
                nn[3][i]   = nn[2][i];
                r_nn[2][i] = r_nn[1][i];
                nn[2][i]   = nn[1][i];
                r_nn[1][i] = r2;
                nn[1][i]   = j;
            }
            else if (r2 < r_nn[2][i])
            {
                r_nn[3][i] = r_nn[2][i];
                nn[3][i]   = nn[2][i];
                r_nn[2][i] = r2;
                nn[2][i]   = j;
            }
            else if (r2 < r_nn[3][i])
            {
                r_nn[3][i] = r2;
                nn[3][i]   = j;
            }
        }


        /* calculate mean distance between nearest neighbours */
        rmean = 0;
        for (j = 0; (j < 4); j++)
        {
            r_nn[j][i] = std::sqrt(r_nn[j][i]);
            rmean += r_nn[j][i];
        }
        rmean /= 4;

        sgmol[i] = 0.0;
        skmol[i] = 0.0;

        /* Chau1998a eqn 3 */
        /* angular part tetrahedrality order parameter per atom */
        for (j = 0; (j < 3); j++)
        {
            for (k = j + 1; (k < 4); k++)
            {
                pbc_dx(&pbc, x[ix], x[index[nn[k][i]]], rk);
                pbc_dx(&pbc, x[ix], x[index[nn[j][i]]], rj);

                unitv(rk, urk);
                unitv(rj, urj);

                cost  = iprod(urk, urj) + onethird;
                cost2 = cost * cost;

                sgmol[i] += cost2;
            }
        }
        /* normalize sgmol between 0.0 and 1.0 */
        sgmol[i] = 3 * sgmol[i] / 32;
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

        /* Compute sliced stuff in x y z*/
        slindex_x = static_cast<int>(std::round((1 + x[i][XX] / box[XX][XX]) * nslicex)) % nslicex;
        slindex_y = static_cast<int>(std::round((1 + x[i][YY] / box[YY][YY]) * nslicey)) % nslicey;
        slindex_z = static_cast<int>(std::round((1 + x[i][ZZ] / box[ZZ][ZZ]) * nslicez)) % nslicez;
        sggrid[slindex_x][slindex_y][slindex_z] += sgmol[i];
        skgrid[slindex_x][slindex_y][slindex_z] += skmol[i];
        (sl_count[slindex_x][slindex_y][slindex_z])++;
    } /* loop over entries in index file */

    *sgmean /= maxidx;
    *skmean /= maxidx;

    for (i = 0; (i < nslicex); i++)
    {
        for (j = 0; j < nslicey; j++)
        {
            for (k = 0; k < nslicez; k++)
            {
                if (sl_count[i][j][k] > 0)
                {
                    sggrid[i][j][k] /= sl_count[i][j][k];
                    skgrid[i][j][k] /= sl_count[i][j][k];
                }
            }
        }
    }

    sfree(sl_count);
    sfree(sgmol);
    sfree(skmol);
    for (i = 0; (i < 4); i++)
    {
        sfree(r_nn[i]);
        sfree(nn[i]);
    }
}

/*Determines interface from tetrahedral order parameter in box with specified binwidth.  */
/*Outputs interface positions(bins), the number of timeframes, and the number of surface-mesh points in xy*/

static void calc_tetra_order_interface(const char*       fnNDX,
                                       const char*       fnTPS,
                                       const char*       fnTRX,
                                       real              binw,
                                       int               tblock,
                                       int*              nframes,
                                       int*              nslicex,
                                       int*              nslicey,
                                       real              sgang1,
                                       real              sgang2,
                                       real****          intfpos,
                                       gmx_output_env_t* oenv)
{
    FILE *       fpsg = nullptr, *fpsk = nullptr;
    t_topology   top;
    PbcType      pbcType;
    t_trxstatus* status;
    int          natoms;
    real         t;
    rvec *       xtop, *x;
    matrix       box;
    real         sg, sk, sgintf;
    int**        index   = nullptr;
    char**       grpname = nullptr;
    int          i, j, k, n, *isize, ng, nslicez, framenr;
    real ***sg_grid = nullptr, ***sk_grid = nullptr, ***sg_fravg = nullptr, ***sk_fravg = nullptr,
         ****sk_4d = nullptr, ****sg_4d = nullptr;
    int*       perm;
    int        ndx1, ndx2;
    int        bins;
    const real onehalf = 1.0 / 2.0;
    /* real   ***intfpos[2]; pointers to arrays of two interface positions zcoord(framenr,xbin,ybin): intfpos[interface_index][t][nslicey*x+y]
     * i.e 1D Row-major order in (t,x,y) */


    read_tps_conf(fnTPS, &top, &pbcType, &xtop, nullptr, box, FALSE);

    *nslicex = static_cast<int>(box[XX][XX] / binw + onehalf); /*Calculate slicenr from binwidth*/
    *nslicey = static_cast<int>(box[YY][YY] / binw + onehalf);
    nslicez  = static_cast<int>(box[ZZ][ZZ] / binw + onehalf);


    ng = 1;
    /* get index groups */
    printf("Select the group that contains the atoms you want to use for the tetrahedrality order "
           "parameter calculation:\n");
    snew(grpname, ng);
    snew(index, ng);
    snew(isize, ng);
    get_index(&top.atoms, fnNDX, ng, isize, index, grpname);

    /* Analyze trajectory */
    natoms = read_first_x(oenv, &status, fnTRX, &t, &x, box);
    if (natoms > top.atoms.nr)
    {
        gmx_fatal(FARGS, "Topology (%d atoms) does not match trajectory (%d atoms)", top.atoms.nr, natoms);
    }
    check_index(nullptr, ng, index[0], nullptr, natoms);


    /*Prepare structures for temporary storage of frame info*/
    snew(sg_grid, *nslicex);
    snew(sk_grid, *nslicex);
    for (i = 0; i < *nslicex; i++)
    {
        snew(sg_grid[i], *nslicey);
        snew(sk_grid[i], *nslicey);
        for (j = 0; j < *nslicey; j++)
        {
            snew(sg_grid[i][j], nslicez);
            snew(sk_grid[i][j], nslicez);
        }
    }

    sg_4d    = nullptr;
    sk_4d    = nullptr;
    *nframes = 0;
    framenr  = 0;

    /* Loop over frames*/
    do
    {
        /*Initialize box meshes (temporary storage for each tblock frame -reinitialise every tblock steps */
        if (framenr % tblock == 0)
        {
            srenew(sk_4d, *nframes + 1);
            srenew(sg_4d, *nframes + 1);
            snew(sg_fravg, *nslicex);
            snew(sk_fravg, *nslicex);
            for (i = 0; i < *nslicex; i++)
            {
                snew(sg_fravg[i], *nslicey);
                snew(sk_fravg[i], *nslicey);
                for (j = 0; j < *nslicey; j++)
                {
                    snew(sg_fravg[i][j], nslicez);
                    snew(sk_fravg[i][j], nslicez);
                }
            }
        }

        find_tetra_order_grid(
                top, pbcType, natoms, box, x, isize[0], index[0], &sg, &sk, *nslicex, *nslicey, nslicez, sg_grid, sk_grid);
        GMX_RELEASE_ASSERT(sk_fravg != nullptr, "Trying to dereference NULL sk_fravg pointer");
        for (i = 0; i < *nslicex; i++)
        {
            for (j = 0; j < *nslicey; j++)
            {
                for (k = 0; k < nslicez; k++)
                {
                    sk_fravg[i][j][k] += sk_grid[i][j][k] / tblock;
                    sg_fravg[i][j][k] += sg_grid[i][j][k] / tblock;
                }
            }
        }

        framenr++;

        if (framenr % tblock == 0)
        {
            GMX_RELEASE_ASSERT(sk_4d != nullptr, "Trying to dereference NULL sk_4d pointer");
            sk_4d[*nframes] = sk_fravg;
            sg_4d[*nframes] = sg_fravg;
            (*nframes)++;
        }

    } while (read_next_x(oenv, status, &t, x, box));
    close_trx(status);

    sfree(grpname);
    sfree(index);
    sfree(isize);

    /*Debugging for printing out the entire order parameter meshes.*/
    if (debug)
    {
        fpsg = xvgropen(
                "sg_ang_mesh", "S\\sg\\N Angle Order Parameter / Meshpoint", "(nm)", "S\\sg\\N", oenv);
        fpsk = xvgropen(
                "sk_dist_mesh", "S\\sk\\N Distance Order Parameter / Meshpoint", "(nm)", "S\\sk\\N", oenv);
        for (n = 0; n < (*nframes); n++)
        {
            fprintf(fpsg, "%i\n", n);
            fprintf(fpsk, "%i\n", n);
            for (i = 0; (i < *nslicex); i++)
            {
                for (j = 0; j < *nslicey; j++)
                {
                    for (k = 0; k < nslicez; k++)
                    {
                        fprintf(fpsg,
                                "%4f  %4f  %4f  %8f\n",
                                (i + 0.5) * box[XX][XX] / (*nslicex),
                                (j + 0.5) * box[YY][YY] / (*nslicey),
                                (k + 0.5) * box[ZZ][ZZ] / nslicez,
                                sg_4d[n][i][j][k]);
                        fprintf(fpsk,
                                "%4f  %4f  %4f  %8f\n",
                                (i + 0.5) * box[XX][XX] / (*nslicex),
                                (j + 0.5) * box[YY][YY] / (*nslicey),
                                (k + 0.5) * box[ZZ][ZZ] / nslicez,
                                sk_4d[n][i][j][k]);
                    }
                }
            }
        }
        xvgrclose(fpsg);
        xvgrclose(fpsk);
    }


    /* Find positions of interface z by scanning orderparam for each frame and for each xy-mesh cylinder along z*/

    /*Simple trial: assume interface is in the middle of -sgang1 and sgang2*/
    sgintf = 0.5 * (sgang1 + sgang2);


    /*Allocate memory for interface arrays; */
    snew((*intfpos), 2);
    snew((*intfpos)[0], *nframes);
    snew((*intfpos)[1], *nframes);

    bins = (*nslicex) * (*nslicey);


    snew(perm, nslicez); /*permutation array for sorting along normal coordinate*/


    for (n = 0; n < *nframes; n++)
    {
        snew((*intfpos)[0][n], bins);
        snew((*intfpos)[1][n], bins);
        for (i = 0; i < *nslicex; i++)
        {
            for (j = 0; j < *nslicey; j++)
            {
                rangeArray(perm, nslicez); /*reset permutation array to identity*/
                /*Binsearch returns 2 bin-numbers where the order param is  <= setpoint sgintf*/
                ndx1 = start_binsearch(sg_4d[n][i][j], perm, 0, nslicez / 2 - 1, sgintf, 1);
                ndx2 = start_binsearch(sg_4d[n][i][j], perm, nslicez / 2, nslicez - 1, sgintf, -1);
                /*Use linear interpolation to smooth out the interface position*/

                /*left interface (0)*/
                /*if((sg_4d[n][i][j][perm[ndx1+1]]-sg_4d[n][i][j][perm[ndx1]])/sg_4d[n][i][j][perm[ndx1]] > 0.01){
                   pos=( (sgintf-sg_4d[n][i][j][perm[ndx1]])*perm[ndx1+1]+(sg_4d[n][i][j][perm[ndx1+1]]-sgintf)*perm[ndx1 ])*/
                (*intfpos)[0][n][j + *nslicey * i] = (perm[ndx1] + onehalf) * binw;
                /*right interface (1)*/
                /*alpha=(sgintf-sg_4d[n][i][j][perm[ndx2]])/(sg_4d[n][i][j][perm[ndx2]+1]-sg_4d[n][i][j][perm[ndx2]]);*/
                /*(*intfpos)[1][n][j+*nslicey*i]=((1-alpha)*perm[ndx2]+alpha*(perm[ndx2]+1)+onehalf)*box[ZZ][ZZ]/nslicez;*/
                (*intfpos)[1][n][j + *nslicey * i] = (perm[ndx2] + onehalf) * binw;
            }
        }
    }


    sfree(sk_4d);
    sfree(sg_4d);
}

static void writesurftoxpms(real***                          surf,
                            int                              tblocks,
                            int                              xbins,
                            int                              ybins,
                            real                             bw,
                            gmx::ArrayRef<const std::string> outfiles,
                            int                              maplevels)
{

    char   numbuf[STRLEN];
    int    n, i, j;
    real **profile1, **profile2;
    real   max1, max2, min1, min2, *xticks, *yticks;
    t_rgb  lo = { 1, 1, 1 };
    t_rgb  hi = { 0, 0, 0 };
    FILE * xpmfile1, *xpmfile2;

    /*Prepare xpm structures for output*/

    /*Allocate memory to tick's and matrices*/
    snew(xticks, xbins + 1);
    snew(yticks, ybins + 1);

    profile1 = mk_matrix(xbins, ybins, FALSE);
    profile2 = mk_matrix(xbins, ybins, FALSE);

    for (i = 0; i < xbins + 1; i++)
    {
        xticks[i] += bw;
    }
    for (j = 0; j < ybins + 1; j++)
    {
        yticks[j] += bw;
    }

    xpmfile1 = gmx_ffopen(outfiles[0], "w");
    xpmfile2 = gmx_ffopen(outfiles[1], "w");

    max1 = max2 = 0.0;
    min1 = min2 = 1000.00;

    for (n = 0; n < tblocks; n++)
    {
        sprintf(numbuf, "%5d", n);
        /*Filling matrices for inclusion in xpm-files*/
        for (i = 0; i < xbins; i++)
        {
            for (j = 0; j < ybins; j++)
            {
                profile1[i][j] = (surf[0][n][j + ybins * i]);
                profile2[i][j] = (surf[1][n][j + ybins * i]);
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


static void writeraw(real*** surf, int tblocks, int xbins, int ybins, gmx::ArrayRef<const std::string> fnms)
{
    FILE *raw1, *raw2;
    int   i, j, n;

    raw1 = gmx_ffopen(fnms[0], "w");
    raw2 = gmx_ffopen(fnms[1], "w");
    fprintf(raw1, "#Legend\n#TBlock\n#Xbin Ybin Z t\n");
    fprintf(raw2, "#Legend\n#TBlock\n#Xbin Ybin Z t\n");
    for (n = 0; n < tblocks; n++)
    {
        fprintf(raw1, "%5d\n", n);
        fprintf(raw2, "%5d\n", n);
        for (i = 0; i < xbins; i++)
        {
            for (j = 0; j < ybins; j++)
            {
                fprintf(raw1, "%i  %i  %8.5f\n", i, j, (surf[0][n][j + ybins * i]));
                fprintf(raw2, "%i  %i  %8.5f\n", i, j, (surf[1][n][j + ybins * i]));
            }
        }
    }

    gmx_ffclose(raw1);
    gmx_ffclose(raw2);
}


int gmx_hydorder(int argc, char* argv[])
{
    static const char* desc[] = {
        "[THISMODULE] computes the tetrahedrality order parameters around a ",
        "given atom. Both angle an distance order parameters are calculated. See",
        "P.-L. Chau and A.J. Hardwick, Mol. Phys., 93, (1998), 511-518.",
        "for more details.[PAR]",
        "[THISMODULE] calculates the order parameter in a 3d-mesh in the box, and",
        "with 2 phases in the box gives the user the option to define a 2D interface in time",
        "separating the faces by specifying parameters [TT]-sgang1[tt] and",
        "[TT]-sgang2[tt] (it is important to select these judiciously)."
    };

    int                axis      = 0;
    static int         nsttblock = 1;
    static int         nlevels   = 100;
    static real        binwidth  = 1.0; /* binwidth in mesh           */
    static real        sg1       = 1;
    static real        sg2       = 1; /* order parameters for bulk phases */
    static gmx_bool    bFourier  = FALSE;
    static gmx_bool    bRawOut   = FALSE;
    int                frames, xslices, yslices; /* Dimensions of interface arrays*/
    real***            intfpos; /* Interface arrays (intfnr,t,xy) -potentially large */
    static const char* normal_axis[] = { nullptr, "z", "x", "y", nullptr };

    t_pargs pa[] = {
        { "-d", FALSE, etENUM, { normal_axis }, "Direction of the normal on the membrane" },
        { "-bw", FALSE, etREAL, { &binwidth }, "Binwidth of box mesh" },
        { "-sgang1", FALSE, etREAL, { &sg1 }, "tetrahedral angle parameter in Phase 1 (bulk)" },
        { "-sgang2", FALSE, etREAL, { &sg2 }, "tetrahedral angle parameter in Phase 2 (bulk)" },
        { "-tblock", FALSE, etINT, { &nsttblock }, "Number of frames in one time-block average" },
        { "-nlevel", FALSE, etINT, { &nlevels }, "Number of Height levels in 2D - XPixMaps" }
    };

    t_filenm fnm[] = {
        /* files for gmx order    */
        { efTRX, "-f", nullptr, ffREAD },              /* trajectory file              */
        { efNDX, "-n", nullptr, ffREAD },              /* index file           */
        { efTPR, "-s", nullptr, ffREAD },              /* topology file                */
        { efXPM, "-o", "intf", ffWRMULT },             /* XPM- surface maps	*/
        { efOUT, "-or", "raw", ffOPTWRMULT },          /* xvgr output file           */
        { efOUT, "-Spect", "intfspect", ffOPTWRMULT }, /* Fourier spectrum interfaces */
    };
#define NFILE asize(fnm)

    /*Filenames*/
    const char *      ndxfnm, *tpsfnm, *trxfnm;
    gmx_output_env_t* oenv;

    if (!parse_common_args(
                &argc, argv, PCA_CAN_VIEW | PCA_CAN_TIME, NFILE, fnm, asize(pa), pa, asize(desc), desc, 0, nullptr, &oenv))
    {
        return 0;
    }
    bFourier = opt2bSet("-Spect", NFILE, fnm);
    bRawOut  = opt2bSet("-or", NFILE, fnm);

    if (binwidth < 0.0)
    {
        gmx_fatal(FARGS, "Can not have binwidth < 0");
    }

    ndxfnm = ftp2fn(efNDX, NFILE, fnm);
    tpsfnm = ftp2fn(efTPR, NFILE, fnm);
    trxfnm = ftp2fn(efTRX, NFILE, fnm);

    /* Calculate axis */
    GMX_RELEASE_ASSERT(normal_axis[0] != nullptr,
                       "Option setting inconsistency; normal_axis[0] is NULL");
    if (std::strcmp(normal_axis[0], "x") == 0)
    {
        axis = XX;
    }
    else if (std::strcmp(normal_axis[0], "y") == 0)
    {
        axis = YY;
    }
    else if (std::strcmp(normal_axis[0], "z") == 0)
    {
        axis = ZZ;
    }
    else
    {
        gmx_fatal(FARGS, "Invalid axis, use x, y or z");
    }

    switch (axis)
    {
        case 0: fprintf(stderr, "Taking x axis as normal to the membrane\n"); break;
        case 1: fprintf(stderr, "Taking y axis as normal to the membrane\n"); break;
        case 2: fprintf(stderr, "Taking z axis as normal to the membrane\n"); break;
    }

    /* tetraheder order parameter */
    /* If either of the options is set we compute both */
    gmx::ArrayRef<const std::string> intfn = opt2fns("-o", NFILE, fnm);
    if (intfn.size() != 2)
    {
        gmx_fatal(FARGS, "No or not correct number (2) of output-files: %td", intfn.ssize());
    }
    calc_tetra_order_interface(
            ndxfnm, tpsfnm, trxfnm, binwidth, nsttblock, &frames, &xslices, &yslices, sg1, sg2, &intfpos, oenv);
    writesurftoxpms(intfpos, frames, xslices, yslices, binwidth, intfn, nlevels);

    if (bFourier)
    {
        gmx::ArrayRef<const std::string> spectra = opt2fns("-Spect", NFILE, fnm);
        if (spectra.size() != 2)
        {
            gmx_fatal(FARGS, "No or not correct number (2) of output-files: %td", spectra.ssize());
        }
        powerspectavg(intfpos, frames, xslices, yslices, spectra);
    }

    if (bRawOut)
    {
        gmx::ArrayRef<const std::string> raw = opt2fns("-or", NFILE, fnm);
        if (raw.size() != 2)
        {
            gmx_fatal(FARGS, "No or not correct number (2) of output-files: %td", raw.ssize());
        }
        writeraw(intfpos, frames, xslices, yslices, raw);
    }

    return 0;
}
