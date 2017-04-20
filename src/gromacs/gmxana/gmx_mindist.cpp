/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017, by the GROMACS development team, led by
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

#include <cmath>
#include <cstdlib>
#include <cstring>

#include <algorithm>

#include "gromacs/commandline/pargs.h"
#include "gromacs/commandline/viewit.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pbcutil/rmpbc.h"
#include "gromacs/topology/index.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"


static void periodic_dist(int ePBC,
                          matrix box, rvec x[], int n, int index[],
                          real *rmin, real *rmax, int *min_ind)
{
#define NSHIFT_MAX 26
    int  nsz, nshift, sx, sy, sz, i, j, s;
    real sqr_box, r2min, r2max, r2;
    rvec shift[NSHIFT_MAX], d0, d;

    sqr_box = std::min(norm2(box[XX]), norm2(box[YY]));
    if (ePBC == epbcXYZ)
    {
        sqr_box = std::min(sqr_box, norm2(box[ZZ]));
        nsz     = 1;
    }
    else if (ePBC == epbcXY)
    {
        nsz = 0;
    }
    else
    {
        gmx_fatal(FARGS, "pbc = %s is not supported by g_mindist",
                  epbc_names[ePBC]);
        nsz = 0; /* Keep compilers quiet */
    }

    nshift = 0;
    for (sz = -nsz; sz <= nsz; sz++)
    {
        for (sy = -1; sy <= 1; sy++)
        {
            for (sx = -1; sx <= 1; sx++)
            {
                if (sx != 0 || sy != 0 || sz != 0)
                {
                    for (i = 0; i < DIM; i++)
                    {
                        shift[nshift][i] =
                            sx*box[XX][i] + sy*box[YY][i] + sz*box[ZZ][i];
                    }
                    nshift++;
                }
            }
        }
    }

    r2min = sqr_box;
    r2max = 0;

    for (i = 0; i < n; i++)
    {
        for (j = i+1; j < n; j++)
        {
            rvec_sub(x[index[i]], x[index[j]], d0);
            r2 = norm2(d0);
            if (r2 > r2max)
            {
                r2max = r2;
            }
            for (s = 0; s < nshift; s++)
            {
                rvec_add(d0, shift[s], d);
                r2 = norm2(d);
                if (r2 < r2min)
                {
                    r2min      = r2;
                    min_ind[0] = i;
                    min_ind[1] = j;
                }
            }
        }
    }

    *rmin = std::sqrt(r2min);
    *rmax = std::sqrt(r2max);
}

static void periodic_mindist_plot(const char *trxfn, const char *outfn,
                                  const t_topology *top, int ePBC,
                                  int n, int index[], gmx_bool bSplit,
                                  const gmx_output_env_t *oenv)
{
    FILE        *out;
    const char  *leg[5] = { "min per.", "max int.", "box1", "box2", "box3" };
    t_trxstatus *status;
    real         t;
    rvec        *x;
    matrix       box;
    int          natoms, ind_min[2] = {0, 0}, ind_mini = 0, ind_minj = 0;
    real         rmin, rmax, rmint, tmint;
    gmx_bool     bFirst;
    gmx_rmpbc_t  gpbc = nullptr;

    natoms = read_first_x(oenv, &status, trxfn, &t, &x, box);

    check_index(nullptr, n, index, nullptr, natoms);

    out = xvgropen(outfn, "Minimum distance to periodic image",
                   output_env_get_time_label(oenv), "Distance (nm)", oenv);
    if (output_env_get_print_xvgr_codes(oenv))
    {
        fprintf(out, "@ subtitle \"and maximum internal distance\"\n");
    }
    xvgr_legend(out, 5, leg, oenv);

    rmint = box[XX][XX];
    tmint = 0;

    if (nullptr != top)
    {
        gpbc = gmx_rmpbc_init(&top->idef, ePBC, natoms);
    }

    bFirst = TRUE;
    do
    {
        if (nullptr != top)
        {
            gmx_rmpbc(gpbc, natoms, box, x);
        }

        periodic_dist(ePBC, box, x, n, index, &rmin, &rmax, ind_min);
        if (rmin < rmint)
        {
            rmint    = rmin;
            tmint    = t;
            ind_mini = ind_min[0];
            ind_minj = ind_min[1];
        }
        if (bSplit && !bFirst && std::abs(t/output_env_get_time_factor(oenv)) < 1e-5)
        {
            fprintf(out, "%s\n", output_env_get_print_xvgr_codes(oenv) ? "&" : "");
        }
        fprintf(out, "\t%g\t%6.3f %6.3f %6.3f %6.3f %6.3f\n",
                output_env_conv_time(oenv, t), rmin, rmax, norm(box[0]), norm(box[1]), norm(box[2]));
        bFirst = FALSE;
    }
    while (read_next_x(oenv, status, &t, x, box));

    if (nullptr != top)
    {
        gmx_rmpbc_done(gpbc);
    }

    xvgrclose(out);

    fprintf(stdout,
            "\nThe shortest periodic distance is %g (nm) at time %g (%s),\n"
            "between atoms %d and %d\n",
            rmint, output_env_conv_time(oenv, tmint), output_env_get_time_unit(oenv).c_str(),
            index[ind_mini]+1, index[ind_minj]+1);
}

static void calc_dist(real rcut, gmx_bool bPBC, int ePBC, matrix box, rvec x[],
                      int nx1, int nx2, int index1[], int index2[],
                      gmx_bool bGroup,
                      real *rmin, real *rmax, int *nmin, int *nmax,
                      int *ixmin, int *jxmin, int *ixmax, int *jxmax)
{
    int      i, j, i0 = 0, j1;
    int      ix, jx;
    int     *index3;
    rvec     dx;
    real     r2, rmin2, rmax2, rcut2;
    t_pbc    pbc;
    int      nmin_j, nmax_j;

    *ixmin = -1;
    *jxmin = -1;
    *ixmax = -1;
    *jxmax = -1;
    *nmin  = 0;
    *nmax  = 0;

    rcut2 = gmx::square(rcut);

    /* Must init pbc every step because of pressure coupling */
    if (bPBC)
    {
        set_pbc(&pbc, ePBC, box);
    }
    if (index2)
    {
        i0     = 0;
        j1     = nx2;
        index3 = index2;
    }
    else
    {
        j1     = nx1;
        index3 = index1;
    }

    rmin2 = 1e12;
    rmax2 = -1e12;

    for (j = 0; (j < j1); j++)
    {
        jx = index3[j];
        if (index2 == nullptr)
        {
            i0 = j + 1;
        }
        nmin_j = 0;
        nmax_j = 0;
        for (i = i0; (i < nx1); i++)
        {
            ix = index1[i];
            if (ix != jx)
            {
                if (bPBC)
                {
                    pbc_dx(&pbc, x[ix], x[jx], dx);
                }
                else
                {
                    rvec_sub(x[ix], x[jx], dx);
                }
                r2 = iprod(dx, dx);
                if (r2 < rmin2)
                {
                    rmin2  = r2;
                    *ixmin = ix;
                    *jxmin = jx;
                }
                if (r2 > rmax2)
                {
                    rmax2  = r2;
                    *ixmax = ix;
                    *jxmax = jx;
                }
                if (r2 <= rcut2)
                {
                    nmin_j++;
                }
                else
                {
                    nmax_j++;
                }
            }
        }
        if (bGroup)
        {
            if (nmin_j > 0)
            {
                (*nmin)++;
            }
            if (nmax_j > 0)
            {
                (*nmax)++;
            }
        }
        else
        {
            *nmin += nmin_j;
            *nmax += nmax_j;
        }
    }
    *rmin = std::sqrt(rmin2);
    *rmax = std::sqrt(rmax2);
}

void dist_plot(const char *fn, const char *afile, const char *dfile,
               const char *nfile, const char *rfile, const char *xfile,
               real rcut, gmx_bool bMat, const t_atoms *atoms,
               int ng, int *index[], int gnx[], char *grpn[], gmx_bool bSplit,
               gmx_bool bMin, int nres, int *residue, gmx_bool bPBC, int ePBC,
               gmx_bool bGroup, gmx_bool bEachResEachTime, gmx_bool bPrintResName,
               const gmx_output_env_t *oenv)
{
    FILE            *atm, *dist, *num;
    t_trxstatus     *trxout;
    char             buf[256];
    char           **leg;
    real             t, dmin, dmax, **mindres = nullptr, **maxdres = nullptr;
    int              nmin, nmax;
    t_trxstatus     *status;
    int              i = -1, j, k;
    int              min2, max2, min1r, min2r, max1r, max2r;
    int              min1 = 0;
    int              max1 = 0;
    int              oindex[2];
    rvec            *x0;
    matrix           box;
    gmx_bool         bFirst;
    FILE            *respertime = nullptr;

    if (read_first_x(oenv, &status, fn, &t, &x0, box) == 0)
    {
        gmx_fatal(FARGS, "Could not read coordinates from statusfile\n");
    }

    sprintf(buf, "%simum Distance", bMin ? "Min" : "Max");
    dist = xvgropen(dfile, buf, output_env_get_time_label(oenv), "Distance (nm)", oenv);
    sprintf(buf, "Number of Contacts %s %g nm", bMin ? "<" : ">", rcut);
    num    = nfile ? xvgropen(nfile, buf, output_env_get_time_label(oenv), "Number", oenv) : nullptr;
    atm    = afile ? gmx_ffopen(afile, "w") : nullptr;
    trxout = xfile ? open_trx(xfile, "w") : nullptr;

    if (bMat)
    {
        if (ng == 1)
        {
            snew(leg, 1);
            sprintf(buf, "Internal in %s", grpn[0]);
            leg[0] = gmx_strdup(buf);
            xvgr_legend(dist, 0, (const char**)leg, oenv);
            if (num)
            {
                xvgr_legend(num, 0, (const char**)leg, oenv);
            }
        }
        else
        {
            snew(leg, (ng*(ng-1))/2);
            for (i = j = 0; (i < ng-1); i++)
            {
                for (k = i+1; (k < ng); k++, j++)
                {
                    sprintf(buf, "%s-%s", grpn[i], grpn[k]);
                    leg[j] = gmx_strdup(buf);
                }
            }
            xvgr_legend(dist, j, (const char**)leg, oenv);
            if (num)
            {
                xvgr_legend(num, j, (const char**)leg, oenv);
            }
        }
    }
    else
    {
        snew(leg, ng-1);
        for (i = 0; (i < ng-1); i++)
        {
            sprintf(buf, "%s-%s", grpn[0], grpn[i+1]);
            leg[i] = gmx_strdup(buf);
        }
        xvgr_legend(dist, ng-1, (const char**)leg, oenv);
        if (num)
        {
            xvgr_legend(num, ng-1, (const char**)leg, oenv);
        }
    }

    if (bEachResEachTime)
    {
        sprintf(buf, "%simum Distance", bMin ? "Min" : "Max");
        respertime = xvgropen(rfile, buf, output_env_get_time_label(oenv), "Distance (nm)", oenv);
        xvgr_legend(respertime, ng-1, (const char**)leg, oenv);
        if (bPrintResName && output_env_get_print_xvgr_codes(oenv) )
        {
            fprintf(respertime, "# ");

            for (j = 0; j < nres; j++)
            {
                fprintf(respertime, "%s%d ", *(atoms->resinfo[atoms->atom[index[0][residue[j]]].resind].name), atoms->atom[index[0][residue[j]]].resind);
            }
            fprintf(respertime, "\n");
        }

    }

    if (nres)
    {
        snew(mindres, ng-1);
        snew(maxdres, ng-1);
        for (i = 1; i < ng; i++)
        {
            snew(mindres[i-1], nres);
            snew(maxdres[i-1], nres);
            for (j = 0; j < nres; j++)
            {
                mindres[i-1][j] = 1e6;
            }
            /* maxdres[*][*] is already 0 */
        }
    }
    bFirst = TRUE;
    do
    {
        if (bSplit && !bFirst && std::abs(t/output_env_get_time_factor(oenv)) < 1e-5)
        {
            fprintf(dist, "%s\n", output_env_get_print_xvgr_codes(oenv) ? "&" : "");
            if (num)
            {
                fprintf(num, "%s\n", output_env_get_print_xvgr_codes(oenv) ? "&" : "");
            }
            if (atm)
            {
                fprintf(atm, "%s\n", output_env_get_print_xvgr_codes(oenv) ? "&" : "");
            }
        }
        fprintf(dist, "%12e", output_env_conv_time(oenv, t));
        if (num)
        {
            fprintf(num, "%12e", output_env_conv_time(oenv, t));
        }

        if (bMat)
        {
            if (ng == 1)
            {
                calc_dist(rcut, bPBC, ePBC, box, x0, gnx[0], gnx[0], index[0], index[0], bGroup,
                          &dmin, &dmax, &nmin, &nmax, &min1, &min2, &max1, &max2);
                fprintf(dist, "  %12e", bMin ? dmin : dmax);
                if (num)
                {
                    fprintf(num, "  %8d", bMin ? nmin : nmax);
                }
            }
            else
            {
                for (i = 0; (i < ng-1); i++)
                {
                    for (k = i+1; (k < ng); k++)
                    {
                        calc_dist(rcut, bPBC, ePBC, box, x0, gnx[i], gnx[k], index[i], index[k],
                                  bGroup, &dmin, &dmax, &nmin, &nmax, &min1, &min2, &max1, &max2);
                        fprintf(dist, "  %12e", bMin ? dmin : dmax);
                        if (num)
                        {
                            fprintf(num, "  %8d", bMin ? nmin : nmax);
                        }
                    }
                }
            }
        }
        else
        {
            for (i = 1; (i < ng); i++)
            {
                calc_dist(rcut, bPBC, ePBC, box, x0, gnx[0], gnx[i], index[0], index[i], bGroup,
                          &dmin, &dmax, &nmin, &nmax, &min1, &min2, &max1, &max2);
                fprintf(dist, "  %12e", bMin ? dmin : dmax);
                if (num)
                {
                    fprintf(num, "  %8d", bMin ? nmin : nmax);
                }
                if (nres)
                {
                    for (j = 0; j < nres; j++)
                    {
                        calc_dist(rcut, bPBC, ePBC, box, x0, residue[j+1]-residue[j], gnx[i],
                                  &(index[0][residue[j]]), index[i], bGroup,
                                  &dmin, &dmax, &nmin, &nmax, &min1r, &min2r, &max1r, &max2r);
                        mindres[i-1][j] = std::min(mindres[i-1][j], dmin);
                        maxdres[i-1][j] = std::max(maxdres[i-1][j], dmax);
                    }
                }
            }
        }
        fprintf(dist, "\n");
        if (num)
        {
            fprintf(num, "\n");
        }
        if ( (bMin ? min1 : max1) != -1)
        {
            if (atm)
            {
                fprintf(atm, "%12e  %12d  %12d\n",
                        output_env_conv_time(oenv, t), 1+(bMin ? min1 : max1),
                        1+(bMin ? min2 : max2));
            }
        }

        if (trxout)
        {
            oindex[0] = bMin ? min1 : max1;
            oindex[1] = bMin ? min2 : max2;
            write_trx(trxout, 2, oindex, atoms, i, t, box, x0, nullptr, nullptr);
        }
        bFirst = FALSE;
        /*dmin should be minimum distance for residue and group*/
        if (bEachResEachTime)
        {
            fprintf(respertime, "%12e", t);
            for (i = 1; i < ng; i++)
            {
                for (j = 0; j < nres; j++)
                {
                    fprintf(respertime, " %7g", bMin ? mindres[i-1][j] : maxdres[i-1][j]);
                    /*reset distances for next time point*/
                    mindres[i-1][j] = 1e6;
                    maxdres[i-1][j] = 0;
                }
            }
            fprintf(respertime, "\n");
        }
    }
    while (read_next_x(oenv, status, &t, x0, box));

    close_trj(status);
    xvgrclose(dist);
    if (num)
    {
        xvgrclose(num);
    }
    if (atm)
    {
        gmx_ffclose(atm);
    }
    if (trxout)
    {
        close_trx(trxout);
    }
    if (respertime)
    {
        xvgrclose(respertime);
    }

    if (nres && !bEachResEachTime)
    {
        FILE *res;

        sprintf(buf, "%simum Distance", bMin ? "Min" : "Max");
        res = xvgropen(rfile, buf, "Residue (#)", "Distance (nm)", oenv);
        xvgr_legend(res, ng-1, (const char**)leg, oenv);
        for (j = 0; j < nres; j++)
        {
            fprintf(res, "%4d", j+1);
            for (i = 1; i < ng; i++)
            {
                fprintf(res, " %7g", bMin ? mindres[i-1][j] : maxdres[i-1][j]);
            }
            fprintf(res, "\n");
        }
        xvgrclose(res);
    }

    sfree(x0);
}

int find_residues(const t_atoms *atoms, int n, int index[], int **resindex)
{
    int  i;
    int  nres      = 0, resnr, presnr = 0;
    bool presFound = false;
    int *residx;

    /* build index of first atom numbers for each residue */
    snew(residx, atoms->nres+1);
    for (i = 0; i < n; i++)
    {
        resnr = atoms->atom[index[i]].resind;
        if (!presFound || resnr != presnr)
        {
            residx[nres] = i;
            nres++;
            presnr    = resnr;
            presFound = true;
        }
    }
    if (debug)
    {
        printf("Found %d residues out of %d (%d/%d atoms)\n",
               nres, atoms->nres, atoms->nr, n);
    }
    srenew(residx, nres+1);
    /* mark end of last residue */
    residx[nres] = n;
    *resindex    = residx;
    return nres;
}

void dump_res(FILE *out, int nres, int *resindex, int index[])
{
    int i, j;

    for (i = 0; i < nres-1; i++)
    {
        fprintf(out, "Res %d (%d):", i, resindex[i+1]-resindex[i]);
        for (j = resindex[i]; j < resindex[i+1]; j++)
        {
            fprintf(out, " %d(%d)", j, index[j]);
        }
        fprintf(out, "\n");
    }
}

int gmx_mindist(int argc, char *argv[])
{
    const char       *desc[] = {
        "[THISMODULE] computes the distance between one group and a number of",
        "other groups. Both the minimum distance",
        "(between any pair of atoms from the respective groups)",
        "and the number of contacts within a given",
        "distance are written to two separate output files.",
        "With the [TT]-group[tt] option a contact of an atom in another group",
        "with multiple atoms in the first group is counted as one contact",
        "instead of as multiple contacts.",
        "With [TT]-or[tt], minimum distances to each residue in the first",
        "group are determined and plotted as a function of residue number.[PAR]",
        "With option [TT]-pi[tt] the minimum distance of a group to its",
        "periodic image is plotted. This is useful for checking if a protein",
        "has seen its periodic image during a simulation. Only one shift in",
        "each direction is considered, giving a total of 26 shifts.",
        "It also plots the maximum distance within the group and the lengths",
        "of the three box vectors.[PAR]",
        "Also [gmx-distance] and [gmx-pairdist] calculate distances."
    };

    static gmx_bool   bMat             = FALSE, bPI = FALSE, bSplit = FALSE, bMax = FALSE, bPBC = TRUE;
    static gmx_bool   bGroup           = FALSE;
    static real       rcutoff          = 0.6;
    static int        ng               = 1;
    static gmx_bool   bEachResEachTime = FALSE, bPrintResName = FALSE;
    t_pargs           pa[]             = {
        { "-matrix", FALSE, etBOOL, {&bMat},
          "Calculate half a matrix of group-group distances" },
        { "-max",    FALSE, etBOOL, {&bMax},
          "Calculate *maximum* distance instead of minimum" },
        { "-d",      FALSE, etREAL, {&rcutoff},
          "Distance for contacts" },
        { "-group",      FALSE, etBOOL, {&bGroup},
          "Count contacts with multiple atoms in the first group as one" },
        { "-pi",     FALSE, etBOOL, {&bPI},
          "Calculate minimum distance with periodic images" },
        { "-split",  FALSE, etBOOL, {&bSplit},
          "Split graph where time is zero" },
        { "-ng",       FALSE, etINT, {&ng},
          "Number of secondary groups to compute distance to a central group" },
        { "-pbc",    FALSE, etBOOL, {&bPBC},
          "Take periodic boundary conditions into account" },
        { "-respertime",  FALSE, etBOOL, {&bEachResEachTime},
          "When writing per-residue distances, write distance for each time point" },
        { "-printresname",  FALSE, etBOOL, {&bPrintResName},
          "Write residue names" }
    };
    gmx_output_env_t *oenv;
    t_topology       *top  = nullptr;
    int               ePBC = -1;
    rvec             *x;
    matrix            box;
    gmx_bool          bTop = FALSE;

    int               i, nres = 0;
    const char       *trxfnm, *tpsfnm, *ndxfnm, *distfnm, *numfnm, *atmfnm, *oxfnm, *resfnm;
    char            **grpname;
    int              *gnx;
    int             **index, *residues = nullptr;
    t_filenm          fnm[] = {
        { efTRX, "-f",  nullptr,      ffREAD },
        { efTPS,  nullptr, nullptr,      ffOPTRD },
        { efNDX,  nullptr, nullptr,      ffOPTRD },
        { efXVG, "-od", "mindist",  ffWRITE },
        { efXVG, "-on", "numcont",  ffOPTWR },
        { efOUT, "-o", "atm-pair", ffOPTWR },
        { efTRO, "-ox", "mindist",  ffOPTWR },
        { efXVG, "-or", "mindistres", ffOPTWR }
    };
#define NFILE asize(fnm)

    if (!parse_common_args(&argc, argv,
                           PCA_CAN_VIEW | PCA_CAN_TIME | PCA_TIME_UNIT,
                           NFILE, fnm, asize(pa), pa, asize(desc), desc, 0, nullptr, &oenv))
    {
        return 0;
    }

    trxfnm  = ftp2fn(efTRX, NFILE, fnm);
    ndxfnm  = ftp2fn_null(efNDX, NFILE, fnm);
    distfnm = opt2fn("-od", NFILE, fnm);
    numfnm  = opt2fn_null("-on", NFILE, fnm);
    atmfnm  = ftp2fn_null(efOUT, NFILE, fnm);
    oxfnm   = opt2fn_null("-ox", NFILE, fnm);
    resfnm  = opt2fn_null("-or", NFILE, fnm);
    if (bPI || resfnm != nullptr)
    {
        /* We need a tps file */
        tpsfnm = ftp2fn(efTPS, NFILE, fnm);
    }
    else
    {
        tpsfnm = ftp2fn_null(efTPS, NFILE, fnm);
    }

    if (!tpsfnm && !ndxfnm)
    {
        gmx_fatal(FARGS, "You have to specify either the index file or a tpr file");
    }

    if (bPI)
    {
        ng = 1;
        fprintf(stderr, "Choose a group for distance calculation\n");
    }
    else if (!bMat)
    {
        ng++;
    }

    snew(gnx, ng);
    snew(index, ng);
    snew(grpname, ng);

    if (tpsfnm || resfnm || !ndxfnm)
    {
        snew(top, 1);
        bTop = read_tps_conf(tpsfnm, top, &ePBC, &x, nullptr, box, FALSE);
        if (bPI && !bTop)
        {
            printf("\nWARNING: Without a run input file a trajectory with broken molecules will not give the correct periodic image distance\n\n");
        }
    }
    get_index(top ? &(top->atoms) : nullptr, ndxfnm, ng, gnx, index, grpname);

    if (bMat && (ng == 1))
    {
        ng = gnx[0];
        printf("Special case: making distance matrix between all atoms in group %s\n",
               grpname[0]);
        srenew(gnx, ng);
        srenew(index, ng);
        srenew(grpname, ng);
        for (i = 1; (i < ng); i++)
        {
            gnx[i]      = 1;
            grpname[i]  = grpname[0];
            snew(index[i], 1);
            index[i][0] = index[0][i];
        }
        gnx[0] = 1;
    }

    if (resfnm)
    {
        GMX_RELEASE_ASSERT(top != nullptr, "top pointer cannot be NULL when finding residues");
        nres = find_residues(&(top->atoms), gnx[0], index[0], &residues);

        if (debug)
        {
            dump_res(debug, nres, residues, index[0]);
        }
    }

    if (bPI)
    {
        periodic_mindist_plot(trxfnm, distfnm, top, ePBC, gnx[0], index[0], bSplit, oenv);
    }
    else
    {
        dist_plot(trxfnm, atmfnm, distfnm, numfnm, resfnm, oxfnm,
                  rcutoff, bMat, top ? &(top->atoms) : nullptr,
                  ng, index, gnx, grpname, bSplit, !bMax, nres, residues, bPBC, ePBC,
                  bGroup, bEachResEachTime, bPrintResName, oenv);
    }

    do_view(oenv, distfnm, "-nxy");
    if (!bPI)
    {
        do_view(oenv, numfnm, "-nxy");
    }

    return 0;
}
