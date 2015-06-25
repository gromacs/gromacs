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

#include <math.h>

#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/matio.h"
#include "gromacs/fileio/strdb.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/legacyheaders/names.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/legacyheaders/viewit.h"
#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/index.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"


static void calc_dist(int nind, atom_id index[], rvec x[], int ePBC, matrix box,
                      real **d)
{
    int      i, j;
    real    *xi;
    rvec     dx;
    real     temp2;
    t_pbc    pbc;

    set_pbc(&pbc, ePBC, box);
    for (i = 0; (i < nind-1); i++)
    {
        xi = x[index[i]];
        for (j = i+1; (j < nind); j++)
        {
            pbc_dx(&pbc, xi, x[index[j]], dx);
            temp2   = norm2(dx);
            d[i][j] = sqrt(temp2);
        }
    }
}

static void calc_dist_tot(int nind, atom_id index[], rvec x[],
                          int ePBC, matrix box,
                          real **d, real **dtot, real **dtot2,
                          gmx_bool bNMR, real **dtot1_3, real **dtot1_6)
{
    int      i, j;
    real    *xi;
    real     temp, temp2, temp1_3;
    rvec     dx;
    t_pbc    pbc;

    set_pbc(&pbc, ePBC, box);
    for (i = 0; (i < nind-1); i++)
    {
        xi = x[index[i]];
        for (j = i+1; (j < nind); j++)
        {
            pbc_dx(&pbc, xi, x[index[j]], dx);
            temp2        = norm2(dx);
            temp         = sqrt(temp2);
            d[i][j]      = temp;
            dtot[i][j]  += temp;
            dtot2[i][j] += temp2;
            if (bNMR)
            {
                temp1_3        = 1.0/(temp*temp2);
                dtot1_3[i][j] += temp1_3;
                dtot1_6[i][j] += temp1_3*temp1_3;
            }
        }
    }
}

static void calc_nmr(int nind, int nframes, real **dtot1_3, real **dtot1_6,
                     real *max1_3, real *max1_6)
{
    int     i, j;
    real    temp1_3, temp1_6;

    for (i = 0; (i < nind-1); i++)
    {
        for (j = i+1; (j < nind); j++)
        {
            temp1_3 = pow(dtot1_3[i][j]/nframes, -1.0/3.0);
            temp1_6 = pow(dtot1_6[i][j]/nframes, -1.0/6.0);
            if (temp1_3 > *max1_3)
            {
                *max1_3 = temp1_3;
            }
            if (temp1_6 > *max1_6)
            {
                *max1_6 = temp1_6;
            }
            dtot1_3[i][j] = temp1_3;
            dtot1_6[i][j] = temp1_6;
            dtot1_3[j][i] = temp1_3;
            dtot1_6[j][i] = temp1_6;
        }
    }
}

static char Hnum[] = "123";

typedef struct {
    int  nr;
    real r_3;
    real r_6;
    real i_3;
    real i_6;
} t_noe;

typedef struct {
    int   anr;
    int   ianr;
    int   rnr;
    char *aname;
    char *rname;
} t_noe_gr;

typedef struct {
    int   rnr;
    char *nname;
    char *rname;
    char *aname;
} t_equiv;

static int read_equiv(const char *eq_fn, t_equiv ***equivptr)
{
    FILE     *fp;
    char      line[STRLEN], resname[10], atomname[10], *lp;
    int       neq, na, n, resnr;
    t_equiv **equiv;

    fp    = gmx_ffopen(eq_fn, "r");
    neq   = 0;
    equiv = NULL;
    while (get_a_line(fp, line, STRLEN))
    {
        lp = line;
        /* this is not efficient, but I'm lazy */
        srenew(equiv, neq+1);
        equiv[neq] = NULL;
        na         = 0;
        if (sscanf(lp, "%s %n", atomname, &n) == 1)
        {
            lp += n;
            snew(equiv[neq], 1);
            equiv[neq][0].nname = gmx_strdup(atomname);
            while (sscanf(lp, "%d %s %s %n", &resnr, resname, atomname, &n) == 3)
            {
                /* this is not efficient, but I'm lazy (again) */
                srenew(equiv[neq], na+1);
                equiv[neq][na].rnr   = resnr-1;
                equiv[neq][na].rname = gmx_strdup(resname);
                equiv[neq][na].aname = gmx_strdup(atomname);
                if (na > 0)
                {
                    equiv[neq][na].nname = NULL;
                }
                na++;
                lp += n;
            }
        }
        /* make empty element as flag for end of array */
        srenew(equiv[neq], na+1);
        equiv[neq][na].rnr   = NOTSET;
        equiv[neq][na].rname = NULL;
        equiv[neq][na].aname = NULL;

        /* next */
        neq++;
    }
    gmx_ffclose(fp);

    *equivptr = equiv;

    return neq;
}

static void dump_equiv(FILE *out, int neq, t_equiv **equiv)
{
    int i, j;

    fprintf(out, "Dumping equivalent list\n");
    for (i = 0; i < neq; i++)
    {
        fprintf(out, "%s", equiv[i][0].nname);
        for (j = 0; equiv[i][j].rnr != NOTSET; j++)
        {
            fprintf(out, " %d %s %s",
                    equiv[i][j].rnr, equiv[i][j].rname, equiv[i][j].aname);
        }
        fprintf(out, "\n");
    }
}

static gmx_bool is_equiv(int neq, t_equiv **equiv, char **nname,
                         int rnr1, char *rname1, char *aname1,
                         int rnr2, char *rname2, char *aname2)
{
    int      i, j;
    gmx_bool bFound;

    bFound = FALSE;
    /* we can terminate each loop when bFound is true! */
    for (i = 0; i < neq && !bFound; i++)
    {
        /* find first atom */
        for (j = 0; equiv[i][j].rnr != NOTSET && !bFound; j++)
        {
            bFound = ( equiv[i][j].rnr == rnr1 &&
                       strcmp(equiv[i][j].rname, rname1) == 0 &&
                       strcmp(equiv[i][j].aname, aname1) == 0 );
        }
        if (bFound)
        {
            /* find second atom */
            bFound = FALSE;
            for (j = 0; equiv[i][j].rnr != NOTSET && !bFound; j++)
            {
                bFound = ( equiv[i][j].rnr == rnr2 &&
                           strcmp(equiv[i][j].rname, rname2) == 0 &&
                           strcmp(equiv[i][j].aname, aname2) == 0 );
            }
        }
    }
    if (bFound)
    {
        *nname = gmx_strdup(equiv[i-1][0].nname);
    }

    return bFound;
}

static int analyze_noe_equivalent(const char *eq_fn,
                                  t_atoms *atoms, int isize, atom_id *index,
                                  gmx_bool bSumH,
                                  atom_id *noe_index, t_noe_gr *noe_gr)
{
    FILE     *fp;
    int       i, j, anmil, anmjl, rnri, rnrj, gi, groupnr, neq;
    char     *anmi, *anmj, **nnm;
    gmx_bool  bMatch, bEquiv;
    t_equiv **equiv;

    snew(nnm, isize);
    if (bSumH)
    {
        if (eq_fn)
        {
            neq = read_equiv(eq_fn, &equiv);
            if (debug)
            {
                dump_equiv(debug, neq, equiv);
            }
        }
        else
        {
            neq   = 0;
            equiv = NULL;
        }

        groupnr = 0;
        for (i = 0; i < isize; i++)
        {
            if (equiv && i < isize-1)
            {
                /* check explicit list of equivalent atoms */
                do
                {
                    j      = i+1;
                    rnri   = atoms->atom[index[i]].resind;
                    rnrj   = atoms->atom[index[j]].resind;
                    bEquiv =
                        is_equiv(neq, equiv, &nnm[i],
                                 rnri, *atoms->resinfo[rnri].name, *atoms->atomname[index[i]],
                                 rnrj, *atoms->resinfo[rnrj].name, *atoms->atomname[index[j]]);
                    if (nnm[i] && bEquiv)
                    {
                        nnm[j] = gmx_strdup(nnm[i]);
                    }
                    if (bEquiv)
                    {
                        /* set index for matching atom */
                        noe_index[i] = groupnr;
                        /* skip matching atom */
                        i = j;
                    }
                }
                while (bEquiv && i < isize-1);
            }
            else
            {
                bEquiv = FALSE;
            }
            if (!bEquiv)
            {
                /* look for triplets of consecutive atoms with name XX?,
                   X are any number of letters or digits and ? goes from 1 to 3
                   This is supposed to cover all CH3 groups and the like */
                anmi   = *atoms->atomname[index[i]];
                anmil  = strlen(anmi);
                bMatch = i <= isize-3 && anmi[anmil-1] == '1';
                if (bMatch)
                {
                    for (j = 1; j < 3; j++)
                    {
                        anmj   = *atoms->atomname[index[i+j]];
                        anmjl  = strlen(anmj);
                        bMatch = bMatch && ( anmil == anmjl && anmj[anmjl-1] == Hnum[j] &&
                                             strncmp(anmi, anmj, anmil-1) == 0 );
                    }
                }
                /* set index for this atom */
                noe_index[i] = groupnr;
                if (bMatch)
                {
                    /* set index for next two matching atoms */
                    for (j = 1; j < 3; j++)
                    {
                        noe_index[i+j] = groupnr;
                    }
                    /* skip matching atoms */
                    i += 2;
                }
            }
            groupnr++;
        }
    }
    else
    {
        /* make index without looking for equivalent atoms */
        for (i = 0; i < isize; i++)
        {
            noe_index[i] = i;
        }
        groupnr = isize;
    }
    noe_index[isize] = groupnr;

    if (debug)
    {
        /* dump new names */
        for (i = 0; i < isize; i++)
        {
            rnri = atoms->atom[index[i]].resind;
            fprintf(debug, "%s %s %d -> %s\n", *atoms->atomname[index[i]],
                    *atoms->resinfo[rnri].name, rnri, nnm[i] ? nnm[i] : "");
        }
    }

    for (i = 0; i < isize; i++)
    {
        gi = noe_index[i];
        if (!noe_gr[gi].aname)
        {
            noe_gr[gi].ianr = i;
            noe_gr[gi].anr  = index[i];
            if (nnm[i])
            {
                noe_gr[gi].aname = gmx_strdup(nnm[i]);
            }
            else
            {
                noe_gr[gi].aname = gmx_strdup(*atoms->atomname[index[i]]);
                if (noe_index[i] == noe_index[i+1])
                {
                    noe_gr[gi].aname[strlen(noe_gr[gi].aname)-1] = '*';
                }
            }
            noe_gr[gi].rnr   = atoms->atom[index[i]].resind;
            noe_gr[gi].rname = gmx_strdup(*atoms->resinfo[noe_gr[gi].rnr].name);
            /* dump group definitions */
            if (debug)
            {
                fprintf(debug, "%d %d %d %d %s %s %d\n", i, gi,
                        noe_gr[gi].ianr, noe_gr[gi].anr, noe_gr[gi].aname,
                        noe_gr[gi].rname, noe_gr[gi].rnr);
            }
        }
    }
    for (i = 0; i < isize; i++)
    {
        sfree(nnm[i]);
    }
    sfree(nnm);

    return groupnr;
}

/* #define NSCALE 3 */
/*  static char *noe_scale[NSCALE+1] = { */
/*    "strong", "medium", "weak", "none" */
/*  }; */
#define NSCALE 6

static char *noe2scale(real r3, real r6, real rmax)
{
    int         i, s3, s6;
    static char buf[NSCALE+1];

    /* r goes from 0 to rmax
       NSCALE*r/rmax goes from 0 to NSCALE
       NSCALE - NSCALE*r/rmax goes from NSCALE to 0 */
    s3 = NSCALE - min(NSCALE, (int)(NSCALE*r3/rmax));
    s6 = NSCALE - min(NSCALE, (int)(NSCALE*r6/rmax));

    for (i = 0; i < s3; i++)
    {
        buf[i] = '=';
    }
    for (; i < s6; i++)
    {
        buf[i] = '-';
    }
    buf[i] = '\0';

    return buf;
}

static void calc_noe(int isize, atom_id *noe_index,
                     real **dtot1_3, real **dtot1_6, int gnr, t_noe **noe)
{
    int i, j, gi, gj;

    /* make half matrix of noe-group distances from atom distances */
    for (i = 0; i < isize; i++)
    {
        gi = noe_index[i];
        for (j = i; j < isize; j++)
        {
            gj = noe_index[j];
            noe[gi][gj].nr++;
            noe[gi][gj].i_3 += pow(dtot1_3[i][j], -3);
            noe[gi][gj].i_6 += pow(dtot1_6[i][j], -6);
        }
    }

    /* make averages */
    for (i = 0; i < gnr; i++)
    {
        for (j = i+1; j < gnr; j++)
        {
            noe[i][j].r_3 = pow(noe[i][j].i_3/noe[i][j].nr, -1.0/3.0);
            noe[i][j].r_6 = pow(noe[i][j].i_6/noe[i][j].nr, -1.0/6.0);
            noe[j][i]     = noe[i][j];
        }
    }
}

static void write_noe(FILE *fp, int gnr, t_noe **noe, t_noe_gr *noe_gr, real rmax)
{
    int      i, j;
    real     r3, r6, min3, min6;
    char     buf[10], b3[10], b6[10];
    t_noe_gr gri, grj;

    min3 = min6 = 1e6;
    fprintf(fp,
            ";%4s %3s %4s %4s%3s %4s %4s %4s %4s%3s %5s %5s %8s %2s %2s %s\n",
            "ianr", "anr", "anm", "rnm", "rnr", "ianr", "anr", "anm", "rnm", "rnr",
            "1/r^3", "1/r^6", "intnsty", "Dr", "Da", "scale");
    for (i = 0; i < gnr; i++)
    {
        gri = noe_gr[i];
        for (j = i+1; j < gnr; j++)
        {
            grj  = noe_gr[j];
            r3   = noe[i][j].r_3;
            r6   = noe[i][j].r_6;
            min3 = min(r3, min3);
            min6 = min(r6, min6);
            if (r3 < rmax || r6 < rmax)
            {
                if (grj.rnr == gri.rnr)
                {
                    sprintf(buf, "%2d", grj.anr-gri.anr);
                }
                else
                {
                    buf[0] = '\0';
                }
                if (r3 < rmax)
                {
                    sprintf(b3, "%-5.3f", r3);
                }
                else
                {
                    strcpy(b3, "-");
                }
                if (r6 < rmax)
                {
                    sprintf(b6, "%-5.3f", r6);
                }
                else
                {
                    strcpy(b6, "-");
                }
                fprintf(fp,
                        "%4d %4d %4s %4s%3d %4d %4d %4s %4s%3d %5s %5s %8d %2d %2s %s\n",
                        gri.ianr+1, gri.anr+1, gri.aname, gri.rname, gri.rnr+1,
                        grj.ianr+1, grj.anr+1, grj.aname, grj.rname, grj.rnr+1,
                        b3, b6, (int)(noe[i][j].i_6+0.5), grj.rnr-gri.rnr, buf,
                        noe2scale(r3, r6, rmax));
            }
        }
    }
#define MINI ((i == 3) ? min3 : min6)
    for (i = 3; i <= 6; i += 3)
    {
        if (MINI > rmax)
        {
            fprintf(stdout, "NOTE: no 1/r^%d averaged distances found below %g, "
                    "smallest was %g\n", i, rmax, MINI );
        }
        else
        {
            fprintf(stdout, "Smallest 1/r^%d averaged distance was %g\n", i, MINI );
        }
    }
#undef MINI
}

static void calc_rms(int nind, int nframes,
                     real **dtot, real **dtot2,
                     real **rmsmat,  real *rmsmax,
                     real **rmscmat, real *rmscmax,
                     real **meanmat, real *meanmax)
{
    int     i, j;
    real    mean, mean2, rms, rmsc;
/* N.B. dtot and dtot2 contain the total distance and the total squared
 * distance respectively, BUT they return RMS and the scaled RMS resp.
 */

    *rmsmax  = -1000;
    *rmscmax = -1000;
    *meanmax = -1000;

    for (i = 0; (i < nind-1); i++)
    {
        for (j = i+1; (j < nind); j++)
        {
            mean  = dtot[i][j]/nframes;
            mean2 = dtot2[i][j]/nframes;
            rms   = sqrt(max(0, mean2-mean*mean));
            rmsc  = rms/mean;
            if (mean > *meanmax)
            {
                *meanmax = mean;
            }
            if (rms  > *rmsmax)
            {
                *rmsmax = rms;
            }
            if (rmsc > *rmscmax)
            {
                *rmscmax = rmsc;
            }
            meanmat[i][j] = meanmat[j][i] = mean;
            rmsmat[i][j]  = rmsmat[j][i] = rms;
            rmscmat[i][j] = rmscmat[j][i] = rmsc;
        }
    }
}

real rms_diff(int natom, real **d, real **d_r)
{
    int  i, j;
    real r, r2;

    r2 = 0.0;
    for (i = 0; (i < natom-1); i++)
    {
        for (j = i+1; (j < natom); j++)
        {
            r   = d[i][j]-d_r[i][j];
            r2 += r*r;
        }
    }
    r2 /= (natom*(natom-1))/2;

    return sqrt(r2);
}

int gmx_rmsdist(int argc, char *argv[])
{
    const char     *desc[] = {
        "[THISMODULE] computes the root mean square deviation of atom distances,",
        "which has the advantage that no fit is needed like in standard RMS",
        "deviation as computed by [gmx-rms].",
        "The reference structure is taken from the structure file.",
        "The RMSD at time t is calculated as the RMS",
        "of the differences in distance between atom-pairs in the reference",
        "structure and the structure at time t.[PAR]",
        "[THISMODULE] can also produce matrices of the rms distances, rms distances",
        "scaled with the mean distance and the mean distances and matrices with",
        "NMR averaged distances (1/r^3 and 1/r^6 averaging). Finally, lists",
        "of atom pairs with 1/r^3 and 1/r^6 averaged distance below the",
        "maximum distance ([TT]-max[tt], which will default to 0.6 in this case)",
        "can be generated, by default averaging over equivalent hydrogens",
        "(all triplets of hydrogens named \\*[123]). Additionally a list of",
        "equivalent atoms can be supplied ([TT]-equiv[tt]), each line containing",
        "a set of equivalent atoms specified as residue number and name and",
        "atom name; e.g.:[PAR]",
        "[TT]HB* 3 SER  HB1 3 SER  HB2[tt][PAR]",
        "Residue and atom names must exactly match those in the structure",
        "file, including case. Specifying non-sequential atoms is undefined."

    };

    int             natom, i, j, teller, gi, gj;
    real            t;

    t_topology      top;
    int             ePBC;
    t_atoms        *atoms;
    matrix          box;
    rvec           *x;
    FILE           *fp;

    t_trxstatus    *status;
    int             isize, gnr = 0;
    atom_id        *index, *noe_index;
    char           *grpname;
    real          **d_r, **d, **dtot, **dtot2, **mean, **rms, **rmsc, *resnr;
    real          **dtot1_3 = NULL, **dtot1_6 = NULL;
    real            rmsnow, meanmax, rmsmax, rmscmax;
    real            max1_3, max1_6;
    t_noe_gr       *noe_gr = NULL;
    t_noe         **noe    = NULL;
    t_rgb           rlo, rhi;
    char            buf[255];
    gmx_bool        bRMS, bScale, bMean, bNOE, bNMR3, bNMR6, bNMR;

    static int      nlevels  = 40;
    static real     scalemax = -1.0;
    static gmx_bool bSumH    = TRUE;
    static gmx_bool bPBC     = TRUE;
    output_env_t    oenv;

    t_pargs         pa[] = {
        { "-nlevels",   FALSE, etINT,  {&nlevels},
          "Discretize RMS in this number of levels" },
        { "-max",   FALSE, etREAL, {&scalemax},
          "Maximum level in matrices" },
        { "-sumh",  FALSE, etBOOL, {&bSumH},
          "Average distance over equivalent hydrogens" },
        { "-pbc",   FALSE, etBOOL, {&bPBC},
          "Use periodic boundary conditions when computing distances" }
    };
    t_filenm        fnm[] = {
        { efTRX, "-f",   NULL,       ffREAD },
        { efTPS, NULL,   NULL,       ffREAD },
        { efNDX, NULL,   NULL,       ffOPTRD },
        { efDAT, "-equiv", "equiv",   ffOPTRD },
        { efXVG, NULL,   "distrmsd", ffWRITE },
        { efXPM, "-rms", "rmsdist",  ffOPTWR },
        { efXPM, "-scl", "rmsscale", ffOPTWR },
        { efXPM, "-mean", "rmsmean",  ffOPTWR },
        { efXPM, "-nmr3", "nmr3",     ffOPTWR },
        { efXPM, "-nmr6", "nmr6",     ffOPTWR },
        { efDAT, "-noe", "noe",     ffOPTWR },
    };
#define NFILE asize(fnm)

    if (!parse_common_args(&argc, argv, PCA_CAN_VIEW | PCA_CAN_TIME,
                           NFILE, fnm, asize(pa), pa, asize(desc), desc, 0, NULL, &oenv))
    {
        return 0;
    }

    bRMS   = opt2bSet("-rms", NFILE, fnm);
    bScale = opt2bSet("-scl", NFILE, fnm);
    bMean  = opt2bSet("-mean", NFILE, fnm);
    bNOE   = opt2bSet("-noe", NFILE, fnm);
    bNMR3  = opt2bSet("-nmr3", NFILE, fnm);
    bNMR6  = opt2bSet("-nmr6", NFILE, fnm);
    bNMR   = bNMR3 || bNMR6 || bNOE;

    max1_3 = 0;
    max1_6 = 0;

    /* check input */
    if (bNOE && scalemax < 0)
    {
        scalemax = 0.6;
        fprintf(stderr, "WARNING: using -noe without -max "
                "makes no sense, setting -max to %g\n\n", scalemax);
    }

    /* get topology and index */
    read_tps_conf(ftp2fn(efTPS, NFILE, fnm), buf, &top, &ePBC, &x, NULL, box, FALSE);

    if (!bPBC)
    {
        ePBC = epbcNONE;
    }
    atoms = &(top.atoms);

    get_index(atoms, ftp2fn_null(efNDX, NFILE, fnm), 1, &isize, &index, &grpname);

    /* initialize arrays */
    snew(d, isize);
    snew(dtot, isize);
    snew(dtot2, isize);
    if (bNMR)
    {
        snew(dtot1_3, isize);
        snew(dtot1_6, isize);
    }
    snew(mean, isize);
    snew(rms, isize);
    snew(rmsc, isize);
    snew(d_r, isize);
    snew(resnr, isize);
    for (i = 0; (i < isize); i++)
    {
        snew(d[i], isize);
        snew(dtot[i], isize);
        snew(dtot2[i], isize);
        if (bNMR)
        {
            snew(dtot1_3[i], isize);
            snew(dtot1_6[i], isize);
        }
        snew(mean[i], isize);
        snew(rms[i], isize);
        snew(rmsc[i], isize);
        snew(d_r[i], isize);
        resnr[i] = i+1;
    }

    /*set box type*/
    calc_dist(isize, index, x, ePBC, box, d_r);
    sfree(x);

    /*open output files*/
    fp = xvgropen(ftp2fn(efXVG, NFILE, fnm), "RMS Deviation", "Time (ps)", "RMSD (nm)",
                  oenv);
    if (output_env_get_print_xvgr_codes(oenv))
    {
        fprintf(fp, "@ subtitle \"of distances between %s atoms\"\n", grpname);
    }

    /*do a first step*/
    natom  = read_first_x(oenv, &status, ftp2fn(efTRX, NFILE, fnm), &t, &x, box);
    teller = 0;

    do
    {
        calc_dist_tot(isize, index, x, ePBC, box, d, dtot, dtot2, bNMR, dtot1_3, dtot1_6);

        rmsnow = rms_diff(isize, d, d_r);
        fprintf(fp, "%g  %g\n", t, rmsnow);
        teller++;
    }
    while (read_next_x(oenv, status, &t, x, box));
    fprintf(stderr, "\n");

    xvgrclose(fp);

    close_trj(status);

    calc_rms(isize, teller, dtot, dtot2, rms, &rmsmax, rmsc, &rmscmax, mean, &meanmax);
    fprintf(stderr, "rmsmax = %g, rmscmax = %g\n", rmsmax, rmscmax);

    if (bNMR)
    {
        calc_nmr(isize, teller, dtot1_3, dtot1_6, &max1_3, &max1_6);
    }

    if (scalemax > -1.0)
    {
        rmsmax  = scalemax;
        rmscmax = scalemax;
        meanmax = scalemax;
        max1_3  = scalemax;
        max1_6  = scalemax;
    }

    if (bNOE)
    {
        /* make list of noe atom groups */
        snew(noe_index, isize+1);
        snew(noe_gr, isize);
        gnr = analyze_noe_equivalent(opt2fn_null("-equiv", NFILE, fnm),
                                     atoms, isize, index, bSumH, noe_index, noe_gr);
        fprintf(stdout, "Found %d non-equivalent atom-groups in %d atoms\n",
                gnr, isize);
        /* make half matrix of of noe-group distances from atom distances */
        snew(noe, gnr);
        for (i = 0; i < gnr; i++)
        {
            snew(noe[i], gnr);
        }
        calc_noe(isize, noe_index, dtot1_3, dtot1_6, gnr, noe);
    }

    rlo.r = 1.0, rlo.g = 1.0, rlo.b = 1.0;
    rhi.r = 0.0, rhi.g = 0.0, rhi.b = 0.0;

    if (bRMS)
    {
        write_xpm(opt2FILE("-rms", NFILE, fnm, "w"), 0,
                  "RMS of distance", "RMS (nm)", "Atom Index", "Atom Index",
                  isize, isize, resnr, resnr, rms, 0.0, rmsmax, rlo, rhi, &nlevels);
    }

    if (bScale)
    {
        write_xpm(opt2FILE("-scl", NFILE, fnm, "w"), 0,
                  "Relative RMS", "RMS", "Atom Index", "Atom Index",
                  isize, isize, resnr, resnr, rmsc, 0.0, rmscmax, rlo, rhi, &nlevels);
    }

    if (bMean)
    {
        write_xpm(opt2FILE("-mean", NFILE, fnm, "w"), 0,
                  "Mean Distance", "Distance (nm)", "Atom Index", "Atom Index",
                  isize, isize, resnr, resnr, mean, 0.0, meanmax, rlo, rhi, &nlevels);
    }

    if (bNMR3)
    {
        write_xpm(opt2FILE("-nmr3", NFILE, fnm, "w"), 0, "1/r^3 averaged distances",
                  "Distance (nm)", "Atom Index", "Atom Index",
                  isize, isize, resnr, resnr, dtot1_3, 0.0, max1_3, rlo, rhi, &nlevels);
    }
    if (bNMR6)
    {
        write_xpm(opt2FILE("-nmr6", NFILE, fnm, "w"), 0, "1/r^6 averaged distances",
                  "Distance (nm)", "Atom Index", "Atom Index",
                  isize, isize, resnr, resnr, dtot1_6, 0.0, max1_6, rlo, rhi, &nlevels);
    }

    if (bNOE)
    {
        write_noe(opt2FILE("-noe", NFILE, fnm, "w"), gnr, noe, noe_gr, scalemax);
    }

    do_view(oenv, ftp2fn(efXVG, NFILE, fnm), NULL);

    return 0;
}
