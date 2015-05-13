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
#include <stdio.h>
#include <string.h>

#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/gmxana/gstat.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/legacyheaders/txtdump.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/listed-forces/bonded.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

void print_one(const output_env_t oenv, const char *base, const char *name,
               const char *title, const char *ylabel, int nf, real time[],
               real data[])
{
    FILE *fp;
    char  buf[256], t2[256];
    int   k;

    sprintf(buf, "%s%s.xvg", base, name);
    fprintf(stderr, "\rPrinting %s  ", buf);
    sprintf(t2, "%s %s", title, name);
    fp = xvgropen(buf, t2, "Time (ps)", ylabel, oenv);
    for (k = 0; (k < nf); k++)
    {
        fprintf(fp, "%10g  %10g\n", time[k], data[k]);
    }
    xvgrclose(fp);
}

static int calc_RBbin(real phi, int gmx_unused multiplicity, real gmx_unused core_frac)
{
    /* multiplicity and core_frac NOT used,
     * just given to enable use of pt-to-fn in caller low_ana_dih_trans*/
    static const real r30  = M_PI/6.0;
    static const real r90  = M_PI/2.0;
    static const real r150 = M_PI*5.0/6.0;

    if ((phi < r30) && (phi > -r30))
    {
        return 1;
    }
    else if ((phi > -r150) && (phi < -r90))
    {
        return 2;
    }
    else if ((phi < r150) && (phi > r90))
    {
        return 3;
    }
    return 0;
}

static int calc_Nbin(real phi, int multiplicity, real core_frac)
{
    static const real r360 = 360*DEG2RAD;
    real              rot_width, core_width, core_offset, low, hi;
    int               bin;
    /* with multiplicity 3 and core_frac 0.5
     * 0<g(-)<120, 120<t<240, 240<g(+)<360
     * 0< bin0 < 30, 30<bin1<90, 90<bin0<150, 150<bin2<210, 210<bin0<270, 270<bin3<330, 330<bin0<360
     * so with multiplicity 3, bin1 is core g(-), bin2 is core t, bin3 is
       core g(+), bin0 is between rotamers */
    if (phi < 0)
    {
        phi += r360;
    }

    rot_width   = 360/multiplicity;
    core_width  = core_frac * rot_width;
    core_offset = (rot_width - core_width)/2.0;
    for (bin = 1; bin <= multiplicity; bin++)
    {
        low  = ((bin - 1) * rot_width ) + core_offset;
        hi   = ((bin - 1) * rot_width ) + core_offset + core_width;
        low *= DEG2RAD;
        hi  *= DEG2RAD;
        if ((phi > low) && (phi < hi))
        {
            return bin;
        }
    }
    return 0;
}

void ana_dih_trans(const char *fn_trans, const char *fn_histo,
                   real **dih, int nframes, int nangles,
                   const char *grpname, real *time, gmx_bool bRb,
                   const output_env_t oenv)
{
    /* just a wrapper; declare extra args, then chuck away at end. */
    int      maxchi = 0;
    t_dlist *dlist;
    int     *multiplicity;
    int      nlist = nangles;
    int      k;

    snew(dlist, nlist);
    snew(multiplicity, nangles);
    for (k = 0; (k < nangles); k++)
    {
        multiplicity[k] = 3;
    }

    low_ana_dih_trans(TRUE, fn_trans, TRUE, fn_histo, maxchi,
                      dih, nlist, dlist, nframes,
                      nangles, grpname, multiplicity, time, bRb, 0.5, oenv);
    sfree(dlist);
    sfree(multiplicity);

}

void low_ana_dih_trans(gmx_bool bTrans, const char *fn_trans,
                       gmx_bool bHisto, const char *fn_histo, int maxchi,
                       real **dih, int nlist, t_dlist dlist[], int nframes,
                       int nangles, const char *grpname, int multiplicity[],
                       real *time, gmx_bool bRb, real core_frac,
                       const output_env_t oenv)
{
    FILE *fp;
    int  *tr_f, *tr_h;
    char  title[256];
    int   i, j, k, Dih, ntrans;
    int   cur_bin, new_bin;
    real  ttime, tt;
    real *rot_occ[NROT];
    int   (*calc_bin)(real, int, real);
    real  dt;

    if (1 <= nframes)
    {
        return;
    }
    /* Assumes the frames are equally spaced in time */
    dt = (time[nframes-1]-time[0])/(nframes-1);

    /* Analysis of dihedral transitions */
    fprintf(stderr, "Now calculating transitions...\n");

    if (bRb)
    {
        calc_bin = calc_RBbin;
    }
    else
    {
        calc_bin = calc_Nbin;
    }

    for (k = 0; k < NROT; k++)
    {
        snew(rot_occ[k], nangles);
        for (i = 0; (i < nangles); i++)
        {
            rot_occ[k][i] = 0;
        }
    }
    snew(tr_h, nangles);
    snew(tr_f, nframes);

    /* dih[i][j] is the dihedral angle i in frame j  */
    ntrans = 0;
    for (i = 0; (i < nangles); i++)
    {

        /*#define OLDIE*/
#ifdef OLDIE
        mind = maxd = prev = dih[i][0];
#else
        cur_bin = calc_bin(dih[i][0], multiplicity[i], core_frac);
        rot_occ[cur_bin][i]++;
#endif
        for (j = 1; (j < nframes); j++)
        {
            new_bin = calc_bin(dih[i][j], multiplicity[i], core_frac);
            rot_occ[new_bin][i]++;
#ifndef OLDIE
            if (cur_bin == 0)
            {
                cur_bin = new_bin;
            }
            else if ((new_bin != 0) && (cur_bin != new_bin))
            {
                cur_bin = new_bin;
                tr_f[j]++;
                tr_h[i]++;
                ntrans++;
            }
#else
            /* why is all this md rubbish periodic? Remove 360 degree periodicity */
            if ( (dih[i][j] - prev) > M_PI)
            {
                dih[i][j] -= 2*M_PI;
            }
            else if ( (dih[i][j] - prev) < -M_PI)
            {
                dih[i][j] += 2*M_PI;
            }

            prev = dih[i][j];

            mind = min(mind, dih[i][j]);
            maxd = max(maxd, dih[i][j]);
            if ( (maxd - mind) > 2*M_PI/3) /* or 120 degrees, assuming       */
            {                              /* multiplicity 3. Not so general.*/
                tr_f[j]++;
                tr_h[i]++;
                maxd = mind = dih[i][j]; /* get ready for next transition  */
                ntrans++;
            }
#endif
        } /* end j */
        for (k = 0; k < NROT; k++)
        {
            rot_occ[k][i] /= nframes;
        }
    } /* end i */
    fprintf(stderr, "Total number of transitions: %10d\n", ntrans);
    if (ntrans > 0)
    {
        ttime = (dt*nframes*nangles)/ntrans;
        fprintf(stderr, "Time between transitions:    %10.3f ps\n", ttime);
    }

    /* new by grs - copy transitions from tr_h[] to dlist->ntr[]
     * and rotamer populations from rot_occ to dlist->rot_occ[]
     * based on fn histogramming in g_chi. diff roles for i and j here */

    j = 0;
    for (Dih = 0; (Dih < NONCHI+maxchi); Dih++)
    {
        for (i = 0; (i < nlist); i++)
        {
            if (((Dih  < edOmega) ) ||
                ((Dih == edOmega) && (has_dihedral(edOmega, &(dlist[i])))) ||
                ((Dih  > edOmega) && (dlist[i].atm.Cn[Dih-NONCHI+3] != -1)))
            {
                /* grs debug  printf("Not OK? i %d j %d Dih %d \n", i, j, Dih) ; */
                dlist[i].ntr[Dih] = tr_h[j];
                for (k = 0; k < NROT; k++)
                {
                    dlist[i].rot_occ[Dih][k] = rot_occ[k][j];
                }
                j++;
            }
        }
    }

    /* end addition by grs */

    if (bTrans)
    {
        sprintf(title, "Number of transitions: %s", grpname);
        fp = xvgropen(fn_trans, title, "Time (ps)", "# transitions/timeframe", oenv);
        for (j = 0; (j < nframes); j++)
        {
            fprintf(fp, "%10.3f  %10d\n", time[j], tr_f[j]);
        }
        xvgrclose(fp);
    }

    /* Compute histogram from # transitions per dihedral */
    /* Use old array */
    for (j = 0; (j < nframes); j++)
    {
        tr_f[j] = 0;
    }
    for (i = 0; (i < nangles); i++)
    {
        tr_f[tr_h[i]]++;
    }
    for (j = nframes; ((tr_f[j-1] == 0) && (j > 0)); j--)
    {
        ;
    }

    ttime = dt*nframes;
    if (bHisto)
    {
        sprintf(title, "Transition time: %s", grpname);
        fp = xvgropen(fn_histo, title, "Time (ps)", "#", oenv);
        for (i = j-1; (i > 0); i--)
        {
            if (tr_f[i] != 0)
            {
                fprintf(fp, "%10.3f  %10d\n", ttime/i, tr_f[i]);
            }
        }
        xvgrclose(fp);
    }

    sfree(tr_f);
    sfree(tr_h);
    for (k = 0; k < NROT; k++)
    {
        sfree(rot_occ[k]);
    }

}

void mk_multiplicity_lookup (int *multiplicity, int maxchi,
                             int nlist, t_dlist dlist[], int nangles)
{
    /* new by grs - for dihedral j (as in dih[j]) get multiplicity from dlist
     * and store in multiplicity[j]
     */

    int  j, Dih, i;
    char name[4];

    j = 0;
    for (Dih = 0; (Dih < NONCHI+maxchi); Dih++)
    {
        for (i = 0; (i < nlist); i++)
        {
            strncpy(name, dlist[i].name, 3);
            name[3] = '\0';
            if (((Dih  < edOmega) ) ||
                ((Dih == edOmega) && (has_dihedral(edOmega, &(dlist[i])))) ||
                ((Dih  > edOmega) && (dlist[i].atm.Cn[Dih-NONCHI+3] != -1)))
            {
                /* default - we will correct the rest below */
                multiplicity[j] = 3;

                /* make omegas 2fold, though doesn't make much more sense than 3 */
                if (Dih == edOmega && (has_dihedral(edOmega, &(dlist[i]))))
                {
                    multiplicity[j] = 2;
                }

                /* dihedrals to aromatic rings, COO, CONH2 or guanidinium are 2fold*/
                if (Dih > edOmega && (dlist[i].atm.Cn[Dih-NONCHI+3] != -1))
                {
                    if ( ((strstr(name, "PHE") != NULL) && (Dih == edChi2))  ||
                         ((strstr(name, "TYR") != NULL) && (Dih == edChi2))  ||
                         ((strstr(name, "PTR") != NULL) && (Dih == edChi2))  ||
                         ((strstr(name, "TRP") != NULL) && (Dih == edChi2))  ||
                         ((strstr(name, "HIS") != NULL) && (Dih == edChi2))  ||
                         ((strstr(name, "GLU") != NULL) && (Dih == edChi3))  ||
                         ((strstr(name, "ASP") != NULL) && (Dih == edChi2))  ||
                         ((strstr(name, "GLN") != NULL) && (Dih == edChi3))  ||
                         ((strstr(name, "ASN") != NULL) && (Dih == edChi2))  ||
                         ((strstr(name, "ARG") != NULL) && (Dih == edChi4))  )
                    {
                        multiplicity[j] = 2;
                    }
                }
                j++;
            }
        }
    }
    if (j < nangles)
    {
        fprintf(stderr, "WARNING: not all dihedrals found in topology (only %d out of %d)!\n",
                j, nangles);
    }
    /* Check for remaining dihedrals */
    for (; (j < nangles); j++)
    {
        multiplicity[j] = 3;
    }

}

void mk_chi_lookup (int **lookup, int maxchi,
                    int nlist, t_dlist dlist[])
{

    /* by grs. should rewrite everything to use this. (but haven't,
     * and at mmt only used in get_chi_product_traj
     * returns the dihed number given the residue number (from-0)
     * and chi (from-0) nr. -1 for chi undefined for that res (eg gly, ala..)*/

    int i, j, Dih, Chi;

    j = 0;
    /* NONCHI points to chi1, therefore we have to start counting there. */
    for (Dih = NONCHI; (Dih < NONCHI+maxchi); Dih++)
    {
        for (i = 0; (i < nlist); i++)
        {
            Chi = Dih - NONCHI;
            if (((Dih  < edOmega) ) ||
                ((Dih == edOmega) && (has_dihedral(edOmega, &(dlist[i])))) ||
                ((Dih  > edOmega) && (dlist[i].atm.Cn[Dih-NONCHI+3] != -1)))
            {
                /* grs debug  printf("Not OK? i %d j %d Dih %d \n", i, j, Dih) ; */
                if (Dih > edOmega)
                {
                    lookup[i][Chi] = j;
                }
                j++;
            }
            else
            {
                lookup[i][Chi] = -1;
            }
        }
    }

}


void get_chi_product_traj (real **dih, int nframes, int nlist,
                           int maxchi, t_dlist dlist[], real time[],
                           int **lookup, int *multiplicity, gmx_bool bRb, gmx_bool bNormalize,
                           real core_frac, gmx_bool bAll, const char *fnall,
                           const output_env_t oenv)
{

    gmx_bool bRotZero, bHaveChi = FALSE;
    int      accum = 0, index, i, j, k, Xi, n, b;
    real    *chi_prtrj;
    int     *chi_prhist;
    int      nbin;
    FILE    *fp, *fpall;
    char     hisfile[256], histitle[256], *namept;

    int      (*calc_bin)(real, int, real);

    /* Analysis of dihedral transitions */
    fprintf(stderr, "Now calculating Chi product trajectories...\n");

    if (bRb)
    {
        calc_bin = calc_RBbin;
    }
    else
    {
        calc_bin = calc_Nbin;
    }

    snew(chi_prtrj, nframes);

    /* file for info on all residues */
    if (bNormalize)
    {
        fpall = xvgropen(fnall, "Cumulative Rotamers", "Residue", "Probability", oenv);
    }
    else
    {
        fpall = xvgropen(fnall, "Cumulative Rotamers", "Residue", "# Counts", oenv);
    }

    for (i = 0; (i < nlist); i++)
    {

        /* get nbin, the nr. of cumulative rotamers that need to be considered */
        nbin = 1;
        for (Xi = 0; Xi < maxchi; Xi++)
        {
            index = lookup[i][Xi]; /* chi_(Xi+1) of res i (-1 if off end) */
            if (index >= 0)
            {
                n    = multiplicity[index];
                nbin = n*nbin;
            }
        }
        nbin += 1; /* for the "zero rotamer", outside the core region */

        for (j = 0; (j < nframes); j++)
        {

            bRotZero = FALSE;
            bHaveChi = TRUE;
            index    = lookup[i][0]; /* index into dih of chi1 of res i */
            if (index == -1)
            {
                b        = 0;
                bRotZero = TRUE;
                bHaveChi = FALSE;
            }
            else
            {
                b     = calc_bin(dih[index][j], multiplicity[index], core_frac);
                accum = b - 1;
                if (b == 0)
                {
                    bRotZero = TRUE;
                }
                for (Xi = 1; Xi < maxchi; Xi++)
                {
                    index = lookup[i][Xi]; /* chi_(Xi+1) of res i (-1 if off end) */
                    if (index >= 0)
                    {
                        n     = multiplicity[index];
                        b     = calc_bin(dih[index][j], n, core_frac);
                        accum = n * accum + b - 1;
                        if (b == 0)
                        {
                            bRotZero = TRUE;
                        }
                    }
                }
                accum++;
            }
            if (bRotZero)
            {
                chi_prtrj[j] = 0.0;
            }
            else
            {
                chi_prtrj[j] = accum;
                if (accum+1 > nbin)
                {
                    nbin = accum+1;
                }
            }
        }
        if (bHaveChi)
        {

            if (bAll)
            {
                /* print cuml rotamer vs time */
                print_one(oenv, "chiproduct", dlist[i].name, "chi product for",
                          "cumulative rotamer", nframes, time, chi_prtrj);
            }

            /* make a histogram pf culm. rotamer occupancy too */
            snew(chi_prhist, nbin);
            make_histo(NULL, nframes, chi_prtrj, nbin, chi_prhist, 0, nbin);
            if (bAll)
            {
                sprintf(hisfile, "histo-chiprod%s.xvg", dlist[i].name);
                sprintf(histitle, "cumulative rotamer distribution for %s", dlist[i].name);
                fprintf(stderr, "  and %s  ", hisfile);
                fp = xvgropen(hisfile, histitle, "number", "", oenv);
                if (output_env_get_print_xvgr_codes(oenv))
                {
                    fprintf(fp, "@ xaxis tick on\n");
                    fprintf(fp, "@ xaxis tick major 1\n");
                    fprintf(fp, "@ type xy\n");
                }
                for (k = 0; (k < nbin); k++)
                {
                    if (bNormalize)
                    {
                        fprintf(fp, "%5d  %10g\n", k, (1.0*chi_prhist[k])/nframes);
                    }
                    else
                    {
                        fprintf(fp, "%5d  %10d\n", k, chi_prhist[k]);
                    }
                }
                fprintf(fp, "%s\n", output_env_get_print_xvgr_codes(oenv) ? "&" : "");
                xvgrclose(fp);
            }

            /* and finally print out occupancies to a single file */
            /* get the gmx from-1 res nr by setting a ptr to the number part
             * of dlist[i].name - potential bug for 4-letter res names... */
            namept = dlist[i].name + 3;
            fprintf(fpall, "%5s ", namept);
            for (k = 0; (k < nbin); k++)
            {
                if (bNormalize)
                {
                    fprintf(fpall, "  %10g", (1.0*chi_prhist[k])/nframes);
                }
                else
                {
                    fprintf(fpall, "  %10d", chi_prhist[k]);
                }
            }
            fprintf(fpall, "\n");

            sfree(chi_prhist);
            /* histogram done */
        }
    }

    sfree(chi_prtrj);
    xvgrclose(fpall);
    fprintf(stderr, "\n");

}

void calc_distribution_props(int nh, int histo[], real start,
                             int nkkk, t_karplus kkk[],
                             real *S2)
{
    real d, dc, ds, c1, c2, tdc, tds;
    real fac, ang, invth, Jc;
    int  i, j, th;

    if (nh == 0)
    {
        gmx_fatal(FARGS, "No points in histogram (%s, %d)", __FILE__, __LINE__);
    }
    fac = 2*M_PI/nh;

    /* Compute normalisation factor */
    th = 0;
    for (j = 0; (j < nh); j++)
    {
        th += histo[j];
    }
    invth = 1.0/th;

    for (i = 0; (i < nkkk); i++)
    {
        kkk[i].Jc    = 0;
        kkk[i].Jcsig = 0;
    }
    tdc = 0, tds = 0;
    for (j = 0; (j < nh); j++)
    {
        d    = invth*histo[j];
        ang  = j*fac-start;
        c1   = cos(ang);
        c2   = c1*c1;
        dc   = d*c1;
        ds   = d*sin(ang);
        tdc += dc;
        tds += ds;
        for (i = 0; (i < nkkk); i++)
        {
            c1            = cos(ang+kkk[i].offset);
            c2            = c1*c1;
            Jc            = (kkk[i].A*c2 + kkk[i].B*c1 + kkk[i].C);
            kkk[i].Jc    += histo[j]*Jc;
            kkk[i].Jcsig += histo[j]*sqr(Jc);
        }
    }
    for (i = 0; (i < nkkk); i++)
    {
        kkk[i].Jc    /= th;
        kkk[i].Jcsig  = sqrt(kkk[i].Jcsig/th-sqr(kkk[i].Jc));
    }
    *S2 = tdc*tdc+tds*tds;
}

static void calc_angles(struct t_pbc *pbc,
                        int n3, atom_id index[], real ang[], rvec x_s[])
{
    int  i, ix, t1, t2;
    rvec r_ij, r_kj;
    real costh;

    for (i = ix = 0; (ix < n3); i++, ix += 3)
    {
        ang[i] = bond_angle(x_s[index[ix]], x_s[index[ix+1]], x_s[index[ix+2]],
                            pbc, r_ij, r_kj, &costh, &t1, &t2);
    }
    if (debug)
    {
        fprintf(debug, "Angle[0]=%g, costh=%g, index0 = %d, %d, %d\n",
                ang[0], costh, index[0], index[1], index[2]);
        pr_rvec(debug, 0, "rij", r_ij, DIM, TRUE);
        pr_rvec(debug, 0, "rkj", r_kj, DIM, TRUE);
    }
}

static real calc_fraction(real angles[], int nangles)
{
    int  i;
    real trans = 0, gauche = 0;
    real angle;

    for (i = 0; i < nangles; i++)
    {
        angle = angles[i] * RAD2DEG;

        if (angle > 135 && angle < 225)
        {
            trans += 1.0;
        }
        else if (angle > 270 && angle < 330)
        {
            gauche += 1.0;
        }
        else if (angle < 90 && angle > 30)
        {
            gauche += 1.0;
        }
    }
    if (trans+gauche > 0)
    {
        return trans/(trans+gauche);
    }
    else
    {
        return 0;
    }
}

static void calc_dihs(struct t_pbc *pbc,
                      int n4, atom_id index[], real ang[], rvec x_s[])
{
    int  i, ix, t1, t2, t3;
    rvec r_ij, r_kj, r_kl, m, n;
    real sign, aaa;

    for (i = ix = 0; (ix < n4); i++, ix += 4)
    {
        aaa = dih_angle(x_s[index[ix]], x_s[index[ix+1]], x_s[index[ix+2]],
                        x_s[index[ix+3]], pbc,
                        r_ij, r_kj, r_kl, m, n,
                        &sign, &t1, &t2, &t3);

        ang[i] = aaa; /* not taking into account ryckaert bellemans yet */
    }
}

void make_histo(FILE *log,
                int ndata, real data[], int npoints, int histo[],
                real minx, real maxx)
{
    double dx;
    int    i, ind;

    if (minx == maxx)
    {
        minx = maxx = data[0];
        for (i = 1; (i < ndata); i++)
        {
            minx = min(minx, data[i]);
            maxx = max(maxx, data[i]);
        }
        fprintf(log, "Min data: %10g  Max data: %10g\n", minx, maxx);
    }
    dx = (double)npoints/(maxx-minx);
    if (debug)
    {
        fprintf(debug,
                "Histogramming: ndata=%d, nhisto=%d, minx=%g,maxx=%g,dx=%g\n",
                ndata, npoints, minx, maxx, dx);
    }
    for (i = 0; (i < ndata); i++)
    {
        ind = (data[i]-minx)*dx;
        if ((ind >= 0) && (ind < npoints))
        {
            histo[ind]++;
        }
        else
        {
            fprintf(log, "index = %d, data[%d] = %g\n", ind, i, data[i]);
        }
    }
}

void normalize_histo(int npoints, int histo[], real dx, real normhisto[])
{
    int    i;
    double d, fac;

    d = 0;
    for (i = 0; (i < npoints); i++)
    {
        d += dx*histo[i];
    }
    if (d == 0)
    {
        fprintf(stderr, "Empty histogram!\n");
        return;
    }
    fac = 1.0/d;
    for (i = 0; (i < npoints); i++)
    {
        normhisto[i] = fac*histo[i];
    }
}

void read_ang_dih(const char *trj_fn,
                  gmx_bool bAngles, gmx_bool bSaveAll, gmx_bool bRb, gmx_bool bPBC,
                  int maxangstat, int angstat[],
                  int *nframes, real **time,
                  int isize, atom_id index[],
                  real **trans_frac,
                  real **aver_angle,
                  real *dih[],
                  const output_env_t oenv)
{
    struct t_pbc *pbc;
    t_trxstatus  *status;
    int           i, angind, natoms, total, teller;
    int           nangles, n_alloc;
    real          t, fraction, pifac, aa, angle;
    real         *angles[2];
    matrix        box;
    rvec         *x;
    int           cur = 0;
#define prev (1-cur)

    snew(pbc, 1);
    natoms = read_first_x(oenv, &status, trj_fn, &t, &x, box);

    if (bAngles)
    {
        nangles = isize/3;
        pifac   = M_PI;
    }
    else
    {
        nangles = isize/4;
        pifac   = 2.0*M_PI;
    }
    snew(angles[cur], nangles);
    snew(angles[prev], nangles);

    /* Start the loop over frames */
    total       = 0;
    teller      = 0;
    n_alloc     = 0;
    *time       = NULL;
    *trans_frac = NULL;
    *aver_angle = NULL;

    do
    {
        if (teller >= n_alloc)
        {
            n_alloc += 100;
            if (bSaveAll)
            {
                for (i = 0; (i < nangles); i++)
                {
                    srenew(dih[i], n_alloc);
                }
            }
            srenew(*time, n_alloc);
            srenew(*trans_frac, n_alloc);
            srenew(*aver_angle, n_alloc);
        }

        (*time)[teller] = t;

        if (pbc)
        {
            set_pbc(pbc, -1, box);
        }

        if (bAngles)
        {
            calc_angles(pbc, isize, index, angles[cur], x);
        }
        else
        {
            calc_dihs(pbc, isize, index, angles[cur], x);

            /* Trans fraction */
            fraction              = calc_fraction(angles[cur], nangles);
            (*trans_frac)[teller] = fraction;

            /* Change Ryckaert-Bellemans dihedrals to polymer convention
             * Modified 990913 by Erik:
             * We actually shouldn't change the convention, since it's
             * calculated from polymer above, but we change the intervall
             * from [-180,180] to [0,360].
             */
            if (bRb)
            {
                for (i = 0; (i < nangles); i++)
                {
                    if (angles[cur][i] <= 0.0)
                    {
                        angles[cur][i] += 2*M_PI;
                    }
                }
            }

            /* Periodicity in dihedral space... */
            if (bPBC)
            {
                for (i = 0; (i < nangles); i++)
                {
                    real dd = angles[cur][i];
                    angles[cur][i] = atan2(sin(dd), cos(dd));
                }
            }
            else
            {
                if (teller > 1)
                {
                    for (i = 0; (i < nangles); i++)
                    {
                        while (angles[cur][i] <= angles[prev][i] - M_PI)
                        {
                            angles[cur][i] += 2*M_PI;
                        }
                        while (angles[cur][i] > angles[prev][i] + M_PI)
                        {
                            angles[cur][i] -= 2*M_PI;
                        }
                    }
                }
            }
        }

        /* Average angles */
        aa = 0;
        for (i = 0; (i < nangles); i++)
        {
            aa = aa+angles[cur][i];

            /* angle in rad / 2Pi * max determines bin. bins go from 0 to maxangstat,
               even though scale goes from -pi to pi (dihedral) or -pi/2 to pi/2
               (angle) Basically: translate the x-axis by Pi. Translate it back by
               -Pi when plotting.
             */

            angle = angles[cur][i];
            if (!bAngles)
            {
                while (angle < -M_PI)
                {
                    angle += 2*M_PI;
                }
                while (angle >= M_PI)
                {
                    angle -= 2*M_PI;
                }

                angle += M_PI;
            }

            /* Update the distribution histogram */
            angind = (int) ((angle*maxangstat)/pifac + 0.5);
            if (angind == maxangstat)
            {
                angind = 0;
            }
            if ( (angind < 0) || (angind >= maxangstat) )
            {
                /* this will never happen */
                gmx_fatal(FARGS, "angle (%f) index out of range (0..%d) : %d\n",
                          angle, maxangstat, angind);
            }

            angstat[angind]++;
            if (angind == maxangstat)
            {
                fprintf(stderr, "angle %d fr %d = %g\n", i, cur, angle);
            }

            total++;
        }

        /* average over all angles */
        (*aver_angle)[teller] = (aa/nangles);

        /* this copies all current dih. angles to dih[i], teller is frame */
        if (bSaveAll)
        {
            for (i = 0; i < nangles; i++)
            {
                dih[i][teller] = angles[cur][i];
            }
        }

        /* Swap buffers */
        cur = prev;

        /* Increment loop counter */
        teller++;
    }
    while (read_next_x(oenv, status, &t, x, box));
    close_trj(status);

    sfree(x);
    sfree(angles[cur]);
    sfree(angles[prev]);

    *nframes = teller;
}
