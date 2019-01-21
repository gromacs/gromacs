/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017,2018, by the GROMACS development team, led by
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
#include <cstring>

#include <algorithm>

#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/commandline/viewit.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/groio.h"
#include "gromacs/fileio/pdbio.h"
#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/math/do_fit.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/rmpbc.h"
#include "gromacs/topology/index.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"

static const int NOTSET = -9368163;

static void calc_rm_cm(gmx::ArrayRef<const int> index,
                       gmx::ArrayRef<const AtomInfo> atoms, rvec x[], rvec xcm)
{
    real tm, m;

    /* calculate and remove center of mass of reference structure */
    tm = 0;
    clear_rvec(xcm);
    for (int i = 0; i < index.size(); i++)
    {
        m = atoms[index[i]].m_;
        for (int d = 0; d < DIM; d++)
        {
            xcm[d] += m*x[index[i]][d];
        }
        tm += m;
    }
    svmul(1/tm, xcm, xcm);
    for (int i = 0; i < atoms.size(); i++)
    {
        rvec_dec(x[i], xcm);
    }
}

static int build_res_index(gmx::ArrayRef<const int> index,
                           gmx::ArrayRef<const AtomInfo> atom, gmx::ArrayRef<int> rindex)
{
    int r         = 0;
    rindex[r] = atom[index[0]].resind_;
    r++;
    for (int i = 1; i < index.size(); i++)
    {
        if (atom[index[i]].resind_ != rindex[r-1])
        {
            rindex[r] = atom[index[i]].resind_;
            r++;
        }
    }

    return r;
}

static int find_res_end(int i, gmx::ArrayRef<const int> index, gmx::ArrayRef<const AtomInfo> atoms)
{
    int rnr;

    rnr = atoms[index[i]].resind_;
    while (i < index.size() && atoms[index[i]].resind_ == rnr)
    {
        i++;
    }
    return i;
}

static int debug_strcmp(char s1[], char s2[])
{
    if (debug)
    {
        fprintf(debug, " %s-%s", s1, s2);
    }
    return std::strcmp(s1, s2);
}

static int find_next_match_atoms_in_res(int *i1, gmx::ArrayRef<const int> index1,
                                        int m1, gmx::ArrayRef<const AtomInfo> atoms1,
                                        int *i2, gmx::ArrayRef<const int> index2,
                                        int m2, gmx::ArrayRef<const AtomInfo> atoms2)
{
    int      dx, dy, dmax, cmp;
    gmx_bool bFW = FALSE;

    cmp  = NOTSET;
    dmax = std::max(m1-*i1, m2-*i2);
    for (dx = 0, dy = 0; dx < dmax && cmp != 0; dx++)
    {
        for (dy = dx; dy < dmax && cmp != 0; dy++)
        {
            if (dx || dy)
            {
                if (debug)
                {
                    fprintf(debug, ".");
                }
                cmp = NOTSET;
                if (*i1+dx < m1 && *i2+dy < m2)
                {
                    bFW = TRUE;
                    cmp = debug_strcmp(*atoms1[index1[*i1+dx]].atomname, *atoms2[index2[*i2+dy]].atomname);
                    if (debug)
                    {
                        fprintf(debug, "(%d %d)", *i1+dx, *i2+dy);
                    }
                }
                if (cmp != 0 && *i1+dy < m1 && *i2+dx < m2)
                {
                    bFW = FALSE;
                    cmp = debug_strcmp(*atoms1[index1[*i1+dy]].atomname, *atoms2[index2[*i2+dx]].atomname);
                    if (debug)
                    {
                        fprintf(debug, "(%d %d)", *i1+dy, *i2+dx);
                    }
                }
            }
        }
    }
    /* apparently, dx and dy are incremented one more time
       as the loops terminate: we correct this here */
    dx--;
    dy--;
    if (cmp == 0)
    {
        if (debug)
        {
            fprintf(debug, "{%d %d}", *i1 + (bFW ? dx : dy), *i2 + (bFW ? dy : dx) );
        }
        if (bFW)
        {
            *i1 += dx;
            *i2 += dy;
        }
        else
        {
            *i1 += dy;
            *i2 += dx;
        }
    }

    return cmp;
}

static int find_next_match_res(int *rnr1,
                               gmx::ArrayRef<const int> index1,
                               gmx::ArrayRef<const Residue> resinfo1,
                               int *rnr2,
                               gmx::ArrayRef<const int> index2,
                               gmx::ArrayRef<const Residue> resinfo2)
{
    int      dmax, cmp, rr1, rr2;
    gmx_bool bFW = FALSE, bFF = FALSE;

    rr1 = 0;
    while (rr1 < index1.size() && *rnr1 != index1[rr1])
    {
        rr1++;
    }
    rr2 = 0;
    while (rr2 < index2.size() && *rnr2 != index2[rr2])
    {
        rr2++;
    }

    cmp  = NOTSET;
    dmax = std::max(index1.size()-rr1, index2.size()-rr2);
    if (debug)
    {
        fprintf(debug, " R:%d-%zu:%d-%zu:%d ",
                rr1, index1.size(), rr2, index2.size(), dmax);
    }
    int dx = 0, dy = 0;
    for (; dx < dmax && cmp != 0; dx++)
    {
        for (dy = 0; dy <= dx && cmp != 0; dy++)
        {
            if (dx != dy)
            {
                cmp = NOTSET;
                if (rr1+dx < index1.size() && rr2+dy < index2.size())
                {
                    bFW = TRUE;
                    cmp = debug_strcmp(*resinfo1[index1[rr1+dx]].name_,
                                       *resinfo2[index2[rr2+dy]].name_);
                    if (debug)
                    {
                        fprintf(debug, "(%d %d)", rr1+dx, rr2+dy);
                    }
                }
                if (cmp != 0 && rr1+dy < index1.size() && rr2+dx < index2.size())
                {
                    bFW = FALSE;
                    cmp = debug_strcmp(*resinfo1[index1[rr1+dy]].name_,
                                       *resinfo2[index2[rr2+dx]].name_);
                    if (debug)
                    {
                        fprintf(debug, "(%d %d)", rr1+dy, rr2+dx);
                    }
                }
                if (dx != 0 && cmp != 0 && rr1+dx < index1.size() && rr2+dx < index2.size())
                {
                    bFF = TRUE;
                    cmp = debug_strcmp(*resinfo1[index1[rr1+dx]].name_,
                                       *resinfo2[index2[rr2+dx]].name_);
                    if (debug)
                    {
                        fprintf(debug, "(%d %d)", rr1+dx, rr2+dx);
                    }
                }
                else
                {
                    bFF = FALSE;
                }
            }
        }
        /* apparently, dx and dy are incremented one more time
           as the loops terminate: we correct this here */
        dy--;
    }
    dx--;
    /* if we skipped equal on both sides, only skip one residue
       to allow for single mutations of residues... */
    if (bFF)
    {
        if (debug)
        {
            fprintf(debug, "%d.%d.%dX%sX%s", dx, rr1, rr2,
                    *resinfo1[index1[rr1+1]].name_,
                    *resinfo2[index2[rr2+1]].name_);
        }
        dx = 1;
    }
    if (cmp == 0)
    {
        if (debug)
        {
            fprintf(debug, "!");
        }
        if (bFF)
        {
            rr1 += dx;
            rr2 += dx;
        }
        else if (bFW)
        {
            rr1 += dx;
            rr2 += dy;
        }
        else
        {
            rr1 += dy;
            rr2 += dx;
        }
        *rnr1 = index1[rr1];
        *rnr2 = index2[rr2];
    }

    return cmp;
}

static int find_first_atom_in_res(int rnr,
                                  gmx::ArrayRef<const int> index,
                                  gmx::ArrayRef<const AtomInfo> atom)
{
    int i = 0;
    while (i < index.size() && atom[index[i]].resind_ != rnr)
    {
        i++;
    }

    if (atom[index[i]].resind_ == rnr)
    {
        return i;
    }
    else
    {
        return NOTSET;
    }
}

static void find_matching_names(gmx::ArrayRef<int> index1,
                                gmx::ArrayRef<const AtomInfo> atoms1,
                                gmx::ArrayRef<const Residue> resinfo1,
                                gmx::ArrayRef<int> index2,
                                gmx::ArrayRef<const AtomInfo> atoms2,
                                gmx::ArrayRef<const Residue> resinfo2)
{
    int        i1, i2, ii1, ii2, m1, m2;
    int        atcmp, rescmp;
    int        rnr1, rnr2, prnr1, prnr2;

    /* set some handy shortcuts */

    /* build indexes of selected residues */
    std::vector<int> rindex1(resinfo1.size());
    build_res_index(index1, atoms1, rindex1);
    std::vector<int> rindex2(resinfo2.size());
    build_res_index(index2, atoms2, rindex2);

    i1    = i2 = 0;
    ii1   = ii2 = 0;
    atcmp = rescmp = 0;
    prnr1 = prnr2 = NOTSET;
    if (debug)
    {
        fprintf(debug, "Find matching names: %zu, %zu\n", index1.size(), index2.size());
    }
    while (atcmp == 0 && i1 < index1.size() && i2 < index2.size())
    {
        /* prologue */
        rnr1 = atoms1[index1[i1]].resind_;
        rnr2 = atoms2[index2[i2]].resind_;
        if (rnr1 != prnr1 || rnr2 != prnr2)
        {
            if (debug)
            {
                fprintf(debug, "R: %s%d %s%d\n",
                        *resinfo1[rnr1].name_, rnr1, *resinfo2[rnr2].name_, rnr2);
            }
            rescmp = std::strcmp(*resinfo1[rnr1].name_, *resinfo2[rnr2].name_);
        }
        if (debug)
        {
            fprintf(debug, "comparing %d %d", i1, i2);
        }
        atcmp = debug_strcmp(*atoms1[index1[i1]].atomname, *atoms2[index2[i2]].atomname);

        /* the works */
        if (atcmp != 0) /* no match -> find match within residues */
        {
            m1 = find_res_end(i1, index1, atoms1);
            m2 = find_res_end(i2, index2, atoms2);
            if (debug)
            {
                fprintf(debug, " [%d<%d %d<%d]", i1, m1, i2, m2);
            }
            atcmp = find_next_match_atoms_in_res(&i1, index1, m1, atoms1,
                                                 &i2, index2, m2, atoms2);
            if (debug)
            {
                fprintf(debug, " -> %d %d %s-%s", i1, i2,
                        *atoms1[index1[i1]].atomname, *atoms2[index2[i2]].atomname);
            }

        }
        if (atcmp != 0) /* still no match -> next residue(s) */
        {
            prnr1  = rnr1;
            prnr2  = rnr2;
            rescmp = find_next_match_res(&rnr1, rindex1, resinfo1,
                                         &rnr2, rindex2, resinfo2);
            if (rnr1 != prnr1)
            {
                i1 = find_first_atom_in_res(rnr1, index1, atoms1);
            }
            if (rnr2 != prnr2)
            {
                i2 = find_first_atom_in_res(rnr2, index2, atoms2);
            }
            if (debug)
            {
                fprintf(debug, " -> %s%d-%s%d %s%d-%s%d",
                        *resinfo1[rnr1].name_, rnr1,
                        *atoms1[index1[i1]].atomname, index1[i1],
                        *resinfo2[rnr2].name_, rnr2,
                        *atoms2[index2[i2]].atomname, index2[i2]);
            }
            m1 = find_res_end(i1, index1, atoms1);
            m2 = find_res_end(i2, index2, atoms2);
            if (debug)
            {
                fprintf(debug, " [%d<%d %d<%d]", i1, m1, i2, m2);
            }
            atcmp = find_next_match_atoms_in_res(&i1, index1, m1, atoms1,
                                                 &i2, index2, m2, atoms2);
            if (debug)
            {
                fprintf(debug, " -> %d %d %s-%s", i1, i2,
                        *atoms1[index1[i1]].atomname, *atoms2[index2[i2]].atomname);
            }
        }
        if (debug)
        {
            fprintf(debug, "(%d %d): %d %d\n", ii1, ii2, atcmp, rescmp);
        }
        if (atcmp == 0) /* if match -> store indices */
        {
            index1[ii1++] = index1[i1];
            index2[ii2++] = index2[i2];
        }
        i1++;
        i2++;

        /* epilogue */
        prnr1 = rnr1;
        prnr2 = rnr2;
    }

    if (ii1 != ii2)
    {
        gmx_fatal(FARGS, "DEATH HORROR: non-equal number of matching atoms!\n");
    }
    if (ii1 == i1 && ii2 == i2)
    {
        printf("All atoms in index groups 1 and 2 match\n");
    }
    else
    {
        if (i1 == i2 && ii1 == ii2)
        {
            printf("Both index groups modified from %d to %d atoms\n", i1, ii1);
        }
        else
        {
            if (ii1 != i1)
            {
                printf("Index group 1 modified from %d to %d atoms\n", i1, ii1);
            }
            if (ii2 != i2)
            {
                printf("Index group 2 modified from %d to %d atoms\n", i2, ii2);
            }
        }
    }
}
/* 1 */

int gmx_confrms(int argc, char *argv[])
{
    const char     *desc[] = {
        "[THISMODULE] computes the root mean square deviation (RMSD) of two",
        "structures after least-squares fitting the second structure on the first one.",
        "The two structures do NOT need to have the same number of atoms,",
        "only the two index groups used for the fit need to be identical.",
        "With [TT]-name[tt] only matching atom names from the selected groups",
        "will be used for the fit and RMSD calculation. This can be useful ",
        "when comparing mutants of a protein.",
        "[PAR]",
        "The superimposed structures are written to file. In a [REF].pdb[ref] file",
        "the two structures will be written as separate models",
        "(use [TT]rasmol -nmrpdb[tt]). Also in a [REF].pdb[ref] file, B-factors",
        "calculated from the atomic MSD values can be written with [TT]-bfac[tt].",
    };
    static gmx_bool bOne  = FALSE, bRmpbc = FALSE, bMW = TRUE, bName = FALSE,
                    bBfac = FALSE, bFit = TRUE, bLabel = FALSE;

    t_pargs  pa[] = {
        { "-one", FALSE, etBOOL, {&bOne},   "Only write the fitted structure to file" },
        { "-mw",  FALSE, etBOOL, {&bMW},    "Mass-weighted fitting and RMSD" },
        { "-pbc", FALSE, etBOOL, {&bRmpbc}, "Try to make molecules whole again" },
        { "-fit", FALSE, etBOOL, {&bFit},
          "Do least squares superposition of the target structure to the reference" },
        { "-name", FALSE, etBOOL, {&bName},
          "Only compare matching atom names" },
        { "-label", FALSE, etBOOL, {&bLabel},
          "Added chain labels A for first and B for second structure"},
        { "-bfac", FALSE, etBOOL, {&bBfac},
          "Output B-factors from atomic MSD values" }
    };
    t_filenm fnm[] = {
        { efTPS, "-f1",  "conf1.gro", ffREAD  },
        { efSTX, "-f2",  "conf2",     ffREAD  },
        { efSTO, "-o",   "fit.pdb",   ffWRITE },
        { efNDX, "-n1",  "fit1",      ffOPTRD },
        { efNDX, "-n2",  "fit2",      ffOPTRD },
        { efNDX, "-no",  "match",     ffOPTWR }
    };
#define NFILE asize(fnm)

    /* the two structure files */
    const char       *conf1file, *conf2file, *matchndxfile, *outfile;
    FILE             *fp;
    char             *name1, *name2;
    int               ePBC1, ePBC2;
    int               warn = 0;
    int               at;
    real             *w_rls, mass, totmass;
    rvec             *x1, *v1, *x2, *v2, *fit_x;
    matrix            box1, box2;

    gmx_output_env_t *oenv;

    /* counters */
    int     i, m;

    /* center of mass calculation */
    rvec    xcm1, xcm2;

    /* variables for fit */
    char    *groupnames1, *groupnames2;
    int      isize1, isize2;
    int     *index1, *index2;
    real     rms, msd, minmsd, maxmsd;
    real    *msds;


    if (!parse_common_args(&argc, argv, PCA_CAN_VIEW,
                           NFILE, fnm, asize(pa), pa, asize(desc), desc, 0, nullptr, &oenv))
    {
        return 0;
    }
    matchndxfile = opt2fn_null("-no", NFILE, fnm);
    conf1file    = ftp2fn(efTPS, NFILE, fnm);
    conf2file    = ftp2fn(efSTX, NFILE, fnm);

    /* reading reference structure from first structure file */
    fprintf(stderr, "\nReading first structure file\n");
    t_topology *top1 = new t_topology;
    read_tps_conf(conf1file, top1, &ePBC1, &x1, &v1, box1, TRUE);
    std::vector<AtomInfo> atoms1(top1->atoms.begin(), top1->atoms.end());
    std::vector<Residue> resinfo1(top1->resinfo.begin(), top1->resinfo.end());
    fprintf(stderr, "%s\nContaining %zu atoms in %zu residues\n",
            *top1->name, atoms1.size(), resinfo1.size());

    if (bRmpbc)
    {
        rm_gropbc(atoms1, x1, box1);
    }

    fprintf(stderr, "Select group from first structure\n");
    get_index(atoms1, resinfo1, opt2fn_null("-n1", NFILE, fnm),
              1, &isize1, &index1, &groupnames1);
    printf("\n");

    if (bFit && (isize1 < 3))
    {
        gmx_fatal(FARGS, "Need >= 3 points to fit!\n");
    }

    /* reading second structure file */
    fprintf(stderr, "\nReading second structure file\n");
    t_topology *top2 = new t_topology;
    read_tps_conf(conf2file, top2, &ePBC2, &x2, &v2, box2, TRUE);
    std::vector<AtomInfo> atoms2(top2->atoms.begin(), top2->atoms.end());
    std::vector<Residue> resinfo2(top2->resinfo.begin(), top2->resinfo.end());
    fprintf(stderr, "%s\nContaining %zu atoms in %zu residues\n",
            *top2->name, atoms2.size(), resinfo2.size());

    if (bRmpbc)
    {
        rm_gropbc(atoms2, x2, box2);
    }

    fprintf(stderr, "Select group from second structure\n");
    get_index(atoms2, resinfo2, opt2fn_null("-n2", NFILE, fnm),
              1, &isize2, &index2, &groupnames2);

    if (bName)
    {
        find_matching_names(gmx::arrayRefFromArray(index1, isize1), atoms1, resinfo1,
                            gmx::arrayRefFromArray(index2, isize2), atoms2, resinfo2);
        if (matchndxfile)
        {
            fp = gmx_ffopen(matchndxfile, "w");
            fprintf(fp, "; Matching atoms between %s from %s and %s from %s\n",
                    groupnames1, conf1file, groupnames2, conf2file);
            fprintf(fp, "[ Match_%s_%s ]\n", conf1file, groupnames1);
            for (i = 0; i < isize1; i++)
            {
                fprintf(fp, "%4d%s", index1[i]+1, (i%15 == 14 || i == isize1-1) ? "\n" : " ");
            }
            fprintf(fp, "[ Match_%s_%s ]\n", conf2file, groupnames2);
            for (i = 0; i < isize2; i++)
            {
                fprintf(fp, "%4d%s", index2[i]+1, (i%15 == 14 || i == isize2-1) ? "\n" : " ");
            }
        }
    }

    /* check isizes, must be equal */
    if (isize2 != isize1)
    {
        gmx_fatal(FARGS, "You selected groups with differen number of atoms.\n");
    }

    for (i = 0; i < isize1; i++)
    {
        name1 = *atoms1[index1[i]].atomname;
        name2 = *atoms2[index2[i]].atomname;
        if (std::strcmp(name1, name2) != 0)
        {
            if (warn < 20)
            {
                fprintf(stderr,
                        "Warning: atomnames at index %d don't match: %d %s, %d %s\n",
                        i+1, index1[i]+1, name1, index2[i]+1, name2);
            }
            warn++;
        }
        if (!bMW)
        {
            atoms1[index1[i]].m_ = 1;
            atoms2[index2[i]].m_ = 1;
        }
    }
    if (warn)
    {
        fprintf(stderr, "%d atomname%s did not match\n", warn, (warn == 1) ? "" : "s");
    }

    if (bFit)
    {
        /* calculate and remove center of mass of structures */
        calc_rm_cm(gmx::arrayRefFromArray(index1, isize1), atoms1, x1, xcm1);
        calc_rm_cm(gmx::arrayRefFromArray(index2, isize2), atoms2, x2, xcm2);

        snew(w_rls, atoms2.size());
        snew(fit_x, atoms2.size());
        for (at = 0; (at < isize1); at++)
        {
            w_rls[index2[at]] = atoms1[index1[at]].m_;
            copy_rvec(x1[index1[at]], fit_x[index2[at]]);
        }

        /* do the least squares fit to the reference structure */
        do_fit(atoms2.size(), w_rls, fit_x, x2);

        sfree(fit_x);
        sfree(w_rls);
        w_rls = nullptr;
    }
    else
    {
        clear_rvec(xcm1);
        clear_rvec(xcm2);
        w_rls = nullptr;
    }

    /* calculate the rms deviation */
    rms     = 0;
    totmass = 0;
    maxmsd  = -1e18;
    minmsd  =  1e18;
    snew(msds, isize1);
    for (at = 0; at < isize1; at++)
    {
        mass = atoms1[index1[at]].m_;
        for (m = 0; m < DIM; m++)
        {
            msd       = gmx::square(x1[index1[at]][m] - x2[index2[at]][m]);
            rms      += msd*mass;
            msds[at] += msd;
        }
        maxmsd   = std::max(maxmsd, msds[at]);
        minmsd   = std::min(minmsd, msds[at]);
        totmass += mass;
    }
    rms = std::sqrt(rms/totmass);

    printf("Root mean square deviation after lsq fit = %g nm\n", rms);
    if (bBfac)
    {
        printf("Atomic MSD's range from %g to %g nm^2\n", minmsd, maxmsd);
    }

    if (bFit)
    {
        /* reset coordinates of reference and fitted structure */
        for (i = 0; i < gmx::index(atoms1.size()); i++)
        {
            for (m = 0; m < DIM; m++)
            {
                x1[i][m] += xcm1[m];
            }
        }
        for (i = 0; i < gmx::index(atoms2.size()); i++)
        {
            for (m = 0; m < DIM; m++)
            {
                x2[i][m] += xcm1[m];
            }
        }
    }

    outfile = ftp2fn(efSTO, NFILE, fnm);
    std::vector<PdbEntry> pdb1;
    std::vector<PdbEntry> pdb2;
    switch (fn2ftp(outfile))
    {
        case efPDB:
        case efBRK:
        case efENT:
            if (bBfac || bLabel)
            {
                if(!allAtomsHavePdbInfo(top1->pdb) || top1->pdb.size() != top1->atoms.size())
                {
                    pdb1.resize(atoms1.size());
                }
                else
                {
                    pdb1 = std::vector<PdbEntry>(top1->pdb.begin(), top1->pdb.end());
                }

                /* Avoid segfaults when writing the pdb-file */
                for (int i = 0; i < gmx::index(atoms1.size()); i++)
                {
                    pdb1[i].type_         = eptAtom;
                    pdb1[i].occup_        = 1.00;
                    pdb1[i].haveAnisotropic_ = FALSE;
                    pdb1[i].isSet_ = true;
                    if (bBfac)
                    {
                        pdb1[i].bfac_ = 0;
                    }
                    if (bLabel)
                    {
                        resinfo1[atoms1[i].resind_].chainid_ = 'A';
                    }
                }

                for (int i = 0; i < isize1; i++)
                {
                    if (bBfac)
                    {
                        pdb1[index1[i]].bfac_ = (800*M_PI*M_PI/3.0)*msds[i];
                    }
                }
                if(!allAtomsHavePdbInfo(top2->pdb) || top2->pdb.size() != top2->atoms.size())
                {
                    pdb2.resize(atoms2.size());
                }
                else
                {
                    pdb2 = std::vector<PdbEntry>(top2->pdb.begin(), top2->pdb.end());
                }

                for (int i = 0; i < gmx::index(atoms2.size()); i++)
                {
                    pdb2[i].type_         = eptAtom;
                    pdb2[i].occup_        = 1.00;
                    pdb2[i].haveAnisotropic_ = FALSE;
                    pdb2[i].isSet_ = true;
                    if (bBfac)
                    {
                        pdb2[i].bfac_ = 0;
                    }
                    if (bLabel)
                    {
                        resinfo2[atoms1[i].resind_].chainid_ = 'B';
                    }
                }

                for (i = 0; i < isize2; i++)
                {
                    if (bBfac)
                    {
                        pdb2[index2[i]].bfac_ = (800*M_PI*M_PI/3.0)*msds[i];
                    }
                }
            }
            fp = gmx_ffopen(outfile, "w");
            if (!bOne)
            {
                write_pdbfile(fp, *top1->name, atoms1, resinfo1, pdb1, x1, ePBC1, box1, ' ', 1, nullptr, TRUE);
            }
            write_pdbfile(fp, *top2->name, atoms2, resinfo2, pdb2, x2, ePBC2, box2, ' ', bOne ? -1 : 2, nullptr, TRUE);
            gmx_ffclose(fp);
            break;
        case efGRO:
            if (bBfac)
            {
                fprintf(stderr, "WARNING: cannot write B-factor values to gro file\n");
            }
            fp = gmx_ffopen(outfile, "w");
            if (!bOne)
            {
                write_hconf_p(fp, *top1->name, atoms1, resinfo1, x1, v1, box1);
            }
            write_hconf_p(fp, *top2->name, atoms2, resinfo1, x2, v2, box2);
            gmx_ffclose(fp);
            break;
        default:
            if (bBfac)
            {
                fprintf(stderr, "WARNING: cannot write B-factor values to %s file\n",
                        ftp2ext(fn2ftp(outfile)));
            }
            if (!bOne)
            {
                fprintf(stderr,
                        "WARNING: cannot write the reference structure to %s file\n",
                        ftp2ext(fn2ftp(outfile)));
            }
            write_sto_conf(outfile, *top2->name, atoms2, resinfo2, pdb2, x2, v2, ePBC2, box2);
            break;
    }

    view_all(oenv, NFILE, fnm);

    return 0;
}
