/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014, by the GROMACS development team, led by
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
#include <algorithm>
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/legacyheaders/mtop_util.h"
#include "gromacs/legacyheaders/viewit.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/fileio/enxio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/math/vec.h"

#ifdef DEBUG
static real minthird = -1.0/3.0
#endif
static real minsixth = -1.0/6.0;

//static const char       * drleg[] = {
//    "Running average",
//    "Instantaneous"
//};

static double mypow(double x, double y)
{
    if (x > 0)
    {
        return pow(x, y);
    }
    else
    {
        return 0.0;
    }
}

static void get_orires_parms(const char *topnm,
                             int *nor, int *nex, int **label, real **obs)
{
    gmx_mtop_t      mtop;
    gmx_localtop_t *top;
    t_inputrec      ir;
    t_iparams      *ip;
    int             natoms, i;
    t_iatom        *iatom;
    int             nb;
    matrix          box;

    read_tpx(topnm, &ir, box, &natoms, NULL, NULL, NULL, &mtop);
    top = gmx_mtop_generate_local_top(&mtop, &ir);

    ip       = top->idef.iparams;
    iatom    = top->idef.il[F_ORIRES].iatoms;

    /* Count how many distance restraint there are... */
    nb = top->idef.il[F_ORIRES].nr;
    if (nb == 0)
    {
        gmx_fatal(FARGS, "No orientation restraints in topology!\n");
    }

    *nor = nb/3;
    *nex = 0;
    snew(*label, *nor);
    snew(*obs, *nor);
    for (i = 0; i < nb; i += 3)
    {
        (*label)[i/3] = ip[iatom[i]].orires.label;
        (*obs)[i/3]   = ip[iatom[i]].orires.obs;
        if (ip[iatom[i]].orires.ex >= *nex)
        {
            *nex = ip[iatom[i]].orires.ex+1;
        }
    }
    fprintf(stderr, "Found %d orientation restraints with %d experiments",
            *nor, *nex);
}

static int get_bounds(const char *topnm,
                      real **bounds, int **index, int **dr_pair, int *npairs,
                      gmx_mtop_t *mtop, gmx_localtop_t **ltop, t_inputrec *ir)
{
    gmx_localtop_t *top;
    t_functype     *functype;
    t_iparams      *ip;
    int             natoms, i, j, k, type, ftype, natom;
    t_ilist        *disres;
    t_iatom        *iatom;
    real           *b;
    int            *ind, *pair;
    int             nb, label1;
    matrix          box;

    read_tpx(topnm, ir, box, &natoms, NULL, NULL, NULL, mtop);
    snew(*ltop, 1);
    top   = gmx_mtop_generate_local_top(mtop, ir);
    *ltop = top;

    functype = top->idef.functype;
    ip       = top->idef.iparams;

    /* Count how many distance restraint there are... */
    nb = top->idef.il[F_DISRES].nr;
    if (nb == 0)
    {
        gmx_fatal(FARGS, "No distance restraints in topology!\n");
    }

    /* Allocate memory */
    snew(b, nb);
    snew(ind, nb);
    snew(pair, nb+1);

    /* Fill the bound array */
    nb = 0;
    for (i = 0; (i < top->idef.ntypes); i++)
    {
        ftype = functype[i];
        if (ftype == F_DISRES)
        {

            label1  = ip[i].disres.label;
            b[nb]   = ip[i].disres.up1;
            ind[nb] = label1;
            nb++;
        }
    }
    *bounds = b;

    /* Fill the index array */
    label1  = -1;
    disres  = &(top->idef.il[F_DISRES]);
    iatom   = disres->iatoms;
    for (i = j = k = 0; (i < disres->nr); )
    {
        type  = iatom[i];
        ftype = top->idef.functype[type];
        natom = interaction_function[ftype].nratoms+1;
        if (label1 != top->idef.iparams[type].disres.label)
        {
            pair[j] = k;
            label1  = top->idef.iparams[type].disres.label;
            j++;
        }
        k++;
        i += natom;
    }
    pair[j]  = k;
    *npairs  = k;
    if (j != nb)
    {
        gmx_incons("get_bounds for distance restraints");
    }

    *index   = ind;
    *dr_pair = pair;

    return nb;
}

static void calc_violations(real rt[], real rav3[], int nb, int index[],
                            real bounds[], real *viol, double *st, double *sa)
{
    const   real sixth = 1.0/6.0;
    int          i, j;
    double       rsum, rav, sumaver, sumt;

    sumaver = 0;
    sumt    = 0;
    for (i = 0; (i < nb); i++)
    {
        rsum = 0.0;
        rav  = 0.0;
        for (j = index[i]; (j < index[i+1]); j++)
        {
            if (viol)
            {
                viol[j] += mypow(rt[j], -3.0);
            }
            rav     += sqr(rav3[j]);
            rsum    += mypow(rt[j], -6);
        }
        rsum    = std::max(0.0, mypow(rsum, -sixth)-bounds[i]);
        rav     = std::max(0.0, mypow(rav, -sixth)-bounds[i]);

        sumt    += rsum;
        sumaver += rav;
    }
    *st = sumt;
    *sa = sumaver;
}

static void analyse_disre(const char *voutfn,    int nframes,
                          real violaver[], real bounds[], int index[],
                          int pair[],      int nbounds,
                          const output_env_t oenv)
{
    FILE   *vout;
    double  sum, sumt, sumaver;
    int     i, j;

    /* Subtract bounds from distances, to calculate violations */
    calc_violations(violaver, violaver,
                    nbounds, pair, bounds, NULL, &sumt, &sumaver);

#ifdef DEBUG
    fprintf(stdout, "\nSum of violations averaged over simulation: %g nm\n",
            sumaver);
    fprintf(stdout, "Largest violation averaged over simulation: %g nm\n\n",
            sumt);
#endif
    vout = xvgropen(voutfn, "r\\S-3\\N average violations", "DR Index", "nm",
                    oenv);
    sum  = 0.0;
    sumt = 0.0;
    for (i = 0; (i < nbounds); i++)
    {
        /* Do ensemble averaging */
        sumaver = 0;
        for (j = pair[i]; (j < pair[i+1]); j++)
        {
            sumaver += sqr(violaver[j]/nframes);
        }
        sumaver = std::max(0.0, mypow(sumaver, minsixth)-bounds[i]);

        sumt   += sumaver;
        sum     = std::max(sum, sumaver);
        fprintf(vout, "%10d  %10.5e\n", index[i], sumaver);
    }
#ifdef DEBUG
    for (j = 0; (j < dr.ndr); j++)
    {
        fprintf(vout, "%10d  %10.5e\n", j, mypow(violaver[j]/nframes, minthird));
    }
#endif
    gmx_ffclose(vout);

    fprintf(stdout, "\nSum of violations averaged over simulation: %g nm\n",
            sumt);
    fprintf(stdout, "Largest violation averaged over simulation: %g nm\n\n", sum);

    do_view(oenv, voutfn, "-graphtype bar");
}
