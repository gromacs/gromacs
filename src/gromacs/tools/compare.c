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
/* This file is completely threadsafe - keep it that way! */
#include "gmxpre.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

#include "gromacs/fileio/enxio.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/legacyheaders/names.h"
#include "gromacs/legacyheaders/txtdump.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"

static void cmp_int(FILE *fp, const char *s, int index, int i1, int i2)
{
    if (i1 != i2)
    {
        if (index != -1)
        {
            fprintf(fp, "%s[%d] (%d - %d)\n", s, index, i1, i2);
        }
        else
        {
            fprintf(fp, "%s (%d - %d)\n", s, i1, i2);
        }
    }
}

static void cmp_int64(FILE *fp, const char *s, gmx_int64_t i1, gmx_int64_t i2)
{
    if (i1 != i2)
    {
        fprintf(fp, "%s (", s);
        fprintf(fp, "%"GMX_PRId64, i1);
        fprintf(fp, " - ");
        fprintf(fp, "%"GMX_PRId64, i2);
        fprintf(fp, ")\n");
    }
}

static void cmp_us(FILE *fp, const char *s, int index, unsigned short i1, unsigned short i2)
{
    if (i1 != i2)
    {
        if (index != -1)
        {
            fprintf(fp, "%s[%d] (%hu - %hu)\n", s, index, i1, i2);
        }
        else
        {
            fprintf(fp, "%s (%hu - %hu)\n", s, i1, i2);
        }
    }
}

static void cmp_uc(FILE *fp, const char *s, int index, unsigned char i1, unsigned char i2)
{
    if (i1 != i2)
    {
        if (index != -1)
        {
            fprintf(fp, "%s[%d] (%d - %d)\n", s, index, i1, i2);
        }
        else
        {
            fprintf(fp, "%s (%d - %d)\n", s, i1, i2);
        }
    }
}

static gmx_bool cmp_bool(FILE *fp, const char *s, int index, gmx_bool b1, gmx_bool b2)
{
    if (b1)
    {
        b1 = 1;
    }
    else
    {
        b1 = 0;
    }
    if (b2)
    {
        b2 = 1;
    }
    else
    {
        b2 = 0;
    }
    if (b1 != b2)
    {
        if (index != -1)
        {
            fprintf(fp, "%s[%d] (%s - %s)\n", s, index,
                    bool_names[b1], bool_names[b2]);
        }
        else
        {
            fprintf(fp, "%s (%s - %s)\n", s,
                    bool_names[b1], bool_names[b2]);
        }
    }
    return b1 && b2;
}

static void cmp_str(FILE *fp, const char *s, int index,
                    const char *s1, const char *s2)
{
    if (strcmp(s1, s2) != 0)
    {
        if (index != -1)
        {
            fprintf(fp, "%s[%d] (%s - %s)\n", s, index, s1, s2);
        }
        else
        {
            fprintf(fp, "%s (%s - %s)\n", s, s1, s2);
        }
    }
}

static gmx_bool equal_real(real i1, real i2, real ftol, real abstol)
{
    return ( ( 2*fabs(i1 - i2) <= (fabs(i1) + fabs(i2))*ftol ) || fabs(i1-i2) <= abstol );
}

static gmx_bool equal_float(float i1, float i2, float ftol, float abstol)
{
    return ( ( 2*fabs(i1 - i2) <= (fabs(i1) + fabs(i2))*ftol ) || fabs(i1-i2) <= abstol );
}

static gmx_bool equal_double(double i1, double i2, real ftol, real abstol)
{
    return ( ( 2*fabs(i1 - i2) <= (fabs(i1) + fabs(i2))*ftol ) || fabs(i1-i2) <= abstol );
}

static void
cmp_real(FILE *fp, const char *s, int index, real i1, real i2, real ftol, real abstol)
{
    if (!equal_real(i1, i2, ftol, abstol))
    {
        if (index != -1)
        {
            fprintf(fp, "%s[%2d] (%e - %e)\n", s, index, i1, i2);
        }
        else
        {
            fprintf(fp, "%s (%e - %e)\n", s, i1, i2);
        }
    }
}

static void
cmp_float(FILE *fp, const char *s, int index, float i1, float i2, float ftol, float abstol)
{
    if (!equal_float(i1, i2, ftol, abstol))
    {
        if (index != -1)
        {
            fprintf(fp, "%s[%2d] (%e - %e)\n", s, index, i1, i2);
        }
        else
        {
            fprintf(fp, "%s (%e - %e)\n", s, i1, i2);
        }
    }
}



static void
cmp_double(FILE *fp, const char *s, int index, double i1, double i2, double ftol, double abstol)
{
    if (!equal_double(i1, i2, ftol, abstol))
    {
        if (index != -1)
        {
            fprintf(fp, "%s[%2d] (%16.9e - %16.9e)\n", s, index, i1, i2);
        }
        else
        {
            fprintf(fp, "%s (%16.9e - %16.9e)\n", s, i1, i2);
        }
    }
}

static void cmp_rvec(FILE *fp, const char *s, int index, rvec i1, rvec i2, real ftol, real abstol)
{
    if (!equal_real(i1[XX], i2[XX], ftol, abstol) ||
        !equal_real(i1[YY], i2[YY], ftol, abstol) ||
        !equal_real(i1[ZZ], i2[ZZ], ftol, abstol))
    {
        if (index != -1)
        {
            fprintf(fp, "%s[%5d] (%12.5e %12.5e %12.5e) - (%12.5e %12.5e %12.5e)\n",
                    s, index, i1[XX], i1[YY], i1[ZZ], i2[XX], i2[YY], i2[ZZ]);
        }
        else
        {
            fprintf(fp, "%s (%12.5e %12.5e %12.5e) - (%12.5e %12.5e %12.5e)\n",
                    s, i1[XX], i1[YY], i1[ZZ], i2[XX], i2[YY], i2[ZZ]);
        }
    }
}

static void cmp_ivec(FILE *fp, const char *s, int index, ivec i1, ivec i2)
{
    if ((i1[XX] != i2[XX]) || (i1[YY] != i2[YY]) || (i1[ZZ] != i2[ZZ]))
    {
        if (index != -1)
        {
            fprintf(fp, "%s[%5d] (%8d,%8d,%8d - %8d,%8d,%8d)\n", s, index,
                    i1[XX], i1[YY], i1[ZZ], i2[XX], i2[YY], i2[ZZ]);
        }
        else
        {
            fprintf(fp, "%s (%8d,%8d,%8d - %8d,%8d,%8d)\n", s,
                    i1[XX], i1[YY], i1[ZZ], i2[XX], i2[YY], i2[ZZ]);
        }
    }
}

static void cmp_ilist(FILE *fp, int ftype, t_ilist *il1, t_ilist *il2)
{
    int  i;
    char buf[256];

    fprintf(fp, "comparing ilist %s\n", interaction_function[ftype].name);
    sprintf(buf, "%s->nr", interaction_function[ftype].name);
    cmp_int(fp, buf, -1, il1->nr, il2->nr);
    sprintf(buf, "%s->iatoms", interaction_function[ftype].name);
    if (((il1->nr > 0) && (!il1->iatoms)) ||
        ((il2->nr > 0) && (!il2->iatoms)) ||
        ((il1->nr != il2->nr)))
    {
        fprintf(fp, "Comparing radically different topologies - %s is different\n",
                buf);
    }
    else
    {
        for (i = 0; (i < il1->nr); i++)
        {
            cmp_int(fp, buf, i, il1->iatoms[i], il2->iatoms[i]);
        }
    }
}

void cmp_iparm(FILE *fp, const char *s, t_functype ft,
               t_iparams ip1, t_iparams ip2, real ftol, real abstol)
{
    int      i;
    gmx_bool bDiff;

    bDiff = FALSE;
    for (i = 0; i < MAXFORCEPARAM && !bDiff; i++)
    {
        bDiff = !equal_real(ip1.generic.buf[i], ip2.generic.buf[i], ftol, abstol);
    }
    if (bDiff)
    {
        fprintf(fp, "%s1: ", s);
        pr_iparams(fp, ft, &ip1);
        fprintf(fp, "%s2: ", s);
        pr_iparams(fp, ft, &ip2);
    }
}

void cmp_iparm_AB(FILE *fp, const char *s, t_functype ft, t_iparams ip1, real ftol, real abstol)
{
    int      nrfpA, nrfpB, p0, i;
    gmx_bool bDiff;

    /* Normally the first parameter is perturbable */
    p0    = 0;
    nrfpA = interaction_function[ft].nrfpA;
    nrfpB = interaction_function[ft].nrfpB;
    if (ft == F_PDIHS)
    {
        nrfpB = 2;
    }
    else if (interaction_function[ft].flags & IF_TABULATED)
    {
        /* For tabulated interactions only the second parameter is perturbable */
        p0    = 1;
        nrfpB = 1;
    }
    bDiff = FALSE;
    for (i = 0; i < nrfpB && !bDiff; i++)
    {
        bDiff = !equal_real(ip1.generic.buf[p0+i], ip1.generic.buf[nrfpA+i], ftol, abstol);
    }
    if (bDiff)
    {
        fprintf(fp, "%s: ", s);
        pr_iparams(fp, ft, &ip1);
    }
}

static void cmp_idef(FILE *fp, t_idef *id1, t_idef *id2, real ftol, real abstol)
{
    int  i;
    char buf1[64], buf2[64];

    fprintf(fp, "comparing idef\n");
    if (id2)
    {
        cmp_int(fp, "idef->ntypes", -1, id1->ntypes, id2->ntypes);
        cmp_int(fp, "idef->atnr",  -1, id1->atnr, id2->atnr);
        for (i = 0; (i < min(id1->ntypes, id2->ntypes)); i++)
        {
            sprintf(buf1, "idef->functype[%d]", i);
            sprintf(buf2, "idef->iparam[%d]", i);
            cmp_int(fp, buf1, i, (int)id1->functype[i], (int)id2->functype[i]);
            cmp_iparm(fp, buf2, id1->functype[i],
                      id1->iparams[i], id2->iparams[i], ftol, abstol);
        }
        cmp_real(fp, "fudgeQQ", -1, id1->fudgeQQ, id2->fudgeQQ, ftol, abstol);
        for (i = 0; (i < F_NRE); i++)
        {
            cmp_ilist(fp, i, &(id1->il[i]), &(id2->il[i]));
        }
    }
    else
    {
        for (i = 0; (i < id1->ntypes); i++)
        {
            cmp_iparm_AB(fp, "idef->iparam", id1->functype[i], id1->iparams[i], ftol, abstol);
        }
    }
}

static void cmp_block(FILE *fp, t_block *b1, t_block *b2, const char *s)
{
    int  i, j, k;
    char buf[32];

    fprintf(fp, "comparing block %s\n", s);
    sprintf(buf, "%s.nr", s);
    cmp_int(fp, buf, -1, b1->nr, b2->nr);
}

static void cmp_blocka(FILE *fp, t_blocka *b1, t_blocka *b2, const char *s)
{
    int  i, j, k;
    char buf[32];

    fprintf(fp, "comparing blocka %s\n", s);
    sprintf(buf, "%s.nr", s);
    cmp_int(fp, buf, -1, b1->nr, b2->nr);
    sprintf(buf, "%s.nra", s);
    cmp_int(fp, buf, -1, b1->nra, b2->nra);
}

static void cmp_atom(FILE *fp, int index, t_atom *a1, t_atom *a2, real ftol, real abstol)
{
    int  i;
    char buf[256];

    if (a2)
    {
        cmp_us(fp, "atom.type", index, a1->type, a2->type);
        cmp_us(fp, "atom.ptype", index, a1->ptype, a2->ptype);
        cmp_int(fp, "atom.resind", index, a1->resind, a2->resind);
        cmp_int(fp, "atom.atomnumber", index, a1->atomnumber, a2->atomnumber);
        cmp_real(fp, "atom.m", index, a1->m, a2->m, ftol, abstol);
        cmp_real(fp, "atom.q", index, a1->q, a2->q, ftol, abstol);
        cmp_us(fp, "atom.typeB", index, a1->typeB, a2->typeB);
        cmp_real(fp, "atom.mB", index, a1->mB, a2->mB, ftol, abstol);
        cmp_real(fp, "atom.qB", index, a1->qB, a2->qB, ftol, abstol);
    }
    else
    {
        cmp_us(fp, "atom.type", index, a1->type, a1->typeB);
        cmp_real(fp, "atom.m", index, a1->m, a1->mB, ftol, abstol);
        cmp_real(fp, "atom.q", index, a1->q, a1->qB, ftol, abstol);
    }
}

static void cmp_atoms(FILE *fp, t_atoms *a1, t_atoms *a2, real ftol, real abstol)
{
    int i;

    fprintf(fp, "comparing atoms\n");

    if (a2)
    {
        cmp_int(fp, "atoms->nr", -1, a1->nr, a2->nr);
        for (i = 0; (i < a1->nr); i++)
        {
            cmp_atom(fp, i, &(a1->atom[i]), &(a2->atom[i]), ftol, abstol);
        }
    }
    else
    {
        for (i = 0; (i < a1->nr); i++)
        {
            cmp_atom(fp, i, &(a1->atom[i]), NULL, ftol, abstol);
        }
    }
}

static void cmp_top(FILE *fp, t_topology *t1, t_topology *t2, real ftol, real abstol)
{
    int i;

    fprintf(fp, "comparing top\n");
    if (t2)
    {
        cmp_idef(fp, &(t1->idef), &(t2->idef), ftol, abstol);
        cmp_atoms(fp, &(t1->atoms), &(t2->atoms), ftol, abstol);
        cmp_block(fp, &t1->cgs, &t2->cgs, "cgs");
        cmp_block(fp, &t1->mols, &t2->mols, "mols");
        cmp_bool(fp, "bIntermolecularInteractions", -1, t1->bIntermolecularInteractions, t2->bIntermolecularInteractions);
        cmp_blocka(fp, &t1->excls, &t2->excls, "excls");
    }
    else
    {
        cmp_idef(fp, &(t1->idef), NULL, ftol, abstol);
        cmp_atoms(fp, &(t1->atoms), NULL, ftol, abstol);
    }
}

static void cmp_groups(FILE *fp, gmx_groups_t *g0, gmx_groups_t *g1,
                       int natoms0, int natoms1)
{
    int  i, j, ndiff;
    char buf[32];

    fprintf(fp, "comparing groups\n");

    for (i = 0; i < egcNR; i++)
    {
        sprintf(buf, "grps[%d].nr", i);
        cmp_int(fp, buf, -1, g0->grps[i].nr, g1->grps[i].nr);
        if (g0->grps[i].nr == g1->grps[i].nr)
        {
            for (j = 0; j < g0->grps[i].nr; j++)
            {
                sprintf(buf, "grps[%d].name[%d]", i, j);
                cmp_str(fp, buf, -1,
                        *g0->grpname[g0->grps[i].nm_ind[j]],
                        *g1->grpname[g1->grps[i].nm_ind[j]]);
            }
        }
        cmp_int(fp, "ngrpnr", i, g0->ngrpnr[i], g1->ngrpnr[i]);
        if (g0->ngrpnr[i] == g1->ngrpnr[i] && natoms0 == natoms1 &&
            (g0->grpnr[i] != NULL || g1->grpnr[i] != NULL))
        {
            for (j = 0; j < natoms0; j++)
            {
                cmp_int(fp, gtypes[i], j, ggrpnr(g0, i, j), ggrpnr(g1, i, j));
            }
        }
    }
    /* We have compared the names in the groups lists,
     * so we can skip the grpname list comparison.
     */
}

static void cmp_rvecs(FILE *fp, const char *title, int n, rvec x1[], rvec x2[],
                      gmx_bool bRMSD, real ftol, real abstol)
{
    int    i, m;
    double d, ssd;

    if (bRMSD)
    {
        ssd = 0;
        for (i = 0; (i < n); i++)
        {
            for (m = 0; m < DIM; m++)
            {
                d    = x1[i][m] - x2[i][m];
                ssd += d*d;
            }
        }
        fprintf(fp, "%s RMSD %g\n", title, sqrt(ssd/n));
    }
    else
    {
        for (i = 0; (i < n); i++)
        {
            cmp_rvec(fp, title, i, x1[i], x2[i], ftol, abstol);
        }
    }
}


/* Similar to cmp_rvecs, but this routine scales the allowed absolute tolerance
 * by the RMS of the force components of x1.
 */
static void cmp_rvecs_rmstol(FILE *fp, const char *title, int n, rvec x1[], rvec x2[],
                             real ftol, real abstol)
{
    int    i, m;
    double d;
    double ave_x1, rms_x1;

    /* It is tricky to compare real values, in particular forces that
     * are sums of lots of terms where the final value might be close to 0.0.
     * To get a reference magnitude we calculate the RMS value of each
     * component in x1, and then set the allowed absolute tolerance to the
     * relative tolerance times this RMS magnitude.
     */
    ave_x1 = 0.0;
    for (i = 0; i < n; i++)
    {
        for (m = 0; m < DIM; m++)
        {
            ave_x1 += x1[i][m];
        }
    }
    ave_x1 /= n*DIM;

    rms_x1 = 0.0;
    for (i = 0; (i < n); i++)
    {
        for (m = 0; m < DIM; m++)
        {
            d       = x1[i][m] - ave_x1;
            rms_x1 += d*d;
        }
    }
    rms_x1 = sqrt(rms_x1/(DIM*n));
    /* And now do the actual comparision with a hopefully realistic abstol. */
    for (i = 0; (i < n); i++)
    {
        cmp_rvec(fp, title, i, x1[i], x2[i], ftol, abstol*rms_x1);
    }
}

static void cmp_grpopts(FILE *fp, t_grpopts *opt1, t_grpopts *opt2, real ftol, real abstol)
{
    int  i, j;
    char buf1[256], buf2[256];

    cmp_int(fp, "inputrec->grpopts.ngtc", -1,  opt1->ngtc, opt2->ngtc);
    cmp_int(fp, "inputrec->grpopts.ngacc", -1, opt1->ngacc, opt2->ngacc);
    cmp_int(fp, "inputrec->grpopts.ngfrz", -1, opt1->ngfrz, opt2->ngfrz);
    cmp_int(fp, "inputrec->grpopts.ngener", -1, opt1->ngener, opt2->ngener);
    for (i = 0; (i < min(opt1->ngtc, opt2->ngtc)); i++)
    {
        cmp_real(fp, "inputrec->grpopts.nrdf", i, opt1->nrdf[i], opt2->nrdf[i], ftol, abstol);
        cmp_real(fp, "inputrec->grpopts.ref_t", i, opt1->ref_t[i], opt2->ref_t[i], ftol, abstol);
        cmp_real(fp, "inputrec->grpopts.tau_t", i, opt1->tau_t[i], opt2->tau_t[i], ftol, abstol);
        cmp_int(fp, "inputrec->grpopts.annealing", i, opt1->annealing[i], opt2->annealing[i]);
        cmp_int(fp, "inputrec->grpopts.anneal_npoints", i,
                opt1->anneal_npoints[i], opt2->anneal_npoints[i]);
        if (opt1->anneal_npoints[i] == opt2->anneal_npoints[i])
        {
            sprintf(buf1, "inputrec->grpopts.anneal_time[%d]", i);
            sprintf(buf2, "inputrec->grpopts.anneal_temp[%d]", i);
            for (j = 0; j < opt1->anneal_npoints[i]; j++)
            {
                cmp_real(fp, buf1, j, opt1->anneal_time[i][j], opt2->anneal_time[i][j], ftol, abstol);
                cmp_real(fp, buf2, j, opt1->anneal_temp[i][j], opt2->anneal_temp[i][j], ftol, abstol);
            }
        }
    }
    if (opt1->ngener == opt2->ngener)
    {
        for (i = 0; i < opt1->ngener; i++)
        {
            for (j = i; j < opt1->ngener; j++)
            {
                sprintf(buf1, "inputrec->grpopts.egp_flags[%d]", i);
                cmp_int(fp, buf1, j,
                        opt1->egp_flags[opt1->ngener*i+j],
                        opt2->egp_flags[opt1->ngener*i+j]);
            }
        }
    }
    for (i = 0; (i < min(opt1->ngacc, opt2->ngacc)); i++)
    {
        cmp_rvec(fp, "inputrec->grpopts.acc", i, opt1->acc[i], opt2->acc[i], ftol, abstol);
    }
    for (i = 0; (i < min(opt1->ngfrz, opt2->ngfrz)); i++)
    {
        cmp_ivec(fp, "inputrec->grpopts.nFreeze", i, opt1->nFreeze[i], opt2->nFreeze[i]);
    }
}

static void cmp_cosines(FILE *fp, const char *s, t_cosines c1[DIM], t_cosines c2[DIM], real ftol, real abstol)
{
    int  i, m;
    char buf[256];

    for (m = 0; (m < DIM); m++)
    {
        sprintf(buf, "inputrec->%s[%d]", s, m);
        cmp_int(fp, buf, 0, c1->n, c2->n);
        for (i = 0; (i < min(c1->n, c2->n)); i++)
        {
            cmp_real(fp, buf, i, c1->a[i], c2->a[i], ftol, abstol);
            cmp_real(fp, buf, i, c1->phi[i], c2->phi[i], ftol, abstol);
        }
    }
}
static void cmp_adress(FILE *fp, t_adress *ad1, t_adress *ad2,
                       real ftol, real abstol)
{
    cmp_int(fp, "ir->adress->type", -1, ad1->type, ad2->type);
    cmp_real(fp, "ir->adress->const_wf", -1, ad1->const_wf, ad2->const_wf, ftol, abstol);
    cmp_real(fp, "ir->adress->ex_width", -1, ad1->ex_width, ad2->ex_width, ftol, abstol);
    cmp_real(fp, "ir->adress->hy_width", -1, ad1->hy_width, ad2->hy_width, ftol, abstol);
    cmp_int(fp, "ir->adress->icor", -1, ad1->icor, ad2->icor);
    cmp_int(fp, "ir->adress->site", -1, ad1->site, ad2->site);
    cmp_rvec(fp, "ir->adress->refs", -1, ad1->refs, ad2->refs, ftol, abstol);
    cmp_real(fp, "ir->adress->ex_forcecap", -1, ad1->ex_forcecap, ad2->ex_forcecap, ftol, abstol);
}

static void cmp_pull(FILE *fp)
{
    fprintf(fp, "WARNING: Both files use COM pulling, but comparing of the pull struct is not implemented (yet). The pull parameters could be the same or different.\n");
}

static void cmp_simtempvals(FILE *fp, t_simtemp *simtemp1, t_simtemp *simtemp2, int n_lambda, real ftol, real abstol)
{
    int i;
    cmp_int(fp, "inputrec->simtempvals->eSimTempScale", -1, simtemp1->eSimTempScale, simtemp2->eSimTempScale);
    cmp_real(fp, "inputrec->simtempvals->simtemp_high", -1, simtemp1->simtemp_high, simtemp2->simtemp_high, ftol, abstol);
    cmp_real(fp, "inputrec->simtempvals->simtemp_low", -1, simtemp1->simtemp_low, simtemp2->simtemp_low, ftol, abstol);
    for (i = 0; i < n_lambda; i++)
    {
        cmp_real(fp, "inputrec->simtempvals->temperatures", -1, simtemp1->temperatures[i], simtemp2->temperatures[i], ftol, abstol);
    }
}

static void cmp_expandedvals(FILE *fp, t_expanded *expand1, t_expanded *expand2, int n_lambda, real ftol, real abstol)
{
    int i;

    cmp_bool(fp, "inputrec->fepvals->bInit_weights", -1, expand1->bInit_weights, expand2->bInit_weights);
    cmp_bool(fp, "inputrec->fepvals->bWLoneovert", -1, expand1->bWLoneovert, expand2->bWLoneovert);

    for (i = 0; i < n_lambda; i++)
    {
        cmp_real(fp, "inputrec->expandedvals->init_lambda_weights", -1,
                 expand1->init_lambda_weights[i], expand2->init_lambda_weights[i], ftol, abstol);
    }

    cmp_int(fp, "inputrec->expandedvals->lambda-stats", -1, expand1->elamstats, expand2->elamstats);
    cmp_int(fp, "inputrec->expandedvals->lambda-mc-move", -1, expand1->elmcmove, expand2->elmcmove);
    cmp_int(fp, "inputrec->expandedvals->lmc-repeats", -1, expand1->lmc_repeats, expand2->lmc_repeats);
    cmp_int(fp, "inputrec->expandedvals->lmc-gibbsdelta", -1, expand1->gibbsdeltalam, expand2->gibbsdeltalam);
    cmp_int(fp, "inputrec->expandedvals->lmc-forced-nstart", -1, expand1->lmc_forced_nstart, expand2->lmc_forced_nstart);
    cmp_int(fp, "inputrec->expandedvals->lambda-weights-equil", -1, expand1->elmceq, expand2->elmceq);
    cmp_int(fp, "inputrec->expandedvals->,weight-equil-number-all-lambda", -1, expand1->equil_n_at_lam, expand2->equil_n_at_lam);
    cmp_int(fp, "inputrec->expandedvals->weight-equil-number-samples", -1, expand1->equil_samples, expand2->equil_samples);
    cmp_int(fp, "inputrec->expandedvals->weight-equil-number-steps", -1, expand1->equil_steps, expand2->equil_steps);
    cmp_real(fp, "inputrec->expandedvals->weight-equil-wl-delta", -1, expand1->equil_wl_delta, expand2->equil_wl_delta, ftol, abstol);
    cmp_real(fp, "inputrec->expandedvals->weight-equil-count-ratio", -1, expand1->equil_ratio, expand2->equil_ratio, ftol, abstol);
    cmp_bool(fp, "inputrec->expandedvals->symmetrized-transition-matrix", -1, expand1->bSymmetrizedTMatrix, expand2->bSymmetrizedTMatrix);
    cmp_int(fp, "inputrec->expandedvals->nstTij", -1, expand1->nstTij, expand2->nstTij);
    cmp_int(fp, "inputrec->expandedvals->mininum-var-min", -1, expand1->minvarmin, expand2->minvarmin); /*default is reasonable */
    cmp_int(fp, "inputrec->expandedvals->weight-c-range", -1, expand1->c_range, expand2->c_range);      /* default is just C=0 */
    cmp_real(fp, "inputrec->expandedvals->wl-scale", -1, expand1->wl_scale, expand2->wl_scale, ftol, abstol);
    cmp_real(fp, "inputrec->expandedvals->init-wl-delta", -1, expand1->init_wl_delta, expand2->init_wl_delta, ftol, abstol);
    cmp_real(fp, "inputrec->expandedvals->wl-ratio", -1, expand1->wl_ratio, expand2->wl_ratio, ftol, abstol);
    cmp_int(fp, "inputrec->expandedvals->nstexpanded", -1, expand1->nstexpanded, expand2->nstexpanded);
    cmp_int(fp, "inputrec->expandedvals->lmc-seed", -1, expand1->lmc_seed, expand2->lmc_seed);
    cmp_real(fp, "inputrec->expandedvals->mc-temperature", -1, expand1->mc_temp, expand2->mc_temp, ftol, abstol);
}

static void cmp_fepvals(FILE *fp, t_lambda *fep1, t_lambda *fep2, real ftol, real abstol)
{
    int i, j;
    cmp_int(fp, "inputrec->nstdhdl", -1, fep1->nstdhdl, fep2->nstdhdl);
    cmp_double(fp, "inputrec->fepvals->init_fep_state", -1, fep1->init_fep_state, fep2->init_fep_state, ftol, abstol);
    cmp_double(fp, "inputrec->fepvals->delta_lambda", -1, fep1->delta_lambda, fep2->delta_lambda, ftol, abstol);
    cmp_int(fp, "inputrec->fepvals->n_lambda", -1, fep1->n_lambda, fep2->n_lambda);
    for (i = 0; i < efptNR; i++)
    {
        for (j = 0; j < min(fep1->n_lambda, fep2->n_lambda); j++)
        {
            cmp_double(fp, "inputrec->fepvals->all_lambda", -1, fep1->all_lambda[i][j], fep2->all_lambda[i][j], ftol, abstol);
        }
    }
    cmp_int(fp, "inputrec->fepvals->lambda_neighbors", 1, fep1->lambda_neighbors,
            fep2->lambda_neighbors);
    cmp_real(fp, "inputrec->fepvals->sc_alpha", -1, fep1->sc_alpha, fep2->sc_alpha, ftol, abstol);
    cmp_int(fp, "inputrec->fepvals->sc_power", -1, fep1->sc_power, fep2->sc_power);
    cmp_real(fp, "inputrec->fepvals->sc_r_power", -1, fep1->sc_r_power, fep2->sc_r_power, ftol, abstol);
    cmp_real(fp, "inputrec->fepvals->sc_sigma", -1, fep1->sc_sigma, fep2->sc_sigma, ftol, abstol);
    cmp_int(fp, "inputrec->fepvals->edHdLPrintEnergy", -1, fep1->edHdLPrintEnergy, fep1->edHdLPrintEnergy);
    cmp_bool(fp, "inputrec->fepvals->bScCoul", -1, fep1->bScCoul, fep1->bScCoul);
    cmp_int(fp, "inputrec->separate_dhdl_file", -1, fep1->separate_dhdl_file, fep2->separate_dhdl_file);
    cmp_int(fp, "inputrec->dhdl_derivatives", -1, fep1->dhdl_derivatives, fep2->dhdl_derivatives);
    cmp_int(fp, "inputrec->dh_hist_size", -1, fep1->dh_hist_size, fep2->dh_hist_size);
    cmp_double(fp, "inputrec->dh_hist_spacing", -1, fep1->dh_hist_spacing, fep2->dh_hist_spacing, ftol, abstol);
}

static void cmp_inputrec(FILE *fp, t_inputrec *ir1, t_inputrec *ir2, real ftol, real abstol)
{
    fprintf(fp, "comparing inputrec\n");

    /* gcc 2.96 doesnt like these defines at all, but issues a huge list
     * of warnings. Maybe it will change in future versions, but for the
     * moment I've spelled them out instead. /EL 000820
     * #define CIB(s) cmp_int(fp,"inputrec->"#s,0,ir1->##s,ir2->##s)
     * #define CII(s) cmp_int(fp,"inputrec->"#s,0,ir1->##s,ir2->##s)
     * #define CIR(s) cmp_real(fp,"inputrec->"#s,0,ir1->##s,ir2->##s,ftol)
     */
    cmp_int(fp, "inputrec->eI", -1, ir1->eI, ir2->eI);
    cmp_int64(fp, "inputrec->nsteps", ir1->nsteps, ir2->nsteps);
    cmp_int64(fp, "inputrec->init_step", ir1->init_step, ir2->init_step);
    cmp_int(fp, "inputrec->simulation_part", -1, ir1->simulation_part, ir2->simulation_part);
    cmp_int(fp, "inputrec->ePBC", -1, ir1->ePBC, ir2->ePBC);
    cmp_int(fp, "inputrec->bPeriodicMols", -1, ir1->bPeriodicMols, ir2->bPeriodicMols);
    cmp_int(fp, "inputrec->cutoff_scheme", -1, ir1->cutoff_scheme, ir2->cutoff_scheme);
    cmp_int(fp, "inputrec->ns_type", -1, ir1->ns_type, ir2->ns_type);
    cmp_int(fp, "inputrec->nstlist", -1, ir1->nstlist, ir2->nstlist);
    cmp_int(fp, "inputrec->nstcomm", -1, ir1->nstcomm, ir2->nstcomm);
    cmp_int(fp, "inputrec->comm_mode", -1, ir1->comm_mode, ir2->comm_mode);
    cmp_int(fp, "inputrec->nstlog", -1, ir1->nstlog, ir2->nstlog);
    cmp_int(fp, "inputrec->nstxout", -1, ir1->nstxout, ir2->nstxout);
    cmp_int(fp, "inputrec->nstvout", -1, ir1->nstvout, ir2->nstvout);
    cmp_int(fp, "inputrec->nstfout", -1, ir1->nstfout, ir2->nstfout);
    cmp_int(fp, "inputrec->nstcalcenergy", -1, ir1->nstcalcenergy, ir2->nstcalcenergy);
    cmp_int(fp, "inputrec->nstenergy", -1, ir1->nstenergy, ir2->nstenergy);
    cmp_int(fp, "inputrec->nstxout_compressed", -1, ir1->nstxout_compressed, ir2->nstxout_compressed);
    cmp_double(fp, "inputrec->init_t", -1, ir1->init_t, ir2->init_t, ftol, abstol);
    cmp_double(fp, "inputrec->delta_t", -1, ir1->delta_t, ir2->delta_t, ftol, abstol);
    cmp_real(fp, "inputrec->x_compression_precision", -1, ir1->x_compression_precision, ir2->x_compression_precision, ftol, abstol);
    cmp_real(fp, "inputrec->fourierspacing", -1, ir1->fourier_spacing, ir2->fourier_spacing, ftol, abstol);
    cmp_int(fp, "inputrec->nkx", -1, ir1->nkx, ir2->nkx);
    cmp_int(fp, "inputrec->nky", -1, ir1->nky, ir2->nky);
    cmp_int(fp, "inputrec->nkz", -1, ir1->nkz, ir2->nkz);
    cmp_int(fp, "inputrec->pme_order", -1, ir1->pme_order, ir2->pme_order);
    cmp_real(fp, "inputrec->ewald_rtol", -1, ir1->ewald_rtol, ir2->ewald_rtol, ftol, abstol);
    cmp_int(fp, "inputrec->ewald_geometry", -1, ir1->ewald_geometry, ir2->ewald_geometry);
    cmp_real(fp, "inputrec->epsilon_surface", -1, ir1->epsilon_surface, ir2->epsilon_surface, ftol, abstol);
    cmp_int(fp, "inputrec->bContinuation", -1, ir1->bContinuation, ir2->bContinuation);
    cmp_int(fp, "inputrec->bShakeSOR", -1, ir1->bShakeSOR, ir2->bShakeSOR);
    cmp_int(fp, "inputrec->etc", -1, ir1->etc, ir2->etc);
    cmp_int(fp, "inputrec->bPrintNHChains", -1, ir1->bPrintNHChains, ir2->bPrintNHChains);
    cmp_int(fp, "inputrec->epc", -1, ir1->epc, ir2->epc);
    cmp_int(fp, "inputrec->epct", -1, ir1->epct, ir2->epct);
    cmp_real(fp, "inputrec->tau_p", -1, ir1->tau_p, ir2->tau_p, ftol, abstol);
    cmp_rvec(fp, "inputrec->ref_p(x)", -1, ir1->ref_p[XX], ir2->ref_p[XX], ftol, abstol);
    cmp_rvec(fp, "inputrec->ref_p(y)", -1, ir1->ref_p[YY], ir2->ref_p[YY], ftol, abstol);
    cmp_rvec(fp, "inputrec->ref_p(z)", -1, ir1->ref_p[ZZ], ir2->ref_p[ZZ], ftol, abstol);
    cmp_rvec(fp, "inputrec->compress(x)", -1, ir1->compress[XX], ir2->compress[XX], ftol, abstol);
    cmp_rvec(fp, "inputrec->compress(y)", -1, ir1->compress[YY], ir2->compress[YY], ftol, abstol);
    cmp_rvec(fp, "inputrec->compress(z)", -1, ir1->compress[ZZ], ir2->compress[ZZ], ftol, abstol);
    cmp_int(fp, "refcoord_scaling", -1, ir1->refcoord_scaling, ir2->refcoord_scaling);
    cmp_rvec(fp, "inputrec->posres_com", -1, ir1->posres_com, ir2->posres_com, ftol, abstol);
    cmp_rvec(fp, "inputrec->posres_comB", -1, ir1->posres_comB, ir2->posres_comB, ftol, abstol);
    cmp_real(fp, "inputrec->verletbuf_tol", -1, ir1->verletbuf_tol, ir2->verletbuf_tol, ftol, abstol);
    cmp_real(fp, "inputrec->rlist", -1, ir1->rlist, ir2->rlist, ftol, abstol);
    cmp_real(fp, "inputrec->rlistlong", -1, ir1->rlistlong, ir2->rlistlong, ftol, abstol);
    cmp_int(fp, "inputrec->nstcalclr", -1, ir1->nstcalclr, ir2->nstcalclr);
    cmp_real(fp, "inputrec->rtpi", -1, ir1->rtpi, ir2->rtpi, ftol, abstol);
    cmp_int(fp, "inputrec->coulombtype", -1, ir1->coulombtype, ir2->coulombtype);
    cmp_int(fp, "inputrec->coulomb_modifier", -1, ir1->coulomb_modifier, ir2->coulomb_modifier);
    cmp_real(fp, "inputrec->rcoulomb_switch", -1, ir1->rcoulomb_switch, ir2->rcoulomb_switch, ftol, abstol);
    cmp_real(fp, "inputrec->rcoulomb", -1, ir1->rcoulomb, ir2->rcoulomb, ftol, abstol);
    cmp_int(fp, "inputrec->vdwtype", -1, ir1->vdwtype, ir2->vdwtype);
    cmp_int(fp, "inputrec->vdw_modifier", -1, ir1->vdw_modifier, ir2->vdw_modifier);  cmp_real(fp, "inputrec->rvdw_switch", -1, ir1->rvdw_switch, ir2->rvdw_switch, ftol, abstol);
    cmp_real(fp, "inputrec->rvdw", -1, ir1->rvdw, ir2->rvdw, ftol, abstol);
    cmp_real(fp, "inputrec->epsilon_r", -1, ir1->epsilon_r, ir2->epsilon_r, ftol, abstol);
    cmp_real(fp, "inputrec->epsilon_rf", -1, ir1->epsilon_rf, ir2->epsilon_rf, ftol, abstol);
    cmp_real(fp, "inputrec->tabext", -1, ir1->tabext, ir2->tabext, ftol, abstol);
    cmp_int(fp, "inputrec->implicit_solvent", -1, ir1->implicit_solvent, ir2->implicit_solvent);
    cmp_int(fp, "inputrec->gb_algorithm", -1, ir1->gb_algorithm, ir2->gb_algorithm);
    cmp_int(fp, "inputrec->nstgbradii", -1, ir1->nstgbradii, ir2->nstgbradii);
    cmp_real(fp, "inputrec->rgbradii", -1, ir1->rgbradii, ir2->rgbradii, ftol, abstol);
    cmp_real(fp, "inputrec->gb_saltconc", -1, ir1->gb_saltconc, ir2->gb_saltconc, ftol, abstol);
    cmp_real(fp, "inputrec->gb_epsilon_solvent", -1, ir1->gb_epsilon_solvent, ir2->gb_epsilon_solvent, ftol, abstol);
    cmp_real(fp, "inputrec->gb_obc_alpha", -1, ir1->gb_obc_alpha, ir2->gb_obc_alpha, ftol, abstol);
    cmp_real(fp, "inputrec->gb_obc_beta", -1, ir1->gb_obc_beta, ir2->gb_obc_beta, ftol, abstol);
    cmp_real(fp, "inputrec->gb_obc_gamma", -1, ir1->gb_obc_gamma, ir2->gb_obc_gamma, ftol, abstol);
    cmp_real(fp, "inputrec->gb_dielectric_offset", -1, ir1->gb_dielectric_offset, ir2->gb_dielectric_offset, ftol, abstol);
    cmp_int(fp, "inputrec->sa_algorithm", -1, ir1->sa_algorithm, ir2->sa_algorithm);
    cmp_real(fp, "inputrec->sa_surface_tension", -1, ir1->sa_surface_tension, ir2->sa_surface_tension, ftol, abstol);

    cmp_int(fp, "inputrec->eDispCorr", -1, ir1->eDispCorr, ir2->eDispCorr);
    cmp_real(fp, "inputrec->shake_tol", -1, ir1->shake_tol, ir2->shake_tol, ftol, abstol);
    cmp_int(fp, "inputrec->efep", -1, ir1->efep, ir2->efep);
    cmp_fepvals(fp, ir1->fepvals, ir2->fepvals, ftol, abstol);
    cmp_int(fp, "inputrec->bSimTemp", -1, ir1->bSimTemp, ir2->bSimTemp);
    if ((ir1->bSimTemp == ir2->bSimTemp) && (ir1->bSimTemp))
    {
        cmp_simtempvals(fp, ir1->simtempvals, ir2->simtempvals, min(ir1->fepvals->n_lambda, ir2->fepvals->n_lambda), ftol, abstol);
    }
    cmp_int(fp, "inputrec->bExpanded", -1, ir1->bExpanded, ir2->bExpanded);
    if ((ir1->bExpanded == ir2->bExpanded) && (ir1->bExpanded))
    {
        cmp_expandedvals(fp, ir1->expandedvals, ir2->expandedvals, min(ir1->fepvals->n_lambda, ir2->fepvals->n_lambda), ftol, abstol);
    }
    cmp_int(fp, "inputrec->nwall", -1, ir1->nwall, ir2->nwall);
    cmp_int(fp, "inputrec->wall_type", -1, ir1->wall_type, ir2->wall_type);
    cmp_int(fp, "inputrec->wall_atomtype[0]", -1, ir1->wall_atomtype[0], ir2->wall_atomtype[0]);
    cmp_int(fp, "inputrec->wall_atomtype[1]", -1, ir1->wall_atomtype[1], ir2->wall_atomtype[1]);
    cmp_real(fp, "inputrec->wall_density[0]", -1, ir1->wall_density[0], ir2->wall_density[0], ftol, abstol);
    cmp_real(fp, "inputrec->wall_density[1]", -1, ir1->wall_density[1], ir2->wall_density[1], ftol, abstol);
    cmp_real(fp, "inputrec->wall_ewald_zfac", -1, ir1->wall_ewald_zfac, ir2->wall_ewald_zfac, ftol, abstol);

    cmp_bool(fp, "inputrec->bPull", -1, ir1->bPull, ir2->bPull);
    if (ir1->bPull && ir2->bPull)
    {
        cmp_pull(fp);
    }

    cmp_int(fp, "inputrec->eDisre", -1, ir1->eDisre, ir2->eDisre);
    cmp_real(fp, "inputrec->dr_fc", -1, ir1->dr_fc, ir2->dr_fc, ftol, abstol);
    cmp_int(fp, "inputrec->eDisreWeighting", -1, ir1->eDisreWeighting, ir2->eDisreWeighting);
    cmp_int(fp, "inputrec->bDisreMixed", -1, ir1->bDisreMixed, ir2->bDisreMixed);
    cmp_int(fp, "inputrec->nstdisreout", -1, ir1->nstdisreout, ir2->nstdisreout);
    cmp_real(fp, "inputrec->dr_tau", -1, ir1->dr_tau, ir2->dr_tau, ftol, abstol);
    cmp_real(fp, "inputrec->orires_fc", -1, ir1->orires_fc, ir2->orires_fc, ftol, abstol);
    cmp_real(fp, "inputrec->orires_tau", -1, ir1->orires_tau, ir2->orires_tau, ftol, abstol);
    cmp_int(fp, "inputrec->nstorireout", -1, ir1->nstorireout, ir2->nstorireout);
    cmp_real(fp, "inputrec->em_stepsize", -1, ir1->em_stepsize, ir2->em_stepsize, ftol, abstol);
    cmp_real(fp, "inputrec->em_tol", -1, ir1->em_tol, ir2->em_tol, ftol, abstol);
    cmp_int(fp, "inputrec->niter", -1, ir1->niter, ir2->niter);
    cmp_real(fp, "inputrec->fc_stepsize", -1, ir1->fc_stepsize, ir2->fc_stepsize, ftol, abstol);
    cmp_int(fp, "inputrec->nstcgsteep", -1, ir1->nstcgsteep, ir2->nstcgsteep);
    cmp_int(fp, "inputrec->nbfgscorr", 0, ir1->nbfgscorr, ir2->nbfgscorr);
    cmp_int(fp, "inputrec->eConstrAlg", -1, ir1->eConstrAlg, ir2->eConstrAlg);
    cmp_int(fp, "inputrec->nProjOrder", -1, ir1->nProjOrder, ir2->nProjOrder);
    cmp_real(fp, "inputrec->LincsWarnAngle", -1, ir1->LincsWarnAngle, ir2->LincsWarnAngle, ftol, abstol);
    cmp_int(fp, "inputrec->nLincsIter", -1, ir1->nLincsIter, ir2->nLincsIter);
    cmp_real(fp, "inputrec->bd_fric", -1, ir1->bd_fric, ir2->bd_fric, ftol, abstol);
    cmp_int64(fp, "inputrec->ld_seed", ir1->ld_seed, ir2->ld_seed);
    cmp_real(fp, "inputrec->cos_accel", -1, ir1->cos_accel, ir2->cos_accel, ftol, abstol);
    cmp_rvec(fp, "inputrec->deform(a)", -1, ir1->deform[XX], ir2->deform[XX], ftol, abstol);
    cmp_rvec(fp, "inputrec->deform(b)", -1, ir1->deform[YY], ir2->deform[YY], ftol, abstol);
    cmp_rvec(fp, "inputrec->deform(c)", -1, ir1->deform[ZZ], ir2->deform[ZZ], ftol, abstol);


    cmp_bool(fp, "ir->bAdress->type", -1, ir1->bAdress, ir2->bAdress);
    if (ir1->bAdress && ir2->bAdress)
    {
        cmp_adress(fp, ir1->adress, ir2->adress, ftol, abstol);
    }

    cmp_int(fp, "inputrec->userint1", -1, ir1->userint1, ir2->userint1);
    cmp_int(fp, "inputrec->userint2", -1, ir1->userint2, ir2->userint2);
    cmp_int(fp, "inputrec->userint3", -1, ir1->userint3, ir2->userint3);
    cmp_int(fp, "inputrec->userint4", -1, ir1->userint4, ir2->userint4);
    cmp_real(fp, "inputrec->userreal1", -1, ir1->userreal1, ir2->userreal1, ftol, abstol);
    cmp_real(fp, "inputrec->userreal2", -1, ir1->userreal2, ir2->userreal2, ftol, abstol);
    cmp_real(fp, "inputrec->userreal3", -1, ir1->userreal3, ir2->userreal3, ftol, abstol);
    cmp_real(fp, "inputrec->userreal4", -1, ir1->userreal4, ir2->userreal4, ftol, abstol);
    cmp_grpopts(fp, &(ir1->opts), &(ir2->opts), ftol, abstol);
    cmp_cosines(fp, "ex", ir1->ex, ir2->ex, ftol, abstol);
    cmp_cosines(fp, "et", ir1->et, ir2->et, ftol, abstol);
}

static void comp_pull_AB(FILE *fp, pull_params_t *pull, real ftol, real abstol)
{
    int i;

    for (i = 0; i < pull->ncoord; i++)
    {
        fprintf(fp, "comparing pull coord %d\n", i);
        cmp_real(fp, "pull-coord->k", -1, pull->coord[i].k, pull->coord[i].kB, ftol, abstol);
    }
}

static void comp_state(t_state *st1, t_state *st2,
                       gmx_bool bRMSD, real ftol, real abstol)
{
    int i, j, nc;

    fprintf(stdout, "comparing flags\n");
    cmp_int(stdout, "flags", -1, st1->flags, st2->flags);
    fprintf(stdout, "comparing box\n");
    cmp_rvecs(stdout, "box", DIM, st1->box, st2->box, FALSE, ftol, abstol);
    fprintf(stdout, "comparing box_rel\n");
    cmp_rvecs(stdout, "box_rel", DIM, st1->box_rel, st2->box_rel, FALSE, ftol, abstol);
    fprintf(stdout, "comparing boxv\n");
    cmp_rvecs(stdout, "boxv", DIM, st1->boxv, st2->boxv, FALSE, ftol, abstol);
    if (st1->flags & (1<<estSVIR_PREV))
    {
        fprintf(stdout, "comparing shake vir_prev\n");
        cmp_rvecs_rmstol(stdout, "svir_prev", DIM, st1->svir_prev, st2->svir_prev, ftol, abstol);
    }
    if (st1->flags & (1<<estFVIR_PREV))
    {
        fprintf(stdout, "comparing force vir_prev\n");
        cmp_rvecs_rmstol(stdout, "fvir_prev", DIM, st1->fvir_prev, st2->fvir_prev, ftol, abstol);
    }
    if (st1->flags & (1<<estPRES_PREV))
    {
        fprintf(stdout, "comparing prev_pres\n");
        cmp_rvecs_rmstol(stdout, "pres_prev", DIM, st1->pres_prev, st2->pres_prev, ftol, abstol);
    }
    cmp_int(stdout, "ngtc", -1, st1->ngtc, st2->ngtc);
    cmp_int(stdout, "nhchainlength", -1, st1->nhchainlength, st2->nhchainlength);
    if (st1->ngtc == st2->ngtc && st1->nhchainlength == st2->nhchainlength)
    {
        for (i = 0; i < st1->ngtc; i++)
        {
            nc = i*st1->nhchainlength;
            for (j = 0; j < nc; j++)
            {
                cmp_real(stdout, "nosehoover_xi",
                         i, st1->nosehoover_xi[nc+j], st2->nosehoover_xi[nc+j], ftol, abstol);
            }
        }
    }
    cmp_int(stdout, "nnhpres", -1, st1->nnhpres, st2->nnhpres);
    if (st1->nnhpres == st2->nnhpres && st1->nhchainlength == st2->nhchainlength)
    {
        for (i = 0; i < st1->nnhpres; i++)
        {
            nc = i*st1->nhchainlength;
            for (j = 0; j < nc; j++)
            {
                cmp_real(stdout, "nosehoover_xi",
                         i, st1->nhpres_xi[nc+j], st2->nhpres_xi[nc+j], ftol, abstol);
            }
        }
    }

    cmp_int(stdout, "natoms", -1, st1->natoms, st2->natoms);
    if (st1->natoms == st2->natoms)
    {
        if ((st1->flags & (1<<estX)) && (st2->flags & (1<<estX)))
        {
            fprintf(stdout, "comparing x\n");
            cmp_rvecs(stdout, "x", st1->natoms, st1->x, st2->x, bRMSD, ftol, abstol);
        }
        if ((st1->flags & (1<<estV)) && (st2->flags & (1<<estV)))
        {
            fprintf(stdout, "comparing v\n");
            cmp_rvecs(stdout, "v", st1->natoms, st1->v, st2->v, bRMSD, ftol, abstol);
        }
    }
}

void comp_tpx(const char *fn1, const char *fn2,
              gmx_bool bRMSD, real ftol, real abstol)
{
    const char  *ff[2];
    t_tpxheader  sh[2];
    t_inputrec   ir[2];
    t_state      state[2];
    gmx_mtop_t   mtop[2];
    t_topology   top[2];
    int          i;

    ff[0] = fn1;
    ff[1] = fn2;
    for (i = 0; i < (fn2 ? 2 : 1); i++)
    {
        read_tpx_state(ff[i], &(ir[i]), &state[i], NULL, &(mtop[i]));
    }
    if (fn2)
    {
        cmp_inputrec(stdout, &ir[0], &ir[1], ftol, abstol);
        /* Convert gmx_mtop_t to t_topology.
         * We should implement direct mtop comparison,
         * but it might be useful to keep t_topology comparison as an option.
         */
        top[0] = gmx_mtop_t_to_t_topology(&mtop[0]);
        top[1] = gmx_mtop_t_to_t_topology(&mtop[1]);
        cmp_top(stdout, &top[0], &top[1], ftol, abstol);
        cmp_groups(stdout, &mtop[0].groups, &mtop[1].groups,
                   mtop[0].natoms, mtop[1].natoms);
        comp_state(&state[0], &state[1], bRMSD, ftol, abstol);
    }
    else
    {
        if (ir[0].efep == efepNO)
        {
            fprintf(stdout, "inputrec->efep = %s\n", efep_names[ir[0].efep]);
        }
        else
        {
            if (ir[0].bPull)
            {
                comp_pull_AB(stdout, ir->pull, ftol, abstol);
            }
            /* Convert gmx_mtop_t to t_topology.
             * We should implement direct mtop comparison,
             * but it might be useful to keep t_topology comparison as an option.
             */
            top[0] = gmx_mtop_t_to_t_topology(&mtop[0]);
            cmp_top(stdout, &top[0], NULL, ftol, abstol);
        }
    }
}

void comp_frame(FILE *fp, t_trxframe *fr1, t_trxframe *fr2,
                gmx_bool bRMSD, real ftol, real abstol)
{
    fprintf(fp, "\n");
    cmp_int(fp, "flags", -1, fr1->flags, fr2->flags);
    cmp_int(fp, "not_ok", -1, fr1->not_ok, fr2->not_ok);
    cmp_int(fp, "natoms", -1, fr1->natoms, fr2->natoms);
    cmp_real(fp, "t0", -1, fr1->t0, fr2->t0, ftol, abstol);
    if (cmp_bool(fp, "bTitle", -1, fr1->bTitle, fr2->bTitle))
    {
        cmp_str(fp, "title", -1, fr1->title, fr2->title);
    }
    if (cmp_bool(fp, "bStep", -1, fr1->bStep, fr2->bStep))
    {
        cmp_int(fp, "step", -1, fr1->step, fr2->step);
    }
    cmp_int(fp, "step", -1, fr1->step, fr2->step);
    if (cmp_bool(fp, "bTime", -1, fr1->bTime, fr2->bTime))
    {
        cmp_real(fp, "time", -1, fr1->time, fr2->time, ftol, abstol);
    }
    if (cmp_bool(fp, "bLambda", -1, fr1->bLambda, fr2->bLambda))
    {
        cmp_real(fp, "lambda", -1, fr1->lambda, fr2->lambda, ftol, abstol);
    }
    if (cmp_bool(fp, "bAtoms", -1, fr1->bAtoms, fr2->bAtoms))
    {
        cmp_atoms(fp, fr1->atoms, fr2->atoms, ftol, abstol);
    }
    if (cmp_bool(fp, "bPrec", -1, fr1->bPrec, fr2->bPrec))
    {
        cmp_real(fp, "prec", -1, fr1->prec, fr2->prec, ftol, abstol);
    }
    if (cmp_bool(fp, "bX", -1, fr1->bX, fr2->bX))
    {
        cmp_rvecs(fp, "x", min(fr1->natoms, fr2->natoms), fr1->x, fr2->x, bRMSD, ftol, abstol);
    }
    if (cmp_bool(fp, "bV", -1, fr1->bV, fr2->bV))
    {
        cmp_rvecs(fp, "v", min(fr1->natoms, fr2->natoms), fr1->v, fr2->v, bRMSD, ftol, abstol);
    }
    if (cmp_bool(fp, "bF", -1, fr1->bF, fr2->bF))
    {
        if (bRMSD)
        {
            cmp_rvecs(fp, "f", min(fr1->natoms, fr2->natoms), fr1->f, fr2->f, bRMSD, ftol, abstol);
        }
        else
        {
            cmp_rvecs_rmstol(fp, "f", min(fr1->natoms, fr2->natoms), fr1->f, fr2->f, ftol, abstol);
        }
    }
    if (cmp_bool(fp, "bBox", -1, fr1->bBox, fr2->bBox))
    {
        cmp_rvecs(fp, "box", 3, fr1->box, fr2->box, FALSE, ftol, abstol);
    }
}

void comp_trx(const output_env_t oenv, const char *fn1, const char *fn2,
              gmx_bool bRMSD, real ftol, real abstol)
{
    int          i;
    const char  *fn[2];
    t_trxframe   fr[2];
    t_trxstatus *status[2];
    gmx_bool     b[2];

    fn[0] = fn1;
    fn[1] = fn2;
    fprintf(stderr, "Comparing trajectory files %s and %s\n", fn1, fn2);
    for (i = 0; i < 2; i++)
    {
        b[i] = read_first_frame(oenv, &status[i], fn[i], &fr[i], TRX_READ_X|TRX_READ_V|TRX_READ_F);
    }

    if (b[0] && b[1])
    {
        do
        {
            comp_frame(stdout, &(fr[0]), &(fr[1]), bRMSD, ftol, abstol);

            for (i = 0; i < 2; i++)
            {
                b[i] = read_next_frame(oenv, status[i], &fr[i]);
            }
        }
        while (b[0] && b[1]);

        for (i = 0; i < 2; i++)
        {
            if (b[i] && !b[1-i])
            {
                fprintf(stdout, "\nEnd of file on %s but not on %s\n", fn[1-i], fn[i]);
            }
            close_trj(status[i]);
        }
    }
    if (!b[0] && !b[1])
    {
        fprintf(stdout, "\nBoth files read correctly\n");
    }
}

static real ener_tensor_diag(int n, int *ind1, int *ind2,
                             gmx_enxnm_t *enm1,
                             int *tensi, int i,
                             t_energy e1[], t_energy e2[])
{
    int  d1, d2;
    int  len;
    int  j;
    real prod1, prod2;
    int  nfound;

    d1 = tensi[i]/DIM;
    d2 = tensi[i] - d1*DIM;

    /* Find the diagonal elements d1 and d2 */
    len    = strlen(enm1[ind1[i]].name);
    prod1  = 1;
    prod2  = 1;
    nfound = 0;
    for (j = 0; j < n; j++)
    {
        if (tensi[j] >= 0 &&
            strlen(enm1[ind1[j]].name) == len &&
            strncmp(enm1[ind1[i]].name, enm1[ind1[j]].name, len-2) == 0 &&
            (tensi[j] == d1*DIM+d1 || tensi[j] == d2*DIM+d2))
        {
            prod1 *= fabs(e1[ind1[j]].e);
            prod2 *= fabs(e2[ind2[j]].e);
            nfound++;
        }
    }

    if (nfound == 2)
    {
        return 0.5*(sqrt(prod1) + sqrt(prod2));
    }
    else
    {
        return 0;
    }
}

static gmx_bool enernm_equal(const char *nm1, const char *nm2)
{
    int len1, len2;

    len1 = strlen(nm1);
    len2 = strlen(nm2);

    /* Remove " (bar)" at the end of a name */
    if (len1 > 6 && strcmp(nm1+len1-6, " (bar)") == 0)
    {
        len1 -= 6;
    }
    if (len2 > 6 && strcmp(nm2+len2-6, " (bar)") == 0)
    {
        len2 -= 6;
    }

    return (len1 == len2 && gmx_strncasecmp(nm1, nm2, len1) == 0);
}

static void cmp_energies(FILE *fp, int step1, int step2,
                         t_energy e1[], t_energy e2[],
                         gmx_enxnm_t *enm1,
                         real ftol, real abstol,
                         int nre, int *ind1, int *ind2, int maxener)
{
    int   i, ii;
    int  *tensi, len, d1, d2;
    real  ftol_i, abstol_i;

    snew(tensi, maxener);
    /* Check for tensor elements ending on "-XX", "-XY", ... , "-ZZ" */
    for (i = 0; (i < maxener); i++)
    {
        ii       = ind1[i];
        tensi[i] = -1;
        len      = strlen(enm1[ii].name);
        if (len > 3 && enm1[ii].name[len-3] == '-')
        {
            d1 = enm1[ii].name[len-2] - 'X';
            d2 = enm1[ii].name[len-1] - 'X';
            if (d1 >= 0 && d1 < DIM &&
                d2 >= 0 && d2 < DIM)
            {
                tensi[i] = d1*DIM + d2;
            }
        }
    }

    for (i = 0; (i < maxener); i++)
    {
        /* Check if this is an off-diagonal tensor element */
        if (tensi[i] >= 0 && tensi[i] != 0 && tensi[i] != 4 && tensi[i] != 8)
        {
            /* Turn on the relative tolerance check (4 is maximum relative diff.) */
            ftol_i = 5;
            /* Do the relative tolerance through an absolute tolerance times
             * the size of diagonal components of the tensor.
             */
            abstol_i = ftol*ener_tensor_diag(nre, ind1, ind2, enm1, tensi, i, e1, e2);
            if (debug)
            {
                fprintf(debug, "tensor '%s' val %f diag %f\n",
                        enm1[i].name, e1[i].e, abstol_i/ftol);
            }
            if (abstol_i > 0)
            {
                /* We found a diagonal, we need to check with the minimum tolerance */
                abstol_i = min(abstol_i, abstol);
            }
            else
            {
                /* We did not find a diagonal, ignore the relative tolerance check */
                abstol_i = abstol;
            }
        }
        else
        {
            ftol_i   = ftol;
            abstol_i = abstol;
        }
        if (!equal_real(e1[ind1[i]].e, e2[ind2[i]].e, ftol_i, abstol_i))
        {
            fprintf(fp, "%-15s  step %3d:  %12g,  step %3d: %12g\n",
                    enm1[ind1[i]].name,
                    step1, e1[ind1[i]].e,
                    step2, e2[ind2[i]].e);
        }
    }

    sfree(tensi);
}

#if 0
static void cmp_disres(t_enxframe *fr1, t_enxframe *fr2, real ftol, real abstol)
{
    int  i;
    char bav[64], bt[64], bs[22];

    cmp_int(stdout, "ndisre", -1, fr1->ndisre, fr2->ndisre);
    if ((fr1->ndisre == fr2->ndisre) && (fr1->ndisre > 0))
    {
        sprintf(bav, "step %s: disre rav", gmx_step_str(fr1->step, bs));
        sprintf(bt, "step %s: disre  rt", gmx_step_str(fr1->step, bs));
        for (i = 0; (i < fr1->ndisre); i++)
        {
            cmp_real(stdout, bav, i, fr1->disre_rm3tav[i], fr2->disre_rm3tav[i], ftol, abstol);
            cmp_real(stdout, bt, i, fr1->disre_rt[i], fr2->disre_rt[i], ftol, abstol);
        }
    }
}
#endif

static void cmp_eblocks(t_enxframe *fr1, t_enxframe *fr2, real ftol, real abstol)
{
    int  i, j, k;
    char buf[64], bs[22];

    cmp_int(stdout, "nblock", -1, fr1->nblock, fr2->nblock);
    if ((fr1->nblock == fr2->nblock) && (fr1->nblock > 0))
    {
        for (j = 0; (j < fr1->nblock); j++)
        {
            t_enxblock *b1, *b2; /* convenience vars */

            b1 = &(fr1->block[j]);
            b2 = &(fr2->block[j]);

            sprintf(buf, "step %s: block[%d]", gmx_step_str(fr1->step, bs), j);
            cmp_int(stdout, buf, -1, b1->nsub, b2->nsub);
            cmp_int(stdout, buf, -1, b1->id, b2->id);

            if ( (b1->nsub == b2->nsub) && (b1->id == b2->id) )
            {
                for (i = 0; i < b1->nsub; i++)
                {
                    t_enxsubblock *s1, *s2;

                    s1 = &(b1->sub[i]);
                    s2 = &(b2->sub[i]);

                    cmp_int(stdout, buf, -1, (int)s1->type, (int)s2->type);
                    cmp_int64(stdout, buf, s1->nr, s2->nr);

                    if ((s1->type == s2->type) && (s1->nr == s2->nr))
                    {
                        switch (s1->type)
                        {
                            case xdr_datatype_float:
                                for (k = 0; k < s1->nr; k++)
                                {
                                    cmp_float(stdout, buf, i,
                                              s1->fval[k], s2->fval[k],
                                              ftol, abstol);
                                }
                                break;
                            case xdr_datatype_double:
                                for (k = 0; k < s1->nr; k++)
                                {
                                    cmp_double(stdout, buf, i,
                                               s1->dval[k], s2->dval[k],
                                               ftol, abstol);
                                }
                                break;
                            case xdr_datatype_int:
                                for (k = 0; k < s1->nr; k++)
                                {
                                    cmp_int(stdout, buf, i,
                                            s1->ival[k], s2->ival[k]);
                                }
                                break;
                            case xdr_datatype_int64:
                                for (k = 0; k < s1->nr; k++)
                                {
                                    cmp_int64(stdout, buf,
                                              s1->lval[k], s2->lval[k]);
                                }
                                break;
                            case xdr_datatype_char:
                                for (k = 0; k < s1->nr; k++)
                                {
                                    cmp_uc(stdout, buf, i,
                                           s1->cval[k], s2->cval[k]);
                                }
                                break;
                            case xdr_datatype_string:
                                for (k = 0; k < s1->nr; k++)
                                {
                                    cmp_str(stdout, buf, i,
                                            s1->sval[k], s2->sval[k]);
                                }
                                break;
                            default:
                                gmx_incons("Unknown data type!!");
                        }
                    }
                }
            }
        }
    }
}

void comp_enx(const char *fn1, const char *fn2, real ftol, real abstol, const char *lastener)
{
    int            nre, nre1, nre2, block;
    ener_file_t    in1, in2;
    int            i, j, maxener, *ind1, *ind2, *have;
    char           buf[256];
    gmx_enxnm_t   *enm1 = NULL, *enm2 = NULL;
    t_enxframe    *fr1, *fr2;
    gmx_bool       b1, b2;

    fprintf(stdout, "comparing energy file %s and %s\n\n", fn1, fn2);

    in1 = open_enx(fn1, "r");
    in2 = open_enx(fn2, "r");
    do_enxnms(in1, &nre1, &enm1);
    do_enxnms(in2, &nre2, &enm2);
    if (nre1 != nre2)
    {
        fprintf(stdout, "There are %d and %d terms in the energy files\n\n",
                nre1, nre2);
    }
    else
    {
        fprintf(stdout, "There are %d terms in the energy files\n\n", nre1);
    }

    snew(ind1, nre1);
    snew(ind2, nre2);
    snew(have, nre2);
    nre = 0;
    for (i = 0; i < nre1; i++)
    {
        for (j = 0; j < nre2; j++)
        {
            if (enernm_equal(enm1[i].name, enm2[j].name))
            {
                ind1[nre] = i;
                ind2[nre] = j;
                have[j]   = 1;
                nre++;
                break;
            }
        }
        if (nre == 0 || ind1[nre-1] != i)
        {
            cmp_str(stdout, "enm", i, enm1[i].name, "-");
        }
    }
    for (i = 0; i < nre2; i++)
    {
        if (have[i] == 0)
        {
            cmp_str(stdout, "enm", i, "-", enm2[i].name);
        }
    }

    maxener = nre;
    for (i = 0; i < nre; i++)
    {
        if ((lastener != NULL) && (strstr(enm1[i].name, lastener) != NULL))
        {
            maxener = i+1;
            break;
        }
    }

    fprintf(stdout, "There are %d terms to compare in the energy files\n\n",
            maxener);

    for (i = 0; i < maxener; i++)
    {
        cmp_str(stdout, "unit", i, enm1[ind1[i]].unit, enm2[ind2[i]].unit);
    }

    snew(fr1, 1);
    snew(fr2, 1);
    do
    {
        b1 = do_enx(in1, fr1);
        b2 = do_enx(in2, fr2);
        if (b1 && !b2)
        {
            fprintf(stdout, "\nEnd of file on %s but not on %s\n", fn2, fn1);
        }
        else if (!b1 && b2)
        {
            fprintf(stdout, "\nEnd of file on %s but not on %s\n", fn1, fn2);
        }
        else if (!b1 && !b2)
        {
            fprintf(stdout, "\nFiles read successfully\n");
        }
        else
        {
            cmp_real(stdout, "t", -1, fr1->t, fr2->t, ftol, abstol);
            cmp_int(stdout, "step", -1, fr1->step, fr2->step);
            /* We don't want to print the nre mismatch for every frame */
            /* cmp_int(stdout,"nre",-1,fr1->nre,fr2->nre); */
            if ((fr1->nre >= nre) && (fr2->nre >= nre))
            {
                cmp_energies(stdout, fr1->step, fr1->step, fr1->ener, fr2->ener,
                             enm1, ftol, abstol, nre, ind1, ind2, maxener);
            }
            /*cmp_disres(fr1,fr2,ftol,abstol);*/
            cmp_eblocks(fr1, fr2, ftol, abstol);
        }
    }
    while (b1 && b2);

    close_enx(in1);
    close_enx(in2);

    free_enxframe(fr2);
    sfree(fr2);
    free_enxframe(fr1);
    sfree(fr1);
}
