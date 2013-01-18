/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "typedefs.h"
#include "xmdrun.h"
#include "futil.h"
#include "xvgr.h"
#include "macros.h"
#include "physics.h"
#include "network.h"
#include "smalloc.h"
#include "string2.h"
#include "readinp.h"
#include "readir.h"
#include "filenm.h"
#include "names.h"
#include "gmxfio.h"

const char *eoNames[eoNR] = {
    "Pres", "Epot", "Vir", "Dist", "Mu", "Force", "Fx", "Fy", "Fz",
    "Px", "Py", "Pz",
    "Polarizability", "Dipole", "Memory", "UseEinter", "UseVirial",
    "CombinationRule"
};

static int Name2eo(char *s)
{
    int i, res;

    res = -1;

    for (i = 0; (i < eoNR); i++)
    {
        if (gmx_strcasecmp(s, eoNames[i]) == 0)
        {
            res = i;
            fprintf(stderr, "Coupling to observable %d (%s)\n", res, eoNames[res]);
            break;
        }
    }

    return res;
}

#define  block_bc(cr,   d) gmx_bcast(     sizeof(d),     &(d), (cr))
#define nblock_bc(cr, nr, d) gmx_bcast((nr)*sizeof((d)[0]), (d), (cr))
#define   snew_bc(cr, d, nr) { if (!MASTER(cr)) {snew((d), (nr)); }}

static void low_comm_tcr(t_commrec *cr, t_coupl_rec *tcr)
{
    nblock_bc(cr, eoObsNR, tcr->ref_value);

    block_bc(cr, tcr->nLJ);
    snew_bc(cr, tcr->tcLJ, tcr->nLJ);
    nblock_bc(cr, tcr->nLJ, tcr->tcLJ);

    block_bc(cr, tcr->nBU);
    snew_bc(cr, tcr->tcBU, tcr->nBU);
    nblock_bc(cr, tcr->nBU, tcr->tcBU);

    block_bc(cr, tcr->nQ);
    snew_bc(cr, tcr->tcQ, tcr->nQ);
    nblock_bc(cr, tcr->nQ, tcr->tcQ);

    block_bc(cr, tcr->nmemory);
    block_bc(cr, tcr->bInter);
    block_bc(cr, tcr->bVirial);
    block_bc(cr, tcr->combrule);
}

void comm_tcr(FILE *log, t_commrec *cr, t_coupl_rec **tcr)
{
    if (!MASTER(cr))
    {
        snew(*tcr, 1);
    }

    low_comm_tcr(cr, *tcr);
}

static void clear_lj(t_coupl_LJ *tc)
{
    tc->at_i   = 0;
    tc->at_j   = 0;
    tc->eObs   = -1;
    tc->bPrint = TRUE;
    tc->c6     = 0.0;
    tc->c12    = 0.0;
    tc->xi_6   = 0.0;
    tc->xi_12  = 0.0;
}

static void clear_bu(t_coupl_BU *tc)
{
    tc->at_i   = 0;
    tc->at_j   = 0;
    tc->eObs   = -1;
    tc->bPrint = TRUE;
    tc->a      = 0.0;
    tc->b      = 0.0;
    tc->c      = 0.0;
    tc->xi_a   = 0.0;
    tc->xi_b   = 0.0;
    tc->xi_c   = 0.0;
}

static void clear_q(t_coupl_Q *tc)
{
    tc->at_i   = 0;
    tc->eObs   = -1;
    tc->bPrint = TRUE;
    tc->Q      = 0.0;
    tc->xi_Q   = 0.0;
}

void copy_ff(t_coupl_rec *tcr, t_forcerec *fr, t_mdatoms *md, t_idef *idef)
{
    int         i, j, ati, atj, type;
    t_coupl_LJ *tclj;
    t_coupl_BU *tcbu;
    t_coupl_Q  *tcq;

    /* Save values for printing */
    for (i = 0; (i < tcr->nLJ); i++)
    {
        tclj = &(tcr->tcLJ[i]);

        ati  = tclj->at_i;
        atj  = tclj->at_j;
        if (atj == -1)
        {
            atj = ati;
        }
        /* nbfp now includes the 6.0/12.0 derivative prefactors */
        tclj->c6  = C6(fr->nbfp, fr->ntype, ati, atj)/6.0;
        tclj->c12 = C12(fr->nbfp, fr->ntype, ati, atj)/12.0;
    }

    for (i = 0; (i < tcr->nBU); i++)
    {
        tcbu = &(tcr->tcBU[i]);

        ati  = tcbu->at_i;
        atj  = tcbu->at_j;
        if (atj == -1)
        {
            atj = ati;
        }
        /* nbfp now includes the 6.0 derivative prefactor */
        tcbu->a = BHAMA(fr->nbfp, fr->ntype, ati, atj);
        tcbu->b = BHAMB(fr->nbfp, fr->ntype, ati, atj);
        tcbu->c = BHAMC(fr->nbfp, fr->ntype, ati, atj)/6.0;
    }

    for (i = 0; (i < tcr->nQ); i++)
    {
        tcq = &(tcr->tcQ[i]);
        for (j = 0; (j < md->nr); j++)
        {
            if (md->typeA[j] == tcq->at_i)
            {
                tcr->tcQ[i].Q = md->chargeA[j];
                break;
            }
        }
    }
    for (i = 0; (i < tcr->nIP); i++)
    {
        /* Let's just copy the whole struct !*/
        type               = tcr->tIP[i].type;
        tcr->tIP[i].iprint = idef->iparams[type];
    }
}

void write_gct(const char *fn, t_coupl_rec *tcr, t_idef *idef)
{
    FILE *fp;
    int   i, ftype;

    fp = gmx_fio_fopen(fn, "w");
    nice_header(fp, fn);
    fprintf(fp, "%-15s = %12g  ; Reference pressure for coupling\n",
            eoNames[eoPres], tcr->ref_value[eoPres]);
    fprintf(fp, "%-15s = %12g  ; Reference potential energy\n",
            eoNames[eoEpot], tcr->ref_value[eoEpot]);
    fprintf(fp, "%-15s = %12g  ; Reference distance\n",
            eoNames[eoDist], tcr->ref_value[eoDist]);
    fprintf(fp, "%-15s = %12g  ; Reference dipole\n",
            eoNames[eoMu], tcr->ref_value[eoMu]);
    fprintf(fp, "%-15s = %12g  ; Reference force\n",
            eoNames[eoForce], tcr->ref_value[eoForce]);
    fprintf(fp, "%-15s = %12g  ; Reference force in X dir\n",
            eoNames[eoFx], tcr->ref_value[eoFx]);
    fprintf(fp, "%-15s = %12g  ; Reference force in Y dir\n",
            eoNames[eoFy], tcr->ref_value[eoFy]);
    fprintf(fp, "%-15s = %12g  ; Reference force in Z dir\n",
            eoNames[eoFz], tcr->ref_value[eoFz]);
    fprintf(fp, "%-15s = %12g  ; Reference pres in X dir\n",
            eoNames[eoPx], tcr->ref_value[eoPx]);
    fprintf(fp, "%-15s = %12g  ; Reference pres in Y dir\n",
            eoNames[eoPy], tcr->ref_value[eoPy]);
    fprintf(fp, "%-15s = %12g  ; Reference pres in Z dir\n",
            eoNames[eoPz], tcr->ref_value[eoPz]);
    fprintf(fp, "%-15s = %12g  ; Polarizability used for the Epot correction\n",
            eoNames[eoPolarizability], tcr->ref_value[eoPolarizability]);
    fprintf(fp, "%-15s = %12g  ; Gas phase dipole moment used for Epot correction\n",
            eoNames[eoDipole], tcr->ref_value[eoDipole]);
    fprintf(fp, "%-15s = %12d  ; Memory for coupling. Makes it converge faster.\n",
            eoNames[eoMemory], tcr->nmemory);
    fprintf(fp, "%-15s = %12s  ; Use intermolecular Epot only (LJ+Coul)\n",
            eoNames[eoInter], yesno_names[tcr->bInter]);
    fprintf(fp, "%-15s = %12s  ; Use virial iso pressure\n",
            eoNames[eoUseVirial], yesno_names[tcr->bVirial]);
    fprintf(fp, "%-15s = %12d  ; Combination rule, same coding as in grompp.\n",
            eoNames[eoCombRule], tcr->combrule);

    fprintf(fp, "\n; Q-Coupling   %6s  %12s\n", "type", "xi");
    for (i = 0; (i < tcr->nQ); i++)
    {
        fprintf(fp, "%-8s = %8s  %6d  %12g\n",
                "Q", eoNames[tcr->tcQ[i].eObs], tcr->tcQ[i].at_i, tcr->tcQ[i].xi_Q);
    }

    fprintf(fp, "\n; %8s %8s  %6s  %6s  %12s  %12s\n", "Couple", "To",
            "i-type", "j-type", "xi-c6", "xi-c12");
    fprintf(fp, "; j-type == -1 means mixing rules will be applied!\n");
    for (i = 0; (i < tcr->nLJ); i++)
    {
        fprintf(fp, "%-8s = %8s  %6d  %6d  %12g  %12g\n",
                "LJ", eoNames[tcr->tcLJ[i].eObs],
                tcr->tcLJ[i].at_i, tcr->tcLJ[i].at_j,
                tcr->tcLJ[i].xi_6, tcr->tcLJ[i].xi_12);
    }

    fprintf(fp, "\n; %8s %8s  %6s  %6s  %12s  %12s  %12s\n", "Couple", "To",
            "i-type", "j-type", "xi-A", "xi-B", "xi-C");
    fprintf(fp, "; j-type == -1 means mixing rules will be applied!\n");
    for (i = 0; (i < tcr->nBU); i++)
    {
        fprintf(fp, "%-8s = %8s  %6d  %6d  %12g  %12g  %12g\n",
                "BU", eoNames[tcr->tcBU[i].eObs],
                tcr->tcBU[i].at_i, tcr->tcBU[i].at_j,
                tcr->tcBU[i].xi_a, tcr->tcBU[i].xi_b, tcr->tcBU[i].xi_c);
    }

    fprintf(fp, "\n; More Coupling\n");
    for (i = 0; (i < tcr->nIP); i++)
    {
        ftype = idef->functype[tcr->tIP[i].type];
        switch (ftype)
        {
            case F_BONDS:
                fprintf(fp, "%-15s = %-8s  %4d  %12g  %12g\n",
                        "Bonds", eoNames[tcr->tIP[i].eObs], tcr->tIP[i].type,
                        tcr->tIP[i].xi.harmonic.krA,
                        tcr->tIP[i].xi.harmonic.rA);
                break;
            default:
                fprintf(stderr, "ftype %s not supported (yet)\n",
                        interaction_function[ftype].longname);
        }
    }
    gmx_fio_fclose(fp);
}

static gmx_bool add_lj(int *nLJ, t_coupl_LJ **tcLJ, char *s, gmx_bool bObsUsed[])
{
    int       j, ati, atj, eo;
    char      buf[256];
    double    xi6, xi12;

    if (sscanf(s, "%s%d%d%lf%lf", buf, &ati, &atj, &xi6, &xi12) != 5)
    {
        return TRUE;
    }
    if ((eo = Name2eo(buf)) == -1)
    {
        gmx_fatal(FARGS, "Invalid observable for LJ coupling: %s", buf);
    }

    for (j = 0; (j < *nLJ); j++)
    {
        if ((((*tcLJ)[j].at_i == ati) && ((*tcLJ)[j].at_j == atj)) &&
            ((*tcLJ)[j].xi_6 || (*tcLJ)[j].xi_12) &&
            ((*tcLJ)[j].eObs == eo))
        {
            break;
        }
    }
    if (j == *nLJ)
    {
        ++(*nLJ);
        srenew((*tcLJ), *nLJ);
    }
    else
    {
        fprintf(stderr, "\n*** WARNING: overwriting entry for LJ coupling '%s'\n", s);
    }

    clear_lj(&((*tcLJ)[j]));
    if (((*tcLJ)[j].eObs = eo) == -1)
    {
        gmx_fatal(FARGS, "Invalid observable for LJ coupling: %s", buf);
    }
    (*tcLJ)[j].at_i   = ati;
    (*tcLJ)[j].at_j   = atj;
    (*tcLJ)[j].xi_6   = xi6;
    (*tcLJ)[j].xi_12  = xi12;
    bObsUsed[eo]      = TRUE;

    return FALSE;
}

static gmx_bool add_bu(int *nBU, t_coupl_BU **tcBU, char *s, gmx_bool bObsUsed[])
{
    int       j, ati, atj, eo;
    char      buf[256];
    double    xia, xib, xic;

    if (sscanf(s, "%s%d%d%lf%lf%lf", buf, &ati, &atj, &xia, &xib, &xic) != 6)
    {
        return TRUE;
    }
    if ((eo = Name2eo(buf)) == -1)
    {
        gmx_fatal(FARGS, "Invalid observable for BU coupling: %s", buf);
    }

    for (j = 0; (j < *nBU); j++)
    {
        if ((((*tcBU)[j].at_i == ati) && ((*tcBU)[j].at_j == atj)) &&
            ((*tcBU)[j].xi_a || (*tcBU)[j].xi_b || (*tcBU)[j].xi_c ) &&
            ((*tcBU)[j].eObs == eo))
        {
            break;
        }
    }
    if (j == *nBU)
    {
        ++(*nBU);
        srenew((*tcBU), *nBU);
    }
    else
    {
        fprintf(stderr, "\n*** WARNING: overwriting entry for BU coupling '%s'\n", s);
    }

    clear_bu(&((*tcBU)[j]));
    if (((*tcBU)[j].eObs = eo) == -1)
    {
        gmx_fatal(FARGS, "Invalid observable for BU coupling: %s", buf);
    }
    (*tcBU)[j].at_i   = ati;
    (*tcBU)[j].at_j   = atj;
    (*tcBU)[j].xi_a   = xia;
    (*tcBU)[j].xi_b   = xib;
    (*tcBU)[j].xi_c   = xic;
    bObsUsed[eo]      = TRUE;

    return FALSE;
}

static gmx_bool add_ip(int *nIP, t_coupl_iparams **tIP, char *s, int ftype, gmx_bool bObsUsed[])
{
    int    i, eo, type;
    char   buf[256];
    double kb, b0;

    switch (ftype)
    {
        case F_BONDS:
            /* Pick out the type */
            if (sscanf(s, "%s%d", buf, &type) != 2)
            {
                return TRUE;
            }
            if ((eo = Name2eo(buf)) == -1)
            {
                gmx_fatal(FARGS, "Invalid observable for IP coupling: %s", buf);
            }

            /* Check whether this entry is there already */
            for (i = 0; (i < *nIP); i++)
            {
                if ((*tIP)[i].type == type)
                {
                    break;
                }
            }
            if (i < *nIP)
            {
                fprintf(stderr, "*** WARNING: overwriting entry for type %d\n", type);
            }
            else
            {
                i = *nIP;
                srenew((*tIP), i+1);
                (*nIP)++;
            }
            if (sscanf(s, "%s%d%lf%lf", buf, &type, &kb, &b0) != 4)
            {
                return TRUE;
            }
            (*tIP)[i].type            = type;
            (*tIP)[i].eObs            = eo;
            (*tIP)[i].xi.harmonic.krA = kb;
            (*tIP)[i].xi.harmonic.rA  = b0;
            bObsUsed[eo]              = TRUE;
            break;
        default:
            fprintf(stderr, "ftype %s not supported (yet)\n",
                    interaction_function[ftype].longname);
            return TRUE;
    }
    return FALSE;
}

static gmx_bool add_q(int *nQ, t_coupl_Q **tcQ, char *s, gmx_bool bObsUsed[])
{
    int       j, ati, eo;
    char      buf[256];
    double    xiQ;

    if (sscanf(s, "%s%d%lf", buf, &ati, &xiQ) != 3)
    {
        return TRUE;
    }

    for (j = 0; (j < *nQ); j++)
    {
        if ((*tcQ)[j].at_i == ati)
        {
            break;
        }
    }
    if (j == *nQ)
    {
        ++(*nQ);
        srenew((*tcQ), *nQ);
    }
    else
    {
        fprintf(stderr, "\n*** WARNING: overwriting entry for Q coupling '%s'\n", s);
    }

    clear_q(&((*tcQ)[j]));
    eo = (*tcQ)[j].eObs = Name2eo(buf);
    if ((*tcQ)[j].eObs == -1)
    {
        gmx_fatal(FARGS, "Invalid observable for Q coupling: %s", buf);
    }
    (*tcQ)[j].at_i   = ati;
    (*tcQ)[j].xi_Q   = xiQ;
    bObsUsed[eo]     = TRUE;

    return FALSE;
}

void read_gct(const char *fn, t_coupl_rec *tcr)
{
    warninp_t     wi;
    t_inpfile    *inp;
    int           i, j, ninp, nQ, nLJ, nBU, nIP;
    gmx_bool      bWrong;

    wi = init_warning(FALSE, 0);

    inp = read_inpfile(fn, &ninp, NULL, wi);

    for (i = 0; (i < eoObsNR); i++)
    {
        tcr->bObsUsed[i] = FALSE;
        RTYPE (eoNames[i],  tcr->ref_value[i],  0.0);
    }
    ITYPE (eoNames[eoMemory],     tcr->nmemory,   1);
    ETYPE (eoNames[eoInter],      tcr->bInter,    yesno_names);
    ETYPE (eoNames[eoUseVirial],  tcr->bVirial,   yesno_names);
    ITYPE (eoNames[eoCombRule],   tcr->combrule,  1);
    tcr->tcLJ = NULL;
    tcr->tcBU = NULL;
    tcr->tcQ  = NULL;
    tcr->tIP  = NULL;
    nQ        = nLJ = nBU = nIP = 0;

    for (i = 0; (i < ninp); i++)
    {
        bWrong = FALSE;
        if (gmx_strcasecmp(inp[i].name, "LJ") == 0)
        {
            bWrong = add_lj(&nLJ, &(tcr->tcLJ), inp[i].value, tcr->bObsUsed);
        }
        else if (gmx_strcasecmp(inp[i].name, "BU") == 0)
        {
            bWrong = add_bu(&nBU, &(tcr->tcBU), inp[i].value, tcr->bObsUsed);
        }
        else if (gmx_strcasecmp(inp[i].name, "Q") == 0)
        {
            bWrong = add_q(&nQ, &(tcr->tcQ), inp[i].value, tcr->bObsUsed);
        }
        else if (gmx_strcasecmp(inp[i].name, "Bonds") == 0)
        {
            bWrong = add_ip(&nIP, &(tcr->tIP), inp[i].value, F_BONDS, tcr->bObsUsed);
        }

        if (bWrong)
        {
            fprintf(stderr, "Wrong line in %s: '%s = %s'\n",
                    fn, inp[i].name, inp[i].value);
        }
        /*sfree(inp[i].name);
           sfree(inp[i].value);*/
    }
    /* Check which ones have to be printed */
    for (i = 1; (i < nQ); i++)
    {
        for (j = 0; (j < i); j++)
        {
            if (tcr->tcQ[i].at_i == tcr->tcQ[j].at_i)
            {
                tcr->tcQ[j].bPrint = FALSE;
            }
        }
    }
    for (i = 1; (i < nLJ); i++)
    {
        for (j = 0; (j < i); j++)
        {
            if (((tcr->tcLJ[i].at_i == tcr->tcLJ[j].at_i) &&
                 (tcr->tcLJ[i].at_j == tcr->tcLJ[j].at_j)) ||
                ((tcr->tcLJ[i].at_i == tcr->tcLJ[j].at_j) &&
                 (tcr->tcLJ[i].at_j == tcr->tcLJ[j].at_i)))
            {
                tcr->tcLJ[j].bPrint = FALSE;
            }
        }
    }

    for (i = 1; (i < nBU); i++)
    {
        for (j = 0; (j < i); j++)
        {
            if (((tcr->tcBU[i].at_i == tcr->tcBU[j].at_i) &&
                 (tcr->tcBU[i].at_j == tcr->tcBU[j].at_j)) ||
                ((tcr->tcBU[i].at_i == tcr->tcBU[j].at_j) &&
                 (tcr->tcBU[i].at_j == tcr->tcBU[j].at_i)))
            {
                tcr->tcBU[j].bPrint = FALSE;
            }
        }
    }

    tcr->nQ  = nQ;
    tcr->nLJ = nLJ;
    tcr->nBU = nBU;
    tcr->nIP = nIP;

    sfree(inp);

    done_warning(wi, FARGS);
}
