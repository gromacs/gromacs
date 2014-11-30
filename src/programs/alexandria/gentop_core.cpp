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
/*! \internal \brief
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include "gmxpre.h"
#include <ctype.h>
#include <stdlib.h>
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/legacyheaders/copyrite.h"
#include "gromacs/listed-forces/bonded.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/fileio/strdb.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/pdbio.h"
#include "gromacs/utility/futil.h"
#include "gromacs/gmxpreprocess/grompp.h"
#include "gromacs/gmxpreprocess/toputil.h"
#include "gromacs/gmxpreprocess/gen_ad.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/random/random.h"
#include "gromacs/legacyheaders/txtdump.h"
#include "gromacs/legacyheaders/readinp.h"
#include "gromacs/legacyheaders/names.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/math/vec.h"
#include "gromacs/topology/atomprop.h"
#include "gromacs/gmxpreprocess/gpp_atomtype.h"
#include "gromacs/pbcutil/pbc.h"
#include "gentop_core.h"
#include "gentop_vsite.h"
#include "poldata.h"
#include "plistwrapper.h"
#include "stringutil.h"

void calc_angles_dihs(t_params *ang, t_params *dih, rvec x[], gmx_bool bPBC,
                      matrix box)
{
    int           i, ai, aj, ak, al, t1, t2, t3;
    rvec          r_ij, r_kj, r_kl, m, n;
    real          sign, th, costh, ph;
    struct t_pbc  pbc;

    if (bPBC)
    {
        set_pbc(&pbc, -1, box);
    }
    if (debug)
    {
        pr_rvecs(debug, 0, "GENTOP", box, DIM);
    }
    for (i = 0; (i < ang->nr); i++)
    {
        ai = ang->param[i].a[0];
        aj = ang->param[i].a[1];
        ak = ang->param[i].a[2];
        th = RAD2DEG*bond_angle(x[ai], x[aj], x[ak], bPBC ? &pbc : NULL,
                                r_ij, r_kj, &costh, &t1, &t2);
        if (debug)
        {
            fprintf(debug, "GENTOP: ai=%3d aj=%3d ak=%3d r_ij=%8.3f r_kj=%8.3f th=%8.3f\n",
                    ai, aj, ak, norm(r_ij), norm(r_kj), th);
        }
        ang->param[i].c[0] = th;
    }
    for (i = 0; (i < dih->nr); i++)
    {
        ai = dih->param[i].a[0];
        aj = dih->param[i].a[1];
        ak = dih->param[i].a[2];
        al = dih->param[i].a[3];
        ph = RAD2DEG*dih_angle(x[ai], x[aj], x[ak], x[al], bPBC ? &pbc : NULL,
                               r_ij, r_kj, r_kl, m, n, &sign, &t1, &t2, &t3);
        if (debug)
        {
            fprintf(debug, "GENTOP: ai=%3d aj=%3d ak=%3d al=%3d r_ij=%8.3f r_kj=%8.3f r_kl=%8.3f ph=%8.3f\n",
                    ai, aj, ak, al, norm(r_ij), norm(r_kj), norm(r_kl), ph);
        }
        dih->param[i].c[0] = ph;
    }
}

void dump_hybridization(FILE *fp, t_atoms *atoms, int nbonds[])
{
    int i;

    for (i = 0; (i < atoms->nr); i++)
    {
        fprintf(fp, "Atom %5s has %d bonds\n", *atoms->atomname[i], nbonds[i]);
    }
}

static void print_pl(FILE *fp, t_params plist[], int ftp, const char *name,
                     char ***atomname)
{
    int i, j, nral, nrfp;

    if (plist[ftp].nr > 0)
    {
        fprintf(fp, "\n");
        fprintf(fp, "[ %s ]\n", name);
        nral = interaction_function[ftp].nratoms;
        nrfp = interaction_function[ftp].nrfpA;
        for (i = 0; (i < plist[ftp].nr); i++)
        {
            for (j = 0; (j < nral); j++)
            {
                fprintf(fp, "  %5s", *atomname[plist[ftp].param[i].a[j]]);
            }
            for (j = 0; (j < nrfp); j++)
            {
                fprintf(fp, "  %10.3e", plist[ftp].param[i].c[j]);
            }
            fprintf(fp, "\n");
        }
    }
}

void print_rtp(const char *filenm, const char *title, t_atoms *atoms,
               t_params plist[], int cgnr[], int nbts, int bts[])
{
    FILE *fp;
    int   i;

    fp = gmx_ffopen(filenm, "w");
    fprintf(fp, "; %s\n", title);
    fprintf(fp, "\n");
    fprintf(fp, "[ %s ]\n", *atoms->resinfo[0].name);
    fprintf(fp, "\n");
    fprintf(fp, "[ atoms ]\n");
    for (i = 0; (i < atoms->nr); i++)
    {
        fprintf(fp, "%-8s  %12s  %8.4f  %5d\n",
                *atoms->atomname[i], *atoms->atomtype[i],
                atoms->atom[i].q, cgnr[i]);
    }
    for (i = 0; (i < nbts); i++)
    {
        print_pl(fp, plist, bts[i], interaction_function[bts[i]].name, atoms->atomname);
    }
    fclose(fp);
}

static int pcompar(const void *a, const void *b)
{
    t_param *pa, *pb;
    int      d;
    pa = (t_param *)a;
    pb = (t_param *)b;

    d = pa->a[0] - pb->a[0];
    if (d == 0)
    {
        d = pa->a[1] - pb->a[1];
    }
    if (d == 0)
    {
        d = pa->a[2] - pb->a[2];
    }
    if (d == 0)
    {
        d = pa->a[3] - pb->a[3];
    }
    /*if (d == 0)
       return strlen(pb->s) - strlen(pa->s);
       else*/
    return d;
}

static int acomp(const void *a, const void *b)
{
    atom_id *aa = (atom_id *)a;
    atom_id *ab = (atom_id *)b;

    return (*aa - *ab);
}

static void my_clean_excls(int nr, t_excls excls[])
{
    int i, j, k;

    for (i = 0; (i < nr); i++)
    {
        if (excls[i].nr > 0)
        {
            qsort(excls[i].e, excls[i].nr, sizeof(excls[i].e[0]), acomp);
            k = 0;
            for (j = 0; (j < excls[i].nr); j++)
            {
                if (excls[i].e[j] != excls[i].e[k])
                {
                    excls[i].e[++k] = excls[i].e[j];
                }
            }
            excls[i].nr = ++k;
        }
    }
}

static void clean_thole(t_params *ps)
{
    int     i, j;
    atom_id a;

    if (ps->nr > 0)
    {
        /* swap atomnumbers in bond if first larger than second: */
        for (i = 0; (i < ps->nr); i++)
        {
            if (ps->param[i].a[2] < ps->param[i].a[0])
            {
                a                 = ps->param[i].a[0];
                ps->param[i].a[0] = ps->param[i].a[2];
                ps->param[i].a[2] = a;
                a                 = ps->param[i].a[1];
                ps->param[i].a[1] = ps->param[i].a[3];
                ps->param[i].a[3] = a;
            }
        }

        /* Sort bonds */
        qsort(ps->param, ps->nr, (size_t)sizeof(ps->param[0]), pcompar);

        /* remove doubles, keep the first one always. */
        j = 1;
        for (i = 1; (i < ps->nr); i++)
        {
            if ((ps->param[i].a[0] != ps->param[j-1].a[0]) ||
                (ps->param[i].a[1] != ps->param[j-1].a[1]) ||
                (ps->param[i].a[2] != ps->param[j-1].a[2]) ||
                (ps->param[i].a[3] != ps->param[j-1].a[3]) )
            {
                cp_param(&(ps->param[j]), &(ps->param[i]));
                j++;
            }
        }
        fprintf(stderr, "Number of Tholes was %d, now %d\n", ps->nr, j);
        ps->nr = j;
    }
    else
    {
        fprintf(stderr, "No Tholes\n");
    }
}

real calc_dip(t_atoms *atoms, rvec x[])
{
    int  i;
    rvec mu, mm;
    real qq;

    clear_rvec(mu);
    for (i = 0; (i < atoms->nr); i++)
    {
        qq = atoms->atom[i].q;
        svmul(qq, x[i], mm);
        rvec_inc(mu, mm);
    }
    return norm(mu)*ENM2DEBYE;
}

void reset_q(t_atoms *atoms)
{
    int i;

    /* Use values from file */
    for (i = 0; (i < atoms->nr); i++)
    {
        atoms->atom[i].qB = atoms->atom[i].q;
    }
}

void symmetrize_charges(gmx_bool bQsym, t_atoms *atoms,
                        std::vector<alexandria::PlistWrapper>::iterator bonds,
                        gmx_poldata_t pd,
                        gmx_atomprop_t aps, const char *symm_string,
                        std::vector<int> &sym_charges)
{
    char  *central, *attached;
    int    nattached, nh, ai, aj, anri, anrj;
    int    anr_central, anr_attached, nrq;
    int    hs[8];
    double qaver, qsum;
    int    hsmin;

    printf("Fix me: symmetrize charges algorithm is broken in: file %s, line %d\n",
           __FILE__, __LINE__ );
    return;
    for (int i = 0; (i < atoms->nr); i++)
    {
        sym_charges.push_back(i);
    }
    if (bQsym)
    {
        if ((NULL != symm_string) && (strlen(symm_string) > 0))
        {
            std::vector<std::string> ss = split(symm_string, ' ');
            if ((int)ss.size() != atoms->nr)
            {
                gmx_fatal(FARGS, "Wrong number (%d) of atom-numbers in symm_string: expected %d",
                          ss.size(), atoms->nr);
            }
            int ii = 0;
            for (std::vector<std::string>::iterator is = ss.begin();
                 (is < ss.end()); ++is)
            {
                sym_charges[ii] = atoi(is->c_str());
                ii++;
            }
        }
        else
        {
            while (gmx_poldata_get_symcharges(pd, &central,
                                              &attached, &nattached) == 1)
            {
                anr_central  = gmx_atomprop_atomnumber(aps, central);
                anr_attached = gmx_atomprop_atomnumber(aps, attached);
                hsmin        = -1;
                for (int i = 0; (i < atoms->nr); i++)
                {
                    if (atoms->atom[i].atomnumber == anr_central)
                    {
                        nh = 0;
                        for (alexandria::ParamIterator j = bonds->beginParam();
                             (j < bonds->endParam()); ++j)
                        {
                            ai   = j->a[0];
                            aj   = j->a[1];
                            anri = atoms->atom[ai].atomnumber;
                            anrj = atoms->atom[aj].atomnumber;

                            if ((ai == i) && (anrj == anr_attached))
                            {
                                hs[nh++] = aj;
                            }
                            else if ((aj == i) && (anri == anr_attached))
                            {
                                hs[nh++] = ai;
                            }
                            if ((hsmin == -1) || (hs[nh-1] < hsmin))
                            {
                                hsmin = hs[nh-1];
                            }
                        }
                        if ((nh == nattached) && (hsmin != -1))
                        {
                            for (int j = 0; (j < nattached); j++)
                            {
                                sym_charges[hs[j]] = hsmin;
                            }
                        }
                    }
                }
            }
        }

        for (int i = 0; (i < atoms->nr); i++)
        {
            qsum = 0;
            nrq  = 0;
            for (int j = i+1; (j < atoms->nr); j++)
            {
                if (sym_charges[j] == sym_charges[i])
                {
                    qsum += atoms->atom[j].q;
                    nrq++;
                }
            }
            if (0 < nrq)
            {
                qaver = qsum/nrq;
                for (int j = 0; (j < atoms->nr); j++)
                {
                    if (sym_charges[j] == sym_charges[i])
                    {
                        atoms->atom[j].q = atoms->atom[j].qB = qaver;
                    }
                }
            }
        }
    }
}

static int *generate_cg_neutral(t_atoms *atoms, gmx_bool bUsePDBcharge)
{
    int     i, n = 1;
    int    *cgnr;
    double  qt = 0;

    snew(cgnr, atoms->nr);
    for (i = 0; (i < atoms->nr); i++)
    {
        if (atoms->pdbinfo && bUsePDBcharge)
        {
            atoms->atom[i].q = atoms->pdbinfo[i].bfac;
        }
        qt     += atoms->atom[i].q;
        cgnr[i] = n;
        if (is_int(qt))
        {
            n++;
            qt = 0;
        }
    }
    return cgnr;
}

static int *generate_cg_group(t_atoms *atoms, 
                              std::vector<alexandria::PlistWrapper> &plist)
{
    int        i, j, k, atn, ai, aj, ncg = 1;
    int       *cgnr;
    gmx_bool   bMV;
    int        monovalent[] = { 0, 1, 9, 17, 35, 53, 85 };
    int        nmv          = asize(monovalent);

    /* Assume that shells and masses have atomnumber 0 */
    snew(cgnr, atoms->nr);
    for (i = 0; (i < atoms->nr); i++)
    {
        cgnr[i] = NOTSET;
    }

    for (i = 0; (i < atoms->nr); i++)
    {
        atn = atoms->atom[i].atomnumber;
        bMV = FALSE;
        for (j = 0; (j < nmv) && !bMV; j++)
        {
            bMV = (atn == monovalent[j]);
        }
        if (!bMV)
        {
            cgnr[i] = ncg++;
        }
    }
    /* Rely on the notion that all H and other monovalent
       atoms are bound to something */
    alexandria::PlistWrapperIterator bonds = SearchPlist(plist, F_BONDS);
    if (plist.end() != bonds)
    {
        for (alexandria::ParamIterator j = bonds->beginParam();
             (j < bonds->endParam()); ++j)
        {
            ai  = j->a[0];
            aj  = j->a[1];
            bMV = FALSE;
            atn = atoms->atom[ai].atomnumber;
            for (k = 0; (k < nmv) && !bMV; k++)
            {
                bMV = (atn == monovalent[k]);
            }
            if (bMV)
            {
                if (cgnr[aj] != NOTSET)
                {
                    cgnr[ai] = cgnr[aj];
                }
                else
                {
                    cgnr[ai] = cgnr[aj] = 1;
                }
            }
            else
            {
                bMV = FALSE;
                atn = atoms->atom[aj].atomnumber;
                for (k = 0; (k < nmv) && !bMV; k++)
                {
                    bMV = (atn == monovalent[k]);
                }
                if (bMV)
                {
                    cgnr[aj] = cgnr[ai];
                }
            }
        }
    }
    /* Rely on the notion that all shells are bound to something */
    alexandria::PlistWrapperIterator pols = SearchPlist(plist, F_POLARIZATION);
    if (plist.end() != pols)
    {
        for (alexandria::ParamIterator j = pols->beginParam();
             (j < pols->endParam()); ++j)
        {
            ai       = j->a[0];
            aj       = j->a[1];
            cgnr[aj] = cgnr[ai];
        }
        for (i = 0; (i < atoms->nr); i++)
        {
            if (cgnr[i] == NOTSET)
            {
                cgnr[i] = ncg++;
            }
        }
    }
    printf("There are %d charge groups\n", ncg-1);

    return cgnr;
}

static int *generate_cg_atom(int natom)
{
    int i, *cgnr;

    snew(cgnr, natom);
    for (i = 0; (i < natom); i++)
    {
        cgnr[i] = i+1;
    }

    return cgnr;
}

int *generate_charge_groups(eChargeGroup cgtp, t_atoms *atoms,
                            std::vector<alexandria::PlistWrapper> &plist,
                            bool bUsePDBcharge,
                            real *qtot, real *mtot)
{
    int i, *cgnr = NULL;
    std::vector<alexandria::PlistWrapper>::iterator pb, ps;
    pb = alexandria::SearchPlist(plist, F_BONDS);
    ps = alexandria::SearchPlist(plist, F_POLARIZATION);

    switch (cgtp)
    {
        case ecgNeutral:
            cgnr = generate_cg_neutral(atoms, bUsePDBcharge);
            break;
        case ecgGroup:
            cgnr = generate_cg_group(atoms, plist);
            break;
        case ecgAtom:
            cgnr = generate_cg_atom(atoms->nr);
            break;
        default:
            gmx_fatal(FARGS, "Invalid charge group generation type %d", cgtp);
    }
    *qtot = *mtot = 0;
    for (i = 0; (i < atoms->nr); i++)
    {
        *qtot += atoms->atom[i].q;
        *mtot += atoms->atom[i].m;
    }
    return cgnr;
}

static int    *cgnr_copy;
static double *atomnumber;
static int cg_comp(const void *a, const void *b)
{
    int   *aa = (int *)a;
    int   *bb = (int *)b;
    double c;

    int    d = cgnr_copy[*aa] - cgnr_copy[*bb];
    if (d == 0)
    {
        c = atomnumber[*aa] - atomnumber[*bb];
        if (c < 0)
        {
            return -1;
        }
        else if (c > 0)
        {
            return 1;
        }
        else
        {
            return 0;
        }
    }
    else
    {
        return d;
    }
}

void sort_on_charge_groups(int *cgnr, t_atoms *atoms, 
                           std::vector<alexandria::PlistWrapper> &pw,
                           rvec x[], t_excls excls[],
                           const char *ndxout, int nmol)
{
    FILE      *fp;
    int        i, j, j0, k, newi, ri, *cg_renum, *ccgg, *inv_renum;
    rvec      *rx;
    t_atom    *ra;
    t_excls   *newexcls;
    char    ***an, ***smn;

    snew(cg_renum, atoms->nr);
    snew(atomnumber, atoms->nr);
    snew(rx, atoms->nr);
    snew(ra, atoms->nr);
    snew(an, atoms->nr);
    for (i = 0; (i < atoms->nr); i++)
    {
        cg_renum[i]   = i;
        atomnumber[i] = 1+i; /*atoms->atom[i].atomnumber;*/
        if ((atoms->atom[i].ptype == eptShell) && (i > 0))
        {
            atomnumber[i] = atomnumber[i-1]+0.1;
        }
    }
    cgnr_copy = cgnr;
    qsort(cg_renum, atoms->nr, sizeof(cg_renum[0]), cg_comp);
    if (debug)
    {
        for (i = 0; (i < atoms->nr); i++)
        {
            fprintf(debug, "cg_renum[%d] = %d\n", i, cg_renum[i]);
        }
    }
    snew(ccgg, atoms->nr);
    for (i = 0; (i < atoms->nr); i++)
    {
        ri = cg_renum[i];
        copy_rvec(x[ri], rx[i]);
        memcpy(&(ra[i]), &(atoms->atom[ri]), sizeof(t_atom));
        an[i]   = atoms->atomname[ri];
        ccgg[i] = cgnr[ri];
    }
    snew(inv_renum, atoms->nr);
    snew(smn, atoms->nr);
    for (i = 0; (i < atoms->nr); i++)
    {
        copy_rvec(rx[i], x[i]);
        memcpy(&(atoms->atom[i]), &(ra[i]), sizeof(t_atom));
        atoms->atomname[i]     = an[i];
        cgnr[i]                = ccgg[i];
        inv_renum[cg_renum[i]] = i;
        smn[i]                 = atoms->atomtype[i];
    }
    for (i = 0; (i < atoms->nr); i++)
    {
        newi               = cg_renum[i];
        atoms->atomtype[i] = smn[newi];
    }
    for (alexandria::PlistWrapperIterator i = pw.begin();
         (i < pw.end()); ++i)
    {
        for (alexandria::ParamIterator j = i->beginParam();
             (j < i->endParam()); j++)
        {
            for (k = 0; (k < NRAL(i->getFtype())); k++)
            {
                j->a[k] = inv_renum[j->a[k]];
            }
        }
    }
    snew(newexcls, atoms->nr);
    for (i = 0; (i < atoms->nr); i++)
    {
        snew(newexcls[i].e, excls[i].nr);
        newexcls[i].nr = excls[i].nr;
        for (j = 0; (j < excls[i].nr); j++)
        {
            newexcls[i].e[j] = excls[i].e[j];
        }
    }
    for (i = 0; (i < atoms->nr); i++)
    {
        newi = inv_renum[i];
        if (newexcls[i].nr > excls[newi].nr)
        {
            srenew(excls[newi].e, newexcls[i].nr);
        }
        for (j = 0; (j < newexcls[i].nr); j++)
        {
            excls[newi].e[j] = inv_renum[newexcls[i].e[j]];
        }
        excls[newi].nr = newexcls[i].nr;
    }
    if (NULL != ndxout)
    {
        fp = fopen(ndxout, "w");
        fprintf(fp, "[ number_backward ]\n");
        for (j = 0; (j < nmol); j++)
        {
            j0 = j*atoms->nr;
            for (i = 0; (i < atoms->nr); i++)
            {
                if (atoms->atom[inv_renum[i]].ptype == eptShell)
                {
                    k = j0+inv_renum[i-1]+1;
                }
                else
                {
                    k = j0+inv_renum[i]+1;
                }
                fprintf(fp, " %d", k);
                if (j == 0)
                {
                    cg_renum[inv_renum[i]] = i;
                }
            }
            fprintf(fp, "\n");
        }
        for (j = 0; (j < nmol); j++)
        {
            j0 = j*atoms->nr;
            fprintf(fp, "[ number_forward ]\n");
            for (i = 0; (i < atoms->nr); i++)
            {
                if (atoms->atom[cg_renum[i]].ptype == eptShell)
                {
                    k = j0+cg_renum[i-1]+1;
                }
                else
                {
                    k = j0+cg_renum[i]+1;
                }
                fprintf(fp, " %d", k);
            }
            fprintf(fp, "\n");
        }
        fclose(fp);
    }
    sfree(rx);
    sfree(ra);
    sfree(an);
    sfree(cg_renum);
    sfree(inv_renum);
    sfree(ccgg);
}
