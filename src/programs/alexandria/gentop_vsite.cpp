/*
 * This source file is part of the Aleandria project.
 *
 * Copyright (C) 2014 David van der Spoel and Paul J. van Maaren
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */
/*! \internal \brief
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "gromacs/legacyheaders/network.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/gmxpreprocess/toputil.h"
#include "poldata.h"
#include "gentop_vsite.h"

typedef struct {
    int nline; /* Must be 3 or 4 */
    int a[4];
} gv_linear;

typedef struct {
    int a[4];
    int nb[4];
} gv_planar;

typedef struct {
    int natom;
    int a[6];
    int nb[6];
} gv_ringplanar;

typedef struct gentop_vsite
{
    int            egvt;
    int            nlinear;
    gv_linear     *lin;
    int            nplanar;
    gv_planar     *plan;
    int            nringplanar;
    gv_ringplanar *rplan;
} gentop_vsite;

gentop_vsite_t gentop_vsite_init(int egvt)
{
    gentop_vsite_t gvt;

    snew(gvt, 1);
    gvt->egvt = egvt;
    if ((egvt < egvtNO) || (egvt >= egvtNR))
    {
        gvt->egvt = egvtNO;
    }

    return gvt;
}

void gentop_vsite_done(gentop_vsite_t *gvt)
{
    sfree((*gvt)->lin);
    sfree((*gvt)->plan);
    sfree(*gvt);
    *gvt = NULL;
}

void gentop_vsite_add_linear(gentop_vsite_t gvt, int ai, int aj, int ak)
{
    int i;

    for (i = 0; (i < gvt->nlinear); i++)
    {
        if (gvt->lin[i].a[1] == aj)
        {
            if (((gvt->lin[i].a[0] == ai)  && (gvt->lin[i].a[2] == ak)) ||
                ((gvt->lin[i].a[0] == aj)  && (gvt->lin[i].a[2] == ai)))
            {
                break;
            }
        }
    }
    if (i == gvt->nlinear)
    {
        srenew(gvt->lin, ++gvt->nlinear);
        gvt->lin[i].nline = NOTSET;
        gvt->lin[i].a[0]  = ai;
        gvt->lin[i].a[1]  = aj;
        gvt->lin[i].a[2]  = ak;
    }
}

void gentop_vsite_add_planar(gentop_vsite_t gvt, int ai, int aj, int ak, int al, int nbonds[])
{
    int i, j;

    for (i = 0; (i < gvt->nplanar); i++)
    {
        if (((gvt->plan[i].a[0] == ai) && (gvt->plan[i].a[1] == aj)  &&
             (gvt->plan[i].a[2] == ak) && (gvt->plan[i].a[3] == al)) ||
            ((gvt->plan[i].a[0] == al) && (gvt->plan[i].a[1] == ak)  &&
             (gvt->plan[i].a[2] == aj) && (gvt->plan[i].a[3] == ai)))
        {
            break;
        }
    }
    if (i == gvt->nplanar)
    {
        srenew(gvt->plan, ++gvt->nplanar);
        gvt->plan[i].a[0] = ai;
        gvt->plan[i].a[1] = aj;
        gvt->plan[i].a[2] = ak;
        gvt->plan[i].a[3] = al;
        for (j = 0; (j < 4); j++)
        {
            gvt->plan[i].nb[j] = nbonds[gvt->plan[i].a[j]];
        }
    }
}

void gentop_vsite_add_ring_planar(gentop_vsite_t gvt,
                                  int natom, int aa[], int nbonds[])
{
    int      i, j;
    gmx_bool bSame;

    for (i = 0; (i < gvt->nringplanar); i++)
    {
        bSame = TRUE;
        for (j = 0; (j < natom); j++)
        {
            bSame = bSame && (gvt->rplan[i].a[j] == aa[j]);
        }
        if (bSame)
        {
            break;
        }
    }
    if (i == gvt->nringplanar)
    {
        srenew(gvt->rplan, ++gvt->nringplanar);
        for (j = 0; (j < natom); j++)
        {
            gvt->rplan[i].a[j]  = aa[j];
            gvt->rplan[i].nb[j] = nbonds[gvt->plan[i].a[j]];
        }
    }
}

static void delete_params(t_params plist[], int etype, int alist[])
{
    int j, k, l, nra;

    nra = interaction_function[etype].nratoms;
    switch (nra)
    {
        case 2:
            /* Remove bonds, if present */
            for (j = 0; (j < plist[etype].nr); j++)
            {
                if (((plist[etype].param[j].a[0] == alist[0]) &&
                     (plist[etype].param[j].a[1] == alist[1])) ||
                    ((plist[etype].param[j].a[1] == alist[0]) &&
                     (plist[etype].param[j].a[0] == alist[1])))
                {
                    if (NULL != debug)
                    {
                        fprintf(debug, "Removing bond beteen atoms %d %d\n",
                                alist[0], alist[1]);
                    }
                    for (k = j+1; (k < plist[etype].nr); k++)
                    {
                        for (l = 0; (l < MAXATOMLIST); l++)
                        {
                            plist[etype].param[k-1].a[l] =
                                plist[etype].param[k].a[l];
                        }
                        for (l = 0; (l < MAXFORCEPARAM); l++)
                        {
                            plist[etype].param[k-1].c[l] =
                                plist[etype].param[k].c[l];
                        }
                    }
                    plist[etype].nr--;
                    j--;
                    break;
                }
            }
            break;
        case 3:
            /* Remove angle, if present */
            for (j = 0; (j < plist[etype].nr); j++)
            {
                if (plist[etype].param[j].a[1] == alist[1])
                {
                    if (((plist[etype].param[j].a[0] == alist[0]) &&
                         (plist[etype].param[j].a[2] == alist[2])) ||
                        ((plist[etype].param[j].a[2] == alist[0]) &&
                         (plist[etype].param[j].a[0] == alist[2])))
                    {
                        if (NULL != debug)
                        {
                            fprintf(debug, "Removing angle beteen atoms %d %d %d\n",
                                    alist[0], alist[1], alist[2]);
                        }
                        for (k = j+1; (k < plist[etype].nr); k++)
                        {
                            for (l = 0; (l < MAXATOMLIST); l++)
                            {
                                plist[etype].param[k-1].a[l] =
                                    plist[etype].param[k].a[l];
                            }
                            for (l = 0; (l < MAXFORCEPARAM); l++)
                            {
                                plist[etype].param[k-1].c[l] =
                                    plist[etype].param[k].c[l];
                            }
                        }
                        plist[etype].nr--;
                        j--;
                        break;
                    }
                }
            }
            break;
        case 4:
            /* Remove dihedral, if present. Allow wildcard in alist[3] (specified as -1) */
            for (j = 0; (j < plist[etype].nr); j++)
            {
                if (((plist[etype].param[j].a[0] == alist[0]) &&
                     (plist[etype].param[j].a[1] == alist[1]) &&
                     (plist[etype].param[j].a[2] == alist[2]) &&
                     ((alist[3] == -1) || (plist[etype].param[j].a[3] == alist[3]))) ||
                    ((plist[etype].param[j].a[3] == alist[0]) &&
                     (plist[etype].param[j].a[2] == alist[1]) &&
                     (plist[etype].param[j].a[1] == alist[2]) &&
                     ((alist[3] == -1) || (plist[etype].param[j].a[0] == alist[3]))) ||
                    ((plist[etype].param[j].a[1] == alist[0]) &&
                     (plist[etype].param[j].a[2] == alist[1]) &&
                     (plist[etype].param[j].a[3] == alist[2]) &&
                     (alist[3] == -1)))
                {
                    if (NULL != debug)
                    {
                        fprintf(debug, "Removing dihedral beteen atoms %d %d %d %d\n",
                                alist[0], alist[1], alist[2], alist[3]);
                    }
                    for (k = j+1; (k < plist[etype].nr); k++)
                    {
                        for (l = 0; (l < MAXATOMLIST); l++)
                        {
                            plist[etype].param[k-1].a[l] =
                                plist[etype].param[k].a[l];
                        }
                        for (l = 0; (l < MAXFORCEPARAM); l++)
                        {
                            plist[etype].param[k-1].c[l] =
                                plist[etype].param[k].c[l];
                        }
                    }
                    plist[etype].nr--;
                    j--;
                }
            }
            break;
        default:
            fprintf(stderr, "Don't know how to remove params from type %s\n",
                    interaction_function[etype].name);
    }
}

static void calc_vsite2parm(t_atoms *atoms, t_params plist[], rvec **x,
                            gv_linear *gvl, t_symtab *symtab,
                            gpp_atomtype_t atype)
{
    int              i, j, natoms, mt;
    const   char    *ml = "ML";
    double           mI, mJ, mK, mL, mT, com, I;
    double           rB, rC, rD, rVV, mV, ac[4];
    rvec             dx;
    t_param          pp, nbml;
    t_atom           aml;

    if (gvl->nline <= 0)
    {
        return;
    }
    memset(&nbml, 0, sizeof(nbml));
    rvec_sub((*x)[gvl->a[0]], (*x)[gvl->a[1]], dx);
    rB    = norm(dx);
    rvec_sub((*x)[gvl->a[1]], (*x)[gvl->a[2]], dx);
    rC    = rB+norm(dx);
    mI    = atoms->atom[gvl->a[0]].m;
    mJ    = atoms->atom[gvl->a[1]].m;
    mK    = atoms->atom[gvl->a[2]].m;
    if (gvl->nline == 4)
    {
        rvec_sub((*x)[gvl->a[2]], (*x)[gvl->a[3]], dx);
        rD = rC+norm(dx);
        mL = atoms->atom[gvl->a[3]].m;
    }
    else
    {
        mL = 0;
        rD = 0;
    }
    mT    = mI+mJ+mK+mL;
    /* We need to keep the COM at the same position and the moment of inertia.
     * In order to do this we have two variables, the position of the dummy
     * and the relative masses (we also need to keep the total mass constant).
     * The atom I should be the one connecting to the remainder of the molecule.
     * We put the first atom I at coordinate 0.
     */
    com   = (mJ*rB+mK*rC+mL*rD)/(mT);
    I     = mI*sqr(com) + mJ*sqr(rB-com) + mK*sqr(rC-com) + mL*sqr(rD-com);
    rVV   = com+I/(com*mT);
    mV    = com*mT/rVV;
    if (NULL != debug)
    {
        fprintf(debug, "com = %g, I = %g, rVV = %g mV = %g rB = %g rC = %g rD = %g\n",
                com, I, rVV, mV, rB, rC, rD);
    }
    mI    = (mJ+mK+mL-mV);
    if (mI <= 0)
    {
        gmx_fatal(FARGS, "Zero or negative mass %f in virtual site construction", mI);
    }
    ac[0] = 0;
    ac[1] = (rB/rVV);
    ac[2] = (rC/rVV);
    ac[3] = (rD/rVV);

    natoms = atoms->nr;
    add_t_atoms(atoms, 1, 0);
    srenew(*x, natoms+1);

    /* Add coordinates for mass-particles */
    clear_rvec((*x)[natoms]);
    for (i = 1; (i < 3); i++)
    {
        for (j = 0; (j < DIM); j++)
        {
            (*x)[natoms][j] += (atoms->atom[gvl->a[i]].m/mV)*(*x)[gvl->a[i]][j];
        }
    }
    /* Update t_atoms for atoms that change/lose their masses */
    atoms->atom[gvl->a[0]].m  = mI;
    atoms->atom[gvl->a[0]].mB = mI;
    for (i = 1; (i < 3); i++)
    {
        atoms->atom[gvl->a[i]].m  = 0;
        atoms->atom[gvl->a[i]].mB = 0;
    }
    /* Set information in atoms structure for the new mass-particles */
    memset(&aml, 0, sizeof(aml));
    mt = add_atomtype(atype, symtab, &aml, ml, &nbml, NOTSET, 0, 0, 0, 0, 0, 0);
    for (i = natoms; (i <= natoms); i++)
    {
        atoms->atom[i].m          = mV;
        atoms->atom[i].mB         = mV;
        atoms->atom[i].atomnumber = 0;
        if (NULL != atoms->atomname)
        {
            atoms->atomname[i]        = put_symtab(symtab, ml);
        }
        if (NULL != atoms->atomtype)
        {
            atoms->atomtype[i]        = put_symtab(symtab, ml);
        }
        if (NULL != atoms->atomtypeB)
        {
            atoms->atomtypeB[i]       = put_symtab(symtab, ml);
        }
        atoms->atom[i].type       = mt;
        atoms->atom[i].typeB      = mt;
    }

    /* Add constraint between masses */
    memset(&pp, 0, sizeof(pp));
    pp.a[0] = gvl->a[0];
    pp.a[1] = natoms;
    pp.c[0] = rVV;
    add_param_to_list(&(plist[F_CONSTR]), &pp);

    /* Add vsites */
    for (i = 1; (i < gvl->nline); i++)
    {
        memset(&pp, 0, sizeof(pp));
        pp.a[0] = gvl->a[i];
        pp.a[1] = gvl->a[0];
        pp.a[2] = natoms;
        pp.c[0] = ac[i];
        add_param_to_list(&(plist[F_VSITE2]), &pp);
    }
}

void gentop_vsite_merge_linear(gentop_vsite_t gvt, gmx_bool bGenVsites)
{
    int i, j, k, l, ai, aj, ndbl;

    for (i = 0; (i < gvt->nlinear); i++)
    {
        gvt->lin[i].nline = 3;
    }

    if (!bGenVsites)
    {
        return;
    }

    for (i = 0; (i < gvt->nlinear); i++)
    {
        for (j = i+1; (j < gvt->nlinear); j++)
        {
            ndbl = 0;
            for (k = 0; (k < gvt->lin[i].nline); k++)
            {
                ai = gvt->lin[i].a[k];
                for (l = 0; (l < gvt->lin[j].nline); l++)
                {
                    aj = gvt->lin[j].a[l];
                    if (ai == aj)
                    {
                        ndbl++;
                    }
                }
            }

            ndbl /= 2;
            if (ndbl > 0)
            {
                fprintf(stderr, "WARNING: merging two linear vsites into one. Please check result.\n");
                if (NULL != debug)
                {
                    fprintf(debug, "Linear group j");
                    for (l = 0; (l < gvt->lin[j].nline); l++)
                    {
                        fprintf(debug, " %d", gvt->lin[j].a[l]);
                    }
                    fprintf(debug, "\n");
                    fprintf(debug, "Linear group i");
                    for (k = 0; (k < gvt->lin[i].nline); k++)
                    {
                        fprintf(debug, " %d", gvt->lin[i].a[k]);
                    }
                    fprintf(debug, "\n");
                }
                if ((gvt->lin[j].a[0] == gvt->lin[i].a[1]) &&
                    (gvt->lin[j].a[1] == gvt->lin[i].a[2]))
                {
                    gvt->lin[i].a[3]  = gvt->lin[j].a[2];
                    gvt->lin[i].nline = 4;
                    gvt->lin[j].nline = 0;
                }
                else if ((gvt->lin[i].a[0] == gvt->lin[j].a[1]) &&
                         (gvt->lin[i].a[1] == gvt->lin[j].a[2]))
                {
                    gvt->lin[j].a[3]  = gvt->lin[i].a[2];
                    gvt->lin[j].nline = 4;
                    gvt->lin[i].nline = 0;
                }
                else
                {
                    gmx_fatal(FARGS, "Atoms in strange order in linear vsites. Check debug file.");
                }
                if (NULL != debug)
                {
                    fprintf(debug, "Linear group j");
                    for (l = 0; (l < gvt->lin[j].nline); l++)
                    {
                        fprintf(debug, " %d", gvt->lin[j].a[l]);
                    }
                    fprintf(debug, "\n");
                    fprintf(debug, "Linear group i");
                    for (k = 0; (k < gvt->lin[i].nline); k++)
                    {
                        fprintf(debug, " %d", gvt->lin[i].a[k]);
                    }
                    fprintf(debug, "\n");
                }
            }
        }
    }
}

void gentop_vsite_generate_special(gentop_vsite_t gvt, gmx_bool bGenVsites,
                                   t_atoms *atoms, rvec **x,
                                   t_params plist[], t_symtab *symtab,
                                   gpp_atomtype_t atype, t_excls **excls,
                                   gmx_poldata_t pd)
{
    int     i, j, nlin_at;
    int     a[MAXATOMLIST], aa[2];
    t_param pp;
    int     ftb, fta, ftp, fti;

    ftb = gmx_poldata_get_bond_ftype(pd);
    fta = gmx_poldata_get_angle_ftype(pd);
    ftp = gmx_poldata_get_dihedral_ftype(pd, egdPDIHS);
    fti = gmx_poldata_get_dihedral_ftype(pd, egdIDIHS);
    gentop_vsite_merge_linear(gvt, bGenVsites);
    nlin_at = 0;
    for (i = 0; (i < gvt->nlinear); i++)
    {
        nlin_at += gvt->lin[i].nline;
    }

    if (NULL != debug)
    {
        fprintf(debug, "Generating %d linear ", nlin_at);
        if (bGenVsites)
        {
            fprintf(debug, "vsites");
        }
        else
        {
            fprintf(debug, "angles");
        }
        fprintf(debug, " and %d impropers\n", gvt->nplanar);
    }
    if ((gvt->egvt == egvtLINEAR) || (gvt->egvt == egvtALL))
    {
        /* If we use vsites (discouraged) each triplet of atoms in a linear arrangement
         * is described by
         * two dummy masses connected by a constraint, that maintain the
         * moment of inertia. The three atoms are then generated from these two
         * positions using a virtual_site2 construction. We need to add 2 extra
         * particles for each linear group.
         * In case we use the special linear angle terms, life gets a lot easier!
         */
        for (i = 0; (i < gvt->nlinear); i++)
        {
            for (j = 0; (j < gvt->lin[i].nline); j++)
            {
                a[j] = gvt->lin[i].a[j];
            }
            for (; (j < MAXATOMLIST); j++)
            {
                a[j] = -1;
            }
            delete_params(plist, fta, a);
            delete_params(plist, ftp, a);
            delete_params(plist, fti, a);

            if (bGenVsites)
            {
                /* Complicated algorithm, watch out */
                for (j = 0; (j < gvt->lin[i].nline-1); j++)
                {
                    aa[0] = a[j+0];
                    aa[1] = a[j+1];
                    delete_params(plist, ftb, aa);
                }

                /* Compute details for the new masses and vsites, and update everything */
                calc_vsite2parm(atoms, plist, x, &gvt->lin[i], symtab, atype);
                srenew((*excls), atoms->nr);
                for (j = atoms->nr-2; (j <= atoms->nr-1); j++)
                {
                    (*excls)[j].nr = 1;
                    snew((*excls)[j].e, 1);
                    (*excls)[j].e[0] = j;
                }
            }
            else
            {
                /* Linear angles */
                memset(&pp, 0, sizeof(pp));
                for (j = 0; (j < 3); j++)
                {
                    pp.a[j] = a[j];
                }
                add_param_to_list(&(plist[F_LINEAR_ANGLES]), &pp);
            }
        }
    }
    if ((gvt->egvt == egvtPLANAR) || (gvt->egvt == egvtALL))
    {
        for (i = 0; (i < gvt->nplanar); i++)
        {
            /* First delete superfluous dihedrals */
            a[1] = gvt->plan[i].a[0];
            a[3] = -1;
            for (j = 1; (j < 4); j++)
            {
                if (gvt->plan[i].nb[j] == 1)
                {
                    if (j == 1)
                    {
                        a[0] = gvt->plan[i].a[1];
                        a[2] = gvt->plan[i].a[2];
                        delete_params(plist, ftp, a);
                        delete_params(plist, fti, a);
                        a[2] = gvt->plan[i].a[3];
                        delete_params(plist, ftp, a);
                        delete_params(plist, fti, a);
                    }
                    else if (j == 2)
                    {
                        a[0] = gvt->plan[i].a[2];
                        a[2] = gvt->plan[i].a[1];
                        delete_params(plist, ftp, a);
                        delete_params(plist, fti, a);
                        a[2] = gvt->plan[i].a[3];
                        delete_params(plist, ftp, a);
                        delete_params(plist, fti, a);
                    }
                    else if (j == 3)
                    {
                        a[0] = gvt->plan[i].a[3];
                        a[2] = gvt->plan[i].a[1];
                        delete_params(plist, ftp, a);
                        delete_params(plist, fti, a);
                        a[2] = gvt->plan[i].a[2];
                        delete_params(plist, ftp, a);
                        delete_params(plist, fti, a);
                    }
                }
            }
            /* Now add impropers! */
            memset(&pp, 0, sizeof(pp));
            for (j = 0; (j < 4); j++)
            {
                pp.a[j] = gvt->plan[i].a[j];
            }
            add_param_to_list(&(plist[fti]), &pp);
        }
    }
}
