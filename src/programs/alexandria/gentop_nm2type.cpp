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
#include "gmxpre.h"
#include <assert.h>
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/legacyheaders/copyrite.h"
#include "gromacs/bonded/bonded.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/math/units.h"
#include "gromacs/legacyheaders/names.h"
#include "gromacs/math/vec.h"
#include "gromacs/random/random.h"
#include "gromacs/legacyheaders/txtdump.h"
#include "gromacs/legacyheaders/readinp.h"
#include "gromacs/legacyheaders/names.h"
#include "gromacs/topology/atomprop.h"
#include "poldata.h"
#include "gentop_nm2type.h"
#include "gentop_vsite.h"

enum {
    egmLinear, egmPlanar, egmTetrahedral, egmNR
};
static const char *geoms[egmNR] = { "linear", "planar", "tetrahedral" };

enum {
    ecolWhite, ecolGrey, ecolBlack
};

typedef struct {
    int        nb, geom, conf, atomnr, nat, myat, nn, ecol;
    int       *bbb, *bondindex;
    int        iAromatic;
    double     nbo, max_valence;
    double    *valence;
    char      *elem, *aname;
    char     **nb_hybrid, **at_ptr, **tp;
} t_atsel;

typedef struct {
    int       ai, aj, bcolor;
    int       nbbo;
    double   *bbo;
} t_bondsel;

static gmx_bool is_planar(rvec xi, rvec xj, rvec xk, rvec xl, t_pbc *pbc,
                          real phi_toler)
{
    int  t1, t2, t3;
    rvec r_ij, r_kj, r_kl, m, n;
    real sign, phi;

    phi = RAD2DEG*dih_angle(xi, xj, xk, xl, pbc, r_ij, r_kj, r_kl, m, n, &sign, &t1, &t2, &t3);

    return (fabs(phi) < phi_toler);
}

static gmx_bool is_linear(rvec xi, rvec xj, rvec xk, t_pbc *pbc,
                          real th_toler)
{
    int  t1, t2;
    rvec r_ij, r_kj;
    real costh, th;

    th = fabs(RAD2DEG*bond_angle(xi, xj, xk, pbc, r_ij, r_kj, &costh, &t1, &t2));

    return (th > th_toler) || (th < 180-th_toler);
}

static void calc_bondorder(FILE *fp,
                           gmx_poldata_t pd, t_atsel ats[], t_params *bonds,
                           double bondorder[])
{
    int     ai, aj, j, lu;
    double  toler, tol, dist, bo;
    char   *unit;

    unit = gmx_poldata_get_length_unit(pd);
    lu   = string2unit(unit);
    /* Tolerance is 0.002 nm */
    toler = gmx2convert(0.002, lu);

    for (j = 0; (j < bonds->nr); j++)
    {
        ai = bonds->param[j].a[0];
        aj = bonds->param[j].a[1];

        if ((NULL != ats[ai].tp[ats[ai].myat]) && (NULL != ats[aj].tp[ats[aj].myat]))
        {
            dist = gmx2convert(bonds->param[j].c[0], lu);
            tol  = toler;
            bo   = 0;
            while ((0 == bo) && (tol < 10*toler))
            {
                bo = gmx_poldata_atype_bondorder(pd,
                                                 ats[ai].tp[ats[ai].myat],
                                                 ats[aj].tp[ats[aj].myat], dist, tol);
                tol *= 2;
            }
            bondorder[j] = bo;
        }
        if (NULL != fp)
        {
            fprintf(fp, "Bond %s%d-%s%d order %g dist %g %s\n",
                    ats[ai].tp[ats[ai].myat], ai+1,
                    ats[aj].tp[ats[aj].myat], aj+1,
                    bondorder[j], dist, unit);
        }
    }
}

static void check_rings(FILE *fp, const char *molname, int natoms, t_atsel ats[])
{
    int iter, nra, nrp, i, k, bbk;

    for (iter = 0; (iter < 4); iter++)
    {
        nra = 0;
        for (i = 0; (i < natoms); i++)
        {
            if ((ats[i].nb == 2) && (ats[i].iAromatic > 0))
            {
                ats[i].iAromatic = ((0 < ats[ats[i].bbb[0]].iAromatic) &&
                                    (0 < ats[ats[i].bbb[1]].iAromatic));
            }
            else if ((ats[i].nb == 3) && (ats[i].geom == egmPlanar))
            {
                /* If two of three neighbors are in an aromatic ring, I am too */
                nrp = 0;
                for (k = 0; (k < 3); k++)
                {
                    bbk = ats[i].bbb[k];
                    if ((ats[bbk].geom == egmPlanar) &&
                        (0 < ats[bbk].iAromatic))
                    {
                        nrp++;
                    }
                }
                if (nrp < 2)
                {
                    ats[i].geom      = egmPlanar;
                    ats[i].iAromatic = 0;
                }
            }
            if ((ats[i].geom == egmPlanar)  && (ats[i].iAromatic > 0))
            {
                nra++;
            }
        }
    }
    if (fp)
    {
        fprintf(fp, "There are %d atoms in ring structures in %s.\n",
                nra, molname);
        for (i = 0; (i < natoms); i++)
        {
            fprintf(fp, "%d - geometry = %s (atom %s with %d bonds) iAromatic = %d\n",
                    i+1, geoms[ats[i].geom], ats[i].aname, ats[i].nb,
                    ats[i].iAromatic);
        }
    }
}

static int first_atom_colour(int fC, int natom, int Col, t_atsel ats[])
/* Return the first node with colour Col starting at fC.
 * return -1 if none found.
 */
{
    int i;

    for (i = fC; (i < natom); i++)
    {
        if ((ats[i].nb > 0) && (ats[i].ecol == Col))
        {
            return i;
        }
    }

    return -1;
}

static int colour_inner(int natoms, t_atsel ats[], int fW, int naG,
                        int nbonds, t_bondsel bondsel[],
                        double bondorder[])
{
    double *bo;
    int     bi, i, j, result;
    int     nbG = 0, nbW = 0, nbWi = 0;

    snew(bo, nbonds);
    for (i = 0; (i < nbonds); i++)
    {
        bo[i] = bondorder[i];
    }

    /* Now check all the possible bond orders for all the bonds outgoing from fW */
    for (i = 0; (i < ats[fW].nb); i++)
    {
        bi = ats[fW].bondindex[i];
        switch (bondsel[bi].bcolor)
        {
            case ecolWhite:
                bondsel[bi].bcolor = ecolGrey;
                result             = 0;
                for (j = 0; (j < bondsel[bi].nbbo) && (0 == result); j++)
                {
                    bo[bi] = bondsel[bi].bbo[j];
                    result = colour_inner(natoms, ats, fW, naG, nbonds, bondsel, bo);
                }
                if (0 != result)
                {
                    nbG++;
                    nbW--;
                    nbWi--;
                }
                else
                {
                    bondsel[bi].bcolor = ecolWhite;
                }
                break;
            case ecolGrey:

                break;
            case ecolBlack:
                if (ecolBlack != ats[bondsel[bi].ai].ecol)
                {
                    fprintf(stderr, "Can not modify atom %d\n", bondsel[bi].ai);
                    ats[bondsel[bi].ai].nbo += bo[bi];
                    if (ecolBlack == ats[bondsel[bi].aj].ecol)
                    {
                        gmx_fatal(FARGS, "Can not modify atom %d", bondsel[bi].aj);
                    }
                    ats[bondsel[bi].aj].nbo += bo[bi];
                    nbWi--;
                }
                break;
        }
    }
    if (0 != result)
    {
        for (i = 0; (i < nbonds); i++)
        {
            bondorder[i] = bo[i];
        }
    }
    sfree(bo);

    return result;
}

static void get_bondorders(FILE *fp, const char *molname, rvec x[],
                           int natoms, t_atsel ats[], t_params *bonds,
                           gmx_poldata_t pd, double bondorder[])
{
    int        lu;
    int        i, ai, aj, nbtot, nbonds;
    t_bondsel *bts;
    int        naW, naG, naB, nbW, nbG, nbB;
    int        fW, fG, *ibo_grey;
    char      *elem_i, *elem_j, *unit;
    rvec       dx;
    double     dx1, toler, *bo_grey;

    nbonds = bonds->nr;
    snew(bts, nbonds);
    for (i = 0; (i < nbonds); i++)
    {
        bts[i].ai     = bonds->param[i].a[0];
        bts[i].aj     = bonds->param[i].a[1];
        bts[i].bcolor = ecolWhite;
    }
    unit = gmx_poldata_get_length_unit(pd);
    lu   = string2unit(unit);
    snew(bo_grey, nbonds);
    snew(ibo_grey, nbonds);

    naW   = naG = naB = 0;
    nbW   = nbG = nbB = 0;
    nbtot = 0;
    for (i = 0; (i < natoms); i++)
    {
        ats[i].ecol = ecolWhite;
        naW++;
        nbtot += ats[i].nb;
        snew(ats[i].bondindex, ats[i].nb);
        ats[i].nb          = 0;
        ats[i].max_valence = gmx_poldata_elem_get_max_valence(pd, ats[i].elem);
    }
    /* Bonds are counted twice, so halve their number */
    nbtot /= 2;
    if (nbtot != nbonds)
    {
        gmx_fatal(FARGS, "Inconsistency in number of bonds");
    }
    /* Get the possible bondorders for each bond */
    for (i = 0; (i < nbonds); i++)
    {
        ai     = bts[i].ai;
        aj     = bts[i].aj;
        elem_i = ats[ai].elem;
        elem_j = ats[aj].elem;
        /* Store the index of this bond in the atoms structure */
        ats[ai].bondindex[ats[ai].nb++] = i;
        ats[aj].bondindex[ats[aj].nb++] = i;

        bts[i].bcolor = ecolWhite;
        rvec_sub(x[ai], x[aj], dx);
        dx1   = gmx2convert(norm(dx), lu);
        toler = 0.2*dx1;
        if (NULL == (bts[i].bbo = gmx_poldata_elem_get_bondorders(pd, elem_i, elem_j,
                                                                  dx1, toler)))
        {
            gmx_fatal(FARGS, "Can not find bond orders for %s - %s at distance %f",
                      elem_i, elem_j, dx1);
        }
        /* Count how many bonds */
        bts[i].nbbo = 0;
        while (0 != bts[i].bbo[bts[i].nbbo])
        {
            bts[i].nbbo++;
        }
        /* If only one possible make the bond black */
        if (1 == bts[i].nbbo)
        {
            bts[i].bcolor = ecolBlack;
            bondorder[i]  = bts[i].bbo[0];
            ats[ai].nbo  += bondorder[i];
            ats[aj].nbo  += bondorder[i];
            nbW--;
            nbB++;
        }
    }
    if (NULL != fp)
    {
        fprintf(fp, "There are %d out of %d bonds for which the bondorder is known in %s\n",
                nbB, nbtot, molname);
    }
    /* Loop over white atoms and bonds */
    fW = fG = 0;
    while (naW > 0)
    {
        /* Make atoms grey */
        if (-1 == (fW = first_atom_colour(fW, natoms, ecolWhite, ats)))
        {
            gmx_fatal(FARGS, "No WHITE atoms found while nW = %d\n", naW);
        }
        /* Make the first white atom grey */
        ats[fW].ecol = ecolGrey;
        naG++;
        naW--;

        while (naG > 0)
        {
            /* Check which valence to use! */
            if (0 < colour_inner(natoms, ats, fW, naG, nbonds, bts, bondorder))
            {
                if ((NULL == ats[fW].valence) ||
                    (ats[fW].nbo == ats[fW].valence[0]))
                {
                    naG--;
                    naB++;
                    ats[fW].ecol = ecolBlack;
                }
            }
        }
    }
}

static double minimize_valence(FILE *fp,
                               int natoms, t_atsel ats[], const char *molname,
                               t_params *bonds, gmx_poldata_t pd,
                               double bondorder[])
{
    int    iter, iter2, maxiter = 10, iAroma;
    int    i, j, ai, aj, ntot, jjj, nnull;
    double valence, ddval, dval, dval_best;
    char   buf[1024], b2[64];

    dval_best = natoms*4;
    for (iter = 0; (iter < maxiter) && (dval_best > 0); iter++)
    {
        check_rings(fp, molname, natoms, ats);
        nnull = 1;
        for (iter2 = 0; (iter2 < maxiter) && (nnull > 0); iter2++)
        {
            nnull = 0;
            for (i = 0; (i < natoms); i++)
            {
                for (j = 0; (j < ats[i].nb); j++)
                {
                    aj = ats[i].bbb[j];
                    if (NULL != ats[aj].at_ptr)
                    {
                        ats[i].nb_hybrid[j] = ats[aj].at_ptr[ats[aj].myat];
                    }
                    else
                    {
                        ats[i].nb_hybrid[j] = ats[aj].elem;
                    }
                }
                iAroma        = (iter == 0) ? -1 : ats[i].iAromatic;
                ats[i].at_ptr = gmx_poldata_get_bonding_rules(pd, ats[i].elem, ats[i].nb,
                                                              ats[i].nb_hybrid,
                                                              geoms[ats[i].geom],
                                                              iAroma);
                if ((NULL == ats[i].at_ptr) || (NULL == ats[i].at_ptr[0]))
                {
                    nnull++;
                }
            }
        }

        if (0 < nnull)
        {
            if (NULL != fp)
            {
                fprintf(fp, "Could not find bonding rules to describe all atoms in %s\n",
                        molname);
            }
            return 1;
        }
        /* Select the first possible type */
        ntot = 1;
        for (i = 0; (i < natoms); i++)
        {
            ats[i].nat = 0;
            while (NULL != ats[i].at_ptr[ats[i].nat])
            {
                srenew(ats[i].tp, ats[i].nat+1);
                ats[i].tp[ats[i].nat] = gmx_poldata_get_type(pd, ats[i].at_ptr[ats[i].nat]);
                if (0 == gmx_poldata_bonding_rule_valence(pd, ats[i].at_ptr[ats[i].nat], &valence))
                {
                    gmx_fatal(FARGS, "Can not find valence for %s",
                              ats[i].tp[ats[i].nat]);
                }
                srenew(ats[i].valence, ats[i].nat+1);
                ats[i].valence[ats[i].nat] = valence;
                ats[i].nat++;
            }
            if (0 == ats[i].nat)
            {
                return 1;
            }
            ats[i].myat = 0;
            ntot       *= ats[i].nat;
            ats[i].nn   = ntot;
        }
        if (NULL != fp)
        {
            fprintf(fp, "There are %d possible combinations of atomtypes in %s\n",
                    ntot, molname);
        }
        if (ntot < 10000)
        {
            /* Loop over all possibilities */
            for (j = 0; (j < ntot) && (dval_best > 0); j++)
            {
                /* Set the atomtypes */
                buf[0] = '\0';
                for (i = 0; (i < natoms); i++)
                {
                    if (i > 0)
                    {
                        jjj = (j/ats[i-1].nn);
                    }
                    else
                    {
                        jjj = j;
                    }
                    ats[i].myat = jjj % ats[i].nat;
                    sprintf(b2, " %s", ats[i].at_ptr[ats[i].myat]);
                    strcat(buf, b2);
                    ats[i].nbo = 0;
                }
                /* Compute the bond orders */
                calc_bondorder(fp, pd, ats, bonds, bondorder);

                /* Compute the valences based on bond order */
                for (i = 0; (i < bonds->nr); i++)
                {
                    ai           = bonds->param[i].a[0];
                    aj           = bonds->param[i].a[1];
                    ats[ai].nbo += bondorder[i];
                    ats[aj].nbo += bondorder[i];
                }
                /* Now check the total deviation from valence */
                dval = 0;
                for (i = 0; (i < natoms); i++)
                {
                    ddval = fabs(ats[i].nbo-ats[i].valence[ats[i].myat]);
                    dval += ddval;
                    if ((ddval > 0) && (NULL != fp))
                    {
                        fprintf(fp, "val for %s%d is %g (should be %g)\n",
                                ats[i].tp[ats[i].myat], i+1, ats[i].nbo,
                                ats[i].valence[ats[i].myat]);
                    }
                }
                if (NULL != fp)
                {
                    fprintf(fp, "j = %d myat = %s dval = %g\n", j, buf, dval);
                }
                if (dval < dval_best)
                {
                    dval_best = dval;
                }
            }
        }
    }
    if (0 < dval_best)
    {
        if (NULL != fp)
        {
            fprintf(fp, "Could not find suitable atomtypes for %s. Valence deviation is %g\n",
                    molname, dval_best);
        }
    }
    return dval_best;
}

int nm2type(FILE *fp, const char *molname, gmx_poldata_t pd, gmx_atomprop_t aps,
            t_symtab *tab, t_atoms *atoms, gmx_bool bRing[], double bondorder[],
            gpp_atomtype_t atype, int *nbonds, t_params *bonds,
            char **gt_atoms, rvec x[], t_pbc *pbc, real th_toler, real phi_toler,
            gentop_vsite_t gvt)
{
    int      i, j, nonebond, nresolved;
    int      ai, aj;
    int      type;
    char    *gt_type;
    char    *envptr;
    char    *empty_string = (char *) "";
    t_atsel *ats;
    double   mm, dval;
    real     value;
    t_atom  *atom;
    t_param *param;

    envptr = getenv("MOLNAME");
    if ((NULL != envptr) && (0 == strcasecmp(molname, envptr)))
    {
        printf("BOE! Molname = %s\n", molname);
    }

    snew(atom, 1);
    snew(param, 1);
    nonebond = 0;
    snew(ats, atoms->nr);
    for (i = 0; (i < atoms->nr); i++)
    {
        /* Atom information */
        atoms->atom[i].type = NOTSET;
        ats[i].aname        = *atoms->atomname[i];
        if ((ats[i].elem = gmx_atomprop_element(aps, atoms->atom[i].atomnumber)) == NULL)
        {
            ats[i].elem = empty_string;
        }
        ats[i].atomnr = atoms->atom[i].atomnumber;
        if ((ats[i].atomnr < 0) && (gmx_atomprop_query(aps, epropElement, "???",
                                                       ats[i].aname, &value)))
        {
            ats[i].atomnr             = gmx_nint(value);
            atoms->atom[i].atomnumber = ats[i].atomnr;
        }
        /* Bond information */
        if (nbonds[i] == 1)
        {
            nonebond++;
        }
        snew(ats[i].bbb, nbonds[i]);
        snew(ats[i].nb_hybrid, nbonds[i]);
        ats[i].nb  = 0;
    }
    /* Fill local bonding arrays */
    for (j = 0; (j < bonds->nr); j++)
    {
        ai = bonds->param[j].a[0];
        aj = bonds->param[j].a[1];
        ats[ai].bbb[ats[ai].nb++] = aj;
        ats[aj].bbb[ats[aj].nb++] = ai;
    }
    /* Check for consistency in number of bonds */
    for (i = 0; (i < atoms->nr); i++)
    {
        if (ats[i].nb != nbonds[i])
        {
            gmx_fatal(FARGS, "Inconsistency between number of nbonds[%d] = %d and bonds structure (%d)", i, nbonds[i], ats[i].nb);
        }
    }
    for (i = 0; (i < atoms->nr); i++)
    {
        /* Now test initial geometry */
        ats[i].geom = egmTetrahedral;
        if ((ats[i].nb == 2) && is_linear(x[i], x[ats[i].bbb[0]], x[ats[i].bbb[1]], pbc, th_toler))
        {
            ats[i].geom = egmLinear;
            if (NULL != gvt)
            {
                gentop_vsite_add_linear(gvt, ats[i].bbb[0], i, ats[i].bbb[1]);
            }
        }
        else if ((ats[i].nb == 3) && is_planar(x[i], x[ats[i].bbb[0]],
                                               x[ats[i].bbb[1]], x[ats[i].bbb[2]],
                                               pbc, phi_toler))
        {
            if (bRing[i])
            {
                ats[i].iAromatic = 1;
            }

            ats[i].geom = egmPlanar;
            if (NULL != gvt)
            {
                gentop_vsite_add_planar(gvt, i, ats[i].bbb[0], ats[i].bbb[1], ats[i].bbb[2],
                                        nbonds);
            }
        }
    }
    /* Iterate geometry test to look for rings. In a six-ring the information
     * spreads in at most three steps around the ring, so four iterations is sufficient.
     */
    if (fp)
    {
        fprintf(fp, "There are %d atoms with one bond in %s.\n",
                nonebond, molname);
    }

    /* Now we will determine for each atom the possible atom types in an iterative
     * fashion until convergence (deviation from valence criteria dval == 0).
     */
    dval = minimize_valence(fp, atoms->nr, ats, molname, bonds, pd, bondorder);
    /*get_bondorders(fp,molname,x,atoms->nr,ats,bonds,pd,bondorder);*/

    nresolved = 0;
    if (0 == dval)
    {
        /* Now set the contents of the atoms structure */
        for (i = 0; (i < atoms->nr); i++)
        {
            gt_atoms[i] = strdup(ats[i].at_ptr[ats[i].myat]);
            gt_type     = gmx_poldata_get_type(pd, gt_atoms[i]);
            mm          = 0;
            if (gmx_atomprop_query(aps, epropMass, "???", ats[i].elem, &value))
            {
                mm = value;
            }
            else
            {
                fprintf(stderr, "Can not find a mass for element %s", ats[i].elem);
            }

            type = get_atomtype_type(gt_type, atype);
            if (type == NOTSET)
            {
                /* Store mass and polarizability in atomtype as well */
                atom->qB = 0;
                atom->m  = mm;
                type     = add_atomtype(atype, tab, atom, gt_type, param,
                                        atom->type, 0, 0, 0, ats[i].atomnr, 0, 0);
            }

            /* atoms->atom[i].q  = 0;*/
            atoms->atom[i].qB    = 0;
            atoms->atom[i].m     = atoms->atom[i].mB = mm;
            atoms->atom[i].type  = type;
            atoms->atom[i].typeB = type;
            atoms->atomtype[i]   = put_symtab(tab, gt_type);

            nresolved++;
        }
    }

    /* fprintf(stderr,"\n");*/
    if (nresolved < atoms->nr)
    {
        char *missing, atom[32];
        snew(missing, 12*(atoms->nr-nresolved));
        for (i = 0; (i < atoms->nr); i++)
        {
            if (gt_atoms[i] == NULL)
            {
                sprintf(atom, " %s-%d", *atoms->atomname[i], i+1);
                strcat(missing, atom);
            }
        }
        if (fp)
        {
            fprintf(fp, "Could not resolve atomtypes %s for: %s.\n",
                    missing, molname);
        }
        sfree(missing);
    }
    sfree(atom);
    sfree(param);

    return nresolved;
}

gpp_atomtype_t set_atom_type(FILE *fp, const char *molname,
                             t_symtab *tab, t_atoms *atoms, t_params *bonds,
                             int nbonds[], gmx_bool bRing[], double bondorder[],
                             char **smnames, gmx_poldata_t pd,
                             gmx_atomprop_t aps, rvec x[], t_pbc *pbc, real th_toler,
                             real ph_toler, gentop_vsite_t gvt)
{
    gpp_atomtype_t atype;
    int            nresolved;
    int            i;

    atype = init_atomtype();
    snew(atoms->atomtype, atoms->nr);
    nresolved = nm2type(fp, molname, pd, aps, tab, atoms, bRing, bondorder, atype, nbonds,
                        bonds, smnames, x, pbc, th_toler, ph_toler, gvt);
    if (nresolved != atoms->nr)
    {
        return NULL;
    }
    else if (debug)
    {
        fprintf(debug, "There are %d different atom types in your sample\n",
                get_atomtype_ntypes(atype));
    }
    if (NULL == atoms->atomtype)
    {
        snew(atoms->atomtype, atoms->nr);
    }
    if (NULL == atoms->atomtypeB)
    {
        snew(atoms->atomtypeB, atoms->nr);
    }
    for (i = 0; (i < atoms->nr); i++)
    {
        atoms->atomtype[i]  = put_symtab(tab, get_atomtype_name(atoms->atom[i].type, atype));
        atoms->atomtypeB[i] = put_symtab(tab, get_atomtype_name(atoms->atom[i].typeB, atype));
    }

    return atype;
}
