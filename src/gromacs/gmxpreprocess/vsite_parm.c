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

#include "vsite_parm.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "gromacs/gmxpreprocess/add_par.h"
#include "gromacs/gmxpreprocess/resall.h"
#include "gromacs/gmxpreprocess/toputil.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/legacyheaders/names.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

typedef struct {
    t_iatom a[4];
    real    c;
} t_mybonded;

typedef struct {
    int      ftype;
    t_param *param;
} vsitebondparam_t;

typedef struct {
    int               nr;
    int               ftype;
    vsitebondparam_t *vsbp;
} at2vsitebond_t;

typedef struct {
    int  nr;
    int *aj;
} at2vsitecon_t;

static int vsite_bond_nrcheck(int ftype)
{
    int nrcheck;

    if ((interaction_function[ftype].flags & (IF_BTYPE | IF_CONSTRAINT | IF_ATYPE)) || (ftype == F_IDIHS))
    {
        nrcheck = NRAL(ftype);
    }
    else
    {
        nrcheck = 0;
    }

    return nrcheck;
}

static void enter_bonded(int nratoms, int *nrbonded, t_mybonded **bondeds,
                         t_param *param)
{
    int j;

    srenew(*bondeds, *nrbonded+1);

    /* copy atom numbers */
    for (j = 0; j < nratoms; j++)
    {
        (*bondeds)[*nrbonded].a[j] = param->a[j];
    }
    /* copy parameter */
    (*bondeds)[*nrbonded].c = param->C0;

    (*nrbonded)++;
}

static void get_bondeds(int nrat, t_iatom atoms[],
                        at2vsitebond_t *at2vb,
                        int *nrbond, t_mybonded **bonds,
                        int *nrang,  t_mybonded **angles,
                        int *nridih, t_mybonded **idihs )
{
    int      k, i, ftype, nrcheck;
    t_param *param;

    for (k = 0; k < nrat; k++)
    {
        for (i = 0; i < at2vb[atoms[k]].nr; i++)
        {
            ftype   = at2vb[atoms[k]].vsbp[i].ftype;
            param   = at2vb[atoms[k]].vsbp[i].param;
            nrcheck = vsite_bond_nrcheck(ftype);
            /* abuse nrcheck to see if we're adding bond, angle or idih */
            switch (nrcheck)
            {
                case 2: enter_bonded(nrcheck, nrbond, bonds, param); break;
                case 3: enter_bonded(nrcheck, nrang, angles, param); break;
                case 4: enter_bonded(nrcheck, nridih, idihs, param); break;
            }
        }
    }
}

static at2vsitebond_t *make_at2vsitebond(int natoms, t_params plist[])
{
    gmx_bool       *bVSI;
    int             ftype, i, j, nrcheck, nr;
    t_iatom        *aa;
    at2vsitebond_t *at2vb;

    snew(at2vb, natoms);

    snew(bVSI, natoms);
    for (ftype = 0; (ftype < F_NRE); ftype++)
    {
        if ((interaction_function[ftype].flags & IF_VSITE) && ftype != F_VSITEN)
        {
            for (i = 0; (i < plist[ftype].nr); i++)
            {
                for (j = 0; j < NRAL(ftype); j++)
                {
                    bVSI[plist[ftype].param[i].a[j]] = TRUE;
                }
            }
        }
    }

    for (ftype = 0; (ftype < F_NRE); ftype++)
    {
        nrcheck = vsite_bond_nrcheck(ftype);
        if (nrcheck > 0)
        {
            for (i = 0; (i < plist[ftype].nr); i++)
            {
                aa = plist[ftype].param[i].a;
                for (j = 0; j < nrcheck; j++)
                {
                    if (bVSI[aa[j]])
                    {
                        nr = at2vb[aa[j]].nr;
                        if (nr % 10 == 0)
                        {
                            srenew(at2vb[aa[j]].vsbp, nr+10);
                        }
                        at2vb[aa[j]].vsbp[nr].ftype = ftype;
                        at2vb[aa[j]].vsbp[nr].param = &plist[ftype].param[i];
                        at2vb[aa[j]].nr++;
                    }
                }
            }
        }
    }
    sfree(bVSI);

    return at2vb;
}

static void done_at2vsitebond(int natoms, at2vsitebond_t *at2vb)
{
    int i;

    for (i = 0; i < natoms; i++)
    {
        if (at2vb[i].nr)
        {
            sfree(at2vb[i].vsbp);
        }
    }
    sfree(at2vb);
}

static at2vsitecon_t *make_at2vsitecon(int natoms, t_params plist[])
{
    gmx_bool      *bVSI;
    int            ftype, i, j, ai, aj, nr;
    at2vsitecon_t *at2vc;

    snew(at2vc, natoms);

    snew(bVSI, natoms);
    for (ftype = 0; (ftype < F_NRE); ftype++)
    {
        if ((interaction_function[ftype].flags & IF_VSITE) && ftype != F_VSITEN)
        {
            for (i = 0; (i < plist[ftype].nr); i++)
            {
                for (j = 0; j < NRAL(ftype); j++)
                {
                    bVSI[plist[ftype].param[i].a[j]] = TRUE;
                }
            }
        }
    }

    for (ftype = 0; (ftype < F_NRE); ftype++)
    {
        if (interaction_function[ftype].flags & IF_CONSTRAINT)
        {
            for (i = 0; (i < plist[ftype].nr); i++)
            {
                ai = plist[ftype].param[i].AI;
                aj = plist[ftype].param[i].AJ;
                if (bVSI[ai] && bVSI[aj])
                {
                    /* Store forward direction */
                    nr = at2vc[ai].nr;
                    if (nr % 10 == 0)
                    {
                        srenew(at2vc[ai].aj, nr+10);
                    }
                    at2vc[ai].aj[nr] = aj;
                    at2vc[ai].nr++;
                    /* Store backward direction */
                    nr = at2vc[aj].nr;
                    if (nr % 10 == 0)
                    {
                        srenew(at2vc[aj].aj, nr+10);
                    }
                    at2vc[aj].aj[nr] = ai;
                    at2vc[aj].nr++;
                }
            }
        }
    }
    sfree(bVSI);

    return at2vc;
}

static void done_at2vsitecon(int natoms, at2vsitecon_t *at2vc)
{
    int i;

    for (i = 0; i < natoms; i++)
    {
        if (at2vc[i].nr)
        {
            sfree(at2vc[i].aj);
        }
    }
    sfree(at2vc);
}

/* for debug */
static void print_bad(FILE *fp,
                      int nrbond, t_mybonded *bonds,
                      int nrang,  t_mybonded *angles,
                      int nridih, t_mybonded *idihs )
{
    int i;

    if (nrbond)
    {
        fprintf(fp, "bonds:");
        for (i = 0; i < nrbond; i++)
        {
            fprintf(fp, " %d-%d (%g)",
                    bonds[i].AI+1, bonds[i].AJ+1, bonds[i].c);
        }
        fprintf(fp, "\n");
    }
    if (nrang)
    {
        fprintf(fp, "angles:");
        for (i = 0; i < nrang; i++)
        {
            fprintf(fp, " %d-%d-%d (%g)",
                    angles[i].AI+1, angles[i].AJ+1,
                    angles[i].AK+1, angles[i].c);
        }
        fprintf(fp, "\n");
    }
    if (nridih)
    {
        fprintf(fp, "idihs:");
        for (i = 0; i < nridih; i++)
        {
            fprintf(fp, " %d-%d-%d-%d (%g)",
                    idihs[i].AI+1, idihs[i].AJ+1,
                    idihs[i].AK+1, idihs[i].AL+1, idihs[i].c);
        }
        fprintf(fp, "\n");
    }
}

static void print_param(FILE *fp, int ftype, int i, t_param *param)
{
    static int pass       = 0;
    static int prev_ftype = NOTSET;
    static int prev_i     = NOTSET;
    int        j;

    if ( (ftype != prev_ftype) || (i != prev_i) )
    {
        pass       = 0;
        prev_ftype = ftype;
        prev_i     = i;
    }
    fprintf(fp, "(%d) plist[%s].param[%d]",
            pass, interaction_function[ftype].name, i);
    for (j = 0; j < NRFP(ftype); j++)
    {
        fprintf(fp, ".c[%d]=%g ", j, param->c[j]);
    }
    fprintf(fp, "\n");
    pass++;
}

static real get_bond_length(int nrbond, t_mybonded bonds[],
                            t_iatom ai, t_iatom aj)
{
    int  i;
    real bondlen;

    bondlen = NOTSET;
    for (i = 0; i < nrbond && (bondlen == NOTSET); i++)
    {
        /* check both ways */
        if ( ( (ai == bonds[i].AI) && (aj == bonds[i].AJ) ) ||
             ( (ai == bonds[i].AJ) && (aj == bonds[i].AI) ) )
        {
            bondlen = bonds[i].c; /* note: bonds[i].c might be NOTSET */
        }
    }
    return bondlen;
}

static real get_angle(int nrang, t_mybonded angles[],
                      t_iatom ai, t_iatom aj, t_iatom ak)
{
    int  i;
    real angle;

    angle = NOTSET;
    for (i = 0; i < nrang && (angle == NOTSET); i++)
    {
        /* check both ways */
        if ( ( (ai == angles[i].AI) && (aj == angles[i].AJ) && (ak == angles[i].AK) ) ||
             ( (ai == angles[i].AK) && (aj == angles[i].AJ) && (ak == angles[i].AI) ) )
        {
            angle = DEG2RAD*angles[i].c;
        }
    }
    return angle;
}

static char *get_atomtype_name_AB(t_atom *atom, gpp_atomtype_t atype)
{
    char *name;

    name = get_atomtype_name(atom->type, atype);

    /* When using the decoupling option, atom types are changed
     * to decoupled for the non-bonded interactions, but the virtual
     * sites constructions should be based on the "normal" interactions.
     * So we return the state B atom type name if the state A atom
     * type is the decoupled one. We should actually check for the atom
     * type number, but that's not passed here. So we check for
     * the decoupled atom type name. This should not cause trouble
     * as this code is only used for topologies with v-sites without
     * parameters generated by pdb2gmx.
     */
    if (strcmp(name, "decoupled") == 0)
    {
        name = get_atomtype_name(atom->typeB, atype);
    }

    return name;
}

static gmx_bool calc_vsite3_param(gpp_atomtype_t atype,
                                  t_param *param, t_atoms *at,
                                  int nrbond, t_mybonded *bonds,
                                  int nrang,  t_mybonded *angles )
{
    /* i = virtual site          |    ,k
     * j = 1st bonded heavy atom | i-j
     * k,l = 2nd bonded atoms    |    `l
     */

    gmx_bool bXH3, bError;
    real     bjk, bjl, a = -1, b = -1;
    /* check if this is part of a NH3 , NH2-umbrella or CH3 group,
     * i.e. if atom k and l are dummy masses (MNH* or MCH3*) */
    if (debug)
    {
        int i;
        for (i = 0; i < 4; i++)
        {
            fprintf(debug, "atom %d type %s ",
                    param->a[i]+1,
                    get_atomtype_name_AB(&at->atom[param->a[i]], atype));
        }
        fprintf(debug, "\n");
    }
    bXH3 =
        ( (gmx_strncasecmp(get_atomtype_name_AB(&at->atom[param->AK], atype), "MNH", 3) == 0) &&
          (gmx_strncasecmp(get_atomtype_name_AB(&at->atom[param->AL], atype), "MNH", 3) == 0) ) ||
        ( (gmx_strncasecmp(get_atomtype_name_AB(&at->atom[param->AK], atype), "MCH3", 4) == 0) &&
          (gmx_strncasecmp(get_atomtype_name_AB(&at->atom[param->AL], atype), "MCH3", 4) == 0) );

    bjk    = get_bond_length(nrbond, bonds, param->AJ, param->AK);
    bjl    = get_bond_length(nrbond, bonds, param->AJ, param->AL);
    bError = (bjk == NOTSET) || (bjl == NOTSET);
    if (bXH3)
    {
        /* now we get some XH2/XH3 group specific construction */
        /* note: we call the heavy atom 'C' and the X atom 'N' */
        real bMM, bCM, bCN, bNH, aCNH, dH, rH, dM, rM;
        int  aN;

        /* check if bonds from heavy atom (j) to dummy masses (k,l) are equal: */
        bError = bError || (bjk != bjl);

        /* the X atom (C or N) in the XH2/XH3 group is the first after the masses: */
        aN = max(param->AK, param->AL)+1;

        /* get common bonds */
        bMM    = get_bond_length(nrbond, bonds, param->AK, param->AL);
        bCM    = bjk;
        bCN    = get_bond_length(nrbond, bonds, param->AJ, aN);
        bError = bError || (bMM == NOTSET) || (bCN == NOTSET);

        /* calculate common things */
        rM  = 0.5*bMM;
        dM  = sqrt( sqr(bCM) - sqr(rM) );

        /* are we dealing with the X atom? */
        if (param->AI == aN)
        {
            /* this is trivial */
            a = b = 0.5 * bCN/dM;

        }
        else
        {
            /* get other bondlengths and angles: */
            bNH    = get_bond_length(nrbond, bonds, aN, param->AI);
            aCNH   = get_angle      (nrang, angles, param->AJ, aN, param->AI);
            bError = bError || (bNH == NOTSET) || (aCNH == NOTSET);

            /* calculate */
            dH  = bCN - bNH * cos(aCNH);
            rH  = bNH * sin(aCNH);

            a = 0.5 * ( dH/dM + rH/rM );
            b = 0.5 * ( dH/dM - rH/rM );
        }
    }
    else
    {
        gmx_fatal(FARGS, "calc_vsite3_param not implemented for the general case "
                  "(atom %d)", param->AI+1);
    }

    param->C0 = a;
    param->C1 = b;

    if (debug)
    {
        fprintf(debug, "params for vsite3 %d: %g %g\n",
                param->AI+1, param->C0, param->C1);
    }

    return bError;
}

static gmx_bool calc_vsite3fd_param(t_param *param,
                                    int nrbond, t_mybonded *bonds,
                                    int nrang,  t_mybonded *angles)
{
    /* i = virtual site          |    ,k
     * j = 1st bonded heavy atom | i-j
     * k,l = 2nd bonded atoms    |    `l
     */

    gmx_bool bError;
    real     bij, bjk, bjl, aijk, aijl, rk, rl;

    bij    = get_bond_length(nrbond, bonds, param->AI, param->AJ);
    bjk    = get_bond_length(nrbond, bonds, param->AJ, param->AK);
    bjl    = get_bond_length(nrbond, bonds, param->AJ, param->AL);
    aijk   = get_angle      (nrang, angles, param->AI, param->AJ, param->AK);
    aijl   = get_angle      (nrang, angles, param->AI, param->AJ, param->AL);
    bError = (bij == NOTSET) || (bjk == NOTSET) || (bjl == NOTSET) ||
        (aijk == NOTSET) || (aijl == NOTSET);

    rk        = bjk * sin(aijk);
    rl        = bjl * sin(aijl);
    param->C0 = rk / (rk + rl);
    param->C1 = -bij; /* 'bond'-length for fixed distance vsite */

    if (debug)
    {
        fprintf(debug, "params for vsite3fd %d: %g %g\n",
                param->AI+1, param->C0, param->C1);
    }
    return bError;
}

static gmx_bool calc_vsite3fad_param(t_param *param,
                                     int nrbond, t_mybonded *bonds,
                                     int nrang,  t_mybonded *angles)
{
    /* i = virtual site          |
     * j = 1st bonded heavy atom | i-j
     * k = 2nd bonded heavy atom |    `k-l
     * l = 3d bonded heavy atom  |
     */

    gmx_bool bSwapParity, bError;
    real     bij, aijk;

    bSwapParity = ( param->C1 == -1 );

    bij    = get_bond_length(nrbond, bonds, param->AI, param->AJ);
    aijk   = get_angle      (nrang, angles, param->AI, param->AJ, param->AK);
    bError = (bij == NOTSET) || (aijk == NOTSET);

    param->C1 = bij;          /* 'bond'-length for fixed distance vsite */
    param->C0 = RAD2DEG*aijk; /* 'bond'-angle for fixed angle vsite */

    if (bSwapParity)
    {
        param->C0 = 360 - param->C0;
    }

    if (debug)
    {
        fprintf(debug, "params for vsite3fad %d: %g %g\n",
                param->AI+1, param->C0, param->C1);
    }
    return bError;
}

static gmx_bool calc_vsite3out_param(gpp_atomtype_t atype,
                                     t_param *param, t_atoms *at,
                                     int nrbond, t_mybonded *bonds,
                                     int nrang,  t_mybonded *angles)
{
    /* i = virtual site          |    ,k
     * j = 1st bonded heavy atom | i-j
     * k,l = 2nd bonded atoms    |    `l
     * NOTE: i is out of the j-k-l plane!
     */

    gmx_bool bXH3, bError, bSwapParity;
    real     bij, bjk, bjl, aijk, aijl, akjl, pijk, pijl, a, b, c;

    /* check if this is part of a NH2-umbrella, NH3 or CH3 group,
     * i.e. if atom k and l are dummy masses (MNH* or MCH3*) */
    if (debug)
    {
        int i;
        for (i = 0; i < 4; i++)
        {
            fprintf(debug, "atom %d type %s ",
                    param->a[i]+1, get_atomtype_name_AB(&at->atom[param->a[i]], atype));
        }
        fprintf(debug, "\n");
    }
    bXH3 =
        ( (gmx_strncasecmp(get_atomtype_name_AB(&at->atom[param->AK], atype), "MNH", 3) == 0) &&
          (gmx_strncasecmp(get_atomtype_name_AB(&at->atom[param->AL], atype), "MNH", 3) == 0) ) ||
        ( (gmx_strncasecmp(get_atomtype_name_AB(&at->atom[param->AK], atype), "MCH3", 4) == 0) &&
          (gmx_strncasecmp(get_atomtype_name_AB(&at->atom[param->AL], atype), "MCH3", 4) == 0) );

    /* check if construction parity must be swapped */
    bSwapParity = ( param->C1 == -1 );

    bjk    = get_bond_length(nrbond, bonds, param->AJ, param->AK);
    bjl    = get_bond_length(nrbond, bonds, param->AJ, param->AL);
    bError = (bjk == NOTSET) || (bjl == NOTSET);
    if (bXH3)
    {
        /* now we get some XH3 group specific construction */
        /* note: we call the heavy atom 'C' and the X atom 'N' */
        real bMM, bCM, bCN, bNH, aCNH, dH, rH, rHx, rHy, dM, rM;
        int  aN;

        /* check if bonds from heavy atom (j) to dummy masses (k,l) are equal: */
        bError = bError || (bjk != bjl);

        /* the X atom (C or N) in the XH3 group is the first after the masses: */
        aN = max(param->AK, param->AL)+1;

        /* get all bondlengths and angles: */
        bMM    = get_bond_length(nrbond, bonds, param->AK, param->AL);
        bCM    = bjk;
        bCN    = get_bond_length(nrbond, bonds, param->AJ, aN);
        bNH    = get_bond_length(nrbond, bonds, aN, param->AI);
        aCNH   = get_angle      (nrang, angles, param->AJ, aN, param->AI);
        bError = bError ||
            (bMM == NOTSET) || (bCN == NOTSET) || (bNH == NOTSET) || (aCNH == NOTSET);

        /* calculate */
        dH  = bCN - bNH * cos(aCNH);
        rH  = bNH * sin(aCNH);
        /* we assume the H's are symmetrically distributed */
        rHx = rH*cos(DEG2RAD*30);
        rHy = rH*sin(DEG2RAD*30);
        rM  = 0.5*bMM;
        dM  = sqrt( sqr(bCM) - sqr(rM) );
        a   = 0.5*( (dH/dM) - (rHy/rM) );
        b   = 0.5*( (dH/dM) + (rHy/rM) );
        c   = rHx / (2*dM*rM);

    }
    else
    {
        /* this is the general construction */

        bij    = get_bond_length(nrbond, bonds, param->AI, param->AJ);
        aijk   = get_angle      (nrang, angles, param->AI, param->AJ, param->AK);
        aijl   = get_angle      (nrang, angles, param->AI, param->AJ, param->AL);
        akjl   = get_angle      (nrang, angles, param->AK, param->AJ, param->AL);
        bError = bError ||
            (bij == NOTSET) || (aijk == NOTSET) || (aijl == NOTSET) || (akjl == NOTSET);

        pijk = cos(aijk)*bij;
        pijl = cos(aijl)*bij;
        a    = ( pijk + (pijk*cos(akjl)-pijl) * cos(akjl) / sqr(sin(akjl)) ) / bjk;
        b    = ( pijl + (pijl*cos(akjl)-pijk) * cos(akjl) / sqr(sin(akjl)) ) / bjl;
        c    = -sqrt( sqr(bij) -
                      ( sqr(pijk) - 2*pijk*pijl*cos(akjl) + sqr(pijl) )
                      / sqr(sin(akjl)) )
            / ( bjk*bjl*sin(akjl) );
    }

    param->C0 = a;
    param->C1 = b;
    if (bSwapParity)
    {
        param->C2 = -c;
    }
    else
    {
        param->C2 =  c;
    }
    if (debug)
    {
        fprintf(debug, "params for vsite3out %d: %g %g %g\n",
                param->AI+1, param->C0, param->C1, param->C2);
    }
    return bError;
}

static gmx_bool calc_vsite4fd_param(t_param *param,
                                    int nrbond, t_mybonded *bonds,
                                    int nrang,  t_mybonded *angles)
{
    /* i = virtual site          |    ,k
     * j = 1st bonded heavy atom | i-j-m
     * k,l,m = 2nd bonded atoms  |    `l
     */

    gmx_bool bError;
    real     bij, bjk, bjl, bjm, aijk, aijl, aijm, akjm, akjl;
    real     pk, pl, pm, cosakl, cosakm, sinakl, sinakm, cl, cm;

    bij    = get_bond_length(nrbond, bonds, param->AI, param->AJ);
    bjk    = get_bond_length(nrbond, bonds, param->AJ, param->AK);
    bjl    = get_bond_length(nrbond, bonds, param->AJ, param->AL);
    bjm    = get_bond_length(nrbond, bonds, param->AJ, param->AM);
    aijk   = get_angle      (nrang, angles, param->AI, param->AJ, param->AK);
    aijl   = get_angle      (nrang, angles, param->AI, param->AJ, param->AL);
    aijm   = get_angle      (nrang, angles, param->AI, param->AJ, param->AM);
    akjm   = get_angle      (nrang, angles, param->AK, param->AJ, param->AM);
    akjl   = get_angle      (nrang, angles, param->AK, param->AJ, param->AL);
    bError = (bij == NOTSET) || (bjk == NOTSET) || (bjl == NOTSET) || (bjm == NOTSET) ||
        (aijk == NOTSET) || (aijl == NOTSET) || (aijm == NOTSET) || (akjm == NOTSET) ||
        (akjl == NOTSET);

    if (!bError)
    {
        pk     = bjk*sin(aijk);
        pl     = bjl*sin(aijl);
        pm     = bjm*sin(aijm);
        cosakl = (cos(akjl) - cos(aijk)*cos(aijl)) / (sin(aijk)*sin(aijl));
        cosakm = (cos(akjm) - cos(aijk)*cos(aijm)) / (sin(aijk)*sin(aijm));
        if (cosakl < -1 || cosakl > 1 || cosakm < -1 || cosakm > 1)
        {
            fprintf(stderr, "virtual site %d: angle ijk = %f, angle ijl = %f, angle ijm = %f\n",
                    param->AI+1, RAD2DEG*aijk, RAD2DEG*aijl, RAD2DEG*aijm);
            gmx_fatal(FARGS, "invalid construction in calc_vsite4fd for atom %d: "
                      "cosakl=%f, cosakm=%f\n", param->AI+1, cosakl, cosakm);
        }
        sinakl = sqrt(1-sqr(cosakl));
        sinakm = sqrt(1-sqr(cosakm));

        /* note: there is a '+' because of the way the sines are calculated */
        cl = -pk / ( pl*cosakl - pk + pl*sinakl*(pm*cosakm-pk)/(pm*sinakm) );
        cm = -pk / ( pm*cosakm - pk + pm*sinakm*(pl*cosakl-pk)/(pl*sinakl) );

        param->C0 = cl;
        param->C1 = cm;
        param->C2 = -bij;
        if (debug)
        {
            fprintf(debug, "params for vsite4fd %d: %g %g %g\n",
                    param->AI+1, param->C0, param->C1, param->C2);
        }
    }

    return bError;
}


static gmx_bool
calc_vsite4fdn_param(t_param *param,
                     int nrbond, t_mybonded *bonds,
                     int nrang,  t_mybonded *angles)
{
    /* i = virtual site          |    ,k
     * j = 1st bonded heavy atom | i-j-m
     * k,l,m = 2nd bonded atoms  |    `l
     */

    gmx_bool bError;
    real     bij, bjk, bjl, bjm, aijk, aijl, aijm;
    real     pk, pl, pm, a, b;

    bij  = get_bond_length(nrbond, bonds, param->AI, param->AJ);
    bjk  = get_bond_length(nrbond, bonds, param->AJ, param->AK);
    bjl  = get_bond_length(nrbond, bonds, param->AJ, param->AL);
    bjm  = get_bond_length(nrbond, bonds, param->AJ, param->AM);
    aijk = get_angle      (nrang, angles, param->AI, param->AJ, param->AK);
    aijl = get_angle      (nrang, angles, param->AI, param->AJ, param->AL);
    aijm = get_angle      (nrang, angles, param->AI, param->AJ, param->AM);

    bError = (bij == NOTSET) || (bjk == NOTSET) || (bjl == NOTSET) || (bjm == NOTSET) ||
        (aijk == NOTSET) || (aijl == NOTSET) || (aijm == NOTSET);

    if (!bError)
    {

        /* Calculate component of bond j-k along the direction i-j */
        pk = -bjk*cos(aijk);

        /* Calculate component of bond j-l along the direction i-j */
        pl = -bjl*cos(aijl);

        /* Calculate component of bond j-m along the direction i-j */
        pm = -bjm*cos(aijm);

        if (fabs(pl) < 1000*GMX_REAL_MIN || fabs(pm) < 1000*GMX_REAL_MIN)
        {
            fprintf(stderr, "virtual site %d: angle ijk = %f, angle ijl = %f, angle ijm = %f\n",
                    param->AI+1, RAD2DEG*aijk, RAD2DEG*aijl, RAD2DEG*aijm);
            gmx_fatal(FARGS, "invalid construction in calc_vsite4fdn for atom %d: "
                      "pl=%f, pm=%f\n", param->AI+1, pl, pm);
        }

        a = pk/pl;
        b = pk/pm;

        param->C0 = a;
        param->C1 = b;
        param->C2 = bij;

        if (debug)
        {
            fprintf(debug, "params for vsite4fdn %d: %g %g %g\n",
                    param->AI+1, param->C0, param->C1, param->C2);
        }
    }

    return bError;
}



int set_vsites(gmx_bool bVerbose, t_atoms *atoms, gpp_atomtype_t atype,
               t_params plist[])
{
    int             i, j, ftype;
    int             nvsite, nrbond, nrang, nridih, nrset;
    gmx_bool        bFirst, bSet, bERROR;
    at2vsitebond_t *at2vb;
    t_mybonded     *bonds;
    t_mybonded     *angles;
    t_mybonded     *idihs;

    bFirst = TRUE;
    bERROR = TRUE;
    nvsite = 0;
    if (debug)
    {
        fprintf(debug, "\nCalculating parameters for virtual sites\n");
    }

    /* Make a reverse list to avoid ninteractions^2 operations */
    at2vb = make_at2vsitebond(atoms->nr, plist);

    for (ftype = 0; (ftype < F_NRE); ftype++)
    {
        if (interaction_function[ftype].flags & IF_VSITE)
        {
            nvsite += plist[ftype].nr;

            if (ftype == F_VSITEN)
            {
                /* We don't calculate parameters for VSITEN */
                continue;
            }

            nrset = 0;
            for (i = 0; (i < plist[ftype].nr); i++)
            {
                /* check if all parameters are set */
                bSet = TRUE;
                for (j = 0; j < NRFP(ftype) && bSet; j++)
                {
                    bSet = plist[ftype].param[i].c[j] != NOTSET;
                }

                if (debug)
                {
                    fprintf(debug, "bSet=%s ", bool_names[bSet]);
                    print_param(debug, ftype, i, &plist[ftype].param[i]);
                }
                if (!bSet)
                {
                    if (bVerbose && bFirst)
                    {
                        fprintf(stderr, "Calculating parameters for virtual sites\n");
                        bFirst = FALSE;
                    }

                    nrbond = nrang = nridih = 0;
                    bonds  = NULL;
                    angles = NULL;
                    idihs  = NULL;
                    nrset++;
                    /* now set the vsite parameters: */
                    get_bondeds(NRAL(ftype), plist[ftype].param[i].a, at2vb,
                                &nrbond, &bonds, &nrang,  &angles, &nridih, &idihs);
                    if (debug)
                    {
                        fprintf(debug, "Found %d bonds, %d angles and %d idihs "
                                "for virtual site %d (%s)\n", nrbond, nrang, nridih,
                                plist[ftype].param[i].AI+1,
                                interaction_function[ftype].longname);
                        print_bad(debug, nrbond, bonds, nrang, angles, nridih, idihs);
                    } /* debug */
                    switch (ftype)
                    {
                        case F_VSITE3:
                            bERROR =
                                calc_vsite3_param(atype, &(plist[ftype].param[i]), atoms,
                                                  nrbond, bonds, nrang, angles);
                            break;
                        case F_VSITE3FD:
                            bERROR =
                                calc_vsite3fd_param(&(plist[ftype].param[i]),
                                                    nrbond, bonds, nrang, angles);
                            break;
                        case F_VSITE3FAD:
                            bERROR =
                                calc_vsite3fad_param(&(plist[ftype].param[i]),
                                                     nrbond, bonds, nrang, angles);
                            break;
                        case F_VSITE3OUT:
                            bERROR =
                                calc_vsite3out_param(atype, &(plist[ftype].param[i]), atoms,
                                                     nrbond, bonds, nrang, angles);
                            break;
                        case F_VSITE4FD:
                            bERROR =
                                calc_vsite4fd_param(&(plist[ftype].param[i]),
                                                    nrbond, bonds, nrang, angles);
                            break;
                        case F_VSITE4FDN:
                            bERROR =
                                calc_vsite4fdn_param(&(plist[ftype].param[i]),
                                                     nrbond, bonds, nrang, angles);
                            break;
                        default:
                            gmx_fatal(FARGS, "Automatic parameter generation not supported "
                                      "for %s atom %d",
                                      interaction_function[ftype].longname,
                                      plist[ftype].param[i].AI+1);
                    } /* switch */
                    if (bERROR)
                    {
                        gmx_fatal(FARGS, "Automatic parameter generation not supported "
                                  "for %s atom %d for this bonding configuration",
                                  interaction_function[ftype].longname,
                                  plist[ftype].param[i].AI+1);
                    }
                    sfree(bonds);
                    sfree(angles);
                    sfree(idihs);
                } /* if bSet */
            }     /* for i */
            if (debug && plist[ftype].nr)
            {
                fprintf(stderr, "Calculated parameters for %d out of %d %s atoms\n",
                        nrset, plist[ftype].nr, interaction_function[ftype].longname);
            }
        } /* if IF_VSITE */

    }
    done_at2vsitebond(atoms->nr, at2vb);

    return nvsite;
}

void set_vsites_ptype(gmx_bool bVerbose, gmx_moltype_t *molt)
{
    int      ftype, i;
    int      nra, nrd;
    t_ilist *il;
    t_iatom *ia, avsite;

    if (bVerbose)
    {
        fprintf(stderr, "Setting particle type to V for virtual sites\n");
    }
    if (debug)
    {
        fprintf(stderr, "checking %d functypes\n", F_NRE);
    }
    for (ftype = 0; ftype < F_NRE; ftype++)
    {
        il = &molt->ilist[ftype];
        if (interaction_function[ftype].flags & IF_VSITE)
        {
            nra    = interaction_function[ftype].nratoms;
            nrd    = il->nr;
            ia     = il->iatoms;

            if (debug && nrd)
            {
                fprintf(stderr, "doing %d %s virtual sites\n",
                        (nrd / (nra+1)), interaction_function[ftype].longname);
            }

            for (i = 0; (i < nrd); )
            {
                /* The virtual site */
                avsite = ia[1];
                molt->atoms.atom[avsite].ptype = eptVSite;

                i  += nra+1;
                ia += nra+1;
            }
        }
    }

}

typedef struct {
    int ftype, parnr;
} t_pindex;

static void check_vsite_constraints(t_params *plist,
                                    int cftype, int vsite_type[])
{
    int       i, k, n;
    atom_id   atom;
    t_params *ps;

    n  = 0;
    ps = &(plist[cftype]);
    for (i = 0; (i < ps->nr); i++)
    {
        for (k = 0; k < 2; k++)
        {
            atom = ps->param[i].a[k];
            if (vsite_type[atom] != NOTSET)
            {
                fprintf(stderr, "ERROR: Cannot have constraint (%d-%d) with virtual site (%d)\n",
                        ps->param[i].AI+1, ps->param[i].AJ+1, atom+1);
                n++;
            }
        }
    }
    if (n)
    {
        gmx_fatal(FARGS, "There were %d virtual sites involved in constraints", n);
    }
}

static void clean_vsite_bonds(t_params *plist, t_pindex pindex[],
                              int cftype, int vsite_type[])
{
    int          ftype, i, j, parnr, k, l, m, n, nvsite, nOut, kept_i, vsitetype;
    int          nconverted, nremoved;
    atom_id      atom, oatom, at1, at2;
    gmx_bool     bKeep, bRemove, bUsed, bPresent, bThisFD, bThisOUT, bAllFD, bFirstTwo;
    t_params    *ps;

    if (cftype == F_CONNBONDS)
    {
        return;
    }

    ps         = &(plist[cftype]);
    kept_i     = 0;
    nconverted = 0;
    nremoved   = 0;
    nOut       = 0;
    for (i = 0; (i < ps->nr); i++) /* for all bonds in the plist */
    {
        int            vsnral      = 0;
        const atom_id *first_atoms = NULL;

        bKeep   = FALSE;
        bRemove = FALSE;
        bAllFD  = TRUE;
        /* check if all virtual sites are constructed from the same atoms */
        nvsite = 0;
        if (debug)
        {
            fprintf(debug, "constr %d %d:", ps->param[i].AI+1, ps->param[i].AJ+1);
        }
        for (k = 0; (k < 2) && !bKeep && !bRemove; k++)
        {
            /* for all atoms in the bond */
            atom = ps->param[i].a[k];
            if (vsite_type[atom] != NOTSET && vsite_type[atom] != F_VSITEN)
            {
                nvsite++;
                bThisFD = ( (pindex[atom].ftype == F_VSITE3FD ) ||
                            (pindex[atom].ftype == F_VSITE3FAD) ||
                            (pindex[atom].ftype == F_VSITE4FD ) ||
                            (pindex[atom].ftype == F_VSITE4FDN ) );
                bThisOUT = ( (pindex[atom].ftype == F_VSITE3OUT) &&
                             (interaction_function[cftype].flags & IF_CONSTRAINT) );
                bAllFD = bAllFD && bThisFD;
                if (bThisFD || bThisOUT)
                {
                    if (debug)
                    {
                        fprintf(debug, " %s", bThisOUT ? "out" : "fd");
                    }
                    oatom = ps->param[i].a[1-k]; /* the other atom */
                    if (vsite_type[oatom] == NOTSET &&
                        vsite_type[oatom] != F_VSITEN &&
                        oatom == plist[pindex[atom].ftype].param[pindex[atom].parnr].AJ)
                    {
                        /* if the other atom isn't a vsite, and it is AI */
                        bRemove = TRUE;
                        if (bThisOUT)
                        {
                            nOut++;
                        }
                        if (debug)
                        {
                            fprintf(debug, " D-AI");
                        }
                    }
                }
                if (!bRemove)
                {
                    /* TODO This fragment, and corresponding logic in
                       clean_vsite_angles and clean_vsite_dihs should
                       be refactored into a common function */
                    if (nvsite == 1)
                    {
                        /* if this is the first vsite we encounter then
                           store construction atoms */
                        /* TODO This would be nicer to implement with
                           a C++ "vector view" class" with an
                           STL-container-like interface. */
                        vsnral      = NRAL(pindex[atom].ftype) - 1;
                        first_atoms = plist[pindex[atom].ftype].param[pindex[atom].parnr].a + 1;
                    }
                    else
                    {
                        assert(vsnral != 0);
                        assert(first_atoms != NULL);
                        /* if it is not the first then
                           check if this vsite is constructed from the same atoms */
                        if (vsnral == NRAL(pindex[atom].ftype)-1)
                        {
                            for (m = 0; (m < vsnral) && !bKeep; m++)
                            {
                                const atom_id *atoms;

                                bPresent = FALSE;
                                atoms    = plist[pindex[atom].ftype].param[pindex[atom].parnr].a + 1;
                                for (n = 0; (n < vsnral) && !bPresent; n++)
                                {
                                    if (atoms[m] == first_atoms[n])
                                    {
                                        bPresent = TRUE;
                                    }
                                }
                                if (!bPresent)
                                {
                                    bKeep = TRUE;
                                    if (debug)
                                    {
                                        fprintf(debug, " !present");
                                    }
                                }
                            }
                        }
                        else
                        {
                            bKeep = TRUE;
                            if (debug)
                            {
                                fprintf(debug, " !same#at");
                            }
                        }
                    }
                }
            }
        }

        if (bRemove)
        {
            bKeep = FALSE;
        }
        else
        {
            /* if we have no virtual sites in this bond, keep it */
            if (nvsite == 0)
            {
                if (debug)
                {
                    fprintf(debug, " no vsite");
                }
                bKeep = TRUE;
            }

            /* TODO This loop and the corresponding loop in
               check_vsite_angles should be refactored into a common
               function */
            /* check if all non-vsite atoms are used in construction: */
            bFirstTwo = TRUE;
            for (k = 0; (k < 2) && !bKeep; k++) /* for all atoms in the bond */
            {
                atom = ps->param[i].a[k];
                if (vsite_type[atom] == NOTSET && vsite_type[atom] != F_VSITEN)
                {
                    bUsed = FALSE;
                    for (m = 0; (m < vsnral) && !bUsed; m++)
                    {
                        assert(first_atoms != NULL);

                        if (atom == first_atoms[m])
                        {
                            bUsed     = TRUE;
                            bFirstTwo = bFirstTwo && m < 2;
                        }
                    }
                    if (!bUsed)
                    {
                        bKeep = TRUE;
                        if (debug)
                        {
                            fprintf(debug, " !used");
                        }
                    }
                }
            }

            if (!( bAllFD && bFirstTwo ) )
            {
                /* Two atom bonded interactions include constraints.
                 * We need to remove constraints between vsite pairs that have
                 * a fixed distance due to being constructed from the same
                 * atoms, since this can be numerically unstable.
                 */
                for (m = 0; m < vsnral && !bKeep; m++) /* all constr. atoms */
                {
                    at1      = first_atoms[m];
                    at2      = first_atoms[(m+1) % vsnral];
                    bPresent = FALSE;
                    for (ftype = 0; ftype < F_NRE; ftype++)
                    {
                        if (interaction_function[ftype].flags & IF_CONSTRAINT)
                        {
                            for (j = 0; (j < plist[ftype].nr) && !bPresent; j++)
                            {
                                /* all constraints until one matches */
                                bPresent = ( ( (plist[ftype].param[j].AI == at1) &&
                                               (plist[ftype].param[j].AJ == at2) ) ||
                                             ( (plist[ftype].param[j].AI == at2) &&
                                               (plist[ftype].param[j].AJ == at1) ) );
                            }
                        }
                    }
                    if (!bPresent)
                    {
                        bKeep = TRUE;
                        if (debug)
                        {
                            fprintf(debug, " !bonded");
                        }
                    }
                }
            }
        }

        if (bKeep)
        {
            if (debug)
            {
                fprintf(debug, " keeping");
            }
            /* now copy the bond to the new array */
            ps->param[kept_i] = ps->param[i];
            kept_i++;
        }
        else if (IS_CHEMBOND(cftype))
        {
            srenew(plist[F_CONNBONDS].param, plist[F_CONNBONDS].nr+1);
            plist[F_CONNBONDS].param[plist[F_CONNBONDS].nr] = ps->param[i];
            plist[F_CONNBONDS].nr++;
            nconverted++;
        }
        else
        {
            nremoved++;
        }
        if (debug)
        {
            fprintf(debug, "\n");
        }
    }

    if (nremoved)
    {
        fprintf(stderr, "Removed   %4d %15ss with virtual sites, %5d left\n",
                nremoved, interaction_function[cftype].longname, kept_i);
    }
    if (nconverted)
    {
        fprintf(stderr, "Converted %4d %15ss with virtual sites to connections, %5d left\n",
                nconverted, interaction_function[cftype].longname, kept_i);
    }
    if (nOut)
    {
        fprintf(stderr, "Warning: removed %d %ss with vsite with %s construction\n"
                "         This vsite construction does not guarantee constant "
                "bond-length\n"
                "         If the constructions were generated by pdb2gmx ignore "
                "this warning\n",
                nOut, interaction_function[cftype].longname,
                interaction_function[F_VSITE3OUT].longname );
    }
    ps->nr = kept_i;
}

static void clean_vsite_angles(t_params *plist, t_pindex pindex[],
                               int cftype, int vsite_type[],
                               at2vsitecon_t *at2vc)
{
    int          i, j, parnr, k, l, m, n, nvsite, kept_i, vsitetype;
    atom_id      atom, at1, at2;
    gmx_bool     bKeep, bUsed, bPresent, bAll3FAD, bFirstTwo;
    t_params    *ps;

    ps     = &(plist[cftype]);
    kept_i = 0;
    for (i = 0; (i < ps->nr); i++) /* for all angles in the plist */
    {
        int            vsnral      = 0;
        const atom_id *first_atoms = NULL;

        bKeep    = FALSE;
        bAll3FAD = TRUE;
        /* check if all virtual sites are constructed from the same atoms */
        nvsite = 0;
        for (k = 0; (k < 3) && !bKeep; k++) /* for all atoms in the angle */
        {
            atom = ps->param[i].a[k];
            if (vsite_type[atom] != NOTSET && vsite_type[atom] != F_VSITEN)
            {
                nvsite++;
                bAll3FAD = bAll3FAD && (pindex[atom].ftype == F_VSITE3FAD);
                if (nvsite == 1)
                {
                    /* store construction atoms of first vsite */
                    vsnral      = NRAL(pindex[atom].ftype) - 1;
                    first_atoms = plist[pindex[atom].ftype].param[pindex[atom].parnr].a + 1;
                }
                else
                {
                    assert(vsnral != 0);
                    assert(first_atoms != NULL);
                    /* check if this vsite is constructed from the same atoms */
                    if (vsnral == NRAL(pindex[atom].ftype)-1)
                    {
                        for (m = 0; (m < vsnral) && !bKeep; m++)
                        {
                            const atom_id *atoms;

                            bPresent = FALSE;
                            atoms    = plist[pindex[atom].ftype].param[pindex[atom].parnr].a + 1;
                            for (n = 0; (n < vsnral) && !bPresent; n++)
                            {
                                if (atoms[m] == first_atoms[n])
                                {
                                    bPresent = TRUE;
                                }
                            }
                            if (!bPresent)
                            {
                                bKeep = TRUE;
                            }
                        }
                    }
                    else
                    {
                        bKeep = TRUE;
                    }
                }
            }
        }

        /* keep all angles with no virtual sites in them or
           with virtual sites with more than 3 constr. atoms */
        if (nvsite == 0 && vsnral > 3)
        {
            bKeep = TRUE;
        }

        /* check if all non-vsite atoms are used in construction: */
        bFirstTwo = TRUE;
        for (k = 0; (k < 3) && !bKeep; k++) /* for all atoms in the angle */
        {
            atom = ps->param[i].a[k];
            if (vsite_type[atom] == NOTSET && vsite_type[atom] != F_VSITEN)
            {
                bUsed = FALSE;
                for (m = 0; (m < vsnral) && !bUsed; m++)
                {
                    assert(first_atoms != NULL);

                    if (atom == first_atoms[m])
                    {
                        bUsed     = TRUE;
                        bFirstTwo = bFirstTwo && m < 2;
                    }
                }
                if (!bUsed)
                {
                    bKeep = TRUE;
                }
            }
        }

        if (!( bAll3FAD && bFirstTwo ) )
        {
            /* check if all constructing atoms are constrained together */
            for (m = 0; m < vsnral && !bKeep; m++) /* all constr. atoms */
            {
                at1      = first_atoms[m];
                at2      = first_atoms[(m+1) % vsnral];
                bPresent = FALSE;
                for (j = 0; j < at2vc[at1].nr; j++)
                {
                    if (at2vc[at1].aj[j] == at2)
                    {
                        bPresent = TRUE;
                    }
                }
                if (!bPresent)
                {
                    bKeep = TRUE;
                }
            }
        }

        if (bKeep)
        {
            /* now copy the angle to the new array */
            ps->param[kept_i] = ps->param[i];
            kept_i++;
        }
    }

    if (ps->nr != kept_i)
    {
        fprintf(stderr, "Removed   %4d %15ss with virtual sites, %5d left\n",
                ps->nr-kept_i, interaction_function[cftype].longname, kept_i);
    }
    ps->nr = kept_i;
}

static void clean_vsite_dihs(t_params *plist, t_pindex pindex[],
                             int cftype, int vsite_type[])
{
    int       i, kept_i;
    t_params *ps;

    ps = &(plist[cftype]);

    kept_i = 0;
    for (i = 0; (i < ps->nr); i++) /* for all dihedrals in the plist */
    {
        int            ftype, parnr, k, l, m, n, nvsite;
        int            vsnral      = 0;
        const atom_id *first_atoms = NULL;
        atom_id        atom;
        gmx_bool       bKeep, bUsed, bPresent;


        bKeep = FALSE;
        /* check if all virtual sites are constructed from the same atoms */
        nvsite = 0;
        for (k = 0; (k < 4) && !bKeep; k++) /* for all atoms in the dihedral */
        {
            atom = ps->param[i].a[k];
            if (vsite_type[atom] != NOTSET && vsite_type[atom] != F_VSITEN)
            {
                if (nvsite == 0)
                {
                    /* store construction atoms of first vsite */
                    vsnral      = NRAL(pindex[atom].ftype) - 1;
                    first_atoms = plist[pindex[atom].ftype].param[pindex[atom].parnr].a + 1;
                    if (debug)
                    {
                        fprintf(debug, "dih w. vsite: %d %d %d %d\n",
                                ps->param[i].AI+1, ps->param[i].AJ+1,
                                ps->param[i].AK+1, ps->param[i].AL+1);
                        fprintf(debug, "vsite %d from: %d %d %d\n",
                                atom+1, first_atoms[0]+1, first_atoms[1]+1, first_atoms[2]+1);
                    }
                }
                else
                {
                    assert(vsnral != 0);
                    assert(first_atoms != NULL);

                    /* check if this vsite is constructed from the same atoms */
                    if (vsnral == NRAL(pindex[atom].ftype)-1)
                    {
                        for (m = 0; (m < vsnral) && !bKeep; m++)
                        {
                            const atom_id *atoms;

                            bPresent = FALSE;
                            atoms    = plist[pindex[atom].ftype].param[pindex[atom].parnr].a + 1;
                            for (n = 0; (n < vsnral) && !bPresent; n++)
                            {
                                if (atoms[m] == first_atoms[n])
                                {
                                    bPresent = TRUE;
                                }
                            }
                            if (!bPresent)
                            {
                                bKeep = TRUE;
                            }
                        }
                    }
                }
                /* TODO clean_site_bonds and _angles do this increment
                   at the top of the loop. Refactor this for
                   consistency */
                nvsite++;
            }
        }

        /* keep all dihedrals with no virtual sites in them */
        if (nvsite == 0)
        {
            bKeep = TRUE;
        }

        /* check if all atoms in dihedral are either virtual sites, or used in
           construction of virtual sites. If so, keep it, if not throw away: */
        for (k = 0; (k < 4) && !bKeep; k++) /* for all atoms in the dihedral */
        {
            assert(vsnral != 0);
            assert(first_atoms != NULL);

            atom = ps->param[i].a[k];
            if (vsite_type[atom] == NOTSET && vsite_type[atom] != F_VSITEN)
            {
                /* vsnral will be set here, we don't get here with nvsite==0 */
                bUsed = FALSE;
                for (m = 0; (m < vsnral) && !bUsed; m++)
                {
                    if (atom == first_atoms[m])
                    {
                        bUsed = TRUE;
                    }
                }
                if (!bUsed)
                {
                    bKeep = TRUE;
                    if (debug)
                    {
                        fprintf(debug, "unused atom in dih: %d\n", atom+1);
                    }
                }
            }
        }

        if (bKeep)
        {
            ps->param[kept_i] = ps->param[i];
            kept_i++;
        }
    }

    if (ps->nr != kept_i)
    {
        fprintf(stderr, "Removed   %4d %15ss with virtual sites, %5d left\n",
                ps->nr-kept_i, interaction_function[cftype].longname, kept_i);
    }
    ps->nr = kept_i;
}

void clean_vsite_bondeds(t_params *plist, int natoms, gmx_bool bRmVSiteBds)
{
    int            i, k, nvsite, ftype, vsite, parnr;
    int           *vsite_type;
    t_pindex      *pindex;
    at2vsitecon_t *at2vc;

    pindex = 0; /* avoid warnings */
    /* make vsite_type array */
    snew(vsite_type, natoms);
    for (i = 0; i < natoms; i++)
    {
        vsite_type[i] = NOTSET;
    }
    nvsite = 0;
    for (ftype = 0; ftype < F_NRE; ftype++)
    {
        if (interaction_function[ftype].flags & IF_VSITE)
        {
            nvsite += plist[ftype].nr;
            i       = 0;
            while (i < plist[ftype].nr)
            {
                vsite = plist[ftype].param[i].AI;
                if (vsite_type[vsite] == NOTSET)
                {
                    vsite_type[vsite] = ftype;
                }
                else
                {
                    gmx_fatal(FARGS, "multiple vsite constructions for atom %d", vsite+1);
                }
                if (ftype == F_VSITEN)
                {
                    while (i < plist[ftype].nr && plist[ftype].param[i].AI == vsite)
                    {
                        i++;
                    }
                }
                else
                {
                    i++;
                }
            }
        }
    }

    /* the rest only if we have virtual sites: */
    if (nvsite)
    {
        fprintf(stderr, "Cleaning up constraints %swith virtual sites\n",
                bRmVSiteBds ? "and constant bonded interactions " : "");

        /* Make a reverse list to avoid ninteractions^2 operations */
        at2vc = make_at2vsitecon(natoms, plist);

        snew(pindex, natoms);
        for (ftype = 0; ftype < F_NRE; ftype++)
        {
            /* Here we skip VSITEN. In neary all practical use cases this
             * is not an issue, since VSITEN is intended for constructing
             * additional interaction sites, not for replacing normal atoms
             * with bonded interactions. Thus we do not expect constant
             * bonded interactions. If VSITEN does get used with constant
             * bonded interactions, these are not removed which only leads
             * to very minor extra computation and constant energy.
             * The only problematic case is onstraints between VSITEN
             * constructions with fixed distance (which is anyhow useless).
             * This will generate a fatal error in check_vsite_constraints.
             */
            if ((interaction_function[ftype].flags & IF_VSITE) &&
                ftype != F_VSITEN)
            {
                for (parnr = 0; (parnr < plist[ftype].nr); parnr++)
                {
                    k               = plist[ftype].param[parnr].AI;
                    pindex[k].ftype = ftype;
                    pindex[k].parnr = parnr;
                }
            }
        }

        if (debug)
        {
            for (i = 0; i < natoms; i++)
            {
                fprintf(debug, "atom %d vsite_type %s\n", i,
                        vsite_type[i] == NOTSET ? "NOTSET" :
                        interaction_function[vsite_type[i]].name);
            }
        }

        /* remove interactions that include virtual sites */
        for (ftype = 0; ftype < F_NRE; ftype++)
        {
            if ( ( ( interaction_function[ftype].flags & IF_BOND ) && bRmVSiteBds ) ||
                 ( interaction_function[ftype].flags & IF_CONSTRAINT ) )
            {
                if (interaction_function[ftype].flags & (IF_BTYPE | IF_CONSTRAINT) )
                {
                    clean_vsite_bonds (plist, pindex, ftype, vsite_type);
                }
                else if (interaction_function[ftype].flags & IF_ATYPE)
                {
                    clean_vsite_angles(plist, pindex, ftype, vsite_type, at2vc);
                }
                else if ( (ftype == F_PDIHS) || (ftype == F_IDIHS) )
                {
                    clean_vsite_dihs  (plist, pindex, ftype, vsite_type);
                }
            }
        }
        /* check that no remaining constraints include virtual sites */
        for (ftype = 0; ftype < F_NRE; ftype++)
        {
            if (interaction_function[ftype].flags & IF_CONSTRAINT)
            {
                check_vsite_constraints(plist, ftype, vsite_type);
            }
        }

        done_at2vsitecon(natoms, at2vc);
    }
    sfree(pindex);
    sfree(vsite_type);
}
