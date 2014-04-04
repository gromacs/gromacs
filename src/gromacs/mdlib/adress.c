/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2009 Christoph Junghans, Brad Lambeth.
 * Copyright (c) 2011,2012,2013,2014, by the GROMACS development team, led by
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

#include "adress.h"
#include "gromacs/math/utilities.h"
#include "pbc.h"
#include "types/simple.h"
#include "typedefs.h"
#include "vec.h"

#include "gmx_fatal.h"

real
adress_weight(rvec                 x,
              int                  adresstype,
              real                 adressr,
              real                 adressw,
              rvec      *          ref,
              t_pbc      *         pbc,
              t_forcerec *         fr )
{
    int  i;
    real l2 = adressr+adressw;
    real sqr_dl, dl;
    real tmp;
    rvec dx;

    sqr_dl = 0.0;

    if (pbc)
    {
        pbc_dx(pbc, (*ref), x, dx);
    }
    else
    {
        rvec_sub((*ref), x, dx);
    }

    switch (adresstype)
    {
        case eAdressOff:
            /* default to explicit simulation */
            return 1;
        case eAdressConst:
            /* constant value for weighting function = adressw */
            return fr->adress_const_wf;
        case eAdressXSplit:
            /* plane through center of ref, varies in x direction */
            sqr_dl         = dx[0]*dx[0];
            break;
        case eAdressSphere:
            /* point at center of ref, assuming cubic geometry */
            for (i = 0; i < 3; i++)
            {
                sqr_dl    += dx[i]*dx[i];
            }
            break;
        default:
            /* default to explicit simulation */
            return 1;
    }

    dl = sqrt(sqr_dl);

    /* molecule is coarse grained */
    if (dl > l2)
    {
        return 0;
    }
    /* molecule is explicit */
    else if (dl < adressr)
    {
        return 1;
    }
    /* hybrid region */
    else
    {
        tmp = cos((dl-adressr)*M_PI/2/adressw);
        return tmp*tmp;
    }
}

void
update_adress_weights_com(FILE gmx_unused    * fplog,
                          int                  cg0,
                          int                  cg1,
                          t_block *            cgs,
                          rvec                 x[],
                          t_forcerec *         fr,
                          t_mdatoms *          mdatoms,
                          t_pbc *              pbc)
{
    int            icg, k, k0, k1, d;
    real           nrcg, inv_ncg, mtot, inv_mtot;
    atom_id *      cgindex;
    rvec           ix;
    int            adresstype;
    real           adressr, adressw;
    rvec *         ref;
    real *         massT;
    real *         wf;


    int n_hyb, n_ex, n_cg;

    n_hyb = 0;
    n_cg  = 0;
    n_ex  = 0;

    adresstype         = fr->adress_type;
    adressr            = fr->adress_ex_width;
    adressw            = fr->adress_hy_width;
    massT              = mdatoms->massT;
    wf                 = mdatoms->wf;
    ref                = &(fr->adress_refs);


    /* Since this is center of mass AdResS, the vsite is not guaranteed
     * to be on the same node as the constructing atoms.  Therefore we
     * loop over the charge groups, calculate their center of mass,
     * then use this to calculate wf for each atom.  This wastes vsite
     * construction, but it's the only way to assure that the explicit
     * atoms have the same wf as their vsite. */

#ifdef DEBUG
    fprintf(fplog, "Calculating center of mass for charge groups %d to %d\n",
            cg0, cg1);
#endif
    cgindex = cgs->index;

    /* Compute the center of mass for all charge groups */
    for (icg = cg0; (icg < cg1); icg++)
    {
        k0      = cgindex[icg];
        k1      = cgindex[icg+1];
        nrcg    = k1-k0;
        if (nrcg == 1)
        {
            wf[k0] = adress_weight(x[k0], adresstype, adressr, adressw, ref, pbc, fr);
            if (wf[k0] == 0)
            {
                n_cg++;
            }
            else if (wf[k0] == 1)
            {
                n_ex++;
            }
            else
            {
                n_hyb++;
            }
        }
        else
        {
            mtot = 0.0;
            for (k = k0; (k < k1); k++)
            {
                mtot += massT[k];
            }
            if (mtot > 0.0)
            {
                inv_mtot = 1.0/mtot;

                clear_rvec(ix);
                for (k = k0; (k < k1); k++)
                {
                    for (d = 0; (d < DIM); d++)
                    {
                        ix[d] += x[k][d]*massT[k];
                    }
                }
                for (d = 0; (d < DIM); d++)
                {
                    ix[d] *= inv_mtot;
                }
            }
            /* Calculate the center of gravity if the charge group mtot=0 (only vsites) */
            else
            {
                inv_ncg = 1.0/nrcg;

                clear_rvec(ix);
                for (k = k0; (k < k1); k++)
                {
                    for (d = 0; (d < DIM); d++)
                    {
                        ix[d] += x[k][d];
                    }
                }
                for (d = 0; (d < DIM); d++)
                {
                    ix[d] *= inv_ncg;
                }
            }

            /* Set wf of all atoms in charge group equal to wf of com */
            wf[k0] = adress_weight(ix, adresstype, adressr, adressw, ref, pbc, fr);

            if (wf[k0] == 0)
            {
                n_cg++;
            }
            else if (wf[k0] == 1)
            {
                n_ex++;
            }
            else
            {
                n_hyb++;
            }

            for (k = (k0+1); (k < k1); k++)
            {
                wf[k] = wf[k0];
            }
        }
    }
}

void update_adress_weights_atom_per_atom(
        int                  cg0,
        int                  cg1,
        t_block *            cgs,
        rvec                 x[],
        t_forcerec *         fr,
        t_mdatoms *          mdatoms,
        t_pbc *              pbc)
{
    int            icg, k, k0, k1, d;
    real           nrcg, inv_ncg, mtot, inv_mtot;
    atom_id *      cgindex;
    rvec           ix;
    int            adresstype;
    real           adressr, adressw;
    rvec *         ref;
    real *         massT;
    real *         wf;


    int n_hyb, n_ex, n_cg;

    n_hyb = 0;
    n_cg  = 0;
    n_ex  = 0;

    adresstype         = fr->adress_type;
    adressr            = fr->adress_ex_width;
    adressw            = fr->adress_hy_width;
    massT              = mdatoms->massT;
    wf                 = mdatoms->wf;
    ref                = &(fr->adress_refs);

    cgindex = cgs->index;

    /* Weighting function is determined for each atom individually.
     * This is an approximation
     * as in the theory requires an interpolation based on the center of masses.
     * Should be used with caution */

    for (icg = cg0; (icg < cg1); icg++)
    {
        k0   = cgindex[icg];
        k1   = cgindex[icg + 1];
        nrcg = k1 - k0;

        for (k = (k0); (k < k1); k++)
        {
            wf[k] = adress_weight(x[k], adresstype, adressr, adressw, ref, pbc, fr);
            if (wf[k] == 0)
            {
                n_cg++;
            }
            else if (wf[k] == 1)
            {
                n_ex++;
            }
            else
            {
                n_hyb++;
            }
        }

    }
}

void
update_adress_weights_cog(t_iparams            ip[],
                          t_ilist              ilist[],
                          rvec                 x[],
                          t_forcerec *         fr,
                          t_mdatoms *          mdatoms,
                          t_pbc *              pbc)
{
    int            i, j, k, nr, nra, inc;
    int            ftype, adresstype;
    t_iatom        avsite, ai, aj, ak, al;
    t_iatom *      ia;
    real           adressr, adressw;
    rvec *         ref;
    real *         wf;
    int            n_hyb, n_ex, n_cg;

    adresstype         = fr->adress_type;
    adressr            = fr->adress_ex_width;
    adressw            = fr->adress_hy_width;
    wf                 = mdatoms->wf;
    ref                = &(fr->adress_refs);


    n_hyb = 0;
    n_cg  = 0;
    n_ex  = 0;


    /* Since this is center of geometry AdResS, we know the vsite
     * is in the same charge group node as the constructing atoms.
     * Loop over vsite types, calculate the weight of the vsite,
     * then assign that weight to the constructing atoms. */

    for (ftype = 0; (ftype < F_NRE); ftype++)
    {
        if (interaction_function[ftype].flags & IF_VSITE)
        {
            nra    = interaction_function[ftype].nratoms;
            nr     = ilist[ftype].nr;
            ia     = ilist[ftype].iatoms;

            for (i = 0; (i < nr); )
            {
                /* The vsite and first constructing atom */
                avsite     = ia[1];
                ai         = ia[2];
                wf[avsite] = adress_weight(x[avsite], adresstype, adressr, adressw, ref, pbc, fr);
                wf[ai]     = wf[avsite];

                if (wf[ai]  == 0)
                {
                    n_cg++;
                }
                else if (wf[ai]  == 1)
                {
                    n_ex++;
                }
                else
                {
                    n_hyb++;
                }

                /* Assign the vsite wf to rest of constructing atoms depending on type */
                inc = nra+1;
                switch (ftype)
                {
                    case F_VSITE2:
                        aj     = ia[3];
                        wf[aj] = wf[avsite];
                        break;
                    case F_VSITE3:
                        aj     = ia[3];
                        wf[aj] = wf[avsite];
                        ak     = ia[4];
                        wf[ak] = wf[avsite];
                        break;
                    case F_VSITE3FD:
                        aj     = ia[3];
                        wf[aj] = wf[avsite];
                        ak     = ia[4];
                        wf[ak] = wf[avsite];
                        break;
                    case F_VSITE3FAD:
                        aj     = ia[3];
                        wf[aj] = wf[avsite];
                        ak     = ia[4];
                        wf[ak] = wf[avsite];
                        break;
                    case F_VSITE3OUT:
                        aj     = ia[3];
                        wf[aj] = wf[avsite];
                        ak     = ia[4];
                        wf[ak] = wf[avsite];
                        break;
                    case F_VSITE4FD:
                        aj     = ia[3];
                        wf[aj] = wf[avsite];
                        ak     = ia[4];
                        wf[ak] = wf[avsite];
                        al     = ia[5];
                        wf[al] = wf[avsite];
                        break;
                    case F_VSITE4FDN:
                        aj     = ia[3];
                        wf[aj] = wf[avsite];
                        ak     = ia[4];
                        wf[ak] = wf[avsite];
                        al     = ia[5];
                        wf[al] = wf[avsite];
                        break;
                    case F_VSITEN:
                        inc    = 3*ip[ia[0]].vsiten.n;
                        for (j = 3; j < inc; j += 3)
                        {
                            ai     = ia[j+2];
                            wf[ai] = wf[avsite];
                        }
                        break;
                    default:
                        gmx_fatal(FARGS, "No such vsite type %d in %s, line %d",
                                  ftype, __FILE__, __LINE__);
                }

                /* Increment loop variables */
                i  += inc;
                ia += inc;
            }
        }
    }
}

void
update_adress_weights_atom(int                  cg0,
                           int                  cg1,
                           t_block *            cgs,
                           rvec                 x[],
                           t_forcerec *         fr,
                           t_mdatoms *          mdatoms,
                           t_pbc *              pbc)
{
    int            icg, k, k0, k1;
    atom_id *      cgindex;
    int            adresstype;
    real           adressr, adressw;
    rvec *         ref;
    real *         massT;
    real *         wf;

    adresstype         = fr->adress_type;
    adressr            = fr->adress_ex_width;
    adressw            = fr->adress_hy_width;
    massT              = mdatoms->massT;
    wf                 = mdatoms->wf;
    ref                = &(fr->adress_refs);
    cgindex            = cgs->index;

    /* Only use first atom in charge group.
     * We still can't be sure that the vsite and constructing
     * atoms are on the same processor, so we must calculate
     * in the same way as com adress. */

    for (icg = cg0; (icg < cg1); icg++)
    {
        k0      = cgindex[icg];
        k1      = cgindex[icg+1];
        wf[k0]  = adress_weight(x[k0], adresstype, adressr, adressw, ref, pbc, fr);

        /* Set wf of all atoms in charge group equal to wf of first atom in charge group*/
        for (k = (k0+1); (k < k1); k++)
        {
            wf[k] = wf[k0];
        }
    }
}

void
adress_thermo_force(int                  start,
                    int                  homenr,
                    t_block *            cgs,
                    rvec                 x[],
                    rvec                 f[],
                    t_forcerec *         fr,
                    t_mdatoms *          mdatoms,
                    t_pbc *              pbc)
{
    int              iatom, n0, nnn, nrcg, i;
    int              adresstype;
    real             adressw, adressr;
    atom_id *        cgindex;
    unsigned short * ptype;
    rvec *           ref;
    real *           wf;
    real             tabscale;
    real *           ATFtab;
    rvec             dr;
    real             w, wsq, wmin1, wmin1sq, wp, wt, rinv, sqr_dl, dl;
    real             eps, eps2, F, Geps, Heps2, Fp, dmu_dwp, dwp_dr, fscal;

    adresstype        = fr->adress_type;
    adressw           = fr->adress_hy_width;
    adressr           = fr->adress_ex_width;
    cgindex           = cgs->index;
    ptype             = mdatoms->ptype;
    ref               = &(fr->adress_refs);
    wf                = mdatoms->wf;

    for (iatom = start; (iatom < start+homenr); iatom++)
    {
        if (egp_coarsegrained(fr, mdatoms->cENER[iatom]))
        {
            if (ptype[iatom] == eptVSite)
            {
                w    = wf[iatom];
                /* is it hybrid or apply the thermodynamics force everywhere?*/
                if (mdatoms->tf_table_index[iatom] != NO_TF_TABLE)
                {
                    if (fr->n_adress_tf_grps > 0)
                    {
                        /* multi component tf is on, select the right table */
                        ATFtab   = fr->atf_tabs[mdatoms->tf_table_index[iatom]].data;
                        tabscale = fr->atf_tabs[mdatoms->tf_table_index[iatom]].scale;
                    }
                    else
                    {
                        /* just on component*/
                        ATFtab   = fr->atf_tabs[DEFAULT_TF_TABLE].data;
                        tabscale = fr->atf_tabs[DEFAULT_TF_TABLE].scale;
                    }

                    fscal            = 0;
                    if (pbc)
                    {
                        pbc_dx(pbc, (*ref), x[iatom], dr);
                    }
                    else
                    {
                        rvec_sub((*ref), x[iatom], dr);
                    }




                    /* calculate distace to adress center again */
                    sqr_dl = 0.0;
                    switch (adresstype)
                    {
                        case eAdressXSplit:
                            /* plane through center of ref, varies in x direction */
                            sqr_dl           = dr[0]*dr[0];
                            rinv             = gmx_invsqrt(dr[0]*dr[0]);
                            break;
                        case eAdressSphere:
                            /* point at center of ref, assuming cubic geometry */
                            for (i = 0; i < 3; i++)
                            {
                                sqr_dl    += dr[i]*dr[i];
                            }
                            rinv             = gmx_invsqrt(sqr_dl);
                            break;
                        default:
                            /* This case should not happen */
                            rinv = 0.0;
                    }

                    dl = sqrt(sqr_dl);
                    /* table origin is adress center */
                    wt               = dl*tabscale;
                    n0               = wt;
                    eps              = wt-n0;
                    eps2             = eps*eps;
                    nnn              = 4*n0;
                    F                = ATFtab[nnn+1];
                    Geps             = eps*ATFtab[nnn+2];
                    Heps2            = eps2*ATFtab[nnn+3];
                    Fp               = F+Geps+Heps2;
                    F                = (Fp+Geps+2.0*Heps2)*tabscale;

                    fscal            = F*rinv;

                    f[iatom][0]        += fscal*dr[0];
                    if (adresstype != eAdressXSplit)
                    {
                        f[iatom][1]    += fscal*dr[1];
                        f[iatom][2]    += fscal*dr[2];
                    }
                }
            }
        }
    }
}

gmx_bool egp_explicit(t_forcerec *   fr, int egp_nr)
{
    return fr->adress_group_explicit[egp_nr];
}

gmx_bool egp_coarsegrained(t_forcerec *   fr, int egp_nr)
{
    return !fr->adress_group_explicit[egp_nr];
}
