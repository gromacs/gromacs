/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2012,2014,2015, by the GROMACS development team, led by
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

#include "nb_generic.h"

#include <cmath>

#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/gmxlib/nonbonded/nb_kernel.h"
#include "gromacs/gmxlib/nonbonded/nonbonded.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/utility/fatalerror.h"

void
gmx_nb_generic_kernel(t_nblist *                nlist,
                      rvec *                    xx,
                      rvec *                    ff,
                      t_forcerec *              fr,
                      t_mdatoms *               mdatoms,
                      nb_kernel_data_t *        kernel_data,
                      t_nrnb *                  nrnb)
{
    int           ntype, table_nelements, ielec, ivdw;
    real          facel;
    int           n, ii, is3, ii3, k, nj0, nj1, jnr, j3, ggid, nnn, n0;
    real          shX, shY, shZ;
    real          fscal, felec, fvdw, velec, vvdw, tx, ty, tz;
    real          rinvsq;
    real          iq;
    real          qq, vctot;
    int           nti, nvdwparam;
    int           tj;
    real          rt, r, eps, eps2, Y, F, Geps, Heps2, VV, FF, Fp, fijD, fijR;
    real          rinvsix;
    real          vvdwtot;
    real          vvdw_rep, vvdw_disp;
    real          ix, iy, iz, fix, fiy, fiz;
    real          jx, jy, jz;
    real          dx, dy, dz, rsq, rinv;
    real          c6, c12, c6grid, cexp1, cexp2, br;
    real *        charge;
    real *        shiftvec;
    real *        vdwparam, *vdwgridparam;
    int *         type;
    real *        fshift;
    real *        velecgrp;
    real *        vvdwgrp;
    real          tabscale;
    real *        VFtab;
    real *        x;
    real *        f;
    int           ewitab;
    real          ewtabscale, eweps, ewrt, ewtabhalfspace;
    real *        ewtab;
    real          rcoulomb2, rvdw, rvdw2, sh_dispersion, sh_repulsion;
    real          rcutoff, rcutoff2;
    real          d, d2, sw, dsw, rinvcorr;
    real          elec_swV3, elec_swV4, elec_swV5, elec_swF2, elec_swF3, elec_swF4;
    real          vdw_swV3, vdw_swV4, vdw_swV5, vdw_swF2, vdw_swF3, vdw_swF4;
    real          ewclj, ewclj2, ewclj6, ewcljrsq, poly, exponent, sh_lj_ewald;
    gmx_bool      bExactElecCutoff, bExactVdwCutoff, bExactCutoff;
    gmx_bool      do_tab;

    x                   = xx[0];
    f                   = ff[0];
    ielec               = nlist->ielec;
    ivdw                = nlist->ivdw;

    fshift              = fr->fshift[0];
    velecgrp            = kernel_data->energygrp_elec;
    vvdwgrp             = kernel_data->energygrp_vdw;

    do_tab = (ielec == GMX_NBKERNEL_ELEC_CUBICSPLINETABLE ||
              ivdw == GMX_NBKERNEL_VDW_CUBICSPLINETABLE);
    if (do_tab)
    {
        tabscale         = kernel_data->table_elec_vdw->scale;
        VFtab            = kernel_data->table_elec_vdw->data;
    }
    else
    {
        tabscale        = 0;
        VFtab           = NULL;
    }
    ewtab               = fr->ic->tabq_coul_FDV0;
    ewtabscale          = fr->ic->tabq_scale;
    ewtabhalfspace      = 0.5/ewtabscale;

    rcoulomb2           = fr->rcoulomb*fr->rcoulomb;
    rvdw                = fr->rvdw;
    rvdw2               = rvdw*rvdw;
    sh_dispersion       = fr->ic->dispersion_shift.cpot;
    sh_repulsion        = fr->ic->repulsion_shift.cpot;
    sh_lj_ewald         = fr->ic->sh_lj_ewald;

    ewclj               = fr->ewaldcoeff_lj;
    ewclj2              = ewclj*ewclj;
    ewclj6              = ewclj2*ewclj2*ewclj2;

    if (fr->coulomb_modifier == eintmodPOTSWITCH)
    {
        d               = fr->rcoulomb-fr->rcoulomb_switch;
        elec_swV3       = -10.0/(d*d*d);
        elec_swV4       =  15.0/(d*d*d*d);
        elec_swV5       =  -6.0/(d*d*d*d*d);
        elec_swF2       = -30.0/(d*d*d);
        elec_swF3       =  60.0/(d*d*d*d);
        elec_swF4       = -30.0/(d*d*d*d*d);
    }
    else
    {
        /* Avoid warnings from stupid compilers (looking at you, Clang!) */
        elec_swV3 = elec_swV4 = elec_swV5 = elec_swF2 = elec_swF3 = elec_swF4 = 0.0;
    }
    if (fr->vdw_modifier == eintmodPOTSWITCH)
    {
        d               = fr->rvdw-fr->rvdw_switch;
        vdw_swV3        = -10.0/(d*d*d);
        vdw_swV4        =  15.0/(d*d*d*d);
        vdw_swV5        =  -6.0/(d*d*d*d*d);
        vdw_swF2        = -30.0/(d*d*d);
        vdw_swF3        =  60.0/(d*d*d*d);
        vdw_swF4        = -30.0/(d*d*d*d*d);
    }
    else
    {
        /* Avoid warnings from stupid compilers (looking at you, Clang!) */
        vdw_swV3 = vdw_swV4 = vdw_swV5 = vdw_swF2 = vdw_swF3 = vdw_swF4 = 0.0;
    }

    bExactElecCutoff    = (fr->coulomb_modifier != eintmodNONE) || fr->eeltype == eelRF_ZERO;
    bExactVdwCutoff     = (fr->vdw_modifier != eintmodNONE);
    bExactCutoff        = bExactElecCutoff && bExactVdwCutoff;

    if (bExactCutoff)
    {
        rcutoff  = ( fr->rcoulomb > fr->rvdw ) ? fr->rcoulomb : fr->rvdw;
        rcutoff2 = rcutoff*rcutoff;
    }
    else
    {
        /* Fix warnings for stupid compilers */
        rcutoff2 = 1e30;
    }

    /* avoid compiler warnings for cases that cannot happen */
    nnn                 = 0;
    eps                 = 0.0;
    eps2                = 0.0;

    /* 3 VdW parameters for Buckingham, otherwise 2 */
    nvdwparam           = (ivdw == GMX_NBKERNEL_VDW_BUCKINGHAM) ? 3 : 2;
    table_nelements     = 12;

    charge              = mdatoms->chargeA;
    type                = mdatoms->typeA;
    facel               = fr->epsfac;
    shiftvec            = fr->shift_vec[0];
    vdwparam            = fr->nbfp;
    ntype               = fr->ntype;
    vdwgridparam        = fr->ljpme_c6grid;

    for (n = 0; (n < nlist->nri); n++)
    {
        is3              = 3*nlist->shift[n];
        shX              = shiftvec[is3];
        shY              = shiftvec[is3+1];
        shZ              = shiftvec[is3+2];
        nj0              = nlist->jindex[n];
        nj1              = nlist->jindex[n+1];
        ii               = nlist->iinr[n];
        ii3              = 3*ii;
        ix               = shX + x[ii3+0];
        iy               = shY + x[ii3+1];
        iz               = shZ + x[ii3+2];
        iq               = facel*charge[ii];
        nti              = nvdwparam*ntype*type[ii];
        vctot            = 0;
        vvdwtot          = 0;
        fix              = 0;
        fiy              = 0;
        fiz              = 0;

        for (k = nj0; (k < nj1); k++)
        {
            jnr              = nlist->jjnr[k];
            j3               = 3*jnr;
            jx               = x[j3+0];
            jy               = x[j3+1];
            jz               = x[j3+2];
            dx               = ix - jx;
            dy               = iy - jy;
            dz               = iz - jz;
            rsq              = dx*dx+dy*dy+dz*dz;
            rinv             = gmx::invsqrt(rsq);
            rinvsq           = rinv*rinv;
            felec            = 0;
            fvdw             = 0;
            velec            = 0;
            vvdw             = 0;

            if (bExactCutoff && rsq >= rcutoff2)
            {
                continue;
            }

            if (ielec == GMX_NBKERNEL_ELEC_CUBICSPLINETABLE || ivdw == GMX_NBKERNEL_VDW_CUBICSPLINETABLE)
            {
                r                = rsq*rinv;
                rt               = r*tabscale;
                n0               = rt;
                eps              = rt-n0;
                eps2             = eps*eps;
                nnn              = table_nelements*n0;
            }

            /* Coulomb interaction. ielec==0 means no interaction */
            if (ielec != GMX_NBKERNEL_ELEC_NONE)
            {
                qq               = iq*charge[jnr];

                switch (ielec)
                {
                    case GMX_NBKERNEL_ELEC_NONE:
                        break;

                    case GMX_NBKERNEL_ELEC_COULOMB:
                        /* Vanilla cutoff coulomb */
                        velec            = qq*rinv;
                        felec            = velec*rinvsq;
                        /* The shift for the Coulomb potential is stored in
                         * the RF parameter c_rf, which is 0 without shift
                         */
                        velec           -= qq*fr->ic->c_rf;
                        break;

                    case GMX_NBKERNEL_ELEC_REACTIONFIELD:
                        /* Reaction-field */
                        velec            = qq*(rinv+fr->k_rf*rsq-fr->c_rf);
                        felec            = qq*(rinv*rinvsq-2.0*fr->k_rf);
                        break;

                    case GMX_NBKERNEL_ELEC_CUBICSPLINETABLE:
                        /* Tabulated coulomb */
                        Y                = VFtab[nnn];
                        F                = VFtab[nnn+1];
                        Geps             = eps*VFtab[nnn+2];
                        Heps2            = eps2*VFtab[nnn+3];
                        Fp               = F+Geps+Heps2;
                        VV               = Y+eps*Fp;
                        FF               = Fp+Geps+2.0*Heps2;
                        velec            = qq*VV;
                        felec            = -qq*FF*tabscale*rinv;
                        break;

                    case GMX_NBKERNEL_ELEC_GENERALIZEDBORN:
                        /* GB */
                        gmx_fatal(FARGS, "Death & horror! GB generic interaction not implemented.\n");
                        break;

                    case GMX_NBKERNEL_ELEC_EWALD:
                        ewrt             = rsq*rinv*ewtabscale;
                        ewitab           = ewrt;
                        eweps            = ewrt-ewitab;
                        ewitab           = 4*ewitab;
                        felec            = ewtab[ewitab]+eweps*ewtab[ewitab+1];
                        rinvcorr         = (fr->coulomb_modifier == eintmodPOTSHIFT) ? rinv-fr->ic->sh_ewald : rinv;
                        velec            = qq*(rinvcorr-(ewtab[ewitab+2]-ewtabhalfspace*eweps*(ewtab[ewitab]+felec)));
                        felec            = qq*rinv*(rinvsq-felec);
                        break;

                    default:
                        gmx_fatal(FARGS, "Death & horror! No generic coulomb interaction for ielec=%d.\n", ielec);
                        break;
                }
                if (fr->coulomb_modifier == eintmodPOTSWITCH)
                {
                    d                = rsq*rinv-fr->rcoulomb_switch;
                    d                = (d > 0.0) ? d : 0.0;
                    d2               = d*d;
                    sw               = 1.0+d2*d*(elec_swV3+d*(elec_swV4+d*elec_swV5));
                    dsw              = d2*(elec_swF2+d*(elec_swF3+d*elec_swF4));
                    /* Apply switch function. Note that felec=f/r since it will be multiplied
                     * by the i-j displacement vector. This means felec'=f'/r=-(v*sw)'/r=
                     * -(v'*sw+v*dsw)/r=-v'*sw/r-v*dsw/r=felec*sw-v*dsw/r
                     */
                    felec            = felec*sw - rinv*velec*dsw;
                    /* Once we have used velec to update felec we can modify velec too */
                    velec           *= sw;
                }
                if (bExactElecCutoff)
                {
                    felec            = (rsq < rcoulomb2) ? felec : 0.0;
                    velec            = (rsq < rcoulomb2) ? velec : 0.0;
                }
                vctot           += velec;
            } /* End of coulomb interactions */


            /* VdW interaction. ivdw==0 means no interaction */
            if (ivdw != GMX_NBKERNEL_VDW_NONE)
            {
                tj               = nti+nvdwparam*type[jnr];

                switch (ivdw)
                {
                    case GMX_NBKERNEL_VDW_NONE:
                        break;

                    case GMX_NBKERNEL_VDW_LENNARDJONES:
                        /* Vanilla Lennard-Jones cutoff */
                        c6               = vdwparam[tj];
                        c12              = vdwparam[tj+1];
                        rinvsix          = rinvsq*rinvsq*rinvsq;
                        vvdw_disp        = c6*rinvsix;
                        vvdw_rep         = c12*rinvsix*rinvsix;
                        fvdw             = (vvdw_rep-vvdw_disp)*rinvsq;
                        if (fr->vdw_modifier == eintmodPOTSHIFT)
                        {
                            vvdw             = (vvdw_rep + c12*sh_repulsion)/12.0 - (vvdw_disp + c6*sh_dispersion)/6.0;
                        }
                        else
                        {
                            vvdw             = vvdw_rep/12.0-vvdw_disp/6.0;
                        }
                        break;

                    case GMX_NBKERNEL_VDW_BUCKINGHAM:
                        /* Buckingham */
                        c6               = vdwparam[tj];
                        cexp1            = vdwparam[tj+1];
                        cexp2            = vdwparam[tj+2];

                        rinvsix          = rinvsq*rinvsq*rinvsq;
                        vvdw_disp        = c6*rinvsix;
                        br               = cexp2*rsq*rinv;
                        vvdw_rep         = cexp1*std::exp(-br);
                        fvdw             = (br*vvdw_rep-vvdw_disp)*rinvsq;
                        if (fr->vdw_modifier == eintmodPOTSHIFT)
                        {
                            vvdw             = (vvdw_rep-cexp1*std::exp(-cexp2*rvdw))-(vvdw_disp + c6*sh_dispersion)/6.0;
                        }
                        else
                        {
                            vvdw             = vvdw_rep-vvdw_disp/6.0;
                        }
                        break;

                    case GMX_NBKERNEL_VDW_CUBICSPLINETABLE:
                        /* Tabulated VdW */
                        c6               = vdwparam[tj];
                        c12              = vdwparam[tj+1];
                        Y                = VFtab[nnn+4];
                        F                = VFtab[nnn+5];
                        Geps             = eps*VFtab[nnn+6];
                        Heps2            = eps2*VFtab[nnn+7];
                        Fp               = F+Geps+Heps2;
                        VV               = Y+eps*Fp;
                        FF               = Fp+Geps+2.0*Heps2;
                        vvdw_disp        = c6*VV;
                        fijD             = c6*FF;
                        Y                = VFtab[nnn+8];
                        F                = VFtab[nnn+9];
                        Geps             = eps*VFtab[nnn+10];
                        Heps2            = eps2*VFtab[nnn+11];
                        Fp               = F+Geps+Heps2;
                        VV               = Y+eps*Fp;
                        FF               = Fp+Geps+2.0*Heps2;
                        vvdw_rep         = c12*VV;
                        fijR             = c12*FF;
                        fvdw             = -(fijD+fijR)*tabscale*rinv;
                        vvdw             = vvdw_disp + vvdw_rep;
                        break;


                    case GMX_NBKERNEL_VDW_LJEWALD:
                        /* LJ-PME */
                        rinvsix          = rinvsq*rinvsq*rinvsq;
                        ewcljrsq         = ewclj2*rsq;
                        exponent         = std::exp(-ewcljrsq);
                        poly             = exponent*(1.0 + ewcljrsq + ewcljrsq*ewcljrsq*0.5);
                        c6               = vdwparam[tj];
                        c12              = vdwparam[tj+1];
                        c6grid           = vdwgridparam[tj];
                        vvdw_disp        = (c6-c6grid*(1.0-poly))*rinvsix;
                        vvdw_rep         = c12*rinvsix*rinvsix;
                        fvdw             = (vvdw_rep - vvdw_disp - c6grid*(1.0/6.0)*exponent*ewclj6)*rinvsq;
                        if (fr->vdw_modifier == eintmodPOTSHIFT)
                        {
                            vvdw             = (vvdw_rep + c12*sh_repulsion)/12.0 - (vvdw_disp + c6*sh_dispersion - c6grid*sh_lj_ewald)/6.0;
                        }
                        else
                        {
                            vvdw             = vvdw_rep/12.0-vvdw_disp/6.0;
                        }
                        break;

                    default:
                        gmx_fatal(FARGS, "Death & horror! No generic VdW interaction for ivdw=%d.\n", ivdw);
                        break;
                }
                if (fr->vdw_modifier == eintmodPOTSWITCH)
                {
                    d                = rsq*rinv-fr->rvdw_switch;
                    d                = (d > 0.0) ? d : 0.0;
                    d2               = d*d;
                    sw               = 1.0+d2*d*(vdw_swV3+d*(vdw_swV4+d*vdw_swV5));
                    dsw              = d2*(vdw_swF2+d*(vdw_swF3+d*vdw_swF4));
                    /* See coulomb interaction for the force-switch formula */
                    fvdw             = fvdw*sw - rinv*vvdw*dsw;
                    vvdw            *= sw;
                }
                if (bExactVdwCutoff)
                {
                    fvdw             = (rsq < rvdw2) ? fvdw : 0.0;
                    vvdw             = (rsq < rvdw2) ? vvdw : 0.0;
                }
                vvdwtot         += vvdw;
            } /* end VdW interactions */

            fscal            = felec+fvdw;

            tx               = fscal*dx;
            ty               = fscal*dy;
            tz               = fscal*dz;
            fix              = fix + tx;
            fiy              = fiy + ty;
            fiz              = fiz + tz;
            f[j3+0]          = f[j3+0] - tx;
            f[j3+1]          = f[j3+1] - ty;
            f[j3+2]          = f[j3+2] - tz;
        }

        f[ii3+0]         = f[ii3+0] + fix;
        f[ii3+1]         = f[ii3+1] + fiy;
        f[ii3+2]         = f[ii3+2] + fiz;
        fshift[is3]      = fshift[is3]+fix;
        fshift[is3+1]    = fshift[is3+1]+fiy;
        fshift[is3+2]    = fshift[is3+2]+fiz;
        ggid             = nlist->gid[n];
        velecgrp[ggid]  += vctot;
        vvdwgrp[ggid]   += vvdwtot;
    }
    /* Estimate flops, average for generic kernel:
     * 12 flops per outer iteration
     * 50 flops per inner iteration
     */
    inc_nrnb(nrnb, eNR_NBKERNEL_GENERIC, nlist->nri*12 + nlist->jindex[n]*50);
}
