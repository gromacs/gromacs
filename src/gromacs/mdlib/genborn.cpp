/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2017, by the GROMACS development team, led by
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

#include "genborn.h"

#include <string.h>

#include <cmath>

#include <algorithm>

#include "gromacs/domdec/domdec.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/fileio/pdbio.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/genborn_allvsall.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/nblist.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/pbcutil/mshift.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxmpi.h"
#include "gromacs/utility/smalloc.h"


typedef struct {
    int  shift;
    int  naj;
    int *aj;
    int  aj_nalloc;
} gbtmpnbl_t;

typedef struct gbtmpnbls {
    int         nlist;
    gbtmpnbl_t *list;
    int         list_nalloc;
} t_gbtmpnbls;

/* This function is exactly the same as the one in listed-forces/bonded.cpp. The reason
 * it is copied here is that the bonded gb-interactions are evaluated
 * not in calc_bonds, but rather in calc_gb_forces
 */
static int pbc_rvec_sub(const t_pbc *pbc, const rvec xi, const rvec xj, rvec dx)
{
    if (pbc)
    {
        return pbc_dx_aiuc(pbc, xi, xj, dx);
    }
    else
    {
        rvec_sub(xi, xj, dx);
        return CENTRAL;
    }
}

static int init_gb_nblist(int natoms, t_nblist *nl)
{
    nl->maxnri      = natoms*4;
    nl->maxnrj      = 0;
    nl->nri         = 0;
    nl->nrj         = 0;
    nl->iinr        = nullptr;
    nl->gid         = nullptr;
    nl->shift       = nullptr;
    nl->jindex      = nullptr;
    nl->jjnr        = nullptr;
    /*nl->nltype      = nltype;*/

    srenew(nl->iinr,   nl->maxnri);
    srenew(nl->gid,    nl->maxnri);
    srenew(nl->shift,  nl->maxnri);
    srenew(nl->jindex, nl->maxnri+1);

    nl->jindex[0] = 0;

    return 0;
}


static int init_gb_still(const t_atomtypes *atype, t_idef *idef, t_atoms *atoms,
                         gmx_genborn_t *born, int natoms)
{

    int   i, j, m, ia, ib;
    real  r, ri, rj, ri2, rj2, r3, r4, ratio, term, h, doffset;

    real *vsol;
    real *gp;

    snew(vsol, natoms);
    snew(gp, natoms);
    snew(born->gpol_still_work, natoms+3);

    doffset = born->gb_doffset;

    for (i = 0; i < natoms; i++)
    {
        born->gpol_globalindex[i]              = born->vsolv_globalindex[i] =
                born->gb_radius_globalindex[i] = 0;
    }

    /* Compute atomic solvation volumes for Still method */
    for (i = 0; i < natoms; i++)
    {
        ri = atype->gb_radius[atoms->atom[i].type];
        born->gb_radius_globalindex[i] = ri;
        r3 = ri*ri*ri;
        born->vsolv_globalindex[i] = (4*M_PI/3)*r3;
    }

    for (j = 0; j < idef->il[F_GB12].nr; j += 3)
    {
        m  = idef->il[F_GB12].iatoms[j];
        ia = idef->il[F_GB12].iatoms[j+1];
        ib = idef->il[F_GB12].iatoms[j+2];

        r = 1.01*idef->iparams[m].gb.st;

        ri   = atype->gb_radius[atoms->atom[ia].type];
        rj   = atype->gb_radius[atoms->atom[ib].type];

        ri2  = ri*ri;
        rj2  = rj*rj;

        ratio  = (rj2-ri2-r*r)/(2*ri*r);
        h      = ri*(1+ratio);
        term   = (M_PI/3.0)*h*h*(3.0*ri-h);

        born->vsolv_globalindex[ia] -= term;

        ratio  = (ri2-rj2-r*r)/(2*rj*r);
        h      = rj*(1+ratio);
        term   = (M_PI/3.0)*h*h*(3.0*rj-h);

        born->vsolv_globalindex[ib] -= term;
    }

    /* Get the self-, 1-2 and 1-3 polarization energies for analytical Still
       method */
    /* Self */
    for (j = 0; j < natoms; j++)
    {
        if (born->use_globalindex[j] == 1)
        {
            born->gpol_globalindex[j] = -0.5*ONE_4PI_EPS0/
                (atype->gb_radius[atoms->atom[j].type]-doffset+STILL_P1);
        }
    }

    /* 1-2 */
    for (j = 0; j < idef->il[F_GB12].nr; j += 3)
    {
        m  = idef->il[F_GB12].iatoms[j];
        ia = idef->il[F_GB12].iatoms[j+1];
        ib = idef->il[F_GB12].iatoms[j+2];

        r = idef->iparams[m].gb.st;

        r4 = r*r*r*r;

        born->gpol_globalindex[ia] = born->gpol_globalindex[ia]+
            STILL_P2*born->vsolv_globalindex[ib]/r4;
        born->gpol_globalindex[ib] = born->gpol_globalindex[ib]+
            STILL_P2*born->vsolv_globalindex[ia]/r4;
    }

    /* 1-3 */
    for (j = 0; j < idef->il[F_GB13].nr; j += 3)
    {
        m  = idef->il[F_GB13].iatoms[j];
        ia = idef->il[F_GB13].iatoms[j+1];
        ib = idef->il[F_GB13].iatoms[j+2];

        r  = idef->iparams[m].gb.st;
        r4 = r*r*r*r;

        born->gpol_globalindex[ia] = born->gpol_globalindex[ia]+
            STILL_P3*born->vsolv_globalindex[ib]/r4;
        born->gpol_globalindex[ib] = born->gpol_globalindex[ib]+
            STILL_P3*born->vsolv_globalindex[ia]/r4;
    }

    sfree(vsol);
    sfree(gp);

    return 0;
}

/* Initialize all GB datastructs and compute polarization energies */
int init_gb(gmx_genborn_t **p_born,
            t_forcerec *fr, const t_inputrec *ir,
            const gmx_mtop_t *mtop, int gb_algorithm)
{
    int             i, jj, natoms;
    real            rai, sk, doffset;

    t_atoms         atoms;
    gmx_genborn_t  *born;
    gmx_localtop_t *localtop;

    natoms   = mtop->natoms;

    atoms    = gmx_mtop_global_atoms(mtop);
    localtop = gmx_mtop_generate_local_top(mtop, ir->efep != efepNO);

    snew(born, 1);
    *p_born = born;

    born->nr  = natoms;

    snew(born->drobc, natoms);
    snew(born->bRad,  natoms);

    /* Allocate memory for the global data arrays */
    snew(born->param_globalindex, natoms+3);
    snew(born->gpol_globalindex,  natoms+3);
    snew(born->vsolv_globalindex, natoms+3);
    snew(born->gb_radius_globalindex, natoms+3);
    snew(born->use_globalindex,    natoms+3);

    snew(fr->invsqrta, natoms);
    snew(fr->dvda,     natoms);

    fr->dadx              = nullptr;
    fr->dadx_rawptr       = nullptr;
    fr->nalloc_dadx       = 0;
    born->gpol_still_work = nullptr;
    born->gpol_hct_work   = nullptr;

    /* snew(born->asurf,natoms); */
    /* snew(born->dasurf,natoms); */

    /* Initialize the gb neighbourlist */
    snew(fr->gblist, 1);
    init_gb_nblist(natoms, fr->gblist);

    /* Do the Vsites exclusions (if any) */
    for (i = 0; i < natoms; i++)
    {
        jj = atoms.atom[i].type;
        if (mtop->atomtypes.gb_radius[atoms.atom[i].type] > 0)
        {
            born->use_globalindex[i] = 1;
        }
        else
        {
            born->use_globalindex[i] = 0;
        }

        /* If we have a Vsite, put vs_globalindex[i]=0 */
        if (C6 (fr->nbfp, fr->ntype, jj, jj) == 0 &&
            C12(fr->nbfp, fr->ntype, jj, jj) == 0 &&
            atoms.atom[i].q == 0)
        {
            born->use_globalindex[i] = 0;
        }
    }

    /* Copy algorithm parameters from inputrecord to local structure */
    born->obc_alpha          = ir->gb_obc_alpha;
    born->obc_beta           = ir->gb_obc_beta;
    born->obc_gamma          = ir->gb_obc_gamma;
    born->gb_doffset         = ir->gb_dielectric_offset;
    born->gb_epsilon_solvent = ir->gb_epsilon_solvent;
    born->epsilon_r          = ir->epsilon_r;

    doffset = born->gb_doffset;

    /* Set the surface tension */
    born->sa_surface_tension = ir->sa_surface_tension;

    /* If Still model, initialise the polarisation energies */
    if (gb_algorithm == egbSTILL)
    {
        init_gb_still(&(mtop->atomtypes), &(localtop->idef), &atoms,
                      born, natoms);
    }


    /* If HCT/OBC,  precalculate the sk*atype->S_hct factors */
    else if (gb_algorithm == egbHCT || gb_algorithm == egbOBC)
    {

        snew(born->gpol_hct_work, natoms+3);

        for (i = 0; i < natoms; i++)
        {
            if (born->use_globalindex[i] == 1)
            {
                rai = mtop->atomtypes.gb_radius[atoms.atom[i].type]-doffset;
                sk  = rai * mtop->atomtypes.S_hct[atoms.atom[i].type];
                born->param_globalindex[i]     = sk;
                born->gb_radius_globalindex[i] = rai;
            }
            else
            {
                born->param_globalindex[i]     = 0;
                born->gb_radius_globalindex[i] = 0;
            }
        }
    }

    /* Allocate memory for work arrays for temporary use */
    snew(born->work, natoms+4);
    snew(born->count, natoms);
    snew(born->nblist_work, natoms);

    /* Domain decomposition specific stuff */
    born->nalloc = 0;

    return 0;
}



static int
calc_gb_rad_still(t_commrec *cr, t_forcerec *fr, gmx_localtop_t *top,
                  rvec x[], t_nblist *nl,
                  gmx_genborn_t *born, t_mdatoms *md)
{
    int  i, k, n, nj0, nj1, ai, aj;
    int  shift;
    real shX, shY, shZ;
    real gpi, dr2, idr4, rvdw, ratio, ccf, theta, term, rai, raj;
    real ix1, iy1, iz1, jx1, jy1, jz1, dx11, dy11, dz11;
    real rinv, idr2, idr6, vaj, dccf, cosq, sinq, prod, gpi2;
    real factor;
    real vai, prod_ai, icf4, icf6;

    factor  = 0.5*ONE_4PI_EPS0;
    n       = 0;

    for (i = 0; i < born->nr; i++)
    {
        born->gpol_still_work[i] = 0;
    }

    for (i = 0; i < nl->nri; i++)
    {
        ai      = nl->iinr[i];

        nj0     = nl->jindex[i];
        nj1     = nl->jindex[i+1];

        /* Load shifts for this list */
        shift   = nl->shift[i];
        shX     = fr->shift_vec[shift][0];
        shY     = fr->shift_vec[shift][1];
        shZ     = fr->shift_vec[shift][2];

        gpi     = 0;

        rai     = top->atomtypes.gb_radius[md->typeA[ai]];
        vai     = born->vsolv[ai];
        prod_ai = STILL_P4*vai;

        /* Load atom i coordinates, add shift vectors */
        ix1     = shX + x[ai][0];
        iy1     = shY + x[ai][1];
        iz1     = shZ + x[ai][2];

        for (k = nj0; k < nj1 && nl->jjnr[k] >= 0; k++)
        {
            aj    = nl->jjnr[k];
            jx1   = x[aj][0];
            jy1   = x[aj][1];
            jz1   = x[aj][2];

            dx11  = ix1-jx1;
            dy11  = iy1-jy1;
            dz11  = iz1-jz1;

            dr2   = dx11*dx11+dy11*dy11+dz11*dz11;
            rinv  = gmx::invsqrt(dr2);
            idr2  = rinv*rinv;
            idr4  = idr2*idr2;
            idr6  = idr4*idr2;

            raj = top->atomtypes.gb_radius[md->typeA[aj]];

            rvdw  = rai + raj;

            ratio = dr2 / (rvdw * rvdw);
            vaj   = born->vsolv[aj];

            if (ratio > STILL_P5INV)
            {
                ccf  = 1.0;
                dccf = 0.0;
            }
            else
            {
                theta = ratio*STILL_PIP5;
                cosq  = cos(theta);
                term  = 0.5*(1.0-cosq);
                ccf   = term*term;
                sinq  = 1.0 - cosq*cosq;
                dccf  = 2.0*term*sinq*gmx::invsqrt(sinq)*theta;
            }

            prod                       = STILL_P4*vaj;
            icf4                       = ccf*idr4;
            icf6                       = (4*ccf-dccf)*idr6;
            born->gpol_still_work[aj] += prod_ai*icf4;
            gpi                        = gpi+prod*icf4;

            /* Save ai->aj and aj->ai chain rule terms */
            fr->dadx[n++]   = prod*icf6;
            fr->dadx[n++]   = prod_ai*icf6;
        }
        born->gpol_still_work[ai] += gpi;
    }

    /* Parallel summations */
    if (DOMAINDECOMP(cr))
    {
        dd_atom_sum_real(cr->dd, born->gpol_still_work);
    }

    /* Calculate the radii */
    for (i = 0; i < fr->natoms_force; i++) /* PELA born->nr */
    {
        if (born->use[i] != 0)
        {
            gpi             = born->gpol[i]+born->gpol_still_work[i];
            gpi2            = gpi * gpi;
            born->bRad[i]   = factor*gmx::invsqrt(gpi2);
            fr->invsqrta[i] = gmx::invsqrt(born->bRad[i]);
        }
    }

    /* Extra communication required for DD */
    if (DOMAINDECOMP(cr))
    {
        dd_atom_spread_real(cr->dd, born->bRad);
        dd_atom_spread_real(cr->dd, fr->invsqrta);
    }

    return 0;

}


static int
calc_gb_rad_hct(t_commrec *cr, t_forcerec *fr, gmx_localtop_t *top,
                rvec x[], t_nblist *nl,
                gmx_genborn_t *born, t_mdatoms *md)
{
    int   i, k, n, ai, aj, nj0, nj1;
    int   shift;
    real  shX, shY, shZ;
    real  rai, raj, dr2, dr, sk, sk_ai, sk2, sk2_ai, lij, uij, diff2, tmp, sum_ai;
    real  rad, min_rad, rinv, rai_inv;
    real  ix1, iy1, iz1, jx1, jy1, jz1, dx11, dy11, dz11;
    real  lij2, uij2, lij3, uij3, t1, t2, t3;
    real  lij_inv, dlij, sk2_rinv, prod, log_term;
    real  doffset, raj_inv, dadx_val;
    real *gb_radius;

    doffset   = born->gb_doffset;
    gb_radius = born->gb_radius;

    for (i = 0; i < born->nr; i++)
    {
        born->gpol_hct_work[i] = 0;
    }

    /* Keep the compiler happy */
    n    = 0;

    for (i = 0; i < nl->nri; i++)
    {
        ai     = nl->iinr[i];

        nj0    = nl->jindex[i];
        nj1    = nl->jindex[i+1];

        /* Load shifts for this list */
        shift   = nl->shift[i];
        shX     = fr->shift_vec[shift][0];
        shY     = fr->shift_vec[shift][1];
        shZ     = fr->shift_vec[shift][2];

        rai     = gb_radius[ai];
        rai_inv = 1.0/rai;

        sk_ai   = born->param[ai];
        sk2_ai  = sk_ai*sk_ai;

        /* Load atom i coordinates, add shift vectors */
        ix1     = shX + x[ai][0];
        iy1     = shY + x[ai][1];
        iz1     = shZ + x[ai][2];

        sum_ai  = 0;

        for (k = nj0; k < nj1 && nl->jjnr[k] >= 0; k++)
        {
            aj    = nl->jjnr[k];

            jx1   = x[aj][0];
            jy1   = x[aj][1];
            jz1   = x[aj][2];

            dx11  = ix1 - jx1;
            dy11  = iy1 - jy1;
            dz11  = iz1 - jz1;

            dr2   = dx11*dx11+dy11*dy11+dz11*dz11;
            rinv  = gmx::invsqrt(dr2);
            dr    = rinv*dr2;

            sk    = born->param[aj];
            raj   = gb_radius[aj];

            /* aj -> ai interaction */
            if (rai < dr+sk)
            {
                lij     = 1.0/(dr-sk);
                dlij    = 1.0;

                if (rai > dr-sk)
                {
                    lij  = rai_inv;
                    dlij = 0.0;
                }

                lij2     = lij*lij;
                lij3     = lij2*lij;

                uij      = 1.0/(dr+sk);
                uij2     = uij*uij;
                uij3     = uij2*uij;

                diff2    = uij2-lij2;

                lij_inv  = gmx::invsqrt(lij2);
                sk2      = sk*sk;
                sk2_rinv = sk2*rinv;
                prod     = 0.25*sk2_rinv;

                log_term = std::log(uij*lij_inv);

                tmp      = lij-uij + 0.25*dr*diff2 + (0.5*rinv)*log_term +
                    prod*(-diff2);

                if (rai < sk-dr)
                {
                    tmp = tmp + 2.0 * (rai_inv-lij);
                }

                t1 = 0.5*lij2 + prod*lij3 - 0.25*(lij*rinv+lij3*dr);
                t2 = -0.5*uij2 - 0.25*sk2_rinv*uij3 + 0.25*(uij*rinv+uij3*dr);
                t3 = 0.125*(1.0+sk2_rinv*rinv)*(-diff2)+0.25*log_term*rinv*rinv;

                dadx_val = (dlij*t1+t2+t3)*rinv; /* rb2 is moved to chainrule */
                /* fr->dadx[n++] = (dlij*t1+duij*t2+t3)*rinv; */
                /* rb2 is moved to chainrule    */

                sum_ai += 0.5*tmp;
            }
            else
            {
                dadx_val = 0.0;
            }
            fr->dadx[n++] = dadx_val;


            /* ai -> aj interaction */
            if (raj < dr + sk_ai)
            {
                lij     = 1.0/(dr-sk_ai);
                dlij    = 1.0;
                raj_inv = 1.0/raj;

                if (raj > dr-sk_ai)
                {
                    lij  = raj_inv;
                    dlij = 0.0;
                }

                lij2     = lij  * lij;
                lij3     = lij2 * lij;

                uij      = 1.0/(dr+sk_ai);
                uij2     = uij  * uij;
                uij3     = uij2 * uij;

                diff2    = uij2-lij2;

                lij_inv  = gmx::invsqrt(lij2);
                sk2      =  sk2_ai; /* sk2_ai = sk_ai * sk_ai in i loop above */
                sk2_rinv = sk2*rinv;
                prod     = 0.25 * sk2_rinv;

                /* log_term = table_log(uij*lij_inv,born->log_table,
                   LOG_TABLE_ACCURACY); */
                log_term = std::log(uij*lij_inv);

                tmp      = lij-uij + 0.25*dr*diff2 + (0.5*rinv)*log_term +
                    prod*(-diff2);

                if (raj < sk_ai-dr)
                {
                    tmp     = tmp + 2.0 * (raj_inv-lij);
                }

                /* duij = 1.0 */
                t1      = 0.5*lij2 + prod*lij3 - 0.25*(lij*rinv+lij3*dr);
                t2      = -0.5*uij2 - 0.25*sk2_rinv*uij3 + 0.25*(uij*rinv+uij3*dr);
                t3      = 0.125*(1.0+sk2_rinv*rinv)*(-diff2)+0.25*log_term*rinv*rinv;

                dadx_val = (dlij*t1+t2+t3)*rinv; /* rb2 is moved to chainrule    */
                /* fr->dadx[n++] = (dlij*t1+duij*t2+t3)*rinv; */ /* rb2 is moved to chainrule    */

                born->gpol_hct_work[aj] += 0.5*tmp;
            }
            else
            {
                dadx_val = 0.0;
            }
            fr->dadx[n++] = dadx_val;
        }

        born->gpol_hct_work[ai] += sum_ai;
    }

    /* Parallel summations */
    if (DOMAINDECOMP(cr))
    {
        dd_atom_sum_real(cr->dd, born->gpol_hct_work);
    }

    for (i = 0; i < fr->natoms_force; i++) /* PELA born->nr */
    {
        if (born->use[i] != 0)
        {
            rai     = top->atomtypes.gb_radius[md->typeA[i]]-doffset;
            sum_ai  = 1.0/rai - born->gpol_hct_work[i];
            min_rad = rai + doffset;
            rad     = 1.0/sum_ai;

            born->bRad[i]   = std::max(rad, min_rad);
            fr->invsqrta[i] = gmx::invsqrt(born->bRad[i]);
        }
    }

    /* Extra communication required for DD */
    if (DOMAINDECOMP(cr))
    {
        dd_atom_spread_real(cr->dd, born->bRad);
        dd_atom_spread_real(cr->dd, fr->invsqrta);
    }


    return 0;
}

static int
calc_gb_rad_obc(t_commrec *cr, t_forcerec *fr, gmx_localtop_t *top,
                rvec x[], t_nblist *nl, gmx_genborn_t *born, t_mdatoms *md)
{
    int   i, k, ai, aj, nj0, nj1, n;
    int   shift;
    real  shX, shY, shZ;
    real  rai, raj, dr2, dr, sk, sk2, lij, uij, diff2, tmp, sum_ai;
    real  sum_ai2, sum_ai3, tsum, tchain, rinv, rai_inv, lij_inv, rai_inv2;
    real  log_term, prod, sk2_rinv, sk_ai, sk2_ai;
    real  ix1, iy1, iz1, jx1, jy1, jz1, dx11, dy11, dz11;
    real  lij2, uij2, lij3, uij3, dlij, t1, t2, t3;
    real  doffset, raj_inv, dadx_val;
    real *gb_radius;

    /* Keep the compiler happy */
    n    = 0;

    doffset   = born->gb_doffset;
    gb_radius = born->gb_radius;

    for (i = 0; i < born->nr; i++)
    {
        born->gpol_hct_work[i] = 0;
    }

    for (i = 0; i < nl->nri; i++)
    {
        ai      = nl->iinr[i];

        nj0     = nl->jindex[i];
        nj1     = nl->jindex[i+1];

        /* Load shifts for this list */
        shift   = nl->shift[i];
        shX     = fr->shift_vec[shift][0];
        shY     = fr->shift_vec[shift][1];
        shZ     = fr->shift_vec[shift][2];

        rai      = gb_radius[ai];
        rai_inv  = 1.0/rai;

        sk_ai    = born->param[ai];
        sk2_ai   = sk_ai*sk_ai;

        /* Load atom i coordinates, add shift vectors */
        ix1      = shX + x[ai][0];
        iy1      = shY + x[ai][1];
        iz1      = shZ + x[ai][2];

        sum_ai   = 0;

        for (k = nj0; k < nj1 && nl->jjnr[k] >= 0; k++)
        {
            aj    = nl->jjnr[k];

            jx1   = x[aj][0];
            jy1   = x[aj][1];
            jz1   = x[aj][2];

            dx11  = ix1 - jx1;
            dy11  = iy1 - jy1;
            dz11  = iz1 - jz1;

            dr2   = dx11*dx11+dy11*dy11+dz11*dz11;
            rinv  = gmx::invsqrt(dr2);
            dr    = dr2*rinv;

            /* sk is precalculated in init_gb() */
            sk    = born->param[aj];
            raj   = gb_radius[aj];

            /* aj -> ai interaction */
            if (rai < dr+sk)
            {
                lij       = 1.0/(dr-sk);
                dlij      = 1.0;

                if (rai > dr-sk)
                {
                    lij  = rai_inv;
                    dlij = 0.0;
                }

                uij      = 1.0/(dr+sk);
                lij2     = lij  * lij;
                lij3     = lij2 * lij;
                uij2     = uij  * uij;
                uij3     = uij2 * uij;

                diff2    = uij2-lij2;

                lij_inv  = gmx::invsqrt(lij2);
                sk2      = sk*sk;
                sk2_rinv = sk2*rinv;
                prod     = 0.25*sk2_rinv;

                log_term = std::log(uij*lij_inv);

                tmp      = lij-uij + 0.25*dr*diff2 + (0.5*rinv)*log_term + prod*(-diff2);

                if (rai < sk-dr)
                {
                    tmp = tmp + 2.0 * (rai_inv-lij);
                }

                /* duij    = 1.0; */
                t1      = 0.5*lij2 + prod*lij3 - 0.25*(lij*rinv+lij3*dr);
                t2      = -0.5*uij2 - 0.25*sk2_rinv*uij3 + 0.25*(uij*rinv+uij3*dr);
                t3      = 0.125*(1.0+sk2_rinv*rinv)*(-diff2)+0.25*log_term*rinv*rinv;

                dadx_val = (dlij*t1+t2+t3)*rinv; /* rb2 is moved to chainrule    */

                sum_ai += 0.5*tmp;
            }
            else
            {
                dadx_val = 0.0;
            }
            fr->dadx[n++] = dadx_val;

            /* ai -> aj interaction */
            if (raj < dr + sk_ai)
            {
                lij     = 1.0/(dr-sk_ai);
                dlij    = 1.0;
                raj_inv = 1.0/raj;

                if (raj > dr-sk_ai)
                {
                    lij  = raj_inv;
                    dlij = 0.0;
                }

                lij2     = lij  * lij;
                lij3     = lij2 * lij;

                uij      = 1.0/(dr+sk_ai);
                uij2     = uij  * uij;
                uij3     = uij2 * uij;

                diff2    = uij2-lij2;

                lij_inv  = gmx::invsqrt(lij2);
                sk2      =  sk2_ai; /* sk2_ai = sk_ai * sk_ai in i loop above */
                sk2_rinv = sk2*rinv;
                prod     = 0.25 * sk2_rinv;

                /* log_term = table_log(uij*lij_inv,born->log_table,LOG_TABLE_ACCURACY); */
                log_term = std::log(uij*lij_inv);

                tmp      = lij-uij + 0.25*dr*diff2 + (0.5*rinv)*log_term + prod*(-diff2);

                if (raj < sk_ai-dr)
                {
                    tmp     = tmp + 2.0 * (raj_inv-lij);
                }

                t1      = 0.5*lij2 + prod*lij3 - 0.25*(lij*rinv+lij3*dr);
                t2      = -0.5*uij2 - 0.25*sk2_rinv*uij3 + 0.25*(uij*rinv+uij3*dr);
                t3      = 0.125*(1.0+sk2_rinv*rinv)*(-diff2)+0.25*log_term*rinv*rinv;

                dadx_val = (dlij*t1+t2+t3)*rinv; /* rb2 is moved to chainrule    */

                born->gpol_hct_work[aj] += 0.5*tmp;

            }
            else
            {
                dadx_val = 0.0;
            }
            fr->dadx[n++] = dadx_val;

        }
        born->gpol_hct_work[ai] += sum_ai;

    }

    /* Parallel summations */
    if (DOMAINDECOMP(cr))
    {
        dd_atom_sum_real(cr->dd, born->gpol_hct_work);
    }

    for (i = 0; i < fr->natoms_force; i++) /* PELA born->nr */
    {
        if (born->use[i] != 0)
        {
            rai        = top->atomtypes.gb_radius[md->typeA[i]];
            rai_inv2   = 1.0/rai;
            rai        = rai-doffset;
            rai_inv    = 1.0/rai;
            sum_ai     = rai * born->gpol_hct_work[i];
            sum_ai2    = sum_ai  * sum_ai;
            sum_ai3    = sum_ai2 * sum_ai;

            tsum          = tanh(born->obc_alpha*sum_ai-born->obc_beta*sum_ai2+born->obc_gamma*sum_ai3);
            born->bRad[i] = rai_inv - tsum*rai_inv2;
            born->bRad[i] = 1.0 / born->bRad[i];

            fr->invsqrta[i] = gmx::invsqrt(born->bRad[i]);

            tchain         = rai * (born->obc_alpha-2*born->obc_beta*sum_ai+3*born->obc_gamma*sum_ai2);
            born->drobc[i] = (1.0-tsum*tsum)*tchain*rai_inv2;
        }
    }

    /* Extra (local) communication required for DD */
    if (DOMAINDECOMP(cr))
    {
        dd_atom_spread_real(cr->dd, born->bRad);
        dd_atom_spread_real(cr->dd, fr->invsqrta);
        dd_atom_spread_real(cr->dd, born->drobc);
    }

    return 0;

}



int calc_gb_rad(t_commrec *cr, t_forcerec *fr, t_inputrec *ir, gmx_localtop_t *top,
                rvec x[], t_nblist *nl, gmx_genborn_t *born, t_mdatoms *md, t_nrnb     *nrnb)
{
    int   cnt;
    int   ndadx;

    if (fr->bAllvsAll && fr->dadx == nullptr)
    {
        /* We might need up to 8 atoms of padding before and after,
         * and another 4 units to guarantee SSE alignment.
         */
        fr->nalloc_dadx = 2*(md->homenr+12)*(md->nr/2+1+12);
        snew(fr->dadx_rawptr, fr->nalloc_dadx);
        fr->dadx = (real *) (((size_t) fr->dadx_rawptr + 16) & (~((size_t) 15)));
    }
    else
    {
        /* In the SSE-enabled gb-loops, when writing to dadx, we
         * always write 2*4 elements at a time, even in the case with only
         * 1-3 j particles, where we only really need to write 2*(1-3)
         * elements. This is because we want dadx to be aligned to a 16-
         * byte boundary, and being able to use _mm_store/load_ps
         */
        ndadx = 2 * (nl->nrj + 3*nl->nri);

        /* First, reallocate the dadx array, we need 3 extra for SSE */
        if (ndadx + 3 > fr->nalloc_dadx)
        {
            fr->nalloc_dadx = over_alloc_large(ndadx) + 3;
            srenew(fr->dadx_rawptr, fr->nalloc_dadx);
            fr->dadx = (real *) (((size_t) fr->dadx_rawptr + 16) & (~((size_t) 15)));
        }
    }

    if (fr->bAllvsAll)
    {
        cnt = md->homenr*(md->nr/2+1);

        if (ir->gb_algorithm == egbSTILL)
        {
            genborn_allvsall_calc_still_radii(fr, md, born, top, x[0], &fr->AllvsAll_workgb);
            /* 13 flops in outer loop, 47 flops in inner loop */
            inc_nrnb(nrnb, eNR_BORN_AVA_RADII_STILL, md->homenr*13+cnt*47);
        }
        else if (ir->gb_algorithm == egbHCT || ir->gb_algorithm == egbOBC)
        {
            genborn_allvsall_calc_hct_obc_radii(fr, md, born, ir->gb_algorithm, top, x[0], &fr->AllvsAll_workgb);
            /* 24 flops in outer loop, 183 in inner */
            inc_nrnb(nrnb, eNR_BORN_AVA_RADII_HCT_OBC, md->homenr*24+cnt*183);
        }
        else
        {
            gmx_fatal(FARGS, "Bad gb algorithm for all-vs-all interactions");
        }
        return 0;
    }

    /* Switch for determining which algorithm to use for Born radii calculation */
#if GMX_DOUBLE

    switch (ir->gb_algorithm)
    {
        case egbSTILL:
            calc_gb_rad_still(cr, fr, top, x, nl, born, md);
            break;
        case egbHCT:
            calc_gb_rad_hct(cr, fr, top, x, nl, born, md);
            break;
        case egbOBC:
            calc_gb_rad_obc(cr, fr, top, x, nl, born, md);
            break;

        default:
            gmx_fatal(FARGS, "Unknown double precision algorithm for Born radii calculation: %d", ir->gb_algorithm);
    }

#else

    switch (ir->gb_algorithm)
    {
        case egbSTILL:
            calc_gb_rad_still(cr, fr, top, x, nl, born, md);
            break;
        case egbHCT:
            calc_gb_rad_hct(cr, fr, top, x, nl, born, md);
            break;
        case egbOBC:
            calc_gb_rad_obc(cr, fr, top, x, nl, born, md);
            break;

        default:
            gmx_fatal(FARGS, "Unknown algorithm for Born radii calculation: %d", ir->gb_algorithm);
    }

#endif /* Double or single precision */

    if (fr->bAllvsAll == FALSE)
    {
        switch (ir->gb_algorithm)
        {
            case egbSTILL:
                /* 17 flops per outer loop iteration, 47 flops per inner loop */
                inc_nrnb(nrnb, eNR_BORN_RADII_STILL, nl->nri*17+nl->nrj*47);
                break;
            case egbHCT:
            case egbOBC:
                /* 61 (assuming 10 for tanh) flops for outer loop iteration, 183 flops per inner loop */
                inc_nrnb(nrnb, eNR_BORN_RADII_HCT_OBC, nl->nri*61+nl->nrj*183);
                break;

            default:
                break;
        }
    }

    return 0;
}



real gb_bonds_tab(rvec x[], rvec f[], rvec fshift[], real *charge, real *p_gbtabscale,
                  real *invsqrta, real *dvda, real *GBtab, t_idef *idef, real epsilon_r,
                  real gb_epsilon_solvent, real facel, const t_pbc *pbc, const t_graph *graph)
{
    int      i, j, n0, m, nnn, ai, aj;
    int      ki;

    real     isai, isaj;
    real     r, rsq11;
    real     rinv11, iq;
    real     isaprod, qq, gbscale, gbtabscale, Y, F, Geps, Heps2, Fp, VV, FF, rt, eps, eps2;
    real     vgb, fgb, fijC, dvdatmp, fscal;
    real     vctot;

    rvec     dx;
    ivec     dt;

    t_iatom *forceatoms;

    /* Scale the electrostatics by gb_epsilon_solvent */
    facel = facel * ((1.0/epsilon_r) - 1.0/gb_epsilon_solvent);

    gbtabscale = *p_gbtabscale;
    vctot      = 0.0;

    for (j = F_GB12; j <= F_GB14; j++)
    {
        forceatoms = idef->il[j].iatoms;

        for (i = 0; i < idef->il[j].nr; )
        {
            /* To avoid reading in the interaction type, we just increment i to pass over
             * the types in the forceatoms array, this saves some memory accesses
             */
            i++;
            ai            = forceatoms[i++];
            aj            = forceatoms[i++];

            ki            = pbc_rvec_sub(pbc, x[ai], x[aj], dx);
            rsq11         = iprod(dx, dx);

            isai          = invsqrta[ai];
            iq            = (-1)*facel*charge[ai];

            rinv11        = gmx::invsqrt(rsq11);
            isaj          = invsqrta[aj];
            isaprod       = isai*isaj;
            qq            = isaprod*iq*charge[aj];
            gbscale       = isaprod*gbtabscale;
            r             = rsq11*rinv11;
            rt            = r*gbscale;
            n0            = static_cast<int>(rt);
            eps           = rt-n0;
            eps2          = eps*eps;
            nnn           = 4*n0;
            Y             = GBtab[nnn];
            F             = GBtab[nnn+1];
            Geps          = eps*GBtab[nnn+2];
            Heps2         = eps2*GBtab[nnn+3];
            Fp            = F+Geps+Heps2;
            VV            = Y+eps*Fp;
            FF            = Fp+Geps+2.0*Heps2;
            vgb           = qq*VV;
            fijC          = qq*FF*gbscale;
            dvdatmp       = -(vgb+fijC*r)*0.5;
            dvda[aj]      = dvda[aj] + dvdatmp*isaj*isaj;
            dvda[ai]      = dvda[ai] + dvdatmp*isai*isai;
            vctot         = vctot + vgb;
            fgb           = -(fijC)*rinv11;

            if (graph)
            {
                ivec_sub(SHIFT_IVEC(graph, ai), SHIFT_IVEC(graph, aj), dt);
                ki = IVEC2IS(dt);
            }

            for (m = 0; (m < DIM); m++)             /*  15		*/
            {
                fscal               = fgb*dx[m];
                f[ai][m]           += fscal;
                f[aj][m]           -= fscal;
                fshift[ki][m]      += fscal;
                fshift[CENTRAL][m] -= fscal;
            }
        }
    }

    return vctot;
}

static real calc_gb_selfcorrections(t_commrec *cr, int natoms,
                                    real *charge, gmx_genborn_t *born, real *dvda, double facel)
{
    int  i, ai, at0, at1;
    real rai, e, derb, q, q2, fi, rai_inv, vtot;

    if (DOMAINDECOMP(cr))
    {
        at0 = 0;
        at1 = cr->dd->nat_home;
    }
    else
    {
        at0 = 0;
        at1 = natoms;

    }

    /* Scale the electrostatics by gb_epsilon_solvent */
    facel = facel * ((1.0/born->epsilon_r) - 1.0/born->gb_epsilon_solvent);

    vtot = 0.0;

    /* Apply self corrections */
    for (i = at0; i < at1; i++)
    {
        ai       = i;

        if (born->use[ai] == 1)
        {
            rai       = born->bRad[ai];
            rai_inv   = 1.0/rai;
            q         = charge[ai];
            q2        = q*q;
            fi        = facel*q2;
            e         = fi*rai_inv;
            derb      = 0.5*e*rai_inv*rai_inv;
            dvda[ai] += derb*rai;
            vtot     -= 0.5*e;
        }
    }

    return vtot;

}

static real calc_gb_nonpolar(t_commrec *cr, t_forcerec *fr, int natoms, gmx_genborn_t *born, gmx_localtop_t *top,
                             real *dvda, t_mdatoms *md)
{
    int  ai, i, at0, at1;
    real e, es, rai, term, probe, tmp, factor;
    real rbi_inv, rbi_inv2;

    if (DOMAINDECOMP(cr))
    {
        at0 = 0;
        at1 = cr->dd->nat_home;
    }
    else
    {
        at0 = 0;
        at1 = natoms;
    }

    /* factor is the surface tension */
    factor = born->sa_surface_tension;

    es    = 0;
    probe = 0.14;
    term  = M_PI*4;

    for (i = at0; i < at1; i++)
    {
        ai        = i;

        if (born->use[ai] == 1)
        {
            rai       = top->atomtypes.gb_radius[md->typeA[ai]];
            rbi_inv   = fr->invsqrta[ai];
            rbi_inv2  = rbi_inv * rbi_inv;
            tmp       = (rai*rbi_inv2)*(rai*rbi_inv2);
            tmp       = tmp*tmp*tmp;
            e         = factor*term*(rai+probe)*(rai+probe)*tmp;
            dvda[ai]  = dvda[ai] - 6*e*rbi_inv2;
            es        = es + e;
        }
    }

    return es;
}



static real calc_gb_chainrule(int natoms, t_nblist *nl, real *dadx, real *dvda, rvec x[], rvec t[], rvec fshift[],
                              rvec shift_vec[], int gb_algorithm, gmx_genborn_t *born)
{
    int          i, k, n, ai, aj, nj0, nj1, n0, n1;
    int          shift;
    real         shX, shY, shZ;
    real         fgb, rbi, fix1, fiy1, fiz1;
    real         ix1, iy1, iz1, jx1, jy1, jz1, dx11, dy11, dz11;
    real         tx, ty, tz, rbai, rbaj, fgb_ai;
    real        *rb;

    n  = 0;
    rb = born->work;

    n0 = 0;
    n1 = natoms;

    if (gb_algorithm == egbSTILL)
    {
        for (i = n0; i < n1; i++)
        {
            rbi   = born->bRad[i];
            rb[i] = (2 * rbi * rbi * dvda[i])/ONE_4PI_EPS0;
        }
    }
    else if (gb_algorithm == egbHCT)
    {
        for (i = n0; i < n1; i++)
        {
            rbi   = born->bRad[i];
            rb[i] = rbi * rbi * dvda[i];
        }
    }
    else if (gb_algorithm == egbOBC)
    {
        for (i = n0; i < n1; i++)
        {
            rbi   = born->bRad[i];
            rb[i] = rbi * rbi * born->drobc[i] * dvda[i];
        }
    }

    for (i = 0; i < nl->nri; i++)
    {
        ai   = nl->iinr[i];

        nj0  = nl->jindex[i];
        nj1  = nl->jindex[i+1];

        /* Load shifts for this list */
        shift   = nl->shift[i];
        shX     = shift_vec[shift][0];
        shY     = shift_vec[shift][1];
        shZ     = shift_vec[shift][2];

        /* Load atom i coordinates, add shift vectors */
        ix1  = shX + x[ai][0];
        iy1  = shY + x[ai][1];
        iz1  = shZ + x[ai][2];

        fix1 = 0;
        fiy1 = 0;
        fiz1 = 0;

        rbai = rb[ai];

        for (k = nj0; k < nj1 && nl->jjnr[k] >= 0; k++)
        {
            aj = nl->jjnr[k];

            jx1     = x[aj][0];
            jy1     = x[aj][1];
            jz1     = x[aj][2];

            dx11    = ix1 - jx1;
            dy11    = iy1 - jy1;
            dz11    = iz1 - jz1;

            rbaj    = rb[aj];

            fgb     = rbai*dadx[n++];
            fgb_ai  = rbaj*dadx[n++];

            /* Total force between ai and aj is the sum of ai->aj and aj->ai */
            fgb     = fgb + fgb_ai;

            tx      = fgb * dx11;
            ty      = fgb * dy11;
            tz      = fgb * dz11;

            fix1    = fix1 + tx;
            fiy1    = fiy1 + ty;
            fiz1    = fiz1 + tz;

            /* Update force on atom aj */
            t[aj][0] = t[aj][0] - tx;
            t[aj][1] = t[aj][1] - ty;
            t[aj][2] = t[aj][2] - tz;
        }

        /* Update force and shift forces on atom ai */
        t[ai][0] = t[ai][0] + fix1;
        t[ai][1] = t[ai][1] + fiy1;
        t[ai][2] = t[ai][2] + fiz1;

        fshift[shift][0] = fshift[shift][0] + fix1;
        fshift[shift][1] = fshift[shift][1] + fiy1;
        fshift[shift][2] = fshift[shift][2] + fiz1;

    }

    return 0;
}


void
calc_gb_forces(t_commrec *cr, t_mdatoms *md, gmx_genborn_t *born, gmx_localtop_t *top,
               rvec x[], rvec f[], t_forcerec *fr, t_idef *idef, int gb_algorithm, int sa_algorithm, t_nrnb *nrnb,
               const t_pbc *pbc, const t_graph *graph, gmx_enerdata_t *enerd)
{
    int  cnt;

    /* PBC or not? */
    const t_pbc *pbc_null;

    if (fr->bMolPBC)
    {
        pbc_null = pbc;
    }
    else
    {
        pbc_null = nullptr;
    }

    if (sa_algorithm == esaAPPROX)
    {
        /* Do a simple ACE type approximation for the non-polar solvation */
        enerd->term[F_NPSOLVATION] += calc_gb_nonpolar(cr, fr, born->nr, born, top, fr->dvda, md);
    }

    /* Calculate the bonded GB-interactions using either table or analytical formula */
    enerd->term[F_GBPOL]       += gb_bonds_tab(x, f, fr->fshift, md->chargeA, &(fr->gbtabscale),
                                               fr->invsqrta, fr->dvda, fr->gbtab->data, idef, born->epsilon_r, born->gb_epsilon_solvent, fr->ic->epsfac, pbc_null, graph);

    /* Calculate self corrections to the GB energies - currently only A state used! (FIXME) */
    enerd->term[F_GBPOL]       += calc_gb_selfcorrections(cr, born->nr, md->chargeA, born, fr->dvda, fr->ic->epsfac);

    /* If parallel, sum the derivative of the potential w.r.t the born radii */
    if (DOMAINDECOMP(cr))
    {
        dd_atom_sum_real(cr->dd, fr->dvda);
        dd_atom_spread_real(cr->dd, fr->dvda);
    }

    if (fr->bAllvsAll)
    {
        genborn_allvsall_calc_chainrule(fr, md, born, x[0], f[0], gb_algorithm, fr->AllvsAll_workgb);
        cnt = md->homenr*(md->nr/2+1);
        /* 9 flops for outer loop, 15 for inner */
        inc_nrnb(nrnb, eNR_BORN_AVA_CHAINRULE, md->homenr*9+cnt*15);
        return;
    }

    calc_gb_chainrule(fr->natoms_force, fr->gblist, fr->dadx, fr->dvda,
                      x, f, fr->fshift, fr->shift_vec, gb_algorithm, born);

    if (!fr->bAllvsAll)
    {
        /* 9 flops for outer loop, 15 for inner */
        inc_nrnb(nrnb, eNR_BORN_CHAINRULE, fr->gblist->nri*9+fr->gblist->nrj*15);
    }
}

static void add_j_to_gblist(gbtmpnbl_t *list, int aj)
{
    if (list->naj >= list->aj_nalloc)
    {
        list->aj_nalloc = over_alloc_large(list->naj+1);
        srenew(list->aj, list->aj_nalloc);
    }

    list->aj[list->naj++] = aj;
}

static gbtmpnbl_t *find_gbtmplist(struct gbtmpnbls *lists, int shift)
{
    int ind, i;

    /* Search the list with the same shift, if there is one */
    ind = 0;
    while (ind < lists->nlist && shift != lists->list[ind].shift)
    {
        ind++;
    }
    if (ind == lists->nlist)
    {
        if (lists->nlist == lists->list_nalloc)
        {
            lists->list_nalloc++;
            srenew(lists->list, lists->list_nalloc);
            for (i = lists->nlist; i < lists->list_nalloc; i++)
            {
                lists->list[i].aj        = nullptr;
                lists->list[i].aj_nalloc = 0;
            }

        }

        lists->list[lists->nlist].shift = shift;
        lists->list[lists->nlist].naj   = 0;
        lists->nlist++;
    }

    return &lists->list[ind];
}

static void add_bondeds_to_gblist(t_ilist *il,
                                  gmx_bool bMolPBC, t_pbc *pbc, t_graph *g, rvec *x,
                                  struct gbtmpnbls *nls)
{
    int         ind, j, ai, aj, found;
    rvec        dx;
    ivec        dt;
    gbtmpnbl_t *list;

    for (ind = 0; ind < il->nr; ind += 3)
    {
        ai = il->iatoms[ind+1];
        aj = il->iatoms[ind+2];

        int shift = CENTRAL;
        if (g != nullptr)
        {
            rvec_sub(x[ai], x[aj], dx);
            ivec_sub(SHIFT_IVEC(g, ai), SHIFT_IVEC(g, aj), dt);
            shift = IVEC2IS(dt);
        }
        else if (bMolPBC)
        {
            shift = pbc_dx_aiuc(pbc, x[ai], x[aj], dx);
        }

        /* Find the list for this shift or create one */
        list = find_gbtmplist(&nls[ai], shift);

        found = 0;

        /* So that we do not add the same bond twice.
         * This happens with some constraints between 1-3 atoms
         * that are in the bond-list but should not be in the GB nb-list */
        for (j = 0; j < list->naj; j++)
        {
            if (list->aj[j] == aj)
            {
                found = 1;
            }
        }

        if (found == 0)
        {
            if (ai == aj)
            {
                gmx_incons("ai == aj");
            }

            add_j_to_gblist(list, aj);
        }
    }
}


int make_gb_nblist(t_commrec *cr, int gb_algorithm,
                   rvec x[], matrix box,
                   t_forcerec *fr, t_idef *idef, t_graph *graph, gmx_genborn_t *born)
{
    int               i, j, k, n, nj0, nj1, ai, shift, s;
    t_nblist         *nblist;
    t_pbc             pbc;

    struct gbtmpnbls *nls;
    gbtmpnbl_t       *list = nullptr;

    set_pbc(&pbc, fr->ePBC, box);
    nls   = born->nblist_work;

    for (i = 0; i < born->nr; i++)
    {
        nls[i].nlist = 0;
    }

    if (fr->bMolPBC)
    {
        set_pbc_dd(&pbc, fr->ePBC, cr->dd->nc, TRUE, box);
    }

    switch (gb_algorithm)
    {
        case egbHCT:
        case egbOBC:
            /* Loop over 1-2, 1-3 and 1-4 interactions */
            for (j = F_GB12; j <= F_GB14; j++)
            {
                add_bondeds_to_gblist(&idef->il[j], fr->bMolPBC, &pbc, graph, x, nls);
            }
            break;
        case egbSTILL:
            /* Loop over 1-4 interactions */
            add_bondeds_to_gblist(&idef->il[F_GB14], fr->bMolPBC, &pbc, graph, x, nls);
            break;
        default:
            gmx_incons("Unknown GB algorithm");
    }

    /* Loop over the VDWQQ and VDW nblists to set up the nonbonded part of the GB list */
    for (n = 0; (n < fr->nnblists); n++)
    {
        for (i = 0; (i < eNL_NR); i++)
        {
            nblist = &(fr->nblists[n].nlist_sr[i]);

            if (nblist->nri > 0 && (i == eNL_VDWQQ || i == eNL_QQ))
            {
                for (j = 0; j < nblist->nri; j++)
                {
                    ai    = nblist->iinr[j];
                    shift = nblist->shift[j];

                    /* Find the list for this shift or create one */
                    list = find_gbtmplist(&nls[ai], shift);

                    nj0 = nblist->jindex[j];
                    nj1 = nblist->jindex[j+1];

                    /* Add all the j-atoms in the non-bonded list to the GB list */
                    for (k = nj0; k < nj1; k++)
                    {
                        add_j_to_gblist(list, nblist->jjnr[k]);
                    }
                }
            }
        }
    }

    /* Zero out some counters */
    fr->gblist->nri = 0;
    fr->gblist->nrj = 0;

    fr->gblist->jindex[0] = fr->gblist->nri;

    for (i = 0; i < fr->natoms_force; i++)
    {
        for (s = 0; s < nls[i].nlist; s++)
        {
            list = &nls[i].list[s];

            /* Only add those atoms that actually have neighbours */
            if (born->use[i] != 0)
            {
                fr->gblist->iinr[fr->gblist->nri]  = i;
                fr->gblist->shift[fr->gblist->nri] = list->shift;
                fr->gblist->nri++;

                for (k = 0; k < list->naj; k++)
                {
                    /* Memory allocation for jjnr */
                    if (fr->gblist->nrj >= fr->gblist->maxnrj)
                    {
                        fr->gblist->maxnrj += over_alloc_large(fr->gblist->maxnrj);

                        if (debug)
                        {
                            fprintf(debug, "Increasing GB neighbourlist j size to %d\n", fr->gblist->maxnrj);
                        }

                        srenew(fr->gblist->jjnr, fr->gblist->maxnrj);
                    }

                    /* Put in list */
                    if (i == list->aj[k])
                    {
                        gmx_incons("i == list->aj[k]");
                    }
                    fr->gblist->jjnr[fr->gblist->nrj++] = list->aj[k];
                }

                fr->gblist->jindex[fr->gblist->nri] = fr->gblist->nrj;
            }
        }
    }

    return 0;
}

void make_local_gb(const t_commrec *cr, gmx_genborn_t *born, int gb_algorithm)
{
    int           i, at0, at1;
    gmx_domdec_t *dd = nullptr;

    if (DOMAINDECOMP(cr))
    {
        dd  = cr->dd;
        at0 = 0;
        at1 = dd->nat_tot;
    }
    else
    {
        /* Single node, just copy pointers and return */
        if (gb_algorithm == egbSTILL)
        {
            born->gpol      = born->gpol_globalindex;
            born->vsolv     = born->vsolv_globalindex;
            born->gb_radius = born->gb_radius_globalindex;
        }
        else
        {
            born->param     = born->param_globalindex;
            born->gb_radius = born->gb_radius_globalindex;
        }

        born->use = born->use_globalindex;

        return;
    }

    /* Reallocation of local arrays if necessary */
    /* fr->natoms_force is equal to dd->nat_tot */
    if (DOMAINDECOMP(cr) && dd->nat_tot > born->nalloc)
    {
        int nalloc;

        nalloc = dd->nat_tot;

        /* Arrays specific to different gb algorithms */
        if (gb_algorithm == egbSTILL)
        {
            srenew(born->gpol,  nalloc+3);
            srenew(born->vsolv, nalloc+3);
            srenew(born->gb_radius, nalloc+3);
            for (i = born->nalloc; (i < nalloc+3); i++)
            {
                born->gpol[i]      = 0;
                born->vsolv[i]     = 0;
                born->gb_radius[i] = 0;
            }
        }
        else
        {
            srenew(born->param, nalloc+3);
            srenew(born->gb_radius, nalloc+3);
            for (i = born->nalloc; (i < nalloc+3); i++)
            {
                born->param[i]     = 0;
                born->gb_radius[i] = 0;
            }
        }

        /* All gb-algorithms use the array for vsites exclusions */
        srenew(born->use,    nalloc+3);
        for (i = born->nalloc; (i < nalloc+3); i++)
        {
            born->use[i] = 0;
        }

        born->nalloc = nalloc;
    }

    /* With dd, copy algorithm specific arrays */
    if (gb_algorithm == egbSTILL)
    {
        for (i = at0; i < at1; i++)
        {
            born->gpol[i]      = born->gpol_globalindex[dd->gatindex[i]];
            born->vsolv[i]     = born->vsolv_globalindex[dd->gatindex[i]];
            born->gb_radius[i] = born->gb_radius_globalindex[dd->gatindex[i]];
            born->use[i]       = born->use_globalindex[dd->gatindex[i]];
        }
    }
    else
    {
        for (i = at0; i < at1; i++)
        {
            born->param[i]     = born->param_globalindex[dd->gatindex[i]];
            born->gb_radius[i] = born->gb_radius_globalindex[dd->gatindex[i]];
            born->use[i]       = born->use_globalindex[dd->gatindex[i]];
        }
    }
}
