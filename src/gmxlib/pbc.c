/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * GROningen Mixture of Alchemy and Childrens' Stories
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>
#include "sysstuff.h"
#include "typedefs.h"
#include "vec.h"
#include "maths.h"
#include "main.h"
#include "pbc.h"
#include "smalloc.h"
#include "txtdump.h"
#include "gmx_fatal.h"
#include "names.h"

/* Skip 0 so we have more chance of detecting if we forgot to call set_pbc. */
enum { epbcdxRECTANGULAR=1, epbcdxTRICLINIC,
    epbcdx2D_RECT,       epbcdx2D_TRIC,
    epbcdx1D_RECT,       epbcdx1D_TRIC,
    epbcdxSCREW_RECT,    epbcdxSCREW_TRIC,
    epbcdxNOPBC,         epbcdxUNSUPPORTED };

/* Margin factor for error message and correction if the box is too skewed */
#define BOX_MARGIN         1.0010
#define BOX_MARGIN_CORRECT 1.0005

int ePBC2npbcdim(int ePBC)
{
    int npbcdim=0;

    switch(ePBC) {
    case epbcXYZ:   npbcdim = 3; break;
    case epbcXY:    npbcdim = 2; break;
    case epbcSCREW: npbcdim = 3; break;
    case epbcNONE:  npbcdim = 0; break;
    default: gmx_fatal(FARGS,"Unknown ePBC=%d in ePBC2npbcdim",ePBC);
    }

    return npbcdim;
}

int inputrec2nboundeddim(t_inputrec *ir)
{
    if (ir->nwall == 2 && ir->ePBC == epbcXY)
    {
        return 3;
    }
    else
    {
        return ePBC2npbcdim(ir->ePBC);
    }
}

void dump_pbc(FILE *fp,t_pbc *pbc) 
{
    rvec sum_box;

    fprintf(fp,"ePBCDX = %d\n",pbc->ePBCDX);
    pr_rvecs(fp,0,"box",pbc->box,DIM);
    pr_rvecs(fp,0,"fbox_diag",&pbc->fbox_diag,1);
    pr_rvecs(fp,0,"hbox_diag",&pbc->hbox_diag,1);
    pr_rvecs(fp,0,"mhbox_diag",&pbc->mhbox_diag,1);
    rvec_add(pbc->hbox_diag,pbc->mhbox_diag,sum_box);
    pr_rvecs(fp,0,"sum of the above two",&sum_box,1);
    fprintf(fp,"max_cutoff2 = %g\n",pbc->max_cutoff2);
    fprintf(fp,"bLimitDistance = %s\n",BOOL(pbc->bLimitDistance));
    fprintf(fp,"limit_distance2 = %g\n",pbc->limit_distance2);
    fprintf(fp,"ntric_vec = %d\n",pbc->ntric_vec);
    if (pbc->ntric_vec > 0) {
        pr_ivecs(fp,0,"tric_shift",pbc->tric_shift,pbc->ntric_vec,FALSE);
        pr_rvecs(fp,0,"tric_vec",pbc->tric_vec,pbc->ntric_vec);
    }
}

const char *check_box(int ePBC,matrix box)
{
    const char *ptr;

    if (ePBC == -1)
        ePBC = guess_ePBC(box);

    if (ePBC == epbcNONE)
        return NULL;

    if ((box[XX][YY] != 0) || (box[XX][ZZ] != 0) || (box[YY][ZZ] != 0)) {
        ptr = "Only triclinic boxes with the first vector parallel to the x-axis and the second vector in the xy-plane are supported.";
    } else if (ePBC == epbcSCREW && (box[YY][XX] != 0 || box[ZZ][XX] != 0)) {
        ptr = "The unit cell can not have off-diagonal x-components with screw pbc";
    } else if (fabs(box[YY][XX]) > BOX_MARGIN*0.5*box[XX][XX] ||
        (ePBC != epbcXY &&
            (fabs(box[ZZ][XX]) > BOX_MARGIN*0.5*box[XX][XX] ||
                fabs(box[ZZ][YY]) > BOX_MARGIN*0.5*box[YY][YY]))) {
        ptr = "Triclinic box is too skewed.";
    } else {
        ptr = NULL;
    }

    return ptr;
}

real max_cutoff2(int ePBC,matrix box)
{
    real min_hv2,min_ss;

    /* Physical limitation of the cut-off
     * by half the length of the shortest box vector.
     */
    min_hv2 = min(0.25*norm2(box[XX]),0.25*norm2(box[YY]));
    if (ePBC != epbcXY)
        min_hv2 = min(min_hv2,0.25*norm2(box[ZZ]));

    /* Limitation to the smallest diagonal element due to optimizations:
     * checking only linear combinations of single box-vectors (2 in x)
     * in the grid search and pbc_dx is a lot faster
     * than checking all possible combinations.
     */
    if (ePBC == epbcXY) {
        min_ss = min(box[XX][XX],box[YY][YY]);
    } else {
        min_ss = min(box[XX][XX],min(box[YY][YY]-fabs(box[ZZ][YY]),box[ZZ][ZZ]));
    }

    return min(min_hv2,min_ss*min_ss);
}

/* this one is mostly harmless... */
static gmx_bool bWarnedGuess=FALSE;

int guess_ePBC(matrix box)
{
    int ePBC;

    if (box[XX][XX]>0 && box[YY][YY]>0 && box[ZZ][ZZ]>0) {
        ePBC = epbcXYZ;
    } else if (box[XX][XX]>0 && box[YY][YY]>0 && box[ZZ][ZZ]==0) {
        ePBC = epbcXY;
    } else if (box[XX][XX]==0 && box[YY][YY]==0 && box[ZZ][ZZ]==0) {
        ePBC = epbcNONE;
    } else {
        if (!bWarnedGuess) {
            fprintf(stderr,"WARNING: Unsupported box diagonal %f %f %f, "
                    "will not use periodic boundary conditions\n\n",
                    box[XX][XX],box[YY][YY],box[ZZ][ZZ]);
            bWarnedGuess = TRUE;
        }
        ePBC = epbcNONE;
    }

    if (debug)
        fprintf(debug,"Guessed pbc = %s from the box matrix\n",epbc_names[ePBC]);

    return ePBC;
}

static int correct_box_elem(FILE *fplog,int step,tensor box,int v,int d)
{
    int shift,maxshift=10;

    shift = 0;

    /* correct elem d of vector v with vector d */
    while (box[v][d] > BOX_MARGIN_CORRECT*0.5*box[d][d]) {
        if (fplog) {
            fprintf(fplog,"Step %d: correcting invalid box:\n",step);
            pr_rvecs(fplog,0,"old box",box,DIM);
        }
        rvec_dec(box[v],box[d]);
        shift--;
        if (fplog) {
            pr_rvecs(fplog,0,"new box",box,DIM);
        }
        if (shift <= -maxshift)
            gmx_fatal(FARGS,
                      "Box was shifted at least %d times. Please see log-file.",
                      maxshift);
    } 
    while (box[v][d] < -BOX_MARGIN_CORRECT*0.5*box[d][d]) {
        if (fplog) {
            fprintf(fplog,"Step %d: correcting invalid box:\n",step);
            pr_rvecs(fplog,0,"old box",box,DIM);
        }
        rvec_inc(box[v],box[d]);
        shift++;
        if (fplog) {
            pr_rvecs(fplog,0,"new box",box,DIM);
        }
        if (shift >= maxshift)
            gmx_fatal(FARGS,
                      "Box was shifted at least %d times. Please see log-file.",
                      maxshift);
    }

    return shift;
}

gmx_bool correct_box(FILE *fplog,int step,tensor box,t_graph *graph)
{
    int  zy,zx,yx,i;
    gmx_bool bCorrected;

    /* check if the box still obeys the restrictions, if not, correct it */
    zy = correct_box_elem(fplog,step,box,ZZ,YY);
    zx = correct_box_elem(fplog,step,box,ZZ,XX);
    yx = correct_box_elem(fplog,step,box,YY,XX);

    bCorrected = (zy || zx || yx);

    if (bCorrected && graph) {
        /* correct the graph */
        for(i=0; i<graph->nnodes; i++) {
            graph->ishift[i][YY] -= graph->ishift[i][ZZ]*zy;
            graph->ishift[i][XX] -= graph->ishift[i][ZZ]*zx;
            graph->ishift[i][XX] -= graph->ishift[i][YY]*yx;
        }
    }

    return bCorrected;
}

int ndof_com(t_inputrec *ir)
{
    int n=0;

    switch (ir->ePBC) {
    case epbcXYZ:
    case epbcNONE:
        n = 3;
        break;
    case epbcXY:
        n = (ir->nwall == 0 ? 3 : 2);
        break;
    case epbcSCREW:
        n = 1;
        break;
    default:
        gmx_incons("Unknown pbc in calc_nrdf");
    }
    
    return n;
}

static void low_set_pbc(t_pbc *pbc,int ePBC,ivec *dd_nc,matrix box)
{
    int  order[5]={0,-1,1,-2,2};
    int  ii,jj,kk,i,j,k,d,dd,jc,kc,npbcdim,shift;
    ivec bPBC;
    real d2old,d2new,d2new_c;
    rvec trial,pos;
    gmx_bool bXY,bUse;
    const char *ptr;

    pbc->ndim_ePBC = ePBC2npbcdim(ePBC);

    copy_mat(box,pbc->box);
    pbc->bLimitDistance = FALSE;
    pbc->max_cutoff2 = 0;
    pbc->dim = -1;

    for(i=0; (i<DIM); i++) {
        pbc->fbox_diag[i]  =  box[i][i];
        pbc->hbox_diag[i]  =  pbc->fbox_diag[i]*0.5;
        pbc->mhbox_diag[i] = -pbc->hbox_diag[i];
    }

    ptr = check_box(ePBC,box);
    if (ePBC == epbcNONE) {
        pbc->ePBCDX = epbcdxNOPBC;
    } else if (ptr) {
        fprintf(stderr,   "Warning: %s\n",ptr);
        pr_rvecs(stderr,0,"         Box",box,DIM);
        fprintf(stderr,   "         Can not fix pbc.\n");
        pbc->ePBCDX = epbcdxUNSUPPORTED;
        pbc->bLimitDistance = TRUE;
        pbc->limit_distance2 = 0;
    } else {
        if (ePBC == epbcSCREW && dd_nc) {
            /* This combinated should never appear here */
            gmx_incons("low_set_pbc called with screw pbc and dd_nc != NULL");
        }

        npbcdim = 0;
        for(i=0; i<DIM; i++) {
            if ((dd_nc && (*dd_nc)[i] > 1) || (ePBC == epbcXY && i == ZZ)) {
                bPBC[i] = 0;
            } else {
                bPBC[i] = 1;
                npbcdim++;
            }
        }
        switch (npbcdim) {
        case 1:
            /* 1D pbc is not an mdp option and it is therefore only used
             * with single shifts.
             */
            pbc->ePBCDX = epbcdx1D_RECT;
            for(i=0; i<DIM; i++)
                if (bPBC[i])
                    pbc->dim = i;
            for(i=0; i<pbc->dim; i++)
                if (pbc->box[pbc->dim][i] != 0)
                    pbc->ePBCDX = epbcdx1D_TRIC;
            break;
        case 2:
            pbc->ePBCDX = epbcdx2D_RECT;
            for(i=0; i<DIM; i++)
                if (!bPBC[i])
                    pbc->dim = i;
            for(i=0; i<DIM; i++)
                if (bPBC[i])
                    for(j=0; j<i; j++)
                        if (pbc->box[i][j] != 0)
                            pbc->ePBCDX = epbcdx2D_TRIC;
            break;
        case 3:
            if (ePBC != epbcSCREW) {
                if (TRICLINIC(box)) {
                    pbc->ePBCDX = epbcdxTRICLINIC;
                } else {
                    pbc->ePBCDX = epbcdxRECTANGULAR;
                }
            } else {
                pbc->ePBCDX = (box[ZZ][YY]==0 ? epbcdxSCREW_RECT : epbcdxSCREW_TRIC);
                if (pbc->ePBCDX == epbcdxSCREW_TRIC) {
                    fprintf(stderr,
                            "Screw pbc is not yet implemented for triclinic boxes.\n"
                            "Can not fix pbc.\n");
                    pbc->ePBCDX = epbcdxUNSUPPORTED;
                }
            }
            break;
        default:
            gmx_fatal(FARGS,"Incorrect number of pbc dimensions with DD: %d",
                      npbcdim);
        }
        pbc->max_cutoff2 = max_cutoff2(ePBC,box);

        if (pbc->ePBCDX == epbcdxTRICLINIC ||
            pbc->ePBCDX == epbcdx2D_TRIC ||
            pbc->ePBCDX == epbcdxSCREW_TRIC) {
            if (debug) {
                pr_rvecs(debug,0,"Box",box,DIM);
                fprintf(debug,"max cutoff %.3f\n",sqrt(pbc->max_cutoff2));
            }
            pbc->ntric_vec = 0;
            /* We will only use single shifts, but we will check a few
             * more shifts to see if there is a limiting distance
             * above which we can not be sure of the correct distance.
             */
            for(kk=0; kk<5; kk++) {
                k = order[kk];
                if (!bPBC[ZZ] && k != 0)
                    continue;
                for(jj=0; jj<5; jj++) {
                    j = order[jj];
                    if (!bPBC[YY] && j != 0)
                        continue;
                    for(ii=0; ii<3; ii++) {
                        i = order[ii];
                        if (!bPBC[XX] && i != 0)
                            continue;
                        /* A shift is only useful when it is trilinic */
                        if (j != 0 || k != 0) {
                            d2old = 0;
                            d2new = 0;
                            for(d=0; d<DIM; d++) {
                                trial[d] = i*box[XX][d] + j*box[YY][d] + k*box[ZZ][d];
                                /* Choose the vector within the brick around 0,0,0 that
                                 * will become the shortest due to shift try.
                                 */
                                if (d == pbc->dim) {
                                    trial[d] = 0;
                                    pos[d] = 0;
                                } else {
                                    if (trial[d] < 0)
                                        pos[d] = min( pbc->hbox_diag[d],-trial[d]);
                                    else
                                        pos[d] = max(-pbc->hbox_diag[d],-trial[d]);
                                }
                                d2old += sqr(pos[d]);
                                d2new += sqr(pos[d] + trial[d]);
                            }
                            if (BOX_MARGIN*d2new < d2old) {
                                if (j < -1 || j > 1 || k < -1 || k > 1) {
                                    /* Check if there is a single shift vector
                                     * that decreases this distance even more.
                                     */
                                    jc = 0;
                                    kc = 0;
                                    if (j < -1 || j > 1)
                                        jc = j/2;
                                    if (k < -1 || k > 1)
                                        kc = k/2;
                                    d2new_c = 0;
                                    for(d=0; d<DIM; d++)
                                        d2new_c += sqr(pos[d] + trial[d] 
                                                                      - jc*box[YY][d] - kc*box[ZZ][d]);
                                    if (d2new_c > BOX_MARGIN*d2new) {
                                        /* Reject this shift vector, as there is no a priori limit
                                         * to the number of shifts that decrease distances.
                                         */
                                        if (!pbc->bLimitDistance || d2new <  pbc->limit_distance2)
                                            pbc->limit_distance2 = d2new;
                                        pbc->bLimitDistance = TRUE;
                                    }
                                } else {
                                    /* Check if shifts with one box vector less do better */
                                    bUse = TRUE;
                                    for(dd=0; dd<DIM; dd++) {
                                        shift = (dd==0 ? i : (dd==1 ? j : k));
                                        if (shift) {
                                            d2new_c = 0;
                                            for(d=0; d<DIM; d++)
                                                d2new_c += sqr(pos[d] + trial[d] - shift*box[dd][d]);
                                            if (d2new_c <= BOX_MARGIN*d2new)
                                                bUse = FALSE;
                                        }
                                    }
                                    if (bUse) {
                                        /* Accept this shift vector. */
                                        if (pbc->ntric_vec >= MAX_NTRICVEC) {
                                            fprintf(stderr,"\nWARNING: Found more than %d triclinic correction vectors, ignoring some.\n"
                                                    "  There is probably something wrong with your box.\n",MAX_NTRICVEC);
                                            pr_rvecs(stderr,0,"         Box",box,DIM);
                                        } else {
                                            copy_rvec(trial,pbc->tric_vec[pbc->ntric_vec]);
                                            pbc->tric_shift[pbc->ntric_vec][XX] = i;
                                            pbc->tric_shift[pbc->ntric_vec][YY] = j;
                                            pbc->tric_shift[pbc->ntric_vec][ZZ] = k;
                                            pbc->ntric_vec++;
                                        }
                                    }
                                }
                                if (debug) {
                                    fprintf(debug,"  tricvec %2d = %2d %2d %2d  %5.2f %5.2f  %5.2f %5.2f %5.2f  %5.2f %5.2f %5.2f\n",
                                            pbc->ntric_vec,i,j,k,
                                            sqrt(d2old),sqrt(d2new),
                                            trial[XX],trial[YY],trial[ZZ],
                                            pos[XX],pos[YY],pos[ZZ]);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void set_pbc(t_pbc *pbc,int ePBC,matrix box)
{
    if (ePBC == -1)
        ePBC = guess_ePBC(box);

    low_set_pbc(pbc,ePBC,NULL,box);
}

t_pbc *set_pbc_dd(t_pbc *pbc,int ePBC,
                  gmx_domdec_t *dd,gmx_bool bSingleDir,matrix box)
{
    ivec nc2;
    int  npbcdim,i;

    if (dd == NULL) {
        npbcdim = DIM;
    } else {
        if (ePBC == epbcSCREW && dd->nc[XX] > 1) {
            /* The rotation has been taken care of during coordinate communication */
            ePBC = epbcXYZ;
        }
        npbcdim = 0;
        for(i=0; i<DIM; i++) {
            if (dd->nc[i] <= (bSingleDir ? 1 : 2)) {
                nc2[i] = 1;
                if (!(ePBC == epbcXY && i == ZZ))
                    npbcdim++;
            } else {
                nc2[i] = dd->nc[i];
            }
        }
    }

    if (npbcdim > 0)
        low_set_pbc(pbc,ePBC,npbcdim<DIM ? &nc2 : NULL,box);

    return (npbcdim > 0 ? pbc : NULL);
}

void pbc_dx(const t_pbc *pbc,const rvec x1, const rvec x2, rvec dx)
{
    int  i,j;
    rvec dx_start,trial;
    real d2min,d2trial;
    gmx_bool bRot;

    rvec_sub(x1,x2,dx);

    switch (pbc->ePBCDX) {
    case epbcdxRECTANGULAR:
        for(i=0; i<DIM; i++) {
            while (dx[i] > pbc->hbox_diag[i]) {
                dx[i] -= pbc->fbox_diag[i];
            }
            while (dx[i] <= pbc->mhbox_diag[i]) {
                dx[i] += pbc->fbox_diag[i];
            }
        }
        break;
    case epbcdxTRICLINIC:
        for(i=DIM-1; i>=0; i--) {
            while (dx[i] > pbc->hbox_diag[i]) {
                for (j=i; j>=0; j--)
                    dx[j] -= pbc->box[i][j];
            }
            while (dx[i] <= pbc->mhbox_diag[i]) {
                for (j=i; j>=0; j--)
                    dx[j] += pbc->box[i][j];
            }
        }
        /* dx is the distance in a rectangular box */
        d2min = norm2(dx);
        if (d2min > pbc->max_cutoff2) {
            copy_rvec(dx,dx_start);
            d2min = norm2(dx);
            /* Now try all possible shifts, when the distance is within max_cutoff
             * it must be the shortest possible distance.
             */
            i = 0;
            while ((d2min > pbc->max_cutoff2) && (i < pbc->ntric_vec)) {
                rvec_add(dx_start,pbc->tric_vec[i],trial);
                d2trial = norm2(trial);
                if (d2trial < d2min) {
                    copy_rvec(trial,dx);
                    d2min = d2trial;
                }
                i++;
            }
        }
        break;
    case epbcdx2D_RECT:
        for(i=0; i<DIM; i++) {
            if (i != pbc->dim) {
                while (dx[i] > pbc->hbox_diag[i]) {
                    dx[i] -= pbc->fbox_diag[i];
                }
                while (dx[i] <= pbc->mhbox_diag[i]) {
                    dx[i] += pbc->fbox_diag[i];
                }
            }
        }
        break;
    case epbcdx2D_TRIC:
        d2min = 0;
        for(i=DIM-1; i>=0; i--) {
            if (i != pbc->dim) {
                while (dx[i] > pbc->hbox_diag[i]) {
                    for (j=i; j>=0; j--)
                        dx[j] -= pbc->box[i][j];
                }
                while (dx[i] <= pbc->mhbox_diag[i]) {
                    for (j=i; j>=0; j--)
                        dx[j] += pbc->box[i][j];
                }
                d2min += dx[i]*dx[i];
            }
        }
        if (d2min > pbc->max_cutoff2) {
            copy_rvec(dx,dx_start);
            d2min = norm2(dx);
            /* Now try all possible shifts, when the distance is within max_cutoff
             * it must be the shortest possible distance.
             */
            i = 0;
            while ((d2min > pbc->max_cutoff2) && (i < pbc->ntric_vec)) {
                rvec_add(dx_start,pbc->tric_vec[i],trial);
                d2trial = 0;
                for(j=0; j<DIM; j++) {
                    if (j != pbc->dim) {
                        d2trial += trial[j]*trial[j];
                    }
                }
                if (d2trial < d2min) {
                    copy_rvec(trial,dx);
                    d2min = d2trial;
                }
                i++;
            }
        }
        break;
    case epbcdxSCREW_RECT:
        /* The shift definition requires x first */
        bRot = FALSE;
        while (dx[XX] > pbc->hbox_diag[XX]) {
            dx[XX] -= pbc->fbox_diag[XX];
            bRot = !bRot;
        }
        while (dx[XX] <= pbc->mhbox_diag[XX]) {
            dx[XX] += pbc->fbox_diag[YY];
            bRot = !bRot;
        }
        if (bRot) {
            /* Rotate around the x-axis in the middle of the box */
            dx[YY] = pbc->box[YY][YY] - x1[YY] - x2[YY];
            dx[ZZ] = pbc->box[ZZ][ZZ] - x1[ZZ] - x2[ZZ];
        }
        /* Normal pbc for y and z */
        for(i=YY; i<=ZZ; i++) {
            while (dx[i] > pbc->hbox_diag[i]) {
                dx[i] -= pbc->fbox_diag[i];
            }
            while (dx[i] <= pbc->mhbox_diag[i]) {
                dx[i] += pbc->fbox_diag[i];
            }
        }
        break;
    case epbcdxNOPBC:
    case epbcdxUNSUPPORTED:
        break;
    default:
        gmx_fatal(FARGS,"Internal error in pbc_dx, set_pbc has not been called");
        break;
    }
}

int pbc_dx_aiuc(const t_pbc *pbc,const rvec x1, const rvec x2, rvec dx)
{
    int  i,j,is;
    rvec dx_start,trial;
    real d2min,d2trial;
    ivec ishift,ishift_start;

    rvec_sub(x1,x2,dx);
    clear_ivec(ishift);

    switch (pbc->ePBCDX) {
    case epbcdxRECTANGULAR:
        for(i=0; i<DIM; i++) {
            if (dx[i] > pbc->hbox_diag[i]) {
                dx[i] -=  pbc->fbox_diag[i];
                ishift[i]--;
            } else if (dx[i] <= pbc->mhbox_diag[i]) {
                dx[i] +=  pbc->fbox_diag[i];
                ishift[i]++;
            }
        }
        break;
    case epbcdxTRICLINIC:
        /* For triclinic boxes the performance difference between
         * if/else and two while loops is negligible.
         * However, the while version can cause extreme delays
         * before a simulation crashes due to large forces which
         * can cause unlimited displacements.
         * Also allowing multiple shifts would index fshift beyond bounds.
         */
        for(i=DIM-1; i>=1; i--) {
            if (dx[i] > pbc->hbox_diag[i]) {
                for (j=i; j>=0; j--)
                    dx[j] -= pbc->box[i][j];
                ishift[i]--;
            } else if (dx[i] <= pbc->mhbox_diag[i]) {
                for (j=i; j>=0; j--)
                    dx[j] += pbc->box[i][j];
                ishift[i]++;
            }
        }
        /* Allow 2 shifts in x */
        if (dx[XX] > pbc->hbox_diag[XX]) {
            dx[XX] -= pbc->fbox_diag[XX];
            ishift[XX]--;
            if (dx[XX] > pbc->hbox_diag[XX]) {
                dx[XX] -= pbc->fbox_diag[XX];
                ishift[XX]--;
            }
        } else if (dx[XX] <= pbc->mhbox_diag[XX]) {
            dx[XX] += pbc->fbox_diag[XX];
            ishift[XX]++;
            if (dx[XX] <= pbc->mhbox_diag[XX]) {
                dx[XX] += pbc->fbox_diag[XX];
                ishift[XX]++;
            }
        }
        /* dx is the distance in a rectangular box */
        d2min = norm2(dx);
        if (d2min > pbc->max_cutoff2) {
            copy_rvec(dx,dx_start);
            copy_ivec(ishift,ishift_start);
            d2min = norm2(dx);
            /* Now try all possible shifts, when the distance is within max_cutoff
             * it must be the shortest possible distance.
             */
            i = 0;
            while ((d2min > pbc->max_cutoff2) && (i < pbc->ntric_vec)) {
                rvec_add(dx_start,pbc->tric_vec[i],trial);
                d2trial = norm2(trial);
                if (d2trial < d2min) {
                    copy_rvec(trial,dx);
                    ivec_add(ishift_start,pbc->tric_shift[i],ishift);
                    d2min = d2trial;
                }
                i++;
            }
        }
        break;
    case epbcdx2D_RECT:
        for(i=0; i<DIM; i++) {
            if (i != pbc->dim) {
                if (dx[i] > pbc->hbox_diag[i]) {
                    dx[i] -= pbc->fbox_diag[i];
                    ishift[i]--;
                } else if (dx[i] <= pbc->mhbox_diag[i]) {
                    dx[i] += pbc->fbox_diag[i];
                    ishift[i]++;
                }
            }
        }
        break;
    case epbcdx2D_TRIC:
        d2min = 0;
        for(i=DIM-1; i>=1; i--) {
            if (i != pbc->dim) {
                if (dx[i] > pbc->hbox_diag[i]) {
                    for (j=i; j>=0; j--)
                        dx[j] -= pbc->box[i][j];
                    ishift[i]--;
                } else if (dx[i] <= pbc->mhbox_diag[i]) {
                    for (j=i; j>=0; j--)
                        dx[j] += pbc->box[i][j];
                    ishift[i]++;
                }
                d2min += dx[i]*dx[i];
            }
        }
        if (pbc->dim != XX) {
            /* Allow 2 shifts in x */
            if (dx[XX] > pbc->hbox_diag[XX]) {
                dx[XX] -= pbc->fbox_diag[XX];
                ishift[XX]--;
                if (dx[XX] > pbc->hbox_diag[XX]) {
                    dx[XX] -= pbc->fbox_diag[XX];
                    ishift[XX]--;
                }
            } else if (dx[XX] <= pbc->mhbox_diag[XX]) {
                dx[XX] += pbc->fbox_diag[XX];
                ishift[XX]++;
                if (dx[XX] <= pbc->mhbox_diag[XX]) {
                    dx[XX] += pbc->fbox_diag[XX];
                    ishift[XX]++;
                }
            }
            d2min += dx[XX]*dx[XX];
        }
        if (d2min > pbc->max_cutoff2) {
            copy_rvec(dx,dx_start);
            copy_ivec(ishift,ishift_start);
            /* Now try all possible shifts, when the distance is within max_cutoff
             * it must be the shortest possible distance.
             */
            i = 0;
            while ((d2min > pbc->max_cutoff2) && (i < pbc->ntric_vec)) {
                rvec_add(dx_start,pbc->tric_vec[i],trial);
                d2trial = 0;
                for(j=0; j<DIM; j++) {
                    if (j != pbc->dim) {
                        d2trial += trial[j]*trial[j];
                    }
                }
                if (d2trial < d2min) {
                    copy_rvec(trial,dx);
                    ivec_add(ishift_start,pbc->tric_shift[i],ishift);
                    d2min = d2trial;
                }
                i++;
            }
        }
        break;
    case epbcdx1D_RECT:
        i = pbc->dim;
        if (dx[i] > pbc->hbox_diag[i]) {
            dx[i] -= pbc->fbox_diag[i];
            ishift[i]--;
        } else if (dx[i] <= pbc->mhbox_diag[i]) {
            dx[i] += pbc->fbox_diag[i];
            ishift[i]++;
        }
        break;
    case epbcdx1D_TRIC:
        i = pbc->dim;
        if (dx[i] > pbc->hbox_diag[i]) {
            rvec_dec(dx,pbc->box[i]);
            ishift[i]--;
        } else if (dx[i] <= pbc->mhbox_diag[i]) {
            rvec_inc(dx,pbc->box[i]);
            ishift[i]++;
        }
        break;
    case epbcdxSCREW_RECT:
        /* The shift definition requires x first */
        if (dx[XX] > pbc->hbox_diag[XX]) {
            dx[XX] -= pbc->fbox_diag[XX];
            ishift[XX]--;
        } else if (dx[XX] <= pbc->mhbox_diag[XX]) {
            dx[XX] += pbc->fbox_diag[XX];
            ishift[XX]++;
        }
        if (ishift[XX] == 1 || ishift[XX] == -1) {
            /* Rotate around the x-axis in the middle of the box */
            dx[YY] = pbc->box[YY][YY] - x1[YY] - x2[YY];
            dx[ZZ] = pbc->box[ZZ][ZZ] - x1[ZZ] - x2[ZZ];
        }
        /* Normal pbc for y and z */
        for(i=YY; i<=ZZ; i++) {
            if (dx[i] > pbc->hbox_diag[i]) {
                dx[i] -= pbc->fbox_diag[i];
                ishift[i]--;
            } else if (dx[i] <= pbc->mhbox_diag[i]) {
                dx[i] += pbc->fbox_diag[i];
                ishift[i]++;
            }
        }
        break;
    case epbcdxNOPBC:
    case epbcdxUNSUPPORTED:
        break;
    default:
        gmx_fatal(FARGS,"Internal error in pbc_dx_aiuc, set_pbc_dd or set_pbc has not been called");
        break;
    }

    is = IVEC2IS(ishift);
    if (debug)
    {
        range_check_mesg(is,0,SHIFTS,"PBC shift vector index range check.");
    }

    return is; 
}

void pbc_dx_d(const t_pbc *pbc,const dvec x1, const dvec x2, dvec dx)
{
    int  i,j;
    dvec dx_start,trial;
    double d2min,d2trial;
    gmx_bool bRot;

    dvec_sub(x1,x2,dx);

    switch (pbc->ePBCDX) {
    case epbcdxRECTANGULAR:
    case epbcdx2D_RECT:
        for(i=0; i<DIM; i++) {
            if (i != pbc->dim) {
                while (dx[i] > pbc->hbox_diag[i]) {
                    dx[i] -= pbc->fbox_diag[i];
                }
                while (dx[i] <= pbc->mhbox_diag[i]) {
                    dx[i] += pbc->fbox_diag[i];
                }
            }
        }
        break;
    case epbcdxTRICLINIC:
    case epbcdx2D_TRIC:
        d2min = 0;
        for(i=DIM-1; i>=0; i--) {
            if (i != pbc->dim) {
                while (dx[i] > pbc->hbox_diag[i]) {
                    for (j=i; j>=0; j--)
                        dx[j] -= pbc->box[i][j];
                }
                while (dx[i] <= pbc->mhbox_diag[i]) {
                    for (j=i; j>=0; j--)
                        dx[j] += pbc->box[i][j];
                }
                d2min += dx[i]*dx[i];
            }
        }
        if (d2min > pbc->max_cutoff2) {
            copy_dvec(dx,dx_start);
            /* Now try all possible shifts, when the distance is within max_cutoff
             * it must be the shortest possible distance.
             */
            i = 0;
            while ((d2min > pbc->max_cutoff2) && (i < pbc->ntric_vec)) {
                for(j=0; j<DIM; j++) {
                    trial[j] = dx_start[j] + pbc->tric_vec[i][j];
                }
                d2trial = 0;
                for(j=0; j<DIM; j++) {
                    if (j != pbc->dim) {
                        d2trial += trial[j]*trial[j];
                    }
                }
                if (d2trial < d2min) {
                    copy_dvec(trial,dx);
                    d2min = d2trial;
                }
                i++;
            }
        }
        break;
    case epbcdxSCREW_RECT:
        /* The shift definition requires x first */
        bRot = FALSE;
        while (dx[XX] > pbc->hbox_diag[XX]) {
            dx[XX] -= pbc->fbox_diag[XX];
            bRot = !bRot;
        }
        while (dx[XX] <= pbc->mhbox_diag[XX]) {
            dx[XX] += pbc->fbox_diag[YY];
            bRot = !bRot;
        }
        if (bRot) {
            /* Rotate around the x-axis in the middle of the box */
            dx[YY] = pbc->box[YY][YY] - x1[YY] - x2[YY];
            dx[ZZ] = pbc->box[ZZ][ZZ] - x1[ZZ] - x2[ZZ];
        }
        /* Normal pbc for y and z */
        for(i=YY; i<=ZZ; i++) {
            while (dx[i] > pbc->hbox_diag[i]) {
                dx[i] -= pbc->fbox_diag[i];
            }
            while (dx[i] <= pbc->mhbox_diag[i]) {
                dx[i] += pbc->fbox_diag[i];
            }
        }
        break;
    case epbcdxNOPBC:
    case epbcdxUNSUPPORTED:
        break;
    default:
        gmx_fatal(FARGS,"Internal error in pbc_dx, set_pbc has not been called");
        break;
    }
}

gmx_bool image_rect(ivec xi,ivec xj,ivec box_size,real rlong2,int *shift,real *r2)
{
    int 	m,t;
    int 	dx,b,b_2;
    real  dxr,rij2;

    rij2=0.0;
    t=0;
    for(m=0; (m<DIM); m++) {
        dx=xi[m]-xj[m];
        t*=DIM;
        b=box_size[m];
        b_2=b/2;
        if (dx < -b_2) {
            t+=2;
            dx+=b;
        }
        else if (dx > b_2)
            dx-=b;
        else
            t+=1;
        dxr=dx;
        rij2+=dxr*dxr;
        if (rij2 >= rlong2) 
            return FALSE;
    }

    *shift = t;
    *r2 = rij2;
    return TRUE;
}

gmx_bool image_cylindric(ivec xi,ivec xj,ivec box_size,real rlong2,
                     int *shift,real *r2)
{
    int 	m,t;
    int 	dx,b,b_2;
    real  dxr,rij2;

    rij2=0.0;
    t=0;
    for(m=0; (m<DIM); m++) {
        dx=xi[m]-xj[m];
        t*=DIM;
        b=box_size[m];
        b_2=b/2;
        if (dx < -b_2) {
            t+=2;
            dx+=b;
        }
        else if (dx > b_2)
            dx-=b;
        else
            t+=1;

        dxr=dx;
        rij2+=dxr*dxr;
        if (m < ZZ) {
            if (rij2 >= rlong2) 
                return FALSE;
        }
    }

    *shift = t;
    *r2 = rij2;
    return TRUE;
}

void calc_shifts(matrix box,rvec shift_vec[])
{
    int k,l,m,d,n,test;

    n=0;
    for(m = -D_BOX_Z; m <= D_BOX_Z; m++)
        for(l = -D_BOX_Y; l <= D_BOX_Y; l++) 
            for(k = -D_BOX_X; k <= D_BOX_X; k++,n++) {
                test = XYZ2IS(k,l,m);
                if (n != test)
                    gmx_incons("inconsistent shift numbering");
                for(d = 0; d < DIM; d++)
                    shift_vec[n][d] = k*box[XX][d] + l*box[YY][d] + m*box[ZZ][d];
            }
}

void calc_box_center(int ecenter,matrix box,rvec box_center)
{
    int d,m;

    clear_rvec(box_center);
    switch (ecenter) {
    case ecenterTRIC:
        for(m=0; (m<DIM); m++)  
            for(d=0; d<DIM; d++)
                box_center[d] += 0.5*box[m][d];
        break;
    case ecenterRECT:
        for(d=0; d<DIM; d++)
            box_center[d] = 0.5*box[d][d];
        break;
    case ecenterZERO:
        break;
    default:
        gmx_fatal(FARGS,"Unsupported value %d for ecenter",ecenter);
    }
}

void calc_triclinic_images(matrix box,rvec img[])
{
    int i;

    /* Calculate 3 adjacent images in the xy-plane */
    copy_rvec(box[0],img[0]);
    copy_rvec(box[1],img[1]);
    if (img[1][XX] < 0)
        svmul(-1,img[1],img[1]);
    rvec_sub(img[1],img[0],img[2]);

    /* Get the next 3 in the xy-plane as mirror images */
    for(i=0; i<3; i++)
        svmul(-1,img[i],img[3+i]);

    /* Calculate the first 4 out of xy-plane images */
    copy_rvec(box[2],img[6]);
    if (img[6][XX] < 0)
        svmul(-1,img[6],img[6]);
    for(i=0; i<3; i++)
        rvec_add(img[6],img[i+1],img[7+i]);

    /* Mirror the last 4 from the previous in opposite rotation */
    for(i=0; i<4; i++)
        svmul(-1,img[6 + (2+i) % 4],img[10+i]);
}

void calc_compact_unitcell_vertices(int ecenter,matrix box,rvec vert[])
{
    rvec img[NTRICIMG],box_center;
    int n,i,j,tmp[4],d;

    calc_triclinic_images(box,img);

    n=0;
    for(i=2; i<=5; i+=3) {
        tmp[0] = i-1;
        if (i==2)
            tmp[1] = 8;
        else 
            tmp[1] = 6;
        tmp[2] = (i+1) % 6;
        tmp[3] = tmp[1]+4;
        for(j=0; j<4; j++) {
            for(d=0; d<DIM; d++)
                vert[n][d] = img[i][d]+img[tmp[j]][d]+img[tmp[(j+1)%4]][d];
            n++;
        }
    }
    for(i=7; i<=13; i+=6) {
        tmp[0] = (i-7)/2;
        tmp[1] = tmp[0]+1;
        if (i==7)
            tmp[2] = 8;
        else
            tmp[2] = 10;
        tmp[3] = i-1;
        for(j=0; j<4; j++) {
            for(d=0; d<DIM; d++)
                vert[n][d] = img[i][d]+img[tmp[j]][d]+img[tmp[(j+1)%4]][d];
            n++;
        }
    }
    for(i=9; i<=11; i+=2) {
        if (i==9)
            tmp[0] = 3;
        else
            tmp[0] = 0;
        tmp[1] = tmp[0]+1;
        if (i==9)
            tmp[2] = 6;
        else
            tmp[2] = 12;
        tmp[3] = i-1;
        for(j=0; j<4; j++) {
            for(d=0; d<DIM; d++)
                vert[n][d] = img[i][d]+img[tmp[j]][d]+img[tmp[(j+1)%4]][d];
            n++;
        }
    }

    calc_box_center(ecenter,box,box_center);
    for(i=0; i<NCUCVERT; i++)
        for(d=0; d<DIM; d++)
            vert[i][d] = vert[i][d]*0.25+box_center[d];
}

int *compact_unitcell_edges()
{
    /* this is an index in vert[] (see calc_box_vertices) */
    /*static int edge[NCUCEDGE*2];*/
    int *edge;
    static const int hexcon[24] = { 0,9, 1,19, 2,15, 3,21, 
        4,17, 5,11, 6,23, 7,13,
        8,20, 10,18, 12,16, 14,22 };
    int e,i,j;
    gmx_bool bFirst = TRUE;

    snew(edge,NCUCEDGE*2);

    if (bFirst) {
        e = 0;
        for(i=0; i<6; i++)
            for(j=0; j<4; j++) {
                edge[e++] = 4*i + j;
                edge[e++] = 4*i + (j+1) % 4;
            }
        for(i=0; i<12*2; i++)
            edge[e++] = hexcon[i];

        bFirst = FALSE;
    }

    return edge;
}

void put_atom_in_box(matrix box,rvec x)
{
    int i,m,d;

    for(m=DIM-1; m>=0; m--) {
        while (x[m] < 0) 
            for(d=0; d<=m; d++)
                x[d] += box[m][d];
        while (x[m] >= box[m][m])
            for(d=0; d<=m; d++)
                x[d] -= box[m][d];
    }
}

void put_atoms_in_box(matrix box,int natoms,rvec x[])
{
    int i,m,d;

    for(i=0; (i<natoms); i++)
        put_atom_in_box(box,x[i]);
}

void put_atoms_in_triclinic_unitcell(int ecenter,matrix box,
                                     int natoms,rvec x[])
{
    rvec   box_center,shift_center;
    real   shm01,shm02,shm12,shift;
    int    i,m,d;

    calc_box_center(ecenter,box,box_center);

    /* The product of matrix shm with a coordinate gives the shift vector
     which is required determine the periodic cell position */
    shm01 = box[1][0]/box[1][1];
    shm02 = (box[1][1]*box[2][0] - box[2][1]*box[1][0])/(box[1][1]*box[2][2]);
    shm12 = box[2][1]/box[2][2];

    clear_rvec(shift_center);
    for(d=0; d<DIM; d++)
        rvec_inc(shift_center,box[d]);
    svmul(0.5,shift_center,shift_center);
    rvec_sub(box_center,shift_center,shift_center);

    shift_center[0] = shm01*shift_center[1] + shm02*shift_center[2];
    shift_center[1] = shm12*shift_center[2];
    shift_center[2] = 0;

    for(i=0; (i<natoms); i++)
        for(m=DIM-1; m>=0; m--) {
            shift = shift_center[m];
            if (m == 0) {
                shift += shm01*x[i][1] + shm02*x[i][2];
            } else if (m == 1) {
                shift += shm12*x[i][2];
            }
            while (x[i][m]-shift < 0)
                for(d=0; d<=m; d++)
                    x[i][d] += box[m][d];
            while (x[i][m]-shift >= box[m][m])
                for(d=0; d<=m; d++)
                    x[i][d] -= box[m][d];
        }
}

const char *
put_atoms_in_compact_unitcell(int ePBC,int ecenter,matrix box,
                              int natoms,rvec x[])
                              {
    t_pbc pbc;
    rvec box_center,dx;
    int  i;

    set_pbc(&pbc,ePBC,box);
    calc_box_center(ecenter,box,box_center);
    for(i=0; i<natoms; i++) {
        pbc_dx(&pbc,x[i],box_center,dx);
        rvec_add(box_center,dx,x[i]);
    }

    return pbc.bLimitDistance ?
        "WARNING: Could not put all atoms in the compact unitcell\n"
        : NULL;
                              }

