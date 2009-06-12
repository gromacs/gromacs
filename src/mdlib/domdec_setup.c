/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
 * 
 * This file is part of Gromacs        Copyright (c) 1991-2008
 * David van der Spoel, Erik Lindahl, Berk Hess, University of Groningen.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org
 * 
 * And Hey:
 * Gnomes, ROck Monsters And Chili Sauce
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include "domdec.h"
#include "network.h"
#include "perf_est.h"
#include "physics.h"
#include "pme.h"
#include "smalloc.h"
#include "typedefs.h"
#include "vec.h"

/* Margin for setting up the DD grid */
#define DD_GRID_MARGIN_PRES_SCALE 1.05

static int factorize(int n,int **fac,int **mfac)
{
	int d,ndiv;

	/* Decompose n in factors */
	snew(*fac,n/2);
	snew(*mfac,n/2);
	d = 2;
	ndiv = 0;
	while (n > 1)
		{
		while (n % d == 0)
			{
			if (ndiv == 0 || (*fac)[ndiv-1] != d)
				{
				ndiv++;
				(*fac)[ndiv-1] = d;
				}
			(*mfac)[ndiv-1]++;
			n /= d;
			}
		d++;
		}
	
	return ndiv;
}

static bool fits_pme_ratio(int nnodes,int npme,float ratio)
{
    return ((double)npme/(double)nnodes > 0.95*ratio); 
}

static bool fits_pme_perf(FILE *fplog,
						  t_inputrec *ir,matrix box,gmx_mtop_t *mtop,
						  int nnodes,int npme,float ratio)
{
    /* Does this division gives a reasonable PME load? */
    return (fits_pme_ratio(nnodes,npme,ratio) &&
			pme_inconvenient_nnodes(ir->nkx,ir->nky,npme) <= 1);
}

static int guess_npme(FILE *fplog,gmx_mtop_t *mtop,t_inputrec *ir,matrix box,
					  int nnodes)
{
	float ratio;
	int  npme,nkx,nky,ndiv,*div,*mdiv,ldiv;
	t_inputrec ir_try;
	
	ratio = pme_load_estimate(mtop,ir,box);
	
	if (fplog)
    {
		fprintf(fplog,"Guess for relative PME load: %.2f\n",ratio);
    }
	
	/* We assume the optimal node ratio is close to the load ratio.
	 * The communication load is neglected,
	 * but (hopefully) this will balance out between PP and PME.
	 */
	
    if (!fits_pme_ratio(nnodes,nnodes/2,ratio))
    {
        /* We would need more than nnodes/2 PME only nodes,
         * which is not possible. Since the PME load is very high,
         * we will not loose much performance when all nodes do PME.
         */

        return 0;
    }

    /* First try to find npme as a factor of nnodes up to nnodes/3 */
	npme = 1;
    while (npme <= nnodes/3) {
        if (nnodes % npme == 0) {
            /* Note that fits_perf might change the PME grid,
             * in the current implementation it does not.
             */
            if (fits_pme_perf(fplog,ir,box,mtop,nnodes,npme,ratio))
			{
				break;
			}
        }
        npme++;
    }
    if (npme > nnodes/3)
    {
        /* Try any possible number for npme */
        npme = 1;
        while (npme <= nnodes/2)
        {
            ndiv = factorize(nnodes-npme,&div,&mdiv);
            ldiv = div[ndiv-1];
            sfree(div);
            sfree(mdiv);
            /* Only use this value if nnodes-npme does not have
             * a large prime factor (5 y, 7 n, 14 n, 15 y).
             */
            if (ldiv <= 3 + (int)(pow(nnodes-npme,1.0/3.0) + 0.5))
            {
                /* Note that fits_perf may change the PME grid */
                if (fits_pme_perf(fplog,ir,box,mtop,nnodes,npme,ratio))
                {
                    break;
                }
            }
            npme++;
        }
    }
    if (npme > nnodes/2)
    {
        gmx_fatal(FARGS,"Could not find an appropriate number of separate PME nodes. i.e. >= %5f*#nodes (%d) and <= #nodes/2 (%d) and reasonable performance wise (grid_x=%d, grid_y=%d).\n"
                  "Use the -npme option of mdrun or change the number of processors or the PME grid dimensions, see the manual for details.",
                  ratio,(int)(0.95*ratio*nnodes+0.5),nnodes/2,ir->nkx,ir->nky);
        /* Keep the compiler happy */
        npme = 0;
    }
    else
    {
        if (fplog)
        {
            fprintf(fplog,
                    "Will use %d particle-particle and %d PME only nodes\n"
                    "This is a guess, check the performance at the end of the log file\n",
                    nnodes-npme,npme);
        }
        fprintf(stderr,"\n"
                "Will use %d particle-particle and %d PME only nodes\n"
                "This is a guess, check the performance at the end of the log file\n",
                nnodes-npme,npme);
    }
    
    return npme;
}

static int lcd(int n1,int n2)
{
    int d,i;
    
    d = 1;
    for(i=2; (i<=n1 && i<=n2); i++)
    {
        if (n1 % i == 0 && n2 % i == 0)
        {
            d = i;
        }
    }
    
  return d;
}

real comm_box_frac(ivec dd_nc,real cutoff,gmx_ddbox_t *ddbox)
{
    int  i,j,k,npp;
    rvec bt,nw;
    real comm_vol;

    for(i=0; i<DIM; i++)
    {
        bt[i] = ddbox->box_size[i]*ddbox->skew_fac[i];
        nw[i] = dd_nc[i]*cutoff/bt[i];
    }

    npp = 1;
    comm_vol = 0;
    for(i=0; i<DIM; i++)
    {
        if (dd_nc[i] > 1)
        {
            npp *= dd_nc[i];
            comm_vol += nw[i];
            for(j=i+1; j<DIM; j++)
            {
                if (dd_nc[j] > 1)
                {
                    comm_vol += nw[i]*nw[j]*M_PI/4;
                    for(k=j+1; k<DIM; k++)
                    {
                        if (dd_nc[k] > 1)
                        {
                            comm_vol += nw[i]*nw[j]*nw[k]*M_PI/6;
                        }
                    }
                }
            }
        }
    }
    /* Normalize by the number of PP nodes */
    comm_vol /= npp;
   
    return comm_vol;
}

static float comm_cost_est(gmx_domdec_t *dd,real limit,real cutoff,
                           matrix box,gmx_ddbox_t *ddbox,t_inputrec *ir,
                           float pbcdxr,
                           int npme,ivec nc)
{
    int  i,j,k,npp;
    rvec bt;
    float comm_vol,comm_vol_pme,cost_pbcdx;
    /* This is the cost of a pbc_dx call relative to the cost
     * of communicating the coordinate and force of an atom.
     * This will be machine dependent.
     * These factors are for x86 with SMP or Infiniband.
     */
    float pbcdx_rect_fac = 0.1;
    float pbcdx_tric_fac = 0.2;
    
    /* Check the DD algorithm restrictions */
    if ((ir->ePBC == epbcXY && ir->nwall < 2 && nc[ZZ] > 1) ||
        (ir->ePBC == epbcSCREW && (nc[XX] == 1 || nc[YY] > 1 || nc[ZZ] > 1)))
    {
        return -1;
    }
    
    /* Check if the triclinic requirements are met */
    for(i=0; i<DIM; i++)
    {
        for(j=i+1; j<ddbox->npbcdim; j++)
        {
            if (box[j][i] != 0 || ir->deform[j][i] != 0 ||
                (ir->epc != epcNO && ir->compress[j][i] != 0))
            {
                if (nc[j] > 1 && nc[i] == 1)
                {
                    return -1;
                }
            }
        }
    }
    
    npp = 1;
    for(i=0; i<DIM; i++)
    {
        npp *= nc[i];
        bt[i] = ddbox->box_size[i]*ddbox->skew_fac[i];
        
        /* Without PBC there are no cell size limits with 2 cells */
        if (!(i >= ddbox->npbcdim && nc[i] <= 2) && bt[i] < nc[i]*limit)
        {
            return -1;
        }
    }
    
    /* When two dimensions are (nearly) equal, use more cells
     * for the smallest index, so the decomposition does not
     * depend sensitively on the rounding of the box elements.
     */
    for(i=0; i<DIM; i++)
    {
        if (npme == 0 || i != XX)
        {
            for(j=i+1; j<DIM; j++)
            {
                if (fabs(bt[j] - bt[i]) < 0.01*bt[i] && nc[j] > nc[i])
                {
                    return -1;
                }
            }
        }
    }
    
    comm_vol = comm_box_frac(nc,cutoff,ddbox);

    /* Determine the largest volume that a PME only needs to communicate */
    comm_vol_pme = 0;
    if ((npme > 0) && (nc[XX] % npme != 0))
    {
        if (nc[XX] > npme)
        {
            comm_vol_pme = (npme==2 ? 1.0/3.0 : 0.5);
        }
        else
        {
            comm_vol_pme = 1.0 - lcd(nc[XX],npme)/(double)npme;
        }
        /* Normalize by the number of PME only nodes */
        comm_vol_pme /= npme;
    }
    
    /* Add cost of pbc_dx for bondeds */
    cost_pbcdx = 0;
    if ((nc[XX] == 1 || nc[YY] == 1) || (nc[ZZ] == 1 && ir->ePBC != epbcXY))
    {
        if ((ddbox->tric_dir[XX] && nc[XX] == 1) ||
            (ddbox->tric_dir[YY] && nc[YY] == 1))
        {
            cost_pbcdx = pbcdxr*pbcdx_tric_fac/npp;
        }
        else
        {
            cost_pbcdx = pbcdxr*pbcdx_rect_fac/npp;
        }
    }
    
    if (debug)
    {
        fprintf(debug,
                "nc %2d %2d %2d vol pp %6.4f pbcdx %6.4f pme %6.4f tot %6.4f\n",
                nc[XX],nc[YY],nc[ZZ],
                comm_vol,cost_pbcdx,comm_vol_pme,
                comm_vol + cost_pbcdx + comm_vol_pme);
    }
    
    return comm_vol + cost_pbcdx + comm_vol_pme;
}

static void assign_factors(gmx_domdec_t *dd,
                           real limit,real cutoff,
                           matrix box,gmx_ddbox_t *ddbox,t_inputrec *ir,
                           float pbcdxr,int npme,
                           int ndiv,int *div,int *mdiv,ivec ir_try,ivec opt)
{
    int x,y,z,i;
    float ce;
    
    if (ndiv == 0)
    {
        ce = comm_cost_est(dd,limit,cutoff,box,ddbox,ir,pbcdxr,npme,ir_try);
        if (ce >= 0 && (opt[XX] == 0 ||
                        ce < comm_cost_est(dd,limit,cutoff,box,ddbox,ir,pbcdxr,
                                           npme,opt)))
        {
            copy_ivec(ir_try,opt);
        }
        
        return;
    }
    
    for(x=mdiv[0]; x>=0; x--)
    {
        for(i=0; i<x; i++)
        {
            ir_try[XX] *= div[0];
        }
        for(y=mdiv[0]-x; y>=0; y--)
        {
            for(i=0; i<y; i++)
            {
                ir_try[YY] *= div[0];
            }
            for(i=0; i<mdiv[0]-x-y; i++)
            {
                ir_try[ZZ] *= div[0];
            }
            
            /* recurse */
            assign_factors(dd,limit,cutoff,box,ddbox,ir,pbcdxr,npme,
                           ndiv-1,div+1,mdiv+1,ir_try,opt);
            
            for(i=0; i<mdiv[0]-x-y; i++)
            {
                ir_try[ZZ] /= div[0];
            }
            for(i=0; i<y; i++)
            {
                ir_try[YY] /= div[0];
            }
        }
        for(i=0; i<x; i++)
        {
            ir_try[XX] /= div[0];
        }
    }
}

static real optimize_ncells(FILE *fplog,
                            int nnodes_tot,int npme_only,
                            bool bDynLoadBal,real dlb_scale,
                            gmx_mtop_t *mtop,matrix box,gmx_ddbox_t *ddbox,
                            t_inputrec *ir,
                            gmx_domdec_t *dd,
                            real cellsize_limit,real cutoff,
                            bool bInterCGBondeds,bool bInterCGMultiBody,
                            ivec nc)
{
    int npp,npme,ndiv,*div,*mdiv,d,nmax;
    bool bExcl_pbcdx;
    float pbcdxr;
    real limit;
    ivec itry;
    
    limit  = cellsize_limit;
    
    dd->nc[XX] = 1;
    dd->nc[YY] = 1;
    dd->nc[ZZ] = 1;

    npp = nnodes_tot - npme_only;
    if (EEL_PME(ir->coulombtype))
    {
        npme = (npme_only > 0 ? npme_only : npp);
    }
    else
    {
        npme = 0;
    }
    
    if (bInterCGBondeds)
    {
        /* For Ewald exclusions pbc_dx is not called */
        bExcl_pbcdx =
            (EEL_EXCL_FORCES(ir->coulombtype) && !EEL_FULL(ir->coulombtype));
        pbcdxr = (double)n_bonded_dx(mtop,bExcl_pbcdx)/(double)mtop->natoms;
    }
    else
    {
        /* Every molecule is a single charge group: no pbc required */
        pbcdxr = 0;
    }
    /* Add a margin for DLB and/or pressure scaling */
    if (bDynLoadBal)
    {
        if (dlb_scale >= 1.0)
        {
            gmx_fatal(FARGS,"The value for option -dds should be smaller than 1");
        }
        if (fplog)
        {
            fprintf(fplog,"Scaling the initial minimum size with 1/%g (option -dds) = %g\n",dlb_scale,1/dlb_scale);
        }
        limit /= dlb_scale;
    }
    else if (ir->epc != epcNO)
    {
        if (fplog)
        {
            fprintf(fplog,"To account for pressure scaling, scaling the initial minimum size with %g\n",DD_GRID_MARGIN_PRES_SCALE);
            limit *= DD_GRID_MARGIN_PRES_SCALE;
        }
    }
    
    if (fplog)
    {
        fprintf(fplog,"Optimizing the DD grid for %d cells with a minimum initial size of %.3f nm\n",npp,limit);

        if (limit > 0)
        {
            fprintf(fplog,"The maximum allowed number of cells is:");
            for(d=0; d<DIM; d++)
            {
                nmax = (int)(ddbox->box_size[d]*ddbox->skew_fac[d]/limit);
                if (d >= ddbox->npbcdim && nmax < 2)
                {
                    nmax = 2;
                }
                fprintf(fplog," %c %d",'X' + d,nmax);
            }
            fprintf(fplog,"\n");
        }
    }
    
    if (debug)
    {
        fprintf(debug,"Average nr of pbc_dx calls per atom %.2f\n",pbcdxr);
    }
    
    /* Decompose npp in factors */
    ndiv = factorize(npp,&div,&mdiv);
    
    itry[XX] = 1;
    itry[YY] = 1;
    itry[ZZ] = 1;
    clear_ivec(nc);
    assign_factors(dd,limit,cutoff,box,ddbox,ir,pbcdxr,
                   npme,ndiv,div,mdiv,itry,nc);
    
    sfree(div);
    sfree(mdiv);
    
    return limit;
}

real dd_choose_grid(FILE *fplog,
                    t_commrec *cr,gmx_domdec_t *dd,t_inputrec *ir,
                    gmx_mtop_t *mtop,matrix box,gmx_ddbox_t *ddbox,
                    bool bDynLoadBal,real dlb_scale,
                    real cellsize_limit,real cutoff_dd,
                    bool bInterCGBondeds,bool bInterCGMultiBody)
{
    int  npme,nkx,nky;
    real limit;
    
    if (MASTER(cr))
    {
        if (EEL_PME(ir->coulombtype))
        {
            if (cr->npmenodes >= 0)
            {
                if (cr->nnodes <= 2 && cr->npmenodes > 0)
                {
                    gmx_fatal(FARGS,
                              "Can not have separate PME nodes with 2 or less nodes");
                }
            }
            else
            {
                if (cr->nnodes < 12 &&
                    pme_inconvenient_nnodes(ir->nkx,ir->nky,cr->nnodes) == 0)
                {
                    cr->npmenodes = 0;
                }
                else
                {
                    cr->npmenodes = guess_npme(fplog,mtop,ir,box,cr->nnodes);
                }
            }
            if (fplog)
            {
                fprintf(fplog,"Using %d separate PME nodes\n",cr->npmenodes);
            }
        }
        else
        {
            if (cr->npmenodes < 0)
            {
                cr->npmenodes = 0;
            }
        }
        
        limit = optimize_ncells(fplog,cr->nnodes,cr->npmenodes,
                                bDynLoadBal,dlb_scale,
                                mtop,box,ddbox,ir,dd,
                                cellsize_limit,cutoff_dd,
                                bInterCGBondeds,bInterCGMultiBody,
                                dd->nc);
    }
    else
    {
        limit = 0;
    }
    /* Communicate the information set by the master to all nodes */
    gmx_bcast(sizeof(dd->nc),dd->nc,cr);
    if (EEL_PME(ir->coulombtype))
    {
        gmx_bcast(sizeof(ir->nkx),&ir->nkx,cr);
        gmx_bcast(sizeof(ir->nky),&ir->nky,cr);
        gmx_bcast(sizeof(cr->npmenodes),&cr->npmenodes,cr);
    }
    else
    {
        cr->npmenodes = 0;
    }
    
    return limit;
}
