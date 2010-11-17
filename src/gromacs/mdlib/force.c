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
 * GROwing Monsters And Cloning Shrimps
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>
#include <string.h>
#include "sysstuff.h"
#include "typedefs.h"
#include "macros.h"
#include "smalloc.h"
#include "macros.h"
#include "physics.h"
#include "force.h"
#include "nonbonded.h"
#include "names.h"
#include "network.h"
#include "pbc.h"
#include "ns.h"
#include "nrnb.h"
#include "bondf.h"
#include "mshift.h"
#include "txtdump.h"
#include "coulomb.h"
#include "pppm.h"
#include "pme.h"
#include "mdrun.h"
#include "domdec.h"
#include "partdec.h"
#include "qmmm.h"
#include "mpelogging.h"


void ns(FILE *fp,
        t_forcerec *fr,
        rvec       x[],
        matrix     box,
        gmx_groups_t *groups,
        t_grpopts  *opts,
        gmx_localtop_t *top,
        t_mdatoms  *md,
        t_commrec  *cr,
        t_nrnb     *nrnb,
        real       lambda,
        real       *dvdlambda,
        gmx_grppairener_t *grppener,
        gmx_bool       bFillGrid,
        gmx_bool       bDoLongRange,
        gmx_bool       bDoForces,
        rvec       *f)
{
  char   *ptr;
  int    nsearch;

  GMX_MPE_LOG(ev_ns_start);
  if (!fr->ns.nblist_initialized)
  {
      init_neighbor_list(fp, fr, md->homenr);
  }
    
  if (fr->bTwinRange) 
    fr->nlr=0;

    nsearch = search_neighbours(fp,fr,x,box,top,groups,cr,nrnb,md,
                                lambda,dvdlambda,grppener,
                                bFillGrid,bDoLongRange,
                                bDoForces,f);
  if (debug)
    fprintf(debug,"nsearch = %d\n",nsearch);
    
  /* Check whether we have to do dynamic load balancing */
  /*if ((nsb->nstDlb > 0) && (mod(step,nsb->nstDlb) == 0))
    count_nb(cr,nsb,&(top->blocks[ebCGS]),nns,fr->nlr,
    &(top->idef),opts->ngener);
  */
  if (fr->ns.dump_nl > 0)
    dump_nblist(fp,cr,fr,fr->ns.dump_nl);

  GMX_MPE_LOG(ev_ns_finish);
}

void do_force_lowlevel(FILE       *fplog,   gmx_large_int_t step,
                       t_forcerec *fr,      t_inputrec *ir,
                       t_idef     *idef,    t_commrec  *cr,
                       t_nrnb     *nrnb,    gmx_wallcycle_t wcycle,
                       t_mdatoms  *md,
                       t_grpopts  *opts,
                       rvec       x[],      history_t  *hist,
                       rvec       f[],
                       gmx_enerdata_t *enerd,
                       t_fcdata   *fcd,
                       gmx_mtop_t     *mtop,
                       gmx_localtop_t *top,
                       gmx_genborn_t *born,
                       t_atomtypes *atype,
                       gmx_bool       bBornRadii,
                       matrix     box,
                       real       lambda,  
                       t_graph    *graph,
                       t_blocka   *excl,    
                       rvec       mu_tot[],
                       int        flags,
                       float      *cycles_pme)
{
    int     i,status;
    int     donb_flags;
    gmx_bool    bDoEpot,bSepDVDL,bSB;
    int     pme_flags;
    matrix  boxs;
    rvec    box_size;
    real    dvdlambda,Vsr,Vlr,Vcorr=0,vdip,vcharge;
    t_pbc   pbc;
    real    dvdgb;
    char    buf[22];
    gmx_enerdata_t ed_lam;
    double  lam_i;
    real    dvdl_dum;

#ifdef GMX_MPI
    double  t0=0.0,t1,t2,t3; /* time measurement for coarse load balancing */
#endif
    
#define PRINT_SEPDVDL(s,v,dvdl) if (bSepDVDL) fprintf(fplog,sepdvdlformat,s,v,dvdl);
    
    GMX_MPE_LOG(ev_force_start);
    set_pbc(&pbc,fr->ePBC,box);
    
    /* Reset box */
    for(i=0; (i<DIM); i++)
    {
        box_size[i]=box[i][i];
    }
    
    bSepDVDL=(fr->bSepDVDL && do_per_step(step,ir->nstlog));
    debug_gmx();
    
    /* do QMMM first if requested */
    if(fr->bQMMM)
    {
        enerd->term[F_EQM] = calculate_QMMM(cr,x,f,fr,md);
    }
    
    if (bSepDVDL)
    {
        fprintf(fplog,"Step %s: non-bonded V and dVdl for node %d:\n",
                gmx_step_str(step,buf),cr->nodeid);
    }
    
    /* Call the short range functions all in one go. */
    GMX_MPE_LOG(ev_do_fnbf_start);
    
    dvdlambda = 0;
    
#ifdef GMX_MPI
    /*#define TAKETIME ((cr->npmenodes) && (fr->timesteps < 12))*/
#define TAKETIME FALSE
    if (TAKETIME)
    {
        MPI_Barrier(cr->mpi_comm_mygroup);
        t0=MPI_Wtime();
    }
#endif
    
    if (ir->nwall)
    {
        dvdlambda = do_walls(ir,fr,box,md,x,f,lambda,
                             enerd->grpp.ener[egLJSR],nrnb);
        PRINT_SEPDVDL("Walls",0.0,dvdlambda);
        enerd->dvdl_lin += dvdlambda;
    }
		
	/* If doing GB, reset dvda and calculate the Born radii */
	if (ir->implicit_solvent)
	{
		/* wallcycle_start(wcycle,ewcGB); */
		
		for(i=0;i<born->nr;i++)
		{
			fr->dvda[i]=0;
		}
		
		if(bBornRadii)
		{
			calc_gb_rad(cr,fr,ir,top,atype,x,&(fr->gblist),born,md,nrnb);
		}
		
		/* wallcycle_stop(wcycle, ewcGB); */
	}
	
    where();
    donb_flags = 0;
    if (flags & GMX_FORCE_FORCES)
    {
        donb_flags |= GMX_DONB_FORCES;
    }
    do_nonbonded(cr,fr,x,f,md,excl,
                 fr->bBHAM ?
                 enerd->grpp.ener[egBHAMSR] :
                 enerd->grpp.ener[egLJSR],
                 enerd->grpp.ener[egCOULSR],
				 enerd->grpp.ener[egGB],box_size,nrnb,
                 lambda,&dvdlambda,-1,-1,donb_flags);
    /* If we do foreign lambda and we have soft-core interactions
     * we have to recalculate the (non-linear) energies contributions.
     */
    if (ir->n_flambda > 0 && (flags & GMX_FORCE_DHDL) && ir->sc_alpha != 0)
    {
        init_enerdata(mtop->groups.grps[egcENER].nr,ir->n_flambda,&ed_lam);
        
        for(i=0; i<enerd->n_lambda; i++)
        {
            lam_i = (i==0 ? lambda : ir->flambda[i-1]);
            dvdl_dum = 0;
            reset_enerdata(&ir->opts,fr,TRUE,&ed_lam,FALSE);
            do_nonbonded(cr,fr,x,f,md,excl,
                         fr->bBHAM ?
                         ed_lam.grpp.ener[egBHAMSR] :
                         ed_lam.grpp.ener[egLJSR],
                         ed_lam.grpp.ener[egCOULSR],
                         enerd->grpp.ener[egGB], box_size,nrnb,
                         lam_i,&dvdl_dum,-1,-1,
                         GMX_DONB_FOREIGNLAMBDA);
            sum_epot(&ir->opts,&ed_lam);
            enerd->enerpart_lambda[i] += ed_lam.term[F_EPOT];
        }
        destroy_enerdata(&ed_lam);
    }
    where();
	
	/* If we are doing GB, calculate bonded forces and apply corrections 
	 * to the solvation forces */
	if (ir->implicit_solvent)  {
		calc_gb_forces(cr,md,born,top,atype,x,f,fr,idef,
                       ir->gb_algorithm,ir->sa_algorithm,nrnb,bBornRadii,&pbc,graph,enerd);
    }

#ifdef GMX_MPI
    if (TAKETIME)
    {
        t1=MPI_Wtime();
        fr->t_fnbf += t1-t0;
    }
#endif
    
    if (ir->sc_alpha != 0)
    {
        enerd->dvdl_nonlin += dvdlambda;
    }
    else
    {
        enerd->dvdl_lin    += dvdlambda;
    }
    Vsr = 0;
    if (bSepDVDL)
    {
        for(i=0; i<enerd->grpp.nener; i++)
        {
            Vsr +=
                (fr->bBHAM ?
                 enerd->grpp.ener[egBHAMSR][i] :
                 enerd->grpp.ener[egLJSR][i])
                + enerd->grpp.ener[egCOULSR][i] + enerd->grpp.ener[egGB][i];
        }
    }
    PRINT_SEPDVDL("VdW and Coulomb SR particle-p.",Vsr,dvdlambda);
    debug_gmx();
    
    GMX_MPE_LOG(ev_do_fnbf_finish);
    
    if (debug)
    {
        pr_rvecs(debug,0,"fshift after SR",fr->fshift,SHIFTS);
    }
    
    /* Shift the coordinates. Must be done before bonded forces and PPPM, 
     * but is also necessary for SHAKE and update, therefore it can NOT 
     * go when no bonded forces have to be evaluated.
     */
    
    /* Here sometimes we would not need to shift with NBFonly,
     * but we do so anyhow for consistency of the returned coordinates.
     */
    if (graph)
    {
        shift_self(graph,box,x);
        if (TRICLINIC(box))
        {
            inc_nrnb(nrnb,eNR_SHIFTX,2*graph->nnodes);
        }
        else
        {
            inc_nrnb(nrnb,eNR_SHIFTX,graph->nnodes);
        }
    }
    /* Check whether we need to do bondeds or correct for exclusions */
    if (fr->bMolPBC &&
        ((flags & GMX_FORCE_BONDED)
         || EEL_RF(fr->eeltype) || EEL_FULL(fr->eeltype)))
    {
        /* Since all atoms are in the rectangular or triclinic unit-cell,
         * only single box vector shifts (2 in x) are required.
         */
        set_pbc_dd(&pbc,fr->ePBC,cr->dd,TRUE,box);
    }
    debug_gmx();
    
    if (flags & GMX_FORCE_BONDED)
    {
        GMX_MPE_LOG(ev_calc_bonds_start);
        calc_bonds(fplog,cr->ms,
                   idef,x,hist,f,fr,&pbc,graph,enerd,nrnb,lambda,md,fcd,
                   DOMAINDECOMP(cr) ? cr->dd->gatindex : NULL, atype, born,
                   fr->bSepDVDL && do_per_step(step,ir->nstlog),step);
        
        /* Check if we have to determine energy differences
         * at foreign lambda's.
         */
        if (ir->n_flambda > 0 && (flags & GMX_FORCE_DHDL) &&
            idef->ilsort != ilsortNO_FE)
        {
            if (idef->ilsort != ilsortFE_SORTED)
            {
                gmx_incons("The bonded interactions are not sorted for free energy");
            }
            init_enerdata(mtop->groups.grps[egcENER].nr,ir->n_flambda,&ed_lam);
            
            for(i=0; i<enerd->n_lambda; i++)
            {
                lam_i = (i==0 ? lambda : ir->flambda[i-1]);
                dvdl_dum = 0;
                reset_enerdata(&ir->opts,fr,TRUE,&ed_lam,FALSE);
                calc_bonds_lambda(fplog,
                                  idef,x,fr,&pbc,graph,&ed_lam,nrnb,lam_i,md,
                                  fcd,
                                  DOMAINDECOMP(cr) ? cr->dd->gatindex : NULL);
                sum_epot(&ir->opts,&ed_lam);
                enerd->enerpart_lambda[i] += ed_lam.term[F_EPOT];
            }
            destroy_enerdata(&ed_lam);
        }
        debug_gmx();
        GMX_MPE_LOG(ev_calc_bonds_finish);
    }

    where();

    *cycles_pme = 0;
    if (EEL_FULL(fr->eeltype))
    {
        bSB = (ir->nwall == 2);
        if (bSB)
        {
            copy_mat(box,boxs);
            svmul(ir->wall_ewald_zfac,boxs[ZZ],boxs[ZZ]);
            box_size[ZZ] *= ir->wall_ewald_zfac;
        }
        
        clear_mat(fr->vir_el_recip);	
        
        if (fr->bEwald)
        {
            if (fr->n_tpi == 0)
            {
                dvdlambda = 0;
                Vcorr = ewald_LRcorrection(fplog,md->start,md->start+md->homenr,
                                           cr,fr,
                                           md->chargeA,
                                           md->nChargePerturbed ? md->chargeB : NULL,
                                           excl,x,bSB ? boxs : box,mu_tot,
                                           ir->ewald_geometry,
                                           ir->epsilon_surface,
                                           lambda,&dvdlambda,&vdip,&vcharge);
                PRINT_SEPDVDL("Ewald excl./charge/dip. corr.",Vcorr,dvdlambda);
                enerd->dvdl_lin += dvdlambda;
            }
            else
            {
                if (ir->ewald_geometry != eewg3D || ir->epsilon_surface != 0)
                {
                    gmx_fatal(FARGS,"TPI with PME currently only works in a 3D geometry with tin-foil boundary conditions");
                }
                /* The TPI molecule does not have exclusions with the rest
                 * of the system and no intra-molecular PME grid contributions
                 * will be calculated in gmx_pme_calc_energy.
                 */
                Vcorr = 0;
            }
        }
        else
        {
            Vcorr = shift_LRcorrection(fplog,md->start,md->homenr,cr,fr,
                                       md->chargeA,excl,x,TRUE,box,
                                       fr->vir_el_recip);
        }
        
        dvdlambda = 0;
        status = 0;
        switch (fr->eeltype)
        {
        case eelPPPM:
            status = gmx_pppm_do(fplog,fr->pmedata,FALSE,x,fr->f_novirsum,
                                 md->chargeA,
                                 box_size,fr->phi,cr,md->start,md->homenr,
                                 nrnb,ir->pme_order,&Vlr);
            break;
        case eelPME:
        case eelPMESWITCH:
        case eelPMEUSER:
        case eelPMEUSERSWITCH:
            if (cr->duty & DUTY_PME)
            {
                if (fr->n_tpi == 0 || (flags & GMX_FORCE_STATECHANGED))
                {
                    pme_flags = GMX_PME_SPREAD_Q | GMX_PME_SOLVE;
                    if (flags & GMX_FORCE_FORCES)
                    {
                        pme_flags |= GMX_PME_CALC_F;
                    }
                    if (flags & GMX_FORCE_VIRIAL)
                    {
                        pme_flags |= GMX_PME_CALC_ENER_VIR;
                    }
                    if (fr->n_tpi > 0)
                    {
                        /* We don't calculate f, but we do want the potential */
                        pme_flags |= GMX_PME_CALC_POT;
                    }
                    wallcycle_start(wcycle,ewcPMEMESH);
                    status = gmx_pme_do(fr->pmedata,
                                        md->start,md->homenr - fr->n_tpi,
                                        x,fr->f_novirsum,
                                        md->chargeA,md->chargeB,
                                        bSB ? boxs : box,cr,
                                        DOMAINDECOMP(cr) ? dd_pme_maxshift_x(cr->dd) : 0,
                                        DOMAINDECOMP(cr) ? dd_pme_maxshift_y(cr->dd) : 0,
                                        nrnb,wcycle,
                                        fr->vir_el_recip,fr->ewaldcoeff,
                                        &Vlr,lambda,&dvdlambda,
                                        pme_flags);
                    *cycles_pme = wallcycle_stop(wcycle,ewcPMEMESH);

                    /* We should try to do as little computation after
                     * this as possible, because parallel PME synchronizes
                     * the nodes, so we want all load imbalance of the rest
                     * of the force calculation to be before the PME call.
                     * DD load balancing is done on the whole time of
                     * the force call (without PME).
                     */
                }
                if (fr->n_tpi > 0)
                {
                    /* Determine the PME grid energy of the test molecule
                     * with the PME grid potential of the other charges.
                     */
                    gmx_pme_calc_energy(fr->pmedata,fr->n_tpi,
                                        x + md->homenr - fr->n_tpi,
                                        md->chargeA + md->homenr - fr->n_tpi,
                                        &Vlr);
                }
                PRINT_SEPDVDL("PME mesh",Vlr,dvdlambda);
            } 
            else
            {
                /* Energies and virial are obtained later from the PME nodes */
                /* but values have to be zeroed out here */
                Vlr=0.0;
            }
            break;
        case eelEWALD:
            Vlr = do_ewald(fplog,FALSE,ir,x,fr->f_novirsum,
                           md->chargeA,md->chargeB,
                           box_size,cr,md->homenr,
                           fr->vir_el_recip,fr->ewaldcoeff,
                           lambda,&dvdlambda,fr->ewald_table);
            PRINT_SEPDVDL("Ewald long-range",Vlr,dvdlambda);
            break;
        default:
            Vlr = 0;
            gmx_fatal(FARGS,"No such electrostatics method implemented %s",
                      eel_names[fr->eeltype]);
        }
        if (status != 0)
        {
            gmx_fatal(FARGS,"Error %d in long range electrostatics routine %s",
                      status,EELTYPE(fr->eeltype));
		}
        enerd->dvdl_lin += dvdlambda;
        enerd->term[F_COUL_RECIP] = Vlr + Vcorr;
        if (debug)
        {
            fprintf(debug,"Vlr = %g, Vcorr = %g, Vlr_corr = %g\n",
                    Vlr,Vcorr,enerd->term[F_COUL_RECIP]);
            pr_rvecs(debug,0,"vir_el_recip after corr",fr->vir_el_recip,DIM);
            pr_rvecs(debug,0,"fshift after LR Corrections",fr->fshift,SHIFTS);
        }
    }
    else
    {
        if (EEL_RF(fr->eeltype))
        {
            dvdlambda = 0;
            
            if (fr->eeltype != eelRF_NEC)
            {
                enerd->term[F_RF_EXCL] =
                    RF_excl_correction(fplog,fr,graph,md,excl,x,f,
                                       fr->fshift,&pbc,lambda,&dvdlambda);
            }
            
            enerd->dvdl_lin += dvdlambda;
            PRINT_SEPDVDL("RF exclusion correction",
                          enerd->term[F_RF_EXCL],dvdlambda);
        }
    }
    where();
    debug_gmx();
	
    if (debug)
    {
        print_nrnb(debug,nrnb); 
    }
    debug_gmx();
    
#ifdef GMX_MPI
    if (TAKETIME)
    {
        t2=MPI_Wtime();
        MPI_Barrier(cr->mpi_comm_mygroup);
        t3=MPI_Wtime();
        fr->t_wait += t3-t2;
        if (fr->timesteps == 11)
        {
            fprintf(stderr,"* PP load balancing info: node %d, step %s, rel wait time=%3.0f%% , load string value: %7.2f\n", 
                    cr->nodeid, gmx_step_str(fr->timesteps,buf), 
                    100*fr->t_wait/(fr->t_wait+fr->t_fnbf), 
                    (fr->t_fnbf+fr->t_wait)/fr->t_fnbf);
        }	  
        fr->timesteps++;
    }
#endif
    
    if (debug)
    {
        pr_rvecs(debug,0,"fshift after bondeds",fr->fshift,SHIFTS);
    }
    
    GMX_MPE_LOG(ev_force_finish);

}

void init_enerdata(int ngener,int n_flambda,gmx_enerdata_t *enerd)
{
    int i,n2;
    
    for(i=0; i<F_NRE; i++)
    {
        enerd->term[i] = 0;
    }
    
    n2=ngener*ngener;
    if (debug)
    {
        fprintf(debug,"Creating %d sized group matrix for energies\n",n2);
    }
    enerd->grpp.nener = n2;
    for(i=0; (i<egNR); i++)
    {
        snew(enerd->grpp.ener[i],n2);
    }

    if (n_flambda)
    {
        enerd->n_lambda = 1 + n_flambda;
        snew(enerd->enerpart_lambda,enerd->n_lambda);
    }
    else
    {
        enerd->n_lambda = 0;
    }
}

void destroy_enerdata(gmx_enerdata_t *enerd)
{
    int i;

    for(i=0; (i<egNR); i++)
    {
        sfree(enerd->grpp.ener[i]);
    }

    if (enerd->n_lambda)
    {
        sfree(enerd->enerpart_lambda);
    }
}

static real sum_v(int n,real v[])
{
  real t;
  int  i;
  
  t = 0.0;
  for(i=0; (i<n); i++)
    t = t + v[i];
    
  return t;
}

void sum_epot(t_grpopts *opts,gmx_enerdata_t *enerd)
{
  gmx_grppairener_t *grpp;
  real *epot;
  int i;
  
  grpp = &enerd->grpp;
  epot = enerd->term;

  /* Accumulate energies */
  epot[F_COUL_SR]  = sum_v(grpp->nener,grpp->ener[egCOULSR]);
  epot[F_LJ]       = sum_v(grpp->nener,grpp->ener[egLJSR]);
  epot[F_LJ14]     = sum_v(grpp->nener,grpp->ener[egLJ14]);
  epot[F_COUL14]   = sum_v(grpp->nener,grpp->ener[egCOUL14]);
  epot[F_COUL_LR]  = sum_v(grpp->nener,grpp->ener[egCOULLR]);
  epot[F_LJ_LR]    = sum_v(grpp->nener,grpp->ener[egLJLR]);
  /* We have already added 1-2,1-3, and 1-4 terms to F_GBPOL */
  epot[F_GBPOL]   += sum_v(grpp->nener,grpp->ener[egGB]);
    
/* lattice part of LR doesnt belong to any group
 * and has been added earlier
 */
  epot[F_BHAM]     = sum_v(grpp->nener,grpp->ener[egBHAMSR]);
  epot[F_BHAM_LR]  = sum_v(grpp->nener,grpp->ener[egBHAMLR]);

  epot[F_EPOT] = 0;
  for(i=0; (i<F_EPOT); i++)
    if (i != F_DISRESVIOL && i != F_ORIRESDEV && i != F_DIHRESVIOL)
      epot[F_EPOT] += epot[i];
}

void sum_dhdl(gmx_enerdata_t *enerd,double lambda,t_inputrec *ir)
{
    int i;
    double dlam,dhdl_lin;

    enerd->term[F_DVDL] = enerd->dvdl_lin + enerd->dvdl_nonlin;

    if (debug)
    {
        fprintf(debug,"dvdl: %f, non-linear %f + linear %f\n",
                enerd->term[F_DVDL],enerd->dvdl_nonlin,enerd->dvdl_lin);
    }

    /* Notes on the foreign lambda free energy difference evaluation:
     * Adding the potential and ekin terms that depend linearly on lambda
     * as delta lam * dvdl to the energy differences is exact.
     * For the constraint dvdl this is not exact, but we have no other option.
     * For the non-bonded LR term we assume that the soft-core (if present)
     * no longer affects the energy beyond the short-range cut-off,
     * which is a very good approximation (except for exotic settings).
     */
    for(i=1; i<enerd->n_lambda; i++)
    {
        dlam = (ir->flambda[i-1] - lambda);
        dhdl_lin =
            enerd->dvdl_lin + enerd->term[F_DKDL] + enerd->term[F_DHDL_CON];
        if (debug)
        {
            fprintf(debug,"enerdiff lam %g: non-linear %f linear %f*%f\n",
                    ir->flambda[i-1],
                    enerd->enerpart_lambda[i] - enerd->enerpart_lambda[0],
                    dlam,dhdl_lin);
        }
        enerd->enerpart_lambda[i] += dlam*dhdl_lin;

    }
}

void reset_enerdata(t_grpopts *opts,
                    t_forcerec *fr,gmx_bool bNS,
                    gmx_enerdata_t *enerd,
                    gmx_bool bMaster)
{
  gmx_bool bKeepLR;
  int  i,j;
  
  /* First reset all energy components, except for the long range terms
   * on the master at non neighbor search steps, since the long range
   * terms have already been summed at the last neighbor search step.
   */
  bKeepLR = (fr->bTwinRange && !bNS);
  for(i=0; (i<egNR); i++) {
    if (!(bKeepLR && bMaster && (i == egCOULLR || i == egLJLR))) {
      for(j=0; (j<enerd->grpp.nener); j++)
	enerd->grpp.ener[i][j] = 0.0;
    }
  }
  enerd->dvdl_lin    = 0.0;
  enerd->dvdl_nonlin = 0.0;

  /* Normal potential energy components */
  for(i=0; (i<=F_EPOT); i++) {
    enerd->term[i] = 0.0;
  }
  /* Initialize the dVdlambda term with the long range contribution */
  enerd->term[F_DVDL]     = 0.0;
  enerd->term[F_DKDL]     = 0.0;
  enerd->term[F_DHDL_CON] = 0.0;
  if (enerd->n_lambda > 0)
  {
      for(i=0; i<enerd->n_lambda; i++)
      {
          enerd->enerpart_lambda[i] = 0.0;
      }
  }
}
