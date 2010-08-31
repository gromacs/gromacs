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

#include "confio.h"
#include "constr.h"
#include "copyrite.h"
#include "invblock.h"
#include "main.h"
#include "mdrun.h"
#include "nrnb.h"
#include "smalloc.h"
#include "vec.h"
#include "physics.h"
#include "names.h"
#include "txtdump.h"
#include "domdec.h"
#include "pdbio.h"
#include "partdec.h"
#include "splitter.h"
#include "mtop_util.h"
#include "gmxfio.h"

typedef struct gmx_constr {
  int              ncon_tot;     /* The total number of constraints    */
  int              nflexcon;     /* The number of flexible constraints */
  int              n_at2con_mt;  /* The size of at2con = #moltypes     */
  t_blocka         *at2con_mt;   /* A list of atoms to constraints     */
  gmx_lincsdata_t  lincsd;       /* LINCS data                         */
  gmx_shakedata_t  shaked;       /* SHAKE data                         */
  gmx_settledata_t settled;      /* SETTLE data                        */
  int              nblocks;      /* The number of SHAKE blocks         */
  int              *sblock;      /* The SHAKE blocks                   */
  int              sblock_nalloc;/* The allocation size of sblock      */
  real             *lagr;        /* Lagrange multipliers for SHAKE     */
  int              lagr_nalloc;  /* The allocation size of lagr        */
  int              maxwarn;      /* The maximum number of warnings     */
  int              warncount_lincs;
  int              warncount_settle;
  gmx_edsam_t      ed;           /* The essential dynamics data        */

  gmx_mtop_t       *warn_mtop;   /* Only used for printing warnings    */
} t_gmx_constr;

typedef struct {
  atom_id iatom[3];
  atom_id blocknr;
} t_sortblock;

static t_vetavars *init_vetavars(real veta,real vetanew, t_inputrec *ir, gmx_ekindata_t *ekind, gmx_bool bPscal) 
{
    t_vetavars *vars;
    double g;
    int i;

    snew(vars,1);
    snew(vars->vscale_nhc,ir->opts.ngtc);
    /* first, set the alpha integrator variable */
    if ((ir->opts.nrdf[0] > 0) && bPscal) 
    {
        vars->alpha = 1.0 + DIM/((double)ir->opts.nrdf[0]);  
    } else {
        vars->alpha = 1.0;
    }
    g = 0.5*veta*ir->delta_t;
    vars->rscale = exp(g)*series_sinhx(g);
    g = -0.25*vars->alpha*veta*ir->delta_t;
    vars->vscale = exp(g)*series_sinhx(g);
    vars->rvscale = vars->vscale*vars->rscale;
    vars->veta = vetanew;
    if ((ekind==NULL) || (!bPscal))
    {
        for (i=0;i<ir->opts.ngtc;i++)
        {
            vars->vscale_nhc[i] = 1;
        }
    } else {
        for (i=0;i<ir->opts.ngtc;i++)
        {
            vars->vscale_nhc[i] = ekind->tcstat[i].vscale_nhc;
        }
    }
    return vars;
}

static void free_vetavars(t_vetavars *vars) 
{
    sfree(vars->vscale_nhc);
    sfree(vars);
}

static int pcomp(const void *p1, const void *p2)
{
  int     db;
  atom_id min1,min2,max1,max2;
  t_sortblock *a1=(t_sortblock *)p1;
  t_sortblock *a2=(t_sortblock *)p2;
  
  db=a1->blocknr-a2->blocknr;
  
  if (db != 0)
    return db;
    
  min1=min(a1->iatom[1],a1->iatom[2]);
  max1=max(a1->iatom[1],a1->iatom[2]);
  min2=min(a2->iatom[1],a2->iatom[2]);
  max2=max(a2->iatom[1],a2->iatom[2]);
  
  if (min1 == min2)
    return max1-max2;
  else
    return min1-min2;
}

static int icomp(const void *p1, const void *p2)
{
  atom_id *a1=(atom_id *)p1;
  atom_id *a2=(atom_id *)p2;

  return (*a1)-(*a2);
}

int n_flexible_constraints(struct gmx_constr *constr)
{
  int nflexcon;

  if (constr)
    nflexcon = constr->nflexcon;
  else
    nflexcon = 0;

  return nflexcon;
}

void too_many_constraint_warnings(int eConstrAlg,int warncount)
{
  const char *abort="- aborting to avoid logfile runaway.\n"
    "This normally happens when your system is not sufficiently equilibrated,"
    "or if you are changing lambda too fast in free energy simulations.\n";
  
  gmx_fatal(FARGS,
	    "Too many %s warnings (%d)\n"
	    "If you know what you are doing you can %s"
	    "set the environment variable GMX_MAXCONSTRWARN to -1,\n"
	    "but normally it is better to fix the problem",
	    (eConstrAlg == econtLINCS) ? "LINCS" : "SETTLE",warncount,
	    (eConstrAlg == econtLINCS) ?
	    "adjust the lincs warning threshold in your mdp file\nor " : "\n");
}

static void write_constr_pdb(const char *fn,const char *title,
                             gmx_mtop_t *mtop,
                             int start,int homenr,t_commrec *cr,
                             rvec x[],matrix box)
{
    char fname[STRLEN],format[STRLEN];
    FILE *out;
    int  dd_ac0=0,dd_ac1=0,i,ii,resnr;
    gmx_domdec_t *dd;
    char *anm,*resnm;
  
    dd = NULL;
    if (PAR(cr))
    {
        sprintf(fname,"%s_n%d.pdb",fn,cr->sim_nodeid);
        if (DOMAINDECOMP(cr))
        {
            dd = cr->dd;
            dd_get_constraint_range(dd,&dd_ac0,&dd_ac1);
            start = 0;
            homenr = dd_ac1;
        }
    }
    else
    {
        sprintf(fname,"%s.pdb",fn);
    }
    sprintf(format,"%s\n",pdbformat);
    
    out = gmx_fio_fopen(fname,"w");
    
    fprintf(out,"TITLE     %s\n",title);
    gmx_write_pdb_box(out,-1,box);
    for(i=start; i<start+homenr; i++)
    {
        if (dd != NULL)
        {
            if (i >= dd->nat_home && i < dd_ac0)
            {
                continue;
            }
            ii = dd->gatindex[i];
        }
        else
        {
            ii = i;
        }
        gmx_mtop_atominfo_global(mtop,ii,&anm,&resnr,&resnm);
        fprintf(out,format,"ATOM",(ii+1)%100000,
                anm,resnm,' ',resnr%10000,' ',
                10*x[i][XX],10*x[i][YY],10*x[i][ZZ]);
    }
    fprintf(out,"TER\n");

    gmx_fio_fclose(out);
}
			     
static void dump_confs(FILE *fplog,gmx_large_int_t step,gmx_mtop_t *mtop,
		       int start,int homenr,t_commrec *cr,
		       rvec x[],rvec xprime[],matrix box)
{
  char buf[256],buf2[22];
 
  char *env=getenv("GMX_SUPPRESS_DUMP");
  if (env)
      return; 
  
  sprintf(buf,"step%sb",gmx_step_str(step,buf2));
  write_constr_pdb(buf,"initial coordinates",
		   mtop,start,homenr,cr,x,box);
  sprintf(buf,"step%sc",gmx_step_str(step,buf2));
  write_constr_pdb(buf,"coordinates after constraining",
		   mtop,start,homenr,cr,xprime,box);
  if (fplog)
  {
      fprintf(fplog,"Wrote pdb files with previous and current coordinates\n");
  }
  fprintf(stderr,"Wrote pdb files with previous and current coordinates\n");
}

static void pr_sortblock(FILE *fp,const char *title,int nsb,t_sortblock sb[])
{
  int i;
  
  fprintf(fp,"%s\n",title);
  for(i=0; (i<nsb); i++)
    fprintf(fp,"i: %5d, iatom: (%5d %5d %5d), blocknr: %5d\n",
	    i,sb[i].iatom[0],sb[i].iatom[1],sb[i].iatom[2],
	    sb[i].blocknr);
}

gmx_bool constrain(FILE *fplog,gmx_bool bLog,gmx_bool bEner,
               struct gmx_constr *constr,
               t_idef *idef,t_inputrec *ir,gmx_ekindata_t *ekind,
               t_commrec *cr,
               gmx_large_int_t step,int delta_step,
               t_mdatoms *md,
               rvec *x,rvec *xprime,rvec *min_proj,matrix box,
               real lambda,real *dvdlambda,
               rvec *v,tensor *vir,
               t_nrnb *nrnb,int econq,gmx_bool bPscal,real veta, real vetanew)
{
    gmx_bool    bOK,bDump;
    int     start,homenr,nrend;
    int     i,j,d;
    int     ncons,error;
    tensor  rmdr;
    rvec    *vstor;
    real    invdt,vir_fac,t;
    t_ilist *settle;
    int     nsettle;
    t_pbc   pbc;
    char    buf[22];
    t_vetavars *vetavar;

    if (econq == econqForceDispl && !EI_ENERGY_MINIMIZATION(ir->eI))
    {
        gmx_incons("constrain called for forces displacements while not doing energy minimization, can not do this while the LINCS and SETTLE constraint connection matrices are mass weighted");
    }
    
    bOK   = TRUE;
    bDump = FALSE;
    
    start  = md->start;
    homenr = md->homenr;
    nrend = start+homenr;

    /* set constants for pressure control integration */ 
    vetavar = init_vetavars(veta,vetanew,ir,ekind,bPscal);

    if (ir->delta_t == 0)
    {
        invdt = 0;
    }
    else
    {
        invdt  = 1/ir->delta_t;
    }

    if (ir->efep != efepNO && EI_DYNAMICS(ir->eI))
    {
        /* Set the constraint lengths for the step at which this configuration
         * is meant to be. The invmasses should not be changed.
         */
        lambda += delta_step*ir->delta_lambda;
    }
    
    if (vir != NULL)
    {
        clear_mat(rmdr);
    }
    
    where();
    if (constr->lincsd)
    {
        bOK = constrain_lincs(fplog,bLog,bEner,ir,step,constr->lincsd,md,cr,
                              x,xprime,min_proj,box,lambda,dvdlambda,
                              invdt,v,vir!=NULL,rmdr,
                              econq,nrnb,
                              constr->maxwarn,&constr->warncount_lincs);
        if (!bOK && constr->maxwarn >= 0)
        {
            if (fplog != NULL)
            {
                fprintf(fplog,"Constraint error in algorithm %s at step %s\n",
                        econstr_names[econtLINCS],gmx_step_str(step,buf));
            }
            bDump = TRUE;
        }
    }	
    
    if (constr->nblocks > 0)
    {
        switch (econq) {
        case (econqCoord):
            bOK = bshakef(fplog,constr->shaked,
                          homenr,md->invmass,constr->nblocks,constr->sblock,
                          idef,ir,box,x,xprime,nrnb,
                          constr->lagr,lambda,dvdlambda,
                          invdt,v,vir!=NULL,rmdr,constr->maxwarn>=0,econq,vetavar);
            break;
        case (econqVeloc):
            bOK = bshakef(fplog,constr->shaked,
                          homenr,md->invmass,constr->nblocks,constr->sblock,
                          idef,ir,box,x,min_proj,nrnb,
                          constr->lagr,lambda,dvdlambda,
                          invdt,NULL,vir!=NULL,rmdr,constr->maxwarn>=0,econq,vetavar);
            break;
        default:
            gmx_fatal(FARGS,"Internal error, SHAKE called for constraining something else than coordinates");
            break;
        }

        if (!bOK && constr->maxwarn >= 0)
        {
            if (fplog != NULL)
            {
                fprintf(fplog,"Constraint error in algorithm %s at step %s\n",
                        econstr_names[econtSHAKE],gmx_step_str(step,buf));
            }
            bDump = TRUE;
        }
    }
        
    settle  = &idef->il[F_SETTLE];
    if (settle->nr > 0)
    {
        nsettle = settle->nr/2;
        
        switch (econq)
        {
        case econqCoord:
            csettle(constr->settled,
                    nsettle,settle->iatoms,x[0],xprime[0],
                    invdt,v[0],vir!=NULL,rmdr,&error,vetavar);
            inc_nrnb(nrnb,eNR_SETTLE,nsettle);
            if (v != NULL)
            {
                inc_nrnb(nrnb,eNR_CONSTR_V,nsettle*3);
            }
            if (vir != NULL)
            {
                inc_nrnb(nrnb,eNR_CONSTR_VIR,nsettle*3);
            }
            
            bOK = (error < 0);
            if (!bOK && constr->maxwarn >= 0)
            {
                char buf[256];
                sprintf(buf,
                        "\nstep " gmx_large_int_pfmt ": Water molecule starting at atom %d can not be "
                        "settled.\nCheck for bad contacts and/or reduce the timestep if appropriate.\n",
                        step,ddglatnr(cr->dd,settle->iatoms[error*2+1]));
                if (fplog)
                {
                    fprintf(fplog,"%s",buf);
                }
                fprintf(stderr,"%s",buf);
                constr->warncount_settle++;
                if (constr->warncount_settle > constr->maxwarn)
                {
                    too_many_constraint_warnings(-1,constr->warncount_settle);
                }
                bDump = TRUE;
                break;
            case econqVeloc:
            case econqDeriv:
            case econqForce:
            case econqForceDispl:
                settle_proj(fplog,constr->settled,econq,
                            nsettle,settle->iatoms,x,
                            xprime,min_proj,vir!=NULL,rmdr,vetavar);
                /* This is an overestimate */
                inc_nrnb(nrnb,eNR_SETTLE,nsettle);
                break;
            case econqDeriv_FlexCon:
                /* Nothing to do, since the are no flexible constraints in settles */
                break;
            default:
                gmx_incons("Unknown constraint quantity for settle");
            }
        }
    }

    free_vetavars(vetavar);
    
    if (vir != NULL)
    {
        switch (econq)
        {
        case econqCoord:
            vir_fac = 0.5/(ir->delta_t*ir->delta_t);
            break;
        case econqVeloc:
            vir_fac = 0.5/ir->delta_t;
            break;
        case econqForce:
        case econqForceDispl:
            vir_fac = 0.5;
            break;
        default:
            vir_fac = 0;
            gmx_incons("Unsupported constraint quantity for virial");
        }
        
        if (EI_VV(ir->eI))
        {
            vir_fac *= 2;  /* only constraining over half the distance here */
        }
        for(i=0; i<DIM; i++)
        {
            for(j=0; j<DIM; j++)
            {
                (*vir)[i][j] = vir_fac*rmdr[i][j];
            }
        }
    }
    
    if (bDump)
    {
        dump_confs(fplog,step,constr->warn_mtop,start,homenr,cr,x,xprime,box);
    }
    
    if (econq == econqCoord)
    {
        if (ir->ePull == epullCONSTRAINT)
        {
            if (EI_DYNAMICS(ir->eI))
            {
                t = ir->init_t + (step + delta_step)*ir->delta_t;
            }
            else
            {
                t = ir->init_t;
            }
            set_pbc(&pbc,ir->ePBC,box);
            pull_constraint(ir->pull,md,&pbc,cr,ir->delta_t,t,x,xprime,v,*vir);
        }
        if (constr->ed && delta_step > 0)
        {
            /* apply the essential dynamcs constraints here */
            do_edsam(ir,step,md,cr,xprime,v,box,constr->ed);
        }
    }
    
    return bOK;
}

real *constr_rmsd_data(struct gmx_constr *constr)
{
  if (constr->lincsd)
    return lincs_rmsd_data(constr->lincsd);
  else
    return NULL;
}

real constr_rmsd(struct gmx_constr *constr,gmx_bool bSD2)
{
  if (constr->lincsd)
    return lincs_rmsd(constr->lincsd,bSD2);
  else
    return 0;
}

static void make_shake_sblock_pd(struct gmx_constr *constr,
				 t_idef *idef,t_mdatoms *md)
{
  int  i,j,m,ncons;
  int  bstart,bnr;
  t_blocka    sblocks;
  t_sortblock *sb;
  t_iatom     *iatom;
  atom_id     *inv_sblock;

  /* Since we are processing the local topology,
   * the F_CONSTRNC ilist has been concatenated to the F_CONSTR ilist.
   */
  ncons = idef->il[F_CONSTR].nr/3;

  init_blocka(&sblocks);
  gen_sblocks(NULL,md->start,md->start+md->homenr,idef,&sblocks,FALSE);
  
  /*
    bstart=(idef->nodeid > 0) ? blocks->multinr[idef->nodeid-1] : 0;
    nblocks=blocks->multinr[idef->nodeid] - bstart;
  */
  bstart  = 0;
  constr->nblocks = sblocks.nr;
  if (debug) 
    fprintf(debug,"ncons: %d, bstart: %d, nblocks: %d\n",
	    ncons,bstart,constr->nblocks);
  
  /* Calculate block number for each atom */
  inv_sblock = make_invblocka(&sblocks,md->nr);
  
  done_blocka(&sblocks);
  
  /* Store the block number in temp array and
   * sort the constraints in order of the sblock number 
   * and the atom numbers, really sorting a segment of the array!
   */
#ifdef DEBUGIDEF 
  pr_idef(fplog,0,"Before Sort",idef);
#endif
  iatom=idef->il[F_CONSTR].iatoms;
  snew(sb,ncons);
  for(i=0; (i<ncons); i++,iatom+=3) {
    for(m=0; (m<3); m++)
      sb[i].iatom[m] = iatom[m];
    sb[i].blocknr = inv_sblock[iatom[1]];
  }
  
  /* Now sort the blocks */
  if (debug) {
    pr_sortblock(debug,"Before sorting",ncons,sb);
    fprintf(debug,"Going to sort constraints\n");
  }
  
  qsort(sb,ncons,(size_t)sizeof(*sb),pcomp);
  
  if (debug) {
    pr_sortblock(debug,"After sorting",ncons,sb);
  }
  
  iatom=idef->il[F_CONSTR].iatoms;
  for(i=0; (i<ncons); i++,iatom+=3) 
    for(m=0; (m<3); m++)
      iatom[m]=sb[i].iatom[m];
#ifdef DEBUGIDEF
  pr_idef(fplog,0,"After Sort",idef);
#endif
  
  j=0;
  snew(constr->sblock,constr->nblocks+1);
  bnr=-2;
  for(i=0; (i<ncons); i++) {
    if (sb[i].blocknr != bnr) {
      bnr=sb[i].blocknr;
      constr->sblock[j++]=3*i;
    }
  }
  /* Last block... */
  constr->sblock[j++] = 3*ncons;
  
  if (j != (constr->nblocks+1)) {
    fprintf(stderr,"bstart: %d\n",bstart);
    fprintf(stderr,"j: %d, nblocks: %d, ncons: %d\n",
	    j,constr->nblocks,ncons);
    for(i=0; (i<ncons); i++)
      fprintf(stderr,"i: %5d  sb[i].blocknr: %5u\n",i,sb[i].blocknr);
    for(j=0; (j<=constr->nblocks); j++)
      fprintf(stderr,"sblock[%3d]=%5d\n",j,(int)constr->sblock[j]);
    gmx_fatal(FARGS,"DEATH HORROR: "
	      "sblocks does not match idef->il[F_CONSTR]");
  }
  sfree(sb);
  sfree(inv_sblock);
}

static void make_shake_sblock_dd(struct gmx_constr *constr,
				 t_ilist *ilcon,t_block *cgs,
				 gmx_domdec_t *dd)
{
  int ncons,c,cg;
  t_iatom *iatom;

  if (dd->ncg_home+1 > constr->sblock_nalloc) {
    constr->sblock_nalloc = over_alloc_dd(dd->ncg_home+1);
    srenew(constr->sblock,constr->sblock_nalloc);
  }
  
  ncons = ilcon->nr/3;
  iatom = ilcon->iatoms;
  constr->nblocks = 0;
  cg = 0;
  for(c=0; c<ncons; c++) {
    if (c == 0 || iatom[1] >= cgs->index[cg+1]) {
      constr->sblock[constr->nblocks++] = 3*c;
      while (iatom[1] >= cgs->index[cg+1])
	cg++;
    }
    iatom += 3;
  }
  constr->sblock[constr->nblocks] = 3*ncons;
}

t_blocka make_at2con(int start,int natoms,
		     t_ilist *ilist,t_iparams *iparams,
		     gmx_bool bDynamics,int *nflexiblecons)
{
  int *count,ncon,con,con_tot,nflexcon,ftype,i,a;
  t_iatom  *ia;
  t_blocka at2con;
  gmx_bool bFlexCon;
  
  snew(count,natoms);
  nflexcon = 0;
  for(ftype=F_CONSTR; ftype<=F_CONSTRNC; ftype++) {
    ncon = ilist[ftype].nr/3;
    ia   = ilist[ftype].iatoms;
    for(con=0; con<ncon; con++) {
      bFlexCon = (iparams[ia[0]].constr.dA == 0 &&
		  iparams[ia[0]].constr.dB == 0);
      if (bFlexCon)
	nflexcon++;
      if (bDynamics || !bFlexCon) {
	for(i=1; i<3; i++) {
	  a = ia[i] - start;
	  count[a]++;
	}
      }
      ia += 3;
    }
  }
  *nflexiblecons = nflexcon;

  at2con.nr = natoms;
  at2con.nalloc_index = at2con.nr+1;
  snew(at2con.index,at2con.nalloc_index);
  at2con.index[0] = 0;
  for(a=0; a<natoms; a++) {
    at2con.index[a+1] = at2con.index[a] + count[a];
    count[a] = 0;
  }
  at2con.nra = at2con.index[natoms];
  at2con.nalloc_a = at2con.nra;
  snew(at2con.a,at2con.nalloc_a);

  /* The F_CONSTRNC constraints have constraint numbers
   * that continue after the last F_CONSTR constraint.
   */
  con_tot = 0;
  for(ftype=F_CONSTR; ftype<=F_CONSTRNC; ftype++) {
    ncon = ilist[ftype].nr/3;
    ia   = ilist[ftype].iatoms;
    for(con=0; con<ncon; con++) {
      bFlexCon = (iparams[ia[0]].constr.dA == 0 &&
		  iparams[ia[0]].constr.dB == 0);
      if (bDynamics || !bFlexCon) {
	for(i=1; i<3; i++) {
	  a = ia[i] - start;
	  at2con.a[at2con.index[a]+count[a]++] = con_tot;
	}
      }
      con_tot++;
      ia += 3;
    }
  }
  
  sfree(count);

  return at2con;
}

void set_constraints(struct gmx_constr *constr,
                     gmx_localtop_t *top,t_inputrec *ir,
                     t_mdatoms *md,t_commrec *cr)
{
    t_idef *idef;
    int    ncons;
    t_ilist *settle;
    int    iO,iH;
    
    idef = &top->idef;
       
    if (constr->ncon_tot > 0)
    {
        /* We are using the local topology,
         * so there are only F_CONSTR constraints.
         */
        ncons = idef->il[F_CONSTR].nr/3;
        
        /* With DD we might also need to call LINCS with ncons=0 for
         * communicating coordinates to other nodes that do have constraints.
         */
        if (ir->eConstrAlg == econtLINCS)
        {
            set_lincs(idef,md,EI_DYNAMICS(ir->eI),cr,constr->lincsd);
        }
        if (ir->eConstrAlg == econtSHAKE)
        {
            if (cr->dd)
            {
                make_shake_sblock_dd(constr,&idef->il[F_CONSTR],&top->cgs,cr->dd);
            }
            else
            {
                make_shake_sblock_pd(constr,idef,md);
            }
            if (ncons > constr->lagr_nalloc)
            {
                constr->lagr_nalloc = over_alloc_dd(ncons);
                srenew(constr->lagr,constr->lagr_nalloc);
            }

            constr->shaked = shake_init();
        }
    }

    if (idef->il[F_SETTLE].nr > 0 && constr->settled == NULL)
    {
        settle = &idef->il[F_SETTLE];
        iO = settle->iatoms[1];
        iH = settle->iatoms[1]+1;
        constr->settled =
            settle_init(md->massT[iO],md->massT[iH],
                        md->invmass[iO],md->invmass[iH],
                        idef->iparams[settle->iatoms[0]].settle.doh,
                        idef->iparams[settle->iatoms[0]].settle.dhh);
    }
    
    /* Make a selection of the local atoms for essential dynamics */
    if (constr->ed && cr->dd)
    {
        dd_make_local_ed_indices(cr->dd,constr->ed);
    }
}

static void constr_recur(t_blocka *at2con,
			 t_ilist *ilist,t_iparams *iparams,gmx_bool bTopB,
			 int at,int depth,int nc,int *path,
			 real r0,real r1,real *r2max,
			 int *count)
{
  int  ncon1;
  t_iatom *ia1,*ia2;
  int  c,con,a1;
  gmx_bool bUse;
  t_iatom *ia;
  real len,rn0,rn1;

  (*count)++;

  ncon1 = ilist[F_CONSTR].nr/3;
  ia1   = ilist[F_CONSTR].iatoms;
  ia2   = ilist[F_CONSTRNC].iatoms;

  /* Loop over all constraints connected to this atom */
  for(c=at2con->index[at]; c<at2con->index[at+1]; c++) {
    con = at2con->a[c];
    /* Do not walk over already used constraints */
    bUse = TRUE;
    for(a1=0; a1<depth; a1++) {
      if (con == path[a1])
	bUse = FALSE;
    }
    if (bUse) {
      ia = constr_iatomptr(ncon1,ia1,ia2,con);
      /* Flexible constraints currently have length 0, which is incorrect */
      if (!bTopB)
	len = iparams[ia[0]].constr.dA;
      else
	len = iparams[ia[0]].constr.dB;
      /* In the worst case the bond directions alternate */
      if (nc % 2 == 0) {
	rn0 = r0 + len;
	rn1 = r1;
      } else {
	rn0 = r0;
	rn1 = r1 + len;
      }
      /* Assume angles of 120 degrees between all bonds */
      if (rn0*rn0 + rn1*rn1 + rn0*rn1 > *r2max) {
	*r2max = rn0*rn0 + rn1*rn1 + r0*rn1;
	if (debug) {
	  fprintf(debug,"Found longer constraint distance: r0 %5.3f r1 %5.3f rmax %5.3f\n", rn0,rn1,sqrt(*r2max));
	  for(a1=0; a1<depth; a1++)
	    fprintf(debug," %d %5.3f",
		    path[a1],
		    iparams[constr_iatomptr(ncon1,ia1,ia2,con)[0]].constr.dA);
	  fprintf(debug," %d %5.3f\n",con,len);
	}
      }
      /* Limit the number of recursions to 1000*nc,
       * so a call does not take more than a second,
       * even for highly connected systems.
       */
      if (depth + 1 < nc && *count < 1000*nc) {
	if (ia[1] == at)
	  a1 = ia[2];
	else
	  a1 = ia[1];
	/* Recursion */
	path[depth] = con;
	constr_recur(at2con,ilist,iparams,
		     bTopB,a1,depth+1,nc,path,rn0,rn1,r2max,count);
	path[depth] = -1;
      }
    }
  }
}

static real constr_r_max_moltype(FILE *fplog,
				 gmx_moltype_t *molt,t_iparams *iparams,
				 t_inputrec *ir)
{
  int natoms,nflexcon,*path,at,count;

  t_blocka at2con;
  real r0,r1,r2maxA,r2maxB,rmax,lam0,lam1;

  if (molt->ilist[F_CONSTR].nr   == 0 &&
      molt->ilist[F_CONSTRNC].nr == 0) {
    return 0;
  }
  
  natoms = molt->atoms.nr;

  at2con = make_at2con(0,natoms,molt->ilist,iparams,
		       EI_DYNAMICS(ir->eI),&nflexcon);
  snew(path,1+ir->nProjOrder);
  for(at=0; at<1+ir->nProjOrder; at++)
    path[at] = -1;

  r2maxA = 0;
  for(at=0; at<natoms; at++) {
    r0 = 0;
    r1 = 0;

    count = 0;
    constr_recur(&at2con,molt->ilist,iparams,
		 FALSE,at,0,1+ir->nProjOrder,path,r0,r1,&r2maxA,&count);
  }
  if (ir->efep == efepNO) {
    rmax = sqrt(r2maxA);
  } else {
    r2maxB = 0;
    for(at=0; at<natoms; at++) {
      r0 = 0;
      r1 = 0;
      count = 0;
      constr_recur(&at2con,molt->ilist,iparams,
		   TRUE,at,0,1+ir->nProjOrder,path,r0,r1,&r2maxB,&count);
    }
    lam0 = ir->init_lambda;
    if (EI_DYNAMICS(ir->eI))
      lam0 += ir->init_step*ir->delta_lambda;
    rmax = (1 - lam0)*sqrt(r2maxA) + lam0*sqrt(r2maxB);
    if (EI_DYNAMICS(ir->eI)) {
      lam1 = ir->init_lambda + (ir->init_step + ir->nsteps)*ir->delta_lambda;
      rmax = max(rmax,(1 - lam1)*sqrt(r2maxA) + lam1*sqrt(r2maxB));
    }
  }

  done_blocka(&at2con);
  sfree(path);

  return rmax;
}

real constr_r_max(FILE *fplog,gmx_mtop_t *mtop,t_inputrec *ir)
{
  int mt;
  real rmax;

  rmax = 0;
  for(mt=0; mt<mtop->nmoltype; mt++) {
    rmax = max(rmax,
	       constr_r_max_moltype(fplog,&mtop->moltype[mt],
				    mtop->ffparams.iparams,ir));
  }
  
  if (fplog)
    fprintf(fplog,"Maximum distance for %d constraints, at 120 deg. angles, all-trans: %.3f nm\n",1+ir->nProjOrder,rmax);

  return rmax;
}

gmx_constr_t init_constraints(FILE *fplog,
                              gmx_mtop_t *mtop,t_inputrec *ir,
                              gmx_edsam_t ed,t_state *state,
                              t_commrec *cr)
{
    int  ncon,nset,nmol,settle_type,i,natoms,mt,nflexcon;
    struct gmx_constr *constr;
    char *env;
    t_ilist *ilist;
    gmx_mtop_ilistloop_t iloop;
    
    ncon =
        gmx_mtop_ftype_count(mtop,F_CONSTR) +
        gmx_mtop_ftype_count(mtop,F_CONSTRNC);
    nset = gmx_mtop_ftype_count(mtop,F_SETTLE);
    
    if (ncon+nset == 0 && ir->ePull != epullCONSTRAINT && ed == NULL) 
    {
        return NULL;
    }
    
    snew(constr,1);
    
    constr->ncon_tot = ncon;
    constr->nflexcon = 0;
    if (ncon > 0) 
    {
        constr->n_at2con_mt = mtop->nmoltype;
        snew(constr->at2con_mt,constr->n_at2con_mt);
        for(mt=0; mt<mtop->nmoltype; mt++) 
        {
            constr->at2con_mt[mt] = make_at2con(0,mtop->moltype[mt].atoms.nr,
                                                mtop->moltype[mt].ilist,
                                                mtop->ffparams.iparams,
                                                EI_DYNAMICS(ir->eI),&nflexcon);
            for(i=0; i<mtop->nmolblock; i++) 
            {
                if (mtop->molblock[i].type == mt) 
                {
                    constr->nflexcon += mtop->molblock[i].nmol*nflexcon;
                }
            }
        }
        
        if (constr->nflexcon > 0) 
        {
            if (fplog) 
            {
                fprintf(fplog,"There are %d flexible constraints\n",
                        constr->nflexcon);
                if (ir->fc_stepsize == 0) 
                {
                    fprintf(fplog,"\n"
                            "WARNING: step size for flexible constraining = 0\n"
                            "         All flexible constraints will be rigid.\n"
                            "         Will try to keep all flexible constraints at their original length,\n"
                            "         but the lengths may exhibit some drift.\n\n");
                    constr->nflexcon = 0;
                }
            }
            if (constr->nflexcon > 0) 
            {
                please_cite(fplog,"Hess2002");
            }
        }
        
        if (ir->eConstrAlg == econtLINCS) 
        {
            constr->lincsd = init_lincs(fplog,mtop,
                                        constr->nflexcon,constr->at2con_mt,
                                        DOMAINDECOMP(cr) && cr->dd->bInterCGcons,
                                        ir->nLincsIter,ir->nProjOrder);
        }
        
        if (ir->eConstrAlg == econtSHAKE) {
            if (DOMAINDECOMP(cr) && cr->dd->bInterCGcons)
            {
                gmx_fatal(FARGS,"SHAKE is not supported with domain decomposition and constraint that cross charge group boundaries, use LINCS");
            }
            if (constr->nflexcon) 
            {
                gmx_fatal(FARGS,"For this system also velocities and/or forces need to be constrained, this can not be done with SHAKE, you should select LINCS");
            }
            please_cite(fplog,"Ryckaert77a");
            if (ir->bShakeSOR) 
            {
                please_cite(fplog,"Barth95a");
            }
        }
    }
  
    if (nset > 0) {
        please_cite(fplog,"Miyamoto92a");
        
        /* Check that we have only one settle type */
        settle_type = -1;
        iloop = gmx_mtop_ilistloop_init(mtop);
        while (gmx_mtop_ilistloop_next(iloop,&ilist,&nmol)) 
        {
            for (i=0; i<ilist[F_SETTLE].nr; i+=2) 
            {
                if (settle_type == -1) 
                {
                    settle_type = ilist[F_SETTLE].iatoms[i];
                } 
                else if (ilist[F_SETTLE].iatoms[i] != settle_type) 
                {
                    gmx_fatal(FARGS,"More than one settle type.\n"
                              "Suggestion: change the least use settle constraints into 3 normal constraints.");
                }
            }
        }
    }
    
    constr->maxwarn = 999;
    env = getenv("GMX_MAXCONSTRWARN");
    if (env) 
    {
        constr->maxwarn = 0;
        sscanf(env,"%d",&constr->maxwarn);
        if (fplog) 
        {
            fprintf(fplog,
                    "Setting the maximum number of constraint warnings to %d\n",
                    constr->maxwarn);
        }
        if (MASTER(cr)) 
        {
            fprintf(stderr,
                    "Setting the maximum number of constraint warnings to %d\n",
                    constr->maxwarn);
        }
    }
    if (constr->maxwarn < 0 && fplog) 
    {
        fprintf(fplog,"maxwarn < 0, will not stop on constraint errors\n");
    }
    constr->warncount_lincs  = 0;
    constr->warncount_settle = 0;
    
    /* Initialize the essential dynamics sampling.
     * Put the pointer to the ED struct in constr */
    constr->ed = ed;
    if (ed != NULL) 
    {
        init_edsam(mtop,ir,cr,ed,state->x,state->box);
    }
    
    constr->warn_mtop = mtop;
    
    return constr;
}

t_blocka *atom2constraints_moltype(gmx_constr_t constr)
{
  return constr->at2con_mt;
}


gmx_bool inter_charge_group_constraints(gmx_mtop_t *mtop)
{
  const gmx_moltype_t *molt;
  const t_block *cgs;
  const t_ilist *il;
  int  mb;
  int  nat,*at2cg,cg,a,ftype,i;
  gmx_bool bInterCG;

  bInterCG = FALSE;
  for(mb=0; mb<mtop->nmolblock && !bInterCG; mb++) {
    molt = &mtop->moltype[mtop->molblock[mb].type];

    if (molt->ilist[F_CONSTR].nr   > 0 ||
	molt->ilist[F_CONSTRNC].nr > 0) {
      cgs  = &molt->cgs;
      snew(at2cg,molt->atoms.nr);
      for(cg=0; cg<cgs->nr; cg++) {
	for(a=cgs->index[cg]; a<cgs->index[cg+1]; a++)
	  at2cg[a] = cg;
      }
      
      for(ftype=F_CONSTR; ftype<=F_CONSTRNC; ftype++) {
	il = &molt->ilist[ftype];
	for(i=0; i<il->nr && !bInterCG; i+=3) {
	  if (at2cg[il->iatoms[i+1]] != at2cg[il->iatoms[i+2]])
	    bInterCG = TRUE;
	}
      }
      sfree(at2cg);
    }
  }

  return bInterCG;
}
