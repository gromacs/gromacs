/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team,
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
 * Gallium Rubidium Oxygen Manganese Argon Carbon Silicon
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>
#include <string.h>

#include "typedefs.h"
#include "smalloc.h"
#include "genborn.h"
#include "vec.h"
#include "grompp.h"
#include "pdbio.h"
#include "names.h"
#include "physics.h"
#include "partdec.h"
#include "domdec.h"
#include "network.h"
#include "gmx_fatal.h"
#include "mtop_util.h"

#ifdef GMX_MPI
#include "mpi.h"
#endif

#if ( defined(GMX_IA32_SSE) || defined(GMX_X86_64_SSE) || defined(GMX_SSE2) )
#ifdef GMX_DOUBLE
#include "genborn_sse2_double.h"
#else
#include "genborn_sse2_single.h"
#endif /* GMX_DOUBLE */
#endif


/* Still parameters - make sure to edit in genborn_sse.c too if you change these! */
#define STILL_P1  0.073*0.1              /* length        */
#define STILL_P2  0.921*0.1*CAL2JOULE    /* energy*length */
#define STILL_P3  6.211*0.1*CAL2JOULE    /* energy*length */
#define STILL_P4  15.236*0.1*CAL2JOULE
#define STILL_P5  1.254 

#define STILL_P5INV (1.0/STILL_P5)
#define STILL_PIP5  (M_PI*STILL_P5)

int init_gb_nblist(int natoms, t_nblist *nl)
{
	nl->maxnri      = natoms*4;
	nl->maxnrj      = 0;
	nl->maxlen      = 0;
	nl->nri         = 0;
	nl->nrj         = 0;
	nl->iinr        = NULL;
	nl->gid         = NULL;
	nl->shift       = NULL;
	nl->jindex      = NULL;
	/*nl->nltype      = nltype;*/
	
	srenew(nl->iinr,   nl->maxnri);
	srenew(nl->gid,    nl->maxnri);
	srenew(nl->shift,  nl->maxnri);
	srenew(nl->jindex, nl->maxnri+1);
	
	nl->jindex[0] = 0;
	
	return 0;
}

int print_nblist(int natoms, t_nblist *nl)
{
	int i,k,ai,aj,nj0,nj1;
	
	printf("genborn.c: print_nblist, natoms=%d\n",nl->nri); 
	
	for(i=0;i<nl->nri;i++)
	{
		ai=nl->iinr[i];
		nj0=nl->jindex[i];
		nj1=nl->jindex[i+1];
	
		for(k=nj0;k<nj1;k++)
		{	
			aj=nl->jjnr[k];
			printf("ai=%d, aj=%d\n",ai,aj);
		}
	}
	
	return 0;	
}

void fill_log_table(const int n, real *table)
{
	real numlog,logfactor;
	int i;
	int *const exp_ptr=((int*)&numlog);
	int x = *exp_ptr;
	
	int incr = 1 << (23-n);
	int p=pow(2,n);

	logfactor = 1.0/log(2.0);
	
	x=0x3F800000;
	*exp_ptr = x;
	
	for(i=0;i<p;++i)
	{
		table[i]=log(numlog)*logfactor; /* log2(numlog)=log(numlog)/log(2.0) */
		x+=incr;
		*exp_ptr=x;
	}
}


real table_log(real val, const real *table, const int n)
{
	int *const exp_ptr = ((int*)&val);
	int x              = *exp_ptr;
	const int log_2    = ((x>>23) & 255) - 127;
	x &= 0x7FFFFF;
	x = x >> (23-n);
	val = table[x];
	return ((val+log_2)*0.69314718);  
}

void gb_pd_send(t_commrec *cr, real *send_data, int nr)
{
#ifdef GMX_MPI	
	int i,cur;
	int *index,*sendc,*disp;
	
	snew(sendc,cr->nnodes);
	snew(disp,cr->nnodes);
	
	index = pd_index(cr);
	cur   = cr->nodeid;
	
	/* Setup count/index arrays */
	for(i=0;i<cr->nnodes;i++)
	{
		sendc[i]  = index[i+1]-index[i];
		disp[i]   = index[i];	
	}
	
	/* Do communication */
	MPI_Gatherv(send_data+index[cur],sendc[cur],GMX_MPI_REAL,send_data,sendc,disp,GMX_MPI_REAL,0,cr->mpi_comm_mygroup);
	MPI_Bcast(send_data,nr,GMX_MPI_REAL,0,cr->mpi_comm_mygroup);
	
#endif
}


int init_gb_plist(t_params *p_list)
{
	p_list->nr    = 0;
    p_list->param = NULL;
	
	return 0;
}


static void assign_gb_param(t_functype ftype,t_iparams *new,
							real old[MAXFORCEPARAM],int comb,real reppow)
{
	int  i,j;
	
  /* Set to zero */
  for(j=0; (j<MAXFORCEPARAM); j++)
    new->generic.buf[j]=0.0;
	
	switch (ftype) {
		case F_GB12:
		case F_GB13:
		case F_GB14:
			new->gb.c6A=old[0];
			new->gb.c12A=old[1];
			new->gb.c6B=old[2];
			new->gb.c12B=old[3];
			new->gb.sar=old[4];
			new->gb.st=old[5];
			new->gb.pi=old[6];
			new->gb.gbr=old[7];
			new->gb.bmlt=old[8];
			break;
		default:
			gmx_fatal(FARGS,"unknown function type %d in %s line %d",
								ftype,__FILE__,__LINE__);		
	}
}

static void 
append_gb_interaction(t_ilist *ilist,
					  int type,int nral,atom_id a[MAXATOMLIST])
{
  int i,where1;
  
  where1     = ilist->nr;
  ilist->nr += nral+1;
	
  ilist->iatoms[where1++]=type;
  for (i=0; (i<nral); i++) 
  {
    ilist->iatoms[where1++]=a[i];
  }
}


static int 
enter_gb_params(gmx_ffparams_t *ffparams, t_functype ftype,
				real forceparams[MAXFORCEPARAM],int comb,real reppow,
				int start,bool bAppend)
{
  t_iparams new;
  int       type;
	
	assign_gb_param(ftype,&new,forceparams,comb,reppow);
  if (!bAppend) {
		for (type=start; (type<ffparams->ntypes); type++) {
      if (ffparams->functype[type]==ftype) {
					if (memcmp(&new,&ffparams->iparams[type],(size_t)sizeof(new)) == 0)
					return type;
      }
    }
  }
  else
    type=ffparams->ntypes;
  if (debug)
    fprintf(debug,"copying new to idef->iparams[%d] (ntypes=%d)\n",
						type,ffparams->ntypes);
  memcpy(&ffparams->iparams[type],&new,(size_t)sizeof(new));
  
  ffparams->ntypes++;
  ffparams->functype[type]=ftype;
	
  return type;
}

int init_gb_still(t_commrec *cr, t_forcerec  *fr, t_atomtypes *atype, t_idef *idef, t_atoms *atoms, gmx_genborn_t *born,int natoms)
{
	
	int i,j,i1,i2,k,m,nbond,nang,ia,ib,ic,id,nb,idx,idx2,at;
	int iam,ibm;
	int at0,at1;
	real length,angle;
	real r,ri,rj,ri2,ri3,rj2,r2,r3,r4,rk,ratio,term,h,doffset;
	real p1,p2,p3,factor,cosine,rab,rbc;
	
	real *vsol;
	real *gp;
	
	snew(vsol,natoms);
	snew(gp,natoms);
	
	if(PAR(cr))
	{
		if(PARTDECOMP(cr))
		{
			pd_at_range(cr,&at0,&at1);
			
			for(i=0;i<natoms;i++)
			{
				vsol[i] = gp[i] = 0;
			}
		}
		else
		{
			at0 = 0;
			at1 = natoms;
		}
	}
	else
	{
		at0 = 0;
		at1 = natoms;
	}
	
	doffset = born->gb_doffset;
	
	for(i=0;i<natoms;i++)
	{
		born->gpol_globalindex[i]=born->vsolv_globalindex[i]=0; 	
	}
	
	/* Compute atomic solvation volumes for Still method */
	for(i=0;i<natoms;i++)
	{	
		ri=atype->gb_radius[atoms->atom[i].type];
		r3=ri*ri*ri;
		born->vsolv_globalindex[i]=(4*M_PI/3)*r3;
	}

	for(j=0;j<idef->il[F_GB12].nr;j+=3)
	{
		m=idef->il[F_GB12].iatoms[j];
		ia=idef->il[F_GB12].iatoms[j+1];
		ib=idef->il[F_GB12].iatoms[j+2];
		
		r=1.01*idef->iparams[m].gb.st;
		
		ri   = atype->gb_radius[atoms->atom[ia].type];
		rj   = atype->gb_radius[atoms->atom[ib].type];
		
		ri2  = ri*ri;
		ri3  = ri2*ri;
		rj2  = rj*rj;
		
		ratio  = (rj2-ri2-r*r)/(2*ri*r);
		h      = ri*(1+ratio);
		term   = (M_PI/3.0)*h*h*(3.0*ri-h);

		if(PARTDECOMP(cr))
		{
			vsol[ia]+=term;
		}
		else
		{
			born->vsolv_globalindex[ia] -= term;
		}
		
		ratio  = (ri2-rj2-r*r)/(2*rj*r);
		h      = rj*(1+ratio);
		term   = (M_PI/3.0)*h*h*(3.0*rj-h);
		
		if(PARTDECOMP(cr))
		{
			vsol[ib]+=term;
		}
		else
		{
			born->vsolv_globalindex[ib] -= term;
		}		
	}
	
	if(PARTDECOMP(cr))
	{
		gmx_sum(natoms,vsol,cr);
		
		for(i=0;i<natoms;i++)
		{
			born->vsolv_globalindex[i]=born->vsolv_globalindex[i]-vsol[i];
		}
	}
	
	/* Get the self-, 1-2 and 1-3 polarization energies for analytical Still method */
	/* Self */
	for(j=0;j<natoms;j++)
	{
		if(born->vs_globalindex[j]==1)
		{
			born->gpol_globalindex[j]=-0.5*ONE_4PI_EPS0/(atype->gb_radius[atoms->atom[j].type]-doffset+STILL_P1);
		}
	}
	
	/* 1-2 */
	for(j=0;j<idef->il[F_GB12].nr;j+=3)
	{
		m=idef->il[F_GB12].iatoms[j];
		ia=idef->il[F_GB12].iatoms[j+1];
		ib=idef->il[F_GB12].iatoms[j+2];
		
		r=idef->iparams[m].gb.st;
		
		r4=r*r*r*r;

		if(PARTDECOMP(cr))
		{
			gp[ia]+=STILL_P2*born->vsolv_globalindex[ib]/r4;
			gp[ib]+=STILL_P2*born->vsolv_globalindex[ia]/r4;
		}
		else
		{
			born->gpol_globalindex[ia]=born->gpol_globalindex[ia]+STILL_P2*born->vsolv_globalindex[ib]/r4;
			born->gpol_globalindex[ib]=born->gpol_globalindex[ib]+STILL_P2*born->vsolv_globalindex[ia]/r4;
		}
	}

	/* 1-3 */
	for(j=0;j<idef->il[F_GB13].nr;j+=3)
	{
		m=idef->il[F_GB13].iatoms[j];
		ia=idef->il[F_GB13].iatoms[j+1];
		ib=idef->il[F_GB13].iatoms[j+2];
	
		r=idef->iparams[m].gb.st;
		r4=r*r*r*r;
	
		if(PARTDECOMP(cr))
		{
			gp[ia]+=STILL_P3*born->vsolv[ib]/r4;
			gp[ib]+=STILL_P3*born->vsolv[ia]/r4;
		}
		else
		{
			born->gpol_globalindex[ia]=born->gpol_globalindex[ia]+STILL_P3*born->vsolv_globalindex[ib]/r4;
			born->gpol_globalindex[ib]=born->gpol_globalindex[ib]+STILL_P3*born->vsolv_globalindex[ia]/r4;
		}		
	}
	
	if(PARTDECOMP(cr))
	{
		gmx_sum(natoms,gp,cr);
		
		for(i=0;i<natoms;i++)
		{
			born->gpol_globalindex[i]=born->gpol_globalindex[i]+gp[i];
		}	
	}
	
	sfree(vsol);
	sfree(gp);
		
	return 0;
}



#define LOG_TABLE_ACCURACY 15 /* Accuracy of the table logarithm */


/* Initialize all GB datastructs and compute polarization energies */
int init_gb(gmx_genborn_t **p_born,t_commrec *cr, t_forcerec *fr, t_inputrec *ir,
			gmx_mtop_t *mtop, rvec x[], real rgbradii, int gb_algorithm)
{
	int i,j,m,ai,aj,jj,natoms,nalloc;
	real rai,sk,p,doffset;
	
	t_atoms        atoms;
	gmx_genborn_t  *born;
	gmx_localtop_t *localtop;

	natoms   = mtop->natoms;
		
	atoms    = gmx_mtop_global_atoms(mtop);
	localtop = gmx_mtop_generate_local_top(mtop,ir);
	
	snew(born,1);
	*p_born = born;

	born->nr  = natoms;
	
	snew(born->drobc, natoms);
	snew(born->bRad,  natoms);
	
	/* Allocate memory for the global data arrays */
	snew(born->param_globalindex, natoms+3);
	snew(born->gpol_globalindex,  natoms+3);
	snew(born->vsolv_globalindex, natoms+3);
	snew(born->vs_globalindex,    natoms+3);
	
	snew(fr->invsqrta, natoms);
	snew(fr->dvda,     natoms);
	
	fr->dadx = NULL;
	fr->nalloc_dadx = 0;

	/* snew(born->asurf,natoms); */
	/* snew(born->dasurf,natoms); */

	/* Initialize the gb neighbourlist */
	init_gb_nblist(natoms, &(fr->gblist));
	
	/* Do the Vsites exclusions (if any) */
	for(i=0;i<natoms;i++)
	{
		jj = atoms.atom[i].type;
		born->vs_globalindex[i]=1;																							
				
		/* If we have a Vsite, put vs_globalindex[i]=0 */
	    if(C6(fr->nbfp,fr->ntype,jj,jj)==0 && C12(fr->nbfp,fr->ntype,jj,jj)==0 && atoms.atom[i].q==0)
		{
			born->vs_globalindex[i]=0;
		}
	}
	
	/* Copy algorithm parameters from inputrecord to local structure */
	born->obc_alpha  = ir->gb_obc_alpha;
	born->obc_beta   = ir->gb_obc_beta;
	born->obc_gamma  = ir->gb_obc_gamma;
	born->gb_doffset = ir->gb_dielectric_offset;
	
	doffset = born->gb_doffset;
	
	/* If Still model, initialise the polarisation energies */
	if(gb_algorithm==egbSTILL)	
	{
		init_gb_still(cr, fr,&(mtop->atomtypes), &(localtop->idef), &atoms, born, natoms);	
	}
	/* If HCT/OBC,  precalculate the sk*atype->S_hct factors */
	else if(gb_algorithm==egbHCT || gb_algorithm==egbOBC)
	{
		for(i=0;i<natoms;i++)
		{	
			if(born->vs_globalindex[i]==1)
			{
				rai            = mtop->atomtypes.gb_radius[atoms.atom[i].type]-doffset; 
				sk             = rai * mtop->atomtypes.S_hct[atoms.atom[i].type];
				born->param_globalindex[i] = sk;
			}
			else
			{
				born->param_globalindex[i] = 0;
			}
		}
	}
		
	/* Init the logarithm table */
	p=pow(2,LOG_TABLE_ACCURACY);
	snew(born->log_table, p);
	
	fill_log_table(LOG_TABLE_ACCURACY, born->log_table);
	
	/* Allocate memory for work arrays for temporary use */
	snew(born->work,natoms+4);
	snew(born->count,natoms);
	snew(born->nblist_work,natoms);
	
	/* If less than 500 atoms (completely arbitrary number), allocate exact */
	if(natoms<500)
	{
		born->nblist_work_nalloc = natoms;
	}
	else
	{
		born->nblist_work_nalloc = 500;
	}
	
	/* Allocate memory for the second dimension in the array used to setup the gb neighbourlist */
	for(i=0;i<natoms;i++)
	{
		snew(born->nblist_work[i],born->nblist_work_nalloc);
	}
	
	/* Domain decomposition specific stuff */
	if(DOMAINDECOMP(cr))
	{
		snew(born->dd_work,natoms);
		born->nlocal = cr->dd->nat_tot; /* cr->dd->nat_tot will be zero here */
	}
	

	
	return 0;
}

int generate_gb_topology(gmx_mtop_t *mtop, t_molinfo *mi)
{
	int i,j,k,type,m,a1,a2,a3,a4,nral,maxtypes,start,comb,mt;
	int natoms,n12,n13,n14,a1type,a2type,a3type;
	double p1,p2,p3,cosine,r2,rab,rbc;
	
	t_params *plist;
	t_params plist_gb12;
	t_params plist_gb13;
	t_params plist_gb14;
	genborn_bonds_t *bonds;

	gmx_ffparams_t *ffp;
	gmx_moltype_t *molt;
	
	/* To keep the compiler happy */
	rab=rbc=0;

	p1=STILL_P1;
	p2=STILL_P2;
	p3=STILL_P3;
	
	plist_gb12.param = NULL;
	plist_gb13.param = NULL;
	plist_gb14.param = NULL;
	bonds             = NULL;
	
	for(mt=0;mt<mtop->nmoltype;mt++)
	{
		plist = mi[mt].plist;
		plist_gb12.nr = plist_gb13.nr = plist_gb14.nr = 0;
		n12 = n13 = n14 = 0;
		
		for(i=0;i<F_NRE;i++)
		{
			if(IS_CHEMBOND(i))
			{
				n12 += plist[i].nr;
			}
			else if(IS_ANGLE(i))
			{
				n13 += plist[i].nr;
			}
			else if(i==F_LJ14)
			{
				n14 += plist[i].nr;
			}
		}
		
		srenew(plist_gb12.param,n12);
		srenew(plist_gb13.param,n13);
		srenew(plist_gb14.param,n14);
		srenew(bonds,mi[mt].atoms.nr);
		
		for(i=0;i<mi[mt].atoms.nr;i++)
		{
			bonds[i].nbonds=0;
		}
			
		/* We need to find all bonds lengths first */
		for(i=0;i<F_NRE;i++)
		{
			if(IS_CHEMBOND(i))
			{
				for(j=0;j<plist[i].nr; j++)
				{
					a1=plist[i].param[j].a[0];
					a2=plist[i].param[j].a[1];

					if(mi[mt].atoms.atom[a1].q!=0 && mi[mt].atoms.atom[a2].q!=0)
					{
						bonds[a1].bond[bonds[a1].nbonds]=a2;
						bonds[a1].length[bonds[a1].nbonds]=plist[i].param[j].c[0];
						bonds[a1].nbonds++;
						bonds[a2].bond[bonds[a2].nbonds]=a1;
						bonds[a2].length[bonds[a2].nbonds]=plist[i].param[j].c[0];
						bonds[a2].nbonds++;
						
						plist_gb12.param[plist_gb12.nr].a[0]=a1;
						plist_gb12.param[plist_gb12.nr].a[1]=a2;
						
						/* LJ parameters */
						plist_gb12.param[plist_gb12.nr].c[0]=-1;
						plist_gb12.param[plist_gb12.nr].c[1]=-1;
						plist_gb12.param[plist_gb12.nr].c[2]=-1;
						plist_gb12.param[plist_gb12.nr].c[3]=-1;
						
						/* GBSA parameters */
						a1type = mi[mt].atoms.atom[a1].type;
						a2type = mi[mt].atoms.atom[a2].type;
						
						plist_gb12.param[plist_gb12.nr].c[4]=mtop->atomtypes.radius[a1type]+mtop->atomtypes.radius[a2type];	
						plist_gb12.param[plist_gb12.nr].c[5]=plist[i].param[j].c[0]; /* bond length */
						plist_gb12.param[plist_gb12.nr].c[6]=p2;
						plist_gb12.param[plist_gb12.nr].c[7]=mtop->atomtypes.gb_radius[a1type]+mtop->atomtypes.gb_radius[a2type];
						plist_gb12.param[plist_gb12.nr].c[8]=0.8875;
						plist_gb12.nr++;
					}
				}
			}
		}
		
		/* Now we detect angles and 1,4 pairs in parallel */
		for(i=0;i<F_NRE;i++)
		{
			if(IS_ANGLE(i))
			{
				if(F_G96ANGLES==i && plist[i].nr>0)
				{
					gmx_fatal(FARGS,"Cannot do GB with Gromos96 angles - yet\n");
				}
				
				for(j=0;j<plist[i].nr; j++)
				{
					a1=plist[i].param[j].a[0];
					a2=plist[i].param[j].a[1];
					a3=plist[i].param[j].a[2];	
					
					plist_gb13.param[plist_gb13.nr].a[0]=a1;
					plist_gb13.param[plist_gb13.nr].a[1]=a3;
					
					/* LJ parameters */	
					plist_gb13.param[plist_gb13.nr].c[0]=-1;
					plist_gb13.param[plist_gb13.nr].c[1]=-1;
					plist_gb13.param[plist_gb13.nr].c[2]=-1;
					plist_gb13.param[plist_gb13.nr].c[3]=-1;
					
					/* GBSA parameters */
					a1type = mi[mt].atoms.atom[a1].type;
					a3type = mi[mt].atoms.atom[a3].type;
					
					plist_gb13.param[plist_gb13.nr].c[4]=mtop->atomtypes.radius[a1type]+mtop->atomtypes.radius[a3type];	
					
					for(k=0;k<bonds[a2].nbonds;k++)
					{
						if(bonds[a2].bond[k]==a1)
						{
							rab = bonds[a2].length[k];
						}
						else if(bonds[a2].bond[k]==a3)
						{
							rbc=bonds[a2].length[k];
						}
					}
					
					cosine=cos(plist[i].param[j].c[0]/RAD2DEG);
					r2=rab*rab+rbc*rbc-(2*rab*rbc*cosine);
					plist_gb13.param[plist_gb13.nr].c[5]=sqrt(r2);
					plist_gb13.param[plist_gb13.nr].c[6]=p3;
					plist_gb13.param[plist_gb13.nr].c[7]=mtop->atomtypes.gb_radius[a1type]+mtop->atomtypes.gb_radius[a3type];
					plist_gb13.param[plist_gb13.nr].c[8]=0.3516;
					plist_gb13.nr++;
				}	
			}
			
			if(F_LJ14 == i)
			{
				for(j=0;j<plist[F_LJ14].nr; j++)
				{
					a1=plist[F_LJ14].param[j].a[0];
					a2=plist[F_LJ14].param[j].a[1];
					
					plist_gb14.param[plist_gb14.nr].a[0]=a1;
					plist_gb14.param[plist_gb14.nr].a[1]=a2;
					
					plist_gb14.param[plist_gb14.nr].c[0]=-1;
					plist_gb14.param[plist_gb14.nr].c[1]=-1;
					plist_gb14.param[plist_gb14.nr].c[2]=-1;
					plist_gb14.param[plist_gb14.nr].c[3]=-1;
					
					/* GBSA parameters */
					a1type = mi[mt].atoms.atom[a1].type;
					a2type = mi[mt].atoms.atom[a2].type;
					
					plist_gb14.param[plist_gb14.nr].c[4]=mtop->atomtypes.radius[a1type]+mtop->atomtypes.radius[a2type];	
					plist_gb14.param[plist_gb14.nr].c[5]=-1;
					plist_gb14.param[plist_gb14.nr].c[6]=p3;
					plist_gb14.param[plist_gb14.nr].c[7]=mtop->atomtypes.gb_radius[a1type]+mtop->atomtypes.gb_radius[a2type];
					plist_gb14.param[plist_gb14.nr].c[8]=0.3516;
					plist_gb14.nr++;
				}	
			}
		}

		/* Put GB 1-2, 1-3, and 1-4 interactions into topology, per moleculetype */
		ffp  = &mtop->ffparams;
		molt = &mtop->moltype[mt];
	
		convert_gb_params(ffp, F_GB12, &plist_gb12,&molt->ilist[F_GB12]);
		convert_gb_params(ffp, F_GB13, &plist_gb13,&molt->ilist[F_GB13]);
		convert_gb_params(ffp, F_GB14, &plist_gb14,&molt->ilist[F_GB14]); 
	}
	
	return 0;
	
}

int convert_gb_params(gmx_ffparams_t *ffparams, t_functype ftype, t_params *gb_plist, t_ilist *il)
{
	int k,nr,nral,delta,maxtypes,comb,type,start;
	real reppow;
	
	nral     = NRAL(ftype);
	start    = ffparams->ntypes;
	maxtypes = ffparams->ntypes; 
	nr	     = gb_plist->nr;
	comb     = 3;
	reppow   = 12;
	
	for (k=0; k<nr; k++) 
	{
		if (maxtypes <= ffparams->ntypes) 
		{
			maxtypes += 1000;
			srenew(ffparams->functype,maxtypes);
			srenew(ffparams->iparams, maxtypes);
			
			if (debug) 
				fprintf(debug,"%s, line %d: srenewed idef->functype and idef->iparams to %d\n",
						__FILE__,__LINE__,maxtypes);
		}
	
		type=enter_gb_params(ffparams,ftype,gb_plist->param[k].c,comb,reppow,start,0);
		
		delta = nr*(nral+1);
		srenew(il->iatoms,il->nr+delta);
		
		append_gb_interaction(il,type,nral,gb_plist->param[k].a);
	}
	
	return 0;

}

static int
calc_gb_rad_still(t_commrec *cr, t_forcerec *fr,int natoms, gmx_localtop_t *top,
				  const t_atomtypes *atype, rvec x[], t_nblist *nl, gmx_genborn_t *born,t_mdatoms *md)
{	
	int i,k,n,nj0,nj1,ai,aj,type;
	real gpi,dr,dr2,dr4,idr4,rvdw,ratio,ccf,theta,term,rai,raj;
	real ix1,iy1,iz1,jx1,jy1,jz1,dx11,dy11,dz11;
	real rinv,idr2,idr6,vaj,dccf,cosq,sinq,prod,gpi2;
	real factor;
	real *sum_gpi;
	
	factor  = 0.5*ONE_4PI_EPS0;
	n       = 0;
	sum_gpi = born->work;
	
	memset(sum_gpi,0,sizeof(real)*born->nr);
	
	for(i=0;i<nl->nri;i++ )
	{
		ai  = nl->iinr[i];
		
		nj0 = nl->jindex[ai];			
		nj1 = nl->jindex[ai+1];
		
		if(DOMAINDECOMP(cr))
		{
			gpi = 0;
		}
		else
		{
			gpi = born->gpol[ai];
		}
		
		rai = top->atomtypes.gb_radius[md->typeA[ai]];
		
		ix1 = x[ai][0];
		iy1 = x[ai][1];
		iz1 = x[ai][2];
		
		for(k=nj0;k<nj1;k++)
		{
			aj    = nl->jjnr[k];
		
			jx1   = x[aj][0];
			jy1   = x[aj][1];
			jz1   = x[aj][2];
			
			dx11  = ix1-jx1;
			dy11  = iy1-jy1;
			dz11  = iz1-jz1;
			
			dr2   = dx11*dx11+dy11*dy11+dz11*dz11; 
			rinv  = invsqrt(dr2);
			idr2  = rinv*rinv;
			idr4  = idr2*idr2;
			idr6  = idr4*idr2;
			
			raj = top->atomtypes.gb_radius[md->typeA[aj]];
			
			rvdw  = rai + raj;
			
			ratio = dr2 / (rvdw * rvdw);
			vaj   = born->vsolv[aj];
			
			if(ratio>STILL_P5INV) 
			{
				ccf=1.0;
				dccf=0.0;
			}
			else
			{
				theta = ratio*STILL_PIP5;
				cosq  = cos(theta);
				term  = 0.5*(1.0-cosq);
				ccf   = term*term;
				sinq  = 1.0 - cosq*cosq;
				dccf  = 2.0*term*sinq*invsqrt(sinq)*STILL_PIP5*ratio;
			}
			
			prod          = STILL_P4*vaj;
			gpi           = gpi+prod*ccf*idr4;
			fr->dadx[n++] = prod*(4*ccf-dccf)*idr6;
		}
		
		if(PAR(cr))
		{
			sum_gpi[ai] = gpi;
		}
		else
		{
			gpi2 = gpi * gpi;
			born->bRad[ai]   = factor*invsqrt(gpi2);
			fr->invsqrta[ai] = invsqrt(born->bRad[ai]);
		}
	}
	
	if(DOMAINDECOMP(cr))
	{
		dd_atom_sum_real(cr->dd,sum_gpi);
		
		for(i=0;i<nl->nri;i++)
		{
			ai   = nl->iinr[i];
			gpi  = born->gpol[ai];
			gpi  = gpi + sum_gpi[ai];
			gpi2 = gpi*gpi;
			
			born->bRad[ai]=factor*invsqrt(gpi2);
			fr->invsqrta[ai]=invsqrt(born->bRad[ai]);
		}
		
		dd_atom_spread_real(cr->dd,born->bRad);
		dd_atom_spread_real(cr->dd,fr->invsqrta);
	}

	return 0;
	
}
	

static int 
calc_gb_rad_hct(t_commrec *cr,t_forcerec *fr,int natoms, gmx_localtop_t *top,
					const t_atomtypes *atype, rvec x[], t_nblist *nl, gmx_genborn_t *born,t_mdatoms *md)
{
	int i,k,n,ai,aj,nj0,nj1,at0,at1;
	real rai,raj,gpi,dr2,dr,sk,sk2,lij,uij,diff2,tmp,sum_ai;
	real rad,min_rad,rinv,rai_inv;
	real ix1,iy1,iz1,jx1,jy1,jz1,dx11,dy11,dz11;
	real lij2, uij2, lij3, uij3, t1,t2,t3;
	real lij_inv,dlij,duij,sk2_rinv,prod,log_term;
	rvec dx;
	real doffset;
	real *sum_tmp;

	doffset = born->gb_doffset;
	sum_tmp = born->work;
	
	/* Keep the compiler happy */
	n    = 0;
	prod = 0;
		
	for(i=0;i<nl->nri;i++)
	{
		ai  = nl->iinr[i];
			
		nj0 = nl->jindex[ai];			
		nj1 = nl->jindex[ai+1];
		
		rai     = top->atomtypes.gb_radius[md->typeA[ai]]-doffset; 
		sum_ai  = 1.0/rai;
		rai_inv = sum_ai;
		
		ix1 = x[ai][0];
		iy1 = x[ai][1];
		iz1 = x[ai][2];
		
		if(PARTDECOMP(cr))
		{
			/* Only have the master node do this, since we only want one value at one time */
			if(MASTER(cr))
			{
				sum_tmp[ai] = sum_ai;
			}
			else
			{
				sum_tmp[ai] = 0;
			}
		}
		else if(DOMAINDECOMP(cr))
		{
			sum_ai = 0;
		}
		
		for(k=nj0;k<nj1;k++)
		{
			aj    = nl->jjnr[k];
				
			jx1   = x[aj][0];
			jy1   = x[aj][1];
			jz1   = x[aj][2];
			
			dx11  = ix1 - jx1;
			dy11  = iy1 - jy1;
			dz11  = iz1 - jz1;
			
			dr2   = dx11*dx11+dy11*dy11+dz11*dz11;
			rinv  = invsqrt(dr2);
			dr    = rinv*dr2;
			
			sk    = born->param[aj];
			
			if(rai < dr+sk)
			{
				lij     = 1.0/(dr-sk);
				dlij    = 1.0;
				
				if(rai>dr-sk) 
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
				
				lij_inv  = invsqrt(lij2);
				sk2      = sk*sk;
				sk2_rinv = sk2*rinv;
				prod     = 0.25*sk2_rinv;
				
				/* log_term = table_log(uij*lij_inv,born->log_table,LOG_TABLE_ACCURACY); */
				log_term = log(uij*lij_inv);
				
				tmp      = lij-uij + 0.25*dr*diff2 + (0.5*rinv)*log_term + prod*(-diff2);
								
				if(rai<sk-dr)
				{
					tmp = tmp + 2.0 * (rai_inv-lij);
				}
					
				/* duij    = 1.0; */
				t1      = 0.5*lij2 + prod*lij3 - 0.25*(lij*rinv+lij3*dr);
				t2      = -0.5*uij2 - 0.25*sk2_rinv*uij3 + 0.25*(uij*rinv+uij3*dr);
				t3      = 0.125*(1.0+sk2_rinv*rinv)*(-diff2)+0.25*log_term*rinv*rinv;
				
				fr->dadx[n++] = (dlij*t1+t2+t3)*rinv; /* rb2 is moved to chainrule	*/
				/* fr->dadx[n++] = (dlij*t1+duij*t2+t3)*rinv; */ /* rb2 is moved to chainrule	*/

				if(PARTDECOMP(cr))
				{
					sum_tmp[ai] -= 0.5*tmp;
				}
				else
				{
					sum_ai -= 0.5*tmp;
				}
			}
		}
	
		if(PAR(cr))
		{
			sum_tmp[ai] = sum_ai;
		}
		else
		{
			min_rad = rai + doffset;
			rad=1.0/sum_ai; 
			
			born->bRad[ai]=rad > min_rad ? rad : min_rad;
			fr->invsqrta[ai]=invsqrt(born->bRad[ai]);
		}
	}

	if(PARTDECOMP(cr))
	{
		pd_at_range(cr,&at0,&at1);
		gmx_sum(natoms,sum_tmp,cr);

		/* Calculate the Born radii so that they are available on all nodes */
		for(i=at0;i<at1;i++)
		{
			ai          = nl->iinr[i];
			min_rad     = top->atomtypes.gb_radius[md->typeA[ai]]; 
			sum_ai      = 1.0/(min_rad-doffset);
			sum_tmp[ai] = sum_ai + sum_tmp[ai];
			rad         = 1.0/sum_tmp[ai];
			
			born->bRad[ai]=rad > min_rad ? rad : min_rad;
			fr->invsqrta[ai]=invsqrt(born->bRad[ai]);
		}
		
		/* Communicate Born radii */
		gb_pd_send(cr,born->bRad,md->nr);
		gb_pd_send(cr,fr->invsqrta,md->nr);
	}
	else if(DOMAINDECOMP(cr))
	{
		dd_atom_sum_real(cr->dd,sum_tmp);
		
		for(i=0;i<nl->nri;i++)
		{
			ai = nl->iinr[i];
			min_rad = top->atomtypes.gb_radius[md->typeA[ai]]; 
			sum_ai = 1.0/(min_rad-doffset);
			sum_tmp[ai] = sum_ai + sum_tmp[ai];
			rad     = 1.0/sum_tmp[ai];
			born->bRad[ai]=rad > min_rad ? rad : min_rad;
			fr->invsqrta[ai]=invsqrt(born->bRad[ai]);
		}
		
		dd_atom_spread_real(cr->dd,born->bRad);
		dd_atom_spread_real(cr->dd,fr->invsqrta);
	}
	
	if(cr->nodeid==0)
	{
		for(i=0;i<nl->nri;i++)
		{
			if(DOMAINDECOMP(cr))
			{
				printf("i=%d, brad=%15.15f\n",cr->dd->gatindex[i],born->bRad[i]);
			}
			else
			{
				printf("i=%d, brad=%15.15f\n",i,born->bRad[i]);
			}
		}
	}
	
	return 0;
}

static int 
calc_gb_rad_obc(t_commrec *cr, t_forcerec *fr, int natoms, gmx_localtop_t *top,
					const t_atomtypes *atype, rvec x[], t_nblist *nl, gmx_genborn_t *born,t_mdatoms *md)
{
	int i,k,ai,aj,nj0,nj1,n,at0,at1;
	real rai,raj,gpi,dr2,dr,sk,sk2,lij,uij,diff2,tmp,sum_ai;
	real rad, min_rad,sum_ai2,sum_ai3,tsum,tchain,rinv,rai_inv,lij_inv,rai_inv2;
	real log_term,prod,sk2_rinv;
	real ix1,iy1,iz1,jx1,jy1,jz1,dx11,dy11,dz11;
	real lij2,uij2,lij3,uij3,dlij,duij,t1,t2,t3;
	real doffset;
	real *sum_tmp;

	/* Keep the compiler happy */
	n=0;
	prod=0;
	
	doffset = born->gb_doffset;
	sum_tmp = born->work;
	
	for(i=0;i<nl->nri;i++)
	{
		ai  = nl->iinr[i];
	
		nj0 = nl->jindex[ai];
		nj1 = nl->jindex[ai+1];
		
		rai      = top->atomtypes.gb_radius[md->typeA[ai]]-doffset;
		sum_ai   = 0;
		rai_inv  = 1.0/rai;
		rai_inv2 = 1.0/top->atomtypes.gb_radius[md->typeA[ai]];
		
		ix1 = x[ai][0];
		iy1 = x[ai][1];
		iz1 = x[ai][2];
				
		for(k=nj0;k<nj1;k++)
		{
			aj    = nl->jjnr[k];
			
			jx1   = x[aj][0];
			jy1   = x[aj][1];
			jz1   = x[aj][2];
			
			dx11  = ix1 - jx1;
			dy11  = iy1 - jy1;
			dz11  = iz1 - jz1;
			
			dr2   = dx11*dx11+dy11*dy11+dz11*dz11;
			rinv  = invsqrt(dr2);
			dr    = dr2*rinv;
		
			/* sk is precalculated in init_gb() */
			sk = born->param[aj];
						
			if(rai < dr+sk)
			{
				lij       = 1.0/(dr-sk);
				dlij      = 1.0; 
								
				if(rai>dr-sk)
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
				
				lij_inv  = invsqrt(lij2);
				sk2      = sk*sk;
				sk2_rinv = sk2*rinv;	
				prod     = 0.25*sk2_rinv;
				
				log_term = log(uij*lij_inv);
				/* log_term = table_log(uij*lij_inv,born->log_table,LOG_TABLE_ACCURACY); */
				tmp      = lij-uij + 0.25*dr*diff2 + (0.5*rinv)*log_term + prod*(-diff2);
				
				if(rai < sk-dr)
				{
					tmp = tmp + 2.0 * (rai_inv-lij);
				}
				
				/* duij    = 1.0; */
				t1      = 0.5*lij2 + prod*lij3 - 0.25*(lij*rinv+lij3*dr); 
				t2      = -0.5*uij2 - 0.25*sk2_rinv*uij3 + 0.25*(uij*rinv+uij3*dr); 
				t3      = 0.125*(1.0+sk2_rinv*rinv)*(-diff2)+0.25*log_term*rinv*rinv; 
					
				fr->dadx[n++] = (dlij*t1+t2+t3)*rinv; /* rb2 is moved to chainrule	*/
				/* fr->dadx[n++] = (dlij*t1+t2+t3)*rinv; */
				
				if(PARTDECOMP(cr))
				{
					sum_tmp[ai] += 0.5*tmp;
				}
				else
				{
					sum_ai += 0.5*tmp;
				
				}
			}	
		}
		
		if(PAR(cr))
		{
			sum_tmp[ai] = sum_ai;
		}
		else
		{
			sum_ai  = rai     * sum_ai;
			sum_ai2 = sum_ai  * sum_ai;
			sum_ai3 = sum_ai2 * sum_ai;
			
			tsum    = tanh(born->obc_alpha*sum_ai-born->obc_beta*sum_ai2+born->obc_gamma*sum_ai3);
			born->bRad[ai] = rai_inv - tsum*rai_inv2;
			born->bRad[ai] = 1.0 / born->bRad[ai];
			
			fr->invsqrta[ai]=invsqrt(born->bRad[ai]);
			
			tchain  = rai * (born->obc_alpha-2*born->obc_beta*sum_ai+3*born->obc_gamma*sum_ai2);
			born->drobc[ai] = (1.0-tsum*tsum)*tchain*rai_inv2;
		}
	}

	if(PARTDECOMP(cr))
	{
		pd_at_range(cr,&at0,&at1);
		gmx_sum(natoms,sum_tmp,cr);

		for(i=at0;i<at1;i++)
		{
			ai      = nl->iinr[i];
			rai     = top->atomtypes.gb_radius[md->typeA[ai]];
			rai_inv = 1.0/rai;
			
			sum_ai  = sum_tmp[ai];
			sum_ai  = rai     * sum_ai;
			sum_ai2 = sum_ai  * sum_ai;
			sum_ai3 = sum_ai2 * sum_ai;
			
			tsum             = tanh(born->obc_alpha*sum_ai-born->obc_beta*sum_ai2+born->obc_gamma*sum_ai3);
			born->bRad[ai]   = rai_inv - tsum/top->atomtypes.gb_radius[md->typeA[ai]];
			born->bRad[ai]   = 1.0 / born->bRad[ai];
			fr->invsqrta[ai] = invsqrt(born->bRad[ai]);
			
			tchain  = rai * (born->obc_alpha-2*born->obc_beta*sum_ai+3*born->obc_gamma*sum_ai2);
			born->drobc[ai] = (1.0-tsum*tsum)*tchain/top->atomtypes.gb_radius[md->typeA[ai]];
		}
		
		gb_pd_send(cr,born->bRad,md->nr);
		gb_pd_send(cr,fr->invsqrta,md->nr);
		gb_pd_send(cr,born->drobc,md->nr);
	}
	else if(DOMAINDECOMP(cr))
	{
		dd_atom_sum_real(cr->dd,sum_tmp);
	
		for(i=0;i<natoms;i++)
		{
			ai      = i;
			rai     = top->atomtypes.gb_radius[md->typeA[ai]]-doffset;
			rai_inv = 1.0/rai;
			
			sum_ai  = sum_tmp[ai];
			sum_ai  = rai     * sum_ai;
			sum_ai2 = sum_ai  * sum_ai;
			sum_ai3 = sum_ai2 * sum_ai;
			
			tsum    = tanh(born->obc_alpha*sum_ai-born->obc_beta*sum_ai2+born->obc_gamma*sum_ai3);
			born->bRad[ai] = rai_inv - tsum/top->atomtypes.gb_radius[md->typeA[ai]];
			
			born->bRad[ai] = 1.0 / born->bRad[ai];
			fr->invsqrta[ai]=invsqrt(born->bRad[ai]);
			
			tchain  = rai * (born->obc_alpha-2*born->obc_beta*sum_ai+3*born->obc_gamma*sum_ai2);
			born->drobc[ai] = (1.0-tsum*tsum)*tchain/top->atomtypes.gb_radius[md->typeA[ai]];
		}
		
		dd_atom_spread_real(cr->dd,born->bRad);
		dd_atom_spread_real(cr->dd,fr->invsqrta);
		dd_atom_spread_real(cr->dd,born->drobc);
	}
	
	return 0;
}



int calc_gb_rad(t_commrec *cr, t_forcerec *fr, t_inputrec *ir,gmx_localtop_t *top,
				const t_atomtypes *atype, rvec x[], rvec f[],t_nblist *nl, gmx_genborn_t *born,t_mdatoms *md)
{	
	/* First check that the allocated size of the dadx array is sufficient */
	if(nl->nrj > fr->nalloc_dadx)
	{
		fr->nalloc_dadx = 1.1*nl->nrj;
		srenew(fr->dadx,fr->nalloc_dadx);
	}
		
#ifdef GMX_DOUBLE
	
#if ( defined(GMX_IA32_SSE2) || defined(GMX_X86_64_SSE2) || defined(GMX_SSE2) )
	/* Currently, the double-precision sse code is disabled, since I have yet to implement
	 * double precision sse versions of sin/cos/log/exp-functions
	 */
	switch(ir->gb_algorithm)
	{
		case egbSTILL:
			/* calc_gb_rad_still_sse2_double(cr,fr,md->nr,mtop, atype, x[0], nl, born, md); */
			calc_gb_rad_still(cr,fr,born->nr,top, atype, x, nl, born, md); 
			break;
		case egbHCT:
			/* calc_gb_rad_hct_sse2_double(cr,fr,md->nr,top, atype, x[0], nl, born, md); */
			calc_gb_rad_hct(cr,fr,born->nr,top, atype, x, nl, born, md); 
			break;
		case egbOBC:
			/* calc_gb_rad_obc_sse2_double(cr,fr,md->nr,top, atype, x[0], nl, born, md); */
			calc_gb_rad_obc(cr,fr,born->nr,top, atype, x, nl, born, md);
			break;
			
		default:
			gmx_fatal(FARGS, "Unknown double precision sse-enabled algorithm for Born radii calculation: %d",ir->gb_algorithm);
	}
	
#else
	/* Switch for determining which algorithm to use for Born radii calculation */
	switch(ir->gb_algorithm)
	{
		case egbSTILL:
			calc_gb_rad_still(cr,fr,born->nr,top,atype,x,nl,born,md); 
			break;
		case egbHCT:
			calc_gb_rad_hct(cr,fr,born->nr,top,atype,x,nl,born,md); 
			break;
		case egbOBC:
			calc_gb_rad_obc(cr,fr,born->nr,top,atype,x,nl,born,md); 
			break;
			
		default:
			gmx_fatal(FARGS, "Unknown double precision algorithm for Born radii calculation: %d",ir->gb_algorithm);
	}
			
#endif
						
#else				
			
#if ( defined(GMX_IA32_SSE) || defined(GMX_X86_64_SSE) || defined(GMX_SSE2) )
	/* x86 or x86-64 with GCC inline assembly and/or SSE intrinsics */
	switch(ir->gb_algorithm)
	{
		case egbSTILL:
			calc_gb_rad_still_sse(cr,fr,born->nr,top, atype, x[0], nl, born, md);
			break;
		case egbHCT:
			calc_gb_rad_hct_sse(cr,fr,born->nr,top, atype, x[0], nl, born, md); 
			break;
		case egbOBC:
			calc_gb_rad_obc_sse(cr,fr,born->nr,top,atype,x[0],nl,born,md); 
			break;
			
		default:
			gmx_fatal(FARGS, "Unknown sse-enabled algorithm for Born radii calculation: %d",ir->gb_algorithm);
	}
	
#else
	
	/* Switch for determining which algorithm to use for Born radii calculation */
	switch(ir->gb_algorithm)
	{
		case egbSTILL:
			calc_gb_rad_still(cr,fr,born->nr,top,atype,x,nl,born,md); 
			break;
		case egbHCT:
			calc_gb_rad_hct(cr,fr,born->nr,top,atype,x,nl,born,md); 
			break;
		case egbOBC:
			calc_gb_rad_obc(cr,fr,born->nr,top,atype,x,nl,born,md); 
			break;
			
		default:
			gmx_fatal(FARGS, "Unknown algorithm for Born radii calculation: %d",ir->gb_algorithm);
	}
	
#endif /* Single precision sse */
			
#endif /* Double or single precision */
	
	return 0;		
}



real gb_bonds_tab(real *x, real *f, real *charge, real *p_gbtabscale,
				  real *invsqrta, real *dvda, real *GBtab, t_idef *idef,
				  real epsilon_r, real facel)
{
	int i,j,n0,nnn,type,ai,aj,ai3,aj3;
	real isai,isaj;
	real r,rsq11,ix1,iy1,iz1,jx1,jy1,jz1;
	real dx11,dy11,dz11,rinv11,iq,facel2;
	real isaprod,qq,gbscale,gbtabscale,Y,F,Geps,Heps2,Fp,VV,FF,rt,eps,eps2;
	real vgb,fgb,vcoul,fijC,dvdatmp,fscal,tx,ty,tz,dvdaj;
	real vctot;
	
	t_iatom *forceatoms;

	gbtabscale=*p_gbtabscale;
	vctot = 0.0;
	
	for(j=F_GB12;j<=F_GB14;j++)
	{
		forceatoms = idef->il[j].iatoms;
		
		for(i=0;i<idef->il[j].nr; )
		{
			type          = forceatoms[i++];
			ai            = forceatoms[i++];
			aj            = forceatoms[i++];
			ai3           = ai*3;
			aj3           = aj*3; 
			isai          = invsqrta[ai];
			ix1           = x[ai3+0];
			iy1           = x[ai3+1];
			iz1           = x[ai3+2];
			iq            = (-1)*facel*charge[ai];
			jx1           = x[aj3+0];
			jy1           = x[aj3+1];
			jz1           = x[aj3+2];
			dx11          = ix1 - jx1;
			dy11          = iy1 - jy1;
			dz11          = iz1 - jz1;
			rsq11         = dx11*dx11+dy11*dy11+dz11*dz11;
			rinv11        = invsqrt(rsq11);
			isaj          = invsqrta[aj];
			isaprod       = isai*isaj;
			qq            = isaprod*iq*charge[aj];
			gbscale       = isaprod*gbtabscale;
			r             = rsq11*rinv11;
			rt            = r*gbscale;
			n0            = rt;
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
			tx            = fgb*dx11;
			ty            = fgb*dy11;
			tz            = fgb*dz11;
	
			f[aj3+0]      = f[aj3+0] - tx;
			f[aj3+1]      = f[aj3+1] - ty;
			f[aj3+2]      = f[aj3+2] - tz;
		
			f[ai3+0]      = f[ai3+0] + tx;
			f[ai3+1]      = f[ai3+1] + ty;
			f[ai3+2]      = f[ai3+2] + tz;
		}
	}
	
	return vctot;
}


real calc_gb_selfcorrections(t_commrec *cr, int natoms, 
			     real *charge, gmx_genborn_t *born, real *dvda, t_mdatoms *md, double facel)
{	
	int i,ai,at0,at1;
	real rai,e,derb,q,q2,fi,rai_inv,vtot;

	if(PARTDECOMP(cr))
	{
		pd_at_range(cr,&at0,&at1);
	}
	else if(DOMAINDECOMP(cr))
	{
		at0=0;
		at1=cr->dd->nat_home;
	}
	else
	{
		at0=0;
		at1=natoms;
		
	}
			
	vtot=0.0;
	
	/* Apply self corrections */
	for(i=at0;i<at1;i++)
	{
		ai       = i;
		
		if(born->vs[ai]==1)
		{
			rai      = born->bRad[ai];
			rai_inv  = 1.0/rai;
			q        = charge[ai];
			q2       = q*q;
			fi       = facel*q2;
			e        = fi*rai_inv;
			derb     = 0.5*e*rai_inv*rai_inv;
			dvda[ai] += derb*rai;
			vtot     -= 0.5*e;
		}
	}
	
   return vtot;	
	
}

real calc_gb_nonpolar(t_commrec *cr, t_forcerec *fr,int natoms,gmx_genborn_t *born, gmx_localtop_t *top, 
					  const t_atomtypes *atype, real *dvda,int gb_algorithm, t_mdatoms *md)
{
	int ai,i,at0,at1;
	real e,es,rai,rbi,term,probe,tmp,factor;
	real rbi_inv,rbi_inv2;
	
	/* To keep the compiler happy */
	factor=0;
	
	if(PARTDECOMP(cr))
	{
		pd_at_range(cr,&at0,&at1);
	}
	else if(DOMAINDECOMP(cr))
	{
		at0 = 0;
		at1 = cr->dd->nat_home;
	}
	else
	{
		at0=0;
		at1=natoms;
	}
	
	/* The surface area factor is 0.0049 for Still model, 0.0054 for HCT/OBC */
	if(gb_algorithm==egbSTILL)
	{
		factor=0.0049*100*CAL2JOULE;
	}
	else	
	{
		factor=0.0054*100*CAL2JOULE;	
	}
	
	/* if(gb_algorithm==egbHCT || gb_algorithm==egbOBC) */
	
	es    = 0;
	probe = 0.14;
	term  = M_PI*4;
	
	for(i=at0;i<at1;i++)
	{
		ai        = i;
		
		if(born->vs[ai]==1)
		{
			rai		  = top->atomtypes.gb_radius[md->typeA[ai]];
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



real calc_gb_chainrule(int natoms, t_nblist *nl, real *dadx, real *dvda, rvec x[], rvec t[], 
					   int gb_algorithm, gmx_genborn_t *born)
{	
	int i,k,n,ai,aj,nj0,nj1;
	real fgb,fij,rb2,rbi,fix1,fiy1,fiz1;
	real ix1,iy1,iz1,jx1,jy1,jz1,dx11,dy11,dz11,rsq11;
	real rinv11,tx,ty,tz,rbai;
	real *rb;
	rvec dx;
	volatile int idx;
	
	n  = 0;	
	rb = born->work;
		
	/* Loop to get the proper form for the Born radius term */
	if(gb_algorithm==egbSTILL) 
	{
		for(i=0;i<natoms;i++)
		{
			rbi   = born->bRad[i];
			rb[i] = (2 * rbi * rbi * dvda[i])/ONE_4PI_EPS0;
		}
	}
	else if(gb_algorithm==egbHCT) 
	{
		for(i=0;i<natoms;i++)
		{
			rbi   = born->bRad[i];
			rb[i] = rbi * rbi * dvda[i];
		}
	}
	else if(gb_algorithm==egbOBC) 
	{
		for(idx=0;idx<natoms;idx++)
		{
			rbi   = born->bRad[idx];
			rb[idx] = rbi * rbi * born->drobc[idx] * dvda[idx];
		}
	}
	
	for(i=0;i<nl->nri;i++)
	{
		ai   = nl->iinr[i];
		nj0	 = nl->jindex[ai];
		nj1  = nl->jindex[ai+1];
		
		ix1  = x[ai][0];
		iy1  = x[ai][1];
		iz1  = x[ai][2];
		
		fix1 = 0;
		fiy1 = 0;
		fiz1 = 0;
		
		rbai = rb[ai];
		
		for(k=nj0;k<nj1;k++)
		{
			aj = nl->jjnr[k];
			
			jx1     = x[aj][0];
			jy1     = x[aj][1];
			jz1     = x[aj][2];
			
			dx11    = ix1 - jx1;
			dy11    = iy1 - jy1;
			dz11    = iz1 - jz1;
			
			fgb     = rbai*dadx[n++]; 
		
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
				
		/* Update force on atom ai */
		t[ai][0] = t[ai][0] + fix1;
		t[ai][1] = t[ai][1] + fiy1;
		t[ai][2] = t[ai][2] + fiz1;
		
	}
	
	return 0;	
}


real calc_gb_forces(t_commrec *cr, t_mdatoms *md, gmx_genborn_t *born, gmx_localtop_t *top, const t_atomtypes *atype, 
                    rvec x[], rvec f[], t_forcerec *fr, t_idef *idef, int gb_algorithm, bool bRad)
{
	real v=0;

	/* Do a simple ACE type approximation for the non-polar solvation */
	v += calc_gb_nonpolar(cr, fr,born->nr, born, top, atype, fr->dvda, gb_algorithm,md);
	
	/* Calculate the bonded GB-interactions */
	v += gb_bonds_tab(x[0],f[0],md->chargeA,&(fr->gbtabscale),
					  fr->invsqrta,fr->dvda,fr->gbtab.tab,idef,fr->epsilon_r, fr->epsfac);
	
	/* Calculate self corrections to the GB energies - currently only A state used! (FIXME) */
	v += calc_gb_selfcorrections(cr,born->nr,md->chargeA, born, fr->dvda, md, fr->epsfac); 		

	/* If parallel, sum the derivative of the potential w.r.t the born radii */
	if(PARTDECOMP(cr))
	{
		gmx_sum(md->nr,fr->dvda, cr);
	}
	else if(DOMAINDECOMP(cr))
	{
		dd_atom_sum_real(cr->dd,fr->dvda);
		dd_atom_spread_real(cr->dd,fr->dvda);
	}
	
#ifdef GMX_DOUBLE	
	
#if ( defined(GMX_IA32_SSE2) || defined(GMX_X86_64_SSE2) || defined(GMX_SSE2) )	
	calc_gb_chainrule_sse2_double(born->nr, &(fr->gblist), fr->dadx, fr->dvda, x[0], f[0], gb_algorithm, born);
#else
	calc_gb_chainrule(born->nr, &(fr->gblist), fr->dadx, fr->dvda, x, f, gb_algorithm, born);
#endif
	
#else
	
#if ( defined(GMX_IA32_SSE) || defined(GMX_X86_64_SSE) || defined(GMX_SSE2) )
	/* x86 or x86-64 with GCC inline assembly and/or SSE intrinsics */
	calc_gb_chainrule_sse(born->nr, &(fr->gblist), fr->dadx, fr->dvda, x[0], f[0], gb_algorithm, born);	
#else
	/* Calculate the forces due to chain rule terms with non sse code */
	calc_gb_chainrule(born->nr, &(fr->gblist), fr->dadx, fr->dvda, x, f, gb_algorithm, born);	
#endif	
#endif

	return v;

}

/* Skriv om den här rutinen så att den får samma format som Charmm SASA
 * Antagligen kan man då beräkna derivatan i samma loop som ytan,
 * vilket är snabbare (görs inte i Charmm, men det borde gå :-) )
 */
int calc_surfStill(t_inputrec *ir,
		  t_idef     *idef,
		  t_atoms    *atoms,
		  rvec       x[],
			rvec       f[],						 
			gmx_genborn_t     *born,
		  t_atomtypes *atype,
			double     *faction,
			int        natoms,
			t_nblist   *nl,
			t_iparams  forceparams[],
			t_iatom    forceatoms[],
			int        nbonds)
{
  int i,j,n,ia,ib,ic;
  
  real pc[3],radp,probd,radnp,s,prob,bmlt;
  real dx,dy,dz,d,rni,dist,bee,ri,rn,asurf,t1ij,t2ij,bij,bji,dbijp;
  real dbjip,tip,tmp,tij,tji,t3ij,t4ij,t5ij,t3ji,t4ji,t5ji,dbij,dbji;
  real dpf,dpx,dpy,dpz,pi,pn,si,sn;
  real *aprob;

  int factor=1;
  
	int k,type,ai,aj,nj0,nj1;
	real dr2,sar,rai,raj,fij;
	rvec dxx;
	
	aprob = born->work;
  /*bonds_t *bonds,*bonds13;*/
  
  /*snew(bonds,natoms);*/
  /*snew(bonds13,natoms);*/
   
  /* Zero out the forces to compare only surface area contribution */
  /*
  printf("Zeroing out forces before surface area..\n");
 
   for(i=0;i<natoms;i++)
    {
      
		if(!(is_hydrogen(*atoms->atomname[i])))
		{
			printf("x=%g, s=%g, p=%g, r=%g\n",
				x[i][0],
				atype->vol[atoms->atom[i].type],
				atype->surftens[atoms->atom[i].type],			
				factor*atype->radius[atoms->atom[i].type]);
	
		}
 
      printf("faction[i]=%g\n",faction[i*3]);
      faction[i*3]=0;
      faction[i*3+1]=0;
      faction[i*3+2]=0;
			
			f[i][0]=f[i][1]=f[i][2]=0;
    }
   exit(1);
	*/
  /* Radius of probe */
  radp=0.14;
  probd=2*radp;
  	
	/*********************************************************
	 *********************************************************
	 *********************************************************
	 *********************************************************
	 *********************************************************
	 Begin SA calculation Gromacs-style
	 *********************************************************
	 *********************************************************
	 *********************************************************
	 *********************************************************
	 *********************************************************/
	
	/* First set up the individual areas */
	for(n=0;n<natoms;n++)
	{	
		rn=atype->radius[atoms->atom[n].type];
		born->asurf[n]=(4*M_PI*(rn+radp)*(rn+radp));
	}
	
	/* Then loop over the bonded interactions */
	for(i=0;i<nbonds; )
	{
		type = forceatoms[i++];
		ai   = forceatoms[i++];
		aj   = forceatoms[i++];
		
		if(!is_hydrogen(*(atoms->atomname[ai])) &&
			 !is_hydrogen(*(atoms->atomname[aj])))
		{
		
			/*printf("genborn.c: xi=%g, xj=%g\n",factor*x[ai][0],factor*x[aj][0]);	*/
		
			rvec_sub(x[ai],x[aj],dxx);	
			dr2  = iprod(dxx,dxx);	
			
			sar  = forceparams[type].gb.sar;
			bmlt = forceparams[type].gb.bmlt;
			rni  = sar+probd;
			
			rn   = atype->radius[atoms->atom[ai].type];
			ri   = atype->radius[atoms->atom[aj].type];
			pn   = atype->surftens[atoms->atom[ai].type];
			pi   = atype->surftens[atoms->atom[aj].type];			
			/*sn   = atype->vol[atoms->atom[ai].type]; */
			/*si   = atype->vol[atoms->atom[aj].type]; */
			
			/*rni  = rn + ri +probd; */
			/*printf("genborn.c: rn=%g, ri=%g\n",ri,rn); */
			if(dr2<rni*rni)
			{
				/*printf("d=%g, s=%g\n",dr2,rni*rni);*/
				dist = dr2*invsqrt(dr2);
				t1ij = M_PI*(rni-dist);
				t2ij = (rn-ri)/dist;
				bij  = (rn+radp)*t1ij*(1-t2ij);
				bji  = (ri+radp)*t1ij*(1+t2ij);
				tij  = pn*bmlt*bij/(4*M_PI*(rn+radp)*(rn+radp));
				tji  = pi*bmlt*bji/(4*M_PI*(ri+radp)*(ri+radp));
				
				born->asurf[ai] = born->asurf[ai]*(1-tij);
				born->asurf[aj] = born->asurf[aj]*(1-tji);
			}
		}
	}
	
	/* Now loop over interactions >= 1-4 */
	bmlt=0.3516;
	
	/*printf("NONBONDED INTERACTIONS\n");*/
	
	for(i=0;i<natoms;i++)
	{
		ai    = i;
		nj0   = nl->jindex[ai];
		nj1   = nl->jindex[ai+1];
		
		rai   = factor*atype->radius[atoms->atom[ai].type];
		pn    = atype->surftens[atoms->atom[ai].type];
		/*sn    = atype->vol[atoms->atom[ai].type];*/
		
		for(k=nj0;k<nj1;k++)
		{
			aj  = nl->jjnr[k];
			
			if(!is_hydrogen(*(atoms->atomname[ai])) &&
					!is_hydrogen(*(atoms->atomname[aj])))
			{
				raj = factor*atype->radius[atoms->atom[aj].type];
				pi  = atype->surftens[atoms->atom[aj].type];
				/*si  = atype->vol[atoms->atom[aj].type];*/
				rvec_sub(x[ai],x[aj],dxx);
			
				dr2 = factor*factor*iprod(dxx,dxx);
				rni = rai + raj + probd;
				/*printf("genborn.c: rn=%g, ri=%g, sar=%g, dr2=%g\n",rai,raj,sar,dr2);	*/
				/*printf("genborn.c: xi=%g, xj=%g\n",factor*x[ai][0],factor*x[aj][0]); */
				
				if(dr2<rni*rni)
				{
					/*printf("d=%g, s=%g\n",dr2,rni*rni);	*/
					dist = dr2*invsqrt(dr2);
					t1ij = M_PI*(rni-dist);
					t2ij = (rai-raj)/dist;
					bij  = (rai+radp)*t1ij*(1-t2ij);
					bji  = (raj+radp)*t1ij*(1+t2ij);
					tij  = pn*bmlt*bij/(4*M_PI*(rai+radp)*(rai+radp));
					tji  = pi*bmlt*bji/(4*M_PI*(raj+radp)*(raj+radp));
									
					born->asurf[ai]=born->asurf[ai]*(1-tij);
					born->asurf[aj]=born->asurf[aj]*(1-tji);
				}
			}
		}
	}
	/*
	printf("AFTER BOTH AREA CALCULATIONS\n");
	n=0;
	for(i=0;i<natoms;i++)
	{
      if(!is_hydrogen(*atoms->atomname[i]))
			{	
				printf("%d, Still=%g, gromacs=%g\n",n,born->asurf[i], born->asurf[i]);
				
				born->as=born->as+born->asurf[i];
				born->as=born->as+born->asurf[i];
			}
		n++;
	}
	*/
	/*printf("Total Still area=%g, Total new area=%g\n",born->as, born->as);*/
	 /*printf("nbonds=%d\n",nbonds);*/
	/* Start to calculate the forces */
	for(i=0;i<nbonds; )
	{
		type = forceatoms[i++];
		ai   = forceatoms[i++];
		aj   = forceatoms[i++];
		
		if(!is_hydrogen(*(atoms->atomname[ai])) &&
			 !is_hydrogen(*(atoms->atomname[aj])))
		{
			rvec_sub(x[ai],x[aj],dxx);	
			
			dr2  = factor*factor*iprod(dxx,dxx);	
		
			sar  = factor*forceparams[type].gb.sar;
			bmlt = forceparams[type].gb.bmlt;
			rni  = sar+probd;
			
			rn   = factor*atype->radius[atoms->atom[ai].type];
			ri   = factor*atype->radius[atoms->atom[aj].type];
			pn   = atype->surftens[atoms->atom[ai].type];
			pi   = atype->surftens[atoms->atom[aj].type];			
			sn   = atype->vol[atoms->atom[ai].type];
			si   = atype->vol[atoms->atom[aj].type];
			
			if(dr2<rni*rni)
			{
				dist = dr2*invsqrt(dr2);
				t1ij = M_PI*(rni-dist);
				t2ij = (rn-ri)/dist;
				bij  = (rn+radp)*t1ij*(1-t2ij);
				bji  = (ri+radp)*t1ij*(1+t2ij);
				
				dbij = M_PI*(rn+radp)*(dr2-(rni*(rn-ri)));
				dbji = M_PI*(ri+radp)*(dr2+(rni*(rn-ri)));
				
				t3ij = sn*born->asurf[ai]*dbij;
				t4ij = (4*M_PI*(rn+radp)*(rn+radp))/(pn*bmlt)-bij;
				t5ij = t3ij/t4ij;
				
				t3ji = si*born->asurf[aj]*dbji;
				t4ji = (4*M_PI*(ri+radp)*(ri+radp))/(pi*bmlt)-bji;
				t5ji = t3ji/t4ji;
				
				dpf  = (t5ij+t5ji)/(dr2*dist);
				/*printf("deriv_cut: ai=%d, xi=%g aj=%d, xj=%g\n",ai,x[ai][0], aj,x[aj][0]);*/
				for(k=0;k<DIM;k++)
				{
					fij = factor*(-dpf)*dxx[k];
					f[ai][k]+=fij;
					f[aj][k]-=fij;
				}
			}
		}
	}
	
	/* Now calculate forces for all interactions >= 1-4 */
	bmlt = 0.3516;
	
	for(i=0;i<natoms;i++)
	{
		ai  = i;
		nj0 = nl->jindex[ai];
		nj1 = nl->jindex[ai+1];
		
		rai   = factor*atype->radius[atoms->atom[ai].type];
		pn    = atype->surftens[atoms->atom[ai].type];
		sn    = atype->vol[atoms->atom[ai].type];
		
		for(k=nj0;k<nj1;k++)
		{
			aj = nl->jjnr[k];
			
			if(!is_hydrogen(*(atoms->atomname[ai])) &&
				 !is_hydrogen(*(atoms->atomname[aj])))
			{
				raj = factor*atype->radius[atoms->atom[aj].type];
				pi  = atype->surftens[atoms->atom[aj].type];
				si  = atype->vol[atoms->atom[aj].type];
				
				rvec_sub(x[ai],x[aj],dxx);
				
				dr2 = factor*factor*iprod(dxx,dxx);
				rni = rai + raj + probd;
				
				if(dr2<rni*rni)
				{
					dist = dr2*invsqrt(dr2);
					t1ij = M_PI*(rni-dist);
					t2ij = (rai-raj)/dist;
					bij  = (rai+radp)*t1ij*(1-t2ij);
					bji  = (raj+radp)*t1ij*(1+t2ij);
					
					dbij = M_PI*(rai+radp)*(dr2-(rni*(rai-raj)));
					dbji = M_PI*(raj+radp)*(dr2+(rni*(rai-raj)));
					
					t3ij = sn*born->asurf[ai]*dbij;
					t4ij = (4*M_PI*(rai+radp)*(rai+radp))/(pn*bmlt)-bij;
					t5ij = t3ij/t4ij;
					
					t3ji = si*born->asurf[aj]*dbji;
					t4ji = (4*M_PI*(raj+radp)*(raj+radp))/(pi*bmlt)-bji;
					t5ji = t3ji/t4ji;
					
					dpf  = (t5ij+t5ji)/(dr2*dist);
					/*printf("deriv_cut: ai=%d, xi=%g aj=%d, xj=%g\n",ai,x[ai][0], aj,x[aj][0]);*/
					for(n=0;n<DIM;n++)
					{
						fij = factor*(-dpf)*dxx[n];
						f[ai][n]+=fij;
						f[aj][n]-=fij;
					}
				}
			}
		}
	}
	/*
	printf("AFTER BOTH FORCES CALCULATIONS\n");
	n=0;
	for(i=0;i<natoms;i++)
	{
		if(!is_hydrogen(*atoms->atomname[i]))
		{	
			printf("%d, gx=%g, gy=%g, gz=%g\n",
						 n,
						 faction[i*3], 
						 faction[i*3+1],
						 faction[i*3+2],
						 f[i][0],
						 f[i][1],
						 f[i][2]);
		}
		n++;
	}
	*/
  return 0;
}

int calc_surfBrooks(t_inputrec *ir,
		    t_idef     *idef,
		    t_atoms    *atoms,
		    rvec       x[],
		    gmx_genborn_t     *born,
		    t_atomtypes *atype,
		    double      *faction,
		    int natoms)
{
  int i,j,k;
  real xi,yi,zi,dx,dy,dz,ri,rj;
  real rho,rho2,rho6,r2,r,aij,aijsum,daij;
  real kappa,sassum,tx,ty,tz,fix1,fiy1,fiz1;

  real ck[5];
  real *Aij;
  real *sasi;
	
	Aij  = born->work;
	sasi = born->work;

  /* Brooks parameter for cutoff between atom pairs
   * Increasing kappa will increase the number of atom pairs
   * included in the calculation, which will also slow the calculation
   */
  kappa=0;

  sassum=0;

  /* Hydrogen atoms are included in via the ck parameters for the
   * heavy atoms
   */
  for(i=0;i<natoms;i++)
    {
      /*if(strncmp(*atoms->atomname[i],"H",1)!=0) */
      /*{ */
	  xi=x[i][0];
	  yi=x[i][1];
	  zi=x[i][2];
	  
	  fix1=0;
	  fiy1=0;
	  fiz1=0;

	  ri=atype->radius[atoms->atom[i].type];
	  aijsum=0;
	  
	  for(j=0;j<natoms;j++)
	    {
	      /*if(strncmp(*atoms->atomname[j],"H",1)!=0 && i!=j) */
	      if(i!=j)
	      {
		  dx=xi-x[j][0];
		  dy=yi-x[j][1];
		  dz=zi-x[j][2];
		  
		  r2=dx*dx+dy*dy+dz*dz;
		  r=sqrt(r2);
		  rj=atype->radius[atoms->atom[j].type];
		  
		  rho=ri+rj+kappa;
		  rho2=rho*rho;
		  rho6=rho2*rho2*rho2;

		    /* Cutoff test */
		    if(r<=rho)
		      {
			aij=pow((rho2-r2)*(rho2+2*r2),2)/rho6;
			daij=((4*r*(rho2-r2)*(rho2-r2))/rho6)-((4*r*(rho2-r2)*(2*r2+rho2))/rho6);
			tx=daij*dx;
			ty=daij*dy;
			tz=daij*dz;
			
			fix1=fix1+tx;
			fiy1=fiy1+ty;
			fiz1=fiz1+tz;

			faction[j*3]=faction[j*3]-tx;
			faction[j*3+1]=faction[j*3+1]-ty;
			faction[j*3+2]=faction[j*3+2]-tz;
			
			aijsum=aijsum+aij;
			printf("xi=%g, xj=%g, fscal=%g\n",xi,x[j][0],daij);
		      }
		}
	    }
	  
	  faction[i*3]=faction[i*3]+fix1;
	  faction[i*3+1]=faction[i*3+1]+fiy1;
	  faction[i*3+2]=faction[i*3+2]+fiz1;
	  
	  /* Calculate surface area coefficient */
	  Aij[i]=pow(aijsum,1/4);
	  
	  for(k=0;k<5;k++)
	    {
	      sasi[i]=sasi[i]+ck[k]*(pow(Aij[i],k));
	    }
	  
	  /* Increase total surface area */
	  sassum=sassum+sasi[i];
	  
	  /*}*/
    }

  printf("Brooks total surface area is: %g\n", sassum);


  return 0;
}

int make_gb_nblist(t_commrec *cr, int natoms, int gb_algorithm, real gbcut, rvec x[], 
				   t_forcerec *fr, t_idef *idef, gmx_genborn_t *born)
{
	int i,l,ii,j,k,n,nj0,nj1,ai,aj,idx,ii_idx,nalloc,at0,at1,found;
	int apa;
	t_nblist *nblist;

	int *count;
	int **atoms;
	
	count = born->count;
	atoms = born->nblist_work;
	
	memset(count,0,sizeof(int)*natoms);
	
	if(gb_algorithm==egbHCT || gb_algorithm==egbOBC)
	{
		/* Loop over 1-2, 1-3 and 1-4 interactions */
		for(j=F_GB12;j<=F_GB14;j++)
		{
			for(k=0;k<idef->il[j].nr;k+=3)
			{
				ai=idef->il[j].iatoms[k+1];
				aj=idef->il[j].iatoms[k+2];
				
				found=0;
				
				/* So that we do not add the same bond twice. This happens with some constraints between 1-3 atoms
				 * that are in the bond-list but should not be in the GB nb-list */
				for(i=0;i<count[ai];i++)
				{
					if(atoms[ai][i]==aj)
					{
						found=1;
					}
				}	
			 
				/* When doing HCT or OBC, we need to add all interactions to the nb-list twice 
				 * since the loop for calculating the Born-radii runs over all vs all atoms */	
				if(found==0)
				{
					atoms[ai][count[ai]]=aj;
					count[ai]++;
			
					atoms[aj][count[aj]]=ai;
					count[aj]++;
				}
			}
		}
	}
	
	if(gb_algorithm==egbSTILL)
	{
		/* Loop over 1-4 interactions */
		for(k=0;k<idef->il[F_GB14].nr;k+=3)
		{
			ai=idef->il[F_GB14].iatoms[k+1];
			aj=idef->il[F_GB14].iatoms[k+2];
		
			found=0;
			
			for(i=0;i<count[ai];i++)
			{
				if(atoms[ai][i]==aj)
				{
					found=1;
				}
			}	
			 
			/* Also for Still, we need to add (1-4) interactions twice */
			if(found==0)
			{
				atoms[ai][count[ai]]=aj;
				count[ai]++;
				
				atoms[aj][count[aj]]=ai;
				count[aj]++;
			}
		}
	}
	
	/* Loop over the VDWQQ and VDW nblists to set up the nonbonded part of the GB list */
	for(n=0; (n<fr->nnblists); n++)
	{
		for(i=0; (i<eNL_NR); i++)
		{
			nblist=&(fr->nblists[n].nlist_sr[i]);
			
			if(nblist->nri>0 && (i==eNL_VDWQQ || i==eNL_QQ))
			{
				for(j=0;j<nblist->nri;j++)
				{
					ai  = nblist->iinr[j];
			
					nj0 = nblist->jindex[j];
					nj1 = nblist->jindex[j+1];
					
					for(k=nj0;k<nj1;k++)
					{
						aj = nblist->jjnr[k];
						
						/* Allocate extra memory if needed */
						if(count[ai]>=born->nblist_work_nalloc || count[aj]>=born->nblist_work_nalloc)
						{
							born->nblist_work_nalloc += 500;
							for(apa=0;apa<natoms;apa++)
							{
								srenew(born->nblist_work[apa],born->nblist_work_nalloc);
							}
						}
												
						if(ai>aj)
						{
							atoms[aj][count[aj]]=ai;
							count[aj]++;
							
							/* We need to add all interactions to the nb-list twice 
							 * since the loop for calculating the Born-radii runs over all vs all atoms 
							 */
							atoms[ai][count[ai]]=aj;
							count[ai]++;
						}
						else
						{
							atoms[ai][count[ai]]=aj;
							count[ai]++;
							
							atoms[aj][count[aj]]=ai;
							count[aj]++;
						}
					}
				}
			}
		}
	}
		
	/* Zero out some counters */
	idx=0;
	ii_idx=0;
	fr->gblist.nri=0;
	fr->gblist.nrj=0;
	
	for(i=0;i<natoms;i++)
	{
		/* Only add those atoms that actually have neighbours (ie. all except vsites) */
		if(!PARTDECOMP(cr))
		{
			if(count[i]>0)
			{
				fr->gblist.iinr[fr->gblist.nri]=i;
				fr->gblist.nri++;
			}
		}
		else
		{
			/* For pd, some atoms can have all their neighbours on one node, and then count[ai] will be zero on the other
			 * nodes, but such ai atoms still need to be added to the nblist on all nodes to get correct summations
			 */
			fr->gblist.iinr[fr->gblist.nri]=i;
			fr->gblist.nri++;
		}
		
		for(k=0;k<count[i];k++)
		{
			/* Memory allocation for jjnr */
			if(fr->gblist.nrj>=fr->gblist.maxnrj)
			{
				fr->gblist.maxnrj += over_alloc_small(fr->gblist.maxnrj);
				
				if(debug)
				{
					fprintf(debug,"Increasing GB neighbourlist j size to %d\n",fr->gblist.maxnrj);
				}
				
				srenew(fr->gblist.jjnr,fr->gblist.maxnrj);
			
			}
			
			/* Put in list */
			fr->gblist.jjnr[idx++]=atoms[i][k];
			fr->gblist.nrj++;
		}
		
		fr->gblist.jindex[ii_idx+1]=idx;	
		ii_idx++;
	}
	
	return 0;
}

void make_local_gb(t_commrec *cr, gmx_genborn_t *born, int gb_algorithm)
{
	int i,at0,at1;
	gmx_domdec_t *dd=NULL;
	
	if(DOMAINDECOMP(cr))
	{
		dd = cr->dd;
		at0 = 0;
		at1 = dd->nat_tot;
	}
	else
	{
		/* Single node or particle decomp (global==local), just copy pointers and return */
		if(gb_algorithm==egbSTILL)
		{
			born->gpol  = born->gpol_globalindex;
			born->vsolv = born->vsolv_globalindex; 
		}
		else
		{
			born->param = born->param_globalindex;
		}
		
		born->vs = born->vs_globalindex;
		
		return;
	}
	
	/* Reallocation of local arrays if necessary */
	if(born->nlocal < dd->nat_tot)
	{
		born->nlocal = dd->nat_tot;
		
		/* Arrays specific to different gb algorithms */
		if(gb_algorithm==egbSTILL)
		{
			srenew(born->gpol,  born->nlocal+3);
			srenew(born->vsolv, born->nlocal+3);
		}
		else
		{
			srenew(born->param, born->nlocal+3);
		}
		
		/* All gb-algorithms use the array for vsites exclusions */
		srenew(born->vs,    born->nlocal+3);
	}
	
	/* With dd, copy algorithm specific arrays */
	if(gb_algorithm==egbSTILL)
	{
		for(i=at0;i<at1;i++)
		{
			born->gpol[i]  = born->gpol_globalindex[dd->gatindex[i]];
			born->vsolv[i] = born->vsolv_globalindex[dd->gatindex[i]];
			born->vs[i]    = born->vs_globalindex[dd->gatindex[i]];
		}
	}
	else
	{
		for(i=at0;i<at1;i++)
		{
			born->param[i] = born->param_globalindex[dd->gatindex[i]];
			born->vs[i]    = born->vs_globalindex[dd->gatindex[i]];
		}
	}
}

