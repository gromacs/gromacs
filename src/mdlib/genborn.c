/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
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
#include "pbc.h"

#ifdef GMX_LIB_MPI
#include <mpi.h>
#endif
#ifdef GMX_THREADS
#include "tmpi.h"
#endif

#if ( defined(GMX_IA32_SSE) || defined(GMX_X86_64_SSE) || defined(GMX_SSE2) )
#ifdef GMX_DOUBLE
#include "genborn_sse2_double.h"
#else
#include "genborn_sse2_single.h"
#endif /* GMX_DOUBLE */
#endif /* GMX_SSE */

/* Still parameters - make sure to edit in genborn_sse.c too if you change these! */
#define STILL_P1  0.073*0.1              /* length        */
#define STILL_P2  0.921*0.1*CAL2JOULE    /* energy*length */
#define STILL_P3  6.211*0.1*CAL2JOULE    /* energy*length */
#define STILL_P4  15.236*0.1*CAL2JOULE
#define STILL_P5  1.254 

#define STILL_P5INV (1.0/STILL_P5)
#define STILL_PIP5  (M_PI*STILL_P5)

typedef struct {
    int shift;
    int naj;
    int *aj;
    int aj_nalloc;
} gbtmpnbl_t;

typedef struct gbtmpnbls {
    int nlist;
    gbtmpnbl_t *list;
    int list_nalloc;
} t_gbtmpnbls;

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


static void assign_gb_param(t_functype ftype,t_iparams *new_par,
                            real old[MAXFORCEPARAM],int comb,real reppow)
{
    int  i,j;
    
    /* Set to zero */
    for(j=0; (j<MAXFORCEPARAM); j++)
    {
        new_par->generic.buf[j] = 0.0;
    }
  
    switch (ftype)
    {
    case F_GB12:
    case F_GB13:
    case F_GB14:
        new_par->gb.c6A=old[0];
        new_par->gb.c12A=old[1];
        new_par->gb.c6B=old[2];
        new_par->gb.c12B=old[3];
        new_par->gb.sar=old[4];
        new_par->gb.st=old[5];
        new_par->gb.pi=old[6];
        new_par->gb.gbr=old[7];
        new_par->gb.bmlt=old[8];
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
  t_iparams new_par;
  int       type;
	
	assign_gb_param(ftype,&new_par,forceparams,comb,reppow);
  if (!bAppend) {
		for (type=start; (type<ffparams->ntypes); type++) {
      if (ffparams->functype[type]==ftype) {
					if (memcmp(&new_par,&ffparams->iparams[type],(size_t)sizeof(new_par)) == 0)
					return type;
      }
    }
  }
  else
    type=ffparams->ntypes;
  if (debug)
    fprintf(debug,"copying new_par to idef->iparams[%d] (ntypes=%d)\n",
						type,ffparams->ntypes);
  memcpy(&ffparams->iparams[type],&new_par,(size_t)sizeof(new_par));
  
  ffparams->ntypes++;
  ffparams->functype[type]=ftype;
	
  return type;
}

int init_gb_still(const t_commrec *cr, t_forcerec  *fr, const t_atomtypes *atype, t_idef *idef, t_atoms *atoms, gmx_genborn_t *born,int natoms)
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
	snew(born->gpol_still_work,natoms+3);
	
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
		if(born->use_globalindex[j]==1)
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
int init_gb(gmx_genborn_t **p_born,
            const t_commrec *cr, t_forcerec *fr, const t_inputrec *ir,
			const gmx_mtop_t *mtop, real rgbradii, int gb_algorithm)
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
	snew(born->use_globalindex,    natoms+3);
	
	snew(fr->invsqrta, natoms);
	snew(fr->dvda,     natoms);
	
	fr->dadx              = NULL;
	fr->nalloc_dadx       = 0;
	born->gpol_still_work = NULL;
	born->gpol_hct_work   = NULL;
	
	/* snew(born->asurf,natoms); */
	/* snew(born->dasurf,natoms); */

	/* Initialize the gb neighbourlist */
	init_gb_nblist(natoms, &(fr->gblist));
	
	/* Do the Vsites exclusions (if any) */
	for(i=0;i<natoms;i++)
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
	    if (C6 (fr->nbfp,fr->ntype,jj,jj) == 0 &&
            C12(fr->nbfp,fr->ntype,jj,jj) == 0 &&
            atoms.atom[i].q == 0)
		{
			born->use_globalindex[i]=0;
		}
	}
	
	/* Copy algorithm parameters from inputrecord to local structure */
	born->obc_alpha  = ir->gb_obc_alpha;
	born->obc_beta   = ir->gb_obc_beta;
	born->obc_gamma  = ir->gb_obc_gamma;
	born->gb_doffset = ir->gb_dielectric_offset;
	born->gb_epsilon_solvent = ir->gb_epsilon_solvent;
	
	doffset = born->gb_doffset;
	
	/* If Still model, initialise the polarisation energies */
	if(gb_algorithm==egbSTILL)	
	{
		init_gb_still(cr, fr,&(mtop->atomtypes), &(localtop->idef), &atoms, born, natoms);	
	}
	/* If HCT/OBC,  precalculate the sk*atype->S_hct factors */
	else if(gb_algorithm==egbHCT || gb_algorithm==egbOBC)
	{
		
		snew(born->gpol_hct_work, natoms+3);
		
		for(i=0;i<natoms;i++)
		{	
			if(born->use_globalindex[i]==1)
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
	int shift;
	real shX,shY,shZ;
	real gpi,dr,dr2,dr4,idr4,rvdw,ratio,ccf,theta,term,rai,raj;
	real ix1,iy1,iz1,jx1,jy1,jz1,dx11,dy11,dz11;
	real rinv,idr2,idr6,vaj,dccf,cosq,sinq,prod,gpi2;
	real factor;
	real vai, prod_ai, icf4,icf6;
	
	factor  = 0.5*ONE_4PI_EPS0;
	n       = 0;
	
	for(i=0;i<born->nr;i++)
	{
		born->gpol_still_work[i]=0;
	}
		
	for(i=0;i<nl->nri;i++ )
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
			rinv  = gmx_invsqrt(dr2);
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
				dccf  = 2.0*term*sinq*gmx_invsqrt(sinq)*theta;
			}
			
			prod          = STILL_P4*vaj;
			icf4          = ccf*idr4;
			icf6          = (4*ccf-dccf)*idr6;

			born->gpol_still_work[aj] += prod_ai*icf4;
			gpi             = gpi+prod*icf4;
			
			/* Save ai->aj and aj->ai chain rule terms */
			fr->dadx[n++]   = prod*icf6;
			fr->dadx[n++]   = prod_ai*icf6;
		}
		
		born->gpol_still_work[ai]+=gpi;
		
	}
	
	/* Parallel summations */
	if(PARTDECOMP(cr))
	{
		gmx_sum(natoms, born->gpol_still_work, cr);
	}
	else if(DOMAINDECOMP(cr))
	{
		dd_atom_sum_real(cr->dd, born->gpol_still_work);
	}
	
	/* Calculate the radii */
	for(i=0;i<nl->nri;i++)
	{
		ai   = nl->iinr[i];
		gpi  = born->gpol[ai]+born->gpol_still_work[ai];
		gpi2 = gpi * gpi;
		born->bRad[ai]   = factor*gmx_invsqrt(gpi2);
		fr->invsqrta[ai] = gmx_invsqrt(born->bRad[ai]);
	}
	
	/* Extra communication required for DD */
	if(DOMAINDECOMP(cr))
	{
		dd_atom_spread_real(cr->dd, born->bRad);
		dd_atom_spread_real(cr->dd, fr->invsqrta);
	}
	
	return 0;
	
}
	

static int 
calc_gb_rad_hct(t_commrec *cr,t_forcerec *fr,int natoms, gmx_localtop_t *top,
				const t_atomtypes *atype, rvec x[], t_nblist *nl, gmx_genborn_t *born,t_mdatoms *md)
{
	int i,k,n,ai,aj,nj0,nj1,at0,at1;
	int shift;
	real shX,shY,shZ;
	real rai,raj,gpi,dr2,dr,sk,sk_ai,sk2,sk2_ai,lij,uij,diff2,tmp,sum_ai;
	real rad,min_rad,rinv,rai_inv;
	real ix1,iy1,iz1,jx1,jy1,jz1,dx11,dy11,dz11;
	real lij2, uij2, lij3, uij3, t1,t2,t3;
	real lij_inv,dlij,duij,sk2_rinv,prod,log_term;
	real doffset,raj_inv;
	
	doffset = born->gb_doffset;
	
	for(i=0;i<born->nr;i++)
	{
		born->gpol_hct_work[i] = 0;
	}
	
	/* Keep the compiler happy */
	n    = 0;
	prod = 0;
		
	for(i=0;i<nl->nri;i++)
	{
		ai     = nl->iinr[i];
			
		nj0    = nl->jindex[ai];			
		nj1    = nl->jindex[ai+1];
		
		/* Load shifts for this list */
		shift   = nl->shift[i];
		shX     = fr->shift_vec[shift][0];
		shY     = fr->shift_vec[shift][1];
		shZ     = fr->shift_vec[shift][2];
		
		rai     = top->atomtypes.gb_radius[md->typeA[ai]]-doffset; 
		rai_inv = 1.0/rai;
		
		sk_ai   = born->param[ai];
		sk2_ai  = sk_ai*sk_ai;
		
		/* Load atom i coordinates, add shift vectors */
		ix1     = shX + x[ai][0];
		iy1     = shY + x[ai][1];
		iz1     = shZ + x[ai][2];
		
		sum_ai  = 0;
		
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
			rinv  = gmx_invsqrt(dr2);
			dr    = rinv*dr2;
			
			sk    = born->param[aj];
			raj   = top->atomtypes.gb_radius[md->typeA[aj]]-doffset; 
			
			/* aj -> ai interaction */
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
				
				lij_inv  = gmx_invsqrt(lij2);
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
					
				t1      = 0.5*lij2 + prod*lij3 - 0.25*(lij*rinv+lij3*dr);
				t2      = -0.5*uij2 - 0.25*sk2_rinv*uij3 + 0.25*(uij*rinv+uij3*dr);
				t3      = 0.125*(1.0+sk2_rinv*rinv)*(-diff2)+0.25*log_term*rinv*rinv;
				
				fr->dadx[n++] = (dlij*t1+t2+t3)*rinv; /* rb2 is moved to chainrule	*/
				/* fr->dadx[n++] = (dlij*t1+duij*t2+t3)*rinv; */ /* rb2 is moved to chainrule	*/

				sum_ai += 0.5*tmp;
			}
			
			/* ai -> aj interaction */
			if(raj < dr + sk_ai)
			{
				lij     = 1.0/(dr-sk_ai);
				dlij    = 1.0;
				raj_inv = 1.0/raj;
				
				if(raj>dr-sk_ai)
				{
					lij = raj_inv;
					dlij = 0.0;
				}
				
				lij2     = lij  * lij;
				lij3     = lij2 * lij;
				
				uij      = 1.0/(dr+sk_ai);
				uij2     = uij  * uij;
				uij3     = uij2 * uij;
				
				diff2    = uij2-lij2;
				
				lij_inv  = gmx_invsqrt(lij2);
				sk2      =  sk2_ai; /* sk2_ai = sk_ai * sk_ai in i loop above */
				sk2_rinv = sk2*rinv;
				prod     = 0.25 * sk2_rinv;
				
				/* log_term = table_log(uij*lij_inv,born->log_table,LOG_TABLE_ACCURACY); */
				log_term = log(uij*lij_inv);
				
				tmp      = lij-uij + 0.25*dr*diff2 + (0.5*rinv)*log_term + prod*(-diff2);
				
				if(raj<sk_ai-dr)
				{
					tmp     = tmp + 2.0 * (raj_inv-lij);
				}
				
				/* duij = 1.0 */
				t1      = 0.5*lij2 + prod*lij3 - 0.25*(lij*rinv+lij3*dr);
				t2      = -0.5*uij2 - 0.25*sk2_rinv*uij3 + 0.25*(uij*rinv+uij3*dr);
				t3      = 0.125*(1.0+sk2_rinv*rinv)*(-diff2)+0.25*log_term*rinv*rinv;
				
				fr->dadx[n++] = (dlij*t1+t2+t3)*rinv; /* rb2 is moved to chainrule	*/
				/* fr->dadx[n++] = (dlij*t1+duij*t2+t3)*rinv; */ /* rb2 is moved to chainrule	*/
				
				born->gpol_hct_work[aj] += 0.5*tmp;
			}
		}
		
		born->gpol_hct_work[ai] += sum_ai;
	}
	
	/* Parallel summations */
	if(PARTDECOMP(cr))
	{
		gmx_sum(natoms, born->gpol_hct_work, cr);
	}
	else if(DOMAINDECOMP(cr))
	{
		dd_atom_sum_real(cr->dd, born->gpol_hct_work);
	}
	
	for(i=0;i<nl->nri;i++)
	{
		ai      = nl->iinr[i];
		rai     = top->atomtypes.gb_radius[md->typeA[ai]]-doffset; 
		sum_ai  = 1.0/rai - born->gpol_hct_work[ai];
		min_rad = rai + doffset;
		rad     = 1.0/sum_ai; 
		
		born->bRad[ai]   = rad > min_rad ? rad : min_rad;
		fr->invsqrta[ai] = gmx_invsqrt(born->bRad[ai]);
	}
	
	/* Extra communication required for DD */
	if(DOMAINDECOMP(cr))
	{
		dd_atom_spread_real(cr->dd, born->bRad);
		dd_atom_spread_real(cr->dd, fr->invsqrta);
	}
	
	
	return 0;
}

static int 
calc_gb_rad_obc(t_commrec *cr, t_forcerec *fr, int natoms, gmx_localtop_t *top,
					const t_atomtypes *atype, rvec x[], t_nblist *nl, gmx_genborn_t *born,t_mdatoms *md)
{
	int i,k,ai,aj,nj0,nj1,n,at0,at1;
	int shift;
	real shX,shY,shZ;
	real rai,raj,gpi,dr2,dr,sk,sk2,lij,uij,diff2,tmp,sum_ai;
	real rad, min_rad,sum_ai2,sum_ai3,tsum,tchain,rinv,rai_inv,lij_inv,rai_inv2;
	real log_term,prod,sk2_rinv,sk_ai,sk2_ai;
	real ix1,iy1,iz1,jx1,jy1,jz1,dx11,dy11,dz11;
	real lij2,uij2,lij3,uij3,dlij,duij,t1,t2,t3;
	real doffset,raj_inv;

	/* Keep the compiler happy */
	n    = 0;
	prod = 0;
	raj  = 0;
	
	doffset = born->gb_doffset;
	
	for(i=0;i<born->nr;i++)
	{
		born->gpol_hct_work[i] = 0;
	}
	
	for(i=0;i<nl->nri;i++)
	{
		ai      = nl->iinr[i];
	
		nj0     = nl->jindex[i];
		nj1     = nl->jindex[i+1];
		
		/* Load shifts for this list */
		shift   = nl->shift[i];
		shX     = fr->shift_vec[shift][0];
		shY     = fr->shift_vec[shift][1];
		shZ     = fr->shift_vec[shift][2];
		
		rai      = top->atomtypes.gb_radius[md->typeA[ai]]-doffset;
		rai_inv  = 1.0/rai;
		rai_inv2 = 1.0/top->atomtypes.gb_radius[md->typeA[ai]];
		
		sk_ai    = born->param[ai];
		sk2_ai   = sk_ai*sk_ai;
		
		/* Load atom i coordinates, add shift vectors */
		ix1      = shX + x[ai][0];
		iy1      = shY + x[ai][1];
		iz1      = shZ + x[ai][2];
		
		sum_ai   = 0;
		
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
			rinv  = gmx_invsqrt(dr2);
			dr    = dr2*rinv;
		
			/* sk is precalculated in init_gb() */
			sk    = born->param[aj];
			raj   = top->atomtypes.gb_radius[md->typeA[aj]]-doffset; 
			
			/* aj -> ai interaction */
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
				
				lij_inv  = gmx_invsqrt(lij2);
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
				
				sum_ai += 0.5*tmp;
			}
				
			/* ai -> aj interaction */
			if(raj < dr + sk_ai)
			{
				lij     = 1.0/(dr-sk_ai);
				dlij    = 1.0;
				raj_inv = 1.0/raj;
				
				if(raj>dr-sk_ai)
				{
					lij = raj_inv;
					dlij = 0.0;
				}
				
				lij2     = lij  * lij;
				lij3     = lij2 * lij;
				
				uij      = 1.0/(dr+sk_ai);
				uij2     = uij  * uij;
				uij3     = uij2 * uij;
				
				diff2    = uij2-lij2;
				
				lij_inv  = gmx_invsqrt(lij2);
				sk2      =  sk2_ai; /* sk2_ai = sk_ai * sk_ai in i loop above */
				sk2_rinv = sk2*rinv;
				prod     = 0.25 * sk2_rinv;
				
				/* log_term = table_log(uij*lij_inv,born->log_table,LOG_TABLE_ACCURACY); */
				log_term = log(uij*lij_inv);
				
				tmp      = lij-uij + 0.25*dr*diff2 + (0.5*rinv)*log_term + prod*(-diff2);
				
				if(raj<sk_ai-dr)
				{
					tmp     = tmp + 2.0 * (raj_inv-lij);
				}
				
				t1      = 0.5*lij2 + prod*lij3 - 0.25*(lij*rinv+lij3*dr);
				t2      = -0.5*uij2 - 0.25*sk2_rinv*uij3 + 0.25*(uij*rinv+uij3*dr);
				t3      = 0.125*(1.0+sk2_rinv*rinv)*(-diff2)+0.25*log_term*rinv*rinv;
				
				fr->dadx[n++] = (dlij*t1+t2+t3)*rinv; /* rb2 is moved to chainrule	*/
				
				born->gpol_hct_work[aj] += 0.5*tmp;
				
			}
		}
		
		born->gpol_hct_work[ai] += sum_ai;
	}
	
	/* Parallel summations */
	if(PARTDECOMP(cr))
	{
		gmx_sum(natoms, born->gpol_hct_work, cr);
	}
	else if(DOMAINDECOMP(cr))
	{
		dd_atom_sum_real(cr->dd, born->gpol_hct_work);
	}
	
	for(i=0;i<nl->nri;i++)
	{
		ai         = nl->iinr[i];
		rai        = top->atomtypes.gb_radius[md->typeA[ai]];
		rai_inv2   = 1.0/rai;
		rai        = rai-doffset; 
		rai_inv    = 1.0/rai;
		sum_ai     = rai * born->gpol_hct_work[ai];
		sum_ai2    = sum_ai  * sum_ai;
		sum_ai3    = sum_ai2 * sum_ai;
		
		tsum    = tanh(born->obc_alpha*sum_ai-born->obc_beta*sum_ai2+born->obc_gamma*sum_ai3);
		born->bRad[ai] = rai_inv - tsum*rai_inv2;
		born->bRad[ai] = 1.0 / born->bRad[ai];
		
		fr->invsqrta[ai]=gmx_invsqrt(born->bRad[ai]);
		
		tchain  = rai * (born->obc_alpha-2*born->obc_beta*sum_ai+3*born->obc_gamma*sum_ai2);
		born->drobc[ai] = (1.0-tsum*tsum)*tchain*rai_inv2;
	}
	
	/* Extra (local) communication required for DD */
	if(DOMAINDECOMP(cr))
	{
		dd_atom_spread_real(cr->dd, born->bRad);
		dd_atom_spread_real(cr->dd, fr->invsqrta);
		dd_atom_spread_real(cr->dd, born->drobc);
	}
	
	return 0;
	
}



int calc_gb_rad(t_commrec *cr, t_forcerec *fr, t_inputrec *ir,gmx_localtop_t *top,
				const t_atomtypes *atype, rvec x[], t_nblist *nl, gmx_genborn_t *born,t_mdatoms *md)
{	
    /* First, reallocate the dadx array, we need at most 3 extra for SSE */
    if (2*nl->nrj + 3 > fr->nalloc_dadx)
    {
        fr->nalloc_dadx = over_alloc_large(2*nl->nrj + 3);
        srenew(fr->dadx,fr->nalloc_dadx);
    }
		
	/* Switch for determining which algorithm to use for Born radii calculation */
#ifdef GMX_DOUBLE
	
#if ( defined(GMX_IA32_SSE2) || defined(GMX_X86_64_SSE2) || defined(GMX_SSE2) )
	/* x86 or x86-64 with GCC inline assembly and/or SSE intrinsics */
	switch(ir->gb_algorithm)
	{
		case egbSTILL:
			 calc_gb_rad_still_sse2_double(cr,fr,md->nr,top, atype, x[0], nl, born, md); 
			break;
		case egbHCT:
			 calc_gb_rad_hct_sse2_double(cr,fr,md->nr,top, atype, x[0], nl, born, md); 
			break;
		case egbOBC:
			calc_gb_rad_obc_sse2_double(cr,fr,md->nr,top, atype, x[0], nl, born, md); 
			break;
			
		default:
			gmx_fatal(FARGS, "Unknown double precision sse-enabled algorithm for Born radii calculation: %d",ir->gb_algorithm);
	}
	
#else
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
				  real gb_epsilon_solvent, real facel)
{
	int i,j,n0,nnn,type,ai,aj,ai3,aj3;
	real isai,isaj;
	real r,rsq11,ix1,iy1,iz1,jx1,jy1,jz1;
	real dx11,dy11,dz11,rinv11,iq,facel2;
	real isaprod,qq,gbscale,gbtabscale,Y,F,Geps,Heps2,Fp,VV,FF,rt,eps,eps2;
	real vgb,fgb,vcoul,fijC,dvdatmp,fscal,tx,ty,tz,dvdaj;
	real vctot;
	
	t_iatom *forceatoms;

	/* Scale the electrostatics by gb_epsilon_solvent */
	facel = facel * (1.0 - 1.0/gb_epsilon_solvent);
	
	gbtabscale=*p_gbtabscale;
	vctot = 0.0;
	
	for(j=F_GB12;j<=F_GB14;j++)
	{
		forceatoms = idef->il[j].iatoms;
		
		for(i=0;i<idef->il[j].nr; )
		{
			/* To avoid reading in the interaction type, we just increment i to pass over
			 * the types in the forceatoms array, this saves some memory accesses
			 */
			i++;
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
			rinv11        = gmx_invsqrt(rsq11);
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
		
	/* Scale the electrostatics by gb_epsilon_solvent */
	facel = facel * (1.0 - 1.0/born->gb_epsilon_solvent);
	
	vtot=0.0;
	
	/* Apply self corrections */
	for(i=at0;i<at1;i++)
	{
		ai       = i;
		
		if(born->use[ai]==1)
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
		
		if(born->use[ai]==1)
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



real calc_gb_chainrule(int natoms, t_nblist *nl, real *dadx, real *dvda, rvec x[], rvec t[], rvec fshift[], 
					   rvec shift_vec[], int gb_algorithm, gmx_genborn_t *born)
{	
	int i,k,n,ai,aj,nj0,nj1;
	int shift;
	real shX,shY,shZ;
	real fgb,fij,rb2,rbi,fix1,fiy1,fiz1;
	real ix1,iy1,iz1,jx1,jy1,jz1,dx11,dy11,dz11,rsq11;
	real rinv11,tx,ty,tz,rbai,rbaj,fgb_ai;
	real *rb;
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
		
		for(k=nj0;k<nj1;k++)
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
		
		fshift[ai][0] = fshift[ai][0] + fix1;
		fshift[ai][1] = fshift[ai][1] + fiy1;
		fshift[ai][2] = fshift[ai][2] + fiz1;
		
	}

	return 0;	
}


real calc_gb_forces(t_commrec *cr, t_mdatoms *md, gmx_genborn_t *born, gmx_localtop_t *top, const t_atomtypes *atype, 
                    rvec x[], rvec f[], t_forcerec *fr, t_idef *idef, int gb_algorithm, bool bRad)
{
	real v=0;

	/* Do a simple ACE type approximation for the non-polar solvation */
	v += calc_gb_nonpolar(cr, fr,born->nr, born, top, atype, fr->dvda, gb_algorithm,md);
		
	/* Calculate the bonded GB-interactions using either table or analytical formula */
#ifdef GMX_DOUBLE	
	v += gb_bonds_tab(x[0],f[0],md->chargeA,&(fr->gbtabscale),
					  fr->invsqrta,fr->dvda,fr->gbtab.tab,idef,born->gb_epsilon_solvent, fr->epsfac);
#else	
#if ( defined(GMX_IA32_SSE2) || defined(GMX_X86_64_SSE2) || defined(GMX_SSE2) )	
	v += gb_bonds_analytic(x[0],f[0],md->chargeA,born->bRad,fr->dvda,idef,born->gb_epsilon_solvent,fr->epsfac);
#else
	v += gb_bonds_tab(x[0],f[0],md->chargeA,&(fr->gbtabscale),
					  fr->invsqrta,fr->dvda,fr->gbtab.tab,idef,born->gb_epsilon_solvent, fr->epsfac);
#endif
#endif
	
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
	 calc_gb_chainrule_sse2_double(born->nr, &(fr->gblist), fr->dadx, fr->dvda, 
								   x[0], f[0], fr->fshift[0],  fr->shift_vec[0],
								   gb_algorithm, born); 
#else
	calc_gb_chainrule(born->nr, &(fr->gblist), fr->dadx, fr->dvda, 
					  x, f, fr->fshift, fr->shift_vec, 
					  gb_algorithm, born);
#endif
	
#else
	
#if ( defined(GMX_IA32_SSE) || defined(GMX_X86_64_SSE) || defined(GMX_SSE2) )
	/* x86 or x86-64 with GCC inline assembly and/or SSE intrinsics */
	calc_gb_chainrule_sse(born->nr, &(fr->gblist), fr->dadx, fr->dvda, 
						  x[0], f[0], fr->fshift[0], fr->shift_vec[0], 
						  gb_algorithm, born);	
#else
	/* Calculate the forces due to chain rule terms with non sse code */
	calc_gb_chainrule(born->nr, &(fr->gblist), fr->dadx, fr->dvda, 
					  x, f, fr->fshift, fr->shift_vec, 
					  gb_algorithm, born);	
#endif	
#endif

	return v;

}

static void add_j_to_gblist(gbtmpnbl_t *list,int aj)
{
    if (list->naj >= list->aj_nalloc)
    {
        list->aj_nalloc = over_alloc_large(list->naj+1);
        srenew(list->aj,list->aj_nalloc);
    }

    list->aj[list->naj++] = aj;
}

static gbtmpnbl_t *find_gbtmplist(struct gbtmpnbls *lists,int shift)
{
    int ind;

    /* Search the list with the same shift, if there is one */
    ind = 0;
    while (ind < lists->nlist && shift != lists->list[ind].shift)
    {
        ind++;
    }
    if (ind == lists->nlist)
    {
        if (lists->nlist >= lists->list_nalloc)
        {
            lists->list_nalloc++;
            srenew(lists->list,lists->list_nalloc);
        }
        
        lists->list[lists->nlist].shift = shift;
        lists->list[lists->nlist].naj   = 0;
        lists->list[lists->nlist].aj    = NULL;
        lists->list[lists->nlist].aj_nalloc = 0;
        lists->nlist++;
    }

    return &lists->list[ind];
}

static void add_bondeds_to_gblist(t_ilist *il,
                                  bool bMolPBC,t_pbc *pbc,t_graph *g,rvec *x,
                                  struct gbtmpnbls *nls)
{
    int  ind,j,ai,aj,shift,found;
    rvec dx;
    ivec dt;
    gbtmpnbl_t *list;

    shift = CENTRAL;
    for(ind=0; ind<il->nr; ind+=3)
    {
        ai = il->iatoms[ind+1];
        aj = il->iatoms[ind+2];
				
        shift = CENTRAL;
        if (g != NULL)
        {
	      rvec_sub(x[ai],x[aj],dx);
	      ivec_sub(SHIFT_IVEC(g,ai),SHIFT_IVEC(g,aj),dt);
	      shift = IVEC2IS(dt);
	    }
        else if (bMolPBC)
        {
	      shift = pbc_dx_aiuc(pbc,x[ai],x[aj],dx);
        }

        /* Find the list for this shift or create one */
        list = find_gbtmplist(&nls[ai],shift);
        
        found=0;
		
        /* So that we do not add the same bond twice.
         * This happens with some constraints between 1-3 atoms
         * that are in the bond-list but should not be in the GB nb-list */
        for(j=0;j<list->naj;j++)
        {
            if (list->aj[j] == aj)
            {
                found = 1;
            }
        }	
		
        if (found == 0)
        {
            add_j_to_gblist(list,aj);
        }
    }
}

int make_gb_nblist(t_commrec *cr, int natoms, int gb_algorithm, real gbcut,
                   rvec x[], matrix box,
				   t_forcerec *fr, t_idef *idef, t_graph *graph, gmx_genborn_t *born)
{
	int i,l,ii,j,k,n,nj0,nj1,ai,aj,ii_idx,nalloc,at0,at1,found,shift,s;
	int apa;
	t_nblist *nblist;
    t_pbc pbc;

	struct gbtmpnbls *nls;
    gbtmpnbl_t *list;
	
	nls   = born->nblist_work;
	
	for(i=0;i<born->nr;i++)
	{
        nls[i].nlist = 0;
	}

    if (fr->bMolPBC)
    {
        set_pbc_dd(&pbc,fr->ePBC,cr->dd,TRUE,box);
    }

    switch (gb_algorithm)
    {
    case egbHCT:
    case egbOBC:
		/* Loop over 1-2, 1-3 and 1-4 interactions */
		for(j=F_GB12;j<=F_GB14;j++)
		{
            add_bondeds_to_gblist(&idef->il[j],fr->bMolPBC,&pbc,graph,x,nls);
		}
        break;
    case egbSTILL:
		/* Loop over 1-4 interactions */
        add_bondeds_to_gblist(&idef->il[F_GB14],fr->bMolPBC,&pbc,graph,x,nls);
        break;
    default:
        gmx_incons("Unknown GB algorithm");
	}
	
	/* Loop over the VDWQQ and VDW nblists to set up the nonbonded part of the GB list */
	for(n=0; (n<fr->nnblists); n++)
	{
		for(i=0; (i<eNL_NR); i++)
		{
			nblist=&(fr->nblists[n].nlist_sr[i]);
			
			if (nblist->nri > 0 && (i==eNL_VDWQQ || i==eNL_QQ))
			{
				for(j=0;j<nblist->nri;j++)
				{
					ai    = nblist->iinr[j];
                    shift = nblist->shift[j];

                    /* Find the list for this shift or create one */
                    list = find_gbtmplist(&nls[ai],shift);

					nj0 = nblist->jindex[j];
					nj1 = nblist->jindex[j+1];
					
                    /* Add all the j-atoms in the non-bonded list to the GB list */
					for(k=nj0;k<nj1;k++)
					{
                        add_j_to_gblist(list,nblist->jjnr[k]);
					}
				}
			}
		}
	}
		
	/* Zero out some counters */
	ii_idx=0;
	fr->gblist.nri=0;
	fr->gblist.nrj=0;
	
	if(DOMAINDECOMP(cr))
	{
		natoms = cr->dd->nat_home;
	}
	
    fr->gblist.jindex[0] = ii_idx;
    for(i=0;i<natoms;i++)
	{
        for(s=0; s<nls[i].nlist; s++)
        {
            list = &nls[i].list[s];

            /* Only add those atoms that actually have neighbours */
            if (born->use[i] != 0)
            {
                fr->gblist.iinr[fr->gblist.nri]  = i;
                fr->gblist.shift[fr->gblist.nri] = list->shift;
                fr->gblist.nri++;
            
                for(k=0; k<list->naj; k++)
                {
                    /* Memory allocation for jjnr */
                    if(fr->gblist.nrj >= fr->gblist.maxnrj)
                    {
                        fr->gblist.maxnrj += over_alloc_large(fr->gblist.maxnrj);
                        
                        if (debug)
                        {
                            fprintf(debug,"Increasing GB neighbourlist j size to %d\n",fr->gblist.maxnrj);
                        }
                        
                        srenew(fr->gblist.jjnr,fr->gblist.maxnrj);
                    }
			
                    /* Put in list */
                    fr->gblist.jjnr[fr->gblist.nrj++] = list->aj[k];
                }
            }
            
            ii_idx++;
            fr->gblist.jindex[ii_idx] = fr->gblist.nrj;	
        }
        if (nls[i].nlist == 0)
        {
            /* Temporary code adding an empty list to make loops work */
            fr->gblist.iinr[fr->gblist.nri]  = i;
            fr->gblist.shift[fr->gblist.nri] = list->shift;
            fr->gblist.nri++;

            ii_idx++;
            fr->gblist.jindex[ii_idx] = fr->gblist.nrj;	
        }
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
		
		born->use = born->use_globalindex;
		
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
		srenew(born->use,    born->nlocal+3);
	}
	
	/* With dd, copy algorithm specific arrays */
	if(gb_algorithm==egbSTILL)
	{
		for(i=at0;i<at1;i++)
		{
			born->gpol[i]  = born->gpol_globalindex[dd->gatindex[i]];
			born->vsolv[i] = born->vsolv_globalindex[dd->gatindex[i]];
			born->use[i]   = born->use_globalindex[dd->gatindex[i]];
		}
	}
	else
	{
		for(i=at0;i<at1;i++)
		{
			born->param[i] = born->param_globalindex[dd->gatindex[i]];
			born->use[i]   = born->use_globalindex[dd->gatindex[i]];
		}
	}
}

