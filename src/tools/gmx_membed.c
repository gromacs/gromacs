/*
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
 * Gallium Rubidium Oxygen Manganese Argon Carbon Silicon
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <signal.h>
#include <stdlib.h>
#include "typedefs.h"
#include "smalloc.h"
#include "sysstuff.h"
#include "vec.h"
#include "statutil.h"
#include "macros.h"
#include "copyrite.h"
#include "main.h"
#include "futil.h"
#include "edsam.h"
#include "checkpoint.h"
#include "vcm.h"
#include "mdebin.h"
#include "nrnb.h"
#include "calcmu.h"
#include "index.h"
#include "vsite.h"
#include "update.h"
#include "ns.h"
#include "trnio.h"
#include "xtcio.h"
#include "mdrun.h"
#include "confio.h"
#include "network.h"
#include "pull.h"
#include "xvgr.h"
#include "physics.h"
#include "names.h"
#include "disre.h"
#include "orires.h"
#include "dihre.h"
#include "pppm.h"
#include "pme.h"
#include "mdatoms.h"
#include "qmmm.h"
#include "mpelogging.h"
#include "domdec.h"
#include "partdec.h"
#include "topsort.h"
#include "coulomb.h"
#include "constr.h"
#include "shellfc.h"
#include "mvdata.h"
#include "checkpoint.h"
#include "mtop_util.h"
#include "tpxio.h"
#include "string2.h"
#include "sighandler.h"
#include "gmx_ana.h"

#ifdef GMX_LIB_MPI
#include <mpi.h>
#endif
#ifdef GMX_THREADS
#include "tmpi.h"
#endif

/* afm stuf */
#include "pull.h"

/* We use the same defines as in mvdata.c here */
#define  block_bc(cr,   d) gmx_bcast(     sizeof(d),     &(d),(cr))
#define nblock_bc(cr,nr,d) gmx_bcast((nr)*sizeof((d)[0]), (d),(cr))
#define   snew_bc(cr,d,nr) { if (!MASTER(cr)) snew((d),(nr)); }

/* The following two variables and the signal_handler function
 * is used from pme.c as well
 */

typedef struct {
	t_state s;
	rvec    *f;
	real    epot;
	real    fnorm;
	real    fmax;
	int     a_fmax;
} em_state_t;

typedef struct {
	int    it_xy;
	int    it_z;
	int    xy_step;
	int    z_step;
	rvec    xmin;
	rvec	xmax;
	rvec	*geom_cent;
	int    pieces;
	int    *nidx;
	atom_id **subindex;
} pos_ins_t;

typedef struct {
	int		id;
	char	*name;
	int 	nr;
	int 	natoms;	    /*nr of atoms per lipid*/
	int	mol1;	    /*id of the first lipid molecule*/
	real 	area;
} lip_t;

typedef struct {
	char	*name;
	t_block mem_at;
	int		*mol_id;
	int		nmol;
	real	lip_area;
	real	zmin;
	real	zmax;
	real	zmed;
} mem_t;

typedef struct {
	int		*mol;
	int		*block;
	int 	nr;
} rmm_t;

int search_string(char *s,int ng,char ***gn)
{
	int i;

	for(i=0; (i<ng); i++)
		if (gmx_strcasecmp(s,*gn[i]) == 0)
			return i;

	gmx_fatal(FARGS,"Group %s not found in indexfile.\nMaybe you have non-default groups in your mdp file, while not using the '-n' option of grompp.\nIn that case use the '-n' option.\n",s);

	return -1;
}

int get_mol_id(int at,int nmblock,gmx_molblock_t *mblock, int *type, int *block)
{
	int mol_id=0;
	int i;

	for(i=0;i<nmblock;i++)
	{
		if(at<(mblock[i].nmol*mblock[i].natoms_mol))
		{
			mol_id+=at/mblock[i].natoms_mol;
			*type = mblock[i].type;
			*block = i;
			return mol_id;
		} else {
			at-= mblock[i].nmol*mblock[i].natoms_mol;
			mol_id+=mblock[i].nmol;
		}
	}

	gmx_fatal(FARGS,"Something is wrong in mol ids, at %d, mol_id %d",at,mol_id);

	return -1;
}

int get_block(int mol_id,int nmblock,gmx_molblock_t *mblock)
{
	int i;
	int nmol=0;

	for(i=0;i<nmblock;i++)
	{
		nmol+=mblock[i].nmol;
		if(mol_id<nmol)
			return i;
	}

	gmx_fatal(FARGS,"mol_id %d larger than total number of molecules %d.\n",mol_id,nmol);

	return -1;
}

int get_tpr_version(const char *infile)
{
	char  	buf[STRLEN];
	gmx_bool  	bDouble;
	int 	precision,fver;
        t_fileio *fio;

	fio = open_tpx(infile,"r");
	gmx_fio_checktype(fio);

	precision = sizeof(real);

	gmx_fio_do_string(fio,buf);
	if (strncmp(buf,"VERSION",7))
		gmx_fatal(FARGS,"Can not read file %s,\n"
				"             this file is from a Gromacs version which is older than 2.0\n"
				"             Make a new one with grompp or use a gro or pdb file, if possible",
				gmx_fio_getname(fio));
	gmx_fio_do_int(fio,precision);
	bDouble = (precision == sizeof(double));
	if ((precision != sizeof(float)) && !bDouble)
		gmx_fatal(FARGS,"Unknown precision in file %s: real is %d bytes "
				"instead of %d or %d",
				gmx_fio_getname(fio),precision,sizeof(float),sizeof(double));
	gmx_fio_setprecision(fio,bDouble);
	fprintf(stderr,"Reading file %s, %s (%s precision)\n",
			gmx_fio_getname(fio),buf,bDouble ? "double" : "single");

	gmx_fio_do_int(fio,fver);

	close_tpx(fio);

	return fver;
}

void set_inbox(int natom, rvec *x)
{
	rvec tmp;
	int  i;

	tmp[XX]=tmp[YY]=tmp[ZZ]=0.0;
	for(i=0;i<natom;i++)
	{
		if(x[i][XX]<tmp[XX])		tmp[XX]=x[i][XX];
		if(x[i][YY]<tmp[YY])		tmp[YY]=x[i][YY];
		if(x[i][ZZ]<tmp[ZZ])		tmp[ZZ]=x[i][ZZ];
	}

	for(i=0;i<natom;i++)
			rvec_inc(x[i],tmp);
}

int get_mtype_list(t_block *at, gmx_mtop_t *mtop, t_block *tlist)
{
	int i,j,nr,mol_id;
        int type=0,block=0;
	gmx_bool bNEW;

	nr=0;
	snew(tlist->index,at->nr);
	for (i=0;i<at->nr;i++)
	{
		bNEW=TRUE;
		mol_id = get_mol_id(at->index[i],mtop->nmolblock,mtop->molblock,&type,&block);
		for(j=0;j<nr;j++)
		{
			if(tlist->index[j]==type)
						bNEW=FALSE;
		}
		if(bNEW==TRUE)
		{
			tlist->index[nr]=type;
			nr++;
		}
	}

	srenew(tlist->index,nr);
	return nr;
}

void check_types(t_block *ins_at,t_block *rest_at,gmx_mtop_t *mtop)
{
	t_block		*ins_mtype,*rest_mtype;
	int			i,j;

	snew(ins_mtype,1);
	snew(rest_mtype,1);
    ins_mtype->nr  = get_mtype_list(ins_at , mtop, ins_mtype );
    rest_mtype->nr = get_mtype_list(rest_at, mtop, rest_mtype);

    for(i=0;i<ins_mtype->nr;i++)
    {
    	for(j=0;j<rest_mtype->nr;j++)
    	{
    		if(ins_mtype->index[i]==rest_mtype->index[j])
    			gmx_fatal(FARGS,"Moleculetype %s is found both in the group to insert and the rest of the system.\n"
    					"Because we need to exclude all interactions between the atoms in the group to\n"
    					"insert, the same moleculetype can not be used in both groups. Change the\n"
    					"moleculetype of the molecules %s in the inserted group. Do not forget to provide\n"
    					"an appropriate *.itp file",*(mtop->moltype[rest_mtype->index[j]].name),
    					*(mtop->moltype[rest_mtype->index[j]].name));
    	}
    }

    sfree(ins_mtype->index);
    sfree(rest_mtype->index);
    sfree(ins_mtype);
    sfree(rest_mtype);
}

int init_ins_at(t_block *ins_at,t_block *rest_at,t_state *state, pos_ins_t *pos_ins,gmx_groups_t *groups,int ins_grp_id, real xy_max)
{
	int i,gid,c=0;
	real x,xmin,xmax,y,ymin,ymax,z,zmin,zmax;

	snew(rest_at->index,state->natoms);

	xmin=xmax=state->x[ins_at->index[0]][XX];
	ymin=ymax=state->x[ins_at->index[0]][YY];
	zmin=zmax=state->x[ins_at->index[0]][ZZ];

	for(i=0;i<state->natoms;i++)
	{
		gid = groups->grpnr[egcFREEZE][i];
		if(groups->grps[egcFREEZE].nm_ind[gid]==ins_grp_id)
		{
			x=state->x[i][XX];
			if (x<xmin) 			xmin=x;
			if (x>xmax)  			xmax=x;
			y=state->x[i][YY];
			if (y<ymin)				ymin=y;
			if (y>ymax)				ymax=y;
			z=state->x[i][ZZ];
			if (z<zmin)				zmin=z;
			if (z>zmax)				zmax=z;
		} else {
			rest_at->index[c]=i;
			c++;
		}
	}

	rest_at->nr=c;
	srenew(rest_at->index,c);

	if(xy_max>1.000001)
	{
		pos_ins->xmin[XX]=xmin-((xmax-xmin)*xy_max-(xmax-xmin))/2;
		pos_ins->xmin[YY]=ymin-((ymax-ymin)*xy_max-(ymax-ymin))/2;

		pos_ins->xmax[XX]=xmax+((xmax-xmin)*xy_max-(xmax-xmin))/2;
		pos_ins->xmax[YY]=ymax+((ymax-ymin)*xy_max-(ymax-ymin))/2;
	} else {
		pos_ins->xmin[XX]=xmin;
		pos_ins->xmin[YY]=ymin;

		pos_ins->xmax[XX]=xmax;
		pos_ins->xmax[YY]=ymax;
	}

	/* 6.0 is estimated thickness of bilayer */
	if( (zmax-zmin) < 6.0 )
	{
		pos_ins->xmin[ZZ]=zmin+(zmax-zmin)/2.0-3.0;
		pos_ins->xmax[ZZ]=zmin+(zmax-zmin)/2.0+3.0;
	} else {
		pos_ins->xmin[ZZ]=zmin;
		pos_ins->xmax[ZZ]=zmax;
	}

	return c;
}

real est_prot_area(pos_ins_t *pos_ins,rvec *r,t_block *ins_at, mem_t *mem_p)
{
	real x,y,dx=0.15,dy=0.15,area=0.0;
	real add;
	int c,at;

	for(x=pos_ins->xmin[XX];x<pos_ins->xmax[XX];x+=dx)
	{
		for(y=pos_ins->xmin[YY];y<pos_ins->xmax[YY];y+=dy)
		{
			c=0;
			add=0.0;
			do
			{
				at=ins_at->index[c];
				if ( (r[at][XX]>=x) && (r[at][XX]<x+dx) &&
						(r[at][YY]>=y) && (r[at][YY]<y+dy) &&
						(r[at][ZZ]>mem_p->zmin+1.0) && (r[at][ZZ]<mem_p->zmax-1.0) )
					add=1.0;
				c++;
			} while ( (c<ins_at->nr) && (add<0.5) );
			area+=add;
		}
	}
	area=area*dx*dy;

	return area;
}

void init_lip(matrix box, gmx_mtop_t *mtop, lip_t *lip)
{
	int i;
	real mem_area;
	int mol1=0;

	mem_area = box[XX][XX]*box[YY][YY]-box[XX][YY]*box[YY][XX];
	for(i=0;i<mtop->nmolblock;i++)
	{
		if(mtop->molblock[i].type == lip->id)
		{
			lip->nr=mtop->molblock[i].nmol;
			lip->natoms=mtop->molblock[i].natoms_mol;
		}
	}
	lip->area=2.0*mem_area/(double)lip->nr;

	for (i=0;i<lip->id;i++)
		mol1+=mtop->molblock[i].nmol;
	lip->mol1=mol1;
}

int init_mem_at(mem_t *mem_p, gmx_mtop_t *mtop, rvec *r, matrix box, pos_ins_t *pos_ins)
{
	int i,j,at,mol,nmol,nmolbox,count;
	t_block *mem_a;
	real z,zmin,zmax,mem_area;
	gmx_bool bNew;
	atom_id *mol_id;
	int type=0,block=0;

	nmol=count=0;
	mem_a=&(mem_p->mem_at);
	snew(mol_id,mem_a->nr);
/*	snew(index,mem_a->nr); */
	zmin=pos_ins->xmax[ZZ];
	zmax=pos_ins->xmin[ZZ];
	for(i=0;i<mem_a->nr;i++)
	{
		at=mem_a->index[i];
		if(	(r[at][XX]>pos_ins->xmin[XX]) && (r[at][XX]<pos_ins->xmax[XX]) &&
			(r[at][YY]>pos_ins->xmin[YY]) && (r[at][YY]<pos_ins->xmax[YY]) &&
			(r[at][ZZ]>pos_ins->xmin[ZZ]) && (r[at][ZZ]<pos_ins->xmax[ZZ]) )
		{
			mol = get_mol_id(at,mtop->nmolblock,mtop->molblock,&type,&block);

			bNew=TRUE;
			for(j=0;j<nmol;j++)
				if(mol == mol_id[j])
					bNew=FALSE;

			if(bNew)
			{
				mol_id[nmol]=mol;
				nmol++;
			}

			z=r[at][ZZ];
			if(z<zmin)					zmin=z;
			if(z>zmax)					zmax=z;

/*			index[count]=at;*/
			count++;
		}
	}

	mem_p->nmol=nmol;
	srenew(mol_id,nmol);
	mem_p->mol_id=mol_id;
/*	srenew(index,count);*/
/*	mem_p->mem_at.nr=count;*/
/*	sfree(mem_p->mem_at.index);*/
/*	mem_p->mem_at.index=index;*/

	if((zmax-zmin)>(box[ZZ][ZZ]-0.5))
		gmx_fatal(FARGS,"Something is wrong with your membrane. Max and min z values are %f and %f.\n"
				"Maybe your membrane is not centered in the box, but located at the box edge in the z-direction,\n"
				"so that one membrane is distributed over two periodic box images. Another possibility is that\n"
				"your water layer is not thick enough.\n",zmax,zmin);
	mem_p->zmin=zmin;
	mem_p->zmax=zmax;
	mem_p->zmed=(zmax-zmin)/2+zmin;

	/*number of membrane molecules in protein box*/
	nmolbox = count/mtop->molblock[block].natoms_mol;
	/*mem_area = box[XX][XX]*box[YY][YY]-box[XX][YY]*box[YY][XX];
	mem_p->lip_area = 2.0*mem_area/(double)mem_p->nmol;*/
	mem_area = (pos_ins->xmax[XX]-pos_ins->xmin[XX])*(pos_ins->xmax[YY]-pos_ins->xmin[YY]);
	mem_p->lip_area = 2.0*mem_area/(double)nmolbox;

	return mem_p->mem_at.nr;
}

void init_resize(t_block *ins_at,rvec *r_ins,pos_ins_t *pos_ins,mem_t *mem_p,rvec *r, gmx_bool bALLOW_ASYMMETRY)
{
	int i,j,at,c,outsidesum,gctr=0;
    int idxsum=0;

    /*sanity check*/
    for (i=0;i<pos_ins->pieces;i++)
          idxsum+=pos_ins->nidx[i];
    if (idxsum!=ins_at->nr)
          gmx_fatal(FARGS,"Piecewise sum of inserted atoms not same as size of group selected to insert.");

    snew(pos_ins->geom_cent,pos_ins->pieces);
    for (i=0;i<pos_ins->pieces;i++)
    {
    	c=0;
    	outsidesum=0;
    	for(j=0;j<DIM;j++)
    		pos_ins->geom_cent[i][j]=0;

    	for(j=0;j<DIM;j++)
    		pos_ins->geom_cent[i][j]=0;
    	for (j=0;j<pos_ins->nidx[i];j++)
    	{
    		at=pos_ins->subindex[i][j];
    		copy_rvec(r[at],r_ins[gctr]);
    		if( (r_ins[gctr][ZZ]<mem_p->zmax) && (r_ins[gctr][ZZ]>mem_p->zmin) )
    		{
    			rvec_inc(pos_ins->geom_cent[i],r_ins[gctr]);
    			c++;
    		}
    		else
    			outsidesum++;
    		gctr++;
    	}
    	if (c>0)
    		svmul(1/(double)c,pos_ins->geom_cent[i],pos_ins->geom_cent[i]);
    	if (!bALLOW_ASYMMETRY)
    		pos_ins->geom_cent[i][ZZ]=mem_p->zmed;

    	fprintf(stderr,"Embedding piece %d with center of geometry: %f %f %f\n",i,pos_ins->geom_cent[i][XX],pos_ins->geom_cent[i][YY],pos_ins->geom_cent[i][ZZ]);
    }
    fprintf(stderr,"\n");
}

void resize(t_block *ins_at, rvec *r_ins, rvec *r, pos_ins_t *pos_ins,rvec fac)
{
	int i,j,k,at,c=0;
	for (k=0;k<pos_ins->pieces;k++)
		for(i=0;i<pos_ins->nidx[k];i++)
		{
			at=pos_ins->subindex[k][i];
			for(j=0;j<DIM;j++)
				r[at][j]=pos_ins->geom_cent[k][j]+fac[j]*(r_ins[c][j]-pos_ins->geom_cent[k][j]);
			c++;
		}
}

int gen_rm_list(rmm_t *rm_p,t_block *ins_at,t_block *rest_at,t_pbc *pbc, gmx_mtop_t *mtop,
		rvec *r, rvec *r_ins, mem_t *mem_p, pos_ins_t *pos_ins, real probe_rad, int low_up_rm, gmx_bool bALLOW_ASYMMETRY)
{
	int i,j,k,l,at,at2,mol_id;
        int type=0,block=0;
	int nrm,nupper,nlower;
	real r_min_rad,z_lip,min_norm;
	gmx_bool bRM;
	rvec dr,dr_tmp;
	real *dist;
	int *order;

	r_min_rad=probe_rad*probe_rad;
	snew(rm_p->mol,mtop->mols.nr);
	snew(rm_p->block,mtop->mols.nr);
	nrm=nupper=0;
	nlower=low_up_rm;
	for(i=0;i<ins_at->nr;i++)
	{
		at=ins_at->index[i];
		for(j=0;j<rest_at->nr;j++)
		{
			at2=rest_at->index[j];
			pbc_dx(pbc,r[at],r[at2],dr);

			if(norm2(dr)<r_min_rad)
			{
				mol_id = get_mol_id(at2,mtop->nmolblock,mtop->molblock,&type,&block);
				bRM=TRUE;
				for(l=0;l<nrm;l++)
					if(rm_p->mol[l]==mol_id)
						bRM=FALSE;
				if(bRM)
				{
					/*fprintf(stderr,"%d wordt toegevoegd\n",mol_id);*/
					rm_p->mol[nrm]=mol_id;
					rm_p->block[nrm]=block;
					nrm++;
					z_lip=0.0;
					for(l=0;l<mem_p->nmol;l++)
					{
						if(mol_id==mem_p->mol_id[l])
						{
							for(k=mtop->mols.index[mol_id];k<mtop->mols.index[mol_id+1];k++)
								z_lip+=r[k][ZZ];
							z_lip/=mtop->molblock[block].natoms_mol;
							if(z_lip<mem_p->zmed)
								nlower++;
							else
								nupper++;
						}
					}
				}
			}
		}
	}

	/*make sure equal number of lipids from upper and lower layer are removed */
	if( (nupper!=nlower) && (!bALLOW_ASYMMETRY) )
	{
		snew(dist,mem_p->nmol);
		snew(order,mem_p->nmol);
		for(i=0;i<mem_p->nmol;i++)
		{
			at = mtop->mols.index[mem_p->mol_id[i]];
			pbc_dx(pbc,r[at],pos_ins->geom_cent[0],dr);
			if (pos_ins->pieces>1)
			{
				/*minimum dr value*/
				min_norm=norm2(dr);
				for (k=1;k<pos_ins->pieces;k++)
				{
					pbc_dx(pbc,r[at],pos_ins->geom_cent[k],dr_tmp);
					if (norm2(dr_tmp) < min_norm)
					{
						min_norm=norm2(dr_tmp);
						copy_rvec(dr_tmp,dr);
					}
				}
			}
			dist[i]=dr[XX]*dr[XX]+dr[YY]*dr[YY];
			j=i-1;
			while (j>=0 && dist[i]<dist[order[j]])
			{
				order[j+1]=order[j];
				j--;
			}
			order[j+1]=i;
		}

		i=0;
		while(nupper!=nlower)
		{
			mol_id=mem_p->mol_id[order[i]];
			block=get_block(mol_id,mtop->nmolblock,mtop->molblock);

			bRM=TRUE;
			for(l=0;l<nrm;l++)
				if(rm_p->mol[l]==mol_id)
					bRM=FALSE;
			if(bRM)
			{
				z_lip=0;
				for(k=mtop->mols.index[mol_id];k<mtop->mols.index[mol_id+1];k++)
					z_lip+=r[k][ZZ];
				z_lip/=mtop->molblock[block].natoms_mol;
				if(nupper>nlower && z_lip<mem_p->zmed)
				{
					rm_p->mol[nrm]=mol_id;
					rm_p->block[nrm]=block;
					nrm++;
					nlower++;
				}
				else if (nupper<nlower && z_lip>mem_p->zmed)
				{
					rm_p->mol[nrm]=mol_id;
					rm_p->block[nrm]=block;
					nrm++;
					nupper++;
				}
			}
			i++;

			if(i>mem_p->nmol)
				gmx_fatal(FARGS,"Trying to remove more lipid molecules than there are in the membrane");
		}
		sfree(dist);
		sfree(order);
	}

	rm_p->nr=nrm;
	srenew(rm_p->mol,nrm);
	srenew(rm_p->block,nrm);

	return nupper+nlower;
}

void rm_group(t_inputrec *ir, gmx_groups_t *groups, gmx_mtop_t *mtop, rmm_t *rm_p, t_state *state, t_block *ins_at, pos_ins_t *pos_ins)
{
	int i,j,k,n,rm,mol_id,at,block;
	rvec *x_tmp,*v_tmp;
	atom_id *list,*new_mols;
	unsigned char  *new_egrp[egcNR];
	gmx_bool bRM;

	snew(list,state->natoms);
	n=0;
	for(i=0;i<rm_p->nr;i++)
	{
		mol_id=rm_p->mol[i];
		at=mtop->mols.index[mol_id];
		block =rm_p->block[i];
		mtop->molblock[block].nmol--;
		for(j=0;j<mtop->molblock[block].natoms_mol;j++)
		{
			list[n]=at+j;
			n++;
		}

		mtop->mols.index[mol_id]=-1;
	}

	mtop->mols.nr-=rm_p->nr;
	mtop->mols.nalloc_index-=rm_p->nr;
	snew(new_mols,mtop->mols.nr);
	for(i=0;i<mtop->mols.nr+rm_p->nr;i++)
	{
		j=0;
		if(mtop->mols.index[i]!=-1)
		{
			new_mols[j]=mtop->mols.index[i];
			j++;
		}
	}
	sfree(mtop->mols.index);
	mtop->mols.index=new_mols;


	mtop->natoms-=n;
	state->natoms-=n;
	state->nalloc=state->natoms;
	snew(x_tmp,state->nalloc);
	snew(v_tmp,state->nalloc);

	for(i=0;i<egcNR;i++)
	{
		if(groups->grpnr[i]!=NULL)
		{
			groups->ngrpnr[i]=state->natoms;
			snew(new_egrp[i],state->natoms);
		}
	}

	rm=0;
	for (i=0;i<state->natoms+n;i++)
	{
		bRM=FALSE;
		for(j=0;j<n;j++)
		{
			if(i==list[j])
			{
				bRM=TRUE;
				rm++;
			}
		}

		if(!bRM)
		{
			for(j=0;j<egcNR;j++)
			{
				if(groups->grpnr[j]!=NULL)
				{
					new_egrp[j][i-rm]=groups->grpnr[j][i];
				}
			}
			copy_rvec(state->x[i],x_tmp[i-rm]);
			copy_rvec(state->v[i],v_tmp[i-rm]);
			for(j=0;j<ins_at->nr;j++)
			{
				if (i==ins_at->index[j])
					ins_at->index[j]=i-rm;
			}
			for(j=0;j<pos_ins->pieces;j++)
			{
				for(k=0;k<pos_ins->nidx[j];k++)
				{
					if (i==pos_ins->subindex[j][k])
						pos_ins->subindex[j][k]=i-rm;
				}
			}
		}
	}
	sfree(state->x);
	state->x=x_tmp;
	sfree(state->v);
	state->v=v_tmp;

	for(i=0;i<egcNR;i++)
	{
		if(groups->grpnr[i]!=NULL)
		{
			sfree(groups->grpnr[i]);
			groups->grpnr[i]=new_egrp[i];
		}
	}
}

int rm_bonded(t_block *ins_at, gmx_mtop_t *mtop)
{
	int i,j,m;
	int type,natom,nmol,at,atom1=0,rm_at=0;
	gmx_bool *bRM,bINS;
	/*this routine lives dangerously by assuming that all molecules of a given type are in order in the structure*/
	/*this routine does not live as dangerously as it seems. There is namely a check in mdrunner_membed to make
         *sure that g_membed exits with a warning when there are molecules of the same type not in the 
	 *ins_at index group. MGWolf 050710 */


	snew(bRM,mtop->nmoltype);
	for (i=0;i<mtop->nmoltype;i++)
	{
		bRM[i]=TRUE;
	}

	for (i=0;i<mtop->nmolblock;i++) 
	{
	    /*loop over molecule blocks*/
		type        =mtop->molblock[i].type;
		natom	    =mtop->molblock[i].natoms_mol;
		nmol		=mtop->molblock[i].nmol;

		for(j=0;j<natom*nmol && bRM[type]==TRUE;j++) 
		{
		    /*loop over atoms in the block*/
			at=j+atom1; /*atom index = block index + offset*/
			bINS=FALSE;

			for (m=0;(m<ins_at->nr) && (bINS==FALSE);m++)
			{
			    /*loop over atoms in insertion index group to determine if we're inserting one*/
				if(at==ins_at->index[m])
				{
					bINS=TRUE;
				}
			}
			bRM[type]=bINS;
		}
		atom1+=natom*nmol; /*update offset*/
		if(bRM[type])
		{
			rm_at+=natom*nmol; /*increment bonded removal counter by # atoms in block*/
		}
	}

	for(i=0;i<mtop->nmoltype;i++)
	{
		if(bRM[i])
		{
			for(j=0;j<F_LJ;j++)
			{
				mtop->moltype[i].ilist[j].nr=0;
			}
			for(j=F_POSRES;j<=F_VSITEN;j++)
			{
				mtop->moltype[i].ilist[j].nr=0;
			}
		}
	}
	sfree(bRM);

	return rm_at;
}

void top_update(const char *topfile, char *ins, rmm_t *rm_p, gmx_mtop_t *mtop)
{
#define TEMP_FILENM "temp.top"
	int	bMolecules=0;
	FILE	*fpin,*fpout;
	char	buf[STRLEN],buf2[STRLEN],*temp;
	int		i,*nmol_rm,nmol,line;

	fpin  = ffopen(topfile,"r");
	fpout = ffopen(TEMP_FILENM,"w");

	snew(nmol_rm,mtop->nmoltype);
	for(i=0;i<rm_p->nr;i++)
		nmol_rm[rm_p->block[i]]++;

	line=0;
	while(fgets(buf,STRLEN,fpin))
	{
		line++;
		if(buf[0]!=';')
		{
			strcpy(buf2,buf);
			if ((temp=strchr(buf2,'\n')) != NULL)
				temp[0]='\0';
			ltrim(buf2);

			if (buf2[0]=='[')
			{
				buf2[0]=' ';
				if ((temp=strchr(buf2,'\n')) != NULL)
					temp[0]='\0';
				rtrim(buf2);
				if (buf2[strlen(buf2)-1]==']')
				{
					buf2[strlen(buf2)-1]='\0';
					ltrim(buf2);
					rtrim(buf2);
					if (gmx_strcasecmp(buf2,"molecules")==0)
						bMolecules=1;
				}
				fprintf(fpout,"%s",buf);
			} else if (bMolecules==1)
			{
				for(i=0;i<mtop->nmolblock;i++)
				{
					nmol=mtop->molblock[i].nmol;
					sprintf(buf,"%-15s %5d\n",*(mtop->moltype[mtop->molblock[i].type].name),nmol);
					fprintf(fpout,"%s",buf);
				}
				bMolecules=2;
			} else if (bMolecules==2)
			{
				/* print nothing */
			} else 
			{
				fprintf(fpout,"%s",buf);
			}
		} else 
		{
			fprintf(fpout,"%s",buf);
		}
	}

	fclose(fpout);
	/* use ffopen to generate backup of topinout */
	fpout=ffopen(topfile,"w");
	fclose(fpout);
	rename(TEMP_FILENM,topfile);
#undef TEMP_FILENM
}

double do_md_membed(FILE *fplog,t_commrec *cr,int nfile,const t_filenm fnm[],
             const output_env_t oenv, gmx_bool bVerbose,gmx_bool bCompact,
             int nstglobalcomm,
             gmx_vsite_t *vsite,gmx_constr_t constr,
             int stepout,t_inputrec *ir,
             gmx_mtop_t *top_global,
             t_fcdata *fcd,
             t_state *state_global,
             t_mdatoms *mdatoms,
             t_nrnb *nrnb,gmx_wallcycle_t wcycle,
             gmx_edsam_t ed,t_forcerec *fr,
             int repl_ex_nst,int repl_ex_seed,
             real cpt_period,real max_hours,
             const char *deviceOptions,
             unsigned long Flags,
             gmx_runtime_t *runtime,
             rvec fac, rvec *r_ins, pos_ins_t *pos_ins, t_block *ins_at,
             real xy_step, real z_step, int it_xy, int it_z)
{
    gmx_mdoutf_t *outf;
    gmx_large_int_t step,step_rel;
    double     run_time;
    double     t,t0,lam0;
    gmx_bool       bGStatEveryStep,bGStat,bNstEner,bCalcEnerPres;
    gmx_bool       bNS,bNStList,bSimAnn,bStopCM,bRerunMD,bNotLastFrame=FALSE,
               bFirstStep,bStateFromTPX,bInitStep,bLastStep,
               bBornRadii,bStartingFromCpt;
    gmx_bool       bDoDHDL=FALSE;
    gmx_bool       do_ene,do_log,do_verbose,bRerunWarnNoV=TRUE,
               bForceUpdate=FALSE,bCPT;
    int        mdof_flags;
    gmx_bool       bMasterState;
    int        force_flags,cglo_flags;
    tensor     force_vir,shake_vir,total_vir,tmp_vir,pres;
    int        i,m;
    t_trxstatus *status;
    rvec       mu_tot;
    t_vcm      *vcm;
    t_state    *bufstate=NULL;
    matrix     *scale_tot,pcoupl_mu,M,ebox;
    gmx_nlheur_t nlh;
    t_trxframe rerun_fr;
/*    gmx_repl_ex_t repl_ex=NULL;*/
    int        nchkpt=1;

    gmx_localtop_t *top;
    t_mdebin *mdebin=NULL;
    t_state    *state=NULL;
    rvec       *f_global=NULL;
    int        n_xtc=-1;
    rvec       *x_xtc=NULL;
    gmx_enerdata_t *enerd;
    rvec       *f=NULL;
    gmx_global_stat_t gstat;
    gmx_update_t upd=NULL;
    t_graph    *graph=NULL;
    globsig_t   gs;

    gmx_bool        bFFscan;
    gmx_groups_t *groups;
    gmx_ekindata_t *ekind, *ekind_save;
    gmx_shellfc_t shellfc;
    int         count,nconverged=0;
    real        timestep=0;
    double      tcount=0;
    gmx_bool        bIonize=FALSE;
    gmx_bool        bTCR=FALSE,bConverged=TRUE,bOK,bSumEkinhOld,bExchanged;
    gmx_bool        bAppend;
    gmx_bool        bResetCountersHalfMaxH=FALSE;
    gmx_bool        bVV,bIterations,bIterate,bFirstIterate,bTemp,bPres,bTrotter;
    real        temp0,dvdl;
    int         a0,a1,ii;
    rvec        *xcopy=NULL,*vcopy=NULL,*cbuf=NULL;
    matrix      boxcopy={{0}},lastbox;
	real        veta_save,pcurr,scalevir,tracevir;
	real        vetanew = 0;
    double      cycles;
	real        last_conserved = 0;
    real        last_ekin = 0;
	t_extmass   MassQ;
    int         **trotter_seq;
    char        sbuf[STEPSTRSIZE],sbuf2[STEPSTRSIZE];
    int         handled_stop_condition=gmx_stop_cond_none; /* compare to get_stop_condition*/
    gmx_iterate_t iterate;
#ifdef GMX_FAHCORE
    /* Temporary addition for FAHCORE checkpointing */
    int chkpt_ret;
#endif

    /* Check for special mdrun options */
    bRerunMD = (Flags & MD_RERUN);
    bIonize  = (Flags & MD_IONIZE);
    bFFscan  = (Flags & MD_FFSCAN);
    bAppend  = (Flags & MD_APPENDFILES);
    bGStatEveryStep = FALSE;
    if (Flags & MD_RESETCOUNTERSHALFWAY)
    {
        if (ir->nsteps > 0)
        {
            /* Signal to reset the counters half the simulation steps. */
            wcycle_set_reset_counters(wcycle,ir->nsteps/2);
        }
        /* Signal to reset the counters halfway the simulation time. */
        bResetCountersHalfMaxH = (max_hours > 0);
    }

    /* md-vv uses averaged full step velocities for T-control
       md-vv-avek uses averaged half step velocities for T-control (but full step ekin for P control)
       md uses averaged half step kinetic energies to determine temperature unless defined otherwise by GMX_EKIN_AVE_VEL; */
    bVV = EI_VV(ir->eI);
    if (bVV) /* to store the initial velocities while computing virial */
    {
        snew(cbuf,top_global->natoms);
    }
    /* all the iteratative cases - only if there are constraints */
    bIterations = ((IR_NPT_TROTTER(ir)) && (constr) && (!bRerunMD));
    bTrotter = (bVV && (IR_NPT_TROTTER(ir) || (IR_NVT_TROTTER(ir))));

    if (bRerunMD)
    {
        /* Since we don't know if the frames read are related in any way,
         * rebuild the neighborlist at every step.
         */
        ir->nstlist       = 1;
        ir->nstcalcenergy = 1;
        nstglobalcomm     = 1;
    }

    check_ir_old_tpx_versions(cr,fplog,ir,top_global);

    nstglobalcomm = check_nstglobalcomm(fplog,cr,nstglobalcomm,ir);
    /*bGStatEveryStep = (nstglobalcomm == 1);*/
    bGStatEveryStep = FALSE;

    if (!bGStatEveryStep && ir->nstlist == -1 && fplog != NULL)
    {
        fprintf(fplog,
                "To reduce the energy communication with nstlist = -1\n"
                "the neighbor list validity should not be checked at every step,\n"
                "this means that exact integration is not guaranteed.\n"
                "The neighbor list validity is checked after:\n"
                "  <n.list life time> - 2*std.dev.(n.list life time)  steps.\n"
                "In most cases this will result in exact integration.\n"
                "This reduces the energy communication by a factor of 2 to 3.\n"
                "If you want less energy communication, set nstlist > 3.\n\n");
    }

    if (bRerunMD || bFFscan)
    {
        ir->nstxtcout = 0;
    }
    groups = &top_global->groups;

    /* Initial values */
    init_md(fplog,cr,ir,oenv,&t,&t0,&state_global->lambda,&lam0,
            nrnb,top_global,&upd,
            nfile,fnm,&outf,&mdebin,
            force_vir,shake_vir,mu_tot,&bSimAnn,&vcm,state_global,Flags);

    clear_mat(total_vir);
    clear_mat(pres);
    /* Energy terms and groups */
    snew(enerd,1);
    init_enerdata(top_global->groups.grps[egcENER].nr,ir->n_flambda,enerd);
    if (DOMAINDECOMP(cr))
    {
        f = NULL;
    }
    else
    {
        snew(f,top_global->natoms);
    }

    /* Kinetic energy data */
    snew(ekind,1);
    init_ekindata(fplog,top_global,&(ir->opts),ekind);
    /* needed for iteration of constraints */
    snew(ekind_save,1);
    init_ekindata(fplog,top_global,&(ir->opts),ekind_save);
    /* Copy the cos acceleration to the groups struct */
    ekind->cosacc.cos_accel = ir->cos_accel;

    gstat = global_stat_init(ir);
    debug_gmx();

    /* Check for polarizable models and flexible constraints */
    shellfc = init_shell_flexcon(fplog,
                                 top_global,n_flexible_constraints(constr),
                                 (ir->bContinuation ||
                                  (DOMAINDECOMP(cr) && !MASTER(cr))) ?
                                 NULL : state_global->x);

/*    if (DEFORM(*ir))
    {
#ifdef GMX_THREADS
        tMPI_Thread_mutex_lock(&deform_init_box_mutex);
#endif
        set_deform_reference_box(upd,
                                 deform_init_init_step_tpx,
                                 deform_init_box_tpx);
#ifdef GMX_THREADS
        tMPI_Thread_mutex_unlock(&deform_init_box_mutex);
#endif
    }*/

/*    {
        double io = compute_io(ir,top_global->natoms,groups,mdebin->ebin->nener,1);
        if ((io > 2000) && MASTER(cr))
            fprintf(stderr,
                    "\nWARNING: This run will generate roughly %.0f Mb of data\n\n",
                    io);
    }*/

    if (DOMAINDECOMP(cr)) {
        top = dd_init_local_top(top_global);

        snew(state,1);
        dd_init_local_state(cr->dd,state_global,state);

        if (DDMASTER(cr->dd) && ir->nstfout) {
            snew(f_global,state_global->natoms);
        }
    } else {
        if (PAR(cr)) {
            /* Initialize the particle decomposition and split the topology */
            top = split_system(fplog,top_global,ir,cr);

            pd_cg_range(cr,&fr->cg0,&fr->hcg);
            pd_at_range(cr,&a0,&a1);
        } else {
            top = gmx_mtop_generate_local_top(top_global,ir);

            a0 = 0;
            a1 = top_global->natoms;
        }

        state = partdec_init_local_state(cr,state_global);
        f_global = f;

        atoms2md(top_global,ir,0,NULL,a0,a1-a0,mdatoms);

        if (vsite) {
            set_vsite_top(vsite,top,mdatoms,cr);
        }

        if (ir->ePBC != epbcNONE && !ir->bPeriodicMols) {
            graph = mk_graph(fplog,&(top->idef),0,top_global->natoms,FALSE,FALSE);
        }

        if (shellfc) {
            make_local_shells(cr,mdatoms,shellfc);
        }

        if (ir->pull && PAR(cr)) {
            dd_make_local_pull_groups(NULL,ir->pull,mdatoms);
        }
    }

    if (DOMAINDECOMP(cr))
    {
        /* Distribute the charge groups over the nodes from the master node */
        dd_partition_system(fplog,ir->init_step,cr,TRUE,1,
                            state_global,top_global,ir,
                            state,&f,mdatoms,top,fr,
                            vsite,shellfc,constr,
                            nrnb,wcycle,FALSE);
    }

    update_mdatoms(mdatoms,state->lambda);

    if (MASTER(cr))
    {
        if (opt2bSet("-cpi",nfile,fnm))
        {
            /* Update mdebin with energy history if appending to output files */
            if ( Flags & MD_APPENDFILES )
            {
                restore_energyhistory_from_state(mdebin,&state_global->enerhist);
            }
            else
            {
                /* We might have read an energy history from checkpoint,
                 * free the allocated memory and reset the counts.
                 */
                done_energyhistory(&state_global->enerhist);
                init_energyhistory(&state_global->enerhist);
            }
        }
        /* Set the initial energy history in state by updating once */
        update_energyhistory(&state_global->enerhist,mdebin);
    }

    if ((state->flags & (1<<estLD_RNG)) && (Flags & MD_READ_RNG)) {
        /* Set the random state if we read a checkpoint file */
        set_stochd_state(upd,state);
    }

    /* Initialize constraints */
    if (constr) {
        if (!DOMAINDECOMP(cr))
            set_constraints(constr,top,ir,mdatoms,cr);
    }

    /* Check whether we have to GCT stuff */
 /*   bTCR = ftp2bSet(efGCT,nfile,fnm);
    if (bTCR) {
        if (MASTER(cr)) {
            fprintf(stderr,"Will do General Coupling Theory!\n");
        }
        gnx = top_global->mols.nr;
        snew(grpindex,gnx);
        for(i=0; (i<gnx); i++) {
            grpindex[i] = i;
        }
    }*/

/*    if (repl_ex_nst > 0 && MASTER(cr))
        repl_ex = init_replica_exchange(fplog,cr->ms,state_global,ir,
                                        repl_ex_nst,repl_ex_seed);*/

    if (!ir->bContinuation && !bRerunMD)
    {
        if (mdatoms->cFREEZE && (state->flags & (1<<estV)))
        {
            /* Set the velocities of frozen particles to zero */
            for(i=mdatoms->start; i<mdatoms->start+mdatoms->homenr; i++)
            {
                for(m=0; m<DIM; m++)
                {
                    if (ir->opts.nFreeze[mdatoms->cFREEZE[i]][m])
                    {
                        state->v[i][m] = 0;
                    }
                }
            }
        }

        if (constr)
        {
            /* Constrain the initial coordinates and velocities */
            do_constrain_first(fplog,constr,ir,mdatoms,state,f,
                               graph,cr,nrnb,fr,top,shake_vir);
        }
        if (vsite)
        {
            /* Construct the virtual sites for the initial configuration */
            construct_vsites(fplog,vsite,state->x,nrnb,ir->delta_t,NULL,
                             top->idef.iparams,top->idef.il,
                             fr->ePBC,fr->bMolPBC,graph,cr,state->box);
        }
    }

    debug_gmx();

    /* I'm assuming we need global communication the first time! MRS */
    cglo_flags = (CGLO_TEMPERATURE | CGLO_GSTAT
                  | (bVV ? CGLO_PRESSURE:0)
                  | (bVV ? CGLO_CONSTRAINT:0)
                  | (bRerunMD ? CGLO_RERUNMD:0)
                  | ((Flags & MD_READ_EKIN) ? CGLO_READEKIN:0));

    bSumEkinhOld = FALSE;
    compute_globals(fplog,gstat,cr,ir,fr,ekind,state,state_global,mdatoms,nrnb,vcm,
                    wcycle,enerd,force_vir,shake_vir,total_vir,pres,mu_tot,
                    constr,NULL,FALSE,state->box,
                    top_global,&pcurr,top_global->natoms,&bSumEkinhOld,cglo_flags);
    if (ir->eI == eiVVAK) {
        /* a second call to get the half step temperature initialized as well */
        /* we do the same call as above, but turn the pressure off -- internally, this
           is recognized as a velocity verlet half-step kinetic energy calculation.
           This minimized excess variables, but perhaps loses some logic?*/

        compute_globals(fplog,gstat,cr,ir,fr,ekind,state,state_global,mdatoms,nrnb,vcm,
                        wcycle,enerd,force_vir,shake_vir,total_vir,pres,mu_tot,
                        constr,NULL,FALSE,state->box,
                        top_global,&pcurr,top_global->natoms,&bSumEkinhOld,
                        cglo_flags &~ CGLO_PRESSURE);
    }

    /* Calculate the initial half step temperature, and save the ekinh_old */
    if (!(Flags & MD_STARTFROMCPT))
    {
        for(i=0; (i<ir->opts.ngtc); i++)
        {
            copy_mat(ekind->tcstat[i].ekinh,ekind->tcstat[i].ekinh_old);
        }
    }
    if (ir->eI != eiVV) 
    {
        enerd->term[F_TEMP] *= 2; /* result of averages being done over previous and current step,
                                     and there is no previous step */
    }
    temp0 = enerd->term[F_TEMP];

    /* if using an iterative algorithm, we need to create a working directory for the state. */
    if (bIterations)
    {
            bufstate = init_bufstate(state);
    }
    if (bFFscan)
    {
        snew(xcopy,state->natoms);
        snew(vcopy,state->natoms);
        copy_rvecn(state->x,xcopy,0,state->natoms);
        copy_rvecn(state->v,vcopy,0,state->natoms);
        copy_mat(state->box,boxcopy);
    }

    /* need to make an initiation call to get the Trotter variables set, as well as other constants for non-trotter
       temperature control */
    trotter_seq = init_npt_vars(ir,state,&MassQ,bTrotter);

    if (MASTER(cr))
    {
        if (constr && !ir->bContinuation && ir->eConstrAlg == econtLINCS)
        {
            fprintf(fplog,
                    "RMS relative constraint deviation after constraining: %.2e\n",
                    constr_rmsd(constr,FALSE));
        }
        fprintf(fplog,"Initial temperature: %g K\n",enerd->term[F_TEMP]);
        if (bRerunMD)
        {
            fprintf(stderr,"starting md rerun '%s', reading coordinates from"
                    " input trajectory '%s'\n\n",
                    *(top_global->name),opt2fn("-rerun",nfile,fnm));
            if (bVerbose)
            {
                fprintf(stderr,"Calculated time to finish depends on nsteps from "
                        "run input file,\nwhich may not correspond to the time "
                        "needed to process input trajectory.\n\n");
            }
        }
        else
        {
            char tbuf[20];
            fprintf(stderr,"starting mdrun '%s'\n",
                    *(top_global->name));
            if (ir->nsteps >= 0)
            {
                sprintf(tbuf,"%8.1f",(ir->init_step+ir->nsteps)*ir->delta_t);
            }
            else
            {
                sprintf(tbuf,"%s","infinite");
            }
            if (ir->init_step > 0)
            {
                fprintf(stderr,"%s steps, %s ps (continuing from step %s, %8.1f ps).\n",
                        gmx_step_str(ir->init_step+ir->nsteps,sbuf),tbuf,
                        gmx_step_str(ir->init_step,sbuf2),
                        ir->init_step*ir->delta_t);
            }
            else
            {
                fprintf(stderr,"%s steps, %s ps.\n",
                        gmx_step_str(ir->nsteps,sbuf),tbuf);
            }
        }
        fprintf(fplog,"\n");
    }

    /* Set and write start time */
    runtime_start(runtime);
    print_date_and_time(fplog,cr->nodeid,"Started mdrun",runtime);
    wallcycle_start(wcycle,ewcRUN);
    if (fplog)
        fprintf(fplog,"\n");

    /* safest point to do file checkpointing is here.  More general point would be immediately before integrator call */
/*#ifdef GMX_FAHCORE
    chkpt_ret=fcCheckPointParallel( cr->nodeid,
                                    NULL,0);
    if ( chkpt_ret == 0 )
        gmx_fatal( 3,__FILE__,__LINE__, "Checkpoint error on step %d\n", 0 );
#endif*/

    debug_gmx();
    /***********************************************************
     *
     *             Loop over MD steps
     *
     ************************************************************/

    /* if rerunMD then read coordinates and velocities from input trajectory */
    if (bRerunMD)
    {
        if (getenv("GMX_FORCE_UPDATE"))
        {
            bForceUpdate = TRUE;
        }

        bNotLastFrame = read_first_frame(oenv,&status,
                                         opt2fn("-rerun",nfile,fnm),
                                         &rerun_fr,TRX_NEED_X | TRX_READ_V);
        if (rerun_fr.natoms != top_global->natoms)
        {
            gmx_fatal(FARGS,
                      "Number of atoms in trajectory (%d) does not match the "
                      "run input file (%d)\n",
                      rerun_fr.natoms,top_global->natoms);
        }
        if (ir->ePBC != epbcNONE)
        {
            if (!rerun_fr.bBox)
            {
                gmx_fatal(FARGS,"Rerun trajectory frame step %d time %f does not contain a box, while pbc is used",rerun_fr.step,rerun_fr.time);
            }
            if (max_cutoff2(ir->ePBC,rerun_fr.box) < sqr(fr->rlistlong))
            {
                gmx_fatal(FARGS,"Rerun trajectory frame step %d time %f has too small box dimensions",rerun_fr.step,rerun_fr.time);
            }

            /* Set the shift vectors.
             * Necessary here when have a static box different from the tpr box.
             */
            calc_shifts(rerun_fr.box,fr->shift_vec);
        }
    }

    /* loop over MD steps or if rerunMD to end of input trajectory */
    bFirstStep = TRUE;
    /* Skip the first Nose-Hoover integration when we get the state from tpx */
    bStateFromTPX = !opt2bSet("-cpi",nfile,fnm);
    bInitStep = bFirstStep && (bStateFromTPX || bVV);
    bStartingFromCpt = (Flags & MD_STARTFROMCPT) && bInitStep;
    bLastStep    = FALSE;
    bSumEkinhOld = FALSE;
    bExchanged   = FALSE;

    init_global_signals(&gs,cr,ir,repl_ex_nst);

    step = ir->init_step;
    step_rel = 0;

    if (ir->nstlist == -1)
    {
        init_nlistheuristics(&nlh,bGStatEveryStep,step);
    }

    bLastStep = (bRerunMD || (ir->nsteps >= 0 && step_rel > ir->nsteps));
    while (!bLastStep || (bRerunMD && bNotLastFrame)) {

        wallcycle_start(wcycle,ewcSTEP);

        GMX_MPE_LOG(ev_timestep1);

        if (bRerunMD) {
            if (rerun_fr.bStep) {
                step = rerun_fr.step;
                step_rel = step - ir->init_step;
            }
            if (rerun_fr.bTime) {
                t = rerun_fr.time;
            }
            else
            {
                t = step;
            }
        }
        else
        {
            bLastStep = (step_rel == ir->nsteps);
            t = t0 + step*ir->delta_t;
        }

        if (ir->efep != efepNO)
        {
            if (bRerunMD && rerun_fr.bLambda && (ir->delta_lambda!=0))
            {
                state_global->lambda = rerun_fr.lambda;
            }
            else
            {
                state_global->lambda = lam0 + step*ir->delta_lambda;
            }
            state->lambda = state_global->lambda;
            bDoDHDL = do_per_step(step,ir->nstdhdl);
        }

        if (bSimAnn)
        {
            update_annealing_target_temp(&(ir->opts),t);
        }

        if (bRerunMD)
        {
            if (!(DOMAINDECOMP(cr) && !MASTER(cr)))
            {
                for(i=0; i<state_global->natoms; i++)
                {
                    copy_rvec(rerun_fr.x[i],state_global->x[i]);
                }
                if (rerun_fr.bV)
                {
                    for(i=0; i<state_global->natoms; i++)
                    {
                        copy_rvec(rerun_fr.v[i],state_global->v[i]);
                    }
                }
                else
                {
                    for(i=0; i<state_global->natoms; i++)
                    {
                        clear_rvec(state_global->v[i]);
                    }
                    if (bRerunWarnNoV)
                    {
                        fprintf(stderr,"\nWARNING: Some frames do not contain velocities.\n"
                                "         Ekin, temperature and pressure are incorrect,\n"
                                "         the virial will be incorrect when constraints are present.\n"
                                "\n");
                        bRerunWarnNoV = FALSE;
                    }
                }
            }
            copy_mat(rerun_fr.box,state_global->box);
            copy_mat(state_global->box,state->box);

            if (vsite && (Flags & MD_RERUN_VSITE))
            {
                if (DOMAINDECOMP(cr))
                {
                    gmx_fatal(FARGS,"Vsite recalculation with -rerun is not implemented for domain decomposition, use particle decomposition");
                }
                if (graph)
                {
                    /* Following is necessary because the graph may get out of sync
                     * with the coordinates if we only have every N'th coordinate set
                     */
                    mk_mshift(fplog,graph,fr->ePBC,state->box,state->x);
                    shift_self(graph,state->box,state->x);
                }
                construct_vsites(fplog,vsite,state->x,nrnb,ir->delta_t,state->v,
                                 top->idef.iparams,top->idef.il,
                                 fr->ePBC,fr->bMolPBC,graph,cr,state->box);
                if (graph)
                {
                    unshift_self(graph,state->box,state->x);
                }
            }
        }

        /* Stop Center of Mass motion */
        bStopCM = (ir->comm_mode != ecmNO && do_per_step(step,ir->nstcomm));

        /* Copy back starting coordinates in case we're doing a forcefield scan */
        if (bFFscan)
        {
            for(ii=0; (ii<state->natoms); ii++)
            {
                copy_rvec(xcopy[ii],state->x[ii]);
                copy_rvec(vcopy[ii],state->v[ii]);
            }
            copy_mat(boxcopy,state->box);
        }

        if (bRerunMD)
        {
            /* for rerun MD always do Neighbour Searching */
            bNS = (bFirstStep || ir->nstlist != 0);
            bNStList = bNS;
        }
        else
        {
            /* Determine whether or not to do Neighbour Searching and LR */
            bNStList = (ir->nstlist > 0  && step % ir->nstlist == 0);

            bNS = (bFirstStep || bExchanged || bNStList ||
                   (ir->nstlist == -1 && nlh.nabnsb > 0));

            if (bNS && ir->nstlist == -1)
            {
                set_nlistheuristics(&nlh,bFirstStep || bExchanged,step);
            }
        }

        /* < 0 means stop at next step, > 0 means stop at next NS step */
        if ( (gs.set[eglsSTOPCOND] < 0 ) ||
             ( (gs.set[eglsSTOPCOND] > 0 ) && ( bNS || ir->nstlist==0)) )
        {
            bLastStep = TRUE;
        }

        /* Determine whether or not to update the Born radii if doing GB */
        bBornRadii=bFirstStep;
        if (ir->implicit_solvent && (step % ir->nstgbradii==0))
        {
            bBornRadii=TRUE;
        }

        do_log = do_per_step(step,ir->nstlog) || bFirstStep || bLastStep;
        do_verbose = bVerbose &&
                  (step % stepout == 0 || bFirstStep || bLastStep);

        if (bNS && !(bFirstStep && ir->bContinuation && !bRerunMD))
        {
            if (bRerunMD)
            {
                bMasterState = TRUE;
            }
            else
            {
                bMasterState = FALSE;
                /* Correct the new box if it is too skewed */
                if (DYNAMIC_BOX(*ir))
                {
                    if (correct_box(fplog,step,state->box,graph))
                    {
                        bMasterState = TRUE;
                    }
                }
                if (DOMAINDECOMP(cr) && bMasterState)
                {
                    dd_collect_state(cr->dd,state,state_global);
                }
            }

            if (DOMAINDECOMP(cr))
            {
                /* Repartition the domain decomposition */
                wallcycle_start(wcycle,ewcDOMDEC);
                dd_partition_system(fplog,step,cr,
                                    bMasterState,nstglobalcomm,
                                    state_global,top_global,ir,
                                    state,&f,mdatoms,top,fr,
                                    vsite,shellfc,constr,
                                    nrnb,wcycle,do_verbose);
                wallcycle_stop(wcycle,ewcDOMDEC);
                /* If using an iterative integrator, reallocate space to match the decomposition */
            }
        }

        if (MASTER(cr) && do_log && !bFFscan)
        {
            print_ebin_header(fplog,step,t,state->lambda);
        }

        if (ir->efep != efepNO)
        {
            update_mdatoms(mdatoms,state->lambda);
        }

        if (bRerunMD && rerun_fr.bV)
        {

            /* We need the kinetic energy at minus the half step for determining
             * the full step kinetic energy and possibly for T-coupling.*/
            /* This may not be quite working correctly yet . . . . */
            compute_globals(fplog,gstat,cr,ir,fr,ekind,state,state_global,mdatoms,nrnb,vcm,
                            wcycle,enerd,NULL,NULL,NULL,NULL,mu_tot,
                            constr,NULL,FALSE,state->box,
                            top_global,&pcurr,top_global->natoms,&bSumEkinhOld,
                            CGLO_RERUNMD | CGLO_GSTAT | CGLO_TEMPERATURE);
        }
        clear_mat(force_vir);

        /* Ionize the atoms if necessary */
/*        if (bIonize)
        {
            ionize(fplog,oenv,mdatoms,top_global,t,ir,state->x,state->v,
                   mdatoms->start,mdatoms->start+mdatoms->homenr,state->box,cr);
        }*/

        /* Update force field in ffscan program */
/*        if (bFFscan)
        {
            if (update_forcefield(fplog,
                                  nfile,fnm,fr,
                                  mdatoms->nr,state->x,state->box)) {
                if (gmx_parallel_env_initialized())
                {
                    gmx_finalize();
                }
                exit(0);
            }
        }*/

        GMX_MPE_LOG(ev_timestep2);

        /* We write a checkpoint at this MD step when:
         * either at an NS step when we signalled through gs,
         * or at the last step (but not when we do not want confout),
         * but never at the first step or with rerun.
         */
/*        bCPT = (((gs.set[eglsCHKPT] && bNS) ||
                 (bLastStep && (Flags & MD_CONFOUT))) &&
                step > ir->init_step && !bRerunMD);
        if (bCPT)
        {
            gs.set[eglsCHKPT] = 0;
        }*/

        /* Determine the energy and pressure:
         * at nstcalcenergy steps and at energy output steps (set below).
         */
        bNstEner = (bGStatEveryStep || do_per_step(step,ir->nstcalcenergy));
        bCalcEnerPres = bNstEner;

        /* Do we need global communication ? */
        bGStat = (bCalcEnerPres || bStopCM ||
                  (ir->nstlist == -1 && !bRerunMD && step >= nlh.step_nscheck));

        do_ene = (do_per_step(step,ir->nstenergy) || bLastStep);

        if (do_ene || do_log)
        {
            bCalcEnerPres = TRUE;
            bGStat    = TRUE;
        }

        /* these CGLO_ options remain the same throughout the iteration */
        cglo_flags = ((bRerunMD ? CGLO_RERUNMD : 0) |
                      (bStopCM ? CGLO_STOPCM : 0) |
                      (bGStat ? CGLO_GSTAT : 0)
            );

        force_flags = (GMX_FORCE_STATECHANGED |
                       ((DYNAMIC_BOX(*ir) || bRerunMD) ? GMX_FORCE_DYNAMICBOX : 0) |
                       GMX_FORCE_ALLFORCES |
                       (bNStList ? GMX_FORCE_DOLR : 0) |
                       GMX_FORCE_SEPLRF |
                       (bCalcEnerPres ? GMX_FORCE_VIRIAL : 0) |
                       (bDoDHDL ? GMX_FORCE_DHDL : 0)
            );

        if (shellfc)
        {
            /* Now is the time to relax the shells */
            count=relax_shell_flexcon(fplog,cr,bVerbose,bFFscan ? step+1 : step,
                                      ir,bNS,force_flags,
                                      bStopCM,top,top_global,
                                      constr,enerd,fcd,
                                      state,f,force_vir,mdatoms,
                                      nrnb,wcycle,graph,groups,
                                      shellfc,fr,bBornRadii,t,mu_tot,
                                      state->natoms,&bConverged,vsite,
                                      outf->fp_field);
            tcount+=count;

            if (bConverged)
            {
                nconverged++;
            }
        }
        else
        {
            /* The coordinates (x) are shifted (to get whole molecules)
             * in do_force.
             * This is parallellized as well, and does communication too.
             * Check comments in sim_util.c
             */

            do_force(fplog,cr,ir,step,nrnb,wcycle,top,top_global,groups,
                     state->box,state->x,&state->hist,
                     f,force_vir,mdatoms,enerd,fcd,
                     state->lambda,graph,
                     fr,vsite,mu_tot,t,outf->fp_field,ed,bBornRadii,
                     (bNS ? GMX_FORCE_NS : 0) | force_flags);
        }

        GMX_BARRIER(cr->mpi_comm_mygroup);

 /*       if (bTCR)
        {
            mu_aver = calc_mu_aver(cr,state->x,mdatoms->chargeA,
                                   mu_tot,&top_global->mols,mdatoms,gnx,grpindex);
        }

        if (bTCR && bFirstStep)
        {
            tcr=init_coupling(fplog,nfile,fnm,cr,fr,mdatoms,&(top->idef));
            fprintf(fplog,"Done init_coupling\n");
            fflush(fplog);
        }*/

        /*  ############### START FIRST UPDATE HALF-STEP ############### */

        if (bVV && !bStartingFromCpt && !bRerunMD)
        {
            if (ir->eI == eiVV)
            {
                if (bInitStep)
                {
                    /* if using velocity verlet with full time step Ekin,
                     * take the first half step only to compute the
                     * virial for the first step. From there,
                     * revert back to the initial coordinates
                     * so that the input is actually the initial step.
                     */
                    copy_rvecn(state->v,cbuf,0,state->natoms); /* should make this better for parallelizing? */
                }

                /* this is for NHC in the Ekin(t+dt/2) version of vv */
                if (!bInitStep)
                {
		  trotter_update(ir,step,ekind,enerd,state,total_vir,mdatoms,&MassQ,trotter_seq,ettTSEQ2);
                }

		if (ir->eI == eiVVAK)
		{
		  update_tcouple(fplog,step,ir,state,ekind,wcycle,upd,&MassQ,mdatoms);
		}

                update_coords(fplog,step,ir,mdatoms,state,
                              f,fr->bTwinRange && bNStList,fr->f_twin,fcd,
                              ekind,M,wcycle,upd,bInitStep,etrtVELOCITY1,
                              cr,nrnb,constr,&top->idef);

                if (bIterations)
                {
                    gmx_iterate_init(&iterate,bIterations && !bInitStep);
                }
                /* for iterations, we save these vectors, as we will be self-consistently iterating
                   the calculations */
                /*#### UPDATE EXTENDED VARIABLES IN TROTTER FORMULATION */

                /* save the state */
                if (bIterations && iterate.bIterate) {
                    copy_coupling_state(state,bufstate,ekind,ekind_save,&(ir->opts));
                }
            }

            bFirstIterate = TRUE;
            while (bFirstIterate || (bIterations && iterate.bIterate))
            {
                if (bIterations && iterate.bIterate)
                {
                    copy_coupling_state(bufstate,state,ekind_save,ekind,&(ir->opts));
                    if (bFirstIterate && bTrotter)
                    {
                        /* The first time through, we need a decent first estimate
                           of veta(t+dt) to compute the constraints.  Do
                           this by computing the box volume part of the
                           trotter integration at this time. Nothing else
                           should be changed by this routine here.  If
                           !(first time), we start with the previous value
                           of veta.  */

                        veta_save = state->veta;
                        trotter_update(ir,step,ekind,enerd,state,total_vir,mdatoms,&MassQ,trotter_seq,ettTSEQ0);
                        vetanew = state->veta;
                        state->veta = veta_save;
                    }
                }

                bOK = TRUE;
                if ( !bRerunMD || rerun_fr.bV || bForceUpdate) {  /* Why is rerun_fr.bV here?  Unclear. */
                    dvdl = 0;

                    update_constraints(fplog,step,&dvdl,ir,ekind,mdatoms,state,graph,f,
                                       &top->idef,shake_vir,NULL,
                                       cr,nrnb,wcycle,upd,constr,
                                       bInitStep,TRUE,bCalcEnerPres,vetanew);

                    if (!bOK && !bFFscan)
                    {
                        gmx_fatal(FARGS,"Constraint error: Shake, Lincs or Settle could not solve the constrains");
                    }

                }
                else if (graph)
                { /* Need to unshift here if a do_force has been
                     called in the previous step */
                    unshift_self(graph,state->box,state->x);
                }


                if (bVV) {
                    /* if VV, compute the pressure and constraints */
                    /* if VV2, the pressure and constraints only if using pressure control.*/
                    bPres = (ir->eI==eiVV || IR_NPT_TROTTER(ir));
                    bTemp = ((ir->eI==eiVV &&(!bInitStep)) || (ir->eI==eiVVAK && IR_NPT_TROTTER(ir)));
                    compute_globals(fplog,gstat,cr,ir,fr,ekind,state,state_global,mdatoms,nrnb,vcm,
                                    wcycle,enerd,force_vir,shake_vir,total_vir,pres,mu_tot,
                                    constr,NULL,FALSE,state->box,
                                    top_global,&pcurr,top_global->natoms,&bSumEkinhOld,
                                    cglo_flags
                                    | CGLO_ENERGY
                                    | (bTemp ? CGLO_TEMPERATURE:0)
                                    | (bPres ? CGLO_PRESSURE : 0)
                                    | (bPres ? CGLO_CONSTRAINT : 0)
                                    | (iterate.bIterate ? CGLO_ITERATE : 0)
                                    | (bFirstIterate ? CGLO_FIRSTITERATE : 0)
                                    | CGLO_SCALEEKIN
                        );
                }
                /* explanation of above:
                   a) We compute Ekin at the full time step
                   if 1) we are using the AveVel Ekin, and it's not the
                   initial step, or 2) if we are using AveEkin, but need the full
                   time step kinetic energy for the pressure.
                   b) If we are using EkinAveEkin for the kinetic energy for the temperture control, we still feed in
                   EkinAveVel because it's needed for the pressure */

                /* temperature scaling and pressure scaling to produce the extended variables at t+dt */
                if (bVV && !bInitStep)
                {
		  trotter_update(ir,step,ekind,enerd,state,total_vir,mdatoms,&MassQ, trotter_seq,ettTSEQ2);
                }

                if (bIterations &&
                    done_iterating(cr,fplog,step,&iterate,bFirstIterate,
                                   state->veta,&vetanew))
                {
                    break;
                }
                bFirstIterate = FALSE;
            }

            if (bTrotter && !bInitStep) {
                copy_mat(shake_vir,state->svir_prev);
                copy_mat(force_vir,state->fvir_prev);
                if (IR_NVT_TROTTER(ir) && ir->eI==eiVV) {
                    /* update temperature and kinetic energy now that step is over - this is the v(t+dt) point */
                    enerd->term[F_TEMP] = sum_ekin(&(ir->opts),ekind,NULL,(ir->eI==eiVV),FALSE,FALSE);
                    enerd->term[F_EKIN] = trace(ekind->ekin);
                }
            }
            /* if it's the initial step, we performed this first step just to get the constraint virial */
            if (bInitStep && ir->eI==eiVV) {
                copy_rvecn(cbuf,state->v,0,state->natoms);
            }

            if (fr->bSepDVDL && fplog && do_log)
            {
                fprintf(fplog,sepdvdlformat,"Constraint",0.0,dvdl);
            }
            enerd->term[F_DHDL_CON] += dvdl;

            GMX_MPE_LOG(ev_timestep1);

        }

        /* MRS -- now done iterating -- compute the conserved quantity */
        if (bVV) {
            last_conserved = 0;
            if (IR_NVT_TROTTER(ir) || IR_NPT_TROTTER(ir))
            {
                last_conserved =
                    NPT_energy(ir,state,&MassQ);
                if ((ir->eDispCorr != edispcEnerPres) && (ir->eDispCorr != edispcAllEnerPres))
                {
                    last_conserved -= enerd->term[F_DISPCORR];
                }
            }
            if (ir->eI==eiVV) {
                last_ekin = enerd->term[F_EKIN]; /* does this get preserved through checkpointing? */
            }
        }

        /* ########  END FIRST UPDATE STEP  ############## */
        /* ########  If doing VV, we now have v(dt) ###### */

        /* ################## START TRAJECTORY OUTPUT ################# */

        /* Now we have the energies and forces corresponding to the
         * coordinates at time t. We must output all of this before
         * the update.
         * for RerunMD t is read from input trajectory
         */
        GMX_MPE_LOG(ev_output_start);

        mdof_flags = 0;
        if (do_per_step(step,ir->nstxout)) { mdof_flags |= MDOF_X; }
        if (do_per_step(step,ir->nstvout)) { mdof_flags |= MDOF_V; }
        if (do_per_step(step,ir->nstfout)) { mdof_flags |= MDOF_F; }
        if (do_per_step(step,ir->nstxtcout)) { mdof_flags |= MDOF_XTC; }
/*        if (bCPT) { mdof_flags |= MDOF_CPT; };*/

#ifdef GMX_FAHCORE
        if (MASTER(cr))
            fcReportProgress( ir->nsteps, step );

        if (bLastStep)
        {
            /* Enforce writing positions and velocities at end of run */
            mdof_flags |= (MDOF_X | MDOF_V);
        }
            /* sync bCPT and fc record-keeping */
/*            if (bCPT && MASTER(cr))
                fcRequestCheckPoint();*/
#endif

        if (mdof_flags != 0)
        {
            wallcycle_start(wcycle,ewcTRAJ);
/*            if (bCPT)
            {
                if (state->flags & (1<<estLD_RNG))
                {
                    get_stochd_state(upd,state);
                }
                if (MASTER(cr))
                {
                    if (bSumEkinhOld)
                    {
                        state_global->ekinstate.bUpToDate = FALSE;
                    }
                    else
                    {
                        update_ekinstate(&state_global->ekinstate,ekind);
                        state_global->ekinstate.bUpToDate = TRUE;
                    }
                    update_energyhistory(&state_global->enerhist,mdebin);
                }
            }*/
            write_traj(fplog,cr,outf,mdof_flags,top_global,
                       step,t,state,state_global,f,f_global,&n_xtc,&x_xtc);
/*            if (bCPT)
            {
                nchkpt++;
                bCPT = FALSE;
            }*/
            debug_gmx();
            if (bLastStep && step_rel == ir->nsteps &&
                (Flags & MD_CONFOUT) && MASTER(cr) &&
                !bRerunMD && !bFFscan)
            {
                /* x and v have been collected in write_traj,
                 * because a checkpoint file will always be written
                 * at the last step.
                 */
                fprintf(stderr,"\nWriting final coordinates.\n");
                if (ir->ePBC != epbcNONE && !ir->bPeriodicMols &&
                    DOMAINDECOMP(cr))
                {
                    /* Make molecules whole only for confout writing */
                    do_pbc_mtop(fplog,ir->ePBC,state->box,top_global,state_global->x);
                }
/*                write_sto_conf_mtop(ftp2fn(efSTO,nfile,fnm),
                                    *top_global->name,top_global,
                                    state_global->x,state_global->v,
                                    ir->ePBC,state->box);*/
                debug_gmx();
            }
            wallcycle_stop(wcycle,ewcTRAJ);
        }
        GMX_MPE_LOG(ev_output_finish);

        /* kludge -- virial is lost with restart for NPT control. Must restart */
        if (bStartingFromCpt && bVV)
        {
            copy_mat(state->svir_prev,shake_vir);
            copy_mat(state->fvir_prev,force_vir);
        }
        /*  ################## END TRAJECTORY OUTPUT ################ */

        /* Determine the pressure:
         * always when we want exact averages in the energy file,
         * at ns steps when we have pressure coupling,
         * otherwise only at energy output steps (set below).
         */

        bNstEner = (bGStatEveryStep || do_per_step(step,ir->nstcalcenergy));
        bCalcEnerPres = bNstEner;

        /* Do we need global communication ? */
        bGStat = (bGStatEveryStep || bStopCM || bNS ||
                  (ir->nstlist == -1 && !bRerunMD && step >= nlh.step_nscheck));

        do_ene = (do_per_step(step,ir->nstenergy) || bLastStep);

        if (do_ene || do_log)
        {
            bCalcEnerPres = TRUE;
            bGStat        = TRUE;
        }

        /* Determine the wallclock run time up till now */
        run_time = gmx_gettime() - (double)runtime->real;

        /* Check whether everything is still allright */
        if (((int)gmx_get_stop_condition() > handled_stop_condition)
#ifdef GMX_THREADS
	    && MASTER(cr)
#endif
	    )
        {
            /* this is just make gs.sig compatible with the hack
               of sending signals around by MPI_Reduce with together with
               other floats */
            if ( gmx_get_stop_condition() == gmx_stop_cond_next_ns )
                gs.sig[eglsSTOPCOND]=1;
            if ( gmx_get_stop_condition() == gmx_stop_cond_next )
                gs.sig[eglsSTOPCOND]=-1;
            /* < 0 means stop at next step, > 0 means stop at next NS step */
            if (fplog)
            {
                fprintf(fplog,
                        "\n\nReceived the %s signal, stopping at the next %sstep\n\n",
                        gmx_get_signal_name(),
                        gs.sig[eglsSTOPCOND]==1 ? "NS " : "");
                fflush(fplog);
            }
            fprintf(stderr,
                    "\n\nReceived the %s signal, stopping at the next %sstep\n\n",
                    gmx_get_signal_name(),
                    gs.sig[eglsSTOPCOND]==1 ? "NS " : "");
            fflush(stderr);
            handled_stop_condition=(int)gmx_get_stop_condition();
        }
        else if (MASTER(cr) && (bNS || ir->nstlist <= 0) &&
                 (max_hours > 0 && run_time > max_hours*60.0*60.0*0.99) &&
                 gs.sig[eglsSTOPCOND] == 0 && gs.set[eglsSTOPCOND] == 0)
        {
            /* Signal to terminate the run */
            gs.sig[eglsSTOPCOND] = 1;
            if (fplog)
            {
                fprintf(fplog,"\nStep %s: Run time exceeded %.3f hours, will terminate the run\n",gmx_step_str(step,sbuf),max_hours*0.99);
            }
            fprintf(stderr, "\nStep %s: Run time exceeded %.3f hours, will terminate the run\n",gmx_step_str(step,sbuf),max_hours*0.99);
        }

        if (bResetCountersHalfMaxH && MASTER(cr) &&
            run_time > max_hours*60.0*60.0*0.495)
        {
            gs.sig[eglsRESETCOUNTERS] = 1;
        }

        if (ir->nstlist == -1 && !bRerunMD)
        {
            /* When bGStatEveryStep=FALSE, global_stat is only called
             * when we check the atom displacements, not at NS steps.
             * This means that also the bonded interaction count check is not
             * performed immediately after NS. Therefore a few MD steps could
             * be performed with missing interactions.
             * But wrong energies are never written to file,
             * since energies are only written after global_stat
             * has been called.
             */
            if (step >= nlh.step_nscheck)
            {
                nlh.nabnsb = natoms_beyond_ns_buffer(ir,fr,&top->cgs,
                                                     nlh.scale_tot,state->x);
            }
            else
            {
                /* This is not necessarily true,
                 * but step_nscheck is determined quite conservatively.
                 */
                nlh.nabnsb = 0;
            }
        }

        /* In parallel we only have to check for checkpointing in steps
         * where we do global communication,
         *  otherwise the other nodes don't know.
         */
        if (MASTER(cr) && ((bGStat || !PAR(cr)) &&
                           cpt_period >= 0 &&
                           (cpt_period == 0 ||
                            run_time >= nchkpt*cpt_period*60.0)) &&
            gs.set[eglsCHKPT] == 0)
        {
            gs.sig[eglsCHKPT] = 1;
        }

        if (bIterations)
        {
            gmx_iterate_init(&iterate,bIterations);
        }

        /* for iterations, we save these vectors, as we will be redoing the calculations */
        if (bIterations && iterate.bIterate)
        {
            copy_coupling_state(state,bufstate,ekind,ekind_save,&(ir->opts));
        }
        bFirstIterate = TRUE;
        while (bFirstIterate || (bIterations && iterate.bIterate))
        {
            /* We now restore these vectors to redo the calculation with improved extended variables */
            if (bIterations)
            {
                copy_coupling_state(bufstate,state,ekind_save,ekind,&(ir->opts));
            }

            /* We make the decision to break or not -after- the calculation of Ekin and Pressure,
               so scroll down for that logic */

            /* #########   START SECOND UPDATE STEP ################# */
            GMX_MPE_LOG(ev_update_start);
            bOK = TRUE;
            if (!bRerunMD || rerun_fr.bV || bForceUpdate)
            {
                wallcycle_start(wcycle,ewcUPDATE);
                dvdl = 0;
                /* Box is changed in update() when we do pressure coupling,
                 * but we should still use the old box for energy corrections and when
                 * writing it to the energy file, so it matches the trajectory files for
                 * the same timestep above. Make a copy in a separate array.
                 */
                copy_mat(state->box,lastbox);
                /* UPDATE PRESSURE VARIABLES IN TROTTER FORMULATION WITH CONSTRAINTS */
                if (bTrotter)
                {
                    if (bIterations && iterate.bIterate)
                    {
                        if (bFirstIterate)
                        {
                            scalevir = 1;
                        }
                        else
                        {
                            /* we use a new value of scalevir to converge the iterations faster */
                            scalevir = tracevir/trace(shake_vir);
                        }
                        msmul(shake_vir,scalevir,shake_vir);
                        m_add(force_vir,shake_vir,total_vir);
                        clear_mat(shake_vir);
                    }
                    trotter_update(ir,step,ekind,enerd,state,total_vir,mdatoms,&MassQ, trotter_seq,ettTSEQ3);
                }
                /* We can only do Berendsen coupling after we have summed
                 * the kinetic energy or virial. Since the happens
                 * in global_state after update, we should only do it at
                 * step % nstlist = 1 with bGStatEveryStep=FALSE.
                 */

		if (ir->eI != eiVVAK)
                {
		  update_tcouple(fplog,step,ir,state,ekind,wcycle,upd,&MassQ,mdatoms);
                }
                update_pcouple(fplog,step,ir,state,pcoupl_mu,M,wcycle,
                                upd,bInitStep);

		if (bVV)
		{
		    /* velocity half-step update */
		    update_coords(fplog,step,ir,mdatoms,state,f,fr->bTwinRange && bNStList,fr->f_twin,fcd,
				  ekind,M,wcycle,upd,FALSE,etrtVELOCITY2,cr,nrnb,constr,&top->idef);
		}

                /* Above, initialize just copies ekinh into ekin,
                 * it doesn't copy position (for VV),
                 * and entire integrator for MD.
                 */

                if (ir->eI==eiVVAK)
                {
                    copy_rvecn(state->x,cbuf,0,state->natoms);
                }

                update_coords(fplog,step,ir,mdatoms,state,f,fr->bTwinRange && bNStList,fr->f_twin,fcd,
                              ekind,M,wcycle,upd,bInitStep,etrtPOSITION,cr,nrnb,constr,&top->idef);
                wallcycle_stop(wcycle,ewcUPDATE);

                update_constraints(fplog,step,&dvdl,ir,ekind,mdatoms,state,graph,f,
                                   &top->idef,shake_vir,force_vir,
                                   cr,nrnb,wcycle,upd,constr,
                                   bInitStep,FALSE,bCalcEnerPres,state->veta);

                if (ir->eI==eiVVAK)
                {
                    /* erase F_EKIN and F_TEMP here? */
                    /* just compute the kinetic energy at the half step to perform a trotter step */
                    compute_globals(fplog,gstat,cr,ir,fr,ekind,state,state_global,mdatoms,nrnb,vcm,
                                    wcycle,enerd,force_vir,shake_vir,total_vir,pres,mu_tot,
                                    constr,NULL,FALSE,lastbox,
                                    top_global,&pcurr,top_global->natoms,&bSumEkinhOld,
                                    cglo_flags | CGLO_TEMPERATURE | CGLO_CONSTRAINT
                        );
                    wallcycle_start(wcycle,ewcUPDATE);
                    trotter_update(ir,step,ekind,enerd,state,total_vir,mdatoms,&MassQ, trotter_seq,ettTSEQ4);
                    /* now we know the scaling, we can compute the positions again again */
                    copy_rvecn(cbuf,state->x,0,state->natoms);

                    update_coords(fplog,step,ir,mdatoms,state,f,fr->bTwinRange && bNStList,fr->f_twin,fcd,
                                  ekind,M,wcycle,upd,bInitStep,etrtPOSITION,cr,nrnb,constr,&top->idef);
                    wallcycle_stop(wcycle,ewcUPDATE);

                    /* do we need an extra constraint here? just need to copy out of state->v to upd->xp? */
                    /* are the small terms in the shake_vir here due
                     * to numerical errors, or are they important
                     * physically? I'm thinking they are just errors, but not completely sure.
                     * For now, will call without actually constraining, constr=NULL*/
                    update_constraints(fplog,step,&dvdl,ir,ekind,mdatoms,state,graph,f,
                                       &top->idef,tmp_vir,force_vir,
                                       cr,nrnb,wcycle,upd,NULL,
                                       bInitStep,FALSE,bCalcEnerPres,state->veta);
                }
                if (!bOK && !bFFscan)
                {
                    gmx_fatal(FARGS,"Constraint error: Shake, Lincs or Settle could not solve the constrains");
                }

                if (fr->bSepDVDL && fplog && do_log)
                {
                    fprintf(fplog,sepdvdlformat,"Constraint",0.0,dvdl);
                }
                enerd->term[F_DHDL_CON] += dvdl;
            }
            else if (graph)
            {
                /* Need to unshift here */
                unshift_self(graph,state->box,state->x);
            }

            GMX_BARRIER(cr->mpi_comm_mygroup);
            GMX_MPE_LOG(ev_update_finish);

            if (vsite != NULL)
            {
                wallcycle_start(wcycle,ewcVSITECONSTR);
                if (graph != NULL)
                {
                    shift_self(graph,state->box,state->x);
                }
                construct_vsites(fplog,vsite,state->x,nrnb,ir->delta_t,state->v,
                                 top->idef.iparams,top->idef.il,
                                 fr->ePBC,fr->bMolPBC,graph,cr,state->box);

                if (graph != NULL)
                {
                    unshift_self(graph,state->box,state->x);
                }
                wallcycle_stop(wcycle,ewcVSITECONSTR);
            }

            /* ############## IF NOT VV, Calculate globals HERE, also iterate constraints ############ */
            if (ir->nstlist == -1 && bFirstIterate)
            {
                gs.sig[eglsNABNSB] = nlh.nabnsb;
            }
            compute_globals(fplog,gstat,cr,ir,fr,ekind,state,state_global,mdatoms,nrnb,vcm,
                            wcycle,enerd,force_vir,shake_vir,total_vir,pres,mu_tot,
                            constr,
                            bFirstIterate ? &gs : NULL,(step % gs.nstms == 0),
                            lastbox,
                            top_global,&pcurr,top_global->natoms,&bSumEkinhOld,
                            cglo_flags
                            | (!EI_VV(ir->eI) ? CGLO_ENERGY : 0)
                            | (!EI_VV(ir->eI) ? CGLO_TEMPERATURE : 0)
                            | (!EI_VV(ir->eI) || bRerunMD ? CGLO_PRESSURE : 0)
                            | (bIterations && iterate.bIterate ? CGLO_ITERATE : 0)
                            | (bFirstIterate ? CGLO_FIRSTITERATE : 0)
                            | CGLO_CONSTRAINT
                );
            if (ir->nstlist == -1 && bFirstIterate)
            {
                nlh.nabnsb = gs.set[eglsNABNSB];
                gs.set[eglsNABNSB] = 0;
            }
            /* bIterate is set to keep it from eliminating the old ekin kinetic energy terms */
            /* #############  END CALC EKIN AND PRESSURE ################# */

            /* Note: this is OK, but there are some numerical precision issues with using the convergence of
               the virial that should probably be addressed eventually. state->veta has better properies,
               but what we actually need entering the new cycle is the new shake_vir value. Ideally, we could
               generate the new shake_vir, but test the veta value for convergence.  This will take some thought. */

            if (bIterations &&
                done_iterating(cr,fplog,step,&iterate,bFirstIterate,
                               trace(shake_vir),&tracevir))
            {
                break;
            }
            bFirstIterate = FALSE;
        }

        update_box(fplog,step,ir,mdatoms,state,graph,f,
                   ir->nstlist==-1 ? &nlh.scale_tot : NULL,pcoupl_mu,nrnb,wcycle,upd,bInitStep,FALSE);

        /* ################# END UPDATE STEP 2 ################# */
        /* #### We now have r(t+dt) and v(t+dt/2)  ############# */

        /* The coordinates (x) were unshifted in update */
/*        if (bFFscan && (shellfc==NULL || bConverged))
        {
            if (print_forcefield(fplog,enerd->term,mdatoms->homenr,
                                 f,NULL,xcopy,
                                 &(top_global->mols),mdatoms->massT,pres))
            {
                if (gmx_parallel_env_initialized())
                {
                    gmx_finalize();
                }
                fprintf(stderr,"\n");
                exit(0);
            }
        }*/
        if (!bGStat)
        {
            /* We will not sum ekinh_old,
             * so signal that we still have to do it.
             */
            bSumEkinhOld = TRUE;
        }

/*        if (bTCR)
        {*/
            /* Only do GCT when the relaxation of shells (minimization) has converged,
             * otherwise we might be coupling to bogus energies.
             * In parallel we must always do this, because the other sims might
             * update the FF.
             */

            /* Since this is called with the new coordinates state->x, I assume
             * we want the new box state->box too. / EL 20040121
             */
/*            do_coupling(fplog,oenv,nfile,fnm,tcr,t,step,enerd->term,fr,
                        ir,MASTER(cr),
                        mdatoms,&(top->idef),mu_aver,
                        top_global->mols.nr,cr,
                        state->box,total_vir,pres,
                        mu_tot,state->x,f,bConverged);
            debug_gmx();
        }*/

        /* #########  BEGIN PREPARING EDR OUTPUT  ###########  */

        sum_dhdl(enerd,state->lambda,ir);
        /* use the directly determined last velocity, not actually the averaged half steps */
        if (bTrotter && ir->eI==eiVV)
        {
            enerd->term[F_EKIN] = last_ekin;
        }
        enerd->term[F_ETOT] = enerd->term[F_EPOT] + enerd->term[F_EKIN];

        switch (ir->etc)
        {
        case etcNO:
            break;
        case etcBERENDSEN:
            break;
        case etcNOSEHOOVER:
            if (IR_NVT_TROTTER(ir)) {
                enerd->term[F_ECONSERVED] = enerd->term[F_ETOT] + last_conserved;
            } else {
                enerd->term[F_ECONSERVED] = enerd->term[F_ETOT] +
                    NPT_energy(ir,state,&MassQ);
            }
            break;
        case etcVRESCALE:
            enerd->term[F_ECONSERVED] =
                enerd->term[F_ETOT] + vrescale_energy(&(ir->opts),
                                                      state->therm_integral);
            break;
        default:
            break;
        }

        /* Check for excessively large energies */
/*        if (bIonize)
        {
#ifdef GMX_DOUBLE
            real etot_max = 1e200;
#else
            real etot_max = 1e30;
#endif
            if (fabs(enerd->term[F_ETOT]) > etot_max)
            {
                fprintf(stderr,"Energy too large (%g), giving up\n",
                        enerd->term[F_ETOT]);
            }
        }*/
        /* #########  END PREPARING EDR OUTPUT  ###########  */

        /* Time for performance */
        if (((step % stepout) == 0) || bLastStep)
        {
            runtime_upd_proc(runtime);
        }

        /* Output stuff */
        if (MASTER(cr))
        {
            gmx_bool do_dr,do_or;

            if (!(bStartingFromCpt && (EI_VV(ir->eI))))
            {
                if (bNstEner)
                {
                    upd_mdebin(mdebin,bDoDHDL,TRUE,
                               t,mdatoms->tmass,enerd,state,lastbox,
                               shake_vir,force_vir,total_vir,pres,
                               ekind,mu_tot,constr);
                }
                else
                {
                    upd_mdebin_step(mdebin);
                }

                do_dr  = do_per_step(step,ir->nstdisreout);
                do_or  = do_per_step(step,ir->nstorireout);

                print_ebin(outf->fp_ene,do_ene,do_dr,do_or,do_log?fplog:NULL,
                           step,t,
                           eprNORMAL,bCompact,mdebin,fcd,groups,&(ir->opts));
            }
            if (ir->ePull != epullNO)
            {
                pull_print_output(ir->pull,step,t);
            }

            if (do_per_step(step,ir->nstlog))
            {
                if(fflush(fplog) != 0)
                {
                    gmx_fatal(FARGS,"Cannot flush logfile - maybe you are out of disk space?");
                }
            }
        }


        /* Remaining runtime */
        if (MULTIMASTER(cr) && (do_verbose || gmx_got_usr_signal() ))
        {
            if (shellfc)
            {
                fprintf(stderr,"\n");
            }
            print_time(stderr,runtime,step,ir,cr);
        }

		/* Set new positions for the group to embed */
		if(!bLastStep){
			if(step_rel<=it_xy)
			{
				fac[0]+=xy_step;
				fac[1]+=xy_step;
			} else if (step_rel<=(it_xy+it_z))
			{
				fac[2]+=z_step;
			}
			resize(ins_at,r_ins,state_global->x,pos_ins,fac);
		}

        /* Replica exchange */
/*        bExchanged = FALSE;
        if ((repl_ex_nst > 0) && (step > 0) && !bLastStep &&
            do_per_step(step,repl_ex_nst))
        {
            bExchanged = replica_exchange(fplog,cr,repl_ex,
                                          state_global,enerd->term,
                                          state,step,t);
        }
        if (bExchanged && PAR(cr))
        {
            if (DOMAINDECOMP(cr))
            {
                dd_partition_system(fplog,step,cr,TRUE,1,
                                    state_global,top_global,ir,
                                    state,&f,mdatoms,top,fr,
                                    vsite,shellfc,constr,
                                    nrnb,wcycle,FALSE);
            }
            else
            {
                bcast_state(cr,state,FALSE);
            }
        }*/

        bFirstStep = FALSE;
        bInitStep = FALSE;
        bStartingFromCpt = FALSE;

        /* #######  SET VARIABLES FOR NEXT ITERATION IF THEY STILL NEED IT ###### */
	/* With all integrators, except VV, we need to retain the pressure
         * at the current step for coupling at the next step.
         */
        if ((state->flags & (1<<estPRES_PREV)) &&
            (bGStatEveryStep ||
             (ir->nstpcouple > 0 && step % ir->nstpcouple == 0)))
        {
            /* Store the pressure in t_state for pressure coupling
             * at the next MD step.
             */
            copy_mat(pres,state->pres_prev);
        }

        /* #######  END SET VARIABLES FOR NEXT ITERATION ###### */

        if (bRerunMD)
        {
            /* read next frame from input trajectory */
            bNotLastFrame = read_next_frame(oenv,status,&rerun_fr);
        }

        if (!bRerunMD || !rerun_fr.bStep)
        {
            /* increase the MD step number */
            step++;
            step_rel++;
        }

        cycles = wallcycle_stop(wcycle,ewcSTEP);
        if (DOMAINDECOMP(cr) && wcycle)
        {
            dd_cycles_add(cr->dd,cycles,ddCyclStep);
        }

        if (step_rel == wcycle_get_reset_counters(wcycle) ||
            gs.set[eglsRESETCOUNTERS] != 0)
        {
            /* Reset all the counters related to performance over the run */
            reset_all_counters(fplog,cr,step,&step_rel,ir,wcycle,nrnb,runtime);
            wcycle_set_reset_counters(wcycle,-1);
            bResetCountersHalfMaxH = FALSE;
            gs.set[eglsRESETCOUNTERS] = 0;
        }
    }
    /* End of main MD loop */
    debug_gmx();
    write_sto_conf_mtop(ftp2fn(efSTO,nfile,fnm),
                                        *top_global->name,top_global,
                                        state_global->x,state_global->v,
                                        ir->ePBC,state->box);

    /* Stop the time */
    runtime_end(runtime);

    if (bRerunMD)
    {
        close_trj(status);
    }

    if (!(cr->duty & DUTY_PME))
    {
        /* Tell the PME only node to finish */
        gmx_pme_finish(cr);
    }

    if (MASTER(cr))
    {
        if (ir->nstcalcenergy > 0 && !bRerunMD)
        {
            print_ebin(outf->fp_ene,FALSE,FALSE,FALSE,fplog,step,t,
                       eprAVER,FALSE,mdebin,fcd,groups,&(ir->opts));
        }
    }

    done_mdoutf(outf);

    debug_gmx();

    if (ir->nstlist == -1 && nlh.nns > 0 && fplog)
    {
        fprintf(fplog,"Average neighborlist lifetime: %.1f steps, std.dev.: %.1f steps\n",nlh.s1/nlh.nns,sqrt(nlh.s2/nlh.nns - sqr(nlh.s1/nlh.nns)));
        fprintf(fplog,"Average number of atoms that crossed the half buffer length: %.1f\n\n",nlh.ab/nlh.nns);
    }

    if (shellfc && fplog)
    {
        fprintf(fplog,"Fraction of iterations that converged:           %.2f %%\n",
                (nconverged*100.0)/step_rel);
        fprintf(fplog,"Average number of force evaluations per MD step: %.2f\n\n",
                tcount/step_rel);
    }

/*    if (repl_ex_nst > 0 && MASTER(cr))
    {
        print_replica_exchange_statistics(fplog,repl_ex);
    }*/

    runtime->nsteps_done = step_rel;

    return 0;
}


int mdrunner_membed(FILE *fplog,t_commrec *cr,int nfile,const t_filenm fnm[],
             const output_env_t oenv, gmx_bool bVerbose,gmx_bool bCompact,
             int nstglobalcomm,
             ivec ddxyz,int dd_node_order,real rdd,real rconstr,
             const char *dddlb_opt,real dlb_scale,
             const char *ddcsx,const char *ddcsy,const char *ddcsz,
             int nstepout,int resetstep,int nmultisim,int repl_ex_nst,int repl_ex_seed,
             real pforce,real cpt_period,real max_hours,
             const char *deviceOptions,
             unsigned long Flags,
             real xy_fac, real xy_max, real z_fac, real z_max,
             int it_xy, int it_z, real probe_rad, int low_up_rm,
             int pieces, gmx_bool bALLOW_ASYMMETRY, int maxwarn)
{
    double     nodetime=0,realtime;
    t_inputrec *inputrec;
    t_state    *state=NULL;
    matrix     box;
    gmx_ddbox_t ddbox;
    int        npme_major,npme_minor;
    real       tmpr1,tmpr2;
    t_nrnb     *nrnb;
    gmx_mtop_t *mtop=NULL;
    t_mdatoms  *mdatoms=NULL;
    t_forcerec *fr=NULL;
    t_fcdata   *fcd=NULL;
    real       ewaldcoeff=0;
    gmx_pme_t  *pmedata=NULL;
    gmx_vsite_t *vsite=NULL;
    gmx_constr_t constr;
    int        i,m,nChargePerturbed=-1,status,nalloc;
    char       *gro;
    gmx_wallcycle_t wcycle;
    gmx_bool       bReadRNG,bReadEkin;
    int        list;
    gmx_runtime_t runtime;
    int        rc;
    gmx_large_int_t reset_counters;
    gmx_edsam_t ed=NULL;
    t_commrec   *cr_old=cr;
    int        nthreads=1,nthreads_requested=1;


	char			*ins;
	int 			rm_bonded_at,fr_id,fr_i=0,tmp_id,warn=0;
	int        		ng,j,max_lip_rm,ins_grp_id,ins_nat,mem_nat,ntype,lip_rm,tpr_version;
	real			xy_step=0,z_step=0;
	real		 	prot_area;
	rvec			*r_ins=NULL,fac;
	t_block 		*ins_at,*rest_at;
	pos_ins_t 		*pos_ins;
	mem_t			*mem_p;
	rmm_t			*rm_p;
	gmx_groups_t 		*groups;
	gmx_bool		 	bExcl=FALSE;
	t_atoms			atoms;
	t_pbc			*pbc;
	char		        **piecename=NULL;

    /* CAUTION: threads may be started later on in this function, so
       cr doesn't reflect the final parallel state right now */
    snew(inputrec,1);
    snew(mtop,1);

    if (bVerbose && SIMMASTER(cr))
    {
        fprintf(stderr,"Getting Loaded...\n");
    }

    if (Flags & MD_APPENDFILES)
    {
        fplog = NULL;
    }

    snew(state,1);
    if (MASTER(cr))
    {
        /* Read (nearly) all data required for the simulation */
        read_tpx_state(ftp2fn(efTPX,nfile,fnm),inputrec,state,NULL,mtop);

        /* NOW the threads will be started: */
#ifdef GMX_THREADS
#endif
    }
    /* END OF CAUTION: cr is now reliable */

    if (PAR(cr))
    {
        /* now broadcast everything to the non-master nodes/threads: */
        init_parallel(fplog, cr, inputrec, mtop);
    }
    /* now make sure the state is initialized and propagated */
    set_state_entries(state,inputrec,cr->nnodes);

    if (can_use_allvsall(inputrec,mtop,TRUE,cr,fplog))
    {
        /* All-vs-all loops do not work with domain decomposition */
        Flags |= MD_PARTDEC;
    }

    if (!EEL_PME(inputrec->coulombtype) || (Flags & MD_PARTDEC))
    {
        cr->npmenodes = 0;
    }

	snew(ins_at,1);
	snew(pos_ins,1);
	if(MASTER(cr))
	{
		tpr_version = get_tpr_version(ftp2fn(efTPX,nfile,fnm));
		if (tpr_version<58)
			gmx_fatal(FARGS,"Version of *.tpr file to old (%d). Rerun grompp with gromacs VERSION 4.0.3 or newer.\n",tpr_version);

		if( inputrec->eI != eiMD )
			gmx_input("Change integrator to md in mdp file.");

		if(PAR(cr))
			gmx_input("Sorry, parallel g_membed is not yet fully functrional.");

		groups=&(mtop->groups);

		atoms=gmx_mtop_global_atoms(mtop);
		snew(mem_p,1);
		fprintf(stderr,"\nSelect a group to embed in the membrane:\n");
		get_index(&atoms,ftp2fn_null(efNDX,nfile,fnm),1,&(ins_at->nr),&(ins_at->index),&ins);
		ins_grp_id = search_string(ins,groups->ngrpname,(groups->grpname));
		fprintf(stderr,"\nSelect a group to embed %s into (e.g. the membrane):\n",ins);
		get_index(&atoms,ftp2fn_null(efNDX,nfile,fnm),1,&(mem_p->mem_at.nr),&(mem_p->mem_at.index),&(mem_p->name));

		pos_ins->pieces=pieces;
		snew(pos_ins->nidx,pieces);
		snew(pos_ins->subindex,pieces);
		snew(piecename,pieces);	
		if (pieces>1)
		{
			fprintf(stderr,"\nSelect pieces to embed:\n");
			get_index(&atoms,ftp2fn_null(efNDX,nfile,fnm),pieces,pos_ins->nidx,pos_ins->subindex,piecename);
		}
		else
		{	
			/*use whole embedded group*/
			snew(pos_ins->nidx,1);
			snew(pos_ins->subindex,1);
			pos_ins->nidx[0]=ins_at->nr;
			pos_ins->subindex[0]=ins_at->index;
		}

		if(probe_rad<0.2199999)
		{
			warn++;
			fprintf(stderr,"\nWarning %d:\nA probe radius (-rad) smaller than 0.2 can result in overlap between waters "
					"and the group to embed, which will result in Lincs errors etc.\nIf you are sure, you can increase maxwarn.\n\n",warn);
		}

		if(xy_fac<0.09999999)
		{
			warn++;
			fprintf(stderr,"\nWarning %d:\nThe initial size of %s is probably too smal.\n"
					"If you are sure, you can increase maxwarn.\n\n",warn,ins);
		}

		if(it_xy<1000)
		{
			warn++;
			fprintf(stderr,"\nWarning %d;\nThe number of steps used to grow the xy-coordinates of %s (%d) is probably too small.\n"
					"Increase -nxy or, if you are sure, you can increase maxwarn.\n\n",warn,ins,it_xy);
		}

		if( (it_z<100) && ( z_fac<0.99999999 || z_fac>1.0000001) )
                {
                        warn++;
                        fprintf(stderr,"\nWarning %d;\nThe number of steps used to grow the z-coordinate of %s (%d) is probably too small.\n"
                                       "Increase -nz or, if you are sure, you can increase maxwarn.\n\n",warn,ins,it_z);
                }

		if(it_xy+it_z>inputrec->nsteps)
		{
			warn++;
			fprintf(stderr,"\nWarning %d:\nThe number of growth steps (-nxy + -nz) is larger than the number of steps in the tpr.\n"
					"If you are sure, you can increase maxwarn.\n\n",warn);
		}

		fr_id=-1;
		if( inputrec->opts.ngfrz==1)
			gmx_fatal(FARGS,"You did not specify \"%s\" as a freezegroup.",ins);
		for(i=0;i<inputrec->opts.ngfrz;i++)
		{
			tmp_id = mtop->groups.grps[egcFREEZE].nm_ind[i];
			if(ins_grp_id==tmp_id)
			{
				fr_id=tmp_id;
				fr_i=i;
			}
		}
		if (fr_id == -1 )
			gmx_fatal(FARGS,"\"%s\" not as freezegroup defined in the mdp-file.",ins);

		for(i=0;i<DIM;i++)
			if( inputrec->opts.nFreeze[fr_i][i] != 1)
				gmx_fatal(FARGS,"freeze dimensions for %s are not Y Y Y\n",ins);

		ng = groups->grps[egcENER].nr;
		if (ng == 1)
			gmx_input("No energy groups defined. This is necessary for energy exclusion in the freeze group");

		for(i=0;i<ng;i++)
		{
			for(j=0;j<ng;j++)
			{
				if (inputrec->opts.egp_flags[ng*i+j] == EGP_EXCL)
				{
					bExcl = TRUE;
					if ( (groups->grps[egcENER].nm_ind[i] != ins_grp_id) || (groups->grps[egcENER].nm_ind[j] != ins_grp_id) )
						gmx_fatal(FARGS,"Energy exclusions \"%s\" and  \"%s\" do not match the group to embed \"%s\"",
								*groups->grpname[groups->grps[egcENER].nm_ind[i]],
								*groups->grpname[groups->grps[egcENER].nm_ind[j]],ins);
				}
			}
		}
		if (!bExcl)
			gmx_input("No energy exclusion groups defined. This is necessary for energy exclusion in the freeze group");

		/* Set all atoms in box*/
		/*set_inbox(state->natoms,state->x);*/

		/* Guess the area the protein will occupy in the membrane plane	 Calculate area per lipid*/
		snew(rest_at,1);
		ins_nat = init_ins_at(ins_at,rest_at,state,pos_ins,groups,ins_grp_id,xy_max);
		/* Check moleculetypes in insertion group */
		check_types(ins_at,rest_at,mtop);

		mem_nat = init_mem_at(mem_p,mtop,state->x,state->box,pos_ins);

		prot_area = est_prot_area(pos_ins,state->x,ins_at,mem_p);
		if ( (prot_area>7.5) && ( (state->box[XX][XX]*state->box[YY][YY]-state->box[XX][YY]*state->box[YY][XX])<50) )
		{
			warn++;
			fprintf(stderr,"\nWarning %d:\nThe xy-area is very small compared to the area of the protein.\n"
					"This might cause pressure problems during the growth phase. Just try with\n"
					"current setup (-maxwarn + 1), but if pressure problems occur, lower the\n"
					"compressibility in the mdp-file or use no pressure coupling at all.\n\n",warn);
		}
		if(warn>maxwarn)
					gmx_fatal(FARGS,"Too many warnings.\n");

		printf("The estimated area of the protein in the membrane is %.3f nm^2\n",prot_area);
		printf("\nThere are %d lipids in the membrane part that overlaps the protein.\nThe area per lipid is %.4f nm^2.\n",mem_p->nmol,mem_p->lip_area);

		/* Maximum number of lipids to be removed*/
		max_lip_rm=(int)(2*prot_area/mem_p->lip_area);
		printf("Maximum number of lipids that will be removed is %d.\n",max_lip_rm);

		printf("\nWill resize the protein by a factor of %.3f in the xy plane and %.3f in the z direction.\n"
				"This resizing will be done with respect to the geometrical center of all protein atoms\n"
				"that span the membrane region, i.e. z between %.3f and %.3f\n\n",xy_fac,z_fac,mem_p->zmin,mem_p->zmax);

		/* resize the protein by xy and by z if necessary*/
		snew(r_ins,ins_at->nr);
		init_resize(ins_at,r_ins,pos_ins,mem_p,state->x,bALLOW_ASYMMETRY);
		fac[0]=fac[1]=xy_fac;
		fac[2]=z_fac;

		xy_step =(xy_max-xy_fac)/(double)(it_xy);
		z_step  =(z_max-z_fac)/(double)(it_z-1);

		resize(ins_at,r_ins,state->x,pos_ins,fac);

		/* remove overlapping lipids and water from the membrane box*/
		/*mark molecules to be removed*/
		snew(pbc,1);
		set_pbc(pbc,inputrec->ePBC,state->box);

		snew(rm_p,1);
		lip_rm = gen_rm_list(rm_p,ins_at,rest_at,pbc,mtop,state->x, r_ins, mem_p,pos_ins,probe_rad,low_up_rm,bALLOW_ASYMMETRY);
        lip_rm -= low_up_rm;

		if(fplog)
			for(i=0;i<rm_p->nr;i++)
				fprintf(fplog,"rm mol %d\n",rm_p->mol[i]);

		for(i=0;i<mtop->nmolblock;i++)
		{
			ntype=0;
			for(j=0;j<rm_p->nr;j++)
				if(rm_p->block[j]==i)
					ntype++;
			printf("Will remove %d %s molecules\n",ntype,*(mtop->moltype[mtop->molblock[i].type].name));
		}

		if(lip_rm>max_lip_rm)
		{
			warn++;
			fprintf(stderr,"\nWarning %d:\nTrying to remove a larger lipid area than the estimated protein area\n"
					"Try making the -xyinit resize factor smaller. If you are sure about this increase maxwarn.\n\n",warn);
		}

		/*remove all lipids and waters overlapping and update all important structures*/
		rm_group(inputrec,groups,mtop,rm_p,state,ins_at,pos_ins);

		rm_bonded_at = rm_bonded(ins_at,mtop);
		if (rm_bonded_at != ins_at->nr)
		{
			fprintf(stderr,"Warning: The number of atoms for which the bonded interactions are removed is %d, "
					"while %d atoms are embedded. Make sure that the atoms to be embedded are not in the same"
					"molecule type as atoms that are not to be embedded.\n",rm_bonded_at,ins_at->nr);
		}

		if(warn>maxwarn)
			gmx_fatal(FARGS,"Too many warnings.\nIf you are sure these warnings are harmless, you can increase -maxwarn");

		if (MASTER(cr))
		{
			if (ftp2bSet(efTOP,nfile,fnm))
				top_update(opt2fn("-p",nfile,fnm),ins,rm_p,mtop);
		}

		sfree(pbc);
		sfree(rest_at);
	}

#ifdef GMX_FAHCORE
    fcRegisterSteps(inputrec->nsteps,inputrec->init_step);
#endif

    /* NMR restraints must be initialized before load_checkpoint,
     * since with time averaging the history is added to t_state.
     * For proper consistency check we therefore need to extend
     * t_state here.
     * So the PME-only nodes (if present) will also initialize
     * the distance restraints.
     */
    snew(fcd,1);

    /* This needs to be called before read_checkpoint to extend the state */
    init_disres(fplog,mtop,inputrec,cr,Flags & MD_PARTDEC,fcd,state,FALSE);

    if (gmx_mtop_ftype_count(mtop,F_ORIRES) > 0)
    {
        if (PAR(cr) && !(Flags & MD_PARTDEC))
        {
            gmx_fatal(FARGS,"Orientation restraints do not work (yet) with domain decomposition, use particle decomposition (mdrun option -pd)");
        }
        /* Orientation restraints */
        if (MASTER(cr))
        {
            init_orires(fplog,mtop,state->x,inputrec,cr->ms,&(fcd->orires),
                        state);
        }
    }

    if (DEFORM(*inputrec))
    {
        /* Store the deform reference box before reading the checkpoint */
        if (SIMMASTER(cr))
        {
            copy_mat(state->box,box);
        }
        if (PAR(cr))
        {
            gmx_bcast(sizeof(box),box,cr);
        }
        /* Because we do not have the update struct available yet
         * in which the reference values should be stored,
         * we store them temporarily in static variables.
         * This should be thread safe, since they are only written once
         * and with identical values.
         */
/*        deform_init_init_step_tpx = inputrec->init_step;*/
/*        copy_mat(box,deform_init_box_tpx);*/
    }

    if (opt2bSet("-cpi",nfile,fnm))
    {
        /* Check if checkpoint file exists before doing continuation.
         * This way we can use identical input options for the first and subsequent runs...
         */
        if( gmx_fexist_master(opt2fn_master("-cpi",nfile,fnm,cr),cr) )
        {
            load_checkpoint(opt2fn_master("-cpi",nfile,fnm,cr),&fplog,
                            cr,Flags & MD_PARTDEC,ddxyz,
                            inputrec,state,&bReadRNG,&bReadEkin,
                            (Flags & MD_APPENDFILES),
			    (Flags & MD_APPENDFILESSET));

            if (bReadRNG)
            {
                Flags |= MD_READ_RNG;
            }
            if (bReadEkin)
            {
                Flags |= MD_READ_EKIN;
            }
        }
    }

    if ((MASTER(cr) || (Flags & MD_SEPPOT)) && (Flags & MD_APPENDFILES))
    {
        gmx_log_open(ftp2fn(efLOG,nfile,fnm),cr,!(Flags & MD_SEPPOT),
                             Flags,&fplog);
    }

    if (SIMMASTER(cr))
    {
        copy_mat(state->box,box);
    }

    if (PAR(cr))
    {
        gmx_bcast(sizeof(box),box,cr);
    }

    if (bVerbose && SIMMASTER(cr))
    {
        fprintf(stderr,"Loaded with Money\n\n");
    }

    if (PAR(cr) && !((Flags & MD_PARTDEC) || EI_TPI(inputrec->eI)))
    {
        cr->dd = init_domain_decomposition(fplog,cr,Flags,ddxyz,rdd,rconstr,
                                           dddlb_opt,dlb_scale,
                                           ddcsx,ddcsy,ddcsz,
                                           mtop,inputrec,
                                           box,state->x,
                                           &ddbox,&npme_major,&npme_minor);

        make_dd_communicators(fplog,cr,dd_node_order);

        /* Set overallocation to avoid frequent reallocation of arrays */
        set_over_alloc_dd(TRUE);
    }
    else
    {
        /* PME, if used, is done on all nodes with 1D decomposition */
        cr->npmenodes = 0;
        cr->duty = (DUTY_PP | DUTY_PME);
        npme_major = cr->nnodes;
        npme_minor = 1;

        if (inputrec->ePBC == epbcSCREW)
        {
            gmx_fatal(FARGS,
                      "pbc=%s is only implemented with domain decomposition",
                      epbc_names[inputrec->ePBC]);
        }
    }

    if (PAR(cr))
    {
        /* After possible communicator splitting in make_dd_communicators.
         * we can set up the intra/inter node communication.
         */
        gmx_setup_nodecomm(fplog,cr);
    }

    wcycle = wallcycle_init(fplog,resetstep,cr);
    if (PAR(cr))
    {
        /* Master synchronizes its value of reset_counters with all nodes
         * including PME only nodes */
        reset_counters = wcycle_get_reset_counters(wcycle);
        gmx_bcast_sim(sizeof(reset_counters),&reset_counters,cr);
        wcycle_set_reset_counters(wcycle, reset_counters);
    }


    snew(nrnb,1);
    if (cr->duty & DUTY_PP)
    {
        /* For domain decomposition we allocate dynamically
         * in dd_partition_system.
         */
        if (DOMAINDECOMP(cr))
        {
            bcast_state_setup(cr,state);
        }
        else
        {
            if (PAR(cr))
            {
                if (!MASTER(cr))
                {
                    snew(state,1);
                }
                bcast_state(cr,state,TRUE);
            }
        }

        /* Dihedral Restraints */
        if (gmx_mtop_ftype_count(mtop,F_DIHRES) > 0)
        {
            init_dihres(fplog,mtop,inputrec,fcd);
        }

        /* Initiate forcerecord */
        fr = mk_forcerec();
        init_forcerec(fplog,oenv,fr,fcd,inputrec,mtop,cr,box,FALSE,
                      opt2fn("-table",nfile,fnm),
                      opt2fn("-tablep",nfile,fnm),
                      opt2fn("-tableb",nfile,fnm),FALSE,pforce);

        /* version for PCA_NOT_READ_NODE (see md.c) */
        /*init_forcerec(fplog,fr,fcd,inputrec,mtop,cr,box,FALSE,
          "nofile","nofile","nofile",FALSE,pforce);
          */
        fr->bSepDVDL = ((Flags & MD_SEPPOT) == MD_SEPPOT);

        /* Initialize QM-MM */
        if(fr->bQMMM)
        {
            init_QMMMrec(cr,box,mtop,inputrec,fr);
        }

        /* Initialize the mdatoms structure.
         * mdatoms is not filled with atom data,
         * as this can not be done now with domain decomposition.
         */
        mdatoms = init_mdatoms(fplog,mtop,inputrec->efep!=efepNO);

        /* Initialize the virtual site communication */
        vsite = init_vsite(mtop,cr);

        calc_shifts(box,fr->shift_vec);

        /* With periodic molecules the charge groups should be whole at start up
         * and the virtual sites should not be far from their proper positions.
         */
        if (!inputrec->bContinuation && MASTER(cr) &&
            !(inputrec->ePBC != epbcNONE && inputrec->bPeriodicMols))
        {
            /* Make molecules whole at start of run */
            if (fr->ePBC != epbcNONE)
            {
                do_pbc_first_mtop(fplog,inputrec->ePBC,box,mtop,state->x);
            }
            if (vsite)
            {
                /* Correct initial vsite positions are required
                 * for the initial distribution in the domain decomposition
                 * and for the initial shell prediction.
                 */
                construct_vsites_mtop(fplog,vsite,mtop,state->x);
            }
        }

        /* Initiate PPPM if necessary */
        if (fr->eeltype == eelPPPM)
        {
            if (mdatoms->nChargePerturbed)
            {
                gmx_fatal(FARGS,"Free energy with %s is not implemented",
                          eel_names[fr->eeltype]);
            }
            status = gmx_pppm_init(fplog,cr,oenv,FALSE,TRUE,box,
                                   getenv("GMXGHAT"),inputrec, (Flags & MD_REPRODUCIBLE));
            if (status != 0)
            {
                gmx_fatal(FARGS,"Error %d initializing PPPM",status);
            }
        }

        if (EEL_PME(fr->eeltype))
        {
            ewaldcoeff = fr->ewaldcoeff;
            pmedata = &fr->pmedata;
        }
        else
        {
            pmedata = NULL;
        }
    }
    else
    {
        /* This is a PME only node */

        /* We don't need the state */
        done_state(state);

        ewaldcoeff = calc_ewaldcoeff(inputrec->rcoulomb, inputrec->ewald_rtol);
        snew(pmedata,1);
    }

    /* Initiate PME if necessary,
     * either on all nodes or on dedicated PME nodes only. */
    if (EEL_PME(inputrec->coulombtype))
    {
        if (mdatoms)
        {
            nChargePerturbed = mdatoms->nChargePerturbed;
        }
        if (cr->npmenodes > 0)
        {
            /* The PME only nodes need to know nChargePerturbed */
            gmx_bcast_sim(sizeof(nChargePerturbed),&nChargePerturbed,cr);
        }
        if (cr->duty & DUTY_PME)
        {
            status = gmx_pme_init(pmedata,cr,npme_major,npme_minor,inputrec,
                                  mtop ? mtop->natoms : 0,nChargePerturbed,
                                  (Flags & MD_REPRODUCIBLE));
            if (status != 0)
            {
                gmx_fatal(FARGS,"Error %d initializing PME",status);
            }
        }
    }


/*    if (integrator[inputrec->eI].func == do_md
#ifdef GMX_OPENMM
        ||
        integrator[inputrec->eI].func == do_md_openmm
#endif
        )
    {*/
        /* Turn on signal handling on all nodes */
        /*
         * (A user signal from the PME nodes (if any)
         * is communicated to the PP nodes.
         */
        signal_handler_install();
/*    }*/

    if (cr->duty & DUTY_PP)
    {
        if (inputrec->ePull != epullNO)
        {
            /* Initialize pull code */
            init_pull(fplog,inputrec,nfile,fnm,mtop,cr,oenv,
                      EI_DYNAMICS(inputrec->eI) && MASTER(cr),Flags);
        }

        constr = init_constraints(fplog,mtop,inputrec,ed,state,cr);

        if (DOMAINDECOMP(cr))
        {
            dd_init_bondeds(fplog,cr->dd,mtop,vsite,constr,inputrec,
                            Flags & MD_DDBONDCHECK,fr->cginfo_mb);

            set_dd_parameters(fplog,cr->dd,dlb_scale,inputrec,fr,&ddbox);

            setup_dd_grid(fplog,cr->dd);
        }

        /* Now do whatever the user wants us to do (how flexible...) */
        do_md_membed(fplog,cr,nfile,fnm,
                                      oenv,bVerbose,bCompact,
                                      nstglobalcomm,
                                      vsite,constr,
                                      nstepout,inputrec,mtop,
                                      fcd,state,
                                      mdatoms,nrnb,wcycle,ed,fr,
                                      repl_ex_nst,repl_ex_seed,
                                      cpt_period,max_hours,
                                      deviceOptions,
                                      Flags,
                                      &runtime,
                                      fac, r_ins, pos_ins, ins_at,
                                      xy_step, z_step, it_xy, it_z);

        if (inputrec->ePull != epullNO)
        {
            finish_pull(fplog,inputrec->pull);
        }
    }
    else
    {
        /* do PME only */
        gmx_pmeonly(*pmedata,cr,nrnb,wcycle,ewaldcoeff,FALSE,inputrec);
    }

    if (EI_DYNAMICS(inputrec->eI) || EI_TPI(inputrec->eI))
    {
        /* Some timing stats */
        if (MASTER(cr))
        {
            if (runtime.proc == 0)
            {
                runtime.proc = runtime.real;
            }
        }
        else
        {
            runtime.real = 0;
        }
    }

    wallcycle_stop(wcycle,ewcRUN);

    /* Finish up, write some stuff
     * if rerunMD, don't write last frame again
     */
    finish_run(fplog,cr,ftp2fn(efSTO,nfile,fnm),
               inputrec,nrnb,wcycle,&runtime,
               EI_DYNAMICS(inputrec->eI) && !MULTISIM(cr));

    /* Does what it says */
    print_date_and_time(fplog,cr->nodeid,"Finished mdrun",&runtime);

    /* Close logfile already here if we were appending to it */
    if (MASTER(cr) && (Flags & MD_APPENDFILES))
    {
        gmx_log_close(fplog);
    }

    if (pieces>1)
    {
    	sfree(piecename);
    }

    rc=(int)gmx_get_stop_condition();

    return rc;
}

int gmx_membed(int argc,char *argv[])
{
	const char *desc[] = {
			"[TT]g_membed[tt] embeds a membrane protein into an equilibrated lipid bilayer at the position",
			"and orientation specified by the user.[PAR]",
			"SHORT MANUAL[BR]------------[BR]",
			"The user should merge the structure files of the protein and membrane (+solvent), creating a",
			"single structure file with the protein overlapping the membrane at the desired position and",
			"orientation. The box size is taken from the membrane structure file. The corresponding topology",
			"files should also be merged. Consecutively, create a [TT].tpr[tt] file (input for [TT]g_membed[tt]) from these files,"
			"with the following options included in the [TT].mdp[tt] file.[BR]",
			" - [TT]integrator      = md[tt][BR]",
			" - [TT]energygrp       = Protein[tt] (or other group that you want to insert)[BR]",
			" - [TT]freezegrps      = Protein[tt][BR]",
			" - [TT]freezedim       = Y Y Y[tt][BR]",
			" - [TT]energygrp_excl  = Protein Protein[tt][BR]",
			"The output is a structure file containing the protein embedded in the membrane. If a topology",
			"file is provided, the number of lipid and ",
			"solvent molecules will be updated to match the new structure file.[BR]",
			"For a more extensive manual see Wolf et al, J Comp Chem 31 (2010) 2169-2174, Appendix.[PAR]",
			"SHORT METHOD DESCRIPTION[BR]",
			"------------------------[BR]",
			"1. The protein is resized around its center of mass by a factor [TT]-xy[tt] in the xy-plane",
			"(the membrane plane) and a factor [TT]-z[tt] in the [IT]z[it]-direction (if the size of the",
			"protein in the z-direction is the same or smaller than the width of the membrane, a",
			"[TT]-z[tt] value larger than 1 can prevent that the protein will be enveloped by the lipids).[BR]",
			"2. All lipid and solvent molecules overlapping with the resized protein are removed. All",
			"intraprotein interactions are turned off to prevent numerical issues for small values of [TT]-xy[tt]",
			" or [TT]-z[tt][BR]",
			"3. One md step is performed.[BR]",
			"4. The resize factor ([TT]-xy[tt] or [TT]-z[tt]) is incremented by a small amount ((1-xy)/nxy or (1-z)/nz) and the",
			"protein is resized again around its center of mass. The resize factor for the xy-plane",
			"is incremented first. The resize factor for the z-direction is not changed until the [TT]-xy[tt] factor",
			"is 1 (thus after [TT]-nxy[tt] iterations).[BR]",
			"5. Repeat step 3 and 4 until the protein reaches its original size ([TT]-nxy[tt] + [TT]-nz[tt] iterations).[BR]",
			"For a more extensive method description see Wolf et al, J Comp Chem, 31 (2010) 2169-2174.[PAR]",
			"NOTE[BR]----[BR]",
			" - Protein can be any molecule you want to insert in the membrane.[BR]",
			" - It is recommended to perform a short equilibration run after the embedding",
			"(see Wolf et al, J Comp Chem 31 (2010) 2169-2174), to re-equilibrate the membrane. Clearly",
			"protein equilibration might require longer.[PAR]"
	};
	t_commrec    *cr;
	t_filenm fnm[] = {
			{ efTPX, "-f",      "into_mem", ffREAD },
			{ efNDX, "-n",      "index",    ffOPTRD },
			{ efTOP, "-p",      "topol",    ffOPTRW },
			{ efTRN, "-o",      NULL,       ffWRITE },
			{ efXTC, "-x",      NULL,       ffOPTWR },
			{ efCPT, "-cpi",    NULL,       ffOPTRD },
			{ efCPT, "-cpo",    NULL,       ffOPTWR },
			{ efSTO, "-c",      "membedded",  ffWRITE },
			{ efEDR, "-e",      "ener",     ffWRITE },
			{ efLOG, "-g",      "md",       ffWRITE },
			{ efEDI, "-ei",     "sam",      ffOPTRD },
			{ efTRX, "-rerun",  "rerun",    ffOPTRD },
			{ efXVG, "-table",  "table",    ffOPTRD },
			{ efXVG, "-tablep", "tablep",   ffOPTRD },
			{ efXVG, "-tableb", "table",    ffOPTRD },
			{ efXVG, "-dhdl",   "dhdl",     ffOPTWR },
			{ efXVG, "-field",  "field",    ffOPTWR },
			{ efXVG, "-table",  "table",    ffOPTRD },
			{ efXVG, "-tablep", "tablep",   ffOPTRD },
			{ efXVG, "-tableb", "table",    ffOPTRD },
			{ efTRX, "-rerun",  "rerun",    ffOPTRD },
			{ efXVG, "-tpi",    "tpi",      ffOPTWR },
			{ efXVG, "-tpid",   "tpidist",  ffOPTWR },
			{ efEDI, "-ei",     "sam",      ffOPTRD },
			{ efEDO, "-eo",     "sam",      ffOPTWR },
			{ efGCT, "-j",      "wham",     ffOPTRD },
			{ efGCT, "-jo",     "bam",      ffOPTWR },
			{ efXVG, "-ffout",  "gct",      ffOPTWR },
			{ efXVG, "-devout", "deviatie", ffOPTWR },
			{ efXVG, "-runav",  "runaver",  ffOPTWR },
			{ efXVG, "-px",     "pullx",    ffOPTWR },
			{ efXVG, "-pf",     "pullf",    ffOPTWR },
			{ efMTX, "-mtx",    "nm",       ffOPTWR },
			{ efNDX, "-dn",     "dipole",   ffOPTWR },
                        { efRND, "-multidir",NULL,      ffOPTRDMULT}
	};
#define NFILE asize(fnm)

	/* Command line options ! */
	gmx_bool bCart        = FALSE;
	gmx_bool bPPPME       = FALSE;
	gmx_bool bPartDec     = FALSE;
	gmx_bool bDDBondCheck = TRUE;
	gmx_bool bDDBondComm  = TRUE;
	gmx_bool bVerbose     = FALSE;
	gmx_bool bCompact     = TRUE;
	gmx_bool bSepPot      = FALSE;
	gmx_bool bRerunVSite  = FALSE;
	gmx_bool bIonize      = FALSE;
	gmx_bool bConfout     = TRUE;
	gmx_bool bReproducible = FALSE;

	int  npme=-1;
	int  nmultisim=0;
	int  nstglobalcomm=-1;
	int  repl_ex_nst=0;
	int  repl_ex_seed=-1;
	int  nstepout=100;
	int  nthreads=0; /* set to determine # of threads automatically */
	int  resetstep=-1;

	rvec realddxyz={0,0,0};
	const char *ddno_opt[ddnoNR+1] =
	{ NULL, "interleave", "pp_pme", "cartesian", NULL };
	const char *dddlb_opt[] =
	{ NULL, "auto", "no", "yes", NULL };
	real rdd=0.0,rconstr=0.0,dlb_scale=0.8,pforce=-1;
	char *ddcsx=NULL,*ddcsy=NULL,*ddcsz=NULL;
	real cpt_period=15.0,max_hours=-1;
	gmx_bool bAppendFiles=TRUE,bAddPart=TRUE;
	gmx_bool bResetCountersHalfWay=FALSE;
	output_env_t oenv=NULL;
	const char *deviceOptions = "";

	real xy_fac = 0.5;
	real xy_max = 1.0;
	real z_fac = 1.0;
	real z_max = 1.0;
	int it_xy = 1000;
	int it_z = 0;
	real probe_rad = 0.22;
	int low_up_rm = 0;
	int maxwarn=0;
	int pieces=1;
    gmx_bool bALLOW_ASYMMETRY=FALSE;


/* arguments relevant to OPENMM only*/
#ifdef GMX_OPENMM
    gmx_input("g_membed not functional in openmm");
#endif

	t_pargs pa[] = {
			{ "-xyinit",   FALSE, etREAL,  {&xy_fac},   	"Resize factor for the protein in the xy dimension before starting embedding" },
			{ "-xyend",   FALSE, etREAL,  {&xy_max},   		"Final resize factor in the xy dimension" },
			{ "-zinit",    FALSE, etREAL,  {&z_fac},  		"Resize factor for the protein in the z dimension before starting embedding" },
			{ "-zend",    FALSE, etREAL,  {&z_max},    		"Final resize faction in the z dimension" },
			{ "-nxy",     FALSE,  etINT,  {&it_xy},         "Number of iteration for the xy dimension" },
			{ "-nz",      FALSE,  etINT,  {&it_z},          "Number of iterations for the z dimension" },
			{ "-rad",     FALSE, etREAL,  {&probe_rad},     "Probe radius to check for overlap between the group to embed and the membrane"},
			{ "-pieces",  FALSE,  etINT,  {&pieces},        "Perform piecewise resize. Select parts of the group to insert and resize these with respect to their own geometrical center." },
            { "-asymmetry",FALSE, etBOOL,{&bALLOW_ASYMMETRY}, "Allow asymmetric insertion, i.e. the number of lipids removed from the upper and lower leaflet will not be checked." },
            { "-ndiff" ,  FALSE, etINT, {&low_up_rm},       "Number of lipids that will additionally be removed from the lower (negative number) or upper (positive number) membrane leaflet." },
			{ "-maxwarn", FALSE, etINT, {&maxwarn},			"Maximum number of warning allowed" },
  { "-pd",      FALSE, etBOOL,{&bPartDec},
    "HIDDENUse particle decompostion" },
  { "-dd",      FALSE, etRVEC,{&realddxyz},
    "HIDDENDomain decomposition grid, 0 is optimize" },
  { "-nt",      FALSE, etINT, {&nthreads},
    "HIDDENNumber of threads to start (0 is guess)" },
  { "-npme",    FALSE, etINT, {&npme},
    "HIDDENNumber of separate nodes to be used for PME, -1 is guess" },
  { "-ddorder", FALSE, etENUM, {ddno_opt},
    "HIDDENDD node order" },
  { "-ddcheck", FALSE, etBOOL, {&bDDBondCheck},
    "HIDDENCheck for all bonded interactions with DD" },
  { "-ddbondcomm", FALSE, etBOOL, {&bDDBondComm},
    "HIDDENUse special bonded atom communication when [TT]-rdd[tt] > cut-off" },
  { "-rdd",     FALSE, etREAL, {&rdd},
    "HIDDENThe maximum distance for bonded interactions with DD (nm), 0 is determine from initial coordinates" },
  { "-rcon",    FALSE, etREAL, {&rconstr},
    "HIDDENMaximum distance for P-LINCS (nm), 0 is estimate" },
  { "-dlb",     FALSE, etENUM, {dddlb_opt},
    "HIDDENDynamic load balancing (with DD)" },
  { "-dds",     FALSE, etREAL, {&dlb_scale},
    "HIDDENMinimum allowed dlb scaling of the DD cell size" },
  { "-ddcsx",   FALSE, etSTR, {&ddcsx},
    "HIDDENThe DD cell sizes in x" },
  { "-ddcsy",   FALSE, etSTR, {&ddcsy},
    "HIDDENThe DD cell sizes in y" },
  { "-ddcsz",   FALSE, etSTR, {&ddcsz},
    "HIDDENThe DD cell sizes in z" },
  { "-gcom",    FALSE, etINT,{&nstglobalcomm},
    "HIDDENGlobal communication frequency" },
  { "-compact", FALSE, etBOOL,{&bCompact},
    "Write a compact log file" },
  { "-seppot",  FALSE, etBOOL, {&bSepPot},
    "HIDDENWrite separate V and dVdl terms for each interaction type and node to the log file(s)" },
  { "-pforce",  FALSE, etREAL, {&pforce},
    "HIDDENPrint all forces larger than this (kJ/mol nm)" },
  { "-reprod",  FALSE, etBOOL,{&bReproducible},
    "HIDDENTry to avoid optimizations that affect binary reproducibility" },
  { "-multi",   FALSE, etINT,{&nmultisim},
    "HIDDENDo multiple simulations in parallel" },
  { "-replex",  FALSE, etINT, {&repl_ex_nst},
    "HIDDENAttempt replica exchange periodically with this period (steps)" },
  { "-reseed",  FALSE, etINT, {&repl_ex_seed},
    "HIDDENSeed for replica exchange, -1 is generate a seed" },
  { "-rerunvsite", FALSE, etBOOL, {&bRerunVSite},
    "HIDDENRecalculate virtual site coordinates with [TT]-rerun[tt]" },
  { "-ionize",  FALSE, etBOOL,{&bIonize},
    "HIDDENDo a simulation including the effect of an X-Ray bombardment on your system" },
  { "-confout", TRUE, etBOOL, {&bConfout},
    "HIDDENWrite the last configuration with [TT]-c[tt] and force checkpointing at the last step" },
  { "-stepout", FALSE, etINT, {&nstepout},
    "HIDDENFrequency of writing the remaining runtime" },
  { "-resetstep", FALSE, etINT, {&resetstep},
    "HIDDENReset cycle counters after these many time steps" },
  { "-resethway", FALSE, etBOOL, {&bResetCountersHalfWay},
    "HIDDENReset the cycle counters after half the number of steps or halfway [TT]-maxh[tt]" },
  { "-v",       FALSE, etBOOL,{&bVerbose},
    "Be loud and noisy" },
  { "-maxh",   FALSE, etREAL, {&max_hours},
    "HIDDENTerminate after 0.99 times this time (hours)" },
  { "-cpt",     FALSE, etREAL, {&cpt_period},
    "HIDDENCheckpoint interval (minutes)" },
  { "-append",  FALSE, etBOOL, {&bAppendFiles},
    "HIDDENAppend to previous output files when continuing from checkpoint" },
  { "-addpart",  FALSE, etBOOL, {&bAddPart},
    "HIDDENAdd the simulation part number to all output files when continuing from checkpoint" },
	};
	gmx_edsam_t  ed;
	unsigned long Flags, PCA_Flags;
	ivec     ddxyz;
	int      dd_node_order;
	gmx_bool     HaveCheckpoint;
	FILE     *fplog,*fptest;
	int      sim_part,sim_part_fn;
	const char *part_suffix=".part";
	char     suffix[STRLEN];
	int      rc;
        char **multidir=NULL;

	cr = init_par(&argc,&argv);

	PCA_Flags = (PCA_KEEP_ARGS | PCA_NOEXIT_ON_ARGS | PCA_CAN_SET_DEFFNM
			| (MASTER(cr) ? 0 : PCA_QUIET));


	/* Comment this in to do fexist calls only on master
	 * works not with rerun or tables at the moment
	 * also comment out the version of init_forcerec in md.c
	 * with NULL instead of opt2fn
	 */
	/*
   if (!MASTER(cr))
   {
   PCA_Flags |= PCA_NOT_READ_NODE;
   }
	 */

	parse_common_args(&argc,argv,PCA_Flags, NFILE,fnm,asize(pa),pa,
			asize(desc),desc,0,NULL, &oenv);

	/* we set these early because they might be used in init_multisystem()
   Note that there is the potential for npme>nnodes until the number of
   threads is set later on, if there's thread parallelization. That shouldn't
   lead to problems. */
	dd_node_order = nenum(ddno_opt);
	cr->npmenodes = npme;

#ifdef GMX_THREADS
	/* now determine the number of threads automatically. The threads are
   only started at mdrunner_threads, though. */
	if (nthreads<1)
	{
		nthreads=tMPI_Thread_get_hw_number();
	}
#else
	nthreads=1;
#endif

        /* now check the -multi and -multidir option */
        if (opt2bSet("-multidir", NFILE, fnm))
        {
            int i;
            if (nmultisim > 0)
            {
                gmx_fatal(FARGS, "mdrun -multi and -multidir options are mutually     exclusive.");
            }
            nmultisim = opt2fns(&multidir, "-multidir", NFILE, fnm);
        }


	if (repl_ex_nst != 0 && nmultisim < 2)
		gmx_fatal(FARGS,"Need at least two replicas for replica exchange (option -multi)");

	if (nmultisim > 1) {
#ifndef GMX_THREADS
                gmx_bool bParFn = (multidir == NULL);
		init_multisystem(cr,nmultisim,multidir,NFILE,fnm,TRUE);
#else
		gmx_fatal(FARGS,"mdrun -multi is not supported with the thread library.Please compile GROMACS with MPI support");
#endif
	}

	/* Check if there is ANY checkpoint file available */
	sim_part    = 1;
	sim_part_fn = sim_part;
	if (opt2bSet("-cpi",NFILE,fnm))
	{
		bAppendFiles =
			read_checkpoint_simulation_part(opt2fn_master("-cpi", NFILE,fnm,cr),
							&sim_part_fn,NULL,cr,
							bAppendFiles,NFILE,fnm,
							part_suffix,&bAddPart);
		if (sim_part_fn==0 && MASTER(cr))
		{
			fprintf(stdout,"No previous checkpoint file present, assuming this is a new run.\n");
		}
		else
		{
			sim_part = sim_part_fn + 1;
		}
	}
	else
	{
		bAppendFiles = FALSE;
	}

	if (!bAppendFiles)
	{
		sim_part_fn = sim_part;
	}

	if (bAddPart && sim_part_fn > 1)
	{
		/* This is a continuation run, rename trajectory output files
       (except checkpoint files) */
		/* create new part name first (zero-filled) */
		sprintf(suffix,"%s%04d",part_suffix,sim_part_fn);

		add_suffix_to_output_names(fnm,NFILE,suffix);
		fprintf(stdout,"Checkpoint file is from part %d, new output files will be suffixed '%s'.\n",sim_part-1,suffix);
	}

	Flags = opt2bSet("-rerun",NFILE,fnm) ? MD_RERUN : 0;
	Flags = Flags | (bSepPot       ? MD_SEPPOT       : 0);
	Flags = Flags | (bIonize       ? MD_IONIZE       : 0);
	Flags = Flags | (bPartDec      ? MD_PARTDEC      : 0);
	Flags = Flags | (bDDBondCheck  ? MD_DDBONDCHECK  : 0);
	Flags = Flags | (bDDBondComm   ? MD_DDBONDCOMM   : 0);
	Flags = Flags | (bConfout      ? MD_CONFOUT      : 0);
	Flags = Flags | (bRerunVSite   ? MD_RERUN_VSITE  : 0);
	Flags = Flags | (bReproducible ? MD_REPRODUCIBLE : 0);
	Flags = Flags | (bAppendFiles  ? MD_APPENDFILES  : 0);
	Flags = Flags | (opt2parg_bSet("-append", asize(pa),pa) ? MD_APPENDFILESSET : 0); 
	Flags = Flags | (sim_part>1    ? MD_STARTFROMCPT : 0);
	Flags = Flags | (bResetCountersHalfWay ? MD_RESETCOUNTERSHALFWAY : 0);


	/* We postpone opening the log file if we are appending, so we can
   first truncate the old log file and append to the correct position
   there instead.  */
	if ((MASTER(cr) || bSepPot) && !bAppendFiles)
	{
		gmx_log_open(ftp2fn(efLOG,NFILE,fnm),cr,!bSepPot,Flags,&fplog);
		CopyRight(fplog,argv[0]);
		please_cite(fplog,"Hess2008b");
		please_cite(fplog,"Spoel2005a");
		please_cite(fplog,"Lindahl2001a");
		please_cite(fplog,"Berendsen95a");
	}
	else
	{
		fplog = NULL;
	}

	ddxyz[XX] = (int)(realddxyz[XX] + 0.5);
	ddxyz[YY] = (int)(realddxyz[YY] + 0.5);
	ddxyz[ZZ] = (int)(realddxyz[ZZ] + 0.5);

	/* even if nthreads = 1, we still call this one */

	rc = mdrunner_membed(fplog, cr, NFILE, fnm, oenv, bVerbose, bCompact,
			nstglobalcomm,
			ddxyz, dd_node_order, rdd, rconstr, dddlb_opt[0], dlb_scale,
			ddcsx, ddcsy, ddcsz, nstepout, resetstep, nmultisim, repl_ex_nst,
			repl_ex_seed, pforce, cpt_period, max_hours, deviceOptions, Flags,
			xy_fac,xy_max,z_fac,z_max,
			it_xy,it_z,probe_rad,low_up_rm,
			pieces,bALLOW_ASYMMETRY,maxwarn);

	if (gmx_parallel_env_initialized())
		gmx_finalize();

	if (MULTIMASTER(cr)) {
		thanx(stderr);
	}

	/* Log file has to be closed in mdrunner if we are appending to it
   (fplog not set here) */
	fprintf(stderr,"Please cite:\nWolf et al, J Comp Chem 31 (2010) 2169-2174.\n");

	if (MASTER(cr) && !bAppendFiles)
	{
		gmx_log_close(fplog);
	}

	return rc;
}
