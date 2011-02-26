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
} rm_t;

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

int gen_rm_list(rm_t *rm_p,t_block *ins_at,t_block *rest_at,t_pbc *pbc, gmx_mtop_t *mtop,
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

void rm_group(t_inputrec *ir, gmx_groups_t *groups, gmx_mtop_t *mtop, rm_t *rm_p, t_state *state, t_block *ins_at, pos_ins_t *pos_ins)
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

void top_update(const char *topfile, char *ins, rm_t *rm_p, gmx_mtop_t *mtop)
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

void md_print_warning(const t_commrec *cr,FILE *fplog,const char *buf)
{
    if (MASTER(cr))
    {
        fprintf(stderr,"\n%s\n",buf);
    }
    if (fplog)
    {
        fprintf(fplog,"\n%s\n",buf);
    }
}

/*  simulation conditions to transmit. Keep in mind that they are
    transmitted to other nodes through an MPI_Reduce after
    casting them to a real (so the signals can be sent together with other
    data). This means that the only meaningful values are positive,
    negative or zero. */
enum { eglsNABNSB, eglsCHKPT, eglsSTOPCOND, eglsRESETCOUNTERS, eglsNR };
/* Is the signal in one simulation independent of other simulations? */
gmx_bool gs_simlocal[eglsNR] = { TRUE, FALSE, FALSE, TRUE };

typedef struct {
    int nstms;       /* The frequency for intersimulation communication */
    int sig[eglsNR]; /* The signal set by one process in do_md */
    int set[eglsNR]; /* The communicated signal, equal for all processes */
} globsig_t;


static int multisim_min(const gmx_multisim_t *ms,int nmin,int n)
{
    int  *buf;
    gmx_bool bPos,bEqual;
    int  s,d;

    snew(buf,ms->nsim);
    buf[ms->sim] = n;
    gmx_sumi_sim(ms->nsim,buf,ms);
    bPos   = TRUE;
    bEqual = TRUE;
    for(s=0; s<ms->nsim; s++)
    {
        bPos   = bPos   && (buf[s] > 0);
        bEqual = bEqual && (buf[s] == buf[0]);
    }
    if (bPos)
    {
        if (bEqual)
        {
            nmin = min(nmin,buf[0]);
        }
        else
        {
            /* Find the least common multiple */
            for(d=2; d<nmin; d++)
            {
                s = 0;
                while (s < ms->nsim && d % buf[s] == 0)
                {
                    s++;
                }
                if (s == ms->nsim)
                {
                    /* We found the LCM and it is less than nmin */
                    nmin = d;
                    break;
                }
            }
        }
    }
    sfree(buf);

    return nmin;
}

static int multisim_nstsimsync(const t_commrec *cr,
                               const t_inputrec *ir,int repl_ex_nst)
{
    int nmin;

    if (MASTER(cr))
    {
        nmin = INT_MAX;
        nmin = multisim_min(cr->ms,nmin,ir->nstlist);
        nmin = multisim_min(cr->ms,nmin,ir->nstcalcenergy);
        nmin = multisim_min(cr->ms,nmin,repl_ex_nst);
        if (nmin == INT_MAX)
        {
            gmx_fatal(FARGS,"Can not find an appropriate interval for inter-simulation communication, since nstlist, nstcalcenergy and -replex are all <= 0");
        }
        /* Avoid inter-simulation communication at every (second) step */
        if (nmin <= 2)
        {
            nmin = 10;
        }
    }

    gmx_bcast(sizeof(int),&nmin,cr);

    return nmin;
}

static void init_global_signals(globsig_t *gs,const t_commrec *cr,
                                const t_inputrec *ir,int repl_ex_nst)
{
    int i;

    if (MULTISIM(cr))
    {
        gs->nstms = multisim_nstsimsync(cr,ir,repl_ex_nst);
        if (debug)
        {
            fprintf(debug,"Syncing simulations for checkpointing and termination every %d steps\n",gs->nstms);
        }
    }
    else
    {
        gs->nstms = 1;
    }

    for(i=0; i<eglsNR; i++)
    {
        gs->sig[i] = 0;
        gs->set[i] = 0;
    }
}

static void copy_coupling_state(t_state *statea,t_state *stateb,
                                gmx_ekindata_t *ekinda,gmx_ekindata_t *ekindb, t_grpopts* opts)
{

    /* MRS note -- might be able to get rid of some of the arguments.  Look over it when it's all debugged */

    int i,j,nc;

    /* Make sure we have enough space for x and v */
    if (statea->nalloc > stateb->nalloc)
    {
        stateb->nalloc = statea->nalloc;
        srenew(stateb->x,stateb->nalloc);
        srenew(stateb->v,stateb->nalloc);
    }

    stateb->natoms     = statea->natoms;
    stateb->ngtc       = statea->ngtc;
    stateb->nnhpres    = statea->nnhpres;
    stateb->veta       = statea->veta;
    if (ekinda)
    {
        copy_mat(ekinda->ekin,ekindb->ekin);
        for (i=0; i<stateb->ngtc; i++)
        {
            ekindb->tcstat[i].T = ekinda->tcstat[i].T;
            ekindb->tcstat[i].Th = ekinda->tcstat[i].Th;
            copy_mat(ekinda->tcstat[i].ekinh,ekindb->tcstat[i].ekinh);
            copy_mat(ekinda->tcstat[i].ekinf,ekindb->tcstat[i].ekinf);
            ekindb->tcstat[i].ekinscalef_nhc =  ekinda->tcstat[i].ekinscalef_nhc;
            ekindb->tcstat[i].ekinscaleh_nhc =  ekinda->tcstat[i].ekinscaleh_nhc;
            ekindb->tcstat[i].vscale_nhc =  ekinda->tcstat[i].vscale_nhc;
        }
    }
    copy_rvecn(statea->x,stateb->x,0,stateb->natoms);
    copy_rvecn(statea->v,stateb->v,0,stateb->natoms);
    copy_mat(statea->box,stateb->box);
    copy_mat(statea->box_rel,stateb->box_rel);
    copy_mat(statea->boxv,stateb->boxv);

    for (i = 0; i<stateb->ngtc; i++)
    {
        nc = i*opts->nhchainlength;
        for (j=0; j<opts->nhchainlength; j++)
        {
            stateb->nosehoover_xi[nc+j]  = statea->nosehoover_xi[nc+j];
            stateb->nosehoover_vxi[nc+j] = statea->nosehoover_vxi[nc+j];
        }
    }
    if (stateb->nhpres_xi != NULL)
    {
        for (i = 0; i<stateb->nnhpres; i++)
        {
            nc = i*opts->nhchainlength;
            for (j=0; j<opts->nhchainlength; j++)
            {
                stateb->nhpres_xi[nc+j]  = statea->nhpres_xi[nc+j];
                stateb->nhpres_vxi[nc+j] = statea->nhpres_vxi[nc+j];
            }
        }
    }
}

static void compute_globals(FILE *fplog, gmx_global_stat_t gstat, t_commrec *cr, t_inputrec *ir,
                            t_forcerec *fr, gmx_ekindata_t *ekind,
                            t_state *state, t_state *state_global, t_mdatoms *mdatoms,
                            t_nrnb *nrnb, t_vcm *vcm, gmx_wallcycle_t wcycle,
                            gmx_enerdata_t *enerd,tensor force_vir, tensor shake_vir, tensor total_vir,
                            tensor pres, rvec mu_tot, gmx_constr_t constr,
                            globsig_t *gs,gmx_bool bInterSimGS,
                            matrix box, gmx_mtop_t *top_global, real *pcurr,
                            int natoms, gmx_bool *bSumEkinhOld, int flags)
{
  return;
}

static void check_nst_param(FILE *fplog,t_commrec *cr,
                            const char *desc_nst,int nst,
                            const char *desc_p,int *p)
{
    char buf[STRLEN];

    if (*p > 0 && *p % nst != 0)
    {
        /* Round up to the next multiple of nst */
        *p = ((*p)/nst + 1)*nst;
        sprintf(buf,"NOTE: %s changes %s to %d\n",desc_nst,desc_p,*p);
        md_print_warning(cr,fplog,buf);
    }
}

static void reset_all_counters(FILE *fplog,t_commrec *cr,
                               gmx_large_int_t step,
                               gmx_large_int_t *step_rel,t_inputrec *ir,
                               gmx_wallcycle_t wcycle,t_nrnb *nrnb,
                               gmx_runtime_t *runtime)
{
    char buf[STRLEN],sbuf[STEPSTRSIZE];

    /* Reset all the counters related to performance over the run */
    sprintf(buf,"Step %s: resetting all time and cycle counters\n",
            gmx_step_str(step,sbuf));
    md_print_warning(cr,fplog,buf);

    wallcycle_stop(wcycle,ewcRUN);
    wallcycle_reset_all(wcycle);
    if (DOMAINDECOMP(cr))
    {
        reset_dd_statistics_counters(cr->dd);
    }
    init_nrnb(nrnb);
    ir->init_step += *step_rel;
    ir->nsteps    -= *step_rel;
    *step_rel = 0;
    wallcycle_start(wcycle,ewcRUN);
    runtime_start(runtime);
    print_date_and_time(fplog,cr->nodeid,"Restarted time",runtime);
}

static int check_nstglobalcomm(FILE *fplog,t_commrec *cr,
                               int nstglobalcomm,t_inputrec *ir)
{
    char buf[STRLEN];

    if (!EI_DYNAMICS(ir->eI))
    {
        nstglobalcomm = 1;
    }

    if (nstglobalcomm == -1)
    {
        if (ir->nstcalcenergy == 0 && ir->nstlist == 0)
        {
            nstglobalcomm = 10;
            if (ir->nstenergy > 0 && ir->nstenergy < nstglobalcomm)
            {
                nstglobalcomm = ir->nstenergy;
            }
        }
        else
        {
            /* We assume that if nstcalcenergy > nstlist,
             * nstcalcenergy is a multiple of nstlist.
             */
            if (ir->nstcalcenergy == 0 ||
                (ir->nstlist > 0 && ir->nstlist < ir->nstcalcenergy))
            {
                nstglobalcomm = ir->nstlist;
            }
            else
            {
                nstglobalcomm = ir->nstcalcenergy;
            }
        }
    }
    else
    {
        if (ir->nstlist > 0 &&
            nstglobalcomm > ir->nstlist && nstglobalcomm % ir->nstlist != 0)
        {
            nstglobalcomm = (nstglobalcomm / ir->nstlist)*ir->nstlist;
            sprintf(buf,"WARNING: nstglobalcomm is larger than nstlist, but not a multiple, setting it to %d\n",nstglobalcomm);
            md_print_warning(cr,fplog,buf);
        }
        if (nstglobalcomm > ir->nstcalcenergy)
        {
            check_nst_param(fplog,cr,"-gcom",nstglobalcomm,
                            "nstcalcenergy",&ir->nstcalcenergy);
        }

        check_nst_param(fplog,cr,"-gcom",nstglobalcomm,
                        "nstenergy",&ir->nstenergy);

        check_nst_param(fplog,cr,"-gcom",nstglobalcomm,
                        "nstlog",&ir->nstlog);
    }

    if (ir->comm_mode != ecmNO && ir->nstcomm < nstglobalcomm)
    {
        sprintf(buf,"WARNING: Changing nstcomm from %d to %d\n",
                ir->nstcomm,nstglobalcomm);
        md_print_warning(cr,fplog,buf);
        ir->nstcomm = nstglobalcomm;
    }

    return nstglobalcomm;
}

void check_ir_old_tpx_versions(t_commrec *cr,FILE *fplog,
                               t_inputrec *ir,gmx_mtop_t *mtop)
{
    /* Check required for old tpx files */
    if (IR_TWINRANGE(*ir) && ir->nstlist > 1 &&
        ir->nstcalcenergy % ir->nstlist != 0)
    {
        md_print_warning(cr,fplog,"Old tpr file with twin-range settings: modifying energy calculation and/or T/P-coupling frequencies");

        if (gmx_mtop_ftype_count(mtop,F_CONSTR) +
            gmx_mtop_ftype_count(mtop,F_CONSTRNC) > 0 &&
            ir->eConstrAlg == econtSHAKE)
        {
            md_print_warning(cr,fplog,"With twin-range cut-off's and SHAKE the virial and pressure are incorrect");
            if (ir->epc != epcNO)
            {
                gmx_fatal(FARGS,"Can not do pressure coupling with twin-range cut-off's and SHAKE");
            }
        }
        check_nst_param(fplog,cr,"nstlist",ir->nstlist,
                        "nstcalcenergy",&ir->nstcalcenergy);
    	check_nst_param(fplog,cr,"nstcalcenergy",ir->nstcalcenergy,
        	        "nstenergy",&ir->nstenergy);
        check_nst_param(fplog,cr,"nstcalcenergy",ir->nstcalcenergy,
                        "nstlog",&ir->nstlog);
        if (ir->efep != efepNO)
        {
            check_nst_param(fplog,cr,"nstcalcenergy",ir->nstcalcenergy,
                            "nstdhdl",&ir->nstdhdl);
        }
    }
}

typedef struct {
    gmx_bool       bGStatEveryStep;
    gmx_large_int_t step_ns;
    gmx_large_int_t step_nscheck;
    gmx_large_int_t nns;
    matrix     scale_tot;
    int        nabnsb;
    double     s1;
    double     s2;
    double     ab;
    double     lt_runav;
    double     lt_runav2;
} gmx_nlheur_t;

static void reset_nlistheuristics(gmx_nlheur_t *nlh,gmx_large_int_t step)
{
    nlh->lt_runav  = 0;
    nlh->lt_runav2 = 0;
    nlh->step_nscheck = step;
}

static void init_nlistheuristics(gmx_nlheur_t *nlh,
                                 gmx_bool bGStatEveryStep,gmx_large_int_t step)
{
    nlh->bGStatEveryStep = bGStatEveryStep;
    nlh->nns       = 0;
    nlh->nabnsb    = 0;
    nlh->s1        = 0;
    nlh->s2        = 0;
    nlh->ab        = 0;

    reset_nlistheuristics(nlh,step);
}

static void update_nliststatistics(gmx_nlheur_t *nlh,gmx_large_int_t step)
{
    gmx_large_int_t nl_lt;
    char sbuf[STEPSTRSIZE],sbuf2[STEPSTRSIZE];

    /* Determine the neighbor list life time */
    nl_lt = step - nlh->step_ns;
    if (debug)
    {
        fprintf(debug,"%d atoms beyond ns buffer, updating neighbor list after %s steps\n",nlh->nabnsb,gmx_step_str(nl_lt,sbuf));
    }
    nlh->nns++;
    nlh->s1 += nl_lt;
    nlh->s2 += nl_lt*nl_lt;
    nlh->ab += nlh->nabnsb;
    if (nlh->lt_runav == 0)
    {
        nlh->lt_runav  = nl_lt;
        /* Initialize the fluctuation average
         * such that at startup we check after 0 steps.
         */
        nlh->lt_runav2 = sqr(nl_lt/2.0);
    }
    /* Running average with 0.9 gives an exp. history of 9.5 */
    nlh->lt_runav2 = 0.9*nlh->lt_runav2 + 0.1*sqr(nlh->lt_runav - nl_lt);
    nlh->lt_runav  = 0.9*nlh->lt_runav  + 0.1*nl_lt;
    if (nlh->bGStatEveryStep)
    {
        /* Always check the nlist validity */
        nlh->step_nscheck = step;
    }
    else
    {
        /* We check after:  <life time> - 2*sigma
         * The factor 2 is quite conservative,
         * but we assume that with nstlist=-1 the user
         * prefers exact integration over performance.
         */
        nlh->step_nscheck = step
                  + (int)(nlh->lt_runav - 2.0*sqrt(nlh->lt_runav2)) - 1;
    }
    if (debug)
    {
        fprintf(debug,"nlist life time %s run av. %4.1f sig %3.1f check %s check with -gcom %d\n",
                gmx_step_str(nl_lt,sbuf),nlh->lt_runav,sqrt(nlh->lt_runav2),
                gmx_step_str(nlh->step_nscheck-step+1,sbuf2),
                (int)(nlh->lt_runav - 2.0*sqrt(nlh->lt_runav2)));
    }
}

static void set_nlistheuristics(gmx_nlheur_t *nlh,gmx_bool bReset,gmx_large_int_t step)
{
    int d;

    if (bReset)
    {
        reset_nlistheuristics(nlh,step);
    }
    else
    {
        update_nliststatistics(nlh,step);
    }

    nlh->step_ns = step;
    /* Initialize the cumulative coordinate scaling matrix */
    clear_mat(nlh->scale_tot);
    for(d=0; d<DIM; d++)
    {
        nlh->scale_tot[d][d] = 1.0;
    }
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
  return 0;
}

int gmx_membed(int argc,char *argv[])
{
	const char *desc[] = {
			"[TT]g_membed[tt] embeds a membrane protein into an equilibrated lipid bilayer at the position",
			"and orientation specified by the user.\n",
			"\n",
			"SHORT MANUAL\n------------\n",
			"The user should merge the structure files of the protein and membrane (+solvent), creating a",
			"single structure file with the protein overlapping the membrane at the desired position and",
			"orientation. Box size should be taken from the membrane structure file. The corresponding topology",
			"files should also be merged. Consecutively, create a [TT].tpr[tt] file (input for [TT]g_membed[tt]) from these files,"
			"with the following options included in the [TT].mdp[tt] file.\n",
			" - [TT]integrator      = md[tt][BR]",
			" - [TT]energygrp       = Protein[tt] (or other group that you want to insert)[BR]",
			" - [TT]freezegrps      = Protein[tt][BR]",
			" - [TT]freezedim       = Y Y Y[tt][BR]",
			" - [TT]energygrp_excl  = Protein Protein[tt][BR]",
			"The output is a structure file containing the protein embedded in the membrane. If a topology",
			"file is provided, the number of lipid and ",
			"solvent molecules will be updated to match the new structure file.\n",
			"For a more extensive manual see Wolf et al, J Comp Chem 31 (2010) 2169-2174, Appendix.\n",
			"\n",
			"SHORT METHOD DESCRIPTION\n",
			"------------------------\n",
			"1. The protein is resized around its center of mass by a factor -xy in the xy-plane",
			"(the membrane plane) and a factor -z in the z-direction (if the size of the",
			"protein in the z-direction is the same or smaller than the width of the membrane, a",
			"-z value larger than 1 can prevent that the protein will be enveloped by the lipids).\n",
			"2. All lipid and solvent molecules overlapping with the resized protein are removed. All",
			"intraprotein interactions are turned off to prevent numerical issues for small values of -xy",
			" or -z\n",
			"3. One md step is performed.\n",
			"4. The resize factor (-xy or -z) is incremented by a small amount ((1-xy)/nxy or (1-z)/nz) and the",
			"protein is resized again around its center of mass. The resize factor for the xy-plane",
			"is incremented first. The resize factor for the z-direction is not changed until the -xy factor",
			"is 1 (thus after -nxy iteration).\n",
			"5. Repeat step 3 and 4 until the protein reaches its original size (-nxy + -nz iterations).\n",
			"For a more extensive method descrition see Wolf et al, J Comp Chem, 31 (2010) 2169-2174.\n",
			"\n",
			"NOTE\n----\n",
			" - Protein can be any molecule you want to insert in the membrane.\n",
			" - It is recommended to perform a short equilibration run after the embedding",
			"(see Wolf et al, J Comp Chem 31 (2010) 2169-2174, to re-equilibrate the membrane. Clearly",
			"protein equilibration might require longer.\n",
			"\n"
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
			{ efNDX, "-dn",     "dipole",   ffOPTWR }
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
    "HIDDENAttempt replica exchange every # steps" },
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
		nthreads=tMPI_Get_recommended_nthreads();
	}
#else
	nthreads=1;
#endif


	if (repl_ex_nst != 0 && nmultisim < 2)
		gmx_fatal(FARGS,"Need at least two replicas for replica exchange (option -multi)");

	if (nmultisim > 1) {
#ifndef GMX_THREADS
		init_multisystem(cr,nmultisim,NFILE,fnm,TRUE);
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
	return 0;
}
