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
 * Green Red Orange Magenta Azure Cyan Skyblue
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <typedefs.h>

#include "smalloc.h"
#include "macros.h"
#include "math.h"
#include "xvgr.h"
#include "copyrite.h"
#include "statutil.h"
#include "string2.h"
#include "vec.h"
#include "index.h"
#include "pbc.h"
#include "gmx_fatal.h"
#include "futil.h"
#include "gstat.h"
#include "gmx_ana.h"


static void add_contact_time(int **ccount,int *ccount_nalloc,int t)
{
  int i;

  if (t+2 >= *ccount_nalloc) {
    srenew(*ccount,t+2);
    for(i=*ccount_nalloc; i<t+2; i++)
      (*ccount)[i] = 0;
    *ccount_nalloc = t+2;
  }
  (*ccount)[t]++;
}

int gmx_dist(int argc,char *argv[])
{
  const char *desc[] = {
    "[TT]g_dist[tt] can calculate the distance between the centers of mass of two",
    "groups of atoms as a function of time. The total distance and its",
    "[IT]x[it]-, [IT]y[it]-, and [IT]z[it]-components are plotted.[PAR]",
    "Or when [TT]-dist[tt] is set, print all the atoms in group 2 that are",
    "closer than a certain distance to the center of mass of group 1.[PAR]",
    "With options [TT]-lt[tt] and [TT]-dist[tt] the number of contacts",
    "of all atoms in group 2 that are closer than a certain distance",
    "to the center of mass of group 1 are plotted as a function of the time",
    "that the contact was continuously present.[PAR]",
    "Other programs that calculate distances are [TT]g_mindist[tt]",
    "and [TT]g_bond[tt]."
  };
  
  t_topology *top=NULL;
  int  ePBC;
  real t,t0,cut2,dist2;
  rvec *x=NULL,*v=NULL,dx;
  matrix box;
  t_trxstatus *status;
  int natoms;

  int g,d,i,j,res,teller=0;
  atom_id aid;

  int     ngrps;     /* the number of index groups */
  atom_id **index, natoms_groups;   /* the index for the atom numbers */
  int     *isize;    /* the size of each group */
  char    **grpname; /* the name of each group */
  rvec    *com;
  real    *mass;
  FILE    *fp=NULL,*fplt=NULL;
  gmx_bool    bCutoff,bPrintDist,bLifeTime;
  t_pbc   *pbc;
  int     *contact_time=NULL,*ccount=NULL,ccount_nalloc=0,sum;
  char    buf[STRLEN];
  output_env_t oenv;
  gmx_rmpbc_t  gpbc=NULL;
  
  const char *leg[4] = { "|d|","d\\sx\\N","d\\sy\\N","d\\sz\\N" };

  static real cut=0;
  
  static t_pargs pa[] = {
    { "-dist",      FALSE, etREAL, {&cut},
      "Print all atoms in group 2 closer than dist to the center of mass of group 1" }
  };
#define NPA asize(pa)

  t_filenm fnm[] = {
    { efTRX, "-f", NULL, ffREAD },
    { efTPX, NULL, NULL, ffREAD },
    { efNDX, NULL, NULL, ffOPTRD },
    { efXVG, NULL, "dist", ffOPTWR },
    { efXVG, "-lt", "lifetime", ffOPTWR },
  };
#define NFILE asize(fnm)


  CopyRight(stderr,argv[0]);

  parse_common_args(&argc,argv,PCA_CAN_TIME | PCA_BE_NICE,
		    NFILE,fnm,NPA,pa,asize(desc),desc,0,NULL,&oenv);
  
  bCutoff = opt2parg_bSet("-dist",NPA,pa);
  cut2 = cut*cut;
  bLifeTime = opt2bSet("-lt",NFILE,fnm);
  bPrintDist = (bCutoff && !bLifeTime);
  
  top=read_top(ftp2fn(efTPX,NFILE,fnm),&ePBC);
  
  /* read index files */
  ngrps = 2;
  snew(com,ngrps);
  snew(grpname,ngrps);
  snew(index,ngrps);
  snew(isize,ngrps);
  get_index(&top->atoms,ftp2fn(efNDX,NFILE,fnm),ngrps,isize,index,grpname);
  
  /* calculate mass */
  natoms_groups=0;
  snew(mass,ngrps);
  for(g=0;(g<ngrps);g++) {
    mass[g]=0;
    for(i=0;(i<isize[g]);i++) {
      if (index[g][i]>natoms_groups)
	natoms_groups=index[g][i];
      if (index[g][i] >= top->atoms.nr)
	gmx_fatal(FARGS,"Atom number %d, item %d of group %d, is larger than number of atoms in the topolgy (%d)\n",index[g][i]+1,i+1,g+1,top->atoms.nr+1);
      mass[g]+=top->atoms.atom[index[g][i]].m;
    }
  }
  /* The number is one more than the highest index */
  natoms_groups++;

  natoms=read_first_x(oenv,&status,ftp2fn(efTRX,NFILE,fnm),&t,&x,box);
  t0 = t;

  if (natoms_groups>natoms)
    gmx_fatal(FARGS,"Atom number %d in an index group is larger than number of atoms in the trajectory (%d)\n",(int)natoms_groups,natoms);

  if (!bCutoff) {
    /* open output file */
    fp = xvgropen(ftp2fn(efXVG,NFILE,fnm),
		  "Distance","Time (ps)","Distance (nm)",oenv);
    xvgr_legend(fp,4,leg,oenv);
  } else {
    ngrps = 1;
    if (bLifeTime)
      snew(contact_time,isize[1]);
  }
  if (ePBC != epbcNONE)
    snew(pbc,1);
  else
    pbc = NULL;
    
  gpbc = gmx_rmpbc_init(&top->idef,ePBC,natoms,box);
  do {
    /* initialisation for correct distance calculations */
    if (pbc) {
      set_pbc(pbc,ePBC,box);
      /* make molecules whole again */
      gmx_rmpbc(gpbc,natoms,box,x);
    }
    /* calculate center of masses */
    for(g=0;(g<ngrps);g++) {
      if (isize[g] == 1) {
	copy_rvec(x[index[g][0]],com[g]);
      } else {
	for(d=0;(d<DIM);d++) {
	  com[g][d]=0;
	  for(i=0;(i<isize[g]);i++) {
	    com[g][d] += x[index[g][i]][d] * top->atoms.atom[index[g][i]].m;
	  }
	  com[g][d] /= mass[g];
	}
      }
    }
    
    if (!bCutoff) {
      /* write to output */
      fprintf(fp,"%12.7f ",t);
      for(g=0;(g<ngrps/2);g++) {
	if (pbc)
	  pbc_dx(pbc,com[2*g],com[2*g+1],dx);
	else
	  rvec_sub(com[2*g],com[2*g+1],dx);
	
	fprintf(fp,"%12.7f %12.7f %12.7f %12.7f",
		norm(dx),dx[XX],dx[YY],dx[ZZ]);
      }
      fprintf(fp,"\n");
    } else {
      for(i=0;(i<isize[1]);i++) { 
	j=index[1][i];
	if (pbc)
	  pbc_dx(pbc,x[j],com[0],dx);
	else
	  rvec_sub(x[j],com[0],dx);
	
	dist2 = norm2(dx);
	if (dist2<cut2) {
	  if (bPrintDist) {
	    res=top->atoms.atom[j].resind;
	    fprintf(stdout,"\rt: %g  %d %s %d %s  %g (nm)\n",
		    t,top->atoms.resinfo[res].nr,*top->atoms.resinfo[res].name,
		    j+1,*top->atoms.atomname[j],sqrt(dist2));
	  }
	  if (bLifeTime)
	    contact_time[i]++;
	} else {
	  if (bLifeTime) {
	    if (contact_time[i]) {
	      add_contact_time(&ccount,&ccount_nalloc,contact_time[i]-1);
	      contact_time[i] = 0;
	    }
	  }
	}
      }
    }
    
    teller++;
  } while (read_next_x(oenv,status,&t,natoms,x,box));
  gmx_rmpbc_done(gpbc);

  if (!bCutoff)
    ffclose(fp);

  close_trj(status);
  
  if (bCutoff && bLifeTime) {
    /* Add the contacts still present in the last frame */
    for(i=0; i<isize[1]; i++)
      if (contact_time[i])
	add_contact_time(&ccount,&ccount_nalloc,contact_time[i]-1);

    sprintf(buf,"%s - %s within %g nm",
	    grpname[0],grpname[1],cut);
    fp = xvgropen(opt2fn("-lt",NFILE,fnm),
		  buf,"Time (ps)","Number of contacts",oenv);
    for(i=0; i<min(ccount_nalloc,teller-1); i++) {
      /* Account for all subintervals of longer intervals */
      sum = 0;
      for(j=i; j<ccount_nalloc; j++)
	sum += (j-i+1)*ccount[j];

      fprintf(fp,"%10.3f %10.3f\n",i*(t-t0)/(teller-1),sum/(double)(teller-i));
    }
    ffclose(fp);
  }
  
  thanx(stderr);
  return 0;
}
