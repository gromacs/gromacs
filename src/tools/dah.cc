/*
 *       @(#) copyrgt.c 1.12 9/30/97
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0b
 * 
 * Copyright (c) 1990-1997,
 * BIOSON Research Institute, Dept. of Biophysical Chemistry,
 * University of Groningen, The Netherlands
 *
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 *
 * Also check out our WWW page:
 * http://rugmd0.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 *
 * And Hey:
 * Good gRace! Old Maple Actually Chews Slate
 */
#include <maths.h>
#include "dah.h"
#include "typedefs.h"
#include "sysstuff.h"
#include "parse.h"
#include "rdgroup.h"
#include "hbond.h"
#include "fatal.h"
#include "vec.h"
#include "pbc.h"
#include "xvgr.h" 
#include "copyrite.h"

int in_list(atom_id selection,int isize,atom_id *index)
{
  int i=0;
  int return_value=1;
  for(i=0;(i<isize);i++)
    if(selection==index[i])
      return_value=0;
  return(return_value);
}

void search_acceptors(t_topology *top,
		      int *nr_a, atom_id **a,
                      bool nitro_acc,
                      int isize,atom_id *index)
{
  int i;
  for (i=0;(i<top->atoms.nr);i++)
    if (*top->atoms.atomname[i][0] == 'O') {
      if (in_list(i,isize,index)==0) {
        (*nr_a)++;
        (*a)[*nr_a-1]=i;
      }
    }
  if (nitro_acc==TRUE)
    for (i=0;(i<top->atoms.nr);i++)
      if (*top->atoms.atomname[i][0] == 'N') {
        if (in_list(i,isize,index)==0) {
          (*nr_a)++;
          (*a)[*nr_a-1]=i;
        }
      }
}

void search_donors(t_topology *top,t_ilist *interaction,int func_type,
		   int *nr_d,atom_id **d,atom_id **h,int isize,atom_id *index)
{
  int i;
  t_functype type;
  atom_id nr_1,nr_2;
  
  for(i=0;(i<interaction->nr);) {
    /* print some values */
    type=top->idef.functype[interaction->iatoms[i]];

    /* check out this functype */
    if (func_type == F_SETTLE) {
      nr_1=interaction->iatoms[i+1];
      
      if (in_list(nr_1,isize,index)==0) {
	(*nr_d)++;
	(*d)[*nr_d-1]=nr_1;
	(*h)[*nr_d-1]=nr_1+1;
	(*nr_d)++;
	(*d)[*nr_d-1]=nr_1;
	(*h)[*nr_d-1]=nr_1+2;
      }
    }
    else if (type==func_type) {      
      nr_1=interaction->iatoms[i+1];
      nr_2=interaction->iatoms[i+2];
      
      if ((*top->atoms.atomname[nr_1][0])=='H') {
        if (((*top->atoms.atomname[nr_2][0])=='O') ||
            ((*top->atoms.atomname[nr_2][0])=='N')) {
          if (in_list(nr_2,isize,index)==0) {
            (*nr_d)++;
            (*d)[*nr_d-1]=nr_2;
            (*h)[*nr_d-1]=nr_1;
          }
        }
      }
      else
        if ((*top->atoms.atomname[nr_2][0])=='H')
          if (((*top->atoms.atomname[nr_1][0])=='O') ||
              ((*top->atoms.atomname[nr_1][0])=='N')) {
            if (in_list(nr_1,isize,index)==0) {
              (*nr_d)++;
              (*d)[*nr_d-1]=nr_1;
              (*h)[*nr_d-1]=nr_2;
            }
          }
    } 

    /* next functype */
    i+=interaction_function[top->idef.functype[interaction->iatoms[i]]].nratoms+1;
  }
}

void gen_list(int& nr_dah,Hbond **dah,int nr_d,atom_id *d,atom_id *h,int nr_a,atom_id *a)
{
  int i,j;

  for(i=0;(i<nr_d);i++)
    for(j=0;(j<nr_a);j++)
      if (d[i]!=a[j]) {
	dah[nr_dah]=new Hbond(d[i],a[j],h[i]);
	nr_dah++;
      }
}

void gen_dah_list(int *isize,atom_id **index,int nr_groups,
                  Hbond ****dah,int& nr_dah,
		  char *trjfile)
{
  int  status;
  FILE *fp;
  int nr_d=0,nr_a=0,*nr_dd,*nr_aa,i,j,total=0;
  atom_id *d,*h,*a;
  atom_id **dd,**hh,**aa;

  switch (mode) {
  case INTER_GRP :
    dd  = (atom_id **)malloc(nr_groups*sizeof(atom_id *));
    hh  = (atom_id **)malloc(nr_groups*sizeof(atom_id *));
    aa  = (atom_id **)malloc(nr_groups*sizeof(atom_id *));
    nr_dd  = (int *)malloc(nr_groups*sizeof(int));
    nr_aa  = (int *)malloc(nr_groups*sizeof(int));
    for(i=0;(i<nr_groups);i++) {
      dd[i] = (atom_id *)malloc(top->atoms.nr*sizeof(atom_id));
      aa[i] = (atom_id *)malloc(top->atoms.nr*sizeof(atom_id));
      hh[i] = (atom_id *)malloc(top->atoms.nr*sizeof(atom_id));
      nr_aa[i]=0;
      nr_dd[i]=0;
      search_donors(top,&top->idef.il[F_SHAKE],
		    F_SHAKE,&nr_dd[i],&dd[i],&hh[i],isize[i],index[i]);
      search_donors(top,&top->idef.il[F_BONDS],
		    F_BONDS,&nr_dd[i],&dd[i],&hh[i],isize[i],index[i]);
      search_donors(top,&top->idef.il[F_SETTLE],
		    F_SETTLE,&nr_dd[i],&dd[i],&hh[i],isize[i],index[i]);
      search_acceptors(top,&nr_aa[i],&aa[i],FALSE,isize[i],index[i]);
      fprintf(stderr,"group %2d : number of donors = %5d  number of acceptors = %5d \n",
	      i,nr_dd[i],nr_aa[i]);
    }
    for(i=0;(i<nr_groups);i++)
      for(j=0;(j<nr_groups);j++)
	if (j!=i)
	  total+=nr_dd[i]*nr_aa[j];
    **dah = (Hbond **)malloc(total*sizeof(Hbond *));
    nr_dah=0;
    for(i=0;(i<nr_groups);i++)
      for(j=0;(j<nr_groups);j++)
	if (j!=i)
	  gen_list(nr_dah,**dah,nr_dd[i],dd[i],hh[i],nr_aa[j],aa[j]);
    for(i=0;(i<nr_groups);i++) {
      free(dd[i]);
      free(aa[i]);
      free(hh[i]);
    }
    free(dd);
    free(hh);
    free(aa);
    free(nr_dd);
    free(nr_aa);
    break;
  case INTRA_GRP :
    d = (atom_id *)malloc(top->atoms.nr*sizeof(atom_id));
    a = (atom_id *)malloc(top->atoms.nr*sizeof(atom_id));
    h = (atom_id *)malloc(top->atoms.nr*sizeof(atom_id));
    search_donors(top,&top->idef.il[F_SHAKE],
		  F_SHAKE,&nr_d,&d,&h,isize[0],index[0]);
    fprintf(stderr,"%5d donors \n",nr_d);
    search_donors(top,&top->idef.il[F_BONDS],
		  F_BONDS,&nr_d,&d,&h,isize[0],index[0]);
    fprintf(stderr,"%5d donors \n",nr_d);
    search_donors(top,&top->idef.il[F_SETTLE],
		  F_SETTLE,&nr_d,&d,&h,isize[0],index[0]);
    fprintf(stderr,"%5d donors \n",nr_d);
    search_acceptors(top,&nr_a,&a,FALSE,isize[0],index[0]);
    fprintf(stderr,"%5d acceptors \n",nr_a);
    nr_dah=0;
    **dah = (Hbond **)malloc(nr_d*nr_a*sizeof(Hbond *));
    gen_list(nr_dah,**dah,nr_d,d,h,nr_a,a);
    free(d);
    free(a);
    free(h);
    break;
  case SELECTED :
  case INSERT   :
    d = (atom_id *)malloc(top->atoms.nr*sizeof(atom_id));
    a = (atom_id *)malloc(top->atoms.nr*sizeof(atom_id));
    h = (atom_id *)malloc(top->atoms.nr*sizeof(atom_id));
    search_donors(top,&top->idef.il[F_SHAKE],
		  F_SHAKE,&nr_d,&d,&h,isize[0],index[0]);
    search_donors(top,&top->idef.il[F_BONDS],
		  F_BONDS,&nr_d,&d,&h,isize[0],index[0]);
    search_donors(top,&top->idef.il[F_SETTLE],
		  F_SETTLE,&nr_d,&d,&h,isize[0],index[0]);
    search_acceptors(top,&nr_a,&a,FALSE,isize[0],index[0]);
    nr_dah=0;
    **dah = (Hbond **)malloc(isize[0]/3 * sizeof(Hbond *));
    for(i=0;(i<isize[0]);) {
      if (((in_list(index[0][i],nr_d,d)==0) &&
           (in_list(index[0][i+1],nr_d,h)==0) &&
           (in_list(index[0][i+2],nr_a,a)==0))) {
        (**dah)[nr_dah]=new Hbond(index[0][i],index[0][i+2],index[0][i+1]);
        nr_dah++;
        i+=3;
      }
      else {
        fprintf(stderr,"Ilegal format in selected group\n");
        fprintf(stderr,"The right format must be :");
        fprintf(stderr,"donor hydrogen acceptor donor hydrogen acceptor ...\n");
        fprintf(stderr,"program terminated\n");
        exit(0);
      }
    }
    free(d);
    free(a);
    free(h);
    break;
  case COUNT:
    

    dd  = (atom_id **)malloc(nr_groups*sizeof(atom_id *));
    hh  = (atom_id **)malloc(nr_groups*sizeof(atom_id *));
    aa  = (atom_id **)malloc(nr_groups*sizeof(atom_id *));
    nr_dd  = (int *)malloc(nr_groups*sizeof(int));
    nr_aa  = (int *)malloc(nr_groups*sizeof(int));
    for(i=0;(i<nr_groups);i++) {
      dd[i] = (atom_id *)malloc(top->atoms.nr*sizeof(atom_id));
      aa[i] = (atom_id *)malloc(top->atoms.nr*sizeof(atom_id));
      hh[i] = (atom_id *)malloc(top->atoms.nr*sizeof(atom_id));
      nr_aa[i]=0;
      nr_dd[i]=0;
      search_donors(top,&top->idef.il[F_SHAKE],
		    F_SHAKE,&nr_dd[i],&dd[i],&hh[i],isize[i],index[i]);
      search_donors(top,&top->idef.il[F_BONDS],
		    F_BONDS,&nr_dd[i],&dd[i],&hh[i],isize[i],index[i]);
      search_donors(top,&top->idef.il[F_SETTLE],
		    F_SETTLE,&nr_dd[i],&dd[i],&hh[i],isize[i],index[i]);
      search_acceptors(top,&nr_aa[i],&aa[i],FALSE,isize[i],index[i]);
      fprintf(stderr,"group %2d : number of donors = %5d  number of acceptors = %5d \n",i,nr_dd[i],nr_aa[i]);
    }

    fp = xvgropen("number_inter.xvg","Intermolecular Hydrogen Bonds","Time (ps)","Number");



    read_first_x(&status,trjfile,&this_time,&x,box);
    do {
      nr_hbonds=0;

      /* search for hydrogen bonds */
      /* donors of group 0 and acceptors of group 1 */
      for(i=0;(i<nr_dd[0]);i++)
	for(j=0;(j<nr_aa[1]);j++) {
	  static rvec dx;
	  pbc_dx(box,x[dd[0][i]],x[aa[1][j]],dx);
	  if ( norm(dx) < rcut ) {
	    static rvec v_1,v_2;
	    pbc_dx(box,x[dd[0][i]],x[hh[0][i]],v_1);
	    pbc_dx(box,x[dd[0][i]],x[aa[1][j]],v_2);
	    if (acos(cos_angle(v_1,v_2)) < alfcut )
	      nr_hbonds++;
	  }
	}

      /* donors of group 1 and acceptors of group 0 */
      for(i=0;(i<nr_dd[1]);i++)
	for(j=0;(j<nr_aa[0]);j++) {
	  static rvec dx;
	  pbc_dx(box,x[dd[1][i]],x[aa[0][j]],dx);
	  if ( norm(dx) < rcut ) {
	    static rvec v_1,v_2;
	    pbc_dx(box,x[dd[1][i]],x[hh[1][i]],v_1);
	    pbc_dx(box,x[dd[1][i]],x[aa[0][j]],v_2);
	    if (acos(cos_angle(v_1,v_2)) < alfcut )
	      nr_hbonds++;
	  }
	}
      
      fprintf(stderr,"%6.3f %5d                         \r",this_time,nr_hbonds);
      fprintf(fp,"%10g  %5d\n",this_time,nr_hbonds);
      fflush(fp);
    } while (read_next_x(status,&this_time,natoms,x,box));           

    close_trj(status);
    
    fflush(fp);
    fclose(fp);


    for(i=0;(i<nr_groups);i++) {
      free(dd[i]);
      free(aa[i]);
      free(hh[i]);
    }
    free(dd);
    free(hh);
    free(aa);
    free(nr_dd);
    free(nr_aa);

    thanx(stdout);
    exit(0);
    break;
  default :
    exit(0);
  }
}

static int compare_dah(void *s1,void *s2)
{
  Hbond *dah_1,*dah_2;
  dah_1 = (Hbond *)s1;
  dah_2 = (Hbond *)s2;
  
  return dah_1->compare(dah_2);

}


void swap(Hbond **dah,int i,int j)
{
  Hbond *temp;
 
  temp   = dah[i];
  dah[i] = dah[j];
  dah[j] = temp;
}

void quicksort(Hbond **dah, int left,int right,int (*comp)(void *,void *))
{
  int i,last;
  if (left >= right)
    return;
  swap(dah,left,(left+right)/2);
  last = left;
  for(i = left+1;i<=right;i++)
    if ((*comp)(dah[i],dah[left])<0)
      swap(dah,++last,i);
  swap(dah,left,last);
  quicksort(dah,left,last-1,comp);
  quicksort(dah,last+1,right,comp);
}

void check_dah_list(Hbond **dah,int& nr_dah)
{
  Hbond **dah_old;
  int number,i;
  number=nr_dah;
  quicksort(dah,0,number-1,compare_dah);

  /* copy dah list to dah_old list */
  dah_old = (Hbond **)malloc(nr_dah * sizeof(Hbond *));
  for(i=0;(i<nr_dah);i++) {
    dah_old[i]= new Hbond(dah[i]);
  }

  /* delete dah list */
  for(i=0;(i<nr_dah);i++)
    delete dah[i];

  /* generate new list */
  number=1;
  dah[0] = new Hbond(dah_old[0]);

  for(i=1;(i<nr_dah);i++) {
    if (dah_old[i]->compare(dah_old[i-1])!=0) {
      dah[number]=new Hbond(dah_old[i]);
      number++;
    }
  }

  /* delete dah old list */
  for(i=0;(i<nr_dah);i++)
    delete dah_old[i];
  free(dah_old);
  nr_dah=number;

}


void init_dah(Hbond ***dah,int& nr_dah,char *ndxfile,char *trjfile)
{
  int           nr_groups=1;
  atom_id       **index;
  int           *isize;
  char          **grpnames;


  user_input(&nr_groups);

  grpnames=(char **)malloc(nr_groups*sizeof(char *));
  index=(atom_id **)malloc(nr_groups*sizeof(atom_id *));
  isize=(int *)malloc(nr_groups*sizeof(int));
  
  rd_index(ndxfile,nr_groups,isize,index,grpnames);

  gen_dah_list(isize,index,nr_groups,&dah,nr_dah,trjfile);

  if (nr_dah==0) {
    fprintf(stderr,"I can't find any hydrogen bonds in the selected ");
    if (nr_groups>1)
      fprintf(stderr,"groups\n");
    else
      fprintf(stderr,"group\n");
    fprintf(stderr,"Check your %s file\n",ndxfile);
    fprintf(stderr,"program terminated...\n");
    exit(0);
  }

  if ((mode!=SELECTED)&&(mode!=INSERT))
    check_dah_list(*dah,nr_dah);

  if (mode==INSERT) 
    water_list_init(index[1],isize[1]);

}













