/*
 *       $Id$
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 1.6
 * 
 * Copyright (c) 1991-1997
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
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
 * GRoups of Organic Molecules in ACtion for Science
 */

#ifndef _list_cc
#define _list_cc

static char *SRCID_list_cc = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) list.cc 1.23 9/30/97"
#endif /* HAVE_IDENT */
#include "list.h"
#include "dah.h"
#include "string2.h"
#include "sysstuff.h"
#include <physics.h>
#include "futil.h"
#include "xvgr.h"
#include <pbc.h>
#include <vec.h>




#define RESOLUTION 200
#define MIN_ENER  -10.0
#define MAX_ENER  100.0

List::List()
{

  hbond = (Hbond **)malloc(sizeof(Hbond *));
  nr_hbonds = 0;
  t      = NULL;


  stat     = NULL;
  total    = NULL;

  /* number of hydrogen bonds */
  number   = NULL;

  /* distance and angle distribution */
  distance = NULL;
  angle    = NULL;

#ifdef ENERGY
  energy   = NULL;
#endif

  helix    = NULL;

  /* matrix */
  matr     = NULL;
}

List::~List()
{
  int i;


  for(i=0;(i<nr_hbonds);i++)
    delete(hbond[i]);
  free(hbond);
  free(t);

  for(i=0;(i<=eNNr);i++)
    free(helix[i]);
  free(helix);

  for(i=0;(i<typeNR);i++) {
    free(number[i]);
    free(distance[i]);
#ifdef ENERGY
    free(energy[i]);
#endif
    free(angle[i]);
  }
  free(number);
#ifdef ENERGY
  free(energy);
#endif
  free(distance);
  free(angle);
  free(stat);
  free(total);


  if (matr!=NULL) {
    for(i=0;(i<nr_hbonds);i++)
      free(matr[i]);
    free(matr);
  }
  
}

void List::checktype()
{
  int i,intra=0,inter=0;
  /* check type assignment */
  for(i=0;(i<nr_hbonds);i++) {
    if ( hbond[i]->type==INTRA )
      intra++;
    else if ( hbond[i]->type==INTER )
      inter++;
    else {
      fprintf(stderr,"error in type assignment\n");
      exit(0);
    }
  }
  fprintf(stderr,"INTRA : %5d\nINTER : %5d\n",intra,inter);
}

void List::print(FILE *output)
{
  int i;
  for(i=0;(i<nr_hbonds);i++) {
    hbond[i]->print(output);
    fprintf(output,"\n");
  }
}

int List::search(Hbond **dah,int& nr_dah)
{
  int i;
  nr_frames++;
  for(i=0;(i<nr_dah);i++) {
    if (dah[i]->exist()) {
      hbond = (Hbond **)realloc(hbond,++nr_hbonds*sizeof(Hbond *));
      hbond[nr_hbonds-1]=dah[i];
      nr_dah--;
      dah[i]=dah[nr_dah];
    }
  } 
  return nr_hbonds;
}

int List::count(int *nr_dd,atom_id **dd,atom_id **hh,int *nr_aa,atom_id **aa)
{
  int i,j;


  int nr_bonds=0;

  /* hbonds between donors of group 0 and acceptors of group 1 */
  for(i=0;(i<nr_dd[0]);i++) {
    for(j=0;(j<nr_aa[1]);j++) {
      /* calculate distance and angle */
      static rvec dx;
      pbc_dx(box,x[dd[0][i]],x[aa[1][j]],dx);
      if ( norm(dx) < rcut ) {
	static rvec v_1,v_2;
	pbc_dx(box,x[dd[0][i]],x[hh[0][i]],v_1);
	pbc_dx(box,x[dd[0][i]],x[aa[1][j]],v_2);
	if (acos(cos_angle(v_1,v_2)) < alfcut )
	  nr_bonds++;
      }
    }
  }
  
  /* hbonds between donors of group 1 and acceptors of group 0 */
  for(i=0;(i<nr_dd[1]);i++) {
    for(j=0;(j<nr_aa[0]);j++) {
      /* calculate distance and angle */
      static rvec dx;
      pbc_dx(box,x[dd[1][i]],x[aa[0][j]],dx);
      if ( norm(dx) < rcut ) {
	static rvec v_1,v_2;
	pbc_dx(box,x[dd[1][i]],x[hh[1][i]],v_1);
	pbc_dx(box,x[dd[1][i]],x[aa[0][j]],v_2);
	if (acos(cos_angle(v_1,v_2)) < alfcut )
	  nr_bonds++;
      }
    }
  } 
  return nr_bonds;
}

void List::nosearch(Hbond **dah,int& nr_dah)
{
  int i;
  
  hbond = (Hbond **)malloc(nr_dah*sizeof(Hbond *));
  
  for(i=0;(i<nr_dah);i++) 
    hbond[i]=dah[i];
  
  nr_hbonds=nr_dah;
  
  fprintf(stderr,"How many coord frames are in your trj file: ");
  do {
    scanf("%d",&nr_frames);
  } while (nr_frames < 0);
}


void List::hydrogenbond_ndx()
{
  int i,sum;

  FILE *fp;

  for(i=0,sum=0;(i<typeNR);i++)
    sum+=(total[i]*3);
  
  fp=ffopen("hydrogen_bonds.ndx","w");
  fprintf(fp,"%5d  %5d\n",typeNR,sum);
  
  /* internal hydrogen bonds */
  fprintf(fp,"internal  %5d \n",total[INTRA]*3);
  for(i=0;(i<nr_hbonds);i++)
    hbond[i]->ndx(INTRA,fp);
  
  /* intermolecular hydrogen bonds */
  fprintf(fp,"intermolecular  %5d \n",total[INTER]*3);
  for(i=0;(i<nr_hbonds);i++)
    hbond[i]->ndx(INTER,fp);

  /* all hydrogen bonds */
  fprintf(fp,"all_bonds  %5d \n",total[typeNR-1]*3);
  for(i=0;(i<nr_hbonds);i++)
    hbond[i]->ndx(ALL,fp);

  /* close file */
  fflush(fp);
  fclose(fp);
  
}

void List::assign_type()
{
  int n,i;
  atom_id left,right;

  /* resetting of type array has already been done
   * in the initialisation of hbond
   */

  /* Select intra bridges */
  for(n=0;(n<top->blocks[ebMOLS].nr);n++) {
    left =top->blocks[ebMOLS].index[n];
    right=top->blocks[ebMOLS].index[n+1];
    for(i=0;(i<nr_hbonds);i++)
      hbond[i]->settype(left,right);
  }
}

void List::analyse_init()
{
  int i,j;

  assign_type();

  total = (int *)malloc(typeNR*sizeof(int));

  for(i=0;(i<typeNR);i++) 
    total[i]=0;
  for(i=0;(i<nr_hbonds);i++) {
    total[hbond[i]->type]++;
    total[ALL]++;
  }

  hydrogenbond_ndx();

  stat = (int *)malloc(nr_hbonds * sizeof(int));
  for(i=0;(i<nr_hbonds);i++) {
    stat[i]=0;
    hbond[i]->analyse_init();
  }

  switch (mode) {
  case INSERT:
  case SELECTED:
    matr = (char **)malloc(nr_hbonds * sizeof(char *));
    for(i=0;(i<nr_hbonds);i++)
      matr[i]=(char *)malloc(nr_frames *sizeof(char));
    break;
  case INTRA_GRP:
  case INTER_GRP:
    break;
  default:
    fprintf(stderr,"Illegal choice in list.cc\n");
    exit(0);
  }
  
  /* number angle and type graphs */
  t = (real *)malloc(nr_frames*sizeof(real)); 
  number   = (int  **)malloc(typeNR*sizeof(int *));
  angle    = (real **)malloc(typeNR*sizeof(real *));
  distance = (real **)malloc(typeNR*sizeof(real *));
#ifdef ENERGY
  energy   = (real **)malloc(typeNR*sizeof(real *));
#endif
  
  for(i=0;(i<typeNR);i++) {
    /* array of the number of hydrogen bonds */
    number[i]  =(int *)malloc(nr_frames*sizeof(int));  
    for(j=0;(j<nr_frames);j++) {
      number[i][j]=0;
    }
    
    /* angle and distance distributions */
    angle[i]   =(real *)malloc(RESOLUTION*sizeof(real));
    distance[i]=(real *)malloc(RESOLUTION*sizeof(real));
#ifdef ENERGY
    energy[i]=(real *)malloc(RESOLUTION*sizeof(real));
#endif
    for(j=0;(j<RESOLUTION);j++) {
      angle[i][j]=0;
      distance[i][j]=0;
#ifdef ENERGY
      energy[i][j]=0;
#endif
    }
  }
  
  /* helix gesuck van van der spoel */
  helix = (int **)malloc((eNNr+1)*sizeof(int *));
  for(i=0;(i<(eNNr+1));i++) {
    helix[i]=(int *)malloc(nr_frames*sizeof(int));
    for(j=0;(j<nr_frames);j++)
      helix[i][j]=0;
  }

  /* open file for frequency ditributions */
  fp_freq_inter = ffopen("hbmap_inter","w");
  fprintf(fp_freq_inter,"%5d %5d\n",total[INTER],nr_frames);

  fp_freq_intra = ffopen("hbmap_intra","w");
  fprintf(fp_freq_intra,"%5d %5d\n",total[INTRA],nr_frames);

  fp_freq_all = ffopen("hbmap_all","w");
  fprintf(fp_freq_all,"%5d %5d\n",total[ALL],nr_frames);
}

void List::helical(int i)
{
  int j;
  for(j=1;(j<eNNr);j++) {
    if (hbond[i]->exist() && (hbond[i]->hhb == j )) {
      helix[j][this_frame]+=1;
      helix[eNNr][this_frame]+=1;
    }
  }
}

int List::analyse()
{

  int ins,i;

  t[this_frame] = this_time;
  
  for(i=0;(i<nr_hbonds);i++) {
    real a,d=-1000,d2;
#ifdef ENERGY
    real e;
#endif
    int  n;
    

    /* a is the angle of this hydrogen bond */
    a = hbond[i]->angle();

    /* n != 0 if this hydrogen bond exists */
    n = hbond[i]->exist();

    /* d is the distance of this hydrogen bond */
    d2 = hbond[i]->distance2();

    /* e is the qoulomb energy of this hbond */
#ifdef ENERGY
    e = hbond[i]->energy();
#endif



    /* hydrogen bond analysis */
    hbond[i]->analyse();

    /* helical analysis */
    helical(i);

    /* number */
    number[hbond[i]->type][this_frame]+=n;
    number[ALL][this_frame]+=n;


    if ( n ) {
      /* calculate the sqrt of the distance */
      d = sqrt(d2);

      /* distance distribution */
      distance[hbond[i]->type][(int)((d/rcut)*RESOLUTION)]++;
      distance[ALL][(int)((d/rcut)*RESOLUTION)]++;
      
      /* angle distribution */
      angle[hbond[i]->type][(int)((a/alfcut)*RESOLUTION)]++;
      angle[ALL][(int)((a/alfcut)*RESOLUTION)]++;

#ifdef ENERGY
      /* ENERGY distribution */
      fprintf(stderr,"change list.cc in line 355 and 356 for index\n");
      exit(0);
      energy[hbond[i]->type][0]++;
      energy[ALL][0]++;
#endif
    }
       


    /* print frequency distribution */
    if (hbond[i]->type==INTRA) {
      fprintf(fp_freq_intra,"%1d",n);
      fprintf(fp_freq_all,"%1d",n);
    }
    if (hbond[i]->type==INTER) {
      fprintf(fp_freq_inter,"%1d",n);
      fprintf(fp_freq_all,"%1d",n);
    }

    
    switch (mode) {
    case SELECTED:
      /* matrix */
      matr[i][this_frame] = n ? '|' : ' ';
      
      break;
    case INSERT:
      /* matrix */
      ins = hbond[i]->insert();    
      matr[i][this_frame] = n ? ins ? '+' : '|' : ins ? '-' : ' ';
    
      break;
    case INTRA_GRP:
    
      
      break;
    case INTER_GRP:
      
      break;
    default:
      break;
    }
  }

  /* give return to the frequency tables */
  fprintf(fp_freq_intra,"\n");
  fprintf(fp_freq_inter,"\n");
  fprintf(fp_freq_all,"\n");


  return number[ALL][this_frame++];
}

/* print output to the files */
void List::dump(t_atoms *atoms)
{
  int i;
  FILE *fp;

  /* close frequency tables */
  fclose(fp_freq_intra);
  fclose(fp_freq_inter);
  fclose(fp_freq_all);


  /*number*/
  fp = xvgropen("number_total.xvg","All Hydrogen Bonds","Time (ps)","Number"); 
  for(i=0;(i<nr_frames);i++) 
    fprintf(fp,"%10g  %5d\n",t[i],number[ALL][i]);
  fflush(fp);
  fclose(fp);
  
  fp = xvgropen("number_inter.xvg","Intermolecular Hydrogen Bonds","Time (ps)","Number");
  for(i=0;(i<nr_frames);i++) 
    fprintf(fp,"%10g  %5d\n",t[i],number[INTER][i]);
  fflush(fp);
  fclose(fp);
  
  fp = (FILE *)xvgropen("number_internal.xvg","Internal Hydrogen Bonds","Time (ps)","Number");
  for(i=0;(i<nr_frames);i++) 
    fprintf(fp,"%10g  %5d\n",t[i],number[INTRA][i]);
  fflush(fp);
  fclose(fp);

  /*distance ditribution */
  fp = (FILE *)xvgropen("distance_total.xvg","All Hydrogen Bonds","Distance (nm)","Freq");
  for(i=0;(i<RESOLUTION);i++) 
    fprintf(fp,"%10g  %10g\n",rcut*((real)i/(real)RESOLUTION),distance[ALL][i]);
  fflush(fp);
  fclose(fp);

  fp = (FILE *)xvgropen("distance_inter.xvg","IntermolecularHydrogen Bonds","Distance (nm)","Freq");
  for(i=0;(i<RESOLUTION);i++) 
    fprintf(fp,"%10g  %10g\n",rcut*((real)i/(real)RESOLUTION),distance[INTER][i]);
  fflush(fp);
  fclose(fp);
  
  fp = (FILE *)xvgropen("distance_internal.xvg","Internal Hydrogen Bonds","Distance (nm)","Freq");
  for(i=0;(i<RESOLUTION);i++) 
    fprintf(fp,"%10g  %10g\n",rcut*((real)i/(real)RESOLUTION),distance[INTRA][i]);
  fflush(fp);
  fclose(fp);

  /*angle distribution */
  fp = xvgropen("angle_total.xvg","All Hydrogen Bonds","Angle (Deg)","Freq");
  for(i=0;(i<RESOLUTION);i++) 
    fprintf(fp,"%10g  %10g\n",RAD2DEG*alfcut*((real)i/(real)RESOLUTION),angle[ALL][i]);
  fflush(fp);
  fclose(fp);

  fp = xvgropen("angle_inter.xvg","Intermolecular Hydrogen Bonds","Angle (Deg)","Freq");
  for(i=0;(i<RESOLUTION);i++) 
    fprintf(fp,"%10g  %10g\n",RAD2DEG*alfcut*((real)i/(real)RESOLUTION),angle[INTER][i]);
  fflush(fp);
  fclose(fp);

  fp = xvgropen("angle_internal.xvg","Internal Hydrogen Bonds","Angle (Deg)","Freq");
  for(i=0;(i<RESOLUTION);i++) 
    fprintf(fp,"%10g  %10g\n",RAD2DEG*alfcut*((real)i/(real)RESOLUTION),angle[INTRA][i]);
  fflush(fp);
  fclose(fp);

#ifdef ENERGY
  /*energy ditribution */
  fp = (FILE *)xvgropen("energy_total.xvg","All Hydrogen Bonds","Energy (nm)","Freq");
  for(i=0;(i<RESOLUTION);i++) 
    fprintf(fp,"%10g  %10g\n",(MAX_ENER*((real)i/(real)RESOLUTION))+MIN_ENER, energy[ALL][i]);
  fflush(fp);
  fclose(fp);

  fp = (FILE *)xvgropen("energy_inter.xvg","IntermolecularHydrogen Bonds","Energy (nm)","Freq");
  for(i=0;(i<RESOLUTION);i++) 
    fprintf(fp,"%10g  %10g\n",(MAX_ENER*((real)i/(real)RESOLUTION))+MIN_ENER,energy[INTER][i]);
  fflush(fp);
  fclose(fp);
  
  fp = (FILE *)xvgropen("energy_internal.xvg","Internal Hydrogen Bonds","Energy (nm)","Freq");
  for(i=0;(i<RESOLUTION);i++) 
    fprintf(fp,"%10g  %10g\n",(MAX_ENER*((real)i/(real)RESOLUTION))+MIN_ENER,energy[INTRA][i]);
  fflush(fp);
  fclose(fp);
#endif

  /*matrix*/
  if ((mode==SELECTED)||(mode==INSERT)) {
    int i,j;
    FILE *fp;
    fp=ffopen("matrix","w");
    for(i=0;(i<nr_frames);i++) {
      fprintf(fp,"%8.3f  (",t[i]);
      for(j=0;(j<nr_hbonds);j++) 
	fprintf(fp,"%c",matr[j][i]);
      fprintf(fp,")\n");
    }
    fflush(fp);
    fclose(fp);
  }
  
  /* IF MODE IS SELECTED */
  if ( mode == SELECTED ) {
    int n;
    for(n=0;(n<nr_hbonds);n++) {
      FILE *fp;
      char line[STRLEN];
      sprintf(line,"HB_%s",hbond[n]->hb_name(atoms));
      fp = xvgropen(line,"Hydrogenbond Plot","Time (ps)","Distance (nm)");
      fprintf(fp,"@ subtitle \"%s\"\n",hbond[n]->hb_name(atoms));
      for(i=0;(i<nr_frames);i++) {
	fprintf(fp,"%10g  ",t[i]);
	hbond[n]->selected_print(i,fp);
	fprintf(fp,"\n");
      }
      fflush(fp);
      fclose(fp);
    }
  }
  
  /* water inserted graphs */
  if ( mode==INSERT ) {
    int n;
    for(n=0;(n<nr_hbonds);n++) {
      FILE *fp;
      char line[STRLEN];
      sprintf(line,"insert_%d.xvg",n);
      fp = xvgropen(line,"Water Insertion Plot","Time (ps)","Distance (nm)");
      /* print legend */
      {
	char *setname[] = { "D-A" , "D-S" , "A-S" , "Solvent Atom Number"};
	xvgr_legend(fp,4,setname);
      }
      
      /* print description of the date */
      fprintf(fp,"# Legend \n");
      fprintf(fp,"# Column | Decsription \n");
      fprintf(fp,"# 0      | Time (ps) \n");
      fprintf(fp,"# 1      | D-A : Hydrogen bond length (nm) \n");
      fprintf(fp,"# 2      | D-S : Distance of hydrogen bond Donor to solvent atom    \n");
      fprintf(fp,"# 3      | A-S : Distance of hydrogen bond Acceptor to solvent atom\n");
      fprintf(fp,"# 4      | Atom number (gromos) of solvent atom\n");
      for(i=0;(i<nr_frames);i++) {
	fprintf(fp,"%10g  ",t[i]);
	hbond[n]->insert_print(i,fp);
	fprintf(fp,"\n");
      }
      fflush(fp);
      fclose(fp);
    }
  }
  
  /* helical */
  fp = xvgropen("n-n+3.xvg","n,n+3 h-bonds","Time (ps)","Number");
  for(i=0;(i<nr_frames);i++) 
    fprintf(fp,"%10g  %5d\n",t[i],helix[NN3][i]);
  fflush(fp);
  fclose(fp);

  fp = xvgropen("n-n+4.xvg","n,n+4 h-bonds","Time (ps)","Number");
  for(i=0;(i<nr_frames);i++) 
    fprintf(fp,"%10g  %5d\n",t[i],helix[NN4][i]);
  fflush(fp);
  fclose(fp);

  fp = xvgropen("n-n+5.xvg","n,n+5 h-bonds","Time (ps)","Number");
  for(i=0;(i<nr_frames);i++) 
    fprintf(fp,"%10g  %5d\n",t[i],helix[NN5][i]);
  fflush(fp);
  fclose(fp);

  fp = xvgropen("helical.xvg","n,n+3.4.5 h-bonds","Time (ps)","Number");
  for(i=0;(i<nr_frames);i++) 
    fprintf(fp,"%10g  %5d\n",t[i],helix[eNNr][i]);
  fflush(fp);
    fclose(fp);

}



#endif	/* _list_cc */








