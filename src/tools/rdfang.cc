/*
 *       $Id$
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0
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
 * GROtesk MACabre and Sinister
 */
static char *SRCID_rdfang_cc = "$Id$";

#include <math.h>
#include "xvgr.h"
#include "sysstuff.h"
#include "do_rdf.h"
	
#define NTH        12      /* Number of slices in theta space      */

typedef class c_rdfang: public c_rdf {
  int  axis;
  real **graph;
public:
  c_rdfang(int n1,char *gnm1,int n2,char *gnm2,int ax);
  ~c_rdfang();
  void init_graph(void);
  void insert(rvec dx,real r2);
  void print(char *outf1,char *outf2);
} c_rdfang;

c_rdfang::c_rdfang(int n1,char *gnm1,int n2,char *gnm2,int ax) 
: c_rdf(n1,gnm1,n2,gnm2,FALSE)
{
  axis=ax;
  graph=NULL;
}

c_rdfang::~c_rdfang()
{
  int i;
  
  for(i=0; (i<ng); i++)
    free(graph[i]);
  free(graph);
}

void c_rdfang::init_graph(void)
{
  int i;
  
  graph=(real **)calloc(ng,sizeof(graph[0]));
  for(i=0; (i<ng); i++)
    graph[i]=(real *)calloc(NTH,sizeof(graph[0][0]));
}



void c_rdfang::insert(rvec dx,real r2)
{
  int  n,m;
  
  n=(int)sqrt(r2*rseg_2);
  m=(int)((NTH/2)*(1-(dx[axis]/sqrt(r2))));
  graph[n][m]+=1.0;
}

static real ring_vol(real r, real dcosth)
{
  return (2.0*M_PI/3.0)*(r*r*r)*dcosth;
}

void c_rdfang::print(char *outf1,char *outf2)
{
  FILE *out;
  int  i,j;
  real dang,dcosth,vol,gem_gr;
  real rho[NTH];

  /* calculate the angles and costheta intervals */
  dang=M_PI/NTH;
  dcosth=2.0/NTH;
  
  /* Calculate the density in pairs/vol */
  vol=ring_vol((segsize)*(ng),dcosth);
  for(j=0;(j<NTH); j++){
    rho[j]=0;
    for(i=0; (i<ng); i++)
      rho[j]+=graph[i][j];
    rho[j]/=vol;
  }

  out=xvgropen(outf1,"Angle Dependent RDF","r (nm)","g(r)");
  for(j=0; (j<NTH); j++)
    fprintf(out,"@    legend string %2d \"theta = %5.3f \\8p\\0\"\n",
	    j,(j*dang));
  /*  fprintf(out,"@    legend on\n"); */
  for(j=0; (j<NTH); j++){
    for(i=0; (i<ng); i++){
      fprintf(out,"%12.5f  %12.5f\n",i*segsize,
	      (graph[i][j]*inv_segvol[i])/rho[j]);
    }
    fprintf(out,"\n");
  }
  fclose(out);

  /* print only one halve of the symmetric slices (theta=0..PI/2) */
  out=xvgropen(outf2,"Angle Dependent RDF","r (nm)","g(r)");
  for(j=0; (j<NTH); j++)
    fprintf(out,"@    legend string %2d \"theta = %5.3f \\8p\\0\"\n",
	    j,(j*dang));
  /*  fprintf(out,"@    legend on\n"); */
  fprintf(out,"@     s0 linestyle 1\n");
  fprintf(out,"@     s%1d linestyle 3\n",((NTH)/2 - 1));
  for(j=0; (j<((NTH)/2)); j++){
    for(i=0; (i<ng); i++){
      gem_gr=0.5*(graph[i][j] + graph[i][NTH-1-j]);
      fprintf(out,"%12.5f  %12.5f\n",i*segsize,
	      (gem_gr*inv_segvol[i])/rho[j]);
    }
    fprintf(out,"\n");
  }
  fclose(out);
}

void do_rdfang(char *fn,char *outf1,char *outf2,
	       real cutoff,t_block *excl,
	       atom_id index1[],int n1,char *gnm1,
	       atom_id index2[],int n2,char *gnm2,
	       int axis)
{
  c_rdfang rdfang(n1,gnm1,n2,gnm2,axis);
  
  rdfang.calc(fn,cutoff,excl,index1,index2);
  rdfang.print(outf1,outf2);
}

