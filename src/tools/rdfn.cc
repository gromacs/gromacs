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
static char *SRCID_rdfn_cc = "$Id$";

#include <math.h>
#include "xvgr.h"
#include "sysstuff.h"
#include "typedefs.h"
#include "do_rdf.h"
#include "maths.h"
#include "vec.h"
	
class c_rdfn: public c_rdf {
private:
  real *graph;
public:
  c_rdfn(int n1,char *gnm1,int n2,char *gnm2,bool bCM);
  ~c_rdfn();
  void init_graph(void);
  void insert(rvec dx,real r2);
  void print(char *outf);
};

c_rdfn::c_rdfn(int n1,char *gnm1,int n2,char *gnm2,bool bCM) 
: c_rdf(n1,gnm1,n2,gnm2,bCM)
{
  graph=NULL;
}

c_rdfn::~c_rdfn()
{
  if (graph)
    free(graph);
}

void c_rdfn::init_graph(void)
{
  graph=(real *)calloc(ng,sizeof(graph[0]));
}

void c_rdfn::insert(rvec dx,real r2)
{
  int  n;
  
  n=(int)sqrt(r2*rseg_2);
  graph[n]+=1.0;
}

void c_rdfn::print(char *outf)
{
  FILE *out;
  int  i;
  real rho;
  real rrr,total;
  
  /* Calculate the density in pairs/vol */
  rho=0.0;
  for(i=0; (i<ng); i++)
    rho+=graph[i];
  rho /= sphere_vol(ng*segsize);

  out=xvgropen(outf,"Gromacs RDF","r (nm)","g(r)");
  fprintf(out,"@ subtitle \"%s-%s\"\n",nm1,nm2);
  total = 0.0;
  for(i=0; (i<ng-1); i++) {
    rrr=(graph[i]*inv_segvol[i])/rho;
    total+=4*M_PI*rrr*sqr(i*segsize)*segsize;
    fprintf(out,"%12.5f  %12.5f  %12.5f\n",i*segsize,
	    rrr,total);
  }
  fclose(out);
}

void do_rdfn(char *fn,char *outfile,real cutoff,t_block *excl,bool bCM,
	     atom_id index1[],int n1,char *gnm1,
	     atom_id index2[],int n2,char *gnm2)
{
  c_rdfn  rdfn(n1,gnm1,n2,gnm2,bCM);
  
  rdfn.calc(fn,cutoff,excl,index1,index2);
  rdfn.print(outfile);
}

