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
 * GROtesk MACabre and Sinister
 */
static char *SRCID_do_rdf_h = "$Id$";

#include "typedefs.h"
#include "vec.h"
		
class c_rdf {
 public:
  int 	nx1,nx2;	/* The number of atoms in both groups 	*/
  char  *nm1,*nm2;      /* The names of the groups              */
  bool  bCoM;           /* Use center of mass as central particl*/
  int 	ng;		/* The number of point in this graph	*/
  real  hbox2;          /* The half-box squared                 */
  real	segsize;	/* The size of the segment in r-space	*/
  real  rseg_2;         /* Inverse square of prev.              */
  real	*inv_segvol;	/* Inverse segvol_size			*/
public:
  c_rdf(int n1,char *gnm1,int n2,char *gnm2,bool bCM);
  virtual ~c_rdf();
  virtual void insert(rvec dx,real r2) = 0;
  virtual void init_graph() = 0;
  void init_box(matrix box);
  void add(real cutoff,t_block *excl,matrix box,rvec x[], 
	   atom_id index1[],atom_id index2[]);
  void calc(char *fn,real cutoff,t_block *excl,
	    atom_id index1[],atom_id index2[]);
};

extern real sphere_vol(real r);

extern void do_rdfn(char *fn,char *outfile,real cutoff,t_block *excl,bool bCM,
		    atom_id index1[],int n1,char *gnm1,
		    atom_id index2[],int n2,char *gnm2);

extern void do_rdfang(char *fn,char *outf1,char *outf2,
		      real cutoff,t_block *excl,
		      atom_id index1[],int n1,char *gnm1,
		      atom_id index2[],int n2,char *gnm2,
		      int axis);
