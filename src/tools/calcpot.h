/*
 * $Id$
 * 
 *       This source code is part of
 * 
 *        G   R   O   M   A   C   S
 * 
 * GROningen MAchine for Chemical Simulations
 * 
 *               VERSION 2.0
 * 
 * Copyright (c) 1991-1999
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 * 
 * Also check out our WWW page:
 * http://md.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 * 
 * And Hey:
 * Great Red Oystrich Makes All Chemists Sane
 */
static char *SRCID_calcpot_h = "$Id$";

	
extern void init_calcpot(int nfile,t_filenm fnm[],t_topology *top,
			 rvec **x,t_parm *parm,t_commrec *cr,
			 t_graph **graph,t_mdatoms **mdatoms,
			 t_nsborder *nsb,t_groups *grps,
			 t_forcerec **fr,real **coulomb,
			 matrix box);

extern void calc_pot(FILE *logf,t_nsborder *nsb,t_commrec *cr,t_groups *grps,
		     t_parm *parm,t_topology *top,rvec x[],t_forcerec *fr,
		     t_mdatoms *mdatoms,real coulomb[]);

extern void write_pdb_coul();

extern void delete_atom(t_topology *top,int inr);
/* Delete an atom from a topology */

extern void replace_atom(t_topology *top,int inr,char *anm,char *resnm,
			 real q,real m,int type);
/* Replace an atom in a topology by someting else */

