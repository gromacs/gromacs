#include "typedefs.h"
	
extern void init_calcpot(int nfile,t_filenm fnm[],t_topology *top,
			 rvec **x,t_parm *parm,t_commrec *cr,
			 t_graph **graph,t_mdatoms **mdatoms,
			 t_nsborder *nsb,t_groups *grps,
			 t_forcerec **fr,real **coulomb);

extern void calc_pot(FILE *logf,t_nsborder *nsb,t_commrec *cr,t_groups *grps,
		     t_parm *parm,t_topology *top,rvec x[],t_forcerec *fr,
		     t_graph *graph,t_mdatoms *mdatoms,real coulomb[]);

extern void write_pdb_coul();

extern void delete_atom(t_topology *top,int inr);
/* Delete an atom from a topology */

extern void replace_atom(t_topology *top,int inr,char *anm,char *resnm,
			 real q,real m,int type);
/* Replace an atom in a topology by someting else */

