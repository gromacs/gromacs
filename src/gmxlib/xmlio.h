#include "typedefs.h"

extern void read_xml(char *fn,int *step,real *t,real *lambda,
		     t_inputrec *ir,rvec *box,int *natoms,
		     rvec **x,rvec **v,rvec **f,t_topology *top);
extern void write_xml(char *fn,char *title,t_inputrec *ir,rvec *box,
		      int natoms,rvec *x,rvec *v,rvec *f,
		      int nmol,t_atoms atoms[],t_idef *idef);

