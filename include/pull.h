#include "vec.h"
#include "typedefs.h"

#define NINT(x) (x<0?((int)((x)-0.5)):((int)((x)+0.5)))
#define DEBUG_START fprintf(stderr,"\n\nDEBUG\n");
#define DEBUG_END   fprintf(stderr,"\nEND DEBUG\n\n");

/* print to output file for the various types of runs */
extern void print_umbrella(t_pull *pull, int step);
extern void print_afm(t_pull *pull, int step);
extern void print_constraint(t_pull *pull,rvec *force,int step,matrix box);
extern void print_start(t_pull *pull, int step);


/* calculate center of mass of index group, making sure it's inside the box */
extern real calc_com(rvec x[], 
		     int gnx, 
		     atom_id *index, 
		     t_mdatoms *md, 
		     rvec com, 
		     matrix box);


/* calculate center of mass of all atoms x[], index needed to get the right
   masses from the atom array */
extern real calc_com2(rvec x[], 
		      int gnx, 
		      atom_id *index, 
		      t_mdatoms *md, 
		      rvec com, 
		      matrix box);


/* calculate a running average for center of mass */
extern void calc_running_com(t_pull *pull);


/* calculate the center of mass from the true coordinates, without
   corrections for pbc */
extern void correct_t0_pbc(t_pull *pull, 
			   rvec x[], 
			   t_mdatoms *md,
			   matrix box);


/* parse a string for 3 numbers and put them in rvec */
extern void string2rvec(char *buf, 
			rvec x);


extern void read_pullparams(t_pull *pull,
			    char *infile); 


/* find all atoms in group pull->idx[pull->n] that are inside a cylinder
   with as origin com[i][x],com[i][y] with radius pull->r and possibly
   a switch function pull->rc. Remember their weight. Now each group i
   has its own reference group (HOW?) with com defined as 
   Sum(wi*mi*ri)/(Sum(wi*mi). Basically, build an index structure with
   the reference groups for the groups i, plus an array with the 
   weight factors for each of the atoms in those index groups? 
   */
extern void make_refgrps(t_pull *pull,
			 matrix box,
			 t_mdatoms *md);


/* write a numbered .gro file in procedure to make starting structures */
extern void dump_conf(t_pull *pull,
		      rvec x[],
		      matrix box,
		      t_topology *top,
		      int nout, 
		      real time);


/* main pull routine that controls all the action */
extern void pull(t_pull *pull,            /* all pull data */
		 rvec *x,                 
		 rvec *f,                 
		 matrix box,              
		 t_topology *top,         
		 real dt, 
		 int step, 
		 int natoms, 
		 t_mdatoms *md);


/* get memory and initialize the fields of pull that still need it, and
   do runtype specific initialization */
extern void init_pull(FILE *log, 
		      int nfile,
		      t_filenm fnm[], 
		      t_pull *pull,
		      rvec *x, 
		      t_mdatoms *md, 
		      rvec boxsize);







