#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

 
typedef struct {
  int atom_id; /* global atom index */
  int res_id; /* global residue index */
  char atomname[6]; /* unfortunately we need this at least for
		       histidine, but possibly for other residues as
		       well....
		    */
  char resname[6];
  bool bdonor; /* 1= donor, 0=acceptor */
  int state; /* 1= donor, 0=acceptor */
  int nr_protons; /* the number or protons that can hop away */
  int *protons;  /* the atom_id of the protons that can hop. These
		     are simply searched upon initialization. It is
		     assumed that we use complete topologies,
		     i.e. with all protons present. Protons are on/off
		     depending on their charges
		  */
  bool bWater; 
  int nr_links;
  int *links;  /* points to atoms that are part of the same residue.
		*/
  int nr_acceptors; /* known after nbsearching. */
  int *acceptors;  /* j particles that fulfil additional geometric
		      criteria */
  /* upon accepting a proton, state becomes 1; */
} t_qhop_atom;


typedef struct {
  int    *atoms; /* globael atom numbers belonging to the residue */
  int    nr_atoms; 
  int    nr_titrating_sites;
  real   **charge_set;
  int    pryte;
  int    *protonated;   /* 1/0 for each titrating atoms. Used to
			   construct a unique pryte, which indexes the
			   charge_set array.
			*/
  int    max_state; /* 2^n */
  int    state; /* 0 is competely deprotonated, 1 is fully protonated.
		   For histidines, we have to decide the state from
		   the atom name.... */
  
  int    res_nr;
} t_qhop_residue;

enum {OO,ON,NN}; /* parameters, we'll probably have to expand for
		    specific resisues
		 */


typedef struct {

  t_qhop_atom    *qhop_atoms;
  t_qhop_residue *qhop_residues;
  int            nr_qhop_residues;
  int            nr_qhop_atoms;
  int            *global_atom_to_qhop_atom;
  int            qhopfreq;  
  real           qhop_rc; /* cut of for the sphere in which veolicites will be scaled */

} t_qhoprec;
