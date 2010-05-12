#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

enum {eQNONE, eQACC, eQDON, eQACCDON, eQNR};
 
typedef struct {
  int atom_id; /* global atom index */
  int res_id; /* global residue index */
  int qres_id;
/*   char atomname[6]; /\* unfortunately we need this at least for */
/* 		       histidine, but possibly for other residues as */
/* 		       well.... */
/* 		    *\/ */
/*   char resname[6]; */

  char *atomname, *resname; /* will point to names in the depths of qhop_resblocks */

  /* We need to consider ampholytic species, such as water(!),
   * which is why we can't have a bool bdonor anymore. */
  /* bool bdonor; /\* 1= donor, 0=acceptor *\/ */
  int state; /* None, acceptor, donor or both? This can perhaps speed up the search for potential hops.
	      * Takes on values eQNONE, eQACC, eQDON, eQACCDON. */
  
  int nr_protons; /* the number or protons that can hop away */
  int *protons;  /* the atom_id of the protons that can hop. These
		    are simply searched upon initialization. It is
	            assumed that we use complete topologies,
	            i.e. with all protons present. Protons are on/off
	            depending on their charges.

		    No they are not anymore. There is an exisence array for all hydrogens in the system
		    in qhop_db (maybe it should move to t_qhoprec?).
		  */
/*   bool bWater;  /\* This is ugly and should probably go. *\/ */
  /* Is the links stuff needed? */
/*   int nr_links; */
/*   int *links;  /\* points to atoms that are part of the same residue. */
/* 		*\/ */
  int nr_acceptors; /* known after nbsearching. */
  int *acceptors;  /* j particles that fulfil additional geometric
		      criteria */
  /* upon accepting a proton, state becomes eQDON or eQACCDON; */
} t_qhop_atom;


typedef struct {
  /* rtype and res are needed to switch between interaction parameters. */
  int    rtype; /* index in qhop_db.restype[]  and qhop_db.res[]. Will not change upon (de)protonation */
  int    res;   /* index in qhop_db.res[rtype][]. Changes upon (de)protonation. */

  int    *atoms; /* global atom numbers belonging to the residue */
  int    nr_atoms;
  int    nr_titrating_sites;
  int    *titrating_sites; /* points back to acceptor donors in the qhop_atoms structure. */
  /* real   **charge_set; 0*/
/*   int    pryte; */
/*   int    *protonated;   /\* 1/0 for each titrating atoms. Used to */
/* 			   construct a unique pryte, which indexes the */
/* 			   charge_set array. */
/* 			*\/ */
  /* Why max state? Can we take this away? */
/*   int    max_state; /\* 2^n *\/ */
  /* State seems obsolete. */
  /* Well, it could speed up things when looking for possible hops. */
  int    state; /* 0 is competely deprotonated, 1 is fully protonated.
		   For histidines, we have to decide the state from
		   the atom name.... */
  /* Let's change to eQACC for deprotonated, eQACCDON for ampholytic and eQDON for fully protonated. */
  
  int    res_nr; /* Global resnr. */
} t_qhop_residue;

/* This enum can probably go, because thins are now explicit in the xml input. */
/* enum {OO,ON,NN}; /\* parameters, we'll probably have to expand for */
/* 		    specific resisues */
/* 		 *\/ */


typedef struct {

  t_qhop_atom    *qhop_atoms;
  t_qhop_residue *qhop_residues;
  int            nr_qhop_residues;
  int            nr_qhop_atoms;
  int            *global_atom_to_qhop_atom;
  int            qhopfreq;
  real           qhop_rc; /* cut of for the sphere in which veolicites will be scaled */

} t_qhoprec;
