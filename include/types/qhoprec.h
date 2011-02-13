#ifndef QHOPREC_H
#define QHOPREC_H

/* #ifdef HAVE_CONFIG_H */
#include "config.h"
/* #include "hackblock.h" */
/* #include "atoms.h" */
#include "idef.h"
#include "../gmx_random.h"
/* #endif */

enum {eQNONE=0, eQACC=1, eQDON=2, eQACCDON=3, eQNR=4};
/* Trivial. But note that eQACCDON = eQACC|eQDON */

enum {eQIDEF, eQBONLIB}; /* use the first to get the pointer to the idef,
			   the other for the bonded stuff stored for a
			   certain residue. */

typedef struct qhop_db *qhop_db_s; /* Symbolic name for qhop_db. We can't include
				    * gmx_qhop_types.h since that breaks compilation. */

typedef struct {
  int donor_id,
    acceptor_id,
    proton_id,
    regime;
  real E12,
    DE_MM,
    Eb,
    rda,
    prob,
    hbo,
    kappa,
    Er,
    El,
    T;
  /* We put temperature T here, partially because one might
     want local temperatures in the future, but mainly to
     reduce the number of function arguments passed around
     in this more "object oriented" approach. The latter
     goes for most of the real data members here. */
  gmx_bool bFlip;
} t_hop;


/* Keeps track of where bonded interactions are in an ilist.
   We use it to quickly change bonded interactions for a residue
   upon (de)protonation. */
typedef struct {
  int     nb;
  t_iatom **ilist_pos; /* dimensions: [F_NRE][#bonded_of this type] */

  int nr[F_NRE]; /* How many interactions in ilist_pos[] arrays? */

  gmx_bool **indexed; /* Are the bonded interactions indexed.
		       * Dimensions: [F_NRE][#bonded_of this type] */

} qhop_bonded_index;

typedef struct {
  int atom_id; /* global atom index */
  int res_id;  /* global residue index */
  int qres_id; /* qhop_residue index */

  char *atomname, *resname; /* will point to names in the depths of qhop_resblocks */

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

  int nr_acceptors; /* known after nbsearching. */
  int *acceptors;   /* j particles that fulfil additional geometric
		      criteria */
  /* upon accepting a proton, state becomes eQDON or eQACCDON; */
} t_qhop_atom;

typedef struct {
  /* rtype and res are needed to switch between interaction parameters. */
  int    rtype; /* index in qhop_db.restype[]  and qhop_db.res[]. Will not change upon (de)protonation */
  int    res;   /* index in qhop_db.res[rtype][]. Changes upon (de)protonation. */

  int    *atoms; /* global atom numbers belonging to the residue */
  int    nr_atoms;
  char   **atomnames;
  int    nr_titrating_sites;
  int    *titrating_sites; /* points back to acceptor donors in the qhop_atoms structure. */

  int    state; /* 0 is competely deprotonated, 1 is fully protonated.
		   For histidines, we have to decide the state from
		   the atom name.... */
  /* Let's change to eQACC for deprotonated, eQACCDON for ampholytic and eQDON for fully protonated. */
  
  int    res_nr; /* Global resnr. */

  qhop_bonded_index bindex; /* Makes ilistPosG and ilistPosL redundant. For now anyway. */
  int nr_indexed;      /* are the ilists indexed? */
} t_qhop_residue;

/* This enum can probably go, because thigns are now explicit in the xml input. */
/* enum {OO,ON,NN}; /\* parameters, we'll probably have to expand for */
/* 		    specific resisues */
/* 		 *\/ */


typedef struct {

  t_hop          *hop;
  t_qhop_atom    *qhop_atoms;
  t_qhop_residue *qhop_residues;
  qhop_db_s      db;
  gmx_rng_t      rng, rng_int;
  gmx_bool       bFreshNlists;
  int nr_hops,
    nr_qhop_residues,
    nr_qhop_atoms,
    *global_atom_to_qhop_atom,
    qhopfreq,
    qhopmode,
    qhop_rc; /* cut of for the sphere in which veolicites will be scaled */

} t_qhoprec;
#endif
