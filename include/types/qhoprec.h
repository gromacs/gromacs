#ifndef QHOPREC_H
#define QHOPREC_H

/* #ifdef HAVE_CONFIG_H */
#include "config.h"
/* #include "hackblock.h" */
/* #include "atoms.h" */
#include "idef.h"
#include "../gmx_random.h"
#include "constr.h"
/* #endif */

enum {eQNONE=0, eQACC=1, eQDON=2, eQACCDON=3, eQNR=4};
/* Trivial. But note that eQACCDON = eQACC | eQDON */

enum {eQIDEF, eQBONLIB}; /* use the first to get the pointer to the idef,
			    the other for the bonded stuff stored for a
			    certain residue. */

/* Symbolic name for qhop_db. We can't include
 * gmx_qhop_types.h since that breaks compilation. */
typedef struct qhop_db *qhop_db_s; 

typedef struct {
  int donor_id, acceptor_id, proton_id,
    regime, primary_d, primary_a;     
  /* These is used for proton tautomerism
     if primary != donor/acceptor then
     we swap the coordinates of the primary
     and donor/acceptor when titrating. */
  real E12_0, E12, DE_MM, Eb,
    rda, ang, prob, hbo, kappa,
    Er, El, T;
  /* We put temperature T here, partially because one might
     want local temperatures in the future, but mainly to
     reduce the number of function arguments passed around
     in this more "object oriented" approach. The latter
     goes for most of the real data members here. */
  gmx_bool bFlip,bDonated;
  rvec xold, /* Where the acceptor proton used to be */
    xnew;    /* Where the acceptor proton may end up */
} t_hop;


/* Keeps track of where bonded interactions are in an ilist.
   We use it to quickly change bonded interactions for a residue
   upon (de)protonation. */
typedef struct {
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
  int    subres;   /* index in qhop_db.res[rtype][]. Changes upon (de)protonation. */

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

typedef struct {
  gmx_constr_t   constr; /* for some reason we can't pass constr as an
			  * argument to do_qhop(). It turns into 0x0. */
  int            nr_hop,max_nr_hop;
  t_hop          *hop;
  int            nr_qhop_atoms;
  t_qhop_atom    *qhop_atoms;
  int            nr_qhop_residues;
  t_qhop_residue *qhop_residues;
  qhop_db_s      db;
  gmx_rng_t      rng, rng_int;
  gmx_bool       bFreshNlists, bHaveEbefore;
  rvec           *f; /* Use this to avoid stupid segfaults that occur,
		      * but shouldn't have to occur, when do_force() is
		      * called with f==NULL */
  int *global_atom_to_qhop_atom;
} t_qhoprec;
#endif
