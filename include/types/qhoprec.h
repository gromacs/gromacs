#ifndef QHOPREC_H
#define QHOPREC_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#include "idef.h"
#endif

enum {eQNONE=0, eQACC=1, eQDON=2, eQACCDON=3, eQNR=4};
/* Trivial. But note that eQACCDON = eQACC|eQDON */

enum {eQIDEF, eQBONLIB}; /* use the first to get the pointer to the idef,
			   the other for the bonded stuff stored for a
			   certain residue. */

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
    El;
  gmx_bool bFlip;
} t_hop;


/* relates position in idef to a position in a bondeds library for a residue. */
typedef struct {
  t_iatom*  ilib;
  t_iatom** reslib; /* points to entry in resblocks.
		       Use it to find bonded interaction
		       for a specific protonation state. */
} qhop_bonswap;

/* Keeps track of where bonded interactions are in an ilist.
   We use it to quickly change bonded interactions for a residue
   upon (de)protonation. */
typedef struct {
  int     nb;
  qhop_bonswap *bondeds; /* nb elements long. */
} qhop_bonded_index;

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
  int    **ilistPosG; /* Where in the global topology are the
		       * bonded interactions found?
		       * Indexes a t_ilist[F_*][interaction] */
  int    **ilistPosL; /* Where in the local topology are the
		       * bonded interactions found?
		       * Indexes a t_ilist[F_*][interaction] */
  qhop_bonded_index bindex; /* Makes ilistPosG and ilistPosL redundant. For now anyway. */
} t_qhop_residue;

/* This enum can probably go, because thigns are now explicit in the xml input. */
/* enum {OO,ON,NN}; /\* parameters, we'll probably have to expand for */
/* 		    specific resisues */
/* 		 *\/ */


typedef struct {

  t_hop          *hop;
  t_qhop_atom    *qhop_atoms;
  t_qhop_residue *qhop_residues;
  
  int nr_hops,
    nr_qhop_residues,
    nr_qhop_atoms,
    *global_atom_to_qhop_atom,
    qhopfreq,
    qhopmode,
    qhop_rc; /* cut of for the sphere in which veolicites will be scaled */

} t_qhoprec;
#endif
