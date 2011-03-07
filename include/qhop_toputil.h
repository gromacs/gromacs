#ifndef QHOP_TOPUTIL_H
#define QHOP_TOPUTIL_H

#include "types/mdatom.h"
#include "types/topology.h"
#include "types/qhoprec.h"
#include "hackblock.h"
#include "types/gmx_qhop_types.h"
#include "types/constr.h"
#include "types/inputrec.h"
#ifdef __cplusplus
extern "C" {
#endif
#if 0
}
#endif

/* Exchanges the qhopatoms indexed by prim and sec
 * when sec is being titrated. Since prim holds the
 * titrating hydrogen the two must be exchanged.
 * If v==NULL then velocities will remain unchanged */
extern void qhop_tautomer_swap(const t_qhoprec *qr,
			       rvec x[], rvec v[],
			       int prim, int sec);

/* Find the inert hydrogen (zero vdw params).
 * Finding it is trivial (it's the last atomtype),
 * But this function also checks that the vdw interactions
 * are zero and bails otherwise.*/
extern int find_inert_atomtype(const gmx_mtop_t *mtop, const t_forcerec *fr);

/* Flags that the proton is present in the existence function Hext */
extern void set_proton_presence(qhop_H_exist *Hext, const int atomid, const gmx_bool present);

/* Reads the existence function Hext to see wether the proton is present or not */
extern gmx_bool get_proton_presence(const qhop_H_exist *Hext, const int atomid);

  /* Goes through the ilists to locate the bonded
   * interactions of qres. Stores the ilist locations
   * for future hacking of bonded interactions. */
extern void index_ilists(t_qhop_residue *qres,
			 const qhop_db *db,
			 const gmx_localtop_t *top,
			 const t_commrec *cr);

  /* Clears the indexing. Under the hood it does not
   * wipe the ilist locations, but only marks them as obsolete. */
extern void unindex_ilists(t_qhop_residue *qres);

/* Attaches the ilib in db->rb to the end of top->idef.iparams
 * and extends top->idef.functype accordingly. */
extern void qhop_attach_ilib(gmx_localtop_t *top, const qhop_db *db);

  /* Returns the index in top->idef.il[?].iatoms where the
   * parameters for the bond involving proton_id are found. */
extern int qhop_get_proton_bond_params(const qhop_db *db, const t_qhoprec *qr,
				       t_qhop_atom *qatom, gmx_localtop_t *top,
				       int proton_id, const t_commrec *cr);

  /* Adds a constrain between the hydrogen (proton_id) and the heavy atom it's connected to.
   * Use it when activating a proton. */
extern void qhop_constrain(t_qhop_residue *qres, t_qhoprec *qr, const qhop_db *db, gmx_localtop_t *top, t_mdatoms *md, int proton_id, gmx_constr_t constr, const t_inputrec *ir, const t_commrec *cr);

  /* Adds a constrain between the hydrogen (proton_id) and the heavy atom it's connected to.
   * Use it when deactivating a proton. */
extern void qhop_deconstrain(t_qhop_residue *qres, const qhop_db *db, gmx_localtop_t *top, t_mdatoms *md, int proton_id, gmx_constr_t constr, const t_inputrec *ir, const t_commrec *cr);

/*   Exchanges the current bonded interactions*/
/*   for the ones defined by prod. */
extern void qhop_swap_bondeds(t_qhop_residue *swapres,
			      qhop_res *prod,
			      qhop_db *db,
			      gmx_localtop_t *top,
			      const t_commrec *cr);

/*   Exchanges the current vdw interactions
 *   for the ones defined by prod by changing the atomtypes. */
extern void qhop_swap_vdws(const t_qhop_residue *swapres,
			   const qhop_res *prod,
			   t_mdatoms *md,
			   const qhop_db *db);

/*   Exchanges the current masses, invmasses and charges */
/*   for the ones defined by prod. */
extern void qhop_swap_m_and_q(const t_qhop_residue *swapres,
			      const qhop_res *prod,
			      t_mdatoms *md,
			      const qhop_db *db, t_qhoprec *qr);

/* Set the interaction parameters and
 * determine whether titraing sites are
 * acceptors, donors, or both.
 * This function is only called when 
 * initializing qhop. */
extern void set_interactions(t_qhoprec *qr,
			     qhop_db *qdb,
			     gmx_localtop_t *top,
			     t_mdatoms *md,
			     t_qhop_residue *qres,
			     t_commrec *cr);

/* Sets the bqhopdonor[] and bqhopacceptor[] arrays in a t_mdatoms. */
/* extern void qhop_atoms2md(t_mdatoms *md, */
/* 			  const t_qhoprec *qr); */

/* Reads the info in the db->qhop_H_exist and finds the corresponding
 * residue subtype in the qhop_db.
 * Returns the index in db->res[rt][]
 * This function is used to set the correct protonation states
 * at the start of the simulation.
 * Right now matching is based on counting the protons on
 * titrating sites in the qhop_res and in the rtp. Thus,
 * the existence map may need slight reshuffling to really
 * represent the global protonation state. */
extern int which_subRes(const t_qhoprec *qr,
			qhop_db *db,
			const int resnr);

#if 0
{
#endif
#ifdef __cplusplus
}
#endif

#endif /* QHOP_TOPUTIL_H */
