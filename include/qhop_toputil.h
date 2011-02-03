#ifndef QHOP_TOPUTIL_H
#define QHOP_TOPUTIL_H

#include "types/mdatom.h"
#include "types/topology.h"
#include "types/qhoprec.h"
#include "types/gmx_qhop_types.h"

#ifdef __cplusplus
extern "C" {
#endif

extern void set_proton_presence(qhop_H_exist *Hext, const int atomid, const gmx_bool present);

extern gmx_bool get_proton_presence(const qhop_H_exist *Hext, const int atomid);

  /* Goes through the ilists to locate the bonded
   * interactions of qres. Stores the ilist locations
   * for future hacking of bonded interactions. */
extern void index_ilists(t_qhop_residue *qres,
			 qhop_db *db,
			 const gmx_localtop_t *top,
			 const t_commrec *cr);

  /* Clears the indexing. Under the hood it does not
   * wipe the ilist locations, but only marks them as obsolete. */
extern void unindex_ilists(t_qhop_residue *qres);

/*   Exchanges the current bonded interactions*/
/*   for the ones defined by prod. */
extern void qhop_swap_bondeds(t_qhop_residue *swapres, qhop_res *prod);

/*   Exchanges the current vdw interactions
/*   for the ones defined by prod by changing the atomtypes. */
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
 * acceptors, donors, or both. */
extern void set_interactions(t_qhoprec *qr,
			     const qhop_db *qdb,
			     t_mdatoms *md,
			     t_qhop_atom *QA,
			     t_qhop_residue *qres);

/* Sets the bqhopdonor[] and bqhopacceptor[] arrays in a t_mdatoms. */
extern void qhop_atoms2md(t_mdatoms *md,
			  const t_qhoprec *qr);

/* Reads the info in the db->qhop_H_exist and finds the corresponding
 * residue subtype in the qhop_db.
 * Returns the index in db->res[rt][]
 * This function is used to set the correct protonation states
 * at the start of the simulation.
 * Right now matching is based on counting the protons on
 * titrating sites in the qhop_res and in the rtp. Thus,
 * the existence map may need slight reshuffling to really
 * represent the global protonation state. */
extern int which_subRes(const gmx_mtop_t *top,
			const t_qhoprec *qr,
			qhop_db *db,
			const int resnr);

#ifdef __cplusplus
}
#endif

#endif /* QHOP_TOPUTIL_H */
