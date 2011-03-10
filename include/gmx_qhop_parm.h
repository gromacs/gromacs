#ifndef _GMX_QHOP_PARM_H
#define _GMX_QHOP_PARM_H

#include "types/gmx_qhop_types.h"

/* npbon[][][] contains the number of parameters for a certain type of bonded interacion.
 * The outmost dimension denotes the A- and B-state respectively. */
static const int npbon[3][9][2] = {
  /*  bond    G96     morse   cubic   conn    harm    fene    table  tabnc */
  {{2, 4}, {2, 4}, {3, 3}, {3, 3}, {0, 0}, {2, 4}, {2, 2}, {2, 3}, {2, 3}}, /*  <-------------- Bonds */
  /* harm    G96     cr_bb    cr_ba    UB     qang     -     table     -   */
  {{2, 4}, {2, 4}, {3, 3}, {4, 4}, {4, 4}, {6, 6}, {0, 0}, {2, 3}, {0, 0}}, /*  <-------------- Angles */
  /*  prop    imp     RB                        -      fou     -     -     table   prop9 */
  {{3, 5}, {2, 4}, {NR_RBDIHS, 2*NR_RBDIHS}, {0, 0}, {4, 8}, {0,0},{0,0}, {2, 3}, {3, 5}} /* <- Dihedrals */
  };

/* These enums serve to find the correct element in npbon */
enum {ebtypeBOND, ebtypeANGLE, ebtypeDIHEDRAL};
enum {eAstate, eBstate};
  


/* *********************
 * qhop resblock stuff *
 * ********************* */

/* Returns a pointer to fresh qhop_resblock */
extern qhop_resblocks_t qhop_resblock_init();

extern void qhop_add_restype(qhop_resblocks_t rb, char *name, int nres, 
			     qhop_subres_t res);

extern void qhop_add_res(qhop_resblocks_t rb, int resblocknr, 
			 qhop_subres_t res, int nres);

/* Sets the interactions of qres according to the (de)protonation at hydrogen H.
 * If H is not present, it's protonation, otherwise deprotonaton.
 * db->Hmap.H is updated.*/
extern void qhop_set_protonation(const qhop_db *db, t_qhop_residue *qres,
				 const atom_id H);

extern void make_ilib(qhop_db *db);

/* **********************
 * qhop parameter stuff *
 * ********************** */

/* Return a new qhop structure */
extern qhop_t qhop_init();

/* These function get and set the obvious */
extern void qhop_set_donor(qhop_t gqh, const char *donor);

extern void qhop_set_acceptor(qhop_t gqh, const char *acceptor);

extern void qhop_set_don_atom(qhop_t gqh, const char *donor);

extern void qhop_set_acc_atom(qhop_t gqh, const char *acceptor);

extern char *qhop_get_don_atom(const qhop_t gqh);

extern char *qhop_get_acc_atom(const qhop_t gqh);

extern char *qhop_get_donor(const qhop_t gqh);

extern char *qhop_get_acceptor(const qhop_t gqh);

/* Add parameter to gqh, return 1 if OK, 0 if not OK */
extern int qhop_add_param(qhop_t gqh,char *name,char *value,char *unit);

/* Lists the parameters, one by one on repeatedly calling the
   function. Returns 1 if OK, 0 if not OK */
extern int qhop_get_param(qhop_t gqh,char **name,char **value,char **unit);

/* Return a value corresponding to name in *x. Return 1 of OK, 0 if
   not OK */
extern int qhop_get_value(qhop_t gqh,char *name,double *x);

/* Liberate memory */
extern void qhop_done(qhop_t gqh);

#endif
