#ifndef _GMX_QHOP_PARM_H
#define _GMX_QHOP_PARM_H

#include "types/gmx_qhop_types.h"

/* *********************
 * qhop resblock stuff *
 * ********************* */

/* Returns a pointer to fresh qhop_resblock */
extern qhop_resblocks_t qhop_resblock_init();

extern void qhop_add_restype(qhop_resblocks_t rb, char *name, int nres, qhop_res_t res);

extern void qhop_add_res(qhop_resblocks_t rb, int resblocknr, qhop_res_t res, int nres);

extern void qhop_set_protonation(qhop_resblocks_t rb, qhop_res_t res);

/* **********************
 * qhop parameter stuff *
 * ********************** */

/* Return a new qhop structure */
extern qhop_t qhop_init();

/* These function get and set the obvious */
extern void qhop_set_donor(qhop_t gqh,char *donor);

extern void qhop_set_acceptor(qhop_t gqh,char *acceptor);

extern void qhop_set_don_atom(qhop_t gqh,char *donor);

extern void qhop_set_acc_atom(qhop_t gqh,char *acceptor);

extern char *qhop_get_donor(qhop_t gqh);

extern char *qhop_get_acceptor(qhop_t gqh);

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
