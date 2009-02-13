#ifndef _GMX_QHOP_PARM_H
#define _GMX_QHOP_PARM_H
	
typedef struct gmx_qhop_t *gmx_qhop;

/* Return a new gmx_qhop structure */
extern gmx_qhop gmx_qhop_init();

/* These function get and set the obvious */
extern void gmx_qhop_set_donor(gmx_qhop gqh,char *donor);

extern void gmx_qhop_set_acceptor(gmx_qhop gqh,char *acceptor);

extern char *gmx_qhop_get_donor(gmx_qhop gqh);

extern char *gmx_qhop_get_acceptor(gmx_qhop gqh);

/* Add parameter to gqh, return 1 if OK, 0 if not OK */
extern int gmx_qhop_add_param(gmx_qhop gqh,char *name,char *value,char *unit);

/* Lists the parameters, one by one on repeatedly calling the
   function. Returns 1 if OK, 0 if not OK */
extern int gmx_qhop_get_param(gmx_qhop gqh,char **name,char **value,char **unit);

/* Return a value corresponding to name in *x. Return 1 of OK, 0 if
   not OK */
extern int gmx_qhop_get_value(gmx_qhop gqh,char *name,double *x);

/* Liberate memory */
extern void gmx_qhop_done(gmx_qhop gqh);

#endif
