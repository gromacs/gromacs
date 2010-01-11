#ifndef _GMX_QHOP_DB_H
#define _GMX_QHOP_DB_H

#include "resall.h"
#include "types/gmx_qhop_types.h"
#include "hackblock.h"

//typedef struct gmx_qhop_db_t *gmx_qhop_db;

/*Print to stderr*/
extern void qhop_db_print (qhop_parameters *qhp);

/* Return database if succesfull, or NULL on failure */
extern qhop_db_t qhop_db_read(char *forcefield, gmx_mtop_t *top);
 
/* Write the database to a filename. Return 1 on succes, or 0 for
   failure */
extern int qhop_db_write(char *fn,qhop_db_t qdb);

/* Destroy the internal datastructures to free memory. Return 1 on
   succes, 0 for failure */
extern int qhop_db_done(qhop_db_t qdb);

/* Return the number of states in the database for a given residue
   name: e.g. 1 for alanine, 2 for lysine, 4 for histidine. Returns
   NOTSET when the residue is not present in the database. */
extern int qhop_db_get_nstates(qhop_db_t qdb,char *resname);

/* Return the net charge for a given state for a given
   residue. Returns NOTSET when the residue is not in the database, or
   when the state is invalid for the residue. */
extern int qhop_db_get_qstate(qhop_db_t qdb,char *resname,int state);

/* Return a NULL-terminated list of atomnames of the donors in the
   residue for the indicated state. If NULL there are no donors in the
   residue. This assumes atomnames are unique, which is true for
   proteins and nucleic acids at least. */
extern char **qhop_db_get_donors(qhop_db_t qdb,char *resname,int state);

/* Return a NULL-terminated list of atomnames of the acceptors in the
   residue for the indicated state. If NULL there are no acceptors in
   the residue. This assumes atomnames are unique, which is true for
   proteins and nucleic acids at least. */
extern char **qhop_db_get_acceptors(qhop_db_t qdb,char *resname,int state);

/* Fills the array q (length natoms) with the charges corresponding to
   residue name and state. Return 1 on success, NOTSET if the resname
   is not found or the state is incorrect. */
extern int qhop_db_set_charges(qhop_db_t qdb,char *resname,int state,
			       int natoms,real q[]);

/* Fill the qhop_parameters for a given donor/acceptor pair. Returns 1
   if OK or 0 if either donor or acceptor does not exist. */
extern int qhop_db_get_parameters(qhop_db_t qdb,
				  char *donor,char *acceptor,
				  qhop_parameters *qp);

#endif
