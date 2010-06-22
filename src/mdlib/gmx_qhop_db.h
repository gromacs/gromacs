#ifndef _GMX_QHOP_DB_H
#define _GMX_QHOP_DB_H

typedef struct{
  real  alpha, beta, gamma;
  real  k_1, k_2, k_3, m_1, m_2, m_3;
  real  s_A, t_A, v_A, s_B, s_C, t_C, v_C;
  real  f, g, h;
  real  p_1, q_1, q_2, q_3, r_1, r_2, r_3;
} t_qhop_parameters;
	
typedef struct gmx_qhop_db_t *gmx_qhop_db;

/* Return database if successful, or NULL on failure */
extern gmx_qhop_db gmx_qhop_db_read(char *forcefield);
 
/* Write the database to a filename. Return 1 on success, or 0 for
   failure */
extern int gmx_qhop_db_write(char *fn,gmx_qhop_db qdb);

/* Destroy the internal datastructures to free memory. Return 1 on
   success, 0 for failure */
extern int gmx_qhop_db_done(gmx_qhop_db qdb);

/* Return the number of states in the database for a given residue
   name: e.g. 1 for alanine, 2 for lysine, 4 for histidine. Returns
   NOTSET when the residue is not present in the database. */
extern int gmx_qhop_db_get_nstates(gmx_qhop_db qdb,char *resname);

/* Return the net charge for a given state for a given
   residue. Returns NOTSET when the residue is not in the database, or
   when the state is invalid for the residue. */
extern int gmx_qhop_db_get_qstate(gmx_qhop_db qdb,char *resname,int state);

/* Return a NULL-terminated list of atomnames of the donors in the
   residue for the indicated state. If NULL there are no donors in the
   residue. This assumes atomnames are unique, which is true for
   proteins and nucleic acids at least. */
extern char **gmx_qhop_db_get_donors(gmx_qhop_db qdb,char *resname,int state);

/* Return a NULL-terminated list of atomnames of the acceptors in the
   residue for the indicated state. If NULL there are no acceptors in
   the residue. This assumes atomnames are unique, which is true for
   proteins and nucleic acids at least. */
extern char **gmx_qhop_db_get_acceptors(gmx_qhop_db qdb,char *resname,int state);

/* Fills the array q (length natoms) with the charges corresponding to
   residue name and state. Return 1 on success, NOTSET if the resname
   is not found or the state is incorrect. */
extern int gmx_qhop_db_set_charges(gmx_qhop_db qdb,char *resname,int state,
				   int natoms,real q[]);

/* Fill the qhop_parameters for a given donor/acceptor pair. Returns 1
   if OK or 0 if either donor or acceptor does not exist. */
extern int gmx_qhop_db_get_parameters(gmx_qhop_db qdb,
				      char *donor,char *acceptor,
				      t_qhop_parameters *qp);

#endif
