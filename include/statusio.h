/*
 *       $Id$
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 1.6
 * 
 * Copyright (c) 1991-1997
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 *
 * Also check out our WWW page:
 * http://rugmd0.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 *
 * And Hey:
 * Good gRace! Old Maple Actually Chews Slate
 */

#ifndef _statusio_h
#define _statusio_h

static char *SRCID_statusio_h = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) statusio.h 1.46 2/2/97"
#endif /* HAVE_IDENT */

#include "typedefs.h"
#include "sheader.h"

/*
 * This module handles status file io. All read and write operations from
 * and to a status file should use this functions to be independent of the
 * actual file layout (text versus binary file).
 */
#ifdef CPLUSPLUS
extern "C" { 
#endif

extern size_t wr_status(FILE *fp,int step,real t,real lambda,
			t_inputrec *ir,rvec *box,rvec *vir,rvec *pres,
			int natoms,rvec *x,rvec *v,rvec *f,
			int nre,t_energy *e,t_topology *top);
/*
 * Writes a complete status to the file, specified by fp. NULL pointers
 * indicate that this field should not be written. The function returns
 * the number of bytes written.
 */

extern char *rd_hstatus(FILE *fp,t_statheader *sh,int *step,real *t,
                        real *lambda,t_inputrec *ir,rvec *box,
			rvec *vir,rvec *pres,int *natoms,
                        rvec *x,rvec *v,rvec *f,int *nre,t_energy *e,
                        t_topology *top);
/*
 * Reads a complete status from the file, specified by fp. It uses
 * the status header to find the items in the file, also the file
 * should be positioned right for reading the first item. The function
 * returns the version string from the header.
 */

extern char *rd_status(FILE *fp,int *step,real *t,real *lambda,
                       t_inputrec *ir,rvec *box,rvec *vir,rvec *pres,
		       int *natoms,rvec *x,
                       rvec *v,rvec *f,int *nre,t_energy *e,
                       t_topology *top);
/*
 * Reads a complete status from the file, specified by fp. First it
 * reads the header and then invokes rd_hstatus() to read the rest
 * of the status. It returns the version returned from rd_hstatus().
 */

extern void write_status(char *fn,int step,real t,real lambda,t_inputrec *ir,
                         rvec *box,rvec *vir,rvec *pres,
			 int natoms,rvec *x,rvec *v,rvec *f,
                         int nre,t_energy *e,t_topology *top);
/*
 * Writes a complete status to the file, specified by fn. NULL pointers
 * indicate that this field should not be written.
 */

extern char *read_status(char *fn,int *step,real *t,real *lambda,
                         t_inputrec *ir,rvec *box,rvec *vir,rvec *pres,
			 int *natoms,rvec *x,
                         rvec *v,rvec *f,int *nre,t_energy *e,
                         t_topology *top);
/*
 * Reads a complete status from the file, specified by fn. It returns
 * the version returned from rd_hstatus().
 */

extern void read_status_header(char *fn,t_statheader *header);
/*
 * Reads a (status) header from the file, specified by fn. If
 * available, it returns the version string from the file, else
 * it returns a version string from the statusio module.
 */

#ifdef CPLUSPLUS
}
#endif

#endif	/* _statusio_h */
