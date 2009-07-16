/*
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * Gromacs Runs On Most of All Computer Systems
 */

#ifndef _statusio_h
#define _statusio_h

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "typedefs.h"
#include "sheader.h"

/*
 * This module handles status file io. All read and write operations from
 * and to a status file should use this functions to be independent of the
 * actual file layout (text versus binary file).
 */
#ifdef __cplusplus
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

#ifdef __cplusplus
}
#endif

#endif	/* _statusio_h */
