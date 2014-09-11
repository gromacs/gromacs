/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */

#ifndef _main_h
#define _main_h


#include <stdio.h>

#include "gromacs/fileio/filenm.h"
#include "gromacs/legacyheaders/network.h"
#include "gromacs/legacyheaders/typedefs.h"

#ifdef __cplusplus
extern "C" {
#endif

void gmx_log_open(const char *fn, const t_commrec *cr,
                  gmx_bool bAppendFiles, FILE**);
/* Open the log file, if necessary (nprocs > 1) the logfile name is
 * communicated around the ring.
 */

void gmx_log_close(FILE *fp);
/* Close the log file */

void check_multi_int(FILE *log, const gmx_multisim_t *ms,
                     int val, const char *name,
                     gmx_bool bQuiet);
void check_multi_int64(FILE *log, const gmx_multisim_t *ms,
                       gmx_int64_t val, const char *name,
                       gmx_bool bQuiet);
/* Check if val is the same on all processors for a mdrun -multi run
 * The string name is used to print to the log file and in a fatal error
 * if the val's don't match. If bQuiet is true and the check passes,
 * no output is written.
 */

void init_multisystem(t_commrec *cr, int nsim, char **multidirs,
                      int nfile, const t_filenm fnm[], gmx_bool bParFn);
/* Splits the communication into nsim separate simulations
 * and creates a communication structure between the master
 * these simulations.
 * If bParFn is set, the nodeid is appended to the tpx and each output file.
 */

#ifdef __cplusplus
}
#endif

#endif  /* _main_h */
