/*
 * $Id$
 * 
 *       This source code is part of
 * 
 *        G   R   O   M   A   C   S
 * 
 * GROningen MAchine for Chemical Simulations
 * 
 *               VERSION 2.0
 * 
 * Copyright (c) 1991-1999
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 * 
 * Also check out our WWW page:
 * http://md.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 * 
 * And Hey:
 * Green Red Orange Magenta Azure Cyan Skyblue
 */
static char *SRCID_sched_h = "$Id$";

#ifndef SCHED_H
#define SCHED_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

extern void free_comm_schedule(int **sched, int npes);
extern void empty_comm_schedule(int **sched, int npes);
extern int **make_comm_schedule(int npes);
extern void fill_comm_schedule(int **sched, int npes);
extern int check_comm_schedule(int **sched, int npes);
extern void invert_comm_schedule(int **sched, int npes);
extern void sort_comm_schedule(int **sched, int npes, int sort_pe);
extern void print_comm_schedule(int **sched, int npes);

#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* SCHED_H */
