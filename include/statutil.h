/*
 *       $Id$
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0
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
 * GRowing Old MAkes el Chrono Sweat
 */

#ifndef _statutil_h
#define _statutil_h

static char *SRCID_statutil_h = "$Id$";

#ifdef CPLUSPLUS
extern "C" {
#endif

#include <stdio.h>
#include "typedefs.h"
#include "filenm.h"
#include "readinp.h"
#include "wman.h"

typedef int t_first_x(int *status,char *fn,real *t,rvec **x,matrix box);
typedef bool t_next_x(int status,real *t,int natoms,rvec x[],matrix box);

/* I/O function types */

extern void usage(char *prog,char *arg);
/* Error function */

extern char *Program(void);
/* Return the name of the program */
extern char *ShortProgram(void);
/* Id. without leading directory */

/************************************************
 *             Trajectory functions
 ************************************************/


int write_trx(int fnum,int nind,atom_id *ind,t_atoms *atoms,
	      int step,real time,matrix box,rvec x[],rvec *v);
/* write an indexed frame to a TRX file */ 

int close_trx(int fnum);
/* Close a TRX file */

int open_trx(char *outfile,char *filemode);
/* Open a TRX file and return the file number */

extern int check_times(real t);
/* This routine checkes if the read-in time is correct or not;
 * returns -1 if t<tbegin,
 *         0  if tbegin <= t <=tend,
 *         1  if t>tend
 */

extern int read_first_x(int *status,char *fn,
			real *t,rvec **x,matrix box);
/* These routines read first coordinates and box, and allocates 
 * memory for the coordinates, for a trajectory file.
 * The routine returns the number of atoms, or 0 when something is wrong.
 * The integer in status should be passed to calls of read_next_x
 */

extern bool read_next_x(int status,real *t,int natoms,rvec x[],matrix box);
/* Read coordinates and box from a trajectory file. Return TRUE when all well,
 * or FALSE when end of file (or last frame requested by user).
 * status is the integer set in read_first_x.
 */

extern void close_trj(int status);
/* Close trj file as opened with read_first_x */

extern void rewind_trj(int status);
/* Rewind trj file as opened with read_first_x */

extern int read_first_v(int *status,char *fn,real *t,rvec **v,matrix box);
/* Same as above, but for velocities */

extern  bool read_next_v(int status,real *t,int natoms,rvec v[],matrix box);
/* Idem */

extern int read_first_x_v(int *status,char *fn,real *t,rvec **x,rvec **v,matrix box);
/* Same as above, but for coordinates and velocities */
  
extern  bool read_next_x_v(int status,real *t,int natoms,rvec x[],rvec v[],matrix box);
/* Idem */

extern  bool read_next_x_or_v(int status,real *t,int natoms,rvec x[],rvec v[],matrix box);
/* Same as above, but for coordinates *AND/OR* velocities */

extern  bool read_first_x_or_v(int *status,char *fn,real *t,rvec **x,rvec **v,matrix box);
/* Idem */
 
extern bool next_e(FILE *status, real *t, t_energy e[]);
/* Read energy terms from trajectory file */

extern t_topology *read_top(char *fn);
/* Extract a topology data structure from a topology file */

extern void mk_single_top(t_topology *top);
/* Make the topology file single processor ready */

extern char *status_title(FILE *status);
/* Return a title from a topology file */

extern bool bDoView(void);
/* Return TRUE when user requested viewing of the file */

/*****************************************************
 *         Some command line parsing routines 
 *****************************************************/

#define PCA_CAN_VIEW       (1<<5)
#define PCA_CAN_BEGIN      (1<<6)
#define PCA_CAN_END        (1<<7)
#define PCA_CAN_TIME       (PCA_CAN_BEGIN | PCA_CAN_END)
#define PCA_KEEP_ARGS      (1<<8)
#define PCA_SILENT         (1<<9)
#define PCA_NOGET_PARGS    (1<<10)
#define PCA_NOEXIT_ON_ARGS (1<<11)
#define PCA_QUIET          (1<<12)
#define PCA_SET_NPRI       (1<<13)

extern int iscan(int argc,char *argv[],int *i);
/* Scan an int from the argument at *i. If the argument length
 * is > 2, the int is assumed to be in the remainder of the arg,
 * eg: -p32, else the int is assumed to be in the next argument
 * eg: -p 32. If neither is the case the routine exits with an error,
 * otherwise it returns the value found. If the value is in the next
 * argument *i is incremented. You typically would want to pass
 * a loop variable to this routine.
 */

extern double dscan(int argc,char *argv[],int *i);
/* Routine similar to the above, but working on doubles. */

extern char *sscan(int argc,char *argv[],int *i);
/* Routine similar to the above, but working on strings. The pointer
 * returned is a pointer to the argv field.
 */

extern void vscan(int argc,char *argv[],int *i,rvec *vec);
/* Routine similar to the above, but working on rvecs. */

#ifdef HAVE_MOTIF
extern void gmx_gui(int *argc,char *argv[],
		    int nfile,t_filenm fnm[],int npargs,t_pargs pa[],
		    int ndesc,char *desc[],int nbugs,char *bugs[]);
  /* This function plops up a Motif dialog box in which the command-line options
   * can be changed.
   */
#endif

extern void parse_common_args(int *argc,char *argv[],ulong Flags,bool bNice,
			      int nfile,t_filenm fnm[],int npargs,t_pargs pa[],
			      int ndesc,char *desc[],int nbugs,char *bugs[]);
/* Get arguments from the arg-list. The arguments extracted
 * are removed from the list. If manual is NULL a default message is displayed
 * when errors are encountered. The Flags argument, when non-0 enables
 * some input checks. Using this routine also means that the arguments
 * -b and -e will be used for begin and end time, whether this is 
 * appropriate or not!
 */

#ifdef CPLUSPLUS
}
#endif

#endif
