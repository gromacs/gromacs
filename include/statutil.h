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

extern char *Program(void);
/* Return the name of the program */
extern char *ShortProgram(void);
/* Id. without leading directory */

/************************************************
 *             Trajectory functions
 ************************************************/

extern int prec2ndec(real prec);
/* Convert precision in 1/(nm) to number of decimal places */

extern void clear_trxframe(t_trxframe *fr,bool bFirst);
/* Set all content booleans to FALSE.
 * When bFirst = TRUE, set natoms=-1, all pointers to NULL
 *                     and all data to zero.
 */

extern int nframes_read(void);
/* Returns the number of frames read from the trajectory */

int write_trxframe_indexed(int status,t_trxframe *fr,int nind,atom_id *ind);
/* Write an indexed frame to a TRX file, see write_trxframe */

int write_trxframe(int status,t_trxframe *fr);
/* Write a frame to a TRX file. 
 * Only entries for which the boolean is TRUE will be written,
 * except for step, time, lambda and/or box, which may not be
 * omitted for certain trajectory formats.
 * The precision for .xtc and .gro is fr->prec, when fr->bPrec=FALSE,
 * the precision is set to 1000.
 */

int write_trx(int status,int nind,atom_id *ind,t_atoms *atoms,
	      int step,real time,matrix box,rvec x[],rvec *v);
/* Write an indexed frame to a TRX file.
 * v can be NULL. 
 * atoms can be NULL for file types which don't need atom names.
 */ 

void close_trx(int status);
/* Close trj file as opened with read_first_x, read_frist_frame
 * or open_trx. Identical to close_trj.
 */

int open_trx(char *outfile,char *filemode);
/* Open a TRX file and return the file number */

extern bool bRmod(double a,double b);
/* Returns TRUE when a MOD b = 0, using a margin which is slightly
 * larger than the float/double precision.
 */

extern int check_times(real t,real t0);
/* This routine checkes if the read-in time is correct or not;
 * returns -1 if t<tbegin or t MOD dt = t0,
 *         0  if tbegin <= t <=tend,
 *         1  if t>tend
 */

/* For trxframe.flags, used in trxframe read routines.
 * When a READ flag is set, the field will be read when present,
 * but a frame might be returned which does not contain the field.
 * When a NEED flag is set, frames not containing the field will be skipped.
 */
#define TRX_READ_X    (1<<0)
#define TRX_NEED_X    (1<<1)
#define TRX_READ_V    (1<<2)
#define TRX_NEED_V    (1<<3)
#define TRX_READ_F    (1<<4)
#define TRX_NEED_F    (1<<5)
/* Useful for reading natoms from a trajectory without skipping */
#define TRX_DONT_SKIP (1<<6)

/* For trxframe.not_ok */
#define HEADER_NOT_OK (1<<0)
#define DATA_NOT_OK   (1<<1)
#define FRAME_NOT_OK  (HEADER_NOT_OK | DATA_NOT_OK)

extern bool read_first_frame(int *status,char *fn,t_trxframe *fr,int flags);
  /* Read the first frame which is in accordance with flags, which are
   * defined further up in this file. 
   * Returns natoms when succeeded, 0 otherwise.
   * Memory will be allocated for flagged entries.
   * The flags are copied to fr for subsequent calls to read_next_frame.
   * Returns TRUE when succeeded, FALSE otherwise.
   */

extern bool read_next_frame(int status,t_trxframe *fr);
  /* Reads the next frame which is in accordance with fr->flags.
   * Returns TRUE when succeeded, FALSE otherwise.
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
/* Close trj file as opened with read_first_x, read_frist_frame
 * or open_trx. Identical to close_trx.
 */

extern void rewind_trj(int status);
/* Rewind trj file as opened with read_first_x */

extern t_topology *read_top(char *fn);
/* Extract a topology data structure from a topology file */

extern void mk_single_top(t_topology *top);
/* Make the topology file single processor ready */

extern bool bDoView(void);
/* Return TRUE when user requested viewing of the file */

/*****************************************************
 *         Some command line parsing routines 
 *****************************************************/

#define PCA_CAN_VIEW       (1<<5)
/* add option -w to view output files (must be implemented in program) */
#define PCA_CAN_BEGIN      (1<<6)
#define PCA_CAN_END        (1<<7)
#define PCA_CAN_DT         (1<<14)
#define PCA_CAN_TIME       (PCA_CAN_BEGIN | PCA_CAN_END | PCA_CAN_DT)
/* adds options -b and -e for begin and end time for reading trajectories */
#define PCA_KEEP_ARGS      (1<<8)
/* keep parsed args in argv (doesn't make sense without NOEXIT_ON_ARGS) */
#define PCA_SILENT         (1<<9)
/* don't print options by default */
#define PCA_CAN_SET_DEFFNM (1<<10)
/* does something for non-master mdrun processes */
#define PCA_NOEXIT_ON_ARGS (1<<11)
/* no fatal_error when invalid options are encountered */
#define PCA_QUIET          (1<<12)
/* does something for non-master mdrun processes */
#define PCA_SET_NPRI       (1<<13)
/* set weightless prioriy by default */

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

extern void parse_common_args(int *argc,char *argv[],unsigned long Flags,bool bNice,
			      int nfile,t_filenm fnm[],int npargs,t_pargs *pa,
			      int ndesc,char **desc,int nbugs,char **bugs);
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
