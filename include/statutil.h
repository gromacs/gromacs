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

#ifndef _statutil_h
#define _statutil_h

#include <stdio.h>
#include "typedefs.h"
#include "filenm.h"
#include "readinp.h"
#include "wman.h"
#include "pdbio.h"
#include "oenv.h"
#include "gmxfio.h"


#ifdef __cplusplus
extern "C" {
#endif
#if 0 /* avoid screwing up indentation */
}
#endif


/* The code below is to facilitate controlled begin and end of
 trajectory reading. Corresponding routines in
 src/gmxlib/tcontrol.c
 */
enum { TBEGIN, TEND, TDELTA, TNR };

gmx_bool bTimeSet(int tcontrol);

real rTimeValue(int tcontrol); 

void setTimeValue(int tcontrol,real value);

/* End trajectory time control */

/* a dedicated status type contains fp, etc. */
typedef struct t_trxstatus t_trxstatus;
  
typedef int t_first_x(t_trxstatus **status,const char *fn,real *t,rvec **x,
                      matrix box);

typedef gmx_bool t_next_x(t_trxstatus *status,real *t,int natoms,rvec x[],
                      matrix box);

/* I/O function types */

  
/* LEGACY FUNCTIONS 

   The program names, command lines, etc. are now also set in the output_env
   structure. That is now the preferred location, but the functions here
   are still available as legacy functions. Because they all act on inherently
   global informaion, their existence in a multi-threaded environment is not
   a real problem. */
    
/* Return the name of the program */
const char *command_line(void);
void set_command_line(int argc, char *argv[]);

/* set the program name to the provided string, but note
 * that it must be a real file - we determine the library
 * directory from its location!
 */    
const char *Program(void);
void set_program_name(const char *argvzero);
/* Id. without leading directory */
const char *ShortProgram(void);

/************************************************
 *             Trajectory functions
 ************************************************/

int prec2ndec(real prec);
/* Convert precision in 1/(nm) to number of decimal places */

void clear_trxframe(t_trxframe *fr,gmx_bool bFirst);
/* Set all content gmx_booleans to FALSE.
 * When bFirst = TRUE, set natoms=-1, all pointers to NULL
 *                     and all data to zero.
 */

void set_trxframe_ePBC(t_trxframe *fr,int ePBC);
/* Set the type of periodic boundary conditions, ePBC=-1 is not set */

int nframes_read(t_trxstatus *status);
/* Returns the number of frames read from the trajectory */

int write_trxframe_indexed(t_trxstatus *status,t_trxframe *fr,int nind,
                           atom_id *ind, gmx_conect gc);
/* Write an indexed frame to a TRX file, see write_trxframe. gc may be NULL */

int write_trxframe(t_trxstatus *status,t_trxframe *fr,gmx_conect gc);
/* Write a frame to a TRX file. 
 * Only entries for which the gmx_boolean is TRUE will be written,
 * except for step, time, lambda and/or box, which may not be
 * omitted for certain trajectory formats.
 * The precision for .xtc and .gro is fr->prec, when fr->bPrec=FALSE,
 * the precision is set to 1000.
 * gc is important for pdb file writing only and may be NULL.
 */

int write_trx(t_trxstatus *status,int nind,atom_id *ind,t_atoms *atoms,
              int step,real time,matrix box,rvec x[],rvec *v,
              gmx_conect gc);
/* Write an indexed frame to a TRX file.
 * v can be NULL. 
 * atoms can be NULL for file types which don't need atom names.
 */ 

void close_trx(t_trxstatus *status);
/* Close trj file as opened with read_first_x, read_frist_frame
 * or open_trx. Identical to close_trj.
 */

t_trxstatus *open_trx(const char *outfile,const char *filemode);
/* Open a TRX file and return an allocated status pointer */

/* get a fileio from a trxstatus */
t_fileio *trx_get_fileio(t_trxstatus *status);


gmx_bool bRmod_fd(double a, double b, double c,gmx_bool bDouble);
/* Returns TRUE when (a - b) MOD c = 0, using a margin which is slightly
 * larger than the float/double precision.
 */

#ifdef GMX_DOUBLE
#define bRmod(a,b,c) bRmod_fd(a,b,c,TRUE)
#else
#define bRmod(a,b,c) bRmod_fd(a,b,c,FALSE)
#endif

int check_times2(real t,real t0,real tp,real tpp,gmx_bool bDouble);
/* This routine checkes if the read-in time is correct or not;
 * returns -1 if t<tbegin or t MOD dt = t0,
 *          0 if tbegin <= t <=tend+margin,
 *          1 if t>tend
 * where margin is 0.1*min(t-tp,tp-tpp), if this positive, 0 otherwise.
 * tp and tpp should be the time of the previous frame and the one before.
 * The mod is done with single or double precision accuracy depending
 * on the value of bDouble.
 */

int check_times(real t);
/* This routine checkes if the read-in time is correct or not;
 * returns -1 if t<tbegin,
 *          0 if tbegin <= t <=tend,
 *          1 if t>tend
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

int read_first_frame(const output_env_t oenv,t_trxstatus **status,
                            const char *fn, t_trxframe *fr,int flags);
  /* Read the first frame which is in accordance with flags, which are
   * defined further up in this file. 
   * Returns natoms when succeeded, 0 otherwise.
   * Memory will be allocated for flagged entries.
   * The flags are copied to fr for subsequent calls to read_next_frame.
   * Returns TRUE when succeeded, FALSE otherwise.
   */

gmx_bool read_next_frame(const output_env_t oenv,t_trxstatus *status,
                            t_trxframe *fr);
  /* Reads the next frame which is in accordance with fr->flags.
   * Returns TRUE when succeeded, FALSE otherwise.
   */

int read_first_x(const output_env_t oenv,t_trxstatus **status,
                        const char *fn, real *t,rvec **x,matrix box);
/* These routines read first coordinates and box, and allocates 
 * memory for the coordinates, for a trajectory file.
 * The routine returns the number of atoms, or 0 when something is wrong.
 * The integer in status should be passed to calls of read_next_x
 */

gmx_bool read_next_x(const output_env_t oenv,t_trxstatus *status,real *t,
                        int natoms, rvec x[],matrix box);
/* Read coordinates and box from a trajectory file. Return TRUE when all well,
 * or FALSE when end of file (or last frame requested by user).
 * status is the integer set in read_first_x.
 */

void close_trj(t_trxstatus *status);
/* Close trj file as opened with read_first_x, read_frist_frame
 * or open_trx. Identical to close_trx.
 */

void rewind_trj(t_trxstatus *status);
/* Rewind trj file as opened with read_first_x */

t_topology *read_top(const char *fn,int *ePBC);
/* Extract a topology data structure from a topology file.
 * If ePBC!=NULL *ePBC gives the pbc type.
 */

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
#define PCA_TIME_UNIT      (1<<15)
/* set time unit for output */
#define PCA_KEEP_ARGS      (1<<8)
/* keep parsed args in argv (doesn't make sense without NOEXIT_ON_ARGS) */
#define PCA_SILENT         (1<<9)
/* don't print options by default */
#define PCA_CAN_SET_DEFFNM (1<<10)
/* does something for non-master mdrun nodes */
#define PCA_NOEXIT_ON_ARGS (1<<11)
/* no fatal_error when invalid options are encountered */
#define PCA_QUIET          (1<<12)
/* does something for non-master mdrun nodes */
#define PCA_BE_NICE        (1<<13)
/* Default to low priority, unless configured with --disable-nice */
#define PCA_NOT_READ_NODE  (1<<16)
/* Is this node not reading: for parallel all nodes but the master */

int iscan(int argc,char *argv[],int *i);
/* Scan an int from the argument at *i. If the argument length
 * is > 2, the int is assumed to be in the remainder of the arg,
 * eg: -p32, else the int is assumed to be in the next argument
 * eg: -p 32. If neither is the case the routine exits with an error,
 * otherwise it returns the value found. If the value is in the next
 * argument *i is incremented. You typically would want to pass
 * a loop variable to this routine.
 */
gmx_large_int_t istepscan(int argc,char *argv[],int *i);
/* Same as above, but for large integer values */

double dscan(int argc,char *argv[],int *i);
/* Routine similar to the above, but working on doubles. */

char *sscan(int argc,char *argv[],int *i);
/* Routine similar to the above, but working on strings. The pointer
 * returned is a pointer to the argv field.
 */

void vscan(int argc,char *argv[],int *i,rvec *vec);
/* Routine similar to the above, but working on rvecs. */

int nenum(const char *const enumc[]);
/* returns ordinal number of selected enum from args 
 * depends on enumc[0] pointing to one of the other elements
 * array must be terminated by a NULL pointer 
 */

void parse_common_args(int *argc,char *argv[],unsigned long Flags,
                              int nfile,t_filenm fnm[],int npargs,t_pargs *pa,
                              int ndesc,const char **desc,
                              int nbugs,const char **bugs, 
                              output_env_t *oenv);
/* Get arguments from the arg-list. The arguments extracted
 * are removed from the list. If manual is NULL a default message is displayed
 * when errors are encountered. The Flags argument, when non-0 enables
 * some input checks. Using this routine also means that the arguments
 * -b and -e will be used for begin and end time, whether this is 
 * appropriate or not!
 */

#ifdef __cplusplus
}
#endif

#endif
