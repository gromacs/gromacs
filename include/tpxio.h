#ifndef _tpxio_h
#define _tpxio_h

#ifdef CPLUSPLUS 
extern "C" {
#endif

  /**************************************************************
   *
   * The routines in the corresponding c-file tpxio.c
   * are based on the lower level routines in gmxfio.c
   * The integer file pointer returned from open_tpx
   * can also be used with the routines in gmxfio.h
   *
   **************************************************************/
#include "typedefs.h"
#include "xdrf.h"

typedef struct
{
  int	bIr;		/* Non zero if input_rec is present		*/
  int	bBox;		/* Non zero if a box is present			*/
  int	bTop;		/* Non zero if a topology is present		*/
  int	bX;		/* Non zero if coordinates are present		*/
  int	bV;		/* Non zero if velocities are present		*/
  int	bF;		/* Non zero if forces are present		*/

  int	natoms;		/* The total number of atoms			*/
  int	step;		/* Current step number				*/
  real	t;		/* Current time					*/
  real	lambda;		/* Current value of lambda			*/
} t_tpxheader;

/* 
 * These routines handle reading and writing of preprocessed
 * topology files in any of the following formats:
 * TPR : topology in XDR format, portable accross platforms
 * TPB : binary topology, not portable accross platforms
 * TPA : ascii topology (possibbly huge)
 * TRR : trajectory in XDR format (non compressed)
 * TRJ : trajectory in binary format
 *
 * Files are written in the precision with which the source are compiled,
 * but double and single precision can be read by either.
 */

extern int open_tpx(char *fn,char *mode);
/* Return an integer corresponding to the file you have just opened */
  
extern void close_tpx(int fp);
/*  Close the file corresponding to fp */
  
extern void read_tpxheader(char *fn,t_tpxheader *tpx);
/* Read the header from a tpx file and then close it again */

extern void write_tpx(char *fn,int step,real t,real lambda,
		      t_inputrec *ir,rvec *box,int natoms,
		      rvec *x,rvec *v,rvec *f,t_topology *top);
/* Write a file, and close it again. 
 * If fn == NULL, an efTPA file will be written to stdout (which
 * will not be closed afterwards)
 */

extern void read_tpx(char *fn,int *step,real *t,real *lambda,
		     t_inputrec *ir,rvec *box,int *natoms,
		     rvec *x,rvec *v,rvec *f,t_topology *top);
/* Read a file, and close it again. 
 * If fn == NULL, an efTPA file will be read from stdin (which
 * will not be closed afterwards)
 */

extern void fwrite_tpx(int fp,int step,real t,real lambda,
		       t_inputrec *ir,rvec *box,int natoms,
		       rvec *x,rvec *v,rvec *f,t_topology *top);
/* Write a file, and do not close it */

extern void fread_tpx(int fp,int *step,real *t,real *lambda,
		      t_inputrec *ir,rvec *box,int *natoms,
		      rvec *x,rvec *v,rvec *f,t_topology *top);
/* Read a file, and do not close it */

extern bool fn2bTPX(char *file);
/* return if *file is one of the TPX file types */ 

extern bool read_tps_conf(char *infile,char *title,t_topology *top,
			  rvec **x,rvec **v,matrix box,bool bMass);
/* Read title, top.atoms, x, v (if not NULL) and box from an STX file,
 * memory for atoms, x and v will be allocated.  
 * Return TRUE if a topology was read. 
 * If infile is a TPX file read the whole top,
 * else if bMass=TRUE, read the masses into top.atoms from the mass database.
 */

#ifdef CPLUSPLUS
}
#endif

#endif
