/*
 *       @(#) copyrgt.c 1.12 9/30/97
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0b
 * 
 * Copyright (c) 1990-1997,
 * BIOSON Research Institute, Dept. of Biophysical Chemistry,
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

#ifndef __AMBADMA__
#define __AMBADMA__

/*------------------------------------------------------------------------
 *  edition history
 *------------------------------------------------------------------------
 *
 *  version  date    author  description
 *  ----------------------------------------------------------------------
 *  1.1      28Dec93 mkr     creation
 */

/*------------------------------------------------------------------------
 *  system includes
 *------------------------------------------------------------------------
 */

/*------------------------------------------------------------------------
 *  application includes
 *------------------------------------------------------------------------
 */
#include	"general.h"

/*------------------------------------------------------------------------
 *  library variables
 *------------------------------------------------------------------------
 */

/*------------------------------------------------------------------------
 *  global variables
 *------------------------------------------------------------------------
 */

/*------------------------------------------------------------------------
 *  constant definitions
 *------------------------------------------------------------------------
 */
#define	NCHANS	2					/* # of dma-channels for 1 i860-CPU	*/
#define	NLEDS	4					/* # of user LED's for 1 i860-CPU	*/

#define OFF		0					/* switch a LED off					*/
#define ON		1					/* switch a LED on					*/

#define	NOT_CONF	0				/* channel is not configured		*/
#define	CONFIGURED	1				/* channel is configured			*/

#define	AMBA_DMAREG_BASE	0xfd000000	/* base of Amba DMA registers	*/
#define	DMA_MAX_XFER	(1 << 18)	/* maximum size of a DMA transfer	*/	
#define	DMA_MIN_XFER	(3 * 8)		/* minimum size of a DMA transfer	*/	
#define	DMA_DATA_SIZE	8		/* data size of a DMA transfer	*/	

#define DMA_LOW_SHIFT	3			/* shift for low word of DMA start	*/
#define DMA_HIGH_SHIFT	19			/* shift for high word of DMA start	*/
#define	DMA_LOW_MASK	0xffff		/* mask low word of DMA start addr	*/
#define	DMA_HIGH_MASK	0x3f		/* mask high word of DMA start addr	*/

#define	DMA_CR_INIT_RX	0x02		/* initialize the CR register		*/
#define	DMA_CR_INIT_TX	0x00		/* initialize the CR register		*/
#define	DMA_CR_START	0x40		/* start a DMA cycle				*/
#define	DMA_CR_DISABLE	0x01		/* disable (unconfigure) DMA channel*/
#define	DMA_STS_DMAIP	0x20		/* DMA in progress bit				*/
#define	DMA_STS_EODMA	0x02		/* End Of DMA bit					*/

#define DMA_VALID_MASK	0x07		/* mask of bit A2 - A0				*/

#define	DMA_STS			0x1c		/* mask out status of DMA transfer	*/
#define	DMA_STS_REM		0x0c		/* mask out remote status			*/
#define	DMA_REM_ERROR	0x0c		/* remote status is error			*/
#define	DMA_REM_RX		0x08		/* remote status is input			*/
#define	DMA_REM_TX		0x04		/* remote status is output			*/
#define	DMA_REM_CLR		0x00		/* remote status is not configured	*/

#define	DMA_UNDERRUN	0x10		/* an underrun error has occured	*/
#define	DMA_OVERRUN		0x10		/* an overrun error has occured		*/
#define	DMA_RC_SYNC		0x08		/* a receiver sync error has occured*/
#define	DMA_TX_SYNC		0x04		/* a transm. sync error has occured	*/

#define DMA_STAT_PERR	0x4000		/* a parity error occured			*/
#define DMA_STAT_TIMOUT	0x8000		/* a access time-out occured		*/
#define DMA_TIMOUT_WNR	0x2000		/* determine operation (read/write)	*/

#define	TIMOUT_READ		0			/* it was a read operation			*/
#define	TIMOUT_WRITE	1			/* it was a write operation			*/

/*------------------------------------------------------------------------
 *  structure definitions
 *------------------------------------------------------------------------
 */
typedef struct
	{
	volatile WORD	dummy[3];
	volatile WORD	lAddr0;			/* WRITE: low DMA start-addr. reg. 0*/
									/* READ: low parity error reg.		*/
	volatile WORD	dummy0[3];
	volatile WORD	hAddr0;			/* WRITE: high DMA start-addr. reg.0*/
									/* READ: high parity error reg.		*/
	volatile WORD	dummy1[3];
	volatile short int	cc0;			/* cycle count register	0			*/
	volatile WORD	dummy2[3];
	volatile WORD	cr0;			/* control register 0				*/
	volatile WORD	dummy3[3];
	volatile WORD	lAddr1;			/* WRITE: low DMA start-addr. reg. 1*/
									/* READ: low time-out error reg.	*/
	volatile WORD	dummy4[3];
	volatile WORD	hAddr1;			/* WRITE: high DMA start-addr. reg.1*/
									/* READ: high time-out error reg.	*/
	volatile WORD	dummy5[3];
	volatile short int	cc1;			/* cycle count register	1			*/
	volatile WORD	dummy6[3];
	volatile WORD	cr1;			/* control register 1				*/
	volatile WORD	dummy7[3];
	volatile WORD	led01;			/* LED register for LED 0 and 1		*/
	volatile WORD	dummy8[15];
	volatile WORD	led23;			/* LED register for LED 2 and 3		*/
	} AMBADMAREG;

#endif  /* __AMBADMA__ */
