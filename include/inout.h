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
 * Gyas ROwers Mature At Cryogenic Speed
 */

#ifndef _inout_h
#define _inout_h

static char *SRCID_inout_h = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) inout.h 1.6 11/23/92"
#endif /* HAVE_IDENT */

/***************************************************************************/
/*                                                                         */
/* Memory mapped io for the SPC-i860 by A. Sijbers, (c) 1991               */
/*                                                                         */
/* General description:                                                    */
/*                                                                         */
/* The i860 io buffers (both dual ported rams, idt7130 and hd63310) are    */
/* organised as bytes located at the lowest byte address in a 64 bit word. */
/* So to address consecutive bytes of the io buffers, it is necessary to   */
/* increment the i860 byte address by 8.                                   */
/* This module implements a general method of accessing the io buffers by  */
/* specifying the base address of the buffer and an offset. The base       */
/* address is the physical address of the buffer, the offset is byte count */
/* in the buffers local address space (to address the first byte of an io  */
/* buffer, specify 0, for the second byte 1 etc.).                         */
/* Although it is possible to use any address and buffer size combination  */
/* for the put_io_buf and get_io_buf routines, it is strongly recommended  */
/* to use only long word aligned addresses and sizes which are multiple    */
/* of 4 for maximum speed. The value of offset has no influence on the     */
/* transfer rate. For details, see the implementation.                     */
/*                                                                         */
/***************************************************************************/

#define HD63310		0xb0000000	/* HD63310 dual port ram base */
#define IDT7130		0xd0000000	/* IDT7130 dual port ram base */

extern void poke_io_buf(int iobase,int offset,int byte);
     /*
      * Puts byte at io buffer, pointed to by iobase, location offset.
      */

extern int peek_io_buf(int iobase,int offset);
     /*
      * Returns the byte at offset from io buffer pointed to by iobase.
      */

extern void put_io_buf(int iobase,int offset,void *buf,int bufsize);
     /*
      * Puts bufsize bytes of buffer, pointed to by buf, into io buffer, 
      * pointed to by iobase, starting at offset.
      */

extern void get_io_buf(int iobase,int offset,void *buf,int bufsize);
     /*
      * Puts bufsize bytes from io buffer, pointed to by iobase, location 
      * offset into the buffer pointed to by buf.
      */

#endif	/* _inout_h */
