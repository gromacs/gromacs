/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.1
 * Copyright (c) 1991-2001, University of Groningen, The Netherlands
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
 * Getting the Right Output Means no Artefacts in Calculating Stuff
 */
static char *SRCID_rdklib_h = "$Id$";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* Rolf-Dieter Klein, 900907   1.0   SPC 860     */
/* rev 900917 r 2.01                             */
/* rev 901213 r 2.10                             */
/* rev 910307 r 2.20  new functions              */
/* rev 910308 r 2.21  rdblit                     */

int rdputlink(int channel, int count, void *pointer);
int rdgetlink(int channel, int count, void *pointer);
int rdgramode(int mode);
void rdsetpixel(int x,int y,int color);
void rdhorline(int x1,int x2,int y,char *buffer);
void rdvgapal(int i, int r, int g, int b);
void rdvgaall(void *buffer);
void rdline(int x1,int y1,int x2,int y2,int color);
void rdrect(int x1,int y1,int x2,int y2,int color,int fillmode);
void rdellipse(int x1,int y1,int x2,int y2,int color,int fillmode);
void rddrawmode(int visual,int page);
void rdclip(int x1,int y1,int x2,int y2);
int getch(void);
int rdspcnumber(void);
int rdspccount(void);
int kbhit(void);
void int86(int intno,void *inregs,void *outregs);
void int86x(int intno,void *inregs,void *outregs,void *sregs);
void outp(int port,int byte);
int inp(int port);
void rdbuftopc(void *src,unsigned long destpc,int count);
void rdbuffrompc(unsigned long srcpc,void *dest,int count);
int rdlinda1out(int id,int cnt,void *erg);
int rdlinda1in(int id,int cnt,void *erg);
int rdlinda1rd(int id,int cnt,void *erg);
int rdlindaout(int id,int nocnt,int *cntptr,void **ptrptr);
int rdlindard(int id,int nocnt,int *cntptr,void **ptrptr);
int rdlindain(int id,int nocnt,int *cntptr,void **ptrptr);
void rdtextbin(int mode);    /* mode =0 text, =1 binary for next open */
int rdclrcounter(int index); /* =0 ok, else maxindex */
int rdgetcounter(int index); /* get counter, incremented each time */
int rdclrfield(int maxbits); /* clr field, =0 ok, else =maxbits */
int rdnextfield(int offset); /* get next free field start offset, no free=-1 */
void rdblit(int x1,int x2,int y1,int y2,int xdest,int ydest);  /* grafik genoa */
/* end of header */

