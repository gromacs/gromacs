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
 * GRoups of Organic Molecules in ACtion for Science
 */
static char *SRCID_rdklib_h = "$Id$";

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

