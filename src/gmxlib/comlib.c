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
 * Great Red Owns Many ACres of Sand 
 */
static char *SRCID_comlib_c = "$Id$";

#define DUALC_PORT	0xb0001000	/* Communication area        */
#define DUAL_PORT	0xb0000000	/* Data area                 */
#define RX_CMD		8*0x18		/* Command receive           */
#define TX_CMD		8*0x19		/* Command transmitt         */
#define RX_DAT		8*0x1a		/* Data receive single byte  */
#define TX_DAT		8*0x1b		/* Data transmit single byte */

#include "comlib.h"

volatile static unsigned char *dualcp=(unsigned char *)DUALC_PORT;

void put_serverbyte(unsigned char data)
{
  dualcp[TX_DAT]=data;
  dualcp[TX_CMD]=0x01;
  while(dualcp[TX_CMD]&0x01);
}

unsigned char get_serverbyte()
{
  unsigned char data;
  
  while((dualcp[RX_CMD]&0x01)==0x00);
  data=dualcp[RX_DAT];
  dualcp[RX_CMD]=0x00;
  return data;
}

void get_serverdata(void *data,int size)
{
  char *p;

  for(p=data; size>0; size--) *(p++)=get_serverbyte();
}

void put_serverdata(void *data,int size)
{
  char *p;

  for(p=data; size>0; size--) put_serverbyte(*(p++));
}
