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
