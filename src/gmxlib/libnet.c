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
static char *SRCID_libnet_c = "$Id$";

#include "main.h"
#include "fatal.h"
#include "network.h"
#include "typedefs.h"
#include "ns.h"

/* Network routines */

#define MYFATAL(str) fatal_error(0,"Routine %s called in %s, %d",str,__FILE__,__LINE__)

void gmx_tx(int pid,void *bufptr,int bufsize)	
{ 
  MYFATAL("gmx_tx"); 
}

void gmx_tx_wait(int pid)				
{
  MYFATAL("gmx_tx_wait");
}

void gmx_txs(int pid,void *bufptr,int bufsize)	
{
  MYFATAL("gmx_txs");
}

void gmx_rx(int pid,void *bufptr,int bufsize)	
{
  MYFATAL("gmx_rx");
}

void gmx_rx_wait(int pid)				
{
  MYFATAL("gmx_rx_wait");
}

void gmx_rxs(int pid,void *bufptr,int bufsize)	
{
  MYFATAL("gmx_rxs");
}

int gmx_cpu_num(void)
{
  return 1;
}

int gmx_cpu_id(void)
{
  return 0;
}

int  get_idle_rec(void)  					
{ 
  return 0; 
}

int  get_idle_send(void) 					
{ 
  return 0; 
}

void gmx_left_right(int nprocs,int pid,int *left,int *right)
{
  *left=0;
  *right=0;
}
