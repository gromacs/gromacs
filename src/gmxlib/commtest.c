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
 * Great Red Oystrich Makes All Chemists Sane
 */
#include <stdlib.h>
#include "typedefs.h"
#include "network.h"
#include "main.h"
#include "smalloc.h"

static void test_sum(t_commrec *cr)
{
  real array[4] = { 1,2,3,4 };
  int  iarray[4] = { 1,2,3,4 };
  
  gmx_sum(4,array,cr);
  if (MASTER(cr))
    fprintf(stderr,"gsum: %8f, %8f, %8f, %8f\n\n",
	    array[0],array[1],array[2],array[3]);
  gmx_sumi(4,iarray,cr);
  if (MASTER(cr))
    fprintf(stderr,"gsumi: %8d, %8d, %8d, %8d\n\n",
	    iarray[0],iarray[1],iarray[2],iarray[3]);
}

static void test_sync(t_commrec *cr,char *s,int left,int right)
{
  int *send;
  int *recv;
  int i;
  
  snew(send,cr->nprocs);
  snew(recv,cr->nprocs);
  if (MASTER(cr)) {
    for(i=0; (i<cr->nprocs); i++) 
      send[i]=i+1;
    gmx_txs(left,array(send,cr->nprocs));
    gmx_rxs(right,array(recv,cr->nprocs));
  }
  else {
    gmx_rxs(right,array(recv,cr->nprocs));
    for(i=0; (i<cr->nprocs); i++)
      send[i]=recv[i]*2;
    gmx_txs(left,array(send,cr->nprocs));
  }
  if (MASTER(cr)) {
    printf("sync %s:",s);
    for(i=0; (i<cr->nprocs); i++)
      printf("%6d",recv[i]);
    printf("\n\n");
  }
}

static void test_async(t_commrec *cr,char *s,int left,int right)
{
  int *send;
  int *recv;
  int i,j;
  
  snew(send,cr->nprocs);
  snew(recv,cr->nprocs);
  if (MASTER(cr)) {
    for(i=0; (i<cr->nprocs); i++) 
      send[i]=i+1;
  }
  
  for(i=0; (i<cr->nprocs); i++) {
    gmx_tx(left,array(send,cr->nprocs));
    gmx_rx(right,array(recv,cr->nprocs));
    gmx_wait(left,right);

    for(j=0; (j<cr->nprocs); j++)
      send[j]=recv[j]*2;
  }
  if (MASTER(cr)) {
    printf("async %s:",s);
    for(i=0; (i<cr->nprocs); i++)
      printf("%6d",recv[i]);
    printf("\n\n");
  }
}

int main(int argc,char *argv[])
{
  t_commrec *cr;
  int       nprocs;
  
  if (argc < 2) {
    fprintf(stderr,"Usage: %s nprocs\n",argv[0]);
    exit(1);
  }
  else
    nprocs=atoi(argv[1]);
    
  cr=init_par(nprocs,"test.log",argv);
  test_sum(cr);
  test_sync(cr,"left to right",cr->left,cr->right);
  test_sync(cr,"right to left",cr->right,cr->left);
  test_async(cr,"left to right",cr->left,cr->right);
  test_async(cr,"right to left",cr->right,cr->left);
  
  return 0;
}
