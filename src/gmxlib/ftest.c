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
 * GRowing Old MAkes el Chrono Sweat
 */
#include "stdio.h"
#include "stdlib.h"

void func(char *fn,int n,float f[],int buffer)
{
  FILE   *fp;
  size_t size;
  
  if ((fp=fopen(fn,"w")) == NULL) {
    perror(fn);
    exit(1);
  }
  printf("fd:    %d\n",fp->fd);
  printf("flags: %x\n",fp->flags);
  
      
  if (buffer == 0)
    setbuf(fp,NULL);
  fprintf(fp,"This is a recording\n");
  
  size=(f) ? n*sizeof(f[0]) : 0;
  if (size) fwrite(f,1,size,fp);
  
  fprintf(fp,"Your Mileage may vary\n");
  fclose(fp);
}

int main(int argc,char *argv[])
{
#define MAX 1024
  float *fbuf;
  int   i;
  
  fbuf=(float *)malloc(MAX*sizeof(fbuf[0]));
  for(i=0; (i<MAX); i++)
    fbuf[i]=i;
    
  func("ok.dat",MAX,NULL,0);
  func("sucks.dat",MAX,NULL,1);
  
  return 0;
}
