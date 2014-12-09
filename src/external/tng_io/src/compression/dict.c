/* This code is part of the tng compression routines.
 *
 * Written by Daniel Spangberg and Magnus Lundborg
 * Copyright (c) 2010, 2013-2014 The GROMACS development team.
 *
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the Revised BSD License.
 */


#include <string.h>
#include "../../include/compression/dict.h"

void Ptngc_comp_canonical_dict(unsigned int *dict, int *ndict)
{
  int i;
  for (i=0; i<0x20004; i++)
    dict[i]=i;

  *ndict=0x20004;
}

void Ptngc_comp_make_dict_hist(unsigned int *vals, const int nvals,
                         unsigned int *dict, int *ndict,
                         unsigned int *hist)
{
  int i;
  int j=0;

  memset(hist, 0, sizeof(unsigned int)*0x20004);

  for (i=0; i<nvals; i++)
    hist[vals[i]]++;
  for (i=0; i<0x20004; i++)
    if (hist[i]!=0)
      {
        hist[j]=hist[i];
        dict[j]=i;
        j++;
        if(j==nvals)
          break;
      }
  *ndict=j;
}
