/* This code is part of the tng compression routines.
 *
 * Written by Daniel Spangberg
 * Copyright (c) 2010, 2013, The GROMACS development team.
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

void Ptngc_comp_make_dict_hist(unsigned int *vals, int nvals,
                         unsigned int *dict, int *ndict,
                         unsigned int *hist)
{
  int i;
  int j=0;
  for (i=0; i<0x20004; i++)
    hist[i]=0;
  for (i=0; i<0x20004; i++)
    dict[i]=i;
  for (i=0; i<nvals; i++)
    hist[vals[i]]++;
  for (i=0; i<0x20004; i++)
    if (hist[i]!=0)
      {
        hist[j]=hist[i];
        dict[j]=dict[i];
        j++;
      }
  *ndict=j;
}
