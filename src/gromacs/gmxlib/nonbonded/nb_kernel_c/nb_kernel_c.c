/*
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009, The GROMACS Development Team
 *
 * Gromacs is a library for molecular simulation and trajectory analysis,
 * written by Erik Lindahl, David van der Spoel, Berk Hess, and others - for
 * a full list of developers and information, check out http://www.gromacs.org
 *
 * This program is free software; you can redistribute it and/or modify it under 
 * the terms of the GNU Lesser General Public License as published by the Free 
 * Software Foundation; either version 2 of the License, or (at your option) any 
 * later version.
 * As a special exception, you may use this file as part of a free software
 * library without restriction.  Specifically, if other files instantiate
 * templates or use macros or inline functions from this file, or you compile
 * this file and link it with other files to produce an executable, this
 * file does not by itself cause the resulting executable to be covered by
 * the GNU Lesser General Public License.  
 *
 * In plain-speak: do not worry about classes/macros/templates either - only
 * changes to the library have to be LGPL, not an application linking with it.
 *
 * To help fund GROMACS development, we humbly ask that you cite
 * the papers people have written on it - you can find them on the website!
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>

#include "types/nrnb.h"
#include "nb_kernel_c.h"
#include "../nb_kerneltype.h"


/* Include standard kernel headers in local directory */
#include "nb_kernel010.h"
#include "nb_kernel020.h"
#include "nb_kernel030.h"
#include "nb_kernel100.h"
#include "nb_kernel101.h"
#include "nb_kernel102.h"
#include "nb_kernel103.h"
#include "nb_kernel104.h"
#include "nb_kernel110.h"
#include "nb_kernel111.h"
#include "nb_kernel112.h"
#include "nb_kernel113.h"
#include "nb_kernel114.h"
#include "nb_kernel120.h"
#include "nb_kernel121.h"
#include "nb_kernel122.h"
#include "nb_kernel123.h"
#include "nb_kernel124.h"
#include "nb_kernel130.h"
#include "nb_kernel131.h"
#include "nb_kernel132.h"
#include "nb_kernel133.h"
#include "nb_kernel134.h"
#include "nb_kernel200.h"
#include "nb_kernel201.h"
#include "nb_kernel202.h"
#include "nb_kernel203.h"
#include "nb_kernel204.h"
#include "nb_kernel210.h"
#include "nb_kernel211.h"
#include "nb_kernel212.h"
#include "nb_kernel213.h"
#include "nb_kernel214.h"
#include "nb_kernel220.h"
#include "nb_kernel221.h"
#include "nb_kernel222.h"
#include "nb_kernel223.h"
#include "nb_kernel224.h"
#include "nb_kernel230.h"
#include "nb_kernel231.h"
#include "nb_kernel232.h"
#include "nb_kernel233.h"
#include "nb_kernel234.h"
#include "nb_kernel300.h"
#include "nb_kernel301.h"
#include "nb_kernel302.h"
#include "nb_kernel303.h"
#include "nb_kernel304.h"
#include "nb_kernel310.h"
#include "nb_kernel311.h"
#include "nb_kernel312.h"
#include "nb_kernel313.h"
#include "nb_kernel314.h"
#include "nb_kernel320.h"
#include "nb_kernel321.h"
#include "nb_kernel322.h"
#include "nb_kernel323.h"
#include "nb_kernel324.h"
#include "nb_kernel330.h"
#include "nb_kernel331.h"
#include "nb_kernel332.h"
#include "nb_kernel333.h"
#include "nb_kernel334.h"
#include "nb_kernel400.h"
#include "nb_kernel410.h"
#include "nb_kernel420.h"
#include "nb_kernel430.h"


static nb_kernel_t *
kernellist[eNR_NBKERNEL_NR] = 
{
    nb_kernel010,
    nb_kernel020,
    nb_kernel030,
    nb_kernel100,
    nb_kernel101,
    nb_kernel102,
    nb_kernel103,
    nb_kernel104,
    nb_kernel110,
    nb_kernel111,
    nb_kernel112,
    nb_kernel113,
    nb_kernel114,
    nb_kernel120,
    nb_kernel121,
    nb_kernel122,
    nb_kernel123,
    nb_kernel124,
    nb_kernel130,
    nb_kernel131,
    nb_kernel132,
    nb_kernel133,
    nb_kernel134,
    nb_kernel200,
    nb_kernel201,
    nb_kernel202,
    nb_kernel203,
    nb_kernel204,
    nb_kernel210,
    nb_kernel211,
    nb_kernel212,
    nb_kernel213,
    nb_kernel214,
    nb_kernel220,
    nb_kernel221,
    nb_kernel222,
    nb_kernel223,
    nb_kernel224,
    nb_kernel230,
    nb_kernel231,
    nb_kernel232,
    nb_kernel233,
    nb_kernel234,
    nb_kernel300,
    nb_kernel301,
    nb_kernel302,
    nb_kernel303,
    nb_kernel304,
    nb_kernel310,
    nb_kernel311,
    nb_kernel312,
    nb_kernel313,
    nb_kernel314,
    nb_kernel320,
    nb_kernel321,
    nb_kernel322,
    nb_kernel323,
    nb_kernel324,
    nb_kernel330,
    nb_kernel331,
    nb_kernel332,
    nb_kernel333,
    nb_kernel334,
    nb_kernel400,
    nb_kernel410,
    nb_kernel430,
	nb_kernel010nf,
    nb_kernel020nf,
    nb_kernel030nf,
    nb_kernel100nf,
    nb_kernel101nf,
    nb_kernel102nf,
    nb_kernel103nf,
    nb_kernel104nf,
    nb_kernel110nf,
    nb_kernel111nf,
    nb_kernel112nf,
    nb_kernel113nf,
    nb_kernel114nf,
    nb_kernel120nf,
    nb_kernel121nf,
    nb_kernel122nf,
    nb_kernel123nf,
    nb_kernel124nf,
    nb_kernel130nf,
    nb_kernel131nf,
    nb_kernel132nf,
    nb_kernel133nf,
    nb_kernel134nf,
    nb_kernel200nf,
    nb_kernel201nf,
    nb_kernel202nf,
    nb_kernel203nf,
    nb_kernel204nf,
    nb_kernel210nf,
    nb_kernel211nf,
    nb_kernel212nf,
    nb_kernel213nf,
    nb_kernel214nf,
    nb_kernel220nf,
    nb_kernel221nf,
    nb_kernel222nf,
    nb_kernel223nf,
    nb_kernel224nf,
    nb_kernel230nf,
    nb_kernel231nf,
    nb_kernel232nf,
    nb_kernel233nf,
    nb_kernel234nf,
    nb_kernel300nf,
    nb_kernel301nf,
    nb_kernel302nf,
    nb_kernel303nf,
    nb_kernel304nf,
    nb_kernel310nf,
    nb_kernel311nf,
    nb_kernel312nf,
    nb_kernel313nf,
    nb_kernel314nf,
    nb_kernel320nf,
    nb_kernel321nf,
    nb_kernel322nf,
    nb_kernel323nf,
    nb_kernel324nf,
    nb_kernel330nf,
    nb_kernel331nf,
    nb_kernel332nf,
    nb_kernel333nf,
    nb_kernel334nf,
    nb_kernel400nf,
    nb_kernel410nf,
    nb_kernel430nf,
};


void
nb_kernel_setup(FILE *log,nb_kernel_t **list)
{
  int i;
  nb_kernel_t *p;

  if(NULL != log)
    fprintf(log,"Configuring standard C nonbonded kernels...\n");

  for(i=0;i<eNR_NBKERNEL_NR;i++)
  {
    p = kernellist[i];
    if(p!=NULL)
      list[i] = p; 
  }
}    
