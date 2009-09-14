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
#include "nb_kernel_power6.h"
#include "../nb_kerneltype.h"


/* Include standard kernel headers in local directory */
#include "nb_kernel010_power6.h"
#include "nb_kernel020_power6.h"
#include "nb_kernel030_power6.h"
#include "nb_kernel100_power6.h"
#include "nb_kernel101_power6.h"
#include "nb_kernel102_power6.h"
#include "nb_kernel103_power6.h"
#include "nb_kernel104_power6.h"
#include "nb_kernel110_power6.h"
#include "nb_kernel111_power6.h"
#include "nb_kernel112_power6.h"
#include "nb_kernel113_power6.h"
#include "nb_kernel114_power6.h"
#include "nb_kernel120_power6.h"
#include "nb_kernel121_power6.h"
#include "nb_kernel122_power6.h"
#include "nb_kernel123_power6.h"
#include "nb_kernel124_power6.h"
#include "nb_kernel130_power6.h"
#include "nb_kernel131_power6.h"
#include "nb_kernel132_power6.h"
#include "nb_kernel133_power6.h"
#include "nb_kernel134_power6.h"
#include "nb_kernel200_power6.h"
#include "nb_kernel201_power6.h"
#include "nb_kernel202_power6.h"
#include "nb_kernel203_power6.h"
#include "nb_kernel204_power6.h"
#include "nb_kernel210_power6.h"
#include "nb_kernel211_power6.h"
#include "nb_kernel212_power6.h"
#include "nb_kernel213_power6.h"
#include "nb_kernel214_power6.h"
#include "nb_kernel220_power6.h"
#include "nb_kernel221_power6.h"
#include "nb_kernel222_power6.h"
#include "nb_kernel223_power6.h"
#include "nb_kernel224_power6.h"
#include "nb_kernel230_power6.h"
#include "nb_kernel231_power6.h"
#include "nb_kernel232_power6.h"
#include "nb_kernel233_power6.h"
#include "nb_kernel234_power6.h"
#include "nb_kernel300_power6.h"
#include "nb_kernel301_power6.h"
#include "nb_kernel302_power6.h"
#include "nb_kernel303_power6.h"
#include "nb_kernel304_power6.h"
#include "nb_kernel310_power6.h"
#include "nb_kernel311_power6.h"
#include "nb_kernel312_power6.h"
#include "nb_kernel313_power6.h"
#include "nb_kernel314_power6.h"
#include "nb_kernel320_power6.h"
#include "nb_kernel321_power6.h"
#include "nb_kernel322_power6.h"
#include "nb_kernel323_power6.h"
#include "nb_kernel324_power6.h"
#include "nb_kernel330_power6.h"
#include "nb_kernel331_power6.h"
#include "nb_kernel332_power6.h"
#include "nb_kernel333_power6.h"
#include "nb_kernel334_power6.h"
#include "nb_kernel400_power6.h"
#include "nb_kernel410_power6.h"
#include "nb_kernel420_power6.h"
#include "nb_kernel430_power6.h"


static nb_kernel_t *
kernellist[eNR_NBKERNEL_NR] = 
{
    nb_kernel010_power6,
    nb_kernel020_power6,
    nb_kernel030_power6,
    nb_kernel100_power6,
    nb_kernel101_power6,
    nb_kernel102_power6,
    nb_kernel103_power6,
    nb_kernel104_power6,
    nb_kernel110_power6,
    nb_kernel111_power6,
    nb_kernel112_power6,
    nb_kernel113_power6,
    nb_kernel114_power6,
    nb_kernel120_power6,
    nb_kernel121_power6,
    nb_kernel122_power6,
    nb_kernel123_power6,
    nb_kernel124_power6,
    nb_kernel130_power6,
    nb_kernel131_power6,
    nb_kernel132_power6,
    nb_kernel133_power6,
    nb_kernel134_power6,
    nb_kernel200_power6,
    nb_kernel201_power6,
    nb_kernel202_power6,
    nb_kernel203_power6,
    nb_kernel204_power6,
    nb_kernel210_power6,
    nb_kernel211_power6,
    nb_kernel212_power6,
    nb_kernel213_power6,
    nb_kernel214_power6,
    nb_kernel220_power6,
    nb_kernel221_power6,
    nb_kernel222_power6,
    nb_kernel223_power6,
    nb_kernel224_power6,
    nb_kernel230_power6,
    nb_kernel231_power6,
    nb_kernel232_power6,
    nb_kernel233_power6,
    nb_kernel234_power6,
    nb_kernel300_power6,
    nb_kernel301_power6,
    nb_kernel302_power6,
    nb_kernel303_power6,
    nb_kernel304_power6,
    nb_kernel310_power6,
    nb_kernel311_power6,
    nb_kernel312_power6,
    nb_kernel313_power6,
    nb_kernel314_power6,
    nb_kernel320_power6,
    nb_kernel321_power6,
    nb_kernel322_power6,
    nb_kernel323_power6,
    nb_kernel324_power6,
    nb_kernel330_power6,
    nb_kernel331_power6,
    nb_kernel332_power6,
    nb_kernel333_power6,
    nb_kernel334_power6,
    nb_kernel400_power6,
    nb_kernel410_power6,
    nb_kernel430_power6,
	nb_kernel010nf_power6,
    nb_kernel020nf_power6,
    nb_kernel030nf_power6,
    nb_kernel100nf_power6,
    nb_kernel101nf_power6,
    nb_kernel102nf_power6,
    nb_kernel103nf_power6,
    nb_kernel104nf_power6,
    nb_kernel110nf_power6,
    nb_kernel111nf_power6,
    nb_kernel112nf_power6,
    nb_kernel113nf_power6,
    nb_kernel114nf_power6,
    nb_kernel120nf_power6,
    nb_kernel121nf_power6,
    nb_kernel122nf_power6,
    nb_kernel123nf_power6,
    nb_kernel124nf_power6,
    nb_kernel130nf_power6,
    nb_kernel131nf_power6,
    nb_kernel132nf_power6,
    nb_kernel133nf_power6,
    nb_kernel134nf_power6,
    nb_kernel200nf_power6,
    nb_kernel201nf_power6,
    nb_kernel202nf_power6,
    nb_kernel203nf_power6,
    nb_kernel204nf_power6,
    nb_kernel210nf_power6,
    nb_kernel211nf_power6,
    nb_kernel212nf_power6,
    nb_kernel213nf_power6,
    nb_kernel214nf_power6,
    nb_kernel220nf_power6,
    nb_kernel221nf_power6,
    nb_kernel222nf_power6,
    nb_kernel223nf_power6,
    nb_kernel224nf_power6,
    nb_kernel230nf_power6,
    nb_kernel231nf_power6,
    nb_kernel232nf_power6,
    nb_kernel233nf_power6,
    nb_kernel234nf_power6,
    nb_kernel300nf_power6,
    nb_kernel301nf_power6,
    nb_kernel302nf_power6,
    nb_kernel303nf_power6,
    nb_kernel304nf_power6,
    nb_kernel310nf_power6,
    nb_kernel311nf_power6,
    nb_kernel312nf_power6,
    nb_kernel313nf_power6,
    nb_kernel314nf_power6,
    nb_kernel320nf_power6,
    nb_kernel321nf_power6,
    nb_kernel322nf_power6,
    nb_kernel323nf_power6,
    nb_kernel324nf_power6,
    nb_kernel330nf_power6,
    nb_kernel331nf_power6,
    nb_kernel332nf_power6,
    nb_kernel333nf_power6,
    nb_kernel334nf_power6,
    nb_kernel400nf_power6,
    nb_kernel410nf_power6,
    nb_kernel430nf_power6,
};


void
nb_kernel_setup_power6(FILE *log,nb_kernel_t **list)
{
  int i;
  nb_kernel_t *p;

  if(log)
  {
#ifdef GMX_DOUBLE
    fprintf(log,"Configuring double precision IBM Power6-specific Fortran kernels...\n");
#else
    fprintf(log,"Configuring single precision IBM Power6-specific Fortran kernels...\n");
#endif
  }

  for(i=0;i<eNR_NBKERNEL_NR;i++)
  {
      p = kernellist[i];
      if(p!=NULL)
	list[i] = p; 
  }
}    
