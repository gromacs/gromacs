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
#include "nbkernel010_power6.h"
#include "nbkernel020_power6.h"
#include "nbkernel030_power6.h"
#include "nbkernel100_power6.h"
#include "nbkernel101_power6.h"
#include "nbkernel102_power6.h"
#include "nbkernel103_power6.h"
#include "nbkernel104_power6.h"
#include "nbkernel110_power6.h"
#include "nbkernel111_power6.h"
#include "nbkernel112_power6.h"
#include "nbkernel113_power6.h"
#include "nbkernel114_power6.h"
#include "nbkernel120_power6.h"
#include "nbkernel121_power6.h"
#include "nbkernel122_power6.h"
#include "nbkernel123_power6.h"
#include "nbkernel124_power6.h"
#include "nbkernel130_power6.h"
#include "nbkernel131_power6.h"
#include "nbkernel132_power6.h"
#include "nbkernel133_power6.h"
#include "nbkernel134_power6.h"
#include "nbkernel200_power6.h"
#include "nbkernel201_power6.h"
#include "nbkernel202_power6.h"
#include "nbkernel203_power6.h"
#include "nbkernel204_power6.h"
#include "nbkernel210_power6.h"
#include "nbkernel211_power6.h"
#include "nbkernel212_power6.h"
#include "nbkernel213_power6.h"
#include "nbkernel214_power6.h"
#include "nbkernel220_power6.h"
#include "nbkernel221_power6.h"
#include "nbkernel222_power6.h"
#include "nbkernel223_power6.h"
#include "nbkernel224_power6.h"
#include "nbkernel230_power6.h"
#include "nbkernel231_power6.h"
#include "nbkernel232_power6.h"
#include "nbkernel233_power6.h"
#include "nbkernel234_power6.h"
#include "nbkernel300_power6.h"
#include "nbkernel301_power6.h"
#include "nbkernel302_power6.h"
#include "nbkernel303_power6.h"
#include "nbkernel304_power6.h"
#include "nbkernel310_power6.h"
#include "nbkernel311_power6.h"
#include "nbkernel312_power6.h"
#include "nbkernel313_power6.h"
#include "nbkernel314_power6.h"
#include "nbkernel320_power6.h"
#include "nbkernel321_power6.h"
#include "nbkernel322_power6.h"
#include "nbkernel323_power6.h"
#include "nbkernel324_power6.h"
#include "nbkernel330_power6.h"
#include "nbkernel331_power6.h"
#include "nbkernel332_power6.h"
#include "nbkernel333_power6.h"
#include "nbkernel334_power6.h"
#include "nbkernel400_power6.h"
#include "nbkernel410_power6.h"
#include "nbkernel420_power6.h"
#include "nbkernel430_power6.h"


static nb_kernel_t *
kernellist[eNR_NBKERNEL_NR] = 
{
    nbkernel010_power6,
    nbkernel020_power6,
    nbkernel030_power6,
    nbkernel100_power6,
    nbkernel101_power6,
    nbkernel102_power6,
    nbkernel103_power6,
    nbkernel104_power6,
    nbkernel110_power6,
    nbkernel111_power6,
    nbkernel112_power6,
    nbkernel113_power6,
    nbkernel114_power6,
    nbkernel120_power6,
    nbkernel121_power6,
    nbkernel122_power6,
    nbkernel123_power6,
    nbkernel124_power6,
    nbkernel130_power6,
    nbkernel131_power6,
    nbkernel132_power6,
    nbkernel133_power6,
    nbkernel134_power6,
    nbkernel200_power6,
    nbkernel201_power6,
    nbkernel202_power6,
    nbkernel203_power6,
    nbkernel204_power6,
    nbkernel210_power6,
    nbkernel211_power6,
    nbkernel212_power6,
    nbkernel213_power6,
    nbkernel214_power6,
    nbkernel220_power6,
    nbkernel221_power6,
    nbkernel222_power6,
    nbkernel223_power6,
    nbkernel224_power6,
    nbkernel230_power6,
    nbkernel231_power6,
    nbkernel232_power6,
    nbkernel233_power6,
    nbkernel234_power6,
    nbkernel300_power6,
    nbkernel301_power6,
    nbkernel302_power6,
    nbkernel303_power6,
    nbkernel304_power6,
    nbkernel310_power6,
    nbkernel311_power6,
    nbkernel312_power6,
    nbkernel313_power6,
    nbkernel314_power6,
    nbkernel320_power6,
    nbkernel321_power6,
    nbkernel322_power6,
    nbkernel323_power6,
    nbkernel324_power6,
    nbkernel330_power6,
    nbkernel331_power6,
    nbkernel332_power6,
    nbkernel333_power6,
    nbkernel334_power6,
    nbkernel400_power6,
    nbkernel410_power6,
    nbkernel430_power6,
	nbkernel010nf_power6,
    nbkernel020nf_power6,
    nbkernel030nf_power6,
    nbkernel100nf_power6,
    nbkernel101nf_power6,
    nbkernel102nf_power6,
    nbkernel103nf_power6,
    nbkernel104nf_power6,
    nbkernel110nf_power6,
    nbkernel111nf_power6,
    nbkernel112nf_power6,
    nbkernel113nf_power6,
    nbkernel114nf_power6,
    nbkernel120nf_power6,
    nbkernel121nf_power6,
    nbkernel122nf_power6,
    nbkernel123nf_power6,
    nbkernel124nf_power6,
    nbkernel130nf_power6,
    nbkernel131nf_power6,
    nbkernel132nf_power6,
    nbkernel133nf_power6,
    nbkernel134nf_power6,
    nbkernel200nf_power6,
    nbkernel201nf_power6,
    nbkernel202nf_power6,
    nbkernel203nf_power6,
    nbkernel204nf_power6,
    nbkernel210nf_power6,
    nbkernel211nf_power6,
    nbkernel212nf_power6,
    nbkernel213nf_power6,
    nbkernel214nf_power6,
    nbkernel220nf_power6,
    nbkernel221nf_power6,
    nbkernel222nf_power6,
    nbkernel223nf_power6,
    nbkernel224nf_power6,
    nbkernel230nf_power6,
    nbkernel231nf_power6,
    nbkernel232nf_power6,
    nbkernel233nf_power6,
    nbkernel234nf_power6,
    nbkernel300nf_power6,
    nbkernel301nf_power6,
    nbkernel302nf_power6,
    nbkernel303nf_power6,
    nbkernel304nf_power6,
    nbkernel310nf_power6,
    nbkernel311nf_power6,
    nbkernel312nf_power6,
    nbkernel313nf_power6,
    nbkernel314nf_power6,
    nbkernel320nf_power6,
    nbkernel321nf_power6,
    nbkernel322nf_power6,
    nbkernel323nf_power6,
    nbkernel324nf_power6,
    nbkernel330nf_power6,
    nbkernel331nf_power6,
    nbkernel332nf_power6,
    nbkernel333nf_power6,
    nbkernel334nf_power6,
    nbkernel400nf_power6,
    nbkernel410nf_power6,
    nbkernel430nf_power6,
};


void
nb_kernel_setup_power6(FILE *log,nb_kernel_t **list)
{
  int i;
  nb_kernel_t *p;

  if(log)
  {
#ifdef GMX_DOUBLE
    fprintf(log,"Configuring double precision IBM Power6-specific Fortran kernels...\n",);
#else
    fprintf(log,"Configuring single precision IBM Power6-specific Fortran kernels...\n",);
#endif
  }

  for(i=0;i<eNR_NBKERNEL_NR;i++)
  {
      p = kernellist[i];
      if(p!=NULL)
	list[i] = p; 
  }
}    
