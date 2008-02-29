/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.3.3
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

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
 * Groningen Machine for Chemical Simulation
 */

#ifndef _writeps_h
#define _writeps_h

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include "typedefs.h"

typedef enum { 
  eXCenter, eXLeft, eXRight
} eXPos;

typedef enum { 
  eYCenter, eYTop,  eYBottom
} eYPos; 

enum { 
  efontTIMES, efontTIMESITALIC, efontTIMESBOLD, efontTIMESBOLDITALIC,
  efontHELV,  efontHELVITALIC,  efontHELVBOLD,  efontHELVBOLDITALIC,
  efontCOUR,  efontCOURITALIC,  efontCOURBOLD,  efontCOURBOLDITALIC,
  efontNR };


typedef struct t_int_psdata *t_psdata;
/* Only use t_psdata - it is a pointer to an abstract datatype
 * that maintains the state of the postscript currently written.
 */

extern const char *fontnm[efontNR];

extern t_psdata ps_open(char *fn,real x1,real y1,real x2,real y2);

extern void ps_linewidth(t_psdata ps,int lw);
extern void ps_color(t_psdata ps,real r,real g,real b);
extern void ps_rgb(t_psdata ps,t_rgb *rgb);

extern void ps_rgb_box(t_psdata ps,t_rgb *rgb);
extern void ps_rgb_nbox(t_psdata ps,t_rgb *rgb,real n);
extern void ps_init_rgb_box(t_psdata ps,real xbox, real ybox);
extern void ps_init_rgb_nbox(t_psdata ps,real xbox, real ybox);

extern void ps_lineto(t_psdata ps,real x,real y);
extern void ps_linerel(t_psdata ps,real dx,real dy);

extern void ps_moveto(t_psdata ps,real x,real y);
extern void ps_moverel(t_psdata ps,real dx,real dy);

extern void ps_line(t_psdata ps,real x1,real y1,real x2,real y2);

extern void ps_box(t_psdata ps,real x1,real y1,real x2,real y2);
extern void ps_fillbox(t_psdata ps,real x1,real y1,real x2,real y2);

extern void ps_arc(t_psdata ps,real x1,real y1,real rad,real a0,real a1);
extern void ps_fillarc(t_psdata ps,real x1,real y1,real rad,real a0,real a1);
extern void ps_arcslice(t_psdata ps,real xc,real yc,
			real rad1,real rad2,real a0,real a1);
extern void ps_fillarcslice(t_psdata ps,real xc,real yc,
			    real rad1,real rad2,real a0,real a1);

extern void ps_circle(t_psdata ps,real x1,real y1,real rad);

extern void ps_font(t_psdata ps,int font,real size);
extern void ps_strfont(t_psdata ps,char *font,real size);

extern void ps_text(t_psdata ps,real x1,real y1,char *str);
extern void ps_ctext(t_psdata ps,real x1,real y1,char *str,int expos);

extern void ps_close(t_psdata ps);

extern void ps_flip(t_psdata ps,bool bPlus);
/* Rotate over 90 (bPlus) or -90 (!bPlus) degrees */

extern void ps_rotate(t_psdata ps,real angle);

extern void ps_translate(t_psdata ps,real x,real y);

extern void ps_setorigin(t_psdata ps);
extern void ps_unsetorigin(t_psdata ps);

extern void viewps(char *fn);

extern void ps_comment(t_psdata ps,char *s);

#endif	/* _writeps_h */
