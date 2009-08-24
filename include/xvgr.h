/*
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
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
 * Gromacs Runs On Most of All Computer Systems
 */

#ifndef _xvgr_h
#define _xvgr_h

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

#include "sysstuff.h"
#include "typedefs.h"
#include "viewit.h"

/***************************************************
 *            XVGR   DEFINITIONS
 ***************************************************/
enum {
  elNone, elSolid, elDotted, elDashed, 
  elLongDashed, elDotDashed, elNR
};
/* xvgr line-styles */

enum {  
  ecWhite, ecFrank, ecBlack=ecFrank,
  ecRed, ecGreen, ecBlue, ecYellow, ecBrown, ecGray, 
  ecPurple, ecLightBlue, ecViolet, ecHolland, ecLila, ecDarkGray, 
  ecAquamarine, ecOlive, ecNR
};
/* xvgr line-colors */

enum {
  eppNone, eppColor, eppPattern, eppNR
};
/* xvgr pattern type */

enum {
  evView, evWorld, evNR
};
/* view type */

/***************************************************
 *            XVGR   ROUTINES
 ***************************************************/

extern bool use_xmgr(void);
/* Returns if we use xmgr instead of xmgrace */

extern FILE *xvgropen(const char *fn,const char *title,const char *xaxis,
                      const char *yaxis,const output_env_t oenv);
/* Open a file, and write a title, and axis-labels in Xvgr format */

/* Close xvgr file, and clean up internal file buffers correctly */
extern void xvgrclose(FILE *fp);

extern void xvgr_subtitle(FILE *out,const char *subtitle,
                          const output_env_t oenv);
/* Set the subtitle in xvgr */

extern void xvgr_view(FILE *out,real xmin,real ymin,real xmax,real ymax,        
                      const output_env_t oenv);
/* Set the view in xvgr */

extern void xvgr_world(FILE *out,real xmin,real ymin,real xmax,real ymax,
                       const output_env_t oenv);
/* Set the world in xvgr */

extern void xvgr_legend(FILE *out,int nsets,char **setnames,
                        const output_env_t oenv);
/* Make a legend box, and also modifies the view to make room for the legend */

extern void xvgr_line_props(FILE *out,int NrSet,int LineStyle,int LineColor,
                            const output_env_t oenv);
/* Set xvgr line styles and colors */

extern void xvgr_box(FILE *out,
		     int LocType,
		     real xmin,real ymin,real xmax,real ymax,
		     int LineStyle,int LineWidth,int LineColor,
		     int BoxFill,int BoxColor,int BoxPattern,
                     const output_env_t oenv);
/* Make a box */

extern int read_xvg_legend(const char *fn,double ***y,int *ny,char ***legend);
/* Read an xvg file for post processing. The number of rows is returned
 * fn is the filename, y is a pointer to a 2D array (to be allocated by
 * the routine) ny is the number of columns (including X if appropriate).
 * If legend!=NULL, read the legends for the sets (when present),
 * 0 is the first y legend, the legend string will be NULL when not present.
 */

extern int read_xvg(const char *fn,double ***y,int *ny);
/* As read_xvg_legend, but does not read legends. */
 
extern void write_xvg(const char *fn,const char *title,int nx,int ny,real **y,
                      char **leg, const output_env_t oenv);
/* Write a two D array (y) of dimensions nx rows times
 * ny columns to a file. If leg != NULL it will be written too.
 */


/* This function reads ascii (xvg) files and extracts the data sets to a 
 * two dimensional array which is returned.
 */
extern real **read_xvg_time(const char *fn,
			    bool bHaveT,
			    bool bTB,real tb,
			    bool bTE,real te,
			    int nsets_in,int *nset,int *nval,
			    real *dt,real **t);
#ifdef __cplusplus
}
#endif

#endif	/* _xvgr_h */
