/*
 *       $Id$
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 1.6
 * 
 * Copyright (c) 1991-1997
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
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
 * Green Red Orange Magenta Azure Cyan Skyblue
 */

#ifndef _xvgr_h
#define _xvgr_h

static char *SRCID_xvgr_h = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) xvgr.h 1.8 7/28/97"
#endif /* HAVE_IDENT */

#ifdef CPLUSPLUS
extern "C" {
#endif

#include "sysstuff.h"
#include "typedefs.h"
/***************************************************
 *            XGraph routines
 ***************************************************/
extern FILE *xgopen(char *fn,char *header,char *title,char *xaxis,char *yaxis);
/* Open a file, and write a title, and axis-labels in XGraph format */

extern void xgraph_file(char *fn,char *opts);
/* Starts xgraph with a file fn in the background,
 * opts (options to xgraph) may be NULL
 */

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

extern FILE *xvgropen(char *fn,char *title,char *xaxis,char *yaxis);
/* Open a file, and write a title, and axis-labels in Xvgr format */

extern void xvgr_file(char *fn,char *opts);
/* Starts xvgr with a file fn in the background,
 * opts (options to xvgr) may be NULL
 */

extern void xvgr_view(FILE *out,real xmin,real ymin,real xmax,real ymax);
/* Set the view in xvgr */

extern void xvgr_world(FILE *out,real xmin,real ymin,real xmax,real ymax);
/* Set the world in xvgr */

extern void xvgr_legend(FILE *out,int nsets,char *setname[]);
/* Make a legend box, and also modifies the view to make room for the legend */

extern void xvgr_line_props(FILE *out,int NrSet,int LineStyle,int LineColor);
/* Set xvgr line styles and colors */

extern void xvgr_box(FILE *out,
		     int LocType,
		     real xmin,real ymin,real xmax,real ymax,
		     int LineStyle,int LineWidth,int LineColor,
		     int BoxFill,int BoxColor,int BoxPattern);
/* Make a box */

extern int read_xvg(char *fn,real ***y,int *ny);
/* Read an xvg file for post processing. The number of rows is returned
 * fn is the filename, y is a pointer to a 2D array (to be allocated by
 * the routine) ny is the number of columns (including X if appropriate)
 */
 
extern void dump_xvg(char *fn,char *title,int nx,int ny,real **y);
/* Quicly dump a two D array (y) of dimensions nx rows times
 * ny columns to a file.
 */


/****************************************************
 *           Some statistics utilities 
 ****************************************************/
extern void lsq_y_ax(int n, real x[], real y[], real *a);
/* Fit a straight line y=ax thru the n data points x,y */

extern void lsq_y_ax_b(int n, real x[], real y[], real *a, real *b);
/* Fit a straight line y=ax+b thru the n data points x,y */

#ifdef CPLUSPLUS
}
#endif

#endif	/* _xvgr_h */
