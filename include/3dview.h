/*
 * $Id$
 * 
 *       This source code is part of
 * 
 *        G   R   O   M   A   C   S
 * 
 * GROningen MAchine for Chemical Simulations
 * 
 *               VERSION 2.0
 * 
 * Copyright (c) 1991-1999
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 * 
 * Also check out our WWW page:
 * http://md.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 * 
 * And Hey:
 * Good ROcking Metal Altar for Chronical Sinners
 */

#ifndef _3dview_h
#define _3dview_h

static char *SRCID_3dview_h = "$Id$";

#define WW 3

typedef real vec4[4];

typedef real mat4[4][4];

typedef int  iv2[2];

typedef struct {
  matrix box;
  vec4   eye,origin;	/* The eye and origin position		*/
  mat4   proj;		/* Projection matrix 			*/
  mat4   Rot;           /* Total rotation matrix                */
  real   sc_x,sc_y;	/* Scaling for aspect ratio		*/
} t_3dview;

extern void print_m4(FILE *fp,char *s,mat4 A);

extern void print_v4(FILE *fp,char *s,int dim,real *a);

extern void m4_op(mat4 m,rvec x,vec4 v);

extern void unity_m4(mat4 m);

extern void mult_matrix(mat4 A, mat4 B, mat4 C);

extern void rotate(int axis, real angle, mat4 A);

extern void translate(real tx, real ty, real tz, mat4 A);

extern void m4_op(mat4 m,rvec x,vec4 v);

extern void calculate_view(t_3dview *view);

extern t_3dview *init_view(matrix box);
/* Generate the view matrix from the eye pos and the origin,
 * applying also the scaling for the aspect ration.
 * There is no accompanying done_view routine: the struct can simply
 * be sfree'd.
 */

/* The following options are present on the 3d struct:
 * zoom (scaling)
 * rotate around the center of the box
 * reset the view
 */

extern bool zoom_3d(t_3dview *view,real fac);
/* Zoom in or out with factor fac, returns TRUE when zoom succesful,
 * FALSE otherwise.
 */

extern void rotate_3d(t_3dview *view,int axis,bool bPositive);
/* Rotate the eye around the center of the box, around axis */

extern void translate_view(t_3dview *view,int axis,bool bPositive);
/* Translate the origin at which one is looking */

extern void reset_view(t_3dview *view);
/* Reset the viewing to the initial view */

#endif

