/*
 *       $Id$
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0
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
 * GROup of MAchos and Cynical Suckers
 */

#ifndef _buttons_h
#define _buttons_h

static char *SRCID_buttons_h = "$Id$";

#include <xutil.h>

enum { 
  IDROTX,IDROTY,IDROTZ,IDTRANSX,IDTRANSY,IDTRANSZ,IDZOOM,IDBUTNR,
  IDREWIND,IDSTEP,IDFF,IDSTOP_ANI,IDNR,
  IDDRAWMOL,IDLABEL
  };

#define AIR      3            /* extra space between child windows */
#define BORDER  1

#define EnterLeave (EnterWindowMask | LeaveWindowMask)

typedef struct {
  t_windata wd;
  int       ID;
} t_mwbut;

typedef struct {
  t_windata wd;
  int       nbut;
  t_mwbut   *b;
} t_butbox;

extern t_butbox *init_vbox(t_x11 *x11,Window Parent,Window SendTo,
			   ulong fg,ulong bg);

extern void set_vbtime(t_x11 *x11,t_butbox *vbox,char *text);

extern t_butbox *init_bbox(t_x11 *x11,Window Parent,Window SendTo,
			   int width,ulong fg,ulong bg);

extern void show_but(t_x11 *x11,t_butbox *bbox);

extern void hide_but(t_x11 *x11,t_butbox *bbox);

extern void done_bbox(t_x11 *x11,t_butbox *bbox);

#endif
