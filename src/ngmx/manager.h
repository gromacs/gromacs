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
 * Good gRace! Old Maple Actually Chews Slate
 */

#ifndef _manager_h
#define _manager_h

static char *SRCID_manager_h = "$Id$";

#include <stdio.h>
#include "typedefs.h"
#include "x11.h"
#include "xutil.h"
#include "3dview.h"
#include "nleg.h"
#include "buttons.h"

/* Some window sizes */
#define EWIDTH  	200
#define EHEIGHT 	  0
#define LDHEIGHT  	  0
#define LEGHEIGHT 	 60

typedef enum { eOSingle, eOBond, eOHBond, eONR } eObject;

typedef enum { eVNormal, eVSpecial, eVHidden, evNR } eVisible;

enum { eBThin, eBFat, eBVeryFat, eBSpheres, eBNR };

typedef struct {
  t_windata wd;			/* Mol window structure			*/
  bool      bShowHydrogen;	/* Show Hydrogens?			*/
  int       bond_type;		/* Show one of the above bondtypes      */
  bool      bBoxSelect;		/* Show Box?				*/
} t_molwin;

typedef struct {
  eObject 	eO;		/* The type of object			*/
  eVisible	eV;		/* Visibility status of the object	*/
  ulong   	color;		/* The color (only when eV==evSpecial) 	*/
  atom_id 	ai,aj;		/* The atom_id for i (and j if bond)	*/
  real    	z;		/* The Z-coordinate for depht cueing	*/
} t_object;

typedef struct {
  t_block *grps;		/* Blocks with atom numbers		*/
  char    **grpnames;		/* The names of the groups		*/
  bool    *bShow;		/* Show a group ?			*/
} t_filter;

/*
 * t_manager structure:
 *
 * This structure manages the display area for the gmx program.
 * It reads the status file and sends messages when windows need to
 * be updated.
 *
 */
typedef struct {
  int       status;
  char      *trajfile;
  int       natom;		/* The number of atoms			*/
  t_topology top;               /* topology                             */
  rvec      box_size;
  int       step;               /* The actual step number               */
  real      time;               /* The actual time                      */
  rvec      *x;			/* The coordinates			*/
  iv2       *ix;		/* The coordinates after projection	*/
  real      *zz;                /* Z-coords                             */
  matrix    box;		/* The box				*/
  int       nobj;		/* The number of objects		*/
  t_object  *obj;		/* The objects on screen		*/
  bool      *bHydro;		/* TRUE for hydrogen atoms		*/
  bool      *bLabel;            /* Show a label on atom i?              */
  char      **szLab;            /* Array of pointers to labels          */
  ulong     *col;		/* The colour of the atoms		*/
  int       *size;		/* The size of the atoms		*/
  real      *vdw;		/* The VDWaals radius of the atoms	*/
  bool      *bVis;              /* visibility of atoms                  */
  bool      bPbc;               /* Remove Periodic boundary             */
  bool      bAnimate;		/* Animation going on?			*/
  bool      bEof;               /* End of file reached?                 */
  bool      bStop;              /* Stopped by user?                     */
  bool      bSort;		/* Sort the coordinates			*/
  bool      bPlus;		/* Draw plus for single atom		*/
  int       nSkip;		/* Skip n steps after each frame	*/

  t_windata   wd;               /* The manager subwindow                */
  t_windata   title;		/* Title window				*/
  t_3dview    *view;            /* The 3d struct                        */
  t_molwin    *molw;		/* The molecule window			*/
  t_butbox    *vbox;		/* The video box			*/
  t_butbox    *bbox;		/* The button box			*/
  t_legendwin *legw;		/* The legend window			*/
} t_manager;

extern t_manager *init_man(t_x11 *x11,Window Parent,
			   int x,int y,int width,int height,
			   ulong fg,ulong bg);
/* Initiate the display manager */

extern void move_man(t_x11 *x11,t_manager *man,int width,int height);
/* Set the right size for this window */

extern void step_message(t_x11 *x11,t_manager *man);
/* Send a message to the manager */

extern void set_file(t_x11 *x11,t_manager *man,char *trajectory,char *status,
		     char *fnvdw);
/* Read a new trajectory and topology */

extern void map_man(t_x11 *x11,t_manager *man);

extern void move_man(t_x11 *x11,t_manager *man,int width,int height);

extern bool toggle_animate (t_x11 *x11,t_manager *man);

extern bool toggle_pbc (t_manager *man);

extern void no_labels(t_x11 *x11,t_manager *man);
/* Turn off all labels */

extern void done_man(t_x11 *x11,t_manager *man);
/* Clean up man struct */

extern void draw_mol(t_x11 *x11,t_manager *man);

extern void create_visibility(t_manager *man);

extern void do_filter(t_x11 *x11,t_manager *man,t_filter *filter);

#endif
