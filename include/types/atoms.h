/*
 *       @(#) copyrgt.c 1.12 9/30/97
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0b
 * 
 * Copyright (c) 1990-1997,
 * BIOSON Research Institute, Dept. of Biophysical Chemistry,
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
 * Good ROcking Metal Altar for Chronical Sinners
 */
enum {
  eptAtom, eptNucleus, eptShell, eptBond, eptDummy, eptNR
};
/* The particle type */
 
enum {
  egcTC,    egcENER,   egcACC,   egcFREEZE, 
  egcUser1, egcUser2,  egcUser3, egcXTC,
  egcNR 
};

typedef struct {
  real 		m,q;		/* Mass and charge			*/
  real 		mB,qB;		/* Mass and charge for Free Energy calc */
  ushort	type;		/* Atom type				*/
  ushort	typeB;		/* Atom type for Free Energy calc	*/
  int           ptype;		/* Particle type			*/
  int 		resnr;		/* Residue number			*/
  unsigned char grpnr[egcNR];   /* Group numbers			*/
  unsigned char chain;          /* chain identifier                     */
} t_atom;

typedef struct {
  int  type;                    /* PDB record name                      */
  int  atomnr;                  /* PDB atom number                      */
  char pdbresnr[12];            /* PDB res number                       */
  real occup;                   /* Occupancy                            */
  real bfac;                    /* B-factor                             */
  bool bAnisotropic;            /* (an)isotropic switch                 */
  int  uij[6];                  /* Anisotropic B-factor                 */
} t_pdbinfo;

typedef struct {
  int  nr;			/* Number of different groups		*/
  int  *nm_ind;                 /* Index in the group names             */
} t_grps;

typedef struct {
  int           nr;             /* Nr of atoms                          */
  t_atom	*atom;		/* Array of atoms (dim: nr)		*/
				/* The following entries will not 	*/
				/* allways be used (nres==0)	 	*/
  char		***atomname;	/* Array of pointers to atom name	*/
				/* use: (*(atomname[i]))		*/
  int		nres;		/* Nr of residue names			*/
  char		***resname; 	/* Array of pointers to residue names 	*/
				/* use: (*(resname[i]))			*/
  int           ngrpname;       /* Number of groupnames                 */
  char          ***grpname;	/* Names of the groups		        */
  t_block	excl;		/* Exclusions				*/
  t_grps        grps[egcNR];    /* Groups of things                     */
  t_pdbinfo     *pdbinfo;       /* PDB Information, such as aniso. Bfac */
} t_atoms;

#define PERTURBED(a) (((a).mB != (a).m) || ((a).qB != (a).q) || ((a).typeB != (a).type))
