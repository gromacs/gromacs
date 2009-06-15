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
 * Green Red Orange Magenta Azure Cyan Skyblue
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <math.h>

#include "copyrite.h"
#include "filenm.h"
#include "macros.h"
#include "pbc.h"
#include "physics.h"
#include "smalloc.h"
#include "statutil.h"
#include "vec.h"
#include "xvgr.h"

#include <trajana.h>

typedef struct {
  bool                 bAll;
  bool                 bDumpDist;
  bool                 bSplit1, bSplit2;
  int                  nat1, nat2;
  int                  fgi2;
  const char          *type;
  real                *tmp;
  real                *tmpd;
  FILE                *ofp;
  FILE                *dfp;
  gmx_ana_selection_t *refsel;
} t_angledata;

static void
calc_vec(int nat, rvec x[], t_pbc *pbc, rvec xout, rvec cout)
{
  switch (nat) {
    case 2:
      if (pbc)
        pbc_dx(pbc, x[1], x[0], xout);
      else
        rvec_sub(x[1], x[0], xout);
      svmul(0.5, xout, cout);
      rvec_add(x[0], cout, cout);
      break;
    case 3: {
      rvec v1, v2;
      if (pbc) {
        pbc_dx(pbc, x[1], x[0], v1);
        pbc_dx(pbc, x[2], x[0], v2);
      } else {
        rvec_sub(x[1], x[0], v1);
        rvec_sub(x[2], x[0], v2);
      }
      cprod(v1, v2, xout);
      rvec_add(x[0], x[1], cout);
      rvec_add(cout, x[2], cout);
      svmul(1.0/3.0, cout, cout);
      break;
      }
  }
}

static void write_file_frame(t_angledata *d, real time, int n, real ave)
{
  int i;

  if (n > 0) ave /= n;
  if (d->ofp) {
    fprintf(d->ofp, "%10.5f %8.3f %5d", time, ave, n);
    if (d->bAll) {
      for (i = 0; i < n; ++i)
          fprintf(d->ofp, " %8.3f", d->tmp[i]);
    }
    fprintf(d->ofp, "\n");
  }
  if (d->dfp) {
    for (i = 0; i < n; ++i) {
      if (d->bDumpDist)
        fprintf(d->dfp, "%8.3f ", d->tmpd[i]);
      fprintf(d->dfp, "%8.3f\n", d->tmp[i]);
    }
  }
}

static int
analyze_frame(t_topology *top, t_trxframe *fr, t_pbc *pbc,
              int nr, gmx_ana_selection_t *sel[], void *data)
{
  t_angledata *d = (t_angledata *)data;
  int          i, j, k, n;
  int          incr1, incr2;
  real         ave;
  rvec         v1, v2;
  rvec         c1, c2;
  rvec         x[3];
  real         angle;
  real         dist;

  if (d->type[0] == 'z') {
    v2[XX] = 0;
    v2[YY] = 0;
    v2[ZZ] = 1;
    clear_rvec(c2);
  }
  if (d->bSplit1) {
    n = sel[0]->p.nr;
    for (i = 1; i < d->nat1; ++i)
      if (sel[i]->p.nr != n)
        gmx_fatal(FARGS, "The first %d groups should contain the same number of positions", d->nat1);
  } else {
    if (sel[0]->p.nr % d->nat1 != 0)
        gmx_fatal(FARGS, "Number of positions in the first group not divisible by %d",
                  d->nat1);
  }
  if (d->nat2 > 0) {
    if (d->bSplit2) {
      n = sel[d->fgi2]->p.nr;
      for (i = 1; i < d->nat2; ++i)
        if (sel[d->fgi2+i]->p.nr != n)
          gmx_fatal(FARGS, "The last %d groups should contain the same number of positions", d->nat2);
    } else {
      if (sel[d->fgi2]->p.nr % d->nat2 != 0)
          gmx_fatal(FARGS, "Number of positions in the second group not divisible by %d",
                    d->nat2);
    }
  }
  incr1 = d->bSplit1 ? d->nat1 : 1;
  incr2 = d->bSplit2 ? d->nat2 : 1;
  ave = 0; n = 0;
  dist = 0;
  for (i = j = 0; i < sel[0]->p.nr; i += incr1) {
    if (d->bSplit1) {
      for (k = 0; k < d->nat1; ++k)
        copy_rvec(sel[k]->p.x[i], x[k]);
    } else {
      for (k = 0; k < d->nat1; ++k)
        copy_rvec(sel[0]->p.x[i+k], x[k]);
    }
    if (d->type[0] != 'a')
      calc_vec(d->nat1, x, pbc, v1, c1);
    switch (d->type[0]) {
      case 'a':
        if (pbc) {
          pbc_dx(pbc, x[1], x[0], v1);
          pbc_dx(pbc, x[1], x[2], v2);
        } else {
          rvec_sub(x[1], x[0], v1);
          rvec_sub(x[1], x[2], v2);
        }
        break;
      case 'v':
      case 'p':
        if (d->bSplit2) {
          for (k = 0; k < d->nat2; ++k)
            copy_rvec(sel[d->fgi2+k]->p.x[i], x[k]);
        } else {
          for (k = 0; k < d->nat2; ++k)
            copy_rvec(sel[d->fgi2]->p.x[i+k], x[k]);
        }
        calc_vec(d->nat2, x, pbc, v2, c2);
        j += incr2;
        break;
      case 'z':
        c1[XX] = c1[YY] = 0;
        break;
      case 's':
        copy_rvec(d->refsel->p.x[0], c2);
        if (pbc)
          pbc_dx(pbc, c1, c2, v2);
        else
          rvec_sub(c1, c2, v2);
        break;
    }
    angle = acos(cos_angle(v1, v2)) * RAD2DEG;
    if (d->bDumpDist) {
      if (pbc) {
        pbc_dx(pbc, c2, c1, v1);
        dist = norm(v1);
      } else
        dist = sqrt(distance2(c1, c2));
    }
    d->tmp[n] = angle;
    d->tmpd[n] = dist;
    ave += angle;
    ++n;
  }
  write_file_frame(d, fr->time, n, ave);
  return 0;
}

int gmx_sgangle(int argc,char *argv[])
{
  const char *desc[] = {
    "g_sgangle computes angles between generic vectors.",
    "It supports both vectors defined by two positions and normals of planes",
    "defined by three positions.",
    "The z axis or the local normal of a sphere can also be used as ",
    "one of the vectors. ",
    "There is also a convenience option 'angle' for calculation of an angle ",
    "defined by three positions.[PAR]",
    "The type of the angle is specified with [TT]-g1[tt] and [TT]-g2[tt]. ",
    "If [TT]-g1[tt] is [TT]angle[tt], [TT]-g2[tt] should not be specified. ",
    "In this case, one selection is required, and it should contain",
    "triplets of positions that define the angles to be calculated.",
    "If [TT]-g1[tt] is not [TT]angle[tt], [TT]-g2[tt] should not be ",
    "[TT]none[tt], and the two options define the two vectors for ",
    "the calculation. For vectors ([TT]vec[tt]), a selection with",
    "pairs of positions is required, and for planes ([TT]plane[tt]),",
    "triplets of positions are required.",
    "If both vectors are specified by positions, the number of vectors ",
    "should be the same in both selections. ",
    "[TT]-g1 sphnorm[tt] requires a reference selection that defines",
    "the center of the sphere.",
    "[TT]-g1 z[tt] does not require any selection.[PAR]",
    "With [TT]-split1[tt], the positions for [TT]-g1[tt] are specified",
    "using N separate selections with M positions each, instead of the",
    "default M*N positions in one selection.",
    "[TT]-split2[tt] does the same for [TT]-g2[tt].[PAR]",
    "There are two options for output: ",
    "[TT]-o[tt] writes an xvgr file with the time, the average angle ",
    "and the number of angles calculated for each frame. ",
    "With [TT]-all[tt], also the individual angles are written. ",
    "[TT]-od[tt] can be used to dump all the individual angles, ",
    "each on a separate line. This format is better suited for ",
    "further processing, e.g., if angles from multiple runs are needed. ",
  };

  static const char *g1type[] = {NULL, "angle", "vec", "plane", "z", "sphnorm", NULL };
  static const char *g2type[] = {NULL, "none", "vec", "plane", NULL };
  static bool bAll = FALSE;
  static bool bSplit1 = FALSE, bSplit2 = FALSE;
  static bool bDumpDist = FALSE;
  static t_pargs pa[] = {
    { "-g1",  FALSE, etENUM, {g1type},
      "Type of the first analysis/reference group" },
    { "-g2",  FALSE, etENUM, {g2type},
      "Type of the second analysis group" },
    { "-split1", FALSE, etBOOL, {&bSplit1},
      "Each selection of the first group in a separate selection" },
    { "-split2", FALSE, etBOOL, {&bSplit2},
      "Each selection of the second group in a separate selection" },
    { "-all", FALSE, etBOOL, {&bAll},
      "Print individual angles together with the average" },
    { "-dumpd", FALSE, etBOOL, {&bDumpDist},
      "Dump distances together with the angles with -od" },
  };

  static t_filenm fnm[] = {
    { efXVG,  "-o", "angle", ffOPTWR },
    { efXVG, "-od", "adump", ffOPTWR },
  };

  gmx_ana_traj_t       *trj;
  gmx_ana_selection_t **sel;
  t_angledata           d;
  int                   ngrps;
  int                   g;
  int                   na1, na2;

#define NFILE asize(fnm)

  CopyRight(stderr,argv[0]);
  gmx_ana_traj_create(&trj, ANA_REQUIRE_WHOLE | ANA_USER_SELINIT);
  parse_trjana_args(trj, &argc, argv, PCA_CAN_VIEW,
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL);

  /* Make some validity checks */
  if (g1type[0][0] == 'a' && g2type[0][0] != 'n')
    gmx_fatal(FARGS,"Cannot use a second group (-g2) with -g1 angle");
  if (g1type[0][0] != 'a' && g2type[0][0] == 'n')
    gmx_fatal(FARGS,"Should specify a second group (-g2) if the first group "
              "is not an angle");
  if (g1type[0][0] == 'a' && bDumpDist) {
    fprintf(stderr,"\nCannot calculate distances with -g1 angle");
    bDumpDist = FALSE;
  }
  d.bSplit1 = bSplit1;
  d.bSplit2 = bSplit2;

  /* Read the index file. */
  if (g1type[0][0] == 'a') {
    fprintf(stderr,"\nSelect a group for angle calculation "
            "(should contain triplets of atoms):\n");
    ngrps = 1;
    d.nat1 = 3;
    d.nat2 = 0;
  } else if (g1type[0][0] == 'z') {
    fprintf(stderr,"\nWill use the z axis as reference.\n");
    fprintf(stderr,"\nSelect a group for %s calculation "
            "(should contain %ss of atoms):\n",
            g2type[0][0] == 'v' ? "vector" : "plane normal",
            g2type[0][0] == 'v' ? "pair" : "triplet");
    ngrps = 1;
    d.bSplit1 = d.bSplit2;
    d.nat1 = (g2type[0][0] == 'v') ? 2 : 3;
    d.nat2 = 0;
  } else if (g1type[0][0] == 's') {
    fprintf(stderr,"\nSelect a group for the center of the sphere and "
            "a group for %s calculation:\n",
            g2type[0][0] == 'v' ? "vector" : "plane normal");
    gmx_ana_set_nrefgrps(trj, 1);
    ngrps = 1;
    d.bSplit1 = d.bSplit2;
    d.nat1 = (g2type[0][0] == 'v') ? 2 : 3;
    d.nat2 = 0;
  } else {
    fprintf(stderr,"\nSelect two groups for calculation:\n");
    gmx_ana_set_nanagrps(trj, 2);
    ngrps = 2;
    d.nat1 = (g1type[0][0] == 'v') ? 2 : 3;
    d.nat2 = (g2type[0][0] == 'v') ? 2 : 3;
  }
  if (d.bSplit1)
    ngrps += d.nat1 - 1;
  if (d.bSplit2 && d.nat2 > 0)
    ngrps += d.nat2 - 1;
  d.fgi2 = d.bSplit1 ? d.nat1 : 1;
  gmx_ana_set_nanagrps(trj, ngrps);
  gmx_ana_init_selections(trj);
  if (g1type[0][0] == 's') {
    gmx_ana_get_refsel(trj, 0, &d.refsel);
  } else {
    d.refsel = NULL;
  }
  gmx_ana_get_anagrps(trj, &sel);

  /* Check the sizes of the indexes */
  if (!d.bSplit1 && d.nat1 > 0 && sel[0]->p.nr % d.nat1 != 0)
    gmx_fatal(FARGS, "Number of positions in the first group not divisible by %d",
              d.nat1);
  if (!d.bSplit2 && d.nat2 > 0 && sel[1]->p.nr % d.nat2 != 0)
    gmx_fatal(FARGS, "Number of positions in the second group not divisible by %d",
              d.nat2);
  if (d.bSplit1) {
    na1 = sel[0]->p.nr;
    for (g = 1; g < d.nat1; ++g)
      if (sel[g]->p.nr != na1)
        gmx_fatal(FARGS, "The first %d groups should contain the same number of positions", d.nat1);
  } else
    na1 = (d.nat1 > 0 ? sel[0]->p.nr / d.nat1 : 0);
  if (d.bSplit2 && d.nat2 > 0) {
    na2 = sel[d.fgi2]->p.nr;
    for (g = d.fgi2+1; g < ngrps; ++g)
      if (sel[g]->p.nr != na2)
        gmx_fatal(FARGS, "The last %d groups should contain the same number of positions", d.nat2);
  } else
    na2 = (d.nat2 > 0 ? sel[d.fgi2]->p.nr / d.nat2 : 0);
  if (d.nat1 > 0 && d.nat2 > 0 && na1 != na2)
    gmx_fatal(FARGS, "Number of vectors defined by the two groups are not the same");

  d.bAll = bAll;
  d.bDumpDist = bDumpDist;
  d.type = g1type[0];
  snew(d.tmp, max(na1, na2));
  snew(d.tmpd, max(na1, na2));
  if (opt2bSet("-o",NFILE,fnm) || !opt2bSet("-od",NFILE,fnm)) {
    d.ofp = xvgropen(opt2fn("-o",NFILE,fnm),"","Time (ps)","Angle (degrees)");
    xvgr_selections(d.ofp, trj);
  } else
    d.ofp = NULL;
  if (opt2bSet("-od",NFILE,fnm))
    d.dfp = ffopen(opt2fn("-od",NFILE,fnm),"w");
  else
    d.dfp = NULL;

  gmx_ana_do(trj, 0, &analyze_frame, &d);

  if (d.ofp)
    fclose(d.ofp);
  if (d.dfp)
    fclose(d.dfp);

  thanx(stderr);

  return 0;
}
