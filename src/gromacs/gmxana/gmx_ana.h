/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2017, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */

#ifndef _gmx_ana_h
#define _gmx_ana_h

#ifdef __cplusplus
extern "C" {
#endif

int
gmx_anadock(int argc, char *argv[]);

int
gmx_analyze(int argc, char *argv[]);

int
gmx_anaeig(int argc, char *argv[]);

int
gmx_g_angle(int argc, char *argv[]);

int
gmx_bar(int argc, char *argv[]);

int
gmx_bundle(int argc, char *argv[]);

int
gmx_chi(int argc, char *argv[]);

int
gmx_cluster(int argc, char *argv[]);

int
gmx_confrms(int argc, char *argv[]);

int
gmx_covar(int argc, char *argv[]);

int
gmx_current(int argc, char *argv[]);

int
gmx_density(int argc, char *argv[]);

int
gmx_densmap(int argc, char *argv[]);

int
gmx_densorder(int argc, char *argv[]);

int
gmx_dielectric(int argc, char *argv[]);

int
gmx_dipoles(int argc, char *argv[]);

int
gmx_disre(int argc, char *argv[]);

int
gmx_do_dssp(int argc, char *argv[]);

int
gmx_dos(int argc, char *argv[]);

int
gmx_dyecoupl(int argc, char *argv[]);

int
gmx_dyndom(int argc, char *argv[]);

int
gmx_editconf(int argc, char *argv[]);

int
gmx_eneconv(int argc, char *argv[]);

int
gmx_enemat(int argc, char *argv[]);

int
gmx_energy(int argc, char *argv[]);

int
gmx_lie(int argc, char *argv[]);

int
gmx_filter(int argc, char *argv[]);

int
gmx_genion(int argc, char *argv[]);

int
gmx_genpr(int argc, char *argv[]);

int
gmx_gyrate(int argc, char *argv[]);

int
gmx_h2order(int argc, char *argv[]);

int
gmx_hbond(int argc, char *argv[]);

int
gmx_helix(int argc, char *argv[]);

int
gmx_helixorient(int argc, char *argv[]);

int
gmx_hydorder(int argc, char *argv[]);

int
gmx_make_edi(int argc, char *argv[]);

int
gmx_make_ndx(int argc, char *argv[]);

int
gmx_mindist(int argc, char *argv[]);

int
gmx_mk_angndx(int argc, char *argv[]);

int
gmx_msd(int argc, char *argv[]);

int
gmx_morph(int argc, char *argv[]);

int
gmx_nmeig(int argc, char *argv[]);

int
gmx_nmens(int argc, char *argv[]);

int
gmx_nmr(int argc, char *argv[]);

int
gmx_nmtraj(int argc, char *argv[]);

int
gmx_order(int argc, char *argv[]);

int
gmx_polystat(int argc, char *argv[]);

int
gmx_potential(int argc, char *argv[]);

int
gmx_principal(int argc, char *argv[]);

int
gmx_rama(int argc, char *argv[]);

int
gmx_rotmat(int argc, char *argv[]);

int
gmx_rms(int argc, char *argv[]);

int
gmx_rmsdist(int argc, char *argv[]);

int
gmx_rmsf(int argc, char *argv[]);

int
gmx_rotacf(int argc, char *argv[]);

int
gmx_saltbr(int argc, char *argv[]);

int
gmx_sham(int argc, char *argv[]);

int
gmx_sigeps(int argc, char *argv[]);

int
gmx_sorient(int argc, char *argv[]);

int
gmx_spol(int argc, char *argv[]);

int
gmx_spatial(int argc, char *argv[]);

int
gmx_tcaf(int argc, char *argv[]);

int
gmx_traj(int argc, char *argv[]);

int
gmx_trjcat(int argc, char *argv[]);

int
gmx_trjconv(int argc, char *argv[]);

int
gmx_trjorder(int argc, char *argv[]);

int
gmx_tune_pme(int argc, char *argv[]);

int
gmx_velacc(int argc, char *argv[]);

int
gmx_clustsize(int argc, char *argv[]);

int
gmx_mdmat(int argc, char *argv[]);

int
gmx_vanhove(int argc, char *argv[]);

int
gmx_wham(int argc, char *argv[]);

int
gmx_wheel(int argc, char *argv[]);

int
gmx_xpm2ps(int argc, char *argv[]);

int
gmx_membed(int argc, char *argv[]);

int
gmx_pme_error(int argc, char *argv[]);

int
gmx_sans(int argc, char *argv[]);

int
gmx_saxs(int argc, char *argv[]);

#ifdef BUILD_WITH_FDA
int
gmx_fda_graph(int argc, char *argv[]);

int
gmx_fda_get_stress(int argc, char *argv[]);

int
gmx_fda_view_stress(int argc, char *argv[]);

int
gmx_fda_shortest_path(int argc, char *argv[]);
#endif

#ifdef __cplusplus
}
#endif

#endif
/* _gmx_ana_h */
