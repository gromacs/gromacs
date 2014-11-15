/*
 * This source file is part of the Alexandria project.
 *
 * Copyright (C) 2014 David van der Spoel and Paul J. van Maaren
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */
/*! \internal \brief
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef POLDATA_H
#define POLDATA_H

#include "gromacs/legacyheaders/types/commrec.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/legacyheaders/network.h"

/* This source code file is part of the Alexandria project */

/*! \brief
 * Enumerated type holding the charge generation models used in PolData
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */
enum ChargeGenerationModel {
    eqgNone,
    eqgAXp,
    eqgAXg,
    eqgAXs,
    eqgESP,
    eqgRESP,
    eqgRESPG,
    eqgYang,
    eqgBultinck,
    eqgRappe,
    eqgNR
};

enum DihedralType {
    egdPDIHS,
    egdIDIHS,
    egdNR
};

typedef struct gmx_poldata *gmx_poldata_t;

extern gmx_poldata_t gmx_poldata_init();

extern void gmx_poldata_set_filename(gmx_poldata_t pd, char *fn2);

extern int gmx_poldata_get_natypes(gmx_poldata_t pd);
extern int gmx_poldata_get_nptypes(gmx_poldata_t pd);
extern int gmx_poldata_get_ngt_bond(gmx_poldata_t pd);
extern int gmx_poldata_get_ngt_angle(gmx_poldata_t pd);
extern int gmx_poldata_get_ngt_dihedral(gmx_poldata_t pd, int egd);

extern void gmx_poldata_add_bonding_rule(gmx_poldata_t pd,
                                         char *gt_brule, char *atype,
                                         char *geometry, int numbonds,
                                         double valence, int iAromatic,
                                         char *neighbors);

extern int gmx_poldata_get_bonding_rule(gmx_poldata_t pd,
                                        char **gt_brule, char **atype,
                                        char **geometry, int *numbonds,
                                        double *valence, int *iAromatic,
                                        char **neighbors);

extern void gmx_poldata_add_ptype(gmx_poldata_t pd,
                                  const char   *ptype,
                                  const char   *miller,
                                  const char   *bosque,
                                  double        polarizability,
                                  double        sig_pol);

extern void gmx_poldata_add_atype(gmx_poldata_t pd, const char *elem,
                                  const char *desc,
                                  const char *atype,
                                  const char *ptype,
                                  const char *btype,
                                  const char *vdwparams);

extern void gmx_poldata_set_ptype_polarizability(gmx_poldata_t pd, const char *ptype,
                                                 double polarizability, double sig_pol);

extern void gmx_poldata_set_force_field(gmx_poldata_t pd, const char *forcefield);

extern void gmx_poldata_set_polar_unit(gmx_poldata_t pd, const char *polar_unit);

extern void gmx_poldata_set_polar_ref(gmx_poldata_t pd, const char *polar_ref);

extern void gmx_poldata_set_length_unit(gmx_poldata_t pd, const char *length_unit);

extern void gmx_poldata_set_vdw_function(gmx_poldata_t pd, const char *func);

extern char *gmx_poldata_get_vdw_function(gmx_poldata_t pd);

extern int gmx_poldata_get_vdw_ftype(gmx_poldata_t pd);

extern void gmx_poldata_set_nexcl(gmx_poldata_t pd, int nexcl);

extern int gmx_poldata_get_nexcl(gmx_poldata_t pd);

extern void gmx_poldata_set_fudgeQQ(gmx_poldata_t pd, double fudgeQQ);

extern double gmx_poldata_get_fudgeQQ(gmx_poldata_t pd);

extern void gmx_poldata_set_fudgeLJ(gmx_poldata_t pd, double fudgeLJ);

extern double gmx_poldata_get_fudgeLJ(gmx_poldata_t pd);

extern void gmx_poldata_set_combination_rule(gmx_poldata_t pd, char *func);

extern char *gmx_poldata_get_combination_rule(gmx_poldata_t pd);

extern int gmx_poldata_get_comb_rule(gmx_poldata_t pd);

extern char *gmx_poldata_get_polar_unit(gmx_poldata_t pd);

extern char *gmx_poldata_get_polar_ref(gmx_poldata_t pd);

extern char *gmx_poldata_get_force_field(gmx_poldata_t pd);

extern char *gmx_poldata_get_length_unit(gmx_poldata_t pd);

/* Return array of atomtypes compatible with the bonded neighbors.
   The array should be freed, but not the contents of the elements.
 */
extern char **gmx_poldata_get_bonding_rules(gmx_poldata_t pd, char *elem,
                                            int nbond, char *neighbors[],
                                            const char *geometry,
                                            int iAromatic);

extern char *gmx_poldata_get_geometry(gmx_poldata_t pd, char *gt_brule);

extern char *gmx_poldata_get_desc(gmx_poldata_t pd, char *atype);

/* Get the charge from the gentop.dat file */
extern char *gmx_poldata_get_charge(gmx_poldata_t pd, char *atype);

/* Returns 1 if OK, 0 if last */
extern int gmx_poldata_get_atype(gmx_poldata_t pd,
                                 char        **elem,
                                 char        **desc,
                                 char        **atype,
                                 char        **ptype,
                                 char        **btype,
                                 char        **vdwparams);

extern int gmx_poldata_get_ptype(gmx_poldata_t pd,
                                 char        **ptype,
                                 char        **miller,
                                 char        **bosque,
                                 double       *polarizability,
                                 double       *sig_pol);

extern const char *gmx_poldata_atype_to_ptype(gmx_poldata_t pd, const char *atype);

extern const char *gmx_poldata_atype_to_btype(gmx_poldata_t pd, const char *atype);

/* Return 1 if OK, 0 if not found */
extern int gmx_poldata_search_atype(gmx_poldata_t pd,
                                    char         *key,
                                    char        **elem,
                                    char        **desc,
                                    char        **atype,
                                    char        **ptype,
                                    char        **btype,
                                    char        **vdwparams);

/* Return 1 if OK, 0 if not found */
extern int gmx_poldata_get_ptype_pol(gmx_poldata_t pd, const char *ptype,
                                     double *polarizability, double *sig_pol);
extern int gmx_poldata_get_atype_pol(gmx_poldata_t pd, const char *atype,
                                     double *polarizability, double *sig_pol);
/* Return 1 if OK, 0 if not found */
extern int gmx_poldata_bonding_rule_valence(gmx_poldata_t pd, char *gt_brule, double *valence);

extern void gmx_poldata_add_miller(gmx_poldata_t pd,
                                   char         *miller,
                                   int           atomnumber,
                                   double        tau_ahc,
                                   double        alpha_ahp);

/* Return 1 if "miller" was found */
extern int gmx_poldata_get_miller_pol(gmx_poldata_t pd,
                                      char         *miller,
                                      int          *atomnumber,
                                      double       *tau_ahc,
                                      double       *alpha_ahp);

extern int gmx_poldata_get_miller(gmx_poldata_t pd,
                                  char        **miller,
                                  int          *atomnumber,
                                  double       *tau_ahc,
                                  double       *alpha_ahp);

extern void gmx_poldata_set_miller_units(gmx_poldata_t pd, char *tau_unit,
                                         char *ahp_unit);

extern void gmx_poldata_get_miller_units(gmx_poldata_t pd, char **tau_unit,
                                         char **ahp_unit);

/* Returns miller name or NULL if not found */
extern char *gmx_poldata_ptype_to_miller(gmx_poldata_t pd, const char *ptype);

extern void gmx_poldata_add_bosque(gmx_poldata_t pd,
                                   char         *bosque,
                                   double        polarizability);

extern int gmx_poldata_get_bosque(gmx_poldata_t pd,
                                  char        **bosque,
                                  double       *polarizability);

extern void gmx_poldata_set_bosque_unit(gmx_poldata_t pd, char *polar_unit);

extern char *gmx_poldata_get_bosque_unit(gmx_poldata_t pd);

/* Returns bosque name or NULL if not found */
extern char *gmx_poldata_ptype_to_bosque(gmx_poldata_t pd, const char *ptype);

extern int gmx_poldata_get_bosque_pol(gmx_poldata_t pd, char *bosque, double *polarizability);

/* Return 1 on success or 0 otherwise */
extern int gmx_poldata_add_bond(gmx_poldata_t pd, char *atom1, char *atom2,
                                double length, double sigma, int ntrain,
                                double bondorder, char *params);

extern int gmx_poldata_set_bond_params(gmx_poldata_t pd, char *atom1, char *atom2,
                                       double length, double sigma, int ntrain,
                                       double bondorder, char *params);

/* Return bond-index 1-N or 0 if not found */
extern int gmx_poldata_get_bond(gmx_poldata_t pd, char **atom1, char **atom2,
                                double *length, double *sigma, int *ntrain,
                                double *bondorder, char **params);

extern void gmx_poldata_set_bond_function(gmx_poldata_t pd, char *fn);
extern char *gmx_poldata_get_bond_function(gmx_poldata_t pd);
extern int gmx_poldata_get_bond_ftype(gmx_poldata_t pd);

/* Return bond-index 1-N or 0 if not found */
extern int gmx_poldata_search_bond(gmx_poldata_t pd, char *atom1, char *atom2,
                                   double *length, double *sigma, int *ntrain,
                                   double *bondorder, char **params);

/* Returns 1 if there is a bond, 0 if not. Toler is absolute in length-units. */
extern int gmx_poldata_elem_is_bond(gmx_poldata_t pd, char *elem1, char *elem2,
                                    double distance, double toler);

/* Return maximal valence for a give element */
double gmx_poldata_elem_get_max_valence(gmx_poldata_t pd, char *elem);

/* Return NULL-terminated array of potential bondorders */
extern double *gmx_poldata_elem_get_bondorders(gmx_poldata_t pd, char *elem1, char *elem2,
                                               double distance, double toler);
/* Returns the bondorder. Toler is absolute in length-units. */
extern double gmx_poldata_atype_bondorder(gmx_poldata_t pd, char *atype1, char *atype2,
                                          double distance, double toler);

extern void gmx_poldata_set_angle_function(gmx_poldata_t pd, char *fn);
extern char *gmx_poldata_get_angle_function(gmx_poldata_t pd);
extern int gmx_poldata_get_angle_ftype(gmx_poldata_t pd);

/* Return 1 on success, 0 otherwise */
extern int gmx_poldata_add_angle(gmx_poldata_t pd, char *atom1, char *atom2,
                                 char *atom3, double angle, double sigma,
                                 int ntrain, char *params);

extern int gmx_poldata_set_angle_params(gmx_poldata_t pd, char *atom1, char *atom2,
                                        char *atom3,
                                        double angle, double sigma, int ntrain, char *params);

/* Return angle-index 1-N or 0 if not found */
extern int gmx_poldata_get_angle(gmx_poldata_t pd, char **atom1, char **atom2,
                                 char **atom3, double *angle, double *sigma,
                                 int *ntrain, char **params);

/* Return angle-index 1-N or 0 if not found */
extern int gmx_poldata_search_angle(gmx_poldata_t pd, char *atom1, char *atom2,
                                    char *atom3, double *angle, double *sigma,
                                    int *ntrain, char **params);

extern void gmx_poldata_set_angle_unit(gmx_poldata_t pd, char *angle_unit);

extern char *gmx_poldata_get_angle_unit(gmx_poldata_t pd);

extern void gmx_poldata_set_dihedral_function(gmx_poldata_t pd, int egd, char *fn);
extern char *gmx_poldata_get_dihedral_function(gmx_poldata_t pd, int egd);
extern int gmx_poldata_get_dihedral_ftype(gmx_poldata_t pd, int egd);

/* Return 1 on success or 0 otherwise */
extern int gmx_poldata_add_dihedral(gmx_poldata_t pd, int egd, char *atom1, char *atom2,
                                    char *atom3, char *atom4,
                                    double dihedral, double sigma,
                                    int ntrain, char *params);

extern int gmx_poldata_set_dihedral_params(gmx_poldata_t pd, int egd,
                                           char *atom1, char *atom2,
                                           char *atom3, char *atom4,
                                           double angle, double sigma,
                                           int ntrain, char *params);

/* Return dihedral-index 1-N or 0 if not found */
extern int gmx_poldata_get_dihedral(gmx_poldata_t pd, int egd,
                                    char **atom1, char **atom2,
                                    char **atom3, char **atom4,
                                    double *dihedral, double *sigma,
                                    int *ntrain, char **params);

/* Return dihedral-index 1-N or 0 if not found */
extern int gmx_poldata_search_dihedral(gmx_poldata_t pd, int egd,
                                       char *atom1, char *atom2,
                                       char *atom3, char *atom4,
                                       double *dihedral, double *sigma,
                                       int *ntrain, char **params);

extern void gmx_poldata_set_dihedral_unit(gmx_poldata_t pd, int egd,
                                          char *dihedral_unit);

extern char *gmx_poldata_get_dihedral_unit(gmx_poldata_t pd, int egd);

extern void gmx_poldata_add_symcharges(gmx_poldata_t pd, char *central,
                                       char *attached, int numattach);

extern int gmx_poldata_get_symcharges(gmx_poldata_t pd, char **central,
                                      char **attached, int *numattach);

extern int gmx_poldata_search_symcharges(gmx_poldata_t pd, char *central,
                                         char *attached, int numattach);

extern ChargeGenerationModel name2eemtype(const char *name);

extern const char *get_eemtype_name(ChargeGenerationModel eem);

extern char *gmx_poldata_get_eemref(gmx_poldata_t pd, ChargeGenerationModel eqg_model);

extern int gmx_poldata_get_numprops(gmx_poldata_t pd, ChargeGenerationModel eqg_model);

extern int gmx_poldata_have_pol_support(gmx_poldata_t pd, const char *atype);

extern int gmx_poldata_have_eem_support(gmx_poldata_t pd, ChargeGenerationModel eqg_model,
                                        const char *name,
                                        gmx_bool bAllowZeroParameters);

extern double gmx_poldata_get_j00(gmx_poldata_t pd, ChargeGenerationModel eqg_model, char *name);

extern int gmx_poldata_get_nzeta(gmx_poldata_t pd, ChargeGenerationModel eqg_model, char *name);

extern double gmx_poldata_get_zeta(gmx_poldata_t pd, ChargeGenerationModel eqg_model, char *name, int zz);

extern char *gmx_poldata_get_qstr(gmx_poldata_t pd, ChargeGenerationModel eqg_model, char *name);

extern char *gmx_poldata_get_rowstr(gmx_poldata_t pd, ChargeGenerationModel eqg_model, char *name);

extern double gmx_poldata_get_q(gmx_poldata_t pd, ChargeGenerationModel eqg_model, char *name, int zz);

extern int gmx_poldata_get_row(gmx_poldata_t pd, ChargeGenerationModel eqg_model, char *name, int zz);

extern double gmx_poldata_get_chi0(gmx_poldata_t pd, ChargeGenerationModel eqg_model, char *name);

extern char *gmx_poldata_get_opts(gmx_poldata_t pd, ChargeGenerationModel eqg_model, char *name);

extern void gmx_poldata_set_eemprops(gmx_poldata_t pd,
                                     ChargeGenerationModel eqg_model, char *name,
                                     double J0, double chi0,
                                     char *zeta, char *q, char *row);

extern int gmx_poldata_get_eemprops(gmx_poldata_t pd,
                                    ChargeGenerationModel *eqg_model, char **name,
                                    double *J0, double *chi0,
                                    char **zeta, char **q, char **row);

extern void gmx_poldata_set_epref(gmx_poldata_t pd, ChargeGenerationModel eqg_model, char *epref);

extern char *gmx_poldata_get_epref(gmx_poldata_t pd, ChargeGenerationModel eqg_model);

extern int gmx_poldata_list_epref(gmx_poldata_t pd, ChargeGenerationModel *eqg_model, char **epref);

extern void gmx_poldata_comm_eemprops(gmx_poldata_t pd, t_commrec *cr);

extern void gmx_poldata_comm_force_parameters(gmx_poldata_t pd, t_commrec *cr);

#endif
