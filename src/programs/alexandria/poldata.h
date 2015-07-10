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
#include "gmxpre.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/legacyheaders/network.h"
#include "gromacs/utility/smalloc.h"
#include "stringutil.h"
/* This source code file is part of the Alexandria project */

/*! \brief
 * Enumerated type holding the charge distribution models used in PolData
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */
enum ChargeDistributionModel {
    eqdAXp,
    eqdAXg,
    eqdAXs,
    eqdYang,
    eqdBultinck,
    eqdRappe,
    eqdNR
};

enum DihedralType {
    egdPDIHS,
    egdIDIHS,
    egdNR
};

typedef struct {
    char   *type, *miller, *bosque;
    double  polarizability, sig_pol;
} t_ptype;

typedef struct {
    char   *desc, *type, *ptype, *btype, *elem, *vdwparams;
    double  ref_enthalpy;
} t_ffatype;

typedef struct {
    char                    *elem, *rule, *type, *neighbors, *geometry;
    int                      numbonds, iAromatic;
    double                   valence;
    std::vector<std::string> nb;
} t_brule;

typedef struct {
    char   *atom1, *atom2, *params;
    char    elem1[4], elem2[4];
    double  length, sigma, bondorder;
    int     ntrain;
} t_gt_bond;

typedef struct {
    char   *atom1, *atom2, *atom3, *params;
    double  angle, sigma;
    int     ntrain;
} t_gt_angle;

typedef struct {
    char   *atom1, *atom2, *atom3, *atom4, *params;
    double  dihedral, sigma;
    int     ntrain;
} t_gt_dihedral;

typedef struct {
    char   *bosque;
    double  polarizability;
} t_bosque;

typedef struct {
    char   *miller;
    int     atomnumber;
    double  tau_ahc, alpha_ahp;
} t_miller;

typedef struct {
    char *central;
    char *attached;
    int   numattach;
} t_symcharges;

typedef struct {
    ChargeDistributionModel eqd_model;
    char                   *epref;
} t_epref;

#define EEMBUFSIZE 256
#define MAXZETA    12
typedef struct {
    ChargeDistributionModel eqd_model;
    int                     nzeta, row[MAXZETA];
    char                    name[EEMBUFSIZE], zetastr[EEMBUFSIZE], qstr[EEMBUFSIZE], rowstr[EEMBUFSIZE];
    double                  J0, chi0, q[MAXZETA], zeta[MAXZETA];
} t_eemprops;


namespace alexandria
{
  class Poldata
  {
  public:

    Poldata(); //Constructor
    ~Poldata(){}
    void  set_filename(char *fn2);


    int  get_natypes( );
    int  get_nptypes( );
    int  get_ngt_bond( );
    int  get_ngt_angle( );
    int  get_ngt_dihedral(  int egd);

    void  add_bonding_rule(
			   char *gt_brule, char *atype,
			   char *geometry, int numbonds,
			   double valence, int iAromatic,
			   char *neighbors);

    int  get_bonding_rule(
			  char **gt_brule, char **atype,
			  char **geometry, int *numbonds,
			  double *valence, int *iAromatic,
			  char **neighbors);

    void  add_ptype(
		    const char   *ptype,
		    const char   *miller,
		    const char   *bosque,
		    double        polarizability,
		    double        sig_pol);

    void  add_atype(  const char *elem,
		      const char *desc,
		      const char *atype,
		      const char *ptype,
		      const char *btype,
		      const char *vdwparams,
		      double      ref_enthalpy);

    void  set_ptype_polarizability(  const char *ptype,
				     double polarizability, double sig_pol);

    void  set_force_field(  const char *forcefield);

    void  set_polar_unit(  const char *polar_unit);

    void  set_polar_ref(  const char *polar_ref);

    void  set_length_unit(  const char *length_unit);

    void  set_vdw_function(  const char *func);

    char * get_vdw_function( );

    int  get_vdw_ftype( );

    void  set_nexcl(  int nexcl);

    int  get_nexcl( );

    void  set_fudgeQQ(  double fudgeQQ);

    double  get_fudgeQQ( );

    void  set_fudgeLJ(  double fudgeLJ);

    double  get_fudgeLJ( );

    int get_atype_ref_enthalpy( const char *atype,
				double     *Href);

    void  set_combination_rule(  char *func);

    char * get_combination_rule( );

    int  get_comb_rule( );

    char * get_polar_unit( );

    char * get_polar_ref( );

    char * get_force_field( );

    char * get_length_unit( );


    /* Return array of atomtypes compatible with the bonded neighbors.
       The array should be freed, but not the contents of the elements.
    */
    char ** get_bonding_rules(  char *elem,
				int nbond, char *neighbors[],
				const char *geometry,
				int iAromatic);

    char * get_geometry(  char *gt_brule);

    char * get_desc(  char *atype);

    /* Get the charge from the gentop.dat file */
    char * get_charge(  char *atype);

    /* Returns 1 if OK, 0 if last */
    int  get_atype(
		   char        **elem,
		   char        **desc,
		   char        **atype,
		   char        **ptype,
		   char        **btype,
		   char        **vdwparams,
		   double       *ref_enthalpy);

    int  get_ptype(
		   char        **ptype,
		   char        **miller,
		   char        **bosque,
		   double       *polarizability,
		   double       *sig_pol);

    const char * atype_to_ptype(  const char *atype);

    const char * atype_to_btype(  const char *atype);

    /* Return 1 if OK, 0 if not found */
    int  search_atype(
		      char         *key,
		      char        **elem,
		      char        **desc,
		      char        **atype,
		      char        **ptype,
		      char        **btype,
		      char        **vdwparams);

    /* Return 1 if OK, 0 if not found */
    int  get_ptype_pol(  const char *ptype,
			 double *polarizability, double *sig_pol);
    int  get_atype_pol(  const char *atype,
			 double *polarizability, double *sig_pol);

    int get_atype_ref_enthalpy(Poldata * pd, const char *atype,
			       double *Href);

    /* Return 1 if OK, 0 if not found */
    int  bonding_rule_valence(  char *gt_brule, double *valence);

    void  add_miller(
		     char         *miller,
		     int           atomnumber,
		     double        tau_ahc,
		     double        alpha_ahp);

    /* Return 1 if "miller" was found */
    int  get_miller_pol(
			char         *miller,
			int          *atomnumber,
			double       *tau_ahc,
			double       *alpha_ahp);

    int  get_miller(
		    char        **miller,
		    int          *atomnumber,
		    double       *tau_ahc,
		    double       *alpha_ahp);

    void  set_miller_units(  char *tau_unit,
			     char *ahp_unit);

    void  get_miller_units(  char **tau_unit,
			     char **ahp_unit);

    /* Returns miller name or NULL if not found */
    char * ptype_to_miller(  const char *ptype);

    void  add_bosque(
		     char         *bosque,
		     double        polarizability);

    int  get_bosque(
		    char        **bosque,
		    double       *polarizability);

    void  set_bosque_unit(  char *polar_unit);

    char * get_bosque_unit( );

    /* Returns bosque name or NULL if not found */
    char * ptype_to_bosque(  const char *ptype);

    int  get_bosque_pol(  char *bosque, double *polarizability);

    /* Return 1 on success or 0 otherwise */
    int  add_bond(  char *atom1, char *atom2,
		    double length, double sigma, int ntrain,
		    double bondorder, char *params);

    int  set_bond_params(  char *atom1, char *atom2,
			   double length, double sigma, int ntrain,
			   double bondorder, char *params);

    /* Return bond-index 1-N or 0 if not found */
    int  get_bond(  char **atom1, char **atom2,
		    double *length, double *sigma, int *ntrain,
		    double *bondorder, char **params);

    void  set_bond_function(  char *fn);
    char * get_bond_function( );
    int  get_bond_ftype( );

    /* Return bond-index 1-N or 0 if not found */
    int  search_bond(  char *atom1, char *atom2,
		       double *length, double *sigma, int *ntrain,
		       double *bondorder, char **params);

    /* Returns 1 if there is a bond, 0 if not. Toler is absolute in length-units. */
    int  elem_is_bond(  char *elem1, char *elem2,
			double distance, double toler);

    /* Return maximal valence for a give element */
    double  elem_get_max_valence(  char *elem);

    /* Return NULL-terminated array of potential bondorders */
    double * elem_get_bondorders(  char *elem1, char *elem2,
				   double distance, double toler);
    /* Returns the bondorder. Toler is absolute in length-units. */
    double  atype_bondorder(  char *atype1, char *atype2,
			      double distance, double toler);

    void  set_angle_function(  char *fn);
    char * get_angle_function( );
    int  get_angle_ftype( );

    /* Return 1 on success, 0 otherwise */
    int  add_angle(
		   char *atom1, char *atom2,
		   char *atom3, double angle, double sigma,
		   int ntrain, char *params);

    int  set_angle_params(  char *atom1, char *atom2,
			    char *atom3,
			    double angle, double sigma, int ntrain, char *params);

    /* Return angle-index 1-N or 0 if not found */
    int  get_angle(  char **atom1, char **atom2,
		     char **atom3, double *angle, double *sigma,
		     int *ntrain, char **params);

    /* Return angle-index 1-N or 0 if not found */
    int  search_angle(  char *atom1, char *atom2,
			char *atom3, double *angle, double *sigma,
			int *ntrain, char **params);

    void  set_angle_unit(  char *angle_unit);

    char * get_angle_unit();

    void  set_dihedral_function(  int egd, char *fn);
    char * get_dihedral_function(  int egd);
    int  get_dihedral_ftype(  int egd);

    /* Return 1 on success or 0 otherwise */
    int  add_dihedral(  int egd, char *atom1, char *atom2,
			char *atom3, char *atom4,
			double dihedral, double sigma,
			int ntrain, char *params);

    int  set_dihedral_params(  int egd,
			       char *atom1, char *atom2,
			       char *atom3, char *atom4,
			       double angle, double sigma,
			       int ntrain, char *params);

    /* Return dihedral-index 1-N or 0 if not found */
    int  get_dihedral(  int egd,
			char **atom1, char **atom2,
			char **atom3, char **atom4,
			double *dihedral, double *sigma,
			int *ntrain, char **params);

    /* Return dihedral-index 1-N or 0 if not found */
    int  search_dihedral(  int egd,
			   char *atom1, char *atom2,
			   char *atom3, char *atom4,
			   double *dihedral, double *sigma,
			   int *ntrain, char **params);

    void  set_dihedral_unit(  int   egd,
			      char *dihedral_unit);

    char * get_dihedral_unit(  int egd);

    void  add_symcharges(  char *central,
			   char *attached, int numattach);

    int  get_symcharges(  char **central,
			  char **attached, int *numattach);

    int  search_symcharges(  char *central,
			     char *attached, int numattach);

    static ChargeDistributionModel name2eemtype(const char *name);

    static const char *get_eemtype_name(ChargeDistributionModel eem);

    char * get_eemref(  ChargeDistributionModel eqd_model);

    int  get_numprops(  ChargeDistributionModel eqd_model);

    int  have_pol_support(  const char *atype);

    int  have_eem_support(  ChargeDistributionModel eqd_model,
			    const char             *name,
			    gmx_bool                bAllowZeroParameters);

    double  get_j00(  ChargeDistributionModel eqd_model, char *name);

    int  get_nzeta(  ChargeDistributionModel eqd_model, char *name);

    double  get_zeta(  ChargeDistributionModel eqd_model, char *name, int zz);

    char * get_qstr(  ChargeDistributionModel eqd_model, char *name);

    char * get_rowstr(  ChargeDistributionModel eqd_model, char *name);

    double  get_q(  ChargeDistributionModel eqd_model, char *name, int zz);

    int  get_row(  ChargeDistributionModel eqd_model, char *name, int zz);

    double  get_chi0(  ChargeDistributionModel eqd_model, char *name);

    char * get_opts(  ChargeDistributionModel eqd_model, char *name);

    void  set_eemprops(
		       ChargeDistributionModel eqd_model, char *name,
		       double J0, double chi0,
		       char *zeta, char *q, char *row);

    int  get_eemprops(
		      ChargeDistributionModel *eqd_model, char **name,
		      double *J0, double *chi0,
		      char **zeta, char **q, char **row);

    void  set_epref(  ChargeDistributionModel eqd_model, char *epref);

    char * get_epref(ChargeDistributionModel eqd_model);

    int list_epref(  ChargeDistributionModel *eqd_model, char **epref);

    void  comm_eemprops(  t_commrec *cr);

    void  comm_force_parameters( t_commrec *cr);


  private:

    char          *filename = NULL;
    int            nptype = 0, nptype_c = 0;
    t_ptype       *ptype = NULL;
    int            nalexandria = 0, nalexandria_c = 0;
    t_ffatype     *alexandria = NULL;
    int            nbtype = 0;
    char         **btype = NULL;
    int            nbrule = 0, nbrule_c = 0;
    t_brule       *brule = NULL;
    char          *alexandria_polar_unit = NULL;
    char          *alexandria_polar_ref = NULL;
    char          *alexandria_forcefield = NULL;
    int            nexcl = 0;
    double         fudgeQQ = 0, fudgeLJ = 0;
    char          *gt_vdw_function = NULL, *gt_combination_rule = NULL;
    int            gt_vdw_ftype = 0, gt_comb_rule = 0;
    char          *gt_bond_function = NULL;
    int            ngt_bond = 0, ngt_bond_c = 0, gt_bond_ftype = 0;
    char          *gt_length_unit = NULL;
    t_gt_bond     *gt_bond = NULL;
    char          *gt_angle_function = NULL;
    int            ngt_angle = 0, ngt_angle_c = 0, gt_angle_ftype = 0;
    char          *gt_angle_unit = NULL;
    t_gt_angle    *gt_angle = NULL;
    char          *gt_dihedral_function[egdNR];
    int            ngt_dihedral[egdNR], ngt_dihedral_c[egdNR], gt_dihedral_ftype[egdNR];
    t_gt_dihedral *gt_dihedral[egdNR];
    int            nmiller = 0, nmiller_c = 0;
    t_miller      *miller = NULL;
    char          *miller_tau_unit = NULL, *miller_ahp_unit = NULL;
    int            nbosque = 0, nbosque_c = 0;
    t_bosque      *bosque = NULL;
    char          *bosque_polar_unit = NULL;
    int            nsymcharges = 0, nsymcharges_c = 0;
    t_symcharges  *symcharges = NULL;
    int            nep = 0, nep_c = 0;
    t_eemprops    *eep = NULL;
    int            ner = 0, ner_c = 0;
    t_epref       *epr = NULL;

    void add_btype(const char   *btype);

    gmx_bool strcasestr_start(char *needle, char *haystack);

    int count_neighbors(t_brule *brule, int nbond, char *nbhybrid[], int *score);

    t_gt_bond *search_bond(  char *atom1, char *atom2,
			     double bondorder);

    int search_bondtype(  char *atom);

    t_eemprops *get_eep(ChargeDistributionModel eqd_model,
			const char             *name);

    t_gt_dihedral *search_dihedral( int egd, char *atom1, char *atom2,
				    char *atom3, char *atom4);

    static int gtb_comp(const void *a, const void *b);

  };
}
#endif
