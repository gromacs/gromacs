/*
 * This source file is part of the Aleandria project.
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
#include "gmxpre.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <string.h>
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/legacyheaders/network.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/legacyheaders/names.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/cstringutil.h"
#include "poldata.h"
#include "gmx_simple_comm.h"
#include "stringutil.h"

typedef struct {
    char   *type, *miller, *bosque;
    double  polarizability, sig_pol;
} t_ptype;

typedef struct {
    char   *desc, *type, *ptype, *btype, *elem, *vdwparams;
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

#define EEMBUFSIZE 256
#define MAXZETA    12
typedef struct {
    ChargeGenerationModel eqg_model;
    int                   nzeta, row[MAXZETA];
    char                  name[EEMBUFSIZE], zetastr[EEMBUFSIZE], qstr[EEMBUFSIZE], rowstr[EEMBUFSIZE];
    double                J0, chi0, q[MAXZETA], zeta[MAXZETA];
} t_eemprops;

typedef struct {
    ChargeGenerationModel eqg_model;
    char                 *epref;
} t_epref;

typedef struct gmx_poldata {
    char          *filename;
    int            nptype, nptype_c;
    t_ptype       *ptype;
    int            nalexandria, nalexandria_c;
    t_ffatype     *alexandria;
    int            nbtype;
    char         **btype;
    int            nbrule, nbrule_c;
    t_brule       *brule;
    char          *alexandria_polar_unit;
    char          *alexandria_polar_ref;
    char          *alexandria_forcefield;
    int            nexcl;
    double         fudgeQQ, fudgeLJ;
    char          *gt_vdw_function, *gt_combination_rule;
    int            gt_vdw_ftype, gt_comb_rule;
    char          *gt_bond_function;
    int            ngt_bond, ngt_bond_c, gt_bond_ftype;
    char          *gt_length_unit;
    t_gt_bond     *gt_bond;
    char          *gt_angle_function;
    int            ngt_angle, ngt_angle_c, gt_angle_ftype;
    char          *gt_angle_unit;
    t_gt_angle    *gt_angle;
    char          *gt_dihedral_function[egdNR];
    int            ngt_dihedral[egdNR], ngt_dihedral_c[egdNR], gt_dihedral_ftype[egdNR];
    t_gt_dihedral *gt_dihedral[egdNR];
    int            nmiller, nmiller_c;
    t_miller      *miller;
    char          *miller_tau_unit, *miller_ahp_unit;
    int            nbosque, nbosque_c;
    t_bosque      *bosque;
    char          *bosque_polar_unit;
    int            nsymcharges, nsymcharges_c;
    t_symcharges  *symcharges;
    int            nep, nep_c;
    t_eemprops    *eep;
    int            ner, ner_c;
    t_epref       *epr;
} gmx_poldata;

#define assign_str(dst, src)  if (dst) { if (src) {*dst = strdup(src); }else{*dst = NULL; }}
#define assign_scal(dst, src) if (dst) *dst = src

gmx_poldata_t gmx_poldata_init()
{
    gmx_poldata_t pd;
    int           i;

    snew(pd, 1);

    /* Initiate some crucial variables */
    pd->nexcl          = NOTSET;
    pd->fudgeQQ        = NOTSET;
    pd->fudgeLJ        = NOTSET;
    pd->gt_comb_rule   = NOTSET;
    pd->gt_vdw_ftype   = NOTSET;
    pd->gt_bond_ftype  = NOTSET;
    pd->gt_angle_ftype = NOTSET;
    for (i = 0; (i < egdNR); i++)
    {
        pd->gt_dihedral_ftype[i] = NOTSET;
    }

    return pd;
}

void gmx_poldata_set_filename(gmx_poldata_t pd, char *fn2)
{
    if (NULL == fn2)
    {
        fprintf(stderr, "Trying to set poldata filename to NULL\n");
        return;
    }
    if (NULL != pd->filename)
    {
        fprintf(stderr, "Overwriting poldata filename from %s to %s\n",
                pd->filename, fn2);
        sfree(pd->filename);
    }
    pd->filename = strdup(fn2);
}

void gmx_poldata_set_vdw_function(gmx_poldata_t pd, const char *func)
{
    int i;

    if (NULL != pd->gt_vdw_function)
    {
        sfree(pd->gt_vdw_function);
    }
    for (i = 0; (i < F_NRE); i++)
    {
        if (strcasecmp(interaction_function[i].name, func) == 0)
        {
            break;
        }
    }
    if (i == F_NRE)
    {
        gmx_fatal(FARGS, "Van der Waals function '%s' does not exist in gromacs", func);
    }
    pd->gt_vdw_ftype    = i;
    pd->gt_vdw_function = strdup(func);
}

char *gmx_poldata_get_vdw_function(gmx_poldata_t pd)
{
    return pd->gt_vdw_function;
}

int gmx_poldata_get_vdw_ftype(gmx_poldata_t pd)
{
    if (NOTSET == pd->gt_vdw_ftype)
    {
        gmx_fatal(FARGS, "No valid function type for Van der Waals function in %s",
                  pd->filename);
    }
    return pd->gt_vdw_ftype;
}

void gmx_poldata_set_nexcl(gmx_poldata_t pd, int nexcl)
{
    pd->nexcl = nexcl;
}

int gmx_poldata_get_nexcl(gmx_poldata_t pd)
{
    if (NOTSET == pd->nexcl)
    {
        gmx_fatal(FARGS, "Nexclusions not set in %s", pd->filename);
    }
    return pd->nexcl;
}

void gmx_poldata_set_fudgeQQ(gmx_poldata_t pd, double fudgeQQ)
{
    pd->fudgeQQ = fudgeQQ;
}

double gmx_poldata_get_fudgeQQ(gmx_poldata_t pd)
{
    if (NOTSET == pd->fudgeQQ)
    {
        gmx_fatal(FARGS, "fudgeQQ not set in %s", pd->filename);
    }
    return pd->fudgeQQ;
}

void gmx_poldata_set_fudgeLJ(gmx_poldata_t pd, double fudgeLJ)
{
    pd->fudgeLJ = fudgeLJ;
}

double gmx_poldata_get_fudgeLJ(gmx_poldata_t pd)
{
    if (NOTSET == pd->fudgeLJ)
    {
        gmx_fatal(FARGS, "fudgeLJ not set in %s", pd->filename);
    }
    return pd->fudgeLJ;
}

void gmx_poldata_set_combination_rule(gmx_poldata_t pd, char *func)
{
    int i;

    if (NULL != pd->gt_combination_rule)
    {
        sfree(pd->gt_combination_rule);
    }
    for (i = 0; (i < eCOMB_NR); i++)
    {
        if (strcasecmp(ecomb_names[i], func) == 0)
        {
            break;
        }
    }
    if (i == eCOMB_NR)
    {
        gmx_fatal(FARGS, "Combination rule '%s' does not exist in gromacs", func);
    }
    pd->gt_comb_rule        = i;
    pd->gt_combination_rule = strdup(func);
}

char *gmx_poldata_get_combination_rule(gmx_poldata_t pd)
{
    if (NOTSET == pd->gt_vdw_ftype)
    {
        gmx_fatal(FARGS, "No valid function type for Van der Waals function in %s",
                  pd->filename);
    }
    return pd->gt_combination_rule;
}

int gmx_poldata_get_comb_rule(gmx_poldata_t pd)
{
    if (NOTSET == pd->gt_comb_rule)
    {
        gmx_fatal(FARGS, "No valid combination rule in %s", pd->filename);
    }
    return pd->gt_comb_rule;
}

void gmx_poldata_add_ptype(gmx_poldata_t pd,
                           const char   *ptype,
                           const char   *miller,
                           const char   *bosque,
                           double        polarizability,
                           double        sig_pol)
{
    t_ptype *sp;
    int      i;

    for (i = 0; (i < pd->nptype); i++)
    {
        if (strcmp(pd->ptype[i].type, ptype) == 0)
        {
            break;
        }
    }
    if (i == pd->nptype)
    {
        pd->nptype++;
        srenew(pd->ptype, pd->nptype);

        sp                 = &(pd->ptype[i]);
        sp->type           = strdup(ptype);
        sp->bosque         = strdup(bosque);
        sp->miller         = strdup(miller);
        sp->polarizability = polarizability;
        sp->sig_pol        = sig_pol;
    }
    else
    {
        fprintf(stderr, "Polarizability type %s was already added to poldata record\n", ptype);
    }
}

static void gmx_poldata_add_btype(gmx_poldata_t pd,
                                  const char   *btype)
{
    int i;

    for (i = 0; (i < pd->nbtype); i++)
    {
        if (strcmp(btype, pd->btype[i]) == 0)
        {
            break;
        }
    }
    if (i == pd->nbtype)
    {
        srenew(pd->btype, ++pd->nbtype);
        pd->btype[i] = strdup(btype);
    }
}

void gmx_poldata_add_atype(gmx_poldata_t pd,
                           const char   *elem,
                           const char   *desc,
                           const char   *atype,
                           const char   *ptype,
                           const char   *btype,
                           const char   *vdwparams)
{
    t_ffatype *sp;
    int        i;

    for (i = 0; (i < pd->nalexandria); i++)
    {
        if (strcmp(pd->alexandria[i].type, atype) == 0)
        {
            break;
        }
    }
    if (i == pd->nalexandria)
    {
        pd->nalexandria++;
        srenew(pd->alexandria, pd->nalexandria);

        sp                 = &(pd->alexandria[i]);
        sp->elem           = strdup(elem);
        sp->desc           = strdup(desc);
        sp->type           = strdup(atype);
        sp->ptype          = strdup(ptype);
        sp->btype          = strdup(btype);
        sp->vdwparams      = strdup(vdwparams);
        gmx_poldata_add_btype(pd, btype);
    }
    else
    {
        fprintf(stderr, "Atom type %s was already added to poldata record\n", atype);
    }
}

void gmx_poldata_add_bonding_rule(gmx_poldata_t pd,
                                  char *gt_brule, char *atype,
                                  char *geometry, int numbonds,
                                  double valence, int iAromatic,
                                  char *neighbors)
{
    t_brule *sp;
    int      i, j;

    for (j = 0; (j < pd->nalexandria); j++)
    {
        if (strcmp(pd->alexandria[j].type, atype) == 0)
        {
            break;
        }
    }
    if (j < pd->nalexandria)
    {
        for (i = 0; (i < pd->nbrule); i++)
        {
            if (strcmp(pd->brule[i].rule, gt_brule) == 0)
            {
                break;
            }
        }
        if (i == pd->nbrule)
        {
            pd->nbrule++;
            srenew(pd->brule, pd->nbrule);

            sp                 = &(pd->brule[i]);
            sp->elem           = strdup(pd->alexandria[j].elem);
            sp->rule           = strdup(gt_brule);
            sp->type           = strdup(atype);
            sp->neighbors      = strdup(neighbors);
            sp->valence        = valence;
            sp->iAromatic      = iAromatic;
            sp->nb             = split(neighbors, ' ');
            sp->geometry       = strdup(geometry);
            sp->numbonds       = numbonds;
        }
        else
        {
            fprintf(stderr, "Bonding rule %s was already added to poldata record\n", gt_brule);
        }
    }
    else if (NULL != debug)
    {
        fprintf(debug, "Ignoring bonding rule involving unknown atom type %s\n",
                atype);
    }
}

int gmx_poldata_get_bonding_rule(gmx_poldata_t pd,
                                 char **gt_brule, char **atype,
                                 char **geometry, int *numbonds,
                                 double *valence, int *iAromatic,
                                 char **neighbors)
{
    if (pd->nbrule_c < pd->nbrule)
    {
        assign_str(gt_brule, pd->brule[pd->nbrule_c].rule);
        assign_str(atype, pd->brule[pd->nbrule_c].type);
        assign_str(geometry, pd->brule[pd->nbrule_c].geometry);
        assign_scal(numbonds, pd->brule[pd->nbrule_c].numbonds);
        assign_scal(valence, pd->brule[pd->nbrule_c].valence);
        assign_scal(iAromatic, pd->brule[pd->nbrule_c].iAromatic);
        assign_str(neighbors, pd->brule[pd->nbrule_c].neighbors);

        pd->nbrule_c++;

        return 1;
    }
    return 0;
}

int gmx_poldata_get_natypes(gmx_poldata_t pd)
{
    return pd->nalexandria;
}

int gmx_poldata_get_nptypes(gmx_poldata_t pd)
{
    return pd->nptype;
}

int gmx_poldata_get_ngt_bond(gmx_poldata_t pd)
{
    return pd->ngt_bond;
}

int gmx_poldata_get_ngt_angle(gmx_poldata_t pd)
{
    return pd->ngt_angle;
}

int gmx_poldata_get_ngt_dihedral(gmx_poldata_t pd, int egd)
{
    return pd->ngt_dihedral[egd];
}

void gmx_poldata_set_bond_function(gmx_poldata_t pd, char *fn)
{
    int i;

    if (NULL != pd->gt_bond_function)
    {
        sfree(pd->gt_bond_function);
    }
    for (i = 0; (i < F_NRE); i++)
    {
        if (strcasecmp(interaction_function[i].name, fn) == 0)
        {
            break;
        }
    }
    if (i == F_NRE)
    {
        gmx_fatal(FARGS, "Bond function '%s' does not exist in gromacs", fn);
    }
    pd->gt_bond_ftype    = i;
    pd->gt_bond_function = strdup(fn);
}

char *gmx_poldata_get_bond_function(gmx_poldata_t pd)
{
    return pd->gt_bond_function;
}

void gmx_poldata_set_angle_function(gmx_poldata_t pd, char *fn)
{
    int i;

    if (NULL != pd->gt_angle_function)
    {
        sfree(pd->gt_angle_function);
    }
    for (i = 0; (i < F_NRE); i++)
    {
        if (strcasecmp(interaction_function[i].name, fn) == 0)
        {
            break;
        }
    }
    if (i == F_NRE)
    {
        gmx_fatal(FARGS, "Angle function '%s' does not exist in gromacs", fn);
    }
    pd->gt_angle_ftype    = i;
    pd->gt_angle_function = strdup(fn);
}

char *gmx_poldata_get_angle_function(gmx_poldata_t pd)
{
    return pd->gt_angle_function;
}

int gmx_poldata_get_bond_ftype(gmx_poldata_t pd)
{
    return pd->gt_bond_ftype;
}

int gmx_poldata_get_angle_ftype(gmx_poldata_t pd)
{
    return pd->gt_angle_ftype;
}

int gmx_poldata_get_dihedral_ftype(gmx_poldata_t pd, int egd)
{
    return pd->gt_dihedral_ftype[egd];
}

void gmx_poldata_set_dihedral_function(gmx_poldata_t pd, int egd, char *func)
{
    int i;

    if (NULL != pd->gt_dihedral_function[egd])
    {
        sfree(pd->gt_dihedral_function[egd]);
    }
    for (i = 0; (i < F_NRE); i++)
    {
        if (strcasecmp(interaction_function[i].name, func) == 0)
        {
            break;
        }
    }
    if (i == F_NRE)
    {
        gmx_fatal(FARGS, "Dihedral function '%s' does not exist in gromacs", func);
    }
    pd->gt_dihedral_ftype[egd]    = i;
    pd->gt_dihedral_function[egd] = strdup(func);
}

char *gmx_poldata_get_dihedral_function(gmx_poldata_t pd, int egd)
{
    return pd->gt_dihedral_function[egd];
}

void gmx_poldata_set_ptype_polarizability(gmx_poldata_t pd, const char *ptype,
                                          double polarizability, double sig_pol)
{
    t_ptype *sp;
    int      i;

    for (i = 0; (i < pd->nptype); i++)
    {
        if (strcmp(ptype, pd->ptype[i].type) == 0)
        {
            sp                 = &(pd->ptype[i]);
            sp->polarizability = polarizability;
            sp->sig_pol        = sig_pol;
            break;
        }
    }
    if (i == pd->nptype)
    {
        fprintf(stderr, "No such ptype %s when trying to set the polarizability.\n", ptype);
    }
}

char *gmx_poldata_get_polar_unit(gmx_poldata_t pd)
{
    return pd->alexandria_polar_unit;
}

char *gmx_poldata_get_polar_ref(gmx_poldata_t pd)
{
    return pd->alexandria_polar_ref;
}

char *gmx_poldata_get_length_unit(gmx_poldata_t pd)
{
    if (NULL == pd->gt_length_unit)
    {
        gmx_fatal(FARGS, "No length unit in %s",
                  (NULL != pd->filename) ? pd->filename : "unknown");
    }

    return pd->gt_length_unit;
}

void gmx_poldata_set_polar_unit(gmx_poldata_t pd, const char *polar_unit)
{
    pd->alexandria_polar_unit   = strdup(polar_unit);
}

void gmx_poldata_set_polar_ref(gmx_poldata_t pd, const char *polar_ref)
{
    pd->alexandria_polar_ref   = strdup(polar_ref);
}

char *gmx_poldata_get_force_field(gmx_poldata_t pd)
{
    return pd->alexandria_forcefield;
}

void gmx_poldata_set_force_field(gmx_poldata_t pd, const char *forcefield)
{
    pd->alexandria_forcefield   = strdup(forcefield);
}

void gmx_poldata_set_length_unit(gmx_poldata_t pd, const char *length_unit)
{
    pd->gt_length_unit   = strdup(length_unit);
}

char *gmx_poldata_get_geometry(gmx_poldata_t pd, char *gt_brule)
{
    int i;

    if (gt_brule)
    {
        for (i = 0; (i < pd->nbrule); i++)
        {
            if (strcmp(pd->brule[i].rule, gt_brule) == 0)
            {
                return pd->brule[i].geometry;
            }
        }
    }

    return NULL;
}

char *gmx_poldata_get_desc(gmx_poldata_t pd, char *atype)
{
    int i;

    if (atype)
    {
        for (i = 0; (i < pd->nalexandria); i++)
        {
            if (strcmp(pd->alexandria[i].type, atype) == 0)
            {
                return pd->alexandria[i].desc;
            }
        }
    }

    return NULL;
}

char *gmx_poldata_get_type(gmx_poldata_t pd, char *gt_brule)
{
    int i;

    if (gt_brule)
    {
        for (i = 0; (i < pd->nbrule); i++)
        {
            if (strcmp(pd->brule[i].rule, gt_brule) == 0)
            {
                return pd->brule[i].type;
            }
        }
    }

    return NULL;
}

static gmx_bool strcasestr_start(char *needle, char *haystack)
{
    char *ptr;

    ptr = strcasestr(haystack, needle);

    return (ptr == haystack);
}

static int count_neighbors(t_brule *brule, int nbond, char *nbhybrid[], int *score)
{
    int j, ni = 0, *jj, i_found;

    *score = 0;
    snew(jj, nbond+1);
    for (unsigned int i = 0; (i < brule->nb.size()); i++)
    {
        i_found = 0;
        for (j = 0; (j < nbond); j++)
        {
            if (
                (NULL != nbhybrid[j]) &&
                (jj[j] == 0) &&
                (i_found == 0) &&
                strcasestr_start((char *)brule->nb[i].c_str(), nbhybrid[j])
                )
            {
                i_found = 1;
                jj[j]   = j+1;
                ni++;
                (*score) += 1;
                if (strlen(nbhybrid[j]) > 0)
                {
                    (*score) += 1;
                }
            }
        }
    }
    sfree(jj);

    return ni;
}

char **gmx_poldata_get_bonding_rules(gmx_poldata_t pd, char *elem,
                                     int nbond, char *neighbors[],
                                     const char *geometry,
                                     int iAromatic)
{
    unsigned int    nnb;
    int             i, nptr = 0, best = -1, score;
    char          **ptr = NULL;

    for (i = 0; (i < pd->nbrule); i++)
    {
        nnb = count_neighbors(&(pd->brule[i]), nbond, neighbors, &score);
        if ((strcmp(pd->brule[i].elem, elem) == 0) &&
            (strcasecmp(pd->brule[i].geometry, geometry) == 0) &&
            (nbond == pd->brule[i].numbonds) &&
            ((iAromatic >= 0 && (iAromatic == pd->brule[i].iAromatic)) ||
             (iAromatic < 0)) &&
            (nnb == pd->brule[i].nb.size()))
        {
            if (score > best)
            {
                if (NULL != ptr)
                {
                    sfree(ptr);
                    nptr = 0;
                    ptr  = NULL;
                }
            }
            if (score >= best)
            {
                srenew(ptr, ++nptr);
                ptr[nptr-1] = pd->brule[i].rule;
                best        = score;
            }
        }
    }
    srenew(ptr, ++nptr);
    ptr[nptr-1] = NULL;

    return ptr;
}

int gmx_poldata_bonding_rule_valence(gmx_poldata_t pd, char *gt_brule, double *valence)
{
    int i;

    for (i = 0; (i < pd->nbrule); i++)
    {
        if (strcasecmp(gt_brule, pd->brule[i].rule) == 0)
        {
            *valence = pd->brule[i].valence;
            return 1;
        }
    }
    return 0;
}

int gmx_poldata_get_ptype_pol(gmx_poldata_t pd, const char *ptype,
                              double *polar, double *sig_pol)
{
    int j;

    for (j = 0; (j < pd->nptype); j++)
    {
        if (strcmp(ptype, pd->ptype[j].type) == 0)
        {
            *polar   = pd->ptype[j].polarizability;
            *sig_pol = pd->ptype[j].sig_pol;
            return 1;
        }
    }
    return 0;
}

int gmx_poldata_get_atype_pol(gmx_poldata_t pd, const char *atype,
                              double *polar, double *sig_pol)
{
    int i;

    for (i = 0; (i < pd->nalexandria); i++)
    {
        if (strcmp(atype, pd->alexandria[i].type) == 0)
        {
            return gmx_poldata_get_ptype_pol(pd, pd->alexandria[i].ptype, polar, sig_pol);
        }
    }
    return 0;
}

char *gmx_poldata_ptype_to_miller(gmx_poldata_t pd, const char *ptype)
{
    int i;

    for (i = 0; (i < pd->nptype); i++)
    {
        if (strcmp(ptype, pd->ptype[i].type) == 0)
        {
            return pd->ptype[i].miller;
        }
    }
    return NULL;
}

char *gmx_poldata_ptype_to_bosque(gmx_poldata_t pd, const char *ptype)
{
    int i;

    for (i = 0; (i < pd->nptype); i++)
    {
        if (strcmp(ptype, pd->ptype[i].type) == 0)
        {
            return pd->ptype[i].bosque;
        }
    }
    return NULL;
}

int gmx_poldata_get_ptype(gmx_poldata_t pd,
                          char        **ptype,
                          char        **miller,
                          char        **bosque,
                          double       *polarizability,
                          double       *sig_pol)
{
    t_ptype *sp;

    if (pd->nptype_c < pd->nptype)
    {
        sp = &(pd->ptype[pd->nptype_c]);
        assign_scal(polarizability, sp->polarizability);
        assign_scal(sig_pol, sp->sig_pol);
        assign_str(ptype, sp->type);
        assign_str(miller, sp->miller);
        assign_str(bosque, sp->bosque);
        pd->nptype_c++;
        return 1;
    }
    else
    {
        pd->nptype_c = 0;
    }

    return 0;
}

int gmx_poldata_get_atype(gmx_poldata_t pd,
                          char        **elem,
                          char        **desc,
                          char        **atype,
                          char        **ptype,
                          char        **btype,
                          char        **vdwparams)
{
    t_ffatype *sp;

    if (pd->nalexandria_c < pd->nalexandria)
    {
        sp = &(pd->alexandria[pd->nalexandria_c]);
        assign_str(elem, sp->elem);
        assign_str(desc, sp->desc);
        assign_str(atype, sp->type);
        assign_str(ptype, sp->ptype);
        assign_str(btype, sp->btype);
        assign_str(vdwparams, sp->vdwparams);
        pd->nalexandria_c++;
        return 1;
    }
    else
    {
        pd->nalexandria_c = 0;
    }

    return 0;
}

const char *gmx_poldata_atype_to_ptype(gmx_poldata_t pd, const char *atype)
{
    int i;

    for (i = 0; (i < pd->nalexandria); i++)
    {
        if (strcmp(pd->alexandria[i].type, atype) == 0)
        {
            return pd->alexandria[i].ptype;
        }
    }
    return NULL;
}

const char *gmx_poldata_atype_to_btype(gmx_poldata_t pd, const char *atype)
{
    int i;

    for (i = 0; (i < pd->nalexandria); i++)
    {
        if (strcmp(pd->alexandria[i].type, atype) == 0)
        {
            return pd->alexandria[i].btype;
        }
    }
    return NULL;
}

int gmx_poldata_search_atype(gmx_poldata_t pd,
                             char         *key,
                             char        **elem,
                             char        **desc,
                             char        **atype,
                             char        **ptype,
                             char        **btype,
                             char        **vdwparams)
{
    t_ffatype *sp;
    int        i;

    for (i = 0; (i < pd->nalexandria); i++)
    {
        if (strcmp(key, pd->alexandria[i].type) == 0)
        {
            break;
        }
    }

    if (i < pd->nalexandria)
    {
        sp = &(pd->alexandria[i]);
        assign_str(elem, sp->elem);
        assign_str(desc, sp->desc);
        assign_str(atype, sp->type);
        assign_str(ptype, sp->ptype);
        assign_str(btype, sp->btype);
        assign_str(vdwparams, sp->vdwparams);

        return 1;
    }
    else
    {
        return 0;
    }
}

double gmx_poldata_elem_get_max_valence(gmx_poldata_t pd, char *elem)
{
    double mv = 0;
    int    i;

    for (i = 0; (i < pd->nbrule); i++)
    {
        if ((0 == gmx_strcasecmp(pd->brule[i].elem, elem)) &&
            (mv < pd->brule[i].valence))
        {
            mv = pd->brule[i].valence;
        }
    }
    return mv;
}

double *gmx_poldata_elem_get_bondorders(gmx_poldata_t pd, char *elem1, char *elem2,
                                        double distance, double toler)
{
    double  dev, *bo = NULL;
    char   *ba1, *ba2;
    int     i, j, k, nbo;

    if ((NULL == elem1) || (NULL == elem2))
    {
        return 0;
    }
    nbo = 0;
    for (i = 0; (i < pd->ngt_bond); i++)
    {
        if (0 == strlen(pd->gt_bond[i].elem1))
        {
            for (j = 0; (j < pd->nalexandria); j++)
            {
                if (strcmp(pd->alexandria[j].type, pd->gt_bond[i].atom2) == 0)
                {
                    strcpy(pd->gt_bond[i].elem1, pd->alexandria[j].elem);
                }
            }
        }
        if (0 == strlen(pd->gt_bond[i].elem2))
        {
            for (j = 0; (j < pd->nalexandria); j++)
            {
                if (strcmp(pd->alexandria[j].type, pd->gt_bond[i].atom2) == 0)
                {
                    strcpy(pd->gt_bond[i].elem2, pd->alexandria[j].elem);
                }
            }
        }
        ba1 = pd->gt_bond[i].elem1;
        ba2 = pd->gt_bond[i].elem2;
        if (((strcmp(ba1, elem1) == 0) && (strcmp(ba2, elem2) == 0)) ||
            ((strcmp(ba1, elem2) == 0) && (strcmp(ba2, elem1) == 0)))
        {
            dev = fabs((pd->gt_bond[i].length - distance)/pd->gt_bond[i].length);
            if (dev < toler)
            {
                for (k = 0; (k < nbo); k++)
                {
                    if (pd->gt_bond[i].bondorder == bo[k])
                    {
                        break;
                    }
                }
                if (k == nbo)
                {
                    srenew(bo, nbo+2);
                    bo[nbo]   = pd->gt_bond[i].bondorder;
                    bo[nbo+1] = 0;
                    nbo++;
                }
            }
        }
    }
    return bo;
}

int gmx_poldata_elem_is_bond(gmx_poldata_t pd, char *elem1, char *elem2,
                             double distance, double toler)
{
    char  *ba1, *ba2;
    double dev, dev_best = 100000;
    int    j, i;

    if ((NULL == elem1) || (NULL == elem2))
    {
        return 0;
    }
    for (i = 0; (i < pd->ngt_bond); i++)
    {
        if (0 == strlen(pd->gt_bond[i].elem1))
        {
            for (j = 0; (j < pd->nalexandria); j++)
            {
                if (strcmp(pd->alexandria[j].type, pd->gt_bond[i].atom2) == 0)
                {
                    strcpy(pd->gt_bond[i].elem1, pd->alexandria[j].elem);
                }
            }
        }
        if (0 == strlen(pd->gt_bond[i].elem2))
        {
            for (j = 0; (j < pd->nalexandria); j++)
            {
                if (strcmp(pd->alexandria[j].type, pd->gt_bond[i].atom2) == 0)
                {
                    strcpy(pd->gt_bond[i].elem2, pd->alexandria[j].elem);
                }
            }
        }
        ba1 = pd->gt_bond[i].elem1;
        ba2 = pd->gt_bond[i].elem2;
        if (((strcmp(ba1, elem1) == 0) && (strcmp(ba2, elem2) == 0)) ||
            ((strcmp(ba1, elem2) == 0) && (strcmp(ba2, elem1) == 0)))
        {
            dev = fabs((pd->gt_bond[i].length - distance)/pd->gt_bond[i].length);
            if (dev < dev_best)
            {
                dev_best = dev;
            }
        }
    }
    return (dev_best < toler);
}

static int lo_gtb_comp(t_gt_bond *ba, t_gt_bond *bb)
{
    char *a1, *a2, *b1, *b2;
    int   i;

    if (strcmp(ba->atom1, ba->atom2) <= 0)
    {
        a1 = ba->atom1;
        a2 = ba->atom2;
    }
    else
    {
        a2 = ba->atom1;
        a1 = ba->atom2;
    }
    if (strcmp(bb->atom1, bb->atom2) <= 0)
    {
        b1 = bb->atom1;
        b2 = bb->atom2;
    }
    else
    {
        b2 = bb->atom1;
        b1 = bb->atom2;
    }
    i = strcmp(a1, b1);
    if (0 == i)
    {
        i = strcmp(a2, b2);
    }

    return i;
}

static int gtb_comp(const void *a, const void *b)
{
    t_gt_bond *ba = (t_gt_bond *)a;
    t_gt_bond *bb = (t_gt_bond *)b;
    int        i;

    i = lo_gtb_comp(ba, bb);
    if ((0 == i) && ((0 != ba->bondorder) && (0 != bb->bondorder)))
    {
        if (ba->bondorder < bb->bondorder)
        {
            i = -1;
        }
        else if (ba->bondorder > bb->bondorder)
        {
            i = 1;
        }
    }
    return i;
}

static t_gt_bond *search_bond(gmx_poldata_t pd, char *atom1, char *atom2,
                              double bondorder)
{
    t_gt_bond key, *gt_b;
    int       i;

    key.atom1     = atom1;
    key.atom2     = atom2;
    key.bondorder = bondorder;

    gt_b = (t_gt_bond *) bsearch(&key, pd->gt_bond, pd->ngt_bond, sizeof(key), gtb_comp);
    if (NULL != gt_b)
    {
        i = gt_b - pd->gt_bond;
        while ((i > 0) && (lo_gtb_comp(&(pd->gt_bond[i-1]), &(pd->gt_bond[i])) == 0))
        {
            i--;
        }
        gt_b = &(pd->gt_bond[i]);
    }
    return gt_b;
}

double gmx_poldata_atype_bondorder(gmx_poldata_t pd, char *atype1, char *atype2,
                                   double distance, double toler)
{
    double     dev, dev_best = 100000;
    int        i, i_best = -1;
    t_gt_bond *gt_b;

    if ((NULL == atype1) || (NULL == atype2))
    {
        return 0.0;
    }
    gt_b = search_bond(pd, atype1, atype2, 0);
    if (NULL != gt_b)
    {
        i = gt_b - pd->gt_bond;
        do
        {
            dev = fabs(pd->gt_bond[i].length - distance);
            if (dev < dev_best)
            {
                dev_best = dev;
                i_best   = i;
            }
            i++;
        }
        while ((i < pd->ngt_bond) &&
               (0 == lo_gtb_comp(&(pd->gt_bond[i]), &(pd->gt_bond[i-1]))));
    }
    if (dev_best < toler)
    {
        return pd->gt_bond[i_best].bondorder;
    }

    return 0.0;
}

void gmx_poldata_add_miller(gmx_poldata_t pd,
                            char         *miller,
                            int           atomnumber,
                            double        tau_ahc,
                            double        alpha_ahp)
{
    t_miller *mil;

    pd->nmiller++;
    srenew(pd->miller, pd->nmiller);
    mil             = &(pd->miller[pd->nmiller-1]);
    mil->miller     = strdup(miller);
    mil->atomnumber = atomnumber;
    mil->tau_ahc    = tau_ahc;
    mil->alpha_ahp  = alpha_ahp;
}

void gmx_poldata_set_miller_units(gmx_poldata_t pd, char *tau_unit, char *ahp_unit)
{
    pd->miller_tau_unit = strdup(tau_unit);
    pd->miller_ahp_unit = strdup(ahp_unit);
}

void gmx_poldata_get_miller_units(gmx_poldata_t pd, char **tau_unit,
                                  char **ahp_unit)
{
    assign_str(tau_unit, pd->miller_tau_unit);
    assign_str(ahp_unit, pd->miller_ahp_unit);
}

int gmx_poldata_get_miller(gmx_poldata_t pd,
                           char        **miller,
                           int          *atomnumber,
                           double       *tau_ahc,
                           double       *alpha_ahp)
{
    t_miller *mil;
    int       i;

    i = pd->nmiller_c;

    if (i < pd->nmiller)
    {
        mil = &(pd->miller[i]);
        assign_str(miller, mil->miller);
        assign_scal(atomnumber, mil->atomnumber);
        assign_scal(tau_ahc, mil->tau_ahc);
        assign_scal(alpha_ahp, mil->alpha_ahp);
        pd->nmiller_c++;

        return 1;
    }
    else
    {
        pd->nmiller_c = 0;
    }

    return 0;
}

int gmx_poldata_get_miller_pol(gmx_poldata_t pd,
                               char         *miller,
                               int          *atomnumber,
                               double       *tau_ahc,
                               double       *alpha_ahp)
{
    t_miller *mil;
    int       i;

    for (i = 0; (i < pd->nmiller); i++)
    {
        if (strcmp(miller, pd->miller[i].miller) == 0)
        {
            mil = &(pd->miller[i]);
            assign_scal(atomnumber, mil->atomnumber);
            assign_scal(tau_ahc, mil->tau_ahc);
            assign_scal(alpha_ahp, mil->alpha_ahp);

            return 1;
        }
    }

    return 0;
}

void gmx_poldata_add_bosque(gmx_poldata_t pd,
                            char         *bosque,
                            double        polarizability)
{
    t_bosque *bs;

    pd->nbosque++;
    srenew(pd->bosque, pd->nbosque);
    bs                 = &(pd->bosque[pd->nbosque-1]);
    bs->bosque         = strdup(bosque);
    bs->polarizability = polarizability;
}

int gmx_poldata_get_bosque(gmx_poldata_t pd,
                           char        **bosque,
                           double       *polarizability)
{
    if (pd->nbosque_c < pd->nbosque)
    {
        assign_str(bosque, pd->bosque[pd->nbosque_c].bosque);
        assign_scal(polarizability, pd->bosque[pd->nbosque_c].polarizability);
        pd->nbosque_c++;

        return 1;
    }
    else
    {
        pd->nbosque_c = 0;
    }

    return 0;
}

int gmx_poldata_get_bosque_pol(gmx_poldata_t pd,
                               char         *bosque,
                               double       *polarizability)
{
    int i;

    for (i = 0; (i < pd->nbosque); i++)
    {
        if (strcasecmp(bosque, pd->bosque[i].bosque) == 0)
        {
            *polarizability = pd->bosque[i].polarizability;
            return 1;
        }
    }
    return 0;
}

void gmx_poldata_set_bosque_unit(gmx_poldata_t pd, char *polar_unit)
{
    pd->bosque_polar_unit   = strdup(polar_unit);
}

char *gmx_poldata_get_bosque_unit(gmx_poldata_t pd)
{
    return pd->bosque_polar_unit;
}

static int search_bondtype(gmx_poldata_t pd, char *atom)
{
    int j;

    for (j = 0; (j < pd->nbtype); j++)
    {
        if (strcmp(pd->btype[j], atom) == 0)
        {
            return j;
        }
    }
    return -1;
}

/*
 * gt_bond stuff
 */
int gmx_poldata_set_bond_params(gmx_poldata_t pd, char *atom1, char *atom2,
                                double length, double sigma, int ntrain,
                                double bondorder, char *params)
{
    t_gt_bond *gt_b;
    int        i;

    for (i = 0; (i < pd->ngt_bond); i++)
    {
        gt_b = &(pd->gt_bond[i]);
        if (((((strcmp(gt_b->atom1, atom1) == 0) &&
               (strcmp(gt_b->atom2, atom2) == 0)) ||
              ((strcmp(gt_b->atom1, atom2) == 0) &&
               (strcmp(gt_b->atom2, atom1) == 0)))) &&
            ((bondorder == 0) || (gt_b->bondorder == bondorder)))
        {
            break;
        }
    }
    if (i < pd->ngt_bond)
    {
        if (length > 0)
        {
            gt_b->length = length;
        }
        if (sigma > 0)
        {
            gt_b->sigma = sigma;
        }
        if (ntrain > 0)
        {
            gt_b->ntrain = ntrain;
        }
        if (NULL != gt_b->params)
        {
            sfree(gt_b->params);
        }
        gt_b->params = strdup(params);
        return 1;
    }
    return 0;
}

int gmx_poldata_add_bond(gmx_poldata_t pd, char *atom1, char *atom2,
                         double length, double sigma, int ntrain,
                         double bondorder, char *params)
{
    t_gt_bond *gt_b;
    int        a1, a2;

    if (-1 == (a1 = search_bondtype(pd, atom1)))
    {
        return 0;
    }
    if (-1 == (a2 = search_bondtype(pd, atom2)))
    {
        return 0;
    }
    if (gmx_poldata_set_bond_params(pd, atom1, atom2, length, sigma, ntrain,
                                    bondorder, params) == 0)
    {
        pd->ngt_bond++;
        srenew(pd->gt_bond, pd->ngt_bond);
        gt_b            = &(pd->gt_bond[pd->ngt_bond-1]);
        gt_b->atom1     = strdup(atom1);
        strncpy(gt_b->elem1, pd->alexandria[a1].elem, sizeof(gt_b->elem1)-1);
        gt_b->atom2     = strdup(atom2);
        strncpy(gt_b->elem2, pd->alexandria[a2].elem, sizeof(gt_b->elem2)-1);
        gt_b->length    = length;
        gt_b->sigma     = sigma;
        gt_b->ntrain    = ntrain;
        gt_b->bondorder = bondorder;
        gt_b->params    = strdup(params);
        qsort(pd->gt_bond, pd->ngt_bond, sizeof(pd->gt_bond[0]), gtb_comp);
    }
    return 1;
}

int gmx_poldata_get_bond(gmx_poldata_t pd, char **atom1, char **atom2,
                         double *length, double *sigma, int *ntrain,
                         double *bondorder, char **params)
{
    t_gt_bond *gt_b;

    if (pd->ngt_bond_c < pd->ngt_bond)
    {
        gt_b = &(pd->gt_bond[pd->ngt_bond_c]);
        assign_str(atom1, gt_b->atom1);
        assign_str(atom2, gt_b->atom2);
        assign_scal(length, gt_b->length);
        assign_scal(sigma, gt_b->sigma);
        assign_scal(ntrain, gt_b->ntrain);
        assign_scal(bondorder, gt_b->bondorder);
        assign_str(params, gt_b->params);
        pd->ngt_bond_c++;

        return pd->ngt_bond_c;
    }
    pd->ngt_bond_c = 0;

    return 0;
}

int gmx_poldata_search_bond(gmx_poldata_t pd, char *atom1, char *atom2,
                            double *length, double *sigma, int *ntrain,
                            double *bondorder, char **params)
{
    t_gt_bond *gt_b;

    if ((NULL == atom1) || (NULL == atom2))
    {
        return 0;
    }
    gt_b = search_bond(pd, atom1, atom2, 0);
    if (NULL != gt_b)
    {
        if (((strcmp(gt_b->atom1, atom1) == 0) &&
             (strcmp(gt_b->atom2, atom2) == 0)) ||
            ((strcmp(gt_b->atom1, atom2) == 0) &&
             (strcmp(gt_b->atom2, atom1) == 0)))
        {
            assign_scal(length, gt_b->length);
            assign_scal(sigma, gt_b->sigma);
            assign_scal(ntrain, gt_b->ntrain);
            assign_scal(bondorder, gt_b->bondorder);
            assign_str(params, gt_b->params);

            return 1+(gt_b - pd->gt_bond);
        }
    }
    return 0;
}

/*
 * gt_angle stuff
 */
int gmx_poldata_set_angle_params(gmx_poldata_t pd, char *atom1, char *atom2,
                                 char *atom3, double angle, double sigma, int ntrain,
                                 char *params)
{
    t_gt_angle *gt_b;
    int         i;

    for (i = 0; (i < pd->ngt_angle); i++)
    {
        gt_b = &(pd->gt_angle[i]);
        if ((strcmp(gt_b->atom2, atom2) == 0) &&
            (((strcmp(gt_b->atom1, atom1) == 0) &&
              (strcmp(gt_b->atom3, atom3) == 0)) ||
             ((strcmp(gt_b->atom1, atom3) == 0) &&
              (strcmp(gt_b->atom3, atom1) == 0))))
        {
            break;
        }
    }
    if (i < pd->ngt_angle)
    {
        if (angle > 0)
        {
            gt_b->angle = angle;
        }
        if (sigma > 0)
        {
            gt_b->sigma = sigma;
        }
        if (ntrain > 0)
        {
            gt_b->ntrain = ntrain;
        }
        if (NULL != gt_b->params)
        {
            sfree(gt_b->params);
        }
        gt_b->params = strdup(params);
        return 1;
    }
    return 0;
}

int gmx_poldata_add_angle(gmx_poldata_t pd, char *atom1, char *atom2,
                          char *atom3, double angle, double sigma,
                          int ntrain, char *params)
{
    t_gt_angle *gt_b;

    if ((-1 == search_bondtype(pd, atom1)) ||
        (-1 == search_bondtype(pd, atom2)) ||
        (-1 == search_bondtype(pd, atom3)))
    {
        return 0;
    }

    if (0 == gmx_poldata_set_angle_params(pd, atom1, atom2, atom3, angle, sigma, ntrain, params))
    {
        pd->ngt_angle++;
        srenew(pd->gt_angle, pd->ngt_angle);
        gt_b          = &(pd->gt_angle[pd->ngt_angle-1]);
        gt_b->atom1   = strdup(atom1);
        gt_b->atom2   = strdup(atom2);
        gt_b->atom3   = strdup(atom3);
        gt_b->angle   = angle;
        gt_b->sigma   = sigma;
        gt_b->ntrain  = ntrain;
        gt_b->params  = strdup(params);
    }
    return 1;
}

int gmx_poldata_get_angle(gmx_poldata_t pd, char **atom1, char **atom2,
                          char **atom3, double *angle, double *sigma,
                          int *ntrain, char **params)
{
    t_gt_angle *gt_b;

    if (pd->ngt_angle_c < pd->ngt_angle)
    {
        gt_b = &(pd->gt_angle[pd->ngt_angle_c]);
        assign_str(atom1, gt_b->atom1);
        assign_str(atom2, gt_b->atom2);
        assign_str(atom3, gt_b->atom3);
        assign_scal(angle, gt_b->angle);
        assign_scal(sigma, gt_b->sigma);
        assign_scal(ntrain, gt_b->ntrain);
        assign_str(params, gt_b->params);
        pd->ngt_angle_c++;

        return pd->ngt_angle_c;
    }
    pd->ngt_angle_c = 0;

    return 0;
}

int gmx_poldata_search_angle(gmx_poldata_t pd, char *atom1, char *atom2,
                             char *atom3, double *angle, double *sigma,
                             int *ntrain, char **params)
{
    t_gt_angle *gt_b;
    int         i;

    if ((NULL == atom1) || (NULL == atom2) || (NULL == atom3))
    {
        return 0;
    }
    for (i = 0; (i < pd->ngt_angle); i++)
    {
        gt_b = &(pd->gt_angle[i]);
        if ((strcmp(gt_b->atom2, atom2) == 0) &&
            (((strcmp(gt_b->atom1, atom1) == 0) &&
              (strcmp(gt_b->atom3, atom3) == 0)) ||
             ((strcmp(gt_b->atom1, atom3) == 0) &&
              (strcmp(gt_b->atom3, atom1) == 0))))
        {
            assign_scal(angle, gt_b->angle);
            assign_scal(sigma, gt_b->sigma);
            assign_scal(ntrain, gt_b->ntrain);
            assign_str(params, gt_b->params);

            return i+1;
        }
    }
    return 0;
}

void gmx_poldata_set_angle_unit(gmx_poldata_t pd, char *angle_unit)
{
    pd->gt_angle_unit   = strdup(angle_unit);
}

char *gmx_poldata_get_angle_unit(gmx_poldata_t pd)
{
    return pd->gt_angle_unit;
}

/*
 * gt_dihedral stuff
 */
static int gtd_comp(const void *a, const void *b)
{
    t_gt_dihedral *gt_a = (t_gt_dihedral *)a;
    t_gt_dihedral *gt_b = (t_gt_dihedral *)b;
    int            n;

    if (0 == (n = strcmp(gt_a->atom1, gt_b->atom1)))
    {
        if (0 == (n = strcmp(gt_a->atom2, gt_b->atom2)))
        {
            if (0 == (n = strcmp(gt_a->atom3, gt_b->atom3)))
            {
                n = strcmp(gt_a->atom4, gt_b->atom4);
            }
        }
    }

    return n;
}

static t_gt_dihedral *search_dihedral(gmx_poldata_t pd, int egd,
                                      char *atom1, char *atom2,
                                      char *atom3, char *atom4)
{
    t_gt_dihedral gt_a, *gt_res, *gt_dptr;
    int           nd;

    if ((NULL == atom1) || (NULL == atom2) || (NULL == atom3) || (NULL == atom4))
    {
        return NULL;
    }
    gt_dptr    = pd->gt_dihedral[egd];
    nd         = pd->ngt_dihedral[egd];
    gt_a.atom1 = atom1;
    gt_a.atom2 = atom2;
    gt_a.atom3 = atom3;
    gt_a.atom4 = atom4;
    gt_res     = (t_gt_dihedral *) bsearch(&gt_a, gt_dptr, nd, sizeof(gt_a), &gtd_comp);
    if (NULL == gt_res)
    {
        gt_a.atom1 = atom4;
        gt_a.atom2 = atom3;
        gt_a.atom3 = atom2;
        gt_a.atom4 = atom1;
        gt_res     = (t_gt_dihedral *) bsearch(&gt_a, gt_dptr, nd, sizeof(gt_a), gtd_comp);
    }
    return gt_res;
}

int gmx_poldata_set_dihedral_params(gmx_poldata_t pd, int egd,
                                    char *atom1, char *atom2,
                                    char *atom3, char *atom4,
                                    double dihedral, double sigma, int ntrain,
                                    char *params)
{
    t_gt_dihedral *gt_b;


    gt_b = search_dihedral(pd, egd, atom1, atom2, atom3, atom4);
    if (NULL != gt_b)
    {
        if (dihedral > 0)
        {
            gt_b->dihedral = dihedral;
        }
        if (sigma > 0)
        {
            gt_b->sigma = sigma;
        }
        if (ntrain > 0)
        {
            gt_b->ntrain = ntrain;
        }
        if (NULL != gt_b->params)
        {
            sfree(gt_b->params);
        }
        gt_b->params = strdup(params);
        return 1;
    }
    return 0;
}

int gmx_poldata_add_dihedral(gmx_poldata_t pd, int egd,
                             char *atom1, char *atom2,
                             char *atom3, char *atom4, double dihedral,
                             double sigma, int ntrain, char *params)
{
    t_gt_dihedral *gt_b;

    if ((-1 == search_bondtype(pd, atom1)) ||
        (-1 == search_bondtype(pd, atom2)) ||
        (-1 == search_bondtype(pd, atom3)) ||
        (-1 == search_bondtype(pd, atom4)))
    {
        return 0;
    }

    if (0 == gmx_poldata_set_dihedral_params(pd, egd, atom1, atom2,
                                             atom3, atom4, dihedral,
                                             sigma, ntrain, params))
    {
        pd->ngt_dihedral[egd]++;
        srenew(pd->gt_dihedral[egd], pd->ngt_dihedral[egd]);
        gt_b           = &(pd->gt_dihedral[egd][pd->ngt_dihedral[egd]-1]);
        gt_b->atom1    = strdup(atom1);
        gt_b->atom2    = strdup(atom2);
        gt_b->atom3    = strdup(atom3);
        gt_b->atom4    = strdup(atom4);
        gt_b->dihedral = dihedral;
        gt_b->sigma    = sigma;
        gt_b->ntrain   = ntrain;
        gt_b->params   = strdup(params);
        qsort(pd->gt_dihedral[egd], pd->ngt_dihedral[egd], sizeof(pd->gt_dihedral[egd][0]),
              gtd_comp);
    }
    return 1;
}

int gmx_poldata_get_dihedral(gmx_poldata_t pd, int egd,
                             char **atom1, char **atom2,
                             char **atom3, char **atom4, double *dihedral,
                             double *sigma, int *ntrain, char **params)
{
    t_gt_dihedral *gt_b;

    if (pd->ngt_dihedral_c[egd] < pd->ngt_dihedral[egd])
    {
        gt_b = &(pd->gt_dihedral[egd][pd->ngt_dihedral_c[egd]]);
        assign_str(atom1, gt_b->atom1);
        assign_str(atom2, gt_b->atom2);
        assign_str(atom3, gt_b->atom3);
        assign_str(atom4, gt_b->atom4);
        assign_scal(dihedral, gt_b->dihedral);
        assign_scal(sigma, gt_b->sigma);
        assign_scal(ntrain, gt_b->ntrain);
        assign_str(params, gt_b->params);
        pd->ngt_dihedral_c[egd]++;

        return pd->ngt_dihedral_c[egd];
    }
    pd->ngt_dihedral_c[egd] = 0;

    return 0;
}

int gmx_poldata_search_dihedral(gmx_poldata_t pd, int egd,
                                char *atom1, char *atom2,
                                char *atom3, char *atom4,
                                double *dihedral, double *sigma,
                                int *ntrain, char **params)
{
    t_gt_dihedral *gt_res;

    gt_res = search_dihedral(pd, egd, atom1, atom2, atom3, atom4);
    if (NULL != gt_res)
    {
        assign_scal(dihedral, gt_res->dihedral);
        assign_scal(sigma, gt_res->sigma);
        assign_scal(ntrain, gt_res->ntrain);
        assign_str(params, gt_res->params);

        return 1 + (int) (gt_res - pd->gt_dihedral[egd]);
    }
    return 0;
}

void gmx_poldata_add_symcharges(gmx_poldata_t pd, char *central,
                                char *attached, int numattach)
{
    t_symcharges *sc;
    int           i;

    for (i = 0; (i < pd->nsymcharges); i++)
    {
        sc = &(pd->symcharges[i]);
        if ((strcasecmp(sc->central, central) == 0) &&
            (strcasecmp(sc->attached, attached) == 0) &&
            (sc->numattach == numattach))
        {
            break;
        }
    }
    if (i == pd->nsymcharges)
    {
        pd->nsymcharges++;
        srenew(pd->symcharges, pd->nsymcharges);
        sc              = &(pd->symcharges[i]);
        sc->central     = strdup(central);
        sc->attached    = strdup(attached);
        sc->numattach   = numattach;
    }
}

int gmx_poldata_get_symcharges(gmx_poldata_t pd, char **central,
                               char **attached, int *numattach)
{
    t_symcharges *sc;

    if (pd->nsymcharges_c < pd->nsymcharges)
    {
        sc = &(pd->symcharges[pd->nsymcharges_c]);
        assign_str(central, sc->central);
        assign_str(attached, sc->attached);
        assign_scal(numattach, sc->numattach);
        pd->nsymcharges_c++;

        return 1;
    }
    pd->nsymcharges_c = 0;

    return 0;
}

int gmx_poldata_search_symcharges(gmx_poldata_t pd, char *central,
                                  char *attached, int numattach)
{
    t_symcharges *sc;
    int           i;

    for (i = 0; (i < pd->nsymcharges); i++)
    {
        sc = &(pd->symcharges[i]);
        if ((strcasecmp(sc->central, central) == 0) &&
            (strcasecmp(sc->attached, attached) == 0) &&
            (sc->numattach == numattach))
        {
            return 1;
        }
    }

    return 0;
}

/* Electrostatics properties */
static t_eemprops *get_eep(gmx_poldata *pd, ChargeGenerationModel eqg_model,
                           const char *name)
{
    int i;

    for (i = 0; (i < pd->nep); i++)
    {
        if ((strcasecmp(pd->eep[i].name, name) == 0) &&
            (pd->eep[i].eqg_model == eqg_model))
        {
            return &(pd->eep[i]);
        }
    }
    return NULL;
}

void gmx_poldata_set_eemprops(gmx_poldata_t pd,
                              ChargeGenerationModel eqg_model, char *name,
                              double J0, double chi0, char *zeta, char *q, char *row)
{
    gmx_poldata             *pold = (gmx_poldata *) pd;
    t_eemprops              *eep;
    std::vector<std::string> sz, sq, sr;

    eep = get_eep(pold, eqg_model, name);
    if (NULL == eep)
    {
        srenew(pd->eep, ++pd->nep);
        eep = &(pd->eep[pd->nep-1]);
    }
    eep->eqg_model = eqg_model;
    strncpy(eep->name, name, EEMBUFSIZE-1);
    eep->name[EEMBUFSIZE-1] = '\0';
    eep->J0                 = J0;
    sz = split(zeta, ' ');
    sq = split(q, ' ');
    sr = split(row, ' ');
    strncpy(eep->zetastr, zeta, EEMBUFSIZE-1);
    strncpy(eep->qstr, q, EEMBUFSIZE-1);
    strncpy(eep->rowstr, row, EEMBUFSIZE-1);
    unsigned int nn = std::min(sz.size(), std::min(sq.size(), sr.size()));
    unsigned int n;
    for (n = 0; (n < nn); n++)
    {
        if (n < MAXZETA)
        {
            eep->zeta[n] = atof(sz[n].c_str());
            eep->q[n]    = atof(sq[n].c_str());
            eep->row[n]  = atoi(sr[n].c_str());
        }
    }
    if (sz.size() > nn)
    {
        fprintf(stderr, "Warning: more zeta values than q/row values for %s n = %d\n",
                name, nn);
    }
    if (sq.size() > nn)
    {
        fprintf(stderr, "Warning: more q values than zeta/row values for %s n = %d\n",
                name, nn);
    }
    if (sr.size() > nn)
    {
        fprintf(stderr, "Warning: more row values than q/zeta values for %s n = %d\n",
                name, nn);
    }
    eep->nzeta = nn;
    if (nn >= MAXZETA)
    {
        fprintf(stderr, "More than %d zeta and/or q values for %s\n", MAXZETA, eep->name);
        eep->nzeta = MAXZETA;
    }
    for (; (n < MAXZETA); n++)
    {
        eep->zeta[n] = 0;
        eep->q[n]    = 0;
        eep->row[n]  = 0;
    }
    eep->chi0  = chi0;
}

int gmx_poldata_get_eemprops(gmx_poldata_t pd,
                             ChargeGenerationModel *eqg_model, char **name,
                             double *J0, double *chi0, char **zeta, char **q, char **row)
{
    if (pd->nep_c < pd->nep)
    {
        assign_scal(eqg_model, pd->eep[pd->nep_c].eqg_model);
        assign_str(name, pd->eep[pd->nep_c].name);
        assign_scal(J0, pd->eep[pd->nep_c].J0);
        assign_str(zeta, pd->eep[pd->nep_c].zetastr);
        assign_str(q, pd->eep[pd->nep_c].qstr);
        assign_str(row, pd->eep[pd->nep_c].rowstr);
        assign_scal(chi0, pd->eep[pd->nep_c].chi0);
        pd->nep_c++;
        return 1;
    }
    else
    {
        pd->nep_c = 0;
        return 0;
    }
}

int gmx_poldata_get_numprops(gmx_poldata_t pd, ChargeGenerationModel eqg_model)
{
    int i, n = 0;

    for (i = 0; (i < pd->nep); i++)
    {
        if (pd->eep[i].eqg_model == eqg_model)
        {
            n++;
        }
    }

    return n;
}

int gmx_poldata_have_pol_support(gmx_poldata_t pd, const char *atype)
{
    int i;

    for (i = 0; (i < pd->nalexandria); i++)
    {
        if (strcmp(atype, pd->alexandria[i].type) == 0)
        {
            return 1;
        }
    }
    return 0;
}

int gmx_poldata_have_eem_support(gmx_poldata_t pd, ChargeGenerationModel eqg_model,
                                 const char *name,
                                 gmx_bool bAllowZeroParameters)
{
    gmx_poldata *pold = (gmx_poldata *) pd;
    t_eemprops  *eep  = get_eep(pold, eqg_model, name);

    return (eep && (bAllowZeroParameters || ((eep->J0 > 0) && (eep->chi0 > 0))));
}

double gmx_poldata_get_j00(gmx_poldata_t pd, ChargeGenerationModel eqg_model, char *name)
{
    gmx_poldata *pold = (gmx_poldata *) pd;
    t_eemprops  *eer;

    if ((eer = get_eep(pold, eqg_model, name)) != NULL)
    {
        return eer->J0;
    }
    else
    {
        gmx_fatal(FARGS, "No J0 data for eqg_model %d and name %s",
                  eqg_model, name);
    }
    return -1;
}

char *gmx_poldata_get_qstr(gmx_poldata_t pd, ChargeGenerationModel eqg_model, char *name)
{
    t_eemprops *eer;

    if ((eer = get_eep(pd, eqg_model, name)) != NULL)
    {
        return eer->qstr;
    }
    return NULL;
}

char *gmx_poldata_get_rowstr(gmx_poldata_t pd, ChargeGenerationModel eqg_model, char *name)
{
    t_eemprops *eer;

    if ((eer = get_eep(pd, eqg_model, name)) != NULL)
    {
        return eer->rowstr;
    }
    return NULL;
}

int gmx_poldata_get_row(gmx_poldata_t pd, ChargeGenerationModel eqg_model, char *name, int zz)
{
    t_eemprops *eer;

    if ((eer = get_eep(pd, eqg_model, name)) != NULL)
    {
        range_check(zz, 0, eer->nzeta);
        return eer->row[zz];
    }
    return -1;
}

double gmx_poldata_get_zeta(gmx_poldata_t pd, ChargeGenerationModel eqg_model, char *name, int zz)
{
    t_eemprops *eer;

    if ((eer = get_eep(pd, eqg_model, name)) != NULL)
    {
        if ((zz < 0) || (zz >= eer->nzeta))
        {
            printf("Bleh\n");
        }
        range_check(zz, 0, eer->nzeta);
        return eer->zeta[zz];
    }
    return -1;
}

int gmx_poldata_get_nzeta(gmx_poldata_t pd, ChargeGenerationModel eqg_model, char *name)
{
    t_eemprops *eer;

    if ((eer = get_eep(pd, eqg_model, name)) != NULL)
    {
        return eer->nzeta;
    }
    return 0;
}

double gmx_poldata_get_q(gmx_poldata_t pd, ChargeGenerationModel eqg_model, char *name, int zz)
{
    t_eemprops *eer;

    if ((eer = get_eep(pd, eqg_model, name)) != NULL)
    {
        range_check(zz, 0, eer->nzeta);
        return eer->q[zz];
    }
    return -1;
}

double gmx_poldata_get_chi0(gmx_poldata_t pd, ChargeGenerationModel eqg_model, char *name)
{
    t_eemprops *eer;

    if ((eer = get_eep(pd, eqg_model, name)) != NULL)
    {
        return eer->chi0;
    }
    else
    {
        gmx_fatal(FARGS, "No chi0 data for eqg_model %d and name %s", eqg_model, name);
    }
    return -1;
}

void gmx_poldata_set_epref(gmx_poldata_t pd, ChargeGenerationModel eqg_model, char *epref)
{
    int i;

    for (i = 0; (i < pd->ner); i++)
    {
        if (pd->epr[i].eqg_model == eqg_model)
        {
            if (pd->epr[i].epref)
            {
                sfree(pd->epr[i].epref);
            }
            pd->epr[i].epref = strdup(epref);
            break;
        }
    }
    if (i == pd->ner)
    {
        srenew(pd->epr, ++pd->ner);
        pd->epr[i].eqg_model = eqg_model;
        pd->epr[i].epref     = strdup(epref);
    }
}

char *gmx_poldata_get_epref(gmx_poldata_t pd, ChargeGenerationModel eqg_model)
{
    int i;

    for (i = 0; (i < pd->ner); i++)
    {
        if (pd->epr[i].eqg_model == eqg_model)
        {
            return pd->epr[i].epref;
        }
    }
    return NULL;
}

int gmx_poldata_list_epref(gmx_poldata_t pd, ChargeGenerationModel *eqg_model, char **epref)
{
    if (pd->ner_c < pd->ner)
    {
        assign_scal(eqg_model, pd->epr[pd->ner_c].eqg_model);
        assign_str(epref, pd->epr[pd->ner_c].epref);
        pd->ner_c++;
        return 1;
    }
    pd->ner_c = 0;

    return 0;
}

void gmx_poldata_comm_eemprops(gmx_poldata_t pd, t_commrec *cr)
{
    int         i, j, nep;
    t_eemprops *ep;

    if (NULL != debug)
    {
        fprintf(debug, "Going to update eemprops on node %d\n", cr->nodeid);
    }
    if (MASTER(cr))
    {
        for (i = 1; (i < cr->nnodes); i++)
        {
            gmx_send_int(cr, i, pd->nep);
            gmx_send(cr, i, pd->eep, pd->nep*sizeof(pd->eep[0]));
        }
    }
    else
    {
        nep = gmx_recv_int(cr, 0);
        if (nep != pd->nep)
        {
            gmx_fatal(FARGS, "Inconsistency in number of EEM parameters");
        }
        snew(ep, pd->nep);
        gmx_recv(cr, 0, ep, pd->nep*sizeof(ep[0]));
        for (i = 0; (i < pd->nep); i++)
        {
            pd->eep[i] = ep[i];
        }
        sfree(ep);
    }
    if (NULL != debug)
    {
        fprintf(debug, "  EEP  Atom      Chi      J00     Zeta\n");
        for (i = 0; (i < pd->nep); i++)
        {
            fprintf(debug, "%5s %5s %8.3f %8.3f",
                    get_eemtype_name(pd->eep[i].eqg_model),
                    pd->eep[i].name, pd->eep[i].chi0, pd->eep[i].J0);
            for (j = 0; (j < pd->eep[i].nzeta); j++)
            {
                fprintf(debug, " %8.3f", pd->eep[i].zeta[j]);
            }
            fprintf(debug, "\n");
        }
    }
}

void gmx_poldata_comm_force_parameters(gmx_poldata_t pd, t_commrec *cr)
{
    int         i, j, nep;
    t_eemprops *ep;

    if (NULL != debug)
    {
        fprintf(debug, "Going to update force parameters on node %d\n", cr->nodeid);
    }
    if (MASTER(cr))
    {
        for (i = 1; (i < cr->nnodes); i++)
        {
            gmx_send_int(cr, i, pd->nep);
            gmx_send(cr, i, pd->eep, pd->nep*sizeof(pd->eep[0]));
        }
    }
    else
    {
        nep = gmx_recv_int(cr, 0);
        if (nep != pd->nep)
        {
            gmx_fatal(FARGS, "Inconsistency in number of EEM parameters");
        }
        snew(ep, pd->nep);
        gmx_recv(cr, 0, ep, pd->nep*sizeof(ep[0]));
        for (i = 0; (i < pd->nep); i++)
        {
            pd->eep[i] = ep[i];
        }
        sfree(ep);
    }
    if (NULL != debug)
    {
        fprintf(debug, "  EEP  Atom      Chi      J00     Zeta\n");
        for (i = 0; (i < pd->nep); i++)
        {
            fprintf(debug, "%5s %5s %8.3f %8.3f",
                    get_eemtype_name(pd->eep[i].eqg_model),
                    pd->eep[i].name, pd->eep[i].chi0, pd->eep[i].J0);
            for (j = 0; (j < pd->eep[i].nzeta); j++)
            {
                fprintf(debug, " %8.3f", pd->eep[i].zeta[j]);
            }
            fprintf(debug, "\n");
        }
    }
}

typedef struct {
    ChargeGenerationModel eqg;
    const char           *name, *ref;
    gmx_bool              bWeight;
} t_eemtype_props;

static t_eemtype_props eemtype_props[eqgNR] = {
    { eqgNone,     "None",     "None",          FALSE },
    { eqgAXp,      "AXp",      "Maaren2014a",   FALSE },
    { eqgAXg,      "AXg",      "Maaren2014a",   TRUE },
    { eqgAXs,      "AXs",      "Maaren2014a",   TRUE },
    { eqgESP,      "ESP",      "Kollman1991a",  FALSE },
    { eqgRESP,     "RESP",     "Kollman1991a",  FALSE },
    { eqgRESPG,    "RESPG",    "Maaren2014a",   FALSE },
    { eqgYang,     "Yang",     "Yang2006b",     TRUE },
    { eqgBultinck, "Bultinck", "Bultinck2002a", FALSE },
    { eqgRappe,    "Rappe",    "Rappe1991a",    TRUE }
};

ChargeGenerationModel name2eemtype(const char *name)
{
    int i;

    for (i = eqgNone; (i < eqgNR); i++)
    {
        if (strcasecmp(name, eemtype_props[i].name) == 0)
        {
            return eemtype_props[i].eqg;
        }
    }
    return eqgNR;
}

const char *get_eemtype_name(ChargeGenerationModel eem)
{
    int i;

    for (i = eqgNone; (i < eqgNR); i++)
    {
        if (eem == eemtype_props[i].eqg)
        {
            return eemtype_props[i].name;
        }
    }

    return NULL;
}
