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
 * MERCHANTABILITY or FITN
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


int gtb_comp(const void *a, const void *b);
int gtd_comp(const void *a, const void *b);

#define assign_str(dst, src)  if (dst) { if (src) {*dst = strdup(src); }else{*dst = NULL; }}
#define assign_scal(dst, src) if (dst) *dst = src

namespace alexandria
{

Poldata::Poldata()
{
   
    for (int x =0; x < egdNR; x++){
      gt_dihedral_function[x]=NULL;
      ngt_dihedral[x]=0;
      ngt_dihedral_c[x]=0;
      gt_dihedral[x]=NULL;
    }

    /* Initiate some crucial variables */
    this->nexcl          = NOTSET;
    this->fudgeQQ        = NOTSET;
    this->fudgeLJ        = NOTSET;
    this->gt_comb_rule   = NOTSET;
    this->gt_vdw_ftype   = NOTSET;
    this->gt_bond_ftype  = NOTSET;
    this->gt_angle_ftype = NOTSET;
    for (int i = 0; (i < egdNR); i++)
    {
        this->gt_dihedral_ftype[i] = NOTSET;
    }

}

void Poldata::set_filename( char *fn2)
{
    if (NULL == fn2)
    {
        fprintf(stderr, "Trying to set Poldata filename to NULL\n");
        return;
    }
    if (NULL != this->filename)
    {
        fprintf(stderr, "Overwriting Poldata filename from %s to %s\n",
                this->filename, fn2);
        sfree(this->filename);
    }
    this->filename = strdup(fn2);
}

void Poldata::set_vdw_function( const char *func)
{
    int i;

    if (NULL != this->gt_vdw_function)
    {
        sfree(this->gt_vdw_function);
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
    this->gt_vdw_ftype    = i;
    this->gt_vdw_function = strdup(func);
}

char *Poldata::get_vdw_function()
{
    return this->gt_vdw_function;
}

int Poldata::get_vdw_ftype()
{
    if (NOTSET == this->gt_vdw_ftype)
    {
        gmx_fatal(FARGS, "No valid function type for Van der Waals function in %s",
                  this->filename);
    }
    return this->gt_vdw_ftype;
}

void Poldata::set_nexcl( int nexcl)
{
    this->nexcl = nexcl;
}

int Poldata::get_nexcl()
{
    if (NOTSET == this->nexcl)
    {
        gmx_fatal(FARGS, "Nexclusions not set in %s", this->filename);
    }
    return this->nexcl;
}

void Poldata::set_fudgeQQ( double fudgeQQ)
{
    this->fudgeQQ = fudgeQQ;
}

double Poldata::get_fudgeQQ()
{
    if (NOTSET == this->fudgeQQ)
    {
        gmx_fatal(FARGS, "fudgeQQ not set in %s", this->filename);
    }
    return this->fudgeQQ;
}

void Poldata::set_fudgeLJ( double fudgeLJ)
{
    this->fudgeLJ = fudgeLJ;
}

double Poldata::get_fudgeLJ()
{
    if (NOTSET == this->fudgeLJ)
    {
        gmx_fatal(FARGS, "fudgeLJ not set in %s", this->filename);
    }
    return this->fudgeLJ;
}

void Poldata::set_combination_rule( char *func)
{
    int i;

    if (NULL != this->gt_combination_rule)
    {
        sfree(this->gt_combination_rule);
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
    this->gt_comb_rule        = i;
    this->gt_combination_rule = strdup(func);
}

char *Poldata::get_combination_rule()
{
    if (NOTSET == this->gt_vdw_ftype)
    {
        gmx_fatal(FARGS, "No valid function type for Van der Waals function in %s",
                  this->filename);
    }
    return this->gt_combination_rule;
}

int Poldata::get_comb_rule()
{
    if (NOTSET == this->gt_comb_rule)
    {
        gmx_fatal(FARGS, "No valid combination rule in %s", this->filename);
    }
    return this->gt_comb_rule;
}

void Poldata::add_ptype(
                           const char   *ptype,
                           const char   *miller,
                           const char   *bosque,
                           double        polarizability,
                           double        sig_pol)
{
    t_ptype *sp;
    int      i;

    for (i = 0; (i < this->nptype); i++)
    {
        if (strcmp(this->ptype[i].type, ptype) == 0)
        {
            break;
        }
    }
    if (i == this->nptype)
    {
        this->nptype++;
        srenew(this->ptype, this->nptype);

        sp                 = &(this->ptype[i]);
        sp->type           = strdup(ptype);
        sp->bosque         = strdup(bosque);
        sp->miller         = strdup(miller);
        sp->polarizability = polarizability;
        sp->sig_pol        = sig_pol;
    }
    else
    {
        fprintf(stderr, "Polarizability type %s was already added to Poldata record\n", ptype);
    }
}

 void Poldata::add_btype(
                                  const char   *btype)
{
    int i;

    for (i = 0; (i < this->nbtype); i++)
    {
        if (strcmp(btype, this->btype[i]) == 0)
        {
            break;
        }
    }
    if (i == this->nbtype)
    {
        srenew(this->btype, ++this->nbtype);
        this->btype[i] = strdup(btype);
    }
}

void Poldata::add_atype(
                           const char   *elem,
                           const char   *desc,
                           const char   *atype,
                           const char   *ptype,
                           const char   *btype,
                           const char   *vdwparams,
                           double        ref_enthalpy)
{
    t_ffatype *sp;
    int        i;

    for (i = 0; (i < this->nalexandria); i++)
    {
        if (strcmp(this->alexandria[i].type, atype) == 0)
        {
            break;
        }
    }
    if (i == this->nalexandria)
    {
        this->nalexandria++;
        srenew(this->alexandria, this->nalexandria);

        sp                 = &(this->alexandria[i]);
        sp->elem           = strdup(elem);
        sp->desc           = strdup(desc);
        sp->type           = strdup(atype);
        sp->ptype          = strdup(ptype);
        sp->btype          = strdup(btype);
        sp->vdwparams      = strdup(vdwparams);
        sp->ref_enthalpy   = ref_enthalpy;
        add_btype(btype);
    }
    else
    {
        fprintf(stderr, "Atom type %s was already added to Poldata record\n", atype);
    }
}

void Poldata::add_bonding_rule(
                                  char *gt_brule, char *atype,
                                  char *geometry, int numbonds,
                                  double valence, int iAromatic,
                                  char *neighbors)
{
    t_brule *sp;
    int      i, j;

    for (j = 0; (j < this->nalexandria); j++)
    {
        if (strcmp(this->alexandria[j].type, atype) == 0)
        {
            break;
        }
    }
    if (j < this->nalexandria)
    {
        for (i = 0; (i < this->nbrule); i++)
        {
            if (strcmp(this->brule[i].rule, gt_brule) == 0)
            {
                break;
            }
        }
        if (i == this->nbrule)
        {
            this->nbrule++;
            srenew(this->brule, this->nbrule);

            sp                 = &(this->brule[i]);
            sp->elem           = strdup(this->alexandria[j].elem);
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
            fprintf(stderr, "Bonding rule %s was already added to Poldata record\n", gt_brule);
        }
    }
    else if (NULL != debug)
    {
        fprintf(debug, "Ignoring bonding rule involving unknown atom type %s\n",
                atype);
    }
}

int Poldata::get_bonding_rule(
                                 char **gt_brule, char **atype,
                                 char **geometry, int *numbonds,
                                 double *valence, int *iAromatic,
                                 char **neighbors)
{
    if (this->nbrule_c < this->nbrule)
    {
        assign_str(gt_brule, this->brule[this->nbrule_c].rule);
        assign_str(atype, this->brule[this->nbrule_c].type);
        assign_str(geometry, this->brule[this->nbrule_c].geometry);
        assign_scal(numbonds, this->brule[this->nbrule_c].numbonds);
        assign_scal(valence, this->brule[this->nbrule_c].valence);
        assign_scal(iAromatic, this->brule[this->nbrule_c].iAromatic);
        assign_str(neighbors, this->brule[this->nbrule_c].neighbors);

        this->nbrule_c++;

        return 1;
    }
    return 0;
}

int Poldata::get_natypes()
{
    return this->nalexandria;
}

int Poldata::get_nptypes()
{
    return this->nptype;
}

int Poldata::get_ngt_bond()
{
    return this->ngt_bond;
}

int Poldata::get_ngt_angle()
{
    return this->ngt_angle;
}

int Poldata::get_ngt_dihedral( int egd)
{
    return this->ngt_dihedral[egd];
}

void Poldata::set_bond_function( char *fn)
{
    int i;

    if (NULL != this->gt_bond_function)
    {
        sfree(this->gt_bond_function);
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
    this->gt_bond_ftype    = i;
    this->gt_bond_function = strdup(fn);
}

char *Poldata::get_bond_function()
{
    return this->gt_bond_function;
}

void Poldata::set_angle_function( char *fn)
{
    int i;

    if (NULL != this->gt_angle_function)
    {
        sfree(this->gt_angle_function);
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
    this->gt_angle_ftype    = i;
    this->gt_angle_function = strdup(fn);
}

char *Poldata::get_angle_function()
{
    return this->gt_angle_function;
}

int Poldata::get_bond_ftype()
{
    return this->gt_bond_ftype;
}

int Poldata::get_angle_ftype()
{
    return this->gt_angle_ftype;
}

int Poldata::get_dihedral_ftype( int egd)
{
    return this->gt_dihedral_ftype[egd];
}

void Poldata::set_dihedral_function( int egd, char *func)
{
    int i;

    if (NULL != this->gt_dihedral_function[egd])
    {
        sfree(this->gt_dihedral_function[egd]);
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
    this->gt_dihedral_ftype[egd]    = i;
    this->gt_dihedral_function[egd] = strdup(func);
}

char *Poldata::get_dihedral_function( int egd)
{
    return this->gt_dihedral_function[egd];
}

void Poldata::set_ptype_polarizability( const char *ptype,
                                          double polarizability, double sig_pol)
{
    t_ptype *sp;
    int      i;

    for (i = 0; (i < this->nptype); i++)
    {
        if (strcmp(ptype, this->ptype[i].type) == 0)
        {
            sp                 = &(this->ptype[i]);
            sp->polarizability = polarizability;
            sp->sig_pol        = sig_pol;
            break;
        }
    }
    if (i == this->nptype)
    {
        fprintf(stderr, "No such ptype %s when trying to set the polarizability.\n", ptype);
    }
}

char *Poldata::get_polar_unit()
{
    return this->alexandria_polar_unit;
}

char *Poldata::get_polar_ref()
{
    return this->alexandria_polar_ref;
}

char *Poldata::get_length_unit()
{
    if (NULL == this->gt_length_unit)
    {
        gmx_fatal(FARGS, "No length unit in %s",
                  (NULL != this->filename) ? this->filename : "unknown");
    }

    return this->gt_length_unit;
}

void Poldata::set_polar_unit( const char *polar_unit)
{
    this->alexandria_polar_unit   = strdup(polar_unit);
}

void Poldata::set_polar_ref( const char *polar_ref)
{
    this->alexandria_polar_ref   = strdup(polar_ref);
}

char *Poldata::get_force_field()
{
    return this->alexandria_forcefield;
}

void Poldata::set_force_field( const char *forcefield)
{
    this->alexandria_forcefield   = strdup(forcefield);
}

void Poldata::set_length_unit( const char *length_unit)
{
    this->gt_length_unit   = strdup(length_unit);
}

char *Poldata::get_geometry( char *gt_brule)
{
    int i;

    if (gt_brule)
    {
        for (i = 0; (i < this->nbrule); i++)
        {
            if (strcmp(this->brule[i].rule, gt_brule) == 0)
            {
                return this->brule[i].geometry;
            }
        }
    }

    return NULL;
}

char *Poldata::get_desc( char *atype)
{
    int i;

    if (atype)
    {
        for (i = 0; (i < this->nalexandria); i++)
        {
            if (strcmp(this->alexandria[i].type, atype) == 0)
            {
                return this->alexandria[i].desc;
            }
        }
    }

    return NULL;
}

gmx_bool Poldata::strcasestr_start(char *needle, char *haystack)
{
    char *ptr;

    ptr = strcasestr(haystack, needle);

    return (ptr == haystack);
}

int Poldata::count_neighbors(t_brule *brule, int nbond, char *nbhybrid[], int *score)
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

char **Poldata::get_bonding_rules( char *elem,
                                     int nbond, char *neighbors[],
                                     const char *geometry,
                                     int iAromatic)
{
    unsigned int    nnb;
    int             i, nptr = 0, best = -1, score;
    char          **ptr = NULL;

    for (i = 0; (i < this->nbrule); i++)
    {
        nnb = count_neighbors(&(this->brule[i]), nbond, neighbors, &score);
        if ((strcmp(this->brule[i].elem, elem) == 0) &&
            (strcasecmp(this->brule[i].geometry, geometry) == 0) &&
            (nbond == this->brule[i].numbonds) &&
            ((iAromatic >= 0 && (iAromatic == this->brule[i].iAromatic)) ||
             (iAromatic < 0)) &&
            (nnb == this->brule[i].nb.size()))
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
                ptr[nptr-1] = this->brule[i].rule;
                best        = score;
            }
        }
    }
    srenew(ptr, ++nptr);
    ptr[nptr-1] = NULL;

    return ptr;
}

int Poldata::bonding_rule_valence( char *gt_brule, double *valence)
{
    int i;

    for (i = 0; (i < this->nbrule); i++)
    {
        if (strcasecmp(gt_brule, this->brule[i].rule) == 0)
        {
            *valence = this->brule[i].valence;
            return 1;
        }
    }
    return 0;
}

int Poldata::get_ptype_pol( const char *ptype,
                              double *polar, double *sig_pol)
{
    int j;

    for (j = 0; (j < this->nptype); j++)
    {
        if (strcmp(ptype, this->ptype[j].type) == 0)
        {
            if (NULL != polar)
            {
                *polar   = this->ptype[j].polarizability;
            }
            if (NULL != sig_pol)
            {
                *sig_pol = this->ptype[j].sig_pol;
            }
            return 1;
        }
    }
    return 0;
}

int Poldata::get_atype_pol( const char *atype,
                              double *polar, double *sig_pol)
{
    int i;

    for (i = 0; (i < this->nalexandria); i++)
    {
        if (strcmp(atype, this->alexandria[i].type) == 0)
        {
            return get_ptype_pol( this->alexandria[i].ptype, polar, sig_pol);
        }
    }
    return 0;
}


int Poldata::get_atype_ref_enthalpy( const char *atype,
                                       double *Href)
{
    int i;

    for (i = 0; (i < this->nalexandria); i++)
    {
        if (strcmp(atype, this->alexandria[i].type) == 0)
        {
            *Href = this->alexandria[i].ref_enthalpy;
            return 1;
        }
    }
    return 0;
}

char *Poldata::ptype_to_miller( const char *ptype)
{
    int i;

    for (i = 0; (i < this->nptype); i++)
    {
        if (strcmp(ptype, this->ptype[i].type) == 0)
        {
            return this->ptype[i].miller;
        }
    }
    return NULL;
}

char *Poldata::ptype_to_bosque( const char *ptype)
{
    int i;

    for (i = 0; (i < this->nptype); i++)
    {
        if (strcmp(ptype, this->ptype[i].type) == 0)
        {
            return this->ptype[i].bosque;
        }
    }
    return NULL;
}

int Poldata::get_ptype(
                          char        **ptype,
                          char        **miller,
                          char        **bosque,
                          double       *polarizability,
                          double       *sig_pol)
{
    t_ptype *sp;

    if (this->nptype_c < this->nptype)
    {
        sp = &(this->ptype[this->nptype_c]);
        assign_scal(polarizability, sp->polarizability);
        assign_scal(sig_pol, sp->sig_pol);
        assign_str(ptype, sp->type);
        assign_str(miller, sp->miller);
        assign_str(bosque, sp->bosque);
        this->nptype_c++;
        return 1;
    }
    else
    {
        this->nptype_c = 0;
    }

    return 0;
}

int Poldata::get_atype(
                          char        **elem,
                          char        **desc,
                          char        **atype,
                          char        **ptype,
                          char        **btype,
                          char        **vdwparams,
                          double       *ref_enthalpy)
{
    t_ffatype *sp;

    if (this->nalexandria_c < this->nalexandria)
    {
        sp = &(this->alexandria[this->nalexandria_c]);
        assign_str(elem, sp->elem);
        assign_str(desc, sp->desc);
        assign_str(atype, sp->type);
        assign_str(ptype, sp->ptype);
        assign_str(btype, sp->btype);
        assign_str(vdwparams, sp->vdwparams);
        *ref_enthalpy = sp->ref_enthalpy;
        this->nalexandria_c++;
        return 1;
    }
    else
    {
        this->nalexandria_c = 0;
    }

    return 0;
}

const char *Poldata::atype_to_ptype( const char *atype)
{
    int i;

    for (i = 0; (i < this->nalexandria); i++)
    {
        if (strcmp(this->alexandria[i].type, atype) == 0)
        {
            return this->alexandria[i].ptype;
        }
    }
    return NULL;
}

const char *Poldata::atype_to_btype( const char *atype)
{
    int i;

    for (i = 0; (i < this->nalexandria); i++)
    {
        if (strcmp(this->alexandria[i].type, atype) == 0)
        {
            return this->alexandria[i].btype;
        }
    }
    return NULL;
}

int Poldata::search_atype(
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

    for (i = 0; (i < this->nalexandria); i++)
    {
        if (strcmp(key, this->alexandria[i].type) == 0)
        {
            break;
        }
    }

    if (i < this->nalexandria)
    {
        sp = &(this->alexandria[i]);
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

double Poldata::elem_get_max_valence( char *elem)
{
    double mv = 0;
    int    i;

    for (i = 0; (i < this->nbrule); i++)
    {
        if ((0 == gmx_strcasecmp(this->brule[i].elem, elem)) &&
            (mv < this->brule[i].valence))
        {
            mv = this->brule[i].valence;
        }
    }
    return mv;
}

double *Poldata::elem_get_bondorders( char *elem1, char *elem2,
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
    for (i = 0; (i < this->ngt_bond); i++)
    {
        if (0 == strlen(this->gt_bond[i].elem1))
        {
            for (j = 0; (j < this->nalexandria); j++)
            {
                if (strcmp(this->alexandria[j].type, this->gt_bond[i].atom2) == 0)
                {
                    strcpy(this->gt_bond[i].elem1, this->alexandria[j].elem);
                }
            }
        }
        if (0 == strlen(this->gt_bond[i].elem2))
        {
            for (j = 0; (j < this->nalexandria); j++)
            {
                if (strcmp(this->alexandria[j].type, this->gt_bond[i].atom2) == 0)
                {
                    strcpy(this->gt_bond[i].elem2, this->alexandria[j].elem);
                }
            }
        }
        ba1 = this->gt_bond[i].elem1;
        ba2 = this->gt_bond[i].elem2;
        if (((strcmp(ba1, elem1) == 0) && (strcmp(ba2, elem2) == 0)) ||
            ((strcmp(ba1, elem2) == 0) && (strcmp(ba2, elem1) == 0)))
        {
            dev = fabs((this->gt_bond[i].length - distance)/this->gt_bond[i].length);
            if (dev < toler)
            {
                for (k = 0; (k < nbo); k++)
                {
                    if (this->gt_bond[i].bondorder == bo[k])
                    {
                        break;
                    }
                }
                if (k == nbo)
                {
                    srenew(bo, nbo+2);
                    bo[nbo]   = this->gt_bond[i].bondorder;
                    bo[nbo+1] = 0;
                    nbo++;
                }
            }
        }
    }
    return bo;
}

int Poldata::elem_is_bond( char *elem1, char *elem2,
                             double distance, double toler)
{
    char  *ba1, *ba2;
    double dev, dev_best = 100000;
    int    j, i;

    if ((NULL == elem1) || (NULL == elem2))
    {
        return 0;
    }
    for (i = 0; (i < this->ngt_bond); i++)
    {
        if (0 == strlen(this->gt_bond[i].elem1))
        {
            for (j = 0; (j < this->nalexandria); j++)
            {
                if (strcmp(this->alexandria[j].type, this->gt_bond[i].atom2) == 0)
                {
                    strcpy(this->gt_bond[i].elem1, this->alexandria[j].elem);
                }
            }
        }
        if (0 == strlen(this->gt_bond[i].elem2))
        {
            for (j = 0; (j < this->nalexandria); j++)
            {
                if (strcmp(this->alexandria[j].type, this->gt_bond[i].atom2) == 0)
                {
                    strcpy(this->gt_bond[i].elem2, this->alexandria[j].elem);
                }
            }
        }
        ba1 = this->gt_bond[i].elem1;
        ba2 = this->gt_bond[i].elem2;
        if (((strcmp(ba1, elem1) == 0) && (strcmp(ba2, elem2) == 0)) ||
            ((strcmp(ba1, elem2) == 0) && (strcmp(ba2, elem1) == 0)))
        {
            dev = fabs((this->gt_bond[i].length - distance)/this->gt_bond[i].length);
            if (dev < dev_best)
            {
                dev_best = dev;
            }
        }
    }
    return (dev_best < toler);
}

 int lo_gtb_comp(t_gt_bond *ba, t_gt_bond *bb)
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


 int Poldata::gtb_comp(const void *a, const void *b)
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

 t_gt_bond *Poldata::search_bond( char *atom1, char *atom2,
                              double bondorder)
{
    t_gt_bond key, *gt_b;
    int       i;

    key.atom1     = atom1;
    key.atom2     = atom2;
    key.bondorder = bondorder;

    gt_b = (t_gt_bond *) bsearch(&key, this->gt_bond, this->ngt_bond, sizeof(key), gtb_comp);
    if (NULL != gt_b)
    {
        i = gt_b - this->gt_bond;
        while ((i > 0) && (lo_gtb_comp(&(this->gt_bond[i-1]), &(this->gt_bond[i])) == 0))
        {
            i--;
        }
        gt_b = &(this->gt_bond[i]);
    }
    return gt_b;
}

double Poldata::atype_bondorder( char *atype1, char *atype2,
                                   double distance, double toler)
{
    double     dev, dev_best = 100000;
    int        i, i_best = -1;
    t_gt_bond *gt_b;

    if ((NULL == atype1) || (NULL == atype2))
    {
        return 0.0;
    }
    gt_b = search_bond( atype1, atype2, 0);
    if (NULL != gt_b)
    {
        i = gt_b - this->gt_bond;
        do
        {
            dev = fabs(this->gt_bond[i].length - distance);
            if (dev < dev_best)
            {
                dev_best = dev;
                i_best   = i;
            }
            i++;
        }
        while ((i < this->ngt_bond) &&
               (0 == lo_gtb_comp(&(this->gt_bond[i]), &(this->gt_bond[i-1]))));
    }
    if (dev_best < toler)
    {
        return this->gt_bond[i_best].bondorder;
    }

    return 0.0;
}

void Poldata::add_miller(
                            char         *miller,
                            int           atomnumber,
                            double        tau_ahc,
                            double        alpha_ahp)
{
    t_miller *mil;

    this->nmiller++;
    srenew(this->miller, this->nmiller);
    mil             = &(this->miller[this->nmiller-1]);
    mil->miller     = strdup(miller);
    mil->atomnumber = atomnumber;
    mil->tau_ahc    = tau_ahc;
    mil->alpha_ahp  = alpha_ahp;
}

void Poldata::set_miller_units( char *tau_unit, char *ahp_unit)
{
    this->miller_tau_unit = strdup(tau_unit);
    this->miller_ahp_unit = strdup(ahp_unit);
}

void Poldata::get_miller_units( char **tau_unit,
                                  char **ahp_unit)
{
    assign_str(tau_unit, this->miller_tau_unit);
    assign_str(ahp_unit, this->miller_ahp_unit);
}

int Poldata::get_miller(
                           char        **miller,
                           int          *atomnumber,
                           double       *tau_ahc,
                           double       *alpha_ahp)
{
    t_miller *mil;
    int       i;

    i = this->nmiller_c;

    if (i < this->nmiller)
    {
        mil = &(this->miller[i]);
        assign_str(miller, mil->miller);
        assign_scal(atomnumber, mil->atomnumber);
        assign_scal(tau_ahc, mil->tau_ahc);
        assign_scal(alpha_ahp, mil->alpha_ahp);
        this->nmiller_c++;

        return 1;
    }
    else
    {
        this->nmiller_c = 0;
    }

    return 0;
}

int Poldata::get_miller_pol(
                               char         *miller,
                               int          *atomnumber,
                               double       *tau_ahc,
                               double       *alpha_ahp)
{
    t_miller *mil;
    int       i;

    for (i = 0; (i < this->nmiller); i++)
    {
        if (strcmp(miller, this->miller[i].miller) == 0)
        {
            mil = &(this->miller[i]);
            assign_scal(atomnumber, mil->atomnumber);
            assign_scal(tau_ahc, mil->tau_ahc);
            assign_scal(alpha_ahp, mil->alpha_ahp);

            return 1;
        }
    }

    return 0;
}

void Poldata::add_bosque(
                            char         *bosque,
                            double        polarizability)
{
    t_bosque *bs;

    this->nbosque++;
    srenew(this->bosque, this->nbosque);
    bs                 = &(this->bosque[this->nbosque-1]);
    bs->bosque         = strdup(bosque);
    bs->polarizability = polarizability;
}

int Poldata::get_bosque(
                           char        **bosque,
                           double       *polarizability)
{
    if (this->nbosque_c < this->nbosque)
    {
        assign_str(bosque, this->bosque[this->nbosque_c].bosque);
        assign_scal(polarizability, this->bosque[this->nbosque_c].polarizability);
        this->nbosque_c++;

        return 1;
    }
    else
    {
        this->nbosque_c = 0;
    }

    return 0;
}

int Poldata::get_bosque_pol(
                               char         *bosque,
                               double       *polarizability)
{
    int i;

    for (i = 0; (i < this->nbosque); i++)
    {
        if (strcasecmp(bosque, this->bosque[i].bosque) == 0)
        {
            *polarizability = this->bosque[i].polarizability;
            return 1;
        }
    }
    return 0;
}

void Poldata::set_bosque_unit( char *polar_unit)
{
    this->bosque_polar_unit   = strdup(polar_unit);
}

char *Poldata::get_bosque_unit()
{
    return this->bosque_polar_unit;
}


int Poldata::search_bondtype( char *atom)
{
    int j;

    for (j = 0; (j < this->nbtype); j++)
    {
        if (strcmp(this->btype[j], atom) == 0)
        {
            return j;
        }
    }
    return -1;
}

/*
 * gt_bond stuff
 */
int Poldata::set_bond_params( char *atom1, char *atom2,
                                double length, double sigma, int ntrain,
                                double bondorder, char *params)
{
    t_gt_bond *gt_b;
    int        i;

    for (i = 0; (i < this->ngt_bond); i++)
    {
        gt_b = &(this->gt_bond[i]);
        if (((((strcmp(gt_b->atom1, atom1) == 0) &&
               (strcmp(gt_b->atom2, atom2) == 0)) ||
              ((strcmp(gt_b->atom1, atom2) == 0) &&
               (strcmp(gt_b->atom2, atom1) == 0)))) &&
            ((bondorder == 0) || (gt_b->bondorder == bondorder)))
        {
            break;
        }
    }
    if (i < this->ngt_bond)
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

int Poldata::add_bond( char *atom1, char *atom2,
                         double length, double sigma, int ntrain,
                         double bondorder, char *params)
{
    t_gt_bond *gt_b;
    int        a1, a2;

    if (-1 == (a1 = search_bondtype( atom1)))
    {
        return 0;
    }
    if (-1 == (a2 = search_bondtype( atom2)))
    {
        return 0;
    }
    if (Poldata::set_bond_params( atom1, atom2, length, sigma, ntrain,
                                    bondorder, params) == 0)
    {
        this->ngt_bond++;
        srenew(this->gt_bond, this->ngt_bond);
        gt_b            = &(this->gt_bond[this->ngt_bond-1]);
        gt_b->atom1     = strdup(atom1);
        strncpy(gt_b->elem1, this->alexandria[a1].elem, sizeof(gt_b->elem1)-1);
        gt_b->atom2     = strdup(atom2);
        strncpy(gt_b->elem2, this->alexandria[a2].elem, sizeof(gt_b->elem2)-1);
        gt_b->length    = length;
        gt_b->sigma     = sigma;
        gt_b->ntrain    = ntrain;
        gt_b->bondorder = bondorder;
        gt_b->params    = strdup(params);
        qsort(this->gt_bond, this->ngt_bond, sizeof(this->gt_bond[0]), gtb_comp);
    }
    return 1;
}

int Poldata::get_bond( char **atom1, char **atom2,
                         double *length, double *sigma, int *ntrain,
                         double *bondorder, char **params)
{
    t_gt_bond *gt_b;

    if (this->ngt_bond_c < this->ngt_bond)
    {
        gt_b = &(this->gt_bond[this->ngt_bond_c]);
        assign_str(atom1, gt_b->atom1);
        assign_str(atom2, gt_b->atom2);
        assign_scal(length, gt_b->length);
        assign_scal(sigma, gt_b->sigma);
        assign_scal(ntrain, gt_b->ntrain);
        assign_scal(bondorder, gt_b->bondorder);
        assign_str(params, gt_b->params);
        this->ngt_bond_c++;

        return this->ngt_bond_c;
    }
    this->ngt_bond_c = 0;

    return 0;
}

int Poldata::search_bond( char *atom1, char *atom2,
                            double *length, double *sigma, int *ntrain,
                            double *bondorder, char **params)
{
    t_gt_bond *gt_b;

    if ((NULL == atom1) || (NULL == atom2))
    {
        return 0;
    }
    gt_b = search_bond( atom1, atom2, 0);
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

            return 1+(gt_b - this->gt_bond);
        }
    }
    return 0;
}

/*
 * gt_angle stuff
 */
int Poldata::set_angle_params( char *atom1, char *atom2,
                                 char *atom3, double angle, double sigma, int ntrain,
                                 char *params)
{
    t_gt_angle *gt_b;
    int         i;

    for (i = 0; (i < this->ngt_angle); i++)
    {
        gt_b = &(this->gt_angle[i]);
        if ((strcmp(gt_b->atom2, atom2) == 0) &&
            (((strcmp(gt_b->atom1, atom1) == 0) &&
              (strcmp(gt_b->atom3, atom3) == 0)) ||
             ((strcmp(gt_b->atom1, atom3) == 0) &&
              (strcmp(gt_b->atom3, atom1) == 0))))
        {
            break;
        }
    }
    if (i < this->ngt_angle)
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

int Poldata::add_angle(
                          char *atom1, char *atom2,
                          char *atom3, double angle, double sigma,
                          int ntrain, char *params)
{
    t_gt_angle *gt_b;

    if ((-1 == search_bondtype( atom1)) ||
        (-1 == search_bondtype( atom2)) ||
        (-1 == search_bondtype( atom3)))
    {
        return 0;
    }

    if (0 == set_angle_params( atom1, atom2, atom3, angle, sigma, ntrain, params))
    {
        this->ngt_angle++;
        srenew(this->gt_angle, this->ngt_angle);
        gt_b          = &(this->gt_angle[this->ngt_angle-1]);
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

int Poldata::get_angle( char **atom1, char **atom2,
                          char **atom3, double *angle, double *sigma,
                          int *ntrain, char **params)
{
    t_gt_angle *gt_b;

    if (this->ngt_angle_c < this->ngt_angle)
    {
        gt_b = &(this->gt_angle[this->ngt_angle_c]);
        assign_str(atom1, gt_b->atom1);
        assign_str(atom2, gt_b->atom2);
        assign_str(atom3, gt_b->atom3);
        assign_scal(angle, gt_b->angle);
        assign_scal(sigma, gt_b->sigma);
        assign_scal(ntrain, gt_b->ntrain);
        assign_str(params, gt_b->params);
        this->ngt_angle_c++;

        return this->ngt_angle_c;
    }
    this->ngt_angle_c = 0;

    return 0;
}

int Poldata::search_angle( char *atom1, char *atom2,
                             char *atom3, double *angle, double *sigma,
                             int *ntrain, char **params)
{
    t_gt_angle *gt_b;
    int         i;

    if ((NULL == atom1) || (NULL == atom2) || (NULL == atom3))
    {
        return 0;
    }
    for (i = 0; (i < this->ngt_angle); i++)
    {
        gt_b = &(this->gt_angle[i]);
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

void Poldata::set_angle_unit( char *angle_unit)
{
    this->gt_angle_unit   = strdup(angle_unit);
}

char *Poldata::get_angle_unit()
{
    return this->gt_angle_unit;
}

/*
 * gt_dihedral stuff
 */
int gtd_comp(const void *a, const void *b)
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

t_gt_dihedral *Poldata::search_dihedral( int egd,
                                      char *atom1, char *atom2,
                                      char *atom3, char *atom4)
{
    t_gt_dihedral gt_a, *gt_res, *gt_dptr;
    int           nd;

    if ((NULL == atom1) || (NULL == atom2) || (NULL == atom3) || (NULL == atom4))
    {
        return NULL;
    }
    gt_dptr    = this->gt_dihedral[egd];
    nd         = this->ngt_dihedral[egd];
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

int Poldata::set_dihedral_params( int egd,
                                    char *atom1, char *atom2,
                                    char *atom3, char *atom4,
                                    double dihedral, double sigma, int ntrain,
                                    char *params)
{
    t_gt_dihedral *gt_b;


    gt_b = search_dihedral( egd, atom1, atom2, atom3, atom4);
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

int Poldata::add_dihedral( int egd,
                             char *atom1, char *atom2,
                             char *atom3, char *atom4, double dihedral,
                             double sigma, int ntrain, char *params)
{
    t_gt_dihedral *gt_b;

    if ((-1 == search_bondtype( atom1)) ||
        (-1 == search_bondtype( atom2)) ||
        (-1 == search_bondtype( atom3)) ||
        (-1 == search_bondtype( atom4)))
    {
        return 0;
    }

    if (0 == Poldata::set_dihedral_params( egd, atom1, atom2,
                                             atom3, atom4, dihedral,
                                             sigma, ntrain, params))
    {
        this->ngt_dihedral[egd]++;
        srenew(this->gt_dihedral[egd], this->ngt_dihedral[egd]);
        gt_b           = &(this->gt_dihedral[egd][this->ngt_dihedral[egd]-1]);
        gt_b->atom1    = strdup(atom1);
        gt_b->atom2    = strdup(atom2);
        gt_b->atom3    = strdup(atom3);
        gt_b->atom4    = strdup(atom4);
        gt_b->dihedral = dihedral;
        gt_b->sigma    = sigma;
        gt_b->ntrain   = ntrain;
        gt_b->params   = strdup(params);
        qsort(this->gt_dihedral[egd], this->ngt_dihedral[egd], sizeof(this->gt_dihedral[egd][0]),
              gtd_comp);
    }
    return 1;
}

int Poldata::get_dihedral( int egd,
                             char **atom1, char **atom2,
                             char **atom3, char **atom4, double *dihedral,
                             double *sigma, int *ntrain, char **params)
{
    t_gt_dihedral *gt_b;

    if (this->ngt_dihedral_c[egd] < this->ngt_dihedral[egd])
    {
        gt_b = &(this->gt_dihedral[egd][this->ngt_dihedral_c[egd]]);
        assign_str(atom1, gt_b->atom1);
        assign_str(atom2, gt_b->atom2);
        assign_str(atom3, gt_b->atom3);
        assign_str(atom4, gt_b->atom4);
        assign_scal(dihedral, gt_b->dihedral);
        assign_scal(sigma, gt_b->sigma);
        assign_scal(ntrain, gt_b->ntrain);
        assign_str(params, gt_b->params);
        this->ngt_dihedral_c[egd]++;

        return this->ngt_dihedral_c[egd];
    }
    this->ngt_dihedral_c[egd] = 0;

    return 0;
}

int Poldata::search_dihedral( int egd,
                                char *atom1, char *atom2,
                                char *atom3, char *atom4,
                                double *dihedral, double *sigma,
                                int *ntrain, char **params)
{
    t_gt_dihedral *gt_res;

    gt_res = search_dihedral( egd, atom1, atom2, atom3, atom4);
    if (NULL != gt_res)
    {
        assign_scal(dihedral, gt_res->dihedral);
        assign_scal(sigma, gt_res->sigma);
        assign_scal(ntrain, gt_res->ntrain);
        assign_str(params, gt_res->params);

        return 1 + (int) (gt_res - this->gt_dihedral[egd]);
    }
    return 0;
}

void Poldata::add_symcharges( char *central,
                                char *attached, int numattach)
{
    t_symcharges *sc;
    int           i;

    for (i = 0; (i < this->nsymcharges); i++)
    {
        sc = &(this->symcharges[i]);
        if ((strcasecmp(sc->central, central) == 0) &&
            (strcasecmp(sc->attached, attached) == 0) &&
            (sc->numattach == numattach))
        {
            break;
        }
    }
    if (i == this->nsymcharges)
    {
        this->nsymcharges++;
        srenew(this->symcharges, this->nsymcharges);
        sc              = &(this->symcharges[i]);
        sc->central     = strdup(central);
        sc->attached    = strdup(attached);
        sc->numattach   = numattach;
    }
}

int Poldata::get_symcharges( char **central,
                               char **attached, int *numattach)
{
    t_symcharges *sc;

    if (this->nsymcharges_c < this->nsymcharges)
    {
        sc = &(this->symcharges[this->nsymcharges_c]);
        assign_str(central, sc->central);
        assign_str(attached, sc->attached);
        assign_scal(numattach, sc->numattach);
        this->nsymcharges_c++;

        return 1;
    }
    this->nsymcharges_c = 0;

    return 0;
}

int Poldata::search_symcharges( char *central,
                                  char *attached, int numattach)
{
    t_symcharges *sc;
    int           i;

    for (i = 0; (i < this->nsymcharges); i++)
    {
        sc = &(this->symcharges[i]);
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
t_eemprops *Poldata::get_eep(ChargeDistributionModel eqd_model,
                           const char *name)
{
    int i;

    for (i = 0; (i < this->nep); i++)
    {
        if ((strcasecmp(this->eep[i].name, name) == 0) &&
            (this->eep[i].eqd_model == eqd_model))
        {
            return &(this->eep[i]);
        }
    }
    return NULL;
}

void Poldata::set_eemprops(ChargeDistributionModel eqd_model, char *name,
                              double J0, double chi0, char *zeta, char *q, char *row)
{
   
    t_eemprops              *eep;
    std::vector<std::string> sz, sq, sr;

    eep = this->get_eep(eqd_model, name);
    if (NULL == eep)
    {
      this->nep++;
        srenew(this->eep,this->nep );
        eep = &(this->eep[this->nep-1]);
    }
    eep->eqd_model = eqd_model;
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

int Poldata::get_eemprops(
                             ChargeDistributionModel *eqd_model, char **name,
                             double *J0, double *chi0, char **zeta, char **q, char **row)
{
    if (this->nep_c < this->nep)
    {
        assign_scal(eqd_model, this->eep[this->nep_c].eqd_model);
        assign_str(name, this->eep[this->nep_c].name);
        assign_scal(J0, this->eep[this->nep_c].J0);
        assign_str(zeta, this->eep[this->nep_c].zetastr);
        assign_str(q, this->eep[this->nep_c].qstr);
        assign_str(row, this->eep[this->nep_c].rowstr);
        assign_scal(chi0, this->eep[this->nep_c].chi0);
        this->nep_c++;
        return 1;
    }
    else
    {
        this->nep_c = 0;
        return 0;
    }
}

int Poldata::get_numprops( ChargeDistributionModel eqd_model)
{
    int i, n = 0;

    for (i = 0; (i < this->nep); i++)
    {
        if (this->eep[i].eqd_model == eqd_model)
        {
            n++;
        }
    }

    return n;
}

int Poldata::have_pol_support( const char *atype)
{
    int i;

    for (i = 0; (i < this->nalexandria); i++)
    {
        if (strcmp(atype, this->alexandria[i].type) == 0)
        {
            return 1;
        }
    }
    return 0;
}

int Poldata::have_eem_support( ChargeDistributionModel eqd_model,
                                 const char *name,
                                 gmx_bool bAllowZeroParameters)
{
   
  t_eemprops  *eep  = get_eep(eqd_model, name);

    return (eep && (bAllowZeroParameters || ((eep->J0 > 0) && (eep->chi0 > 0))));
}

double Poldata::get_j00( ChargeDistributionModel eqd_model, char *name)
{
    t_eemprops  *eer;

    if ((eer = get_eep(eqd_model, name)) != NULL)
    {
        return eer->J0;
    }
    else
    {
        gmx_fatal(FARGS, "No J0 data for eqd_model %d and name %s",
                  eqd_model, name);
    }
    return -1;
}

char *Poldata::get_qstr( ChargeDistributionModel eqd_model, char *name)
{
    t_eemprops *eer;

    if ((eer = get_eep( eqd_model, name)) != NULL)
    {
        return eer->qstr;
    }
    return NULL;
}

char *Poldata::get_rowstr( ChargeDistributionModel eqd_model, char *name)
{
    t_eemprops *eer;

    if ((eer = get_eep( eqd_model, name)) != NULL)
    {
        return eer->rowstr;
    }
    return NULL;
}

int Poldata::get_row( ChargeDistributionModel eqd_model, char *name, int zz)
{
    t_eemprops *eer;

    if ((eer = get_eep( eqd_model, name)) != NULL)
    {
        range_check(zz, 0, eer->nzeta);
        return eer->row[zz];
    }
    return -1;
}

double Poldata::get_zeta( ChargeDistributionModel eqd_model, char *name, int zz)
{
    t_eemprops *eer;

    if ((eer = get_eep( eqd_model, name)) != NULL)
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

int Poldata::get_nzeta( ChargeDistributionModel eqd_model, char *name)
{
    t_eemprops *eer;

    if ((eer = get_eep( eqd_model, name)) != NULL)
    {
        return eer->nzeta;
    }
    return 0;
}

double Poldata::get_q( ChargeDistributionModel eqd_model, char *name, int zz)
{
    t_eemprops *eer;

    if ((eer = get_eep( eqd_model, name)) != NULL)
    {
        range_check(zz, 0, eer->nzeta);
        return eer->q[zz];
    }
    return -1;
}

double Poldata::get_chi0( ChargeDistributionModel eqd_model, char *name)
{
    t_eemprops *eer;

    if ((eer = get_eep( eqd_model, name)) != NULL)
    {
        return eer->chi0;
    }
    else
    {
        gmx_fatal(FARGS, "No chi0 data for eqd_model %d and name %s", eqd_model, name);
    }
    return -1;
}

void Poldata::set_epref( ChargeDistributionModel eqd_model, char *epref)
{
    int i;

    for (i = 0; (i < this->ner); i++)
    {
        if (this->epr[i].eqd_model == eqd_model)
        {
            if (this->epr[i].epref)
            {
                sfree(this->epr[i].epref);
            }
            this->epr[i].epref = strdup(epref);
            break;
        }
    }
    if (i == this->ner)
    {
        srenew(this->epr, ++this->ner);
        this->epr[i].eqd_model = eqd_model;
        this->epr[i].epref     = strdup(epref);
    }
}

char *Poldata::get_epref( ChargeDistributionModel eqd_model)
{
    int i;

    for (i = 0; (i < this->ner); i++)
    {
        if (this->epr[i].eqd_model == eqd_model)
        {
            return this->epr[i].epref;
        }
    }
    return NULL;
}

int Poldata::list_epref( ChargeDistributionModel *eqd_model, char **epref)
{
    if (this->ner_c < this->ner)
    {
        assign_scal(eqd_model, this->epr[this->ner_c].eqd_model);
        assign_str(epref, this->epr[this->ner_c].epref);
        this->ner_c++;
        return 1;
    }
    this->ner_c = 0;

    return 0;
}

void Poldata::comm_eemprops( t_commrec *cr)
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
            gmx_send_int(cr, i, this->nep);
            gmx_send(cr, i, this->eep, this->nep*sizeof(this->eep[0]));
        }
    }
    else
    {
        nep = gmx_recv_int(cr, 0);
        if (nep != this->nep)
        {
            gmx_fatal(FARGS, "Inconsistency in number of EEM parameters");
        }
        snew(ep, this->nep);
        gmx_recv(cr, 0, ep, this->nep*sizeof(ep[0]));
        for (i = 0; (i < this->nep); i++)
        {
            this->eep[i] = ep[i];
        }
        sfree(ep);
    }
    if (NULL != debug)
    {
        fprintf(debug, "  EEP  Atom      Chi      J00     Zeta\n");
        for (i = 0; (i < this->nep); i++)
        {
            fprintf(debug, "%5s %5s %8.3f %8.3f",
                    get_eemtype_name(this->eep[i].eqd_model),
                    this->eep[i].name, this->eep[i].chi0, this->eep[i].J0);
            for (j = 0; (j < this->eep[i].nzeta); j++)
            {
                fprintf(debug, " %8.3f", this->eep[i].zeta[j]);
            }
            fprintf(debug, "\n");
        }
    }
}

void Poldata::comm_force_parameters(t_commrec *cr)
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
            gmx_send_int(cr, i, this->nep);
            gmx_send(cr, i, this->eep, this->nep*sizeof(this->eep[0]));
        }
    }
    else
    {
        nep = gmx_recv_int(cr, 0);
        if (nep != this->nep)
        {
            gmx_fatal(FARGS, "Inconsistency in number of EEM parameters");
        }
        snew(ep, this->nep);
        gmx_recv(cr, 0, ep, this->nep*sizeof(ep[0]));
        for (i = 0; (i < this->nep); i++)
        {
            this->eep[i] = ep[i];
        }
        sfree(ep);
    }
    if (NULL != debug)
    {
        fprintf(debug, "  EEP  Atom      Chi      J00     Zeta\n");
        for (i = 0; (i < this->nep); i++)
        {
            fprintf(debug, "%5s %5s %8.3f %8.3f",
                    get_eemtype_name(this->eep[i].eqd_model),
                    this->eep[i].name, this->eep[i].chi0, this->eep[i].J0);
            for (j = 0; (j < this->eep[i].nzeta); j++)
            {
                fprintf(debug, " %8.3f", this->eep[i].zeta[j]);
            }
            fprintf(debug, "\n");
        }
    }
}

typedef struct {
    ChargeDistributionModel eqd;
    const char             *name, *ref;
    gmx_bool                bWeight;
} t_eemtype_props;

 t_eemtype_props eemtype_props[eqdNR] = {
    { eqdAXp,      "AXp",      "Maaren2014a",   FALSE },
    { eqdAXg,      "AXg",      "Maaren2014a",   TRUE },
    { eqdAXs,      "AXs",      "Maaren2014a",   TRUE },
    { eqdYang,     "Yang",     "Yang2006b",     TRUE },
    { eqdBultinck, "Bultinck", "Bultinck2002a", FALSE },
    { eqdRappe,    "Rappe",    "Rappe1991a",    TRUE }
};

ChargeDistributionModel Poldata::name2eemtype(const char *name)
{
    int i;

    for (i = 0; (i < eqdNR); i++)
    {
        if (strcasecmp(name, eemtype_props[i].name) == 0)
        {
            return eemtype_props[i].eqd;
        }
    }
    return eqdNR;
}

const char *Poldata::get_eemtype_name(ChargeDistributionModel eem)
{
    int i;

    for (i = 0; (i < eqdNR); i++)
    {
        if (eem == eemtype_props[i].eqd)
        {
            return eemtype_props[i].name;
        }
    }

    return NULL;
}

}
