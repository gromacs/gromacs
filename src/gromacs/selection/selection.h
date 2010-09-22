/*
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009, The GROMACS development team,
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
 */
/*! \file
 * \brief API for handling selection (the \c gmx_ana_selection_t structure and related functions).
 *
 * There should be no need to use the data structures or call the
 * functions in this file directly unless using the selection routines outside
 * the main trajectory analysis API.
 */
#ifndef GMX_SELECTION_SELECTION_H
#define GMX_SELECTION_SELECTION_H

#include "position.h"
#include "indexutil.h"

struct t_selelem;

/** Defines the type of covered fraction. */
typedef enum
{
    CFRAC_NONE,         /**< No covered fraction (everything covered). */
    CFRAC_SOLIDANGLE    /**< Fraction of a solid (3D) angle covered. */
} e_coverfrac_t;

/*! \brief
 * Describes a single selection.
 */
typedef struct gmx_ana_selection_t
{
    /** Name of the selection. */
    char                   *name;
    /** The actual selection string. */
    char                   *selstr;
    /** Selected positions. */
    gmx_ana_pos_t           p;
    /** Masses associated with the positions. */
    real                   *m;
    /** Charges associated with the positions. */
    real                   *q;
    /** Pointer to the index group that holds the selected atoms. */
    struct gmx_ana_index_t *g;
    /** TRUE if the value can change as a function of time. */
    gmx_bool                    bDynamic;
    /** Type of the covered fraction. */
    e_coverfrac_t           cfractype;
    /** TRUE if the covered fraction depends on the frame. */
    gmx_bool                    bCFracDyn;
    /** Covered fraction of the selection for the current frame. */
    real                    cfrac;
    /** The average covered fraction (over the trajectory). */
    real                    avecfrac;

    /*! \brief
     * Pointer to the root of the selection element tree (internal use only).
     *
     * \internal
     * This field is NULL if the selection has been loaded directly from an
     * index file.
     */
    struct t_selelem       *selelem;
    /** Original masses of all possible positions (internal use only). */
    real                   *orgm;
    /** Original charges of all possible positions (internal use only). */
    real                   *orgq;
} gmx_ana_selection_t;

namespace gmx
{

class SelectionCollection;

/*! \brief
 * Wrapper class for accessing selection information.
 */
class Selection
{
    public:
        Selection(t_selelem *elem, const char *selstr);

        const char *name() const  { return _sel.name; }
        const char *selectionText() const { return _sel.selstr; }
        e_index_t type() const { return _sel.p.m.type; }
        bool isDynamic() const { return _sel.bDynamic; }
        int posCount() const { return _sel.p.nr; }
        const gmx_ana_pos_t *positions() const { return &_sel.p; }
        const rvec &x(int i) const { return _sel.p.x[i]; }
        const rvec &v(int i) const { return _sel.p.v[i]; }
        const rvec &f(int i) const { return _sel.p.f[i]; }
        int refId(int i) const { return _sel.p.m.refid[i]; }
        int mapId(int i) const { return _sel.p.m.mapid[i]; }
        real mass(int i) const { return _sel.m[i]; }
        real charge(int i) const { return _sel.q[i]; }
        int atomCount(int i) const
            { return _sel.p.m.mapb.index[i+1] - _sel.p.m.mapb.index[i]; }
        const int *atomIndices(int i) const
            { return _sel.g ? _sel.g->index + _sel.p.m.mapb.index[i] : NULL; }
        int *mapIds() const { return _sel.p.m.mapid; }
        gmx_ana_index_t *indexGroup() const { return _sel.g; }
        real cfrac() const { return _sel.cfrac; }

        void setOriginalId(int i, int id) { _sel.p.m.orgid[i] = id; }
        bool initCoveredFraction(e_coverfrac_t type);

        void printInfo() const;
        void printDebugInfo(int nmaxind) const;

        gmx_ana_selection_t     _sel;

    private:
        ~Selection();

        friend class SelectionCollection;

        // Disallow copy and assign.
        Selection(const Selection &);
        void operator =(const Selection &);
};

} // namespace gmx

#endif
