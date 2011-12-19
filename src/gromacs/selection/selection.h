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
 * \brief
 * Declares gmx::Selection.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \inpublicapi
 * \ingroup module_selection
 */
#ifndef GMX_SELECTION_SELECTION_H
#define GMX_SELECTION_SELECTION_H

#include "../legacyheaders/typedefs.h"

#include "position.h"
#include "indexutil.h"
#include "selectionenums.h"

struct t_selelem;

/*! \internal \brief
 * Describes a single selection.
 *
 * \ingroup module_selection
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
    /** true if the value can change as a function of time. */
    bool                    bDynamic;
    /** Type of the covered fraction. */
    e_coverfrac_t           cfractype;
    /** true if the covered fraction depends on the frame. */
    bool                    bCFracDyn;
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

class SelectionEvaluator;
class SelectionCollection;
class SelectionCompiler;

/*! \brief
 * Provides access to a single selection.
 *
 * \inpublicapi
 * \ingroup module_selection
 */
class Selection
{
    public:
        /*! \brief
         * Creates a new selection object.
         *
         * \param[in] elem   Root of the evaluation tree for this selection.
         * \param[in] selstr String that was parsed to produce this selection.
         */
        Selection(t_selelem *elem, const char *selstr);

        //! Returns the name of the selection.
        const char *name() const  { return _sel.name; }
        //! Returns the string that was parsed to produce this selection.
        const char *selectionText() const { return _sel.selstr; }
        //! Returns true if the size of the selection (posCount()) is dynamic.
        bool isDynamic() const { return _sel.bDynamic; }
        //! Returns the type of positions in the selection.
        e_index_t type() const { return _sel.p.m.type; }
        //! Number of positions in the selection.
        int posCount() const { return _sel.p.nr; }
        //! Returns the \p i'th position for the selection.
        const rvec &x(int i) const { return _sel.p.x[i]; }
        //! Returns the velocity for the \p i'th position.
        const rvec &v(int i) const { return _sel.p.v[i]; }
        //! Returns the force for the \p i'th position.
        const rvec &f(int i) const { return _sel.p.f[i]; }
        /*! \brief
         * Returns the reference ID for the \p i'th position.
         */
        int refId(int i) const { return _sel.p.m.refid[i]; }
        /*! \brief
         * Returns the mapped ID for the \p i'th position.
         */
        int mapId(int i) const { return _sel.p.m.mapid[i]; }
        //! Returns the mass for the \p i'th position.
        real mass(int i) const { return _sel.m[i]; }
        //! Returns the charge for the \p i'th position.
        real charge(int i) const { return _sel.q[i]; }
        //! Returns the number of atoms contributing to the \p i'th position.
        int atomCount(int i) const
            { return _sel.p.m.mapb.index[i+1] - _sel.p.m.mapb.index[i]; }
        //! Returns the atom indices contributing to the \p i'th position.
        const int *atomIndices(int i) const
            { return _sel.g ? _sel.g->index + _sel.p.m.mapb.index[i] : NULL; }
        //! Returns the covered fraction for the current frame.
        real cfrac() const { return _sel.cfrac; }
        //! Deprecated method for direct access to position data.
        const gmx_ana_pos_t *positions() const { return &_sel.p; }
        //! Deprecated method for direct access to atom index data.
        gmx_ana_index_t *indexGroup() const { return _sel.g; }
        //! Deprecated method for direct access to to mapped ID array.
        int *mapIds() const { return _sel.p.m.mapid; }

        //! Returns true if the given flag is set.
        bool hasFlag(SelectionFlag flag) const { return _flags.test(flag); }
        //! Sets the flags for this selection.
        void setFlags(SelectionFlags flags) { _flags = flags; }
        /*! \brief
         * Sets the ID for the \p i'th position for use with mapId().
         */
        void setOriginalId(int i, int id) { _sel.p.m.orgid[i] = id; }
        /*! \brief
         * Initializes information about covered fractions.
         *
         * \param[in] type Type of covered fraction required.
         * \returns   True if the covered fraction can be calculated for the
         *      selection.
         */
        bool initCoveredFraction(e_coverfrac_t type);

        /*! \brief
         * Prints out one-line description of the selection.
         *
         * \param[in] fp      Where to print the information.
         */
        void printInfo(FILE *fp) const;
        /*! \brief
         * Prints out extended information about the selection for debugging.
         *
         * \param[in] fp      Where to print the information.
         * \param[in] nmaxind Maximum number of values to print in lists
         *      (-1 = print all).
         */
        void printDebugInfo(FILE *fp, int nmaxind) const;

    private:
        ~Selection();

        void initializeMassesAndCharges(const t_topology *top);

        gmx_ana_selection_t     _sel;
        SelectionFlags          _flags;

        friend class SelectionCompiler;
        friend class SelectionCollection;
        friend class SelectionEvaluator;

        // Disallow copy and assign.
        Selection(const Selection &);
        void operator =(const Selection &);
};

} // namespace gmx

#endif
