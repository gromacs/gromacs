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

#include <string>

#include "../legacyheaders/typedefs.h"

#include "position.h"
#include "indexutil.h"
#include "selectionenums.h"

struct t_selelem;

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
        ~Selection();

        //! Returns the name of the selection.
        const char *name() const  { return name_.c_str(); }
        //! Returns the string that was parsed to produce this selection.
        const char *selectionText() const { return selectionText_.c_str(); }
        //! Returns true if the size of the selection (posCount()) is dynamic.
        bool isDynamic() const { return bDynamic_; }
        //! Returns the type of positions in the selection.
        e_index_t type() const { return rawPositions_.m.type; }
        //! Number of positions in the selection.
        int posCount() const { return rawPositions_.nr; }
        //! Total number of atoms in the selection.
        int atomCount() const
        {
            return rawPositions_.g != NULL ? rawPositions_.g->isize : 0;
        }
        //! Returns the \p i'th position for the selection.
        const rvec &x(int i) const { return rawPositions_.x[i]; }
        //! Returns the velocity for the \p i'th position.
        const rvec &v(int i) const { return rawPositions_.v[i]; }
        //! Returns the force for the \p i'th position.
        const rvec &f(int i) const { return rawPositions_.f[i]; }
        /*! \brief
         * Returns the reference ID for the \p i'th position.
         */
        int refId(int i) const { return rawPositions_.m.refid[i]; }
        /*! \brief
         * Returns the mapped ID for the \p i'th position.
         */
        int mapId(int i) const { return rawPositions_.m.mapid[i]; }
        //! Returns the mass for the \p i'th position.
        real mass(int i) const { return mass_[i]; }
        //! Returns the charge for the \p i'th position.
        real charge(int i) const { return charge_[i]; }
        //! Returns the number of atoms contributing to the \p i'th position.
        int atomCount(int i) const
        {
            return rawPositions_.m.mapb.index[i+1] - rawPositions_.m.mapb.index[i];
        }
        //! Returns the atom indices contributing to the \p i'th position.
        const int *atomIndices(int i) const
        {
            return rawPositions_.g != NULL
                ? rawPositions_.g->index + rawPositions_.m.mapb.index[i]
                : NULL;
        }
        //! Returns whether the covered fraction can change between frames.
        bool isCoveredFractionDynamic() const { return bDynamicCoveredFraction_; }
        //! Returns the covered fraction for the current frame.
        real coveredFraction() const { return coveredFraction_; }
        //! Deprecated method for direct access to position data.
        const gmx_ana_pos_t *positions() const { return &rawPositions_; }
        //! Deprecated method for direct access to atom index data.
        gmx_ana_index_t *indexGroup() const { return rawPositions_.g; }
        //! Deprecated method for direct access to to mapped ID array.
        int *mapIds() const { return rawPositions_.m.mapid; }

        //! Returns true if the given flag is set.
        bool hasFlag(SelectionFlag flag) const { return flags_.test(flag); }
        //! Sets the flags for this selection.
        void setFlags(SelectionFlags flags) { flags_ = flags; }
        /*! \brief
         * Sets the ID for the \p i'th position for use with mapId().
         */
        void setOriginalId(int i, int id) { rawPositions_.m.orgid[i] = id; }
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
        /*! \brief
         * Computes total masses and charges for all selection positions.
         *
         * \param[in]  top   Topology information.
         *
         * Computed values are cached, and need to be updated for dynamic
         * selections with refreshMassesAndCharges() after the selection has
         * been evaluated.  This is done by SelectionEvaluator.
         *
         * This function is called by SelectionCompiler.
         */
        void initializeMassesAndCharges(const t_topology *top);
        /*! \brief
         * Updates masses and charges after dynamic selection has been
         * evaluated.
         *
         * Called by SelectionEvaluator.
         */
        void refreshMassesAndCharges();
        /*! \brief
         * Updates the covered fraction after a selection has been evaluated.
         *
         * Called by SelectionEvaluator.
         */
        void updateCoveredFractionForFrame();
        /*! \brief
         * Computes average covered fraction after all frames have been evaluated.
         *
         * \param[in] nframes  Number of frames that have been evaluated.
         *
         * \p nframes should be equal to the number of calls to
         * updateCoveredFractionForFrame().
         * Called by SelectionEvaluator::evaluateFinal().
         */
        void computeAverageCoveredFraction(int nframes);
        /*! \brief
         * Restores position information to state it was in after compilation.
         *
         * Depends on SelectionCompiler storing the original atoms in the
         * \a rootElement_ object.
         * Called by SelectionEvaluator::evaluateFinal().
         */
        void restoreOriginalPositions();

        //! Name of the selection.
        std::string             name_;
        //! The actual selection string.
        std::string             selectionText_;
        //! Low-level representation of selected positions.
        gmx_ana_pos_t           rawPositions_;
        SelectionFlags          flags_;
        //! Masses associated with the current positions.
        real                   *mass_;
        //! Charges associated with the current positions.
        real                   *charge_;
        //! Original masses of all possible positions.
        real                   *originalMass_;
        //! Original charges of all possible positions.
        real                   *originalCharge_;
        //! Pointer to the root of the selection evaluation tree.
        t_selelem              *rootElement_;
        //! Type of the covered fraction.
        e_coverfrac_t           coveredFractionType_;
        //! Covered fraction of the selection for the current frame.
        real                    coveredFraction_;
        //! The average covered fraction (over the trajectory).
        real                    averageCoveredFraction_;
        //! true if the value can change as a function of time.
        bool                    bDynamic_;
        //! true if the covered fraction depends on the frame.
        bool                    bDynamicCoveredFraction_;

        /*! \brief
         * Needed for the compiler to access initializeMassesAndCharges().
         *
         * Currently the compiler also used rootElement_ directly for
         * simplicity, but does not modify it.
         */
        friend class SelectionCompiler;
        /*! \brief
         * Needed for the evaluator to access the private methods.
         */
        friend class SelectionEvaluator;

        // Disallow copy and assign.
        Selection(const Selection &);
        void operator =(const Selection &);
};

} // namespace gmx

#endif
