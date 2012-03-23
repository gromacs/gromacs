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
#include <vector>

#include "../legacyheaders/typedefs.h"

#include "../fatalerror/gmxassert.h"
#include "../utility/arrayref.h"
#include "../utility/common.h"

#include "position.h"
#include "indexutil.h"
#include "selectionenums.h"

struct t_selelem;

namespace gmx
{

class SelectionEvaluator;
class SelectionCollection;
class SelectionCompiler;

class SelectionPosition;

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
        //! Returns atom indices of all atoms in the selection.
        ConstArrayRef<int> atomIndices() const
        {
            if (rawPositions_.g == NULL)
            {
                return ConstArrayRef<int>();
            }
            return ConstArrayRef<int>(rawPositions_.g->isize,
                                      rawPositions_.g->index);
        }
        //! Access a single position.
        SelectionPosition position(int i) const;
        /*! \brief
         * Sets the ID for the \p i'th position for use with
         * SelectionPosition::mappedId().
         *
         * This method is not part of SelectionPosition because that interface
         * only provides access to const data by design.
         */
        void setOriginalId(int i, int id) { rawPositions_.m.orgid[i] = id; }

        //! Returns whether the covered fraction can change between frames.
        bool isCoveredFractionDynamic() const { return bDynamicCoveredFraction_; }
        //! Returns the covered fraction for the current frame.
        real coveredFraction() const { return coveredFraction_; }

        //! Deprecated method for direct access to position data.
        const gmx_ana_pos_t *positions() const { return &rawPositions_; }
        //! Deprecated method for direct access to atom index data.
        gmx_ana_index_t *indexGroup() const { return rawPositions_.g; }

        // TODO: Remove direct access to the flags from the public interface.
        //! Returns true if the given flag is set.
        bool hasFlag(SelectionFlag flag) const { return flags_.test(flag); }
        //! Sets the flags for this selection.
        void setFlags(SelectionFlags flags) { flags_ = flags; }
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
         * Additional information about positions.
         *
         * This structure contains information about positions in the
         * selection that is not stored in ::gmx_ana_pos_t.
         */
        struct PositionInfo
        {
            //! Construct position information with unit mass and no charge.
            PositionInfo() : mass(1.0), charge(0.0) {}
            //! Construct position information with the given information.
            PositionInfo(real mass, real charge) : mass(mass), charge(charge) {}

            //! Total mass of atoms that make up the position.
            real                mass;
            //! Total charge of atoms that make up the position.
            real                charge;
        };

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
        //! Information associated with the current positions.
        std::vector<PositionInfo> posInfo_;
        //! Information for all possible positions.
        std::vector<PositionInfo> originalPosInfo_;
        SelectionFlags          flags_;
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
        /*! \brief
         * Needed for proper access to position information.
         */
        friend class SelectionPosition;

        GMX_DISALLOW_COPY_AND_ASSIGN(Selection);
};

/*! \brief
 * Wrapper object to access information about a single selected position.
 *
 * Default copy constructor and assignment operators are used, and work as
 * intended: the copy references the same position and works identically.
 *
 * Methods in this class do not throw.
 *
 * \inpublicapi
 * \ingroup module_selection
 */
class SelectionPosition
{
    public:
        /*! \brief
         * Constructs a wrapper object for given selection position.
         *
         * \param[in] sel    Selection from which the position is wrapped.
         * \param[in] index  Zero-based index of the position to wrap.
         *
         * Does not throw.  Asserts if \p index is out of range.
         */
        SelectionPosition(const Selection *sel, int index)
            : sel_(sel), i_(index)
        {
            GMX_ASSERT(index >= 0 && index < sel->posCount(),
                       "Invalid selection position index");
        }

        /*! \brief
         * Returns type of this position.
         *
         * Currently returns the same as Selection::type().
         */
        e_index_t type() const { return sel_->rawPositions_.m.type; }
        //! Returns coordinates for this position.
        const rvec &x() const
        {
            return sel_->rawPositions_.x[i_];
        }
        //! Returns whether velocity is available for this position.
        bool hasVelocity() const { return sel_->rawPositions_.v != NULL; }
        /*! \brief
         * Returns velocity for this position.
         *
         * Must not be called if hasVelocity() returns false.
         */
        const rvec &v() const
        {
            GMX_ASSERT(hasVelocity(), "Velocities accessed, but unavailable");
            return sel_->rawPositions_.v[i_];
        }
        //! Returns whether force is available for this position.
        bool hasForce() const { return sel_->rawPositions_.f != NULL; }
        /*! \brief
         * Returns velocity for this position.
         *
         * Must not be called if hasForce() returns false.
         */
        const rvec &f() const
        {
            GMX_ASSERT(hasForce(), "Forces accessed, but unavailable");
            return sel_->rawPositions_.f[i_];
        }
        /*! \brief
         * Returns total mass for this position.
         *
         * Returns the total mass of atoms that make up this position.
         * If there are not atoms associated or masses are not available,
         * returns unity.
         */
        real mass() const
        {
            return sel_->posInfo_[i_].mass;
        }
        /*! \brief
         * Returns total charge for this position.
         *
         * Returns the sum of charges of atoms that make up this position.
         * If there are not atoms associated or masses are not available,
         * returns zero.
         */
        real charge() const
        {
            return sel_->posInfo_[i_].charge;
        }
        //! Returns the number of atoms that make up this position.
        int atomCount() const
        {
            return sel_->rawPositions_.m.mapb.index[i_ + 1]
                 - sel_->rawPositions_.m.mapb.index[i_];
        }
        //! Return atom indices that make up this position.
        ConstArrayRef<int> atomIndices() const
        {
            if (sel_->rawPositions_.g == NULL)
            {
                return ConstArrayRef<int>();
            }
            int first = sel_->rawPositions_.m.mapb.index[i_];
            return ConstArrayRef<int>(atomCount(),
                                      &sel_->rawPositions_.g->index[first]);
        }
        /*! \brief
         * Returns reference ID for this position.
         *
         * For dynamic selections, this provides means to associate positions
         * across frames.  After compilation, these IDs are consequently
         * numbered starting from zero.  For each frame, the ID then reflects
         * the location of the position in the original array of positions.
         * If SelectionOption::dynamicMask() has been set for the parent
         * selection, the IDs for positions not present in the current
         * selection are set to -1, otherwise they are removed completely.
         *
         * Example:
         * If a dynamic selection consists of three positions, after
         * compilation refId() will return 0, 1, 2 for them, respectively.
         * If for a particular frame, only the first and the third are present,
         * refId() will return 0, 2.
         * If SelectionOption::dynamicMask() has been set, all three positions
         * can be accessed also in this case and refId() will return 0, -1, 2.
         */
        int refId() const
        {
            return sel_->rawPositions_.m.refid[i_];
        }
        /*! \brief
         * Returns mapped ID for this position.
         *
         * Returns ID of the position that corresponds to that set with
         * Selection::setOriginalId().
         *
         * If for an array \c id, \c setOriginalId(i, c[i]) has been called
         * for each \c i, then it always holds that
         * \c mappedId()==c[refId()].
         */
        int mappedId() const
        {
            return sel_->rawPositions_.m.mapid[i_];
        }

    private:
        const Selection        *sel_;
        int                     i_;
};


inline SelectionPosition
Selection::position(int i) const
{
    return SelectionPosition(this, i);
}

} // namespace gmx

#endif
