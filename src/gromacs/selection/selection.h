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
 * Declares gmx::Selection and supporting classes.
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

#include "../utility/arrayref.h"
#include "../utility/common.h"
#include "../utility/gmxassert.h"

#include "position.h"
#include "indexutil.h"
#include "selectionenums.h"

namespace gmx
{

class SelectionOptionStorage;
class SelectionTreeElement;

class Selection;
class SelectionPosition;

//! Container of selections used in public selection interfaces.
typedef std::vector<Selection> SelectionList;

namespace internal
{

/*! \internal \brief
 * Internal data for a single selection.
 *
 * This class is internal to the selection module, but resides in a public
 * header because of efficiency reasons: it allows frequently used access
 * methods in \ref Selection to be inlined.
 *
 * Methods in this class do not throw unless otherwise specified.
 *
 * \ingroup module_selection
 */
class SelectionData
{
    public:
        /*! \brief
         * Creates a new selection object.
         *
         * \param[in] elem   Root of the evaluation tree for this selection.
         * \param[in] selstr String that was parsed to produce this selection.
         * \throws    std::bad_alloc if out of memory.
         */
        SelectionData(SelectionTreeElement *elem, const char *selstr);
        ~SelectionData();

        //! Returns the string that was parsed to produce this selection.
        const char *selectionText() const { return selectionText_.c_str(); }
        //! Returns true if the size of the selection (posCount()) is dynamic.
        bool isDynamic() const { return bDynamic_; }
        //! Number of positions in the selection.
        int posCount() const { return rawPositions_.nr; }
        //! Returns the root of the evaluation tree for this selection.
        SelectionTreeElement &rootElement() { return rootElement_; }

        //! Returns whether the covered fraction can change between frames.
        bool isCoveredFractionDynamic() const { return bDynamicCoveredFraction_; }

        //! Returns true if the given flag is set.
        bool hasFlag(SelectionFlag flag) const { return flags_.test(flag); }
        //! Sets the flags for this selection.
        void setFlags(SelectionFlags flags) { flags_ = flags; }

        //! \copydoc Selection::initCoveredFraction()
        bool initCoveredFraction(e_coverfrac_t type);

        /*! \brief
         * Computes total masses and charges for all selection positions.
         *
         * \param[in] top   Topology information.
         * \throws    std::bad_alloc if out of memory.
         *
         * Computed values are cached, and need to be updated for dynamic
         * selections with refreshMassesAndCharges() after the selection has
         * been evaluated.  This is done by SelectionEvaluator.
         *
         * This function is called by SelectionCompiler.
         *
         * Strong exception safety guarantee.
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

    private:
        /*! \brief
         * Additional information about positions.
         *
         * This structure contains information about positions in the
         * selection that is not stored in ::gmx_ana_pos_t.
         *
         * Methods in this class do not throw.
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
        //! Root of the selection evaluation tree.
        SelectionTreeElement   &rootElement_;
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
         * Needed to wrap access to information.
         */
        friend class gmx::Selection;
        /*! \brief
         * Needed for proper access to position information.
         */
        friend class gmx::SelectionPosition;

        GMX_DISALLOW_COPY_AND_ASSIGN(SelectionData);
};

} // namespace internal

/*! \brief
 * Provides access to a single selection.
 *
 * This class provides a public interface for accessing selection information.
 * General information about the selection can be accessed with methods name(),
 * selectionText(), isDynamic(), and type().  The first three can be accessed
 * any time after the selection has been parsed, and type() can be accessed
 * after the selection has been compiled.
 *
 * Each selection is made of a set of positions.  Each position has associated
 * coordinates, and possibly velocities and forces if they have been requested
 * and are available.  It also has a set of atoms associated with it; typically
 * the coordinates are the center-of-mass or center-of-geometry coordinates for
 * that set of atoms.  To access the number of positions in the selection, use
 * posCount().  To access individual positions, use position().
 * See SelectionPosition for details of how to use individual positions.
 * setOriginalId() can be used to adjust the return value of
 * SelectionPosition::mappedId(); see that method for details.
 *
 * It is also possible to access the list of atoms that make up all the
 * positions directly: atomCount() returns the total number of atoms in the
 * selection and atomIndices() an array of their indices.
 *
 * Both positions and atoms can be accessed after the selection has been
 * compiled.  For dynamic selections, the return values of these methods change
 * after each evaluation to reflect the situation for the current frame.
 * Before any frame has been evaluated, these methods return the maximal set
 * to which the selection can evaluate.
 *
 * There are two possible modes for how positions for dynamic selections are
 * handled.  In the default mode, posCount() can change, and for each frame,
 * only the positions that are selected in that frame can be accessed.  In a
 * masked mode, posCount() remains constant, i.e., the positions are always
 * evaluated for the maximal set, and SelectionPosition::selected() is used to
 * determine whether a position is selected for a frame.  The masked mode can
 * be requested with SelectionOption::dynamicMask().
 *
 * The class also provides methods for printing out information: printInfo()
 * and printDebugInfo().  These are mainly for internal use by Gromacs.
 *
 * This class works like a pointer type: copying and assignment is lightweight,
 * and all copies work interchangeably, accessing the same internal data.
 *
 * Methods in this class do not throw.
 *
 * \see SelectionPosition
 *
 * \inpublicapi
 * \ingroup module_selection
 */
class Selection
{
    public:
        /*! \brief
         * Creates a selection wrapper that has no associated selection.
         *
         * Any attempt to call methods in the object before a selection is
         * assigned results in undefined behavior.
         */
        Selection() : sel_(NULL) {}
        /*! \brief
         * Creates a new selection object.
         *
         * \param  sel  Selection data to wrap.
         *
         * Only for internal use by the selection module.
         */
        explicit Selection(internal::SelectionData *sel) : sel_(sel) {}

        //! Returns the name of the selection.
        const char *name() const  { return data().name_.c_str(); }
        //! Returns the string that was parsed to produce this selection.
        const char *selectionText() const { return data().selectionText(); }
        //! Returns true if the size of the selection (posCount()) is dynamic.
        bool isDynamic() const { return data().isDynamic(); }
        //! Returns the type of positions in the selection.
        e_index_t type() const { return data().rawPositions_.m.type; }

        //! Total number of atoms in the selection.
        int atomCount() const
        {
            return data().rawPositions_.g != NULL ? data().rawPositions_.g->isize : 0;
        }
        //! Returns atom indices of all atoms in the selection.
        ConstArrayRef<int> atomIndices() const
        {
            if (data().rawPositions_.g == NULL)
            {
                return ConstArrayRef<int>();
            }
            return ConstArrayRef<int>(data().rawPositions_.g->isize,
                                      data().rawPositions_.g->index);
        }
        //! Number of positions in the selection.
        int posCount() const { return data().posCount(); }
        //! Access a single position.
        SelectionPosition position(int i) const;
        /*! \brief
         * Sets the ID for the \p i'th position for use with
         * SelectionPosition::mappedId().
         *
         * \param[in] i  Zero-based index
         * \param[in] id Identifier to set.
         *
         * This method is not part of SelectionPosition because that interface
         * only provides access to const data by design.
         *
         * This method can only be called after compilation, before the
         * selection has been evaluated for any frame.
         *
         * \see SelectionPosition::mappedId()
         */
        void setOriginalId(int i, int id) { data().rawPositions_.m.orgid[i] = id; }

        //! Deprecated method for direct access to position data.
        const gmx_ana_pos_t *positions() const { return &data().rawPositions_; }

        //! Returns whether the covered fraction can change between frames.
        bool isCoveredFractionDynamic() const { return data().isCoveredFractionDynamic(); }
        //! Returns the covered fraction for the current frame.
        real coveredFraction() const { return data().coveredFraction_; }
        /*! \brief
         * Initializes information about covered fractions.
         *
         * \param[in] type Type of covered fraction required.
         * \returns   true if the covered fraction can be calculated for the
         *      selection.
         */
        bool initCoveredFraction(e_coverfrac_t type)
        {
            return data().initCoveredFraction(type);
        }

        /*! \brief
         * Prints out one-line description of the selection.
         *
         * \param[in] fp      Where to print the information.
         *
         * The output contains the name of the selection, the number of atoms
         * and the number of positions, and indication of whether the selection
         * is dynamic.
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
        internal::SelectionData &data()
        {
            GMX_ASSERT(sel_ != NULL,
                       "Attempted to access uninitialized selection");
            return *sel_;
        }
        const internal::SelectionData &data() const
        {
            GMX_ASSERT(sel_ != NULL,
                       "Attempted to access uninitialized selection");
            return *sel_;
        }

        /*! \brief
         * Pointer to internal data for the selection.
         *
         * The memory for this object is managed by a SelectionCollection
         * object, and the \ref Selection class simply provides a public
         * interface for accessing the data.
         */
        internal::SelectionData *sel_;

        /*! \brief
         * Needed to access the data to adjust flags.
         */
        friend class SelectionOptionStorage;
};

/*! \brief
 * Provides access to information about a single selected position.
 *
 * Each position has associated coordinates, and possibly velocities and forces
 * if they have been requested and are available.  It also has a set of atoms
 * associated with it; typically the coordinates are the center-of-mass or
 * center-of-geometry coordinates for that set of atoms.  It is possible that
 * there are not atoms associated if the selection has been provided as a fixed
 * position.
 *
 * After the selection has been compiled, but not yet evaluated, the contents
 * of the coordinate, velocity and force vectors are undefined.
 *
 * Default copy constructor and assignment operators are used, and work as
 * intended: the copy references the same position and works identically.
 *
 * Methods in this class do not throw.
 *
 * \see Selection
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
         * Asserts if \p index is out of range.
         *
         * Only for internal use of the library.  To obtain a SelectionPosition
         * object in other code, use Selection::position().
         */
        SelectionPosition(const internal::SelectionData &sel, int index)
            : sel_(&sel), i_(index)
        {
            GMX_ASSERT(index >= 0 && index < sel.posCount(),
                       "Invalid selection position index");
        }

        /*! \brief
         * Returns type of this position.
         *
         * Currently always returns the same as Selection::type().
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
         * Returns force for this position.
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
         * If there are no atoms associated or masses are not available,
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
         * If there are no atoms associated or charges are not available,
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
         * Returns whether this position is selected in the current frame.
         *
         * The return value is equivalent to \c refid() == -1.  Returns always
         * true if SelectionOption::dynamicMask() has not been set.
         *
         * \see refId()
         */
        bool selected() const
        {
            return refId() >= 0;
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
         * If a dynamic selection consists of at most three positions, after
         * compilation refId() will return 0, 1, 2 for them, respectively.
         * If for a particular frame, only the first and the third are present,
         * refId() will return 0, 2.
         * If SelectionOption::dynamicMask() has been set, all three positions
         * can be accessed also for that frame and refId() will return 0, -1,
         * 2.
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
         * If for an array \c id, \c setOriginalId(i, id[i]) has been called
         * for each \c i, then it always holds that
         * \c mappedId()==id[refId()].
         *
         * Selection::setOriginalId() has not been called, the default values
         * are dependent on type():
         *  - ::INDEX_ATOM: atom indices
         *  - ::INDEX_RES:  residue numbers
         *  - ::INDEX_MOL:  molecule numbers
         *  .
         * All the default values are zero-based
         */
        int mappedId() const
        {
            return sel_->rawPositions_.m.mapid[i_];
        }

    private:
        const internal::SelectionData  *sel_;
        int                             i_;
};


inline SelectionPosition
Selection::position(int i) const
{
    return SelectionPosition(data(), i);
}

} // namespace gmx

#endif
