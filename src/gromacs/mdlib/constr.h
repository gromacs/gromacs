/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
/*! \libinternal \file
 * \brief Declares interface to constraint code.
 *
 * \author Berk Hess <hess@kth.se>
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_mdlib
 * \inlibraryapi
 */

#ifndef GMX_MDLIB_CONSTR_H
#define GMX_MDLIB_CONSTR_H

#include <cstdint>
#include <cstdio>

#include <memory>
#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/topology/idef.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/real.h"

struct gmx_edsam;
struct gmx_localtop_t;
struct gmx_moltype_t;
struct gmx_mtop_t;
struct gmx_multisim_t;
struct gmx_wallcycle;
struct pull_t;
struct t_commrec;
struct t_ilist;
struct t_inputrec;
struct t_nrnb;
struct t_pbc;
class t_state;
enum class ConstraintAlgorithm : int;
namespace gmx
{
template<typename T>
class ArrayRefWithPadding;
template<typename>
class ListOfLists;
class ObservablesReducerBuilder;

//! Describes supported flavours of constrained updates.
enum class ConstraintVariable : int
{
    Positions,     /* Constrain positions (mass weighted)             */
    Velocities,    /* Constrain velocities (mass weighted)            */
    Derivative,    /* Constrain a derivative (mass weighted),         *
                    * for instance velocity or acceleration,          *
                    * constraint virial can not be calculated.        */
    Deriv_FlexCon, /* As Derivative, but only output flex. con.       */
    Force,         /* Constrain forces (non mass-weighted)            */
    // TODO What does this do? Improve the comment.
    ForceDispl /* Like Force, but free particles will have mass
                * 1 and frozen particles mass 0                   */
};

/*! \libinternal
 * \brief Handles constraints */
class Constraints
{
private:
    /*! \brief Constructor
     *
     * Private to enforce use of makeConstraints() factory
     * function. */
    Constraints(const gmx_mtop_t&          mtop,
                const t_inputrec&          ir,
                pull_t*                    pull_work,
                FILE*                      log,
                const t_commrec*           cr,
                bool                       useUpdateGroups,
                const gmx_multisim_t*      ms,
                t_nrnb*                    nrnb,
                gmx_wallcycle*             wcycle,
                bool                       pbcHandlingRequired,
                ObservablesReducerBuilder* observablesReducerBuilder,
                int                        numConstraints,
                int                        numSettles);

public:
    /*! \brief This member type helps implement a factory
     * function, because its objects can access the private
     * constructor. */
    struct CreationHelper;

    ~Constraints();

    /*! \brief Returns the total number of flexible constraints in the system. */
    int numFlexibleConstraints() const;

    /*! \brief Returns whether the system contains perturbed constraints */
    bool havePerturbedConstraints() const;

    /*! \brief Set up all the local constraints for the domain.
     *
     * \todo Make this a callback that is called automatically
     * once a new domain has been made. */
    void setConstraints(gmx_localtop_t*                     top,
                        int                                 numAtoms,
                        int                                 numHomeAtoms,
                        gmx::ArrayRef<const real>           masses,
                        gmx::ArrayRef<const real>           inverseMasses,
                        bool                                hasMassPerturbedAtoms,
                        real                                lambda,
                        gmx::ArrayRef<const unsigned short> cFREEZE);

    /*! \brief Applies constraints to coordinates.
     *
     * When econq=ConstraintVariable::Positions constrains
     * coordinates xprime using the directions in x, min_proj is
     * not used.
     *
     * When econq=ConstraintVariable::Derivative, calculates the
     * components xprime in the constraint directions and
     * subtracts these components from min_proj.  So when
     * min_proj=xprime, the constraint components are projected
     * out.
     *
     * When econq=ConstraintVariable::Deriv_FlexCon, the same is
     * done as with ConstraintVariable::Derivative, but only the
     * components of the flexible constraints are stored.
     *
     * delta_step is used for determining the constraint reference lengths
     * when lenA != lenB or will the pull code with a pulling rate.
     * step + delta_step is the step at which the final configuration
     * is meant to be; for update delta_step = 1.
     *
     * step_scaling can be used to update coordinates based on the time
     * step multiplied by this factor. Thus, normally 1.0 is passed. The
     * SD1 integrator uses 0.5 in one of its calls, to correct positions
     * for half a step of changed velocities.
     *
     * If v!=NULL also constrain v by adding the constraint corrections / dt.
     *
     * If computeVirial is true, calculate the constraint virial.
     *
     * Return whether the application of constraints succeeded without error.
     *
     * /note x is non-const, because non-local atoms need to be communicated.
     */
    bool apply(bool                      computeRmsd,
               int64_t                   step,
               int                       delta_step,
               real                      step_scaling,
               ArrayRefWithPadding<RVec> x,
               ArrayRefWithPadding<RVec> xprime,
               ArrayRef<RVec>            min_proj,
               const matrix              box,
               real                      lambda,
               real*                     dhdlambda,
               ArrayRefWithPadding<RVec> v,
               bool                      computeVirial,
               tensor                    constraintsVirial,
               ConstraintVariable        econq);
    //! Links the essentialdynamics and constraint code.
    void saveEdsamPointer(gmx_edsam* ed);
    //! Getter for use by domain decomposition.
    ArrayRef<const ListOfLists<int>> atom2constraints_moltype() const;
    //! Getter for use by domain decomposition.
    ArrayRef<const std::vector<int>> atom2settle_moltype() const;

    /*! \brief Return the RMSD of the constraints when available. */
    real rmsd() const;

    /*! \brief Get the total number of constraints.
     *
     * \returns Total number of constraints, including SETTLE and LINCS/SHAKE constraints.
     */
    int numConstraintsTotal();

private:
    //! Implementation type.
    class Impl;
    //! Implementation object.
    std::unique_ptr<Impl> impl_;
};

/*! \brief Generate a fatal error because of too many LINCS/SETTLE warnings. */
[[noreturn]] void too_many_constraint_warnings(ConstraintAlgorithm eConstrAlg, int warncount);

/*! \brief Returns whether constraint with parameter \p iparamsIndex is a flexible constraint */
static inline bool isConstraintFlexible(ArrayRef<const t_iparams> iparams, int iparamsIndex)
{
    return (iparams[iparamsIndex].constr.dA == 0 && iparams[iparamsIndex].constr.dB == 0);
};

/* The at2con t_blocka struct returned by the routines below
 * contains a list of constraints per atom.
 * The F_CONSTRNC constraints in this structure number consecutively
 * after the F_CONSTR constraints.
 */

/*! \brief Tells make_at2con how to treat flexible constraints */
enum class FlexibleConstraintTreatment
{
    Include, //!< Include all flexible constraints
    Exclude  //!< Exclude all flexible constraints
};

/*! \brief Returns the flexible constraint treatment depending on whether the integrator is dynamic */
FlexibleConstraintTreatment flexibleConstraintTreatment(bool haveDynamicsIntegrator);

/*! \brief Returns a ListOfLists object to go from atoms to constraints
 *
 * The object will contain constraint indices with lower indices
 * directly matching the order in F_CONSTR and higher indices matching
 * the order in F_CONSTRNC offset by the number of constraints in F_CONSTR.
 *
 * \param[in]  moltype   The molecule data
 * \param[in]  iparams   Interaction parameters, can be null when
 *                       \p flexibleConstraintTreatment==Include
 * \param[in]  flexibleConstraintTreatment  The flexible constraint treatment,
 *                                          see enum above
 *
 * \returns a ListOfLists object with all constraints for each atom
 */
ListOfLists<int> make_at2con(const gmx_moltype_t&           moltype,
                             gmx::ArrayRef<const t_iparams> iparams,
                             FlexibleConstraintTreatment    flexibleConstraintTreatment);

/*! \brief Returns a ListOfLists object to go from atoms to constraints
 *
 * The object will contain constraint indices with lower indices
 * directly matching the order in F_CONSTR and higher indices matching
 * the order in F_CONSTRNC offset by the number of constraints in F_CONSTR.
 *
 * \param[in]  numAtoms  The number of atoms to construct the list for
 * \param[in]  ilist     Interaction list, size F_NRE
 * \param[in]  iparams   Interaction parameters, can be null when
 *                       \p flexibleConstraintTreatment==Include
 * \param[in]  flexibleConstraintTreatment  The flexible constraint treatment,
 *                                          see enum above
 *
 * \returns a ListOfLists object with all constraints for each atom
 */
ListOfLists<int> make_at2con(int                             numAtoms,
                             ArrayRef<const InteractionList> ilist,
                             ArrayRef<const t_iparams>       iparams,
                             FlexibleConstraintTreatment     flexibleConstraintTreatment);

//! Return the number of flexible constraints in the \c ilist and \c iparams.
int countFlexibleConstraints(ArrayRef<const InteractionList> ilist, ArrayRef<const t_iparams> iparams);

/*! \brief Returns the constraint iatoms for a constraint number con
 * which comes from a list where F_CONSTR and F_CONSTRNC constraints
 * are concatenated. */
inline const int* constr_iatomptr(gmx::ArrayRef<const int> iatom_constr,
                                  gmx::ArrayRef<const int> iatom_constrnc,
                                  int                      con)
{
    if (con * 3 < iatom_constr.ssize())
    {
        return iatom_constr.data() + con * 3;
    }
    else
    {
        return iatom_constrnc.data() + con * 3 - iatom_constr.size();
    }
};

/*! \brief Constrain the initial coordinates and velocities */
void do_constrain_first(FILE*                     log,
                        gmx::Constraints*         constr,
                        const t_inputrec&         inputrec,
                        int                       numAtoms,
                        int                       numHomeAtoms,
                        ArrayRefWithPadding<RVec> x,
                        ArrayRefWithPadding<RVec> v,
                        const matrix              box,
                        real                      lambda);

/*! \brief Constrain the velocities only.
 *
 * The dhdlambda contribution has to be added to the bonded interactions
 *
 * \param[in,out] constr         The constraints object
 * \param[in]     computeRmsd    Tells whether the constraint RMS deviation should be computed
 * \param[in]     step           The integration step index
 * \param[in,out] state          The state, x is read and halo data updated, v is constrained
 * \param[out]    dhdlambda      The dHdlambda contraint contribution is returned in this
 * \param[in]     computeVirial  Whether the constraint virial contribution should be computed
 * \param[out]    constraintsVirial  The constraint virial is returned in this
 */
void constrain_velocities(gmx::Constraints* constr,
                          bool              computeRmsd,
                          int64_t           step,
                          t_state*          state,
                          real*             dhdlambda,
                          bool              computeVirial,
                          tensor            constraintsVirial);

/*! \brief Constrain the coordinates.
 *
 * Constrain the coordinates \p xp using reference coordinates in \p state.
 * When present, the velocities in \p state are also constrained.
 * The dhdlambda contribution has to be added to the bonded interactions.
 *
 * \param[in,out] constr         The constraints object
 * \param[in]     computeRmsd    Tells whether the constraint RMS deviation should be computed
 * \param[in]     step           The integration step index
 * \param[in,out] state          The state, x is read and halo data updated, v is constrained
 * \param[in,out] xp             The coordinates to constrain
 * \param[out]    dhdlambda      The dHdlambda contraint contribution is returned in this
 * \param[in]     computeVirial  Whether the constraint virial contribution should be computed
 * \param[out]    constraintsVirial  The constraint virial is returned in this
 */
void constrain_coordinates(gmx::Constraints*         constr,
                           bool                      computeRmsd,
                           int64_t                   step,
                           t_state*                  state,
                           ArrayRefWithPadding<RVec> xp,
                           real*                     dhdlambda,
                           bool                      computeVirial,
                           tensor                    constraintsVirial);

} // namespace gmx

#endif
