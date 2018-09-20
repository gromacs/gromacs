/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
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
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */

#ifndef GMX_APPLIEDFORCES_DENSITYFITTTING_DENSFITDATA_H
#define GMX_APPLIEDFORCES_DENSITYFITTTING_DENSFITDATA_H

#include <vector>
#include <string>

#include "gromacs/utility/real.h"
#include "gromacs/math/griddata/griddata.h"

#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/utility/cstringutil.h"

#include "measures.h"

struct t_blocka;
struct t_inpfile;
struct t_commrec;
struct t_fileio;
struct warninp;
struct gmx_output_env_t;

namespace gmx
{
class ISerializer;


class DensfitData;
class TimeDependentDensfitParameters
{
    public:
        TimeDependentDensfitParameters() = default;

        TimeDependentDensfitParameters(const char timePointString[STRLEN],
                                       const char forceConstantString[STRLEN],
                                       const char sigmaString[STRLEN],
                                       warninp   *wi);


        void serialize(ISerializer *  serializer);

        real currentForceConstant(real time) const;
        real currentSigma(real time) const;

        void print(FILE *fp, int indent) const;

        bool vary() const;

        // const std::vector<real> &timePoints() const;
        // const std::vector<real> &forceConstants() const;

    private:
        std::vector<real> timePoints_;     /* npoints of time values */
        std::vector<real> forceConstants_; /* force constants  k (kJ/mol) */
        std::vector<real> sigmas_;         /* Same for sigma (nm) */
        /*! Linearly interpolate function at x0 between given values x,y.
         *
         * Assumes sorted x_values and y_values corresponding to x_values.
         * If x0 > max(x_values) or x0 < min(x_values) return last or first y_value
         * respectively.
         */
        real interpolateLinearly(real x0, const std::vector<real> &x,
                                 const std::vector<real> &y) const;
        int str_nelem(const char *str, int maxptr, char *ptr[]);
};

struct DensityFittingGroup
{
    void                 serialize(ISerializer * serializer);
    std::vector<int>     ind_;            /* The global atoms numbers                        */
    std::string          name_;
};

/*! \brief Parameters for density fitting that are known before the simulation.
 *
 * Anything that is known about density fitting and is independent of the course
 * of the simulation belongs in this class.
 * \TODO: grompp should be responsible for all pre-calculation and settings in
 * this classs
 */
class DensfitData
{
    public:
        void serialize(ISerializer *  serializer);

        void do_fio(t_fileio *fio, bool bRead);
        void read_densparams(std::vector<t_inpfile> *inp, warninp *wi);

        void printToOut(FILE *fp, const gmx_output_env_t *oenv) const;
        void print(FILE *fp, int indent) const;

        const GridDataFloat3D &referenceMap() const;
        real dist() const;

        void makeGroups(t_blocka *grps, char **gnames);

        const DensityFittingGroup &fittingGroup() const;

        bool keepAndNumberMaps() const;
        bool outputMapThisStep(int64_t step) const;
        bool fitThisStep(int64_t step) const;
        DensityPotential densityPotential() const;
        bool normalizeMapsToUnity() const;

        const TimeDependentDensfitParameters &timeDependent() const;
        int nStFit() const;

    private:
        GridDataFloat3D                 map_ref_;            /* The map to fit to        */
        TimeDependentDensfitParameters  timeDependent_;
        DensityPotential                densityPotential_;
        bool                            normalizeMapsToUnity_;

        real                            dist_;               /* Cutoff distance (in sigmas) for density spreading spread each
                                                                atom only in the range -dist ... +dist */
        int                             nstfit_;             /* Only recalculate V_fit every nstfit time steps  */

        int                             nstout_;             /* Write diagnostic output every nstout time steps */
        int                             nstmapout_;          /* Keep simulated maps in these intervals          */
        bool                            bKeepAndNumberMaps_; /* ... or just keep the last one ... */

        DensityFittingGroup             group_;
};
}

#endif /* end of include guard: \
          GMX_APPLIEDFORCES_DENSITYFITTTING_DENSFITDATA_H */
