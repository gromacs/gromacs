/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016,2017, by the GROMACS development team, led by
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
/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Mohammad Mehdi Ghahremanpour <m.ghahremanpour@hotmail.com>
 */

#ifndef ALEXANDRIA_OPTPARAM_H
#define ALEXANDRIA_OPTPARAM_H

#include <functional>
#include <random>
#include <vector>

#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/oenv.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/units.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/real.h"

namespace alexandria
{

/*! \brief
 * Does Bayesian Monte Carlo (BMC) simulation to find the best paramater set,
 * which has the lowest chi-squared.
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */


class OptParam
{
    private:
        int                     maxiter_;
        const gmx_output_env_t *oenv_;
        const char             *xvgconv_;
        const char             *xvgepot_;
        gmx_bool                bBound_;
        real                    seed_;
        real                    step_;
        real                    temperature_;
        bool                    anneal_;
        
    public:

        OptParam() 
            : 
                maxiter_(100), 
                oenv_(nullptr), 
                xvgconv_(nullptr), 
                xvgepot_(nullptr), 
                bBound_(false), 
                seed_(1993), 
                step_(0.02), 
                temperature_(5), 
                anneal_(true)
            {}

        ~OptParam() {};

        /*! \brief Add command line arguments
         *
         * \param[in] pargs Vector of pargs
         */
        void add_pargs(std::vector<t_pargs> *pargs);

        void Init(const char             *xvgconv,
                  const char             *xvgepot,
                  const gmx_output_env_t *oenv);

        /*! \brief Compute and return the Boltzmann factor
         *
         * \param[in] iter  The iteration number
         * \return The Boltzmann factor
         */
        double computeBeta(int iter);

        //! \brief Return Max # iterations
        int maxIter() const { return maxiter_; }

        //! \brief Return the step
        real step() const { return step_; }

        //! \brief Return whether or not bounds are used for parameters
        bool bounds() const { return bBound_; }

        //! \brief Return xvg file for convergence information
        const char *xvgConv() const { return xvgconv_; }

        //! \brief Return xvg file for epot information
        const char *xvgEpot() const { return xvgepot_; }

        //! \brief Return output environment
        const gmx_output_env_t *oenv() const { return oenv_; }
};

template <class T> class Bayes : public OptParam
{
    using func_t = std::function<T(T v[])>;
    using parm_t = std::vector<T>;

    private:
        func_t  func_;
        parm_t  param_;
        parm_t  psigma_;
        parm_t  pmean_;
        parm_t  lowerBound_;
        parm_t  upperBound_;
        parm_t  bestParam_;
        T      *minEval_;

    public:

        Bayes() {}

        void setFunc(func_t func_, parm_t param_, parm_t lowerBound_, parm_t upperBound_, T *minEval_);

        void setParam(parm_t param);

        void setLowerBound(parm_t lowerBound);

        void setUpperBound(parm_t upperBound);

        /*! \brief
         * Change parameter j based on a random unmber
         * obtained from a uniform distribution.
         */
        void changeParam(int j, real rand);

        /*! \brief
         * Returns the vecor of best found value for each parameter.
         */
        void getBestParam(parm_t &bestParam);

        /*! \brief
         * Returns the vecor of mean value calculated for each parameter.
         */
        void getPmean(parm_t &pmean);

        /*! \brief
         * Returns the vecor of standard deviation calculated for each parameter.
         */
        void getPsigma(parm_t &psigma);

        /*! \brief
         * Run the Bayesian Monte carlo (BMC) simulation
         *
         */
        void simulate();

        ~Bayes() {};
};

template <class T>
void Bayes<T>::setFunc(func_t func,
                       parm_t param,
                       parm_t lowerBound,
                       parm_t upperBound,
                       T     *minEval)
{
    func_       = func;
    param_      = param;
    lowerBound_ = lowerBound;
    upperBound_ = upperBound;
    bestParam_  = param;
    minEval_    = minEval;
}

template <class T>
void Bayes<T>::setParam(parm_t param)
{
    param_ = param;
}

template <class T>
void Bayes<T>::setLowerBound(parm_t lowerBound)
{
    lowerBound_ = lowerBound;
}

template <class T>
void Bayes<T>::setUpperBound(parm_t upperBound)
{
    upperBound_ = upperBound;
}

template <class T>
void Bayes<T>::getBestParam(parm_t &bestParam)
{
    bestParam = bestParam_;
}

template <class T>
void Bayes<T>::getPmean(parm_t &pmean)
{
    pmean = pmean_;
}

template <class T>
void Bayes<T>::getPsigma(parm_t &psigma)
{
    psigma = psigma_;
}

template <class T>
void Bayes<T>::changeParam(int j, real rand)
{
    real delta = (2*rand-1)*step()*fabs(param_[j]);
    param_[j] += delta;
    if (bounds())
    {
        if (param_[j] < lowerBound_[j])
        {
            param_[j] = lowerBound_[j];
        }
        else if (param_[j] > upperBound_[j])
        {
            param_[j] = upperBound_[j];
        }
    }
}

template <class T>
void Bayes<T>::simulate()
{
    T                                storeParam;
    int                              nsum            = 0;
    int                              nParam          = 0; 
    double                           currEval        = 0;
    double                           prevEval        = 0;
    double                           deltaEval       = 0;
    double                           randProbability = 0;
    double                           mcProbability   = 0;    
    parm_t                           sum, sum_of_sq;
    
    FILE                            *fpc             = nullptr;
    FILE                            *fpe             = nullptr;
    
    std::random_device               rd;
    std::mt19937                     gen(rd());
    std::uniform_real_distribution<> uniform(0, 1);

    if (nullptr != xvgConv())
    {
        fpc = xvgropen(xvgConv(), "Parameter convergence", "iteration", "", oenv());
    }
    if (nullptr != xvgEpot())
    {
        fpe = xvgropen(xvgEpot(), "Parameter energy", "iteration", "kT", oenv());
    }

    nParam = param_.size();
    sum.resize(nParam, 0);
    sum_of_sq.resize(nParam, 0);
    pmean_.resize(nParam, 0);
    psigma_.resize(nParam, 0);
    
    prevEval  = func_(param_.data());
    *minEval_ = prevEval;
    for (int iter = 0; iter < nParam*maxIter(); iter++)
    {
        double beta = computeBeta(iter/nParam);
        int       j = static_cast<int>(std::round((1+uniform(gen))*nParam)) % nParam; // Pick random parameter to change
        
        storeParam = param_[j];
        changeParam(j, uniform(gen));
        currEval        = func_(param_.data());
        deltaEval       = currEval-prevEval;
        randProbability = uniform(gen);
        mcProbability   = exp(-beta*deltaEval);
        
        if ((deltaEval < 0) || (mcProbability > randProbability))
        {
            double xiter = (1.0*iter)/nParam;
            if (nullptr != fpc)
            {
                fprintf(fpc, "%8f", xiter);
                for (auto value : param_)
                {
                    fprintf(fpc, "  %10g", value);
                }
                fprintf(fpc, "\n");
                fflush(fpc);
            }
            if (nullptr != fpe)
            {
                fprintf(fpe, "%8f  %10g\n", xiter, prevEval);
                fflush(fpe);
            }
            if (currEval < *minEval_)
            {
                bestParam_ = param_;
                *minEval_  = currEval;
            }
            prevEval = currEval;
        }
        else
        {
            param_[j] = storeParam;
        }
        if (iter >= maxIter()/2)
        {
            for (auto k = 0; k < nParam; k++)
            {
                sum[k]       += param_[k];
                sum_of_sq[k] += gmx::square(param_[k]);
            }
            nsum++;
        }
    }
    if (nsum > 0)
    {
        for (auto k = 0; k < nParam; k++)
        {
            pmean_[k]     = (sum[k]/nsum);
            sum_of_sq[k] /= nsum;
            psigma_[k]    = sqrt(sum_of_sq[k]-gmx::square(pmean_[k]));
        }
    }
    if (nullptr != fpc)
    {
        xvgrclose(fpc);
    }
    if (nullptr != fpe)
    {
        xvgrclose(fpe);
    }
}
}

#endif
