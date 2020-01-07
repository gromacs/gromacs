/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2020
 *
 * Developers:
 *             Mohammad Mehdi Ghahremanpour,
 *             Paul J. van Maaren,
 *             David van der Spoel (Project leader)
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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor,
 * Boston, MA  02110-1301, USA.
 */

/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#include "optparam.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <functional>
#include <string>
#include <vector>

#include "gromacs/random.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

namespace alexandria
{

void OptParam::add_pargs(std::vector<t_pargs> *pargs)
{
    t_pargs pa[] = {
        { "-maxiter", FALSE, etINT, {&maxiter_},
          "Max number of iterations for optimization" },
        { "-temp",    FALSE, etREAL, {&temperature_},
          "'Temperature' for the Monte Carlo simulation" },
        { "-anneal", FALSE, etBOOL, {&anneal_},
          "Use annealing in Monte Carlo simulation, starting from the second half of the simulation." },
        { "-seed",   FALSE, etINT,  {&seed_},
          "Random number seed. If zero, a seed will be generated." },
        { "-step",  FALSE, etREAL, {&step_},
          //          "Step size in parameter optimization. Is used as a fraction of the starting value, should be less than 10%. At each reinit step the step size is updated." }
          "Step size for the parameter optimization. Is used as fraction of the available range per parameter which depends on the parameter type" }
    };
    for (size_t i = 0; i < asize(pa); i++)
    {
        pargs->push_back(pa[i]);
    }
}

void OptParam::setOutputFiles(const char             *xvgconv,
                              const char             *xvgepot,
                              const gmx_output_env_t *oenv)
{
    xvgconv_     = xvgconv;
    xvgepot_     = xvgepot;
    oenv_        = oenv;
}

double OptParam::computeBeta(int iter)
{
    double temp = temperature_;
    if (anneal_)
    {
        if (iter >= maxiter_)
        {
            temp = 0;
        }
        else
        {
            temp = temperature_*(1.0 - iter/(maxiter_ + 1.0));
        }
    }
    return 1/(BOLTZ*temp);
}

double OptParam::computeBeta(int maxiter, int iter, int ncycle)
{
    double temp = temperature_;
    if (anneal_)
    {
        if (iter >= maxiter_)
        {
            temp = 0;
        }
        else
        {
            temp = (0.5*temperature_)*((exp(-iter/(0.2*(maxiter+1)))) * (1.1 + cos((ncycle*M_PI*iter)/(maxiter+1))));
        }
    }
    return 1/(BOLTZ*temp);
}

void Bayes::setFunc(func_t  func,
                    double *minEval)
{
    func_       = func;
    minEval_    = minEval;
}

void Bayes::addParam(real val,
                     real factor)
{
    GMX_RELEASE_ASSERT(factor > 0, "Scaling factor for bounds should be larger than zero");
    if (factor < 1)
    {
        factor = 1/factor;
    }
    param_.push_back(val);
    prevParam_.push_back(val);
    lowerBound_.push_back(val/factor);
    upperBound_.push_back(val*factor);
}

void Bayes::addParam(real val,
                     real lower,
                     real upper)
{
    param_.push_back(val);
    prevParam_.push_back(val);
    lowerBound_.push_back(lower);
    upperBound_.push_back(upper);
}


void Bayes::addParamName(std::string name)
{
    paramNames_.push_back(name);
}

void Bayes::changeParam(size_t j, real rand)
{
    GMX_RELEASE_ASSERT(j < param_.size(), "Parameter out of range");
    //real delta = (2*rand-1)*step()*fabs(param_[j]);
    real delta = (2*rand-1)*step()*(upperBound_[j]-lowerBound_[j]);
    param_[j] += delta;
    if (boxConstraint())
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

double Bayes::objFunction(const double v[])
{
    auto              np = param_.size();
    std::vector<bool> changed(np, false);
    for (size_t i = 0; i < np; i++)
    {
        if (prevParam_[i] != v[i])
        {
            changed[i] = true;
        }
    }
    toPolData(changed);
    return calcDeviation();
}

void Bayes::MCMC()
{
    double                           storeParam;
    int                              nsum            = 0;
    int                              nParam          = 0; 
    double                           currEval        = 0;
    double                           prevEval        = 0;
    double                           deltaEval       = 0;
    double                           randProbability = 0;
    double                           mcProbability   = 0; 
    double                           halfIter        = maxIter()/2;   
    parm_t                           sum, sum_of_sq;
    
    FILE                            *fpc             = nullptr;
    FILE                            *fpe             = nullptr;
    
    std::random_device               rd;
    std::mt19937                     gen(rd());
    std::uniform_real_distribution<> uniform(0, 1);

    if (nullptr != xvgConv())
    {
        fpc = xvgropen(xvgConv(), "Parameter convergence", "iteration", "", oenv());
        if (!paramNames_.empty())
        {
            std::vector<const char*> paramNames;
            for (const auto &paramName : paramNames_)
            {   
                paramNames.push_back(paramName.c_str());
            }
            xvgr_legend(fpc, paramNames.size(), paramNames.data(), oenv());   
        }
    }   
    if (nullptr != xvgEpot())
    {
        fpe = xvgropen(xvgEpot(), "Parameter energy", "iteration", "\\f{12}c\\S2\\f{4}", oenv());
    }
    nParam = param_.size();
    sum.resize(nParam, 0);
    sum_of_sq.resize(nParam, 0);
    pmean_.resize(nParam, 0);
    psigma_.resize(nParam, 0);
    attemptedMoves_.resize(nParam, 0);
    acceptedMoves_.resize(nParam, 0);

    prevEval  = func_(param_.data());
    *minEval_ = prevEval;
    if (debug)
    {
        fprintf(debug, "Initial chi2 value = %g\n", prevEval);
    }
    for (int iter = 0; iter < nParam*maxIter(); iter++)
    {
        double beta = computeBeta(iter/nParam);
        int       j = static_cast<int>(std::round((1+uniform(gen))*nParam)) % nParam; // Pick random parameter to change
        attemptedMoves_[j] = attemptedMoves_[j] + 1;
        prevParam_ = param_;
        storeParam = param_[j];
        changeParam(j, uniform(gen));
        currEval        = func_(param_.data());
        deltaEval       = currEval-prevEval;
        randProbability = uniform(gen);
        mcProbability   = exp(-beta*deltaEval);
        
        if ((deltaEval < 0) || (mcProbability > randProbability))
        {
            if (currEval < *minEval_)
            {
                bestParam_ = param_;
                *minEval_  = currEval;
                if (debug)
                {
                    fprintf(debug, "New minimum at %g", currEval);
                    for(int k = 0; k < nParam; k++)
                    {
                        fprintf(debug, " %g", bestParam_[k]);
                    }
                    fprintf(debug, "\n");
                }
            }
            prevEval = currEval;
            acceptedMoves_[j] = acceptedMoves_[j] + 1;
        }
        else
        {
            param_[j] = storeParam;
        }
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
        if (iter >= halfIter)
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
            double ps2    = std::max(0.0, sum_of_sq[k]-gmx::square(pmean_[k]));
            psigma_[k]    = sqrt(ps2);
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

void Bayes::DRAM()
{
    double                           storeParam;
    int                              nsum            = 0;
    int                              nParam          = 0; 
    double                           currEval        = 0;
    double                           prevEval        = 0;
    double                           deltaEval       = 0;
    double                           randProbability = 0;
    double                           mcProbability   = 0; 
    double                           halfIter        = maxIter()/2;   
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
        fpe = xvgropen(xvgEpot(), "Parameter energy", "iteration", "\\f{12}c\\S2\\f{4}", oenv());
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
        // Pick random parameter to change
        int       j = static_cast<int>(std::round((1+uniform(gen))*nParam)) % nParam; 
        
        storeParam = param_[j];
        changeParam(j, uniform(gen));
        currEval        = func_(param_.data());
        deltaEval       = currEval-prevEval;
        randProbability = uniform(gen);
        mcProbability   = exp(-beta*deltaEval);
        
        if ((deltaEval < 0) || (mcProbability > randProbability))
        {
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
        if (iter >= halfIter)
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

void Bayes::dumpParam(FILE *fp)
{
    if (nullptr != fp)
    {
        fprintf(fp, "Parameters:");
        for (auto &p : param_)
        {
            fprintf(fp, " %.3f", p);
        }
        fprintf(fp, "\n");
    }
}
    
}
