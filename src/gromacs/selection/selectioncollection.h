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
/*! \libinternal \file
 * \brief
 * Declares gmx::SelectionCollection.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \inlibraryapi
 * \ingroup module_selection
 */
#ifndef GMX_SELECTION_SELECTIONCOLLECTION_H
#define GMX_SELECTION_SELECTIONCOLLECTION_H

#include <string>
#include <vector>

#include <typedefs.h>

struct gmx_ana_indexgrps_t;
struct gmx_ana_poscalc_coll_t;

namespace gmx
{

class Options;
class Selection;

class SelectionCollection
{
    public:
        explicit SelectionCollection(gmx_ana_poscalc_coll_t *pcc);
        ~SelectionCollection();

        int init();
        Options *initOptions();

        static int create(SelectionCollection **scp, gmx_ana_poscalc_coll_t *pcc);

        void setReferencePosType(const char *type);
        void setOutputPosType(const char *type);
        void setMaskOnly(bool bMaskOnly);
        void setVelocityOutput(bool bVelOut);
        void setForceOutput(bool bForceOut);
        void setDebugLevel(int debuglevel);

        bool requiresTopology() const;
        int setTopology(t_topology *top, int natoms);
        void setIndexGroups(gmx_ana_indexgrps_t *grps);
        int parseFromStdin(int count, bool bInteractive,
                           std::vector<Selection *> *output);
        int parseFromFile(const std::string &filename,
                          std::vector<Selection *> *output);
        int parseFromString(const std::string &str,
                            std::vector<Selection *> *output);
        int compile();
        int evaluate(t_trxframe *fr, t_pbc *pbc);
        int evaluateFinal(int nframes);

        void printTree(FILE *fp, bool bValues) const;
        void printXvgrInfo(FILE *fp, output_env_t oenv) const;

        class Impl;
        Impl                   *_impl;

    private:
        // Disallow copy and assign.
        SelectionCollection(const SelectionCollection &);
        void operator =(const SelectionCollection &);
};

} // namespace gmx

#endif
