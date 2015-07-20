/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015, by the GROMACS development team, led by
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

#ifndef _SFGRAPH_H_
#define _SFGRAPH_H_

#include <algorithm>
#include <functional>
#include <vector>
#include <tbb/flow_graph.h>

using namespace tbb::flow;

typedef broadcast_node<continue_msg> bnode;
typedef continue_node<continue_msg>  cnode;

/* Simple Flow Graph class */
/* This class allows for easily creating simple TBB flow graphs.
 * It supports any number of independent parallel pipelines but nodes
 * cannot pass data or join. This seems sufficient enough for Gromacs
 * for now and can, of course, be extended as needed.
 */ 
class SFGraph {
public:
    explicit SFGraph(int num_pls) :num_pipelines(num_pls), start(sfgraph), pipelines(num_pls) {}
    ~SFGraph()
    {
         for (int p = 0; p<num_pipelines; p++)
         {
             std::vector<cnode *> pl = pipelines[p];
             for_each(pl.begin(), pl.end(), [](cnode *c){delete c;});
         }  
    }

    template<typename Body>
    void add_node(int plnum, const Body &body)
    {
        if (plnum < 0 || plnum >= num_pipelines)
        {
            return;
        }

        cnode *c = new cnode(sfgraph, body);
        if (pipelines[plnum].empty())
        {
            make_edge(start, *c);
        }
        else
        {
            make_edge(*(pipelines[plnum].back()), *c);
        }
        pipelines[plnum].push_back(c);
    }

    void run()
    {
        start.try_put(continue_msg());
        sfgraph.wait_for_all();
    }

private:
    const int num_pipelines;
    graph sfgraph;
    bnode start;
    std::vector< std::vector<cnode *> > pipelines;
};

#endif /* _SFGRAPH_H_ */
