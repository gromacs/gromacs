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
#include "gmxpre.h"

#include "gromacs/pbcutil/mshift.h"

#include <cmath>
#include <cstdio>
#include <cstring>

#include <algorithm>
#include <array>
#include <filesystem>
#include <string>
#include <vector>

#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/listoflists.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/strconvert.h"
#include "gromacs/utility/stringutil.h"

/************************************************************
 *
 *      S H I F T   U T I L I T I E S
 *
 ************************************************************/


/************************************************************
 *
 *      G R A P H   G E N E R A T I O N    C O D E
 *
 ************************************************************/

using gmx::ArrayRef;
using gmx::IVec;

// Class for generating the edges of the graph
class EdgesGenerator
{
public:
    EdgesGenerator(const int endAtom) : edges_(endAtom) {}

    // Adds an edge, bi-directional
    void addEdge(int a0, int a1);

    // Returns the edges
    const std::vector<std::vector<int>>& edges() const { return edges_; }

private:
    // The edges stored as list (for each atom starting at \p startAtom_) of lists of atoms
    std::vector<std::vector<int>> edges_;
};

void EdgesGenerator::addEdge(const int a0, const int a1)
{
    const auto& edges = edges_[a0];
    if (std::find(edges.begin(), edges.end(), a1) == edges.end())
    {
        edges_[a0].push_back(a1);
        edges_[a1].push_back(a0);
    }
}

/* Make the actual graph for an ilist, returns whether an edge was added.
 *
 * When a non-null part array is supplied with part indices for each atom,
 * edges are only added when atoms have a different part index.
 */
template<typename T>
static bool mk_igraph(EdgesGenerator* edgesG, int ftype, const T& il, int at_end, ArrayRef<const int> part)
{
    int  i, j, np;
    int  end;
    bool addedEdge = false;

    end = il.size();

    i = 0;
    while (i < end)
    {
        np = interaction_function[ftype].nratoms;

        if (np > 1 && il.iatoms[i + 1] < at_end)
        {
            if (il.iatoms[i + np] >= at_end)
            {
                gmx_fatal(FARGS,
                          "Molecule in topology has atom numbers below and "
                          "above natoms (%d).\n"
                          "You are probably trying to use a trajectory which does "
                          "not match the first %d atoms of the run input file.\n"
                          "You can make a matching run input file with gmx convert-tpr.",
                          at_end,
                          at_end);
            }
            if (ftype == F_SETTLE)
            {
                /* Bond all the atoms in the settle */
                edgesG->addEdge(il.iatoms[i + 1], il.iatoms[i + 2]);
                edgesG->addEdge(il.iatoms[i + 1], il.iatoms[i + 3]);
                addedEdge = true;
            }
            else if (part.empty())
            {
                /* Simply add this bond */
                for (j = 1; j < np; j++)
                {
                    edgesG->addEdge(il.iatoms[i + j], il.iatoms[i + j + 1]);
                }
                addedEdge = true;
            }
            else
            {
                /* Add this bond when it connects two unlinked parts of the graph */
                for (j = 1; j < np; j++)
                {
                    if (part[il.iatoms[i + j]] != part[il.iatoms[i + j + 1]])
                    {
                        edgesG->addEdge(il.iatoms[i + j], il.iatoms[i + j + 1]);
                        addedEdge = true;
                    }
                }
            }
        }

        i += np + 1;
    }

    return addedEdge;
}

[[noreturn]] static void g_error(int line, const char* file)
{
    gmx_fatal(FARGS, "Trying to print nonexistent graph (file %s, line %d)", file, line);
}
#define GCHECK(g)    \
    if ((g) == NULL) \
    g_error(__LINE__, __FILE__)

void p_graph(FILE* log, const char* title, const t_graph* g)
{
    int         i;
    const char* cc[egcolNR] = { "W", "G", "B" };

    GCHECK(g);
    fprintf(log, "graph:  %s\n", title);
    fprintf(log, "nnodes: %d\n", g->numNodes());
    fprintf(log, "nbound: %d\n", g->numConnectedAtoms);
    fprintf(log, "start:  %d\n", g->edgeAtomBegin);
    fprintf(log, "end:    %d\n", g->edgeAtomEnd);
    fprintf(log, " atom shiftx shifty shiftz C nedg    e1    e2 etc.\n");
    for (i = 0; i < int(g->edges.size()); i++)
    {
        if (!g->edges[i].empty())
        {
            fprintf(log,
                    "%5d%7d%7d%7d %1s%5zu",
                    g->edgeAtomBegin + i + 1,
                    g->ishift[g->edgeAtomBegin + i][XX],
                    g->ishift[g->edgeAtomBegin + i][YY],
                    g->ishift[g->edgeAtomBegin + i][ZZ],
                    (!g->edgeColor.empty()) ? cc[g->edgeColor[i]] : " ",
                    g->edges[i].size());
            for (const int edge : g->edges[i])
            {
                fprintf(log, " %5d", edge + 1);
            }
            fprintf(log, "\n");
        }
    }
    fflush(log);
}

/* Converts the vector of vector of edges to ListOfLists
 * and removes leading and trailing uncoupled atoms
 */
static gmx::ListOfLists<int> convertGraph(FILE*                 fplog,
                                          const EdgesGenerator& edgesG,
                                          int*                  firstConnectedAtom,
                                          int*                  numConnectedAtoms)
{
    gmx::ListOfLists<int> edgesLists;

    *firstConnectedAtom        = edgesG.edges().size();
    *numConnectedAtoms         = 0;
    int numEmptyEntriesSkipped = 0;
    int max_nedge              = 0;
    for (const auto& edges : edgesG.edges())
    {
        if (edges.empty())
        {
            numEmptyEntriesSkipped++;
        }
        else
        {
            if (edgesLists.empty())
            {
                /* We ignore empty entries before the first connected entry */
                *firstConnectedAtom = numEmptyEntriesSkipped;
            }
            else
            {
                /* Push any empty entries we skipped */
                for (int i = 0; i < numEmptyEntriesSkipped; i++)
                {
                    edgesLists.pushBack({});
                }
            }
            numEmptyEntriesSkipped = 0;

            edgesLists.pushBack(edges);

            (*numConnectedAtoms)++;

            max_nedge = std::max(max_nedge, int(gmx::ssize(edges)));
        }
    }

    if (fplog)
    {
        fprintf(fplog, "Max number of graph edges per atom is %d\n", max_nedge);
        fprintf(fplog, "Total number of graph edges is %d\n", edgesLists.numElements());
    }

    return edgesLists;
}

static gmx_bool determine_graph_parts(const EdgesGenerator& edgesG, ArrayRef<int> partNr)
{
    /* Initialize the part array with all entries different */
    const int endAtom = edgesG.edges().size();
    for (int at_i = 0; at_i < endAtom; at_i++)
    {
        partNr[at_i] = at_i;
    }

    /* Loop over the graph until the part array is fixed */
    bool haveMultipleParts = false;
    int  numAtomsChanged   = 0;
    do
    {
        haveMultipleParts = false;
        numAtomsChanged   = 0;
        for (gmx::Index at_i = 0; at_i < gmx::ssize(edgesG.edges()); at_i++)
        {
            for (const int at_i2 : edgesG.edges()[at_i])
            {
                /* Set part for both nodes to the minimum */
                if (partNr[at_i2] > partNr[at_i])
                {
                    partNr[at_i2] = partNr[at_i];
                    numAtomsChanged++;
                }
                else if (partNr[at_i2] < partNr[at_i])
                {
                    partNr[at_i] = partNr[at_i2];
                    numAtomsChanged++;
                }
            }
            if (partNr[at_i] != partNr[0])
            {
                haveMultipleParts = true;
            }
        }
        if (debug)
        {
            fprintf(debug,
                    "graph partNr[] numAtomsChanged=%d, bMultiPart=%s\n",
                    numAtomsChanged,
                    gmx::boolToString(haveMultipleParts));
        }
    } while (numAtomsChanged > 0);

    return haveMultipleParts;
}

template<typename T>
static t_graph mk_graph_ilist(FILE* fplog, const T* ilist, int at_end, gmx_bool bShakeOnly, gmx_bool bSettle)
{
    EdgesGenerator edgesG(at_end);

    t_graph::BondedParts parts = t_graph::BondedParts::Single;

    if (at_end > 0)
    {
        if (!bShakeOnly)
        {
            /* First add all the real bonds: they should determine the molecular
             * graph.
             */
            for (int i = 0; (i < F_NRE); i++)
            {
                if (interaction_function[i].flags & IF_CHEMBOND)
                {
                    mk_igraph(&edgesG, i, ilist[i], at_end, {});
                }
            }

            /* Determine of which separated parts the IF_CHEMBOND graph consists.
             * Store the part numbers in the partNr array.
             */
            std::vector<int> partNr(at_end);
            const bool       bMultiPart = determine_graph_parts(edgesG, partNr);

            if (bMultiPart)
            {
                /* Then add all the other interactions in fixed lists,
                 * but only when they connect parts of the graph
                 * that are not connected through IF_CHEMBOND interactions.
                 */
                bool addedEdge = false;
                for (int i = 0; (i < F_NRE); i++)
                {
                    if (!(interaction_function[i].flags & IF_CHEMBOND))
                    {
                        bool addedEdgeForType = mk_igraph(&edgesG, i, ilist[i], at_end, partNr);
                        addedEdge             = (addedEdge || addedEdgeForType);
                    }
                }

                if (addedEdge)
                {
                    parts = t_graph::BondedParts::MultipleConnected;
                }
                else
                {
                    parts = t_graph::BondedParts::MultipleDisconnected;
                }
            }
        }
        else
        {
            /* This is a special thing used in splitter.c to generate shake-blocks */
            mk_igraph(&edgesG, F_CONSTR, ilist[F_CONSTR], at_end, {});
            if (bSettle)
            {
                mk_igraph(&edgesG, F_SETTLE, ilist[F_SETTLE], at_end, {});
            }
        }
    }

    t_graph graph;
    /* The naming is somewhat confusing, but we need g->shiftAtomEnd
     * for shifthing coordinates to a new array (not in place) when
     * some atoms are not connected by the graph, which runs from
     * g->edgeAtomBegin (>= 0) to g->edgeAtomEnd (<= g->shiftAtomEnd).
     */
    graph.shiftAtomEnd  = at_end;
    graph.edgeAtomBegin = 0;
    graph.edgeAtomEnd   = at_end;
    graph.parts         = parts;
    if (at_end > 0)
    {
        /* Convert the vector of vector of edges to ListOfLists */
        graph.edges = convertGraph(fplog, edgesG, &graph.edgeAtomBegin, &graph.numConnectedAtoms);

        graph.edgeAtomEnd = graph.edgeAtomBegin + graph.edges.ssize();
    }

    graph.edgeColor.resize(graph.edges.size());
    graph.ishift.resize(graph.shiftAtomEnd);

    if (gmx_debug_at)
    {
        p_graph(debug, "graph", &graph);
    }

    return graph;
}

t_graph mk_graph_moltype(const gmx_moltype_t& moltype)
{
    return mk_graph_ilist(nullptr, moltype.ilist.data(), moltype.atoms.nr, FALSE, FALSE);
}

t_graph mk_graph(const InteractionDefinitions& idef, const int numAtoms)
{
    return mk_graph_ilist(nullptr, idef.il.data(), numAtoms, false, false);
}

t_graph* mk_graph(FILE* fplog, const InteractionDefinitions& idef, int at_end, gmx_bool bShakeOnly, gmx_bool bSettle)
{
    t_graph* g = new (t_graph);

    *g = mk_graph_ilist(fplog, idef.il.data(), at_end, bShakeOnly, bSettle);

    return g;
}

t_graph* mk_graph(FILE* fplog, const t_idef* idef, int at_end, gmx_bool bShakeOnly, gmx_bool bSettle)
{
    t_graph* g = new (t_graph);

    *g = mk_graph_ilist(fplog, idef->il, at_end, bShakeOnly, bSettle);

    return g;
}

void done_graph(t_graph* g)
{
    delete g;
}

/************************************************************
 *
 *      S H I F T   C A L C U L A T I O N   C O D E
 *
 ************************************************************/

static void mk_1shift_tric(int          npbcdim,
                           const matrix box,
                           const rvec   hbox,
                           const rvec   xi,
                           const rvec   xj,
                           const int*   mi,
                           int*         mj)
{
    /* Calculate periodicity for triclinic box... */
    int  m, d;
    rvec dx;

    rvec_sub(xi, xj, dx);

    mj[ZZ] = 0;
    for (m = npbcdim - 1; (m >= 0); m--)
    {
        /* If dx < hbox, then xj will be reduced by box, so that
         * xi - xj will be bigger
         */
        if (dx[m] < -hbox[m])
        {
            mj[m] = mi[m] - 1;
            for (d = m - 1; d >= 0; d--)
            {
                dx[d] += box[m][d];
            }
        }
        else if (dx[m] >= hbox[m])
        {
            mj[m] = mi[m] + 1;
            for (d = m - 1; d >= 0; d--)
            {
                dx[d] -= box[m][d];
            }
        }
        else
        {
            mj[m] = mi[m];
        }
    }
}

static void mk_1shift(int npbcdim, const rvec hbox, const rvec xi, const rvec xj, const int* mi, int* mj)
{
    /* Calculate periodicity for rectangular box... */
    int  m;
    rvec dx;

    rvec_sub(xi, xj, dx);

    mj[ZZ] = 0;
    for (m = 0; (m < npbcdim); m++)
    {
        /* If dx < hbox, then xj will be reduced by box, so that
         * xi - xj will be bigger
         */
        if (dx[m] < -hbox[m])
        {
            mj[m] = mi[m] - 1;
        }
        else if (dx[m] >= hbox[m])
        {
            mj[m] = mi[m] + 1;
        }
        else
        {
            mj[m] = mi[m];
        }
    }
}

static void mk_1shift_screw(const matrix box, const rvec hbox, const rvec xi, const rvec xj, const int* mi, int* mj)
{
    /* Calculate periodicity for rectangular box... */
    int  signi, m;
    rvec dx;

    if ((mi[XX] > 0 && mi[XX] % 2 == 1) || (mi[XX] < 0 && -mi[XX] % 2 == 1))
    {
        signi = -1;
    }
    else
    {
        signi = 1;
    }

    rvec_sub(xi, xj, dx);

    if (dx[XX] < -hbox[XX])
    {
        mj[XX] = mi[XX] - 1;
    }
    else if (dx[XX] >= hbox[XX])
    {
        mj[XX] = mi[XX] + 1;
    }
    else
    {
        mj[XX] = mi[XX];
    }
    if (mj[XX] != mi[XX])
    {
        /* Rotate */
        dx[YY] = xi[YY] - (box[YY][YY] + box[ZZ][YY] - xj[YY]);
        dx[ZZ] = xi[ZZ] - (box[ZZ][ZZ] - xj[ZZ]);
    }
    for (m = 1; (m < DIM); m++)
    {
        /* The signs are taken such that we can first shift x and rotate
         * and then shift y and z.
         */
        if (dx[m] < -hbox[m])
        {
            mj[m] = mi[m] - signi;
        }
        else if (dx[m] >= hbox[m])
        {
            mj[m] = mi[m] + signi;
        }
        else
        {
            mj[m] = mi[m];
        }
    }
}

static int mk_grey(ArrayRef<egCol> edgeColor,
                   t_graph*        g,
                   int*            AtomI,
                   int             npbcdim,
                   const matrix    box,
                   const rvec      x[],
                   int*            nerror)
{
    int      m, ng, ai, g0;
    rvec     dx, hbox;
    gmx_bool bTriclinic;
    ivec     is_aj;
    t_pbc    pbc;

    for (m = 0; (m < DIM); m++)
    {
        hbox[m] = box[m][m] * 0.5;
    }
    bTriclinic = TRICLINIC(box);

    g0 = g->edgeAtomBegin;
    ng = 0;
    ai = g0 + *AtomI;

    /* Loop over all the bonds */
    for (const int aj : g->edges[ai - g0])
    {
        /* If there is a white one, make it grey and set pbc */
        if (g->useScrewPbc)
        {
            mk_1shift_screw(box, hbox, x[ai], x[aj], g->ishift[ai], is_aj);
        }
        else if (bTriclinic)
        {
            mk_1shift_tric(npbcdim, box, hbox, x[ai], x[aj], g->ishift[ai], is_aj);
        }
        else
        {
            mk_1shift(npbcdim, hbox, x[ai], x[aj], g->ishift[ai], is_aj);
        }

        if (edgeColor[aj - g0] == egcolWhite)
        {
            if (aj - g0 < *AtomI)
            {
                *AtomI = aj - g0;
            }
            edgeColor[aj - g0] = egcolGrey;

            copy_ivec(is_aj, g->ishift[aj]);

            ng++;
        }
        else if ((is_aj[XX] != g->ishift[aj][XX]) || (is_aj[YY] != g->ishift[aj][YY])
                 || (is_aj[ZZ] != g->ishift[aj][ZZ]))
        {
            if (gmx_debug_at)
            {
                set_pbc(&pbc, PbcType::Unset, box);
                pbc_dx(&pbc, x[ai], x[aj], dx);
                fprintf(debug,
                        "mk_grey: shifts for atom %d due to atom %d\n"
                        "are (%d,%d,%d), should be (%d,%d,%d)\n"
                        "dx = (%g,%g,%g)\n",
                        aj + 1,
                        ai + 1,
                        is_aj[XX],
                        is_aj[YY],
                        is_aj[ZZ],
                        g->ishift[aj][XX],
                        g->ishift[aj][YY],
                        g->ishift[aj][ZZ],
                        dx[XX],
                        dx[YY],
                        dx[ZZ]);
            }
            (*nerror)++;
        }
    }
    return ng;
}

/* Return the first node/atom with colour Col starting at fC.
 * return -1 if none found.
 */
static gmx::Index first_colour(const int fC, const egCol Col, const t_graph* g, ArrayRef<const egCol> edgeColor)
{
    for (gmx::Index i = fC; i < gmx::ssize(g->edges); i++)
    {
        if (!g->edges[i].empty() && edgeColor[i] == Col)
        {
            return i;
        }
    }

    return -1;
}

/* Returns the maximum length of the graph edges for coordinates x */
static real maxEdgeLength(const t_graph& g, PbcType pbcType, const matrix box, const rvec x[])
{
    t_pbc pbc;

    set_pbc(&pbc, pbcType, box);

    real maxEdgeLength2 = 0;

    for (int node = 0; node < int(g.edges.size()); node++)
    {
        for (const int nodeJ : g.edges[node])
        {
            rvec dx;
            pbc_dx(&pbc, x[node], x[nodeJ], dx);
            maxEdgeLength2 = std::max(maxEdgeLength2, norm2(dx));
        }
    }

    return std::sqrt(maxEdgeLength2);
}

void mk_mshift(FILE* log, t_graph* g, PbcType pbcType, const matrix box, const rvec x[])
{
    static int            nerror_tot = 0;
    int                   npbcdim;
    int                   ng, i;
    int                   nW, nG; /* Number of White and Grey nodes */
    int gmx_used_in_debug nB;     /* Number of Black nodes */
    int                   fW, fG; /* First of each category */
    int                   nerror = 0;

    g->useScrewPbc = (pbcType == PbcType::Screw);

    if (pbcType == PbcType::XY)
    {
        npbcdim = 2;
    }
    else
    {
        npbcdim = 3;
    }

    GCHECK(g);
    /* This puts everything in the central box, that is does not move it
     * at all. If we return without doing this for a system without bonds
     * (i.e. only settles) all water molecules are moved to the opposite octant
     */
    for (i = 0; i < g->shiftAtomEnd; i++)
    {
        g->ishift[i][XX] = g->ishift[i][YY] = g->ishift[i][ZZ] = 0;
    }

    if (!g->numConnectedAtoms)
    {
        return;
    }

    std::fill(g->edgeColor.begin(), g->edgeColor.end(), egcolWhite);

    nW = g->numConnectedAtoms;
    nG = 0;
    nB = 0;

    fW = 0;

    while (nW > 0)
    {
        GMX_ASSERT(nW + nG + nB == g->numConnectedAtoms, "Graph coloring inconsistency");
        /* Find the first white, this will always be a larger
         * number than before, because no nodes are made white
         * in the loop
         */
        if ((fW = first_colour(fW, egcolWhite, g, g->edgeColor)) == -1)
        {
            gmx_fatal(FARGS, "No WHITE nodes found while nW=%d\n", nW);
        }

        /* Make the first white node grey */
        g->edgeColor[fW] = egcolGrey;
        nG++;
        nW--;

        /* Initial value for the first grey */
        fG = fW;
        while (nG > 0)
        {
            if ((fG = first_colour(fG, egcolGrey, g, g->edgeColor)) == -1)
            {
                gmx_fatal(FARGS, "No GREY nodes found while nG=%d\n", nG);
            }

            /* Make the first grey node black */
            g->edgeColor[fG] = egcolBlack;
            nB++;
            nG--;

            /* Make all the neighbours of this black node grey
             * and set their periodicity
             */
            ng = mk_grey(g->edgeColor, g, &fG, npbcdim, box, x, &nerror);
            /* ng is the number of white nodes made grey */
            nG += ng;
            nW -= ng;
        }
    }
    if (nerror > 0)
    {
        /* We use a threshold of 0.25*boxSize for generating a fatal error
         * to allow removing PBC for systems with periodic molecules.
         *
         * TODO: Consider a better solution for systems with periodic
         *       molecules. Ideally analysis tools should only ask to make
         *       non-periodic molecules whole in a system with periodic mols.
         */
        constexpr real c_relativeDistanceThreshold = 0.25;

        int npbcdim = numPbcDimensions(pbcType);
        GMX_RELEASE_ASSERT(npbcdim > 0, "Expect PBC with graph");
        real minBoxSize = norm(box[XX]);
        for (int d = 1; d < npbcdim; d++)
        {
            minBoxSize = std::min(minBoxSize, norm(box[d]));
        }
        real maxDistance = maxEdgeLength(*g, pbcType, box, x);
        if (maxDistance >= c_relativeDistanceThreshold * minBoxSize)
        {
            std::string mesg = gmx::formatString(
                    "There are inconsistent shifts over periodic boundaries in a molecule type "
                    "consisting of %d atoms. The longest distance involved in such interactions is "
                    "%.3f nm which is %s half the box length.",
                    g->shiftAtomEnd,
                    maxDistance,
                    maxDistance >= 0.5 * minBoxSize ? "above" : "close to");

            switch (g->parts)
            {
                case t_graph::BondedParts::MultipleConnected:
                    /* Ideally we should check if the long distances are
                     * actually between the parts, but that would require
                     * a lot of extra code.
                     */
                    mesg += " This molecule type consists of multiple parts, e.g. monomers, that "
                            "are connected by interactions that are not chemical bonds, e.g. "
                            "restraints. Such systems can not be treated. The only solution is "
                            "increasing the box size.";
                    break;
                default:
                    mesg += " Either you have excessively large distances between atoms in bonded "
                            "interactions or your system is exploding.";
            }
            gmx_fatal(FARGS, "%s", mesg.c_str());
        }

        /* The most likely reason for arriving here is a periodic molecule. */

        nerror_tot++;
        if (nerror_tot <= 100)
        {
            fprintf(stderr, "There were %d inconsistent shifts. Check your topology\n", nerror);
            if (log)
            {
                fprintf(log, "There were %d inconsistent shifts. Check your topology\n", nerror);
            }
        }
        if (nerror_tot == 100)
        {
            fprintf(stderr, "Will stop reporting inconsistent shifts\n");
            if (log)
            {
                fprintf(log, "Will stop reporting inconsistent shifts\n");
            }
        }
    }
}

/************************************************************
 *
 *      A C T U A L   S H I F T   C O D E
 *
 ************************************************************/

void shift_x(const t_graph* g, const matrix box, const rvec x[], rvec x_s[])
{
    int j, tx, ty, tz;

    GCHECK(g);
    const int            g0 = g->edgeAtomBegin;
    const int            g1 = g->edgeAtomEnd;
    ArrayRef<const IVec> is = g->ishift;

    for (j = 0; j < g0; j++)
    {
        copy_rvec(x[j], x_s[j]);
    }

    if (g->useScrewPbc)
    {
        for (j = g0; (j < g1); j++)
        {
            tx = is[j][XX];
            ty = is[j][YY];
            tz = is[j][ZZ];

            if ((tx > 0 && tx % 2 == 1) || (tx < 0 && -tx % 2 == 1))
            {
                x_s[j][XX] = x[j][XX] + tx * box[XX][XX];
                x_s[j][YY] = box[YY][YY] + box[ZZ][YY] - x[j][YY];
                x_s[j][ZZ] = box[ZZ][ZZ] - x[j][ZZ];
            }
            else
            {
                x_s[j][XX] = x[j][XX];
            }
            x_s[j][YY] = x[j][YY] + ty * box[YY][YY] + tz * box[ZZ][YY];
            x_s[j][ZZ] = x[j][ZZ] + tz * box[ZZ][ZZ];
        }
    }
    else if (TRICLINIC(box))
    {
        for (j = g0; (j < g1); j++)
        {
            tx = is[j][XX];
            ty = is[j][YY];
            tz = is[j][ZZ];

            x_s[j][XX] = x[j][XX] + tx * box[XX][XX] + ty * box[YY][XX] + tz * box[ZZ][XX];
            x_s[j][YY] = x[j][YY] + ty * box[YY][YY] + tz * box[ZZ][YY];
            x_s[j][ZZ] = x[j][ZZ] + tz * box[ZZ][ZZ];
        }
    }
    else
    {
        for (j = g0; (j < g1); j++)
        {
            tx = is[j][XX];
            ty = is[j][YY];
            tz = is[j][ZZ];

            x_s[j][XX] = x[j][XX] + tx * box[XX][XX];
            x_s[j][YY] = x[j][YY] + ty * box[YY][YY];
            x_s[j][ZZ] = x[j][ZZ] + tz * box[ZZ][ZZ];
        }
    }

    for (j = g1; j < g->shiftAtomEnd; j++)
    {
        copy_rvec(x[j], x_s[j]);
    }
}

void shift_self(const t_graph& g, const matrix box, rvec x[])
{
    int j, tx, ty, tz;

    GMX_RELEASE_ASSERT(!g.useScrewPbc, "screw pbc not implemented for shift_self");

    const int            g0 = g.edgeAtomBegin;
    const int            g1 = g.edgeAtomEnd;
    ArrayRef<const IVec> is = g.ishift;

#ifdef DEBUG
    fprintf(stderr, "Shifting atoms %d to %d\n", g0, g0 + gn);
#endif
    if (TRICLINIC(box))
    {
        for (j = g0; (j < g1); j++)
        {
            tx = is[j][XX];
            ty = is[j][YY];
            tz = is[j][ZZ];

            x[j][XX] = x[j][XX] + tx * box[XX][XX] + ty * box[YY][XX] + tz * box[ZZ][XX];
            x[j][YY] = x[j][YY] + ty * box[YY][YY] + tz * box[ZZ][YY];
            x[j][ZZ] = x[j][ZZ] + tz * box[ZZ][ZZ];
        }
    }
    else
    {
        for (j = g0; (j < g1); j++)
        {
            tx = is[j][XX];
            ty = is[j][YY];
            tz = is[j][ZZ];

            x[j][XX] = x[j][XX] + tx * box[XX][XX];
            x[j][YY] = x[j][YY] + ty * box[YY][YY];
            x[j][ZZ] = x[j][ZZ] + tz * box[ZZ][ZZ];
        }
    }
}

void shift_self(const t_graph* g, const matrix box, rvec x[])
{
    shift_self(*g, box, x);
}

void unshift_x(const t_graph* g, const matrix box, rvec x[], const rvec x_s[])
{
    int j, tx, ty, tz;

    if (g->useScrewPbc)
    {
        gmx_incons("screw pbc not implemented (yet) for unshift_x");
    }

    const int            g0 = g->edgeAtomBegin;
    const int            g1 = g->edgeAtomEnd;
    ArrayRef<const IVec> is = g->ishift;

    for (j = 0; j < g0; j++)
    {
        copy_rvec(x_s[j], x[j]);
    }

    if (TRICLINIC(box))
    {
        for (j = g0; (j < g1); j++)
        {
            tx = is[j][XX];
            ty = is[j][YY];
            tz = is[j][ZZ];

            x[j][XX] = x_s[j][XX] - tx * box[XX][XX] - ty * box[YY][XX] - tz * box[ZZ][XX];
            x[j][YY] = x_s[j][YY] - ty * box[YY][YY] - tz * box[ZZ][YY];
            x[j][ZZ] = x_s[j][ZZ] - tz * box[ZZ][ZZ];
        }
    }
    else
    {
        for (j = g0; (j < g1); j++)
        {
            tx = is[j][XX];
            ty = is[j][YY];
            tz = is[j][ZZ];

            x[j][XX] = x_s[j][XX] - tx * box[XX][XX];
            x[j][YY] = x_s[j][YY] - ty * box[YY][YY];
            x[j][ZZ] = x_s[j][ZZ] - tz * box[ZZ][ZZ];
        }
    }

    for (j = g1; j < g->shiftAtomEnd; j++)
    {
        copy_rvec(x_s[j], x[j]);
    }
}

void unshift_self(const t_graph* g, const matrix box, rvec x[])
{
    int j, tx, ty, tz;

    if (g->useScrewPbc)
    {
        gmx_incons("screw pbc not implemented for unshift_self");
    }

    const int            g0 = g->edgeAtomBegin;
    const int            g1 = g->edgeAtomEnd;
    ArrayRef<const IVec> is = g->ishift;

    if (TRICLINIC(box))
    {
        for (j = g0; (j < g1); j++)
        {
            tx = is[j][XX];
            ty = is[j][YY];
            tz = is[j][ZZ];

            x[j][XX] = x[j][XX] - tx * box[XX][XX] - ty * box[YY][XX] - tz * box[ZZ][XX];
            x[j][YY] = x[j][YY] - ty * box[YY][YY] - tz * box[ZZ][YY];
            x[j][ZZ] = x[j][ZZ] - tz * box[ZZ][ZZ];
        }
    }
    else
    {
        for (j = g0; (j < g1); j++)
        {
            tx = is[j][XX];
            ty = is[j][YY];
            tz = is[j][ZZ];

            x[j][XX] = x[j][XX] - tx * box[XX][XX];
            x[j][YY] = x[j][YY] - ty * box[YY][YY];
            x[j][ZZ] = x[j][ZZ] - tz * box[ZZ][ZZ];
        }
    }
}
#undef GCHECK
