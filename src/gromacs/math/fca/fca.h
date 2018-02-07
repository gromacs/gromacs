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
/*!
 * \brief
 * FCA of a trajectory follows the spirit of
 * PCA in finding important modes of motion through an orthonormal
 * coordinate transformation and subsequent sorting of basis
 * vectors. It has been developed inside the trajectoryanlysis
 * framework. Orignal code from Oliver Lange, then ported to
 * gromacs4.6 by Christian Blau, and current patch by Berenger Bramas.
 *
 * Fixes #2005
 *
 * Implementation according to: Lange, Oliver F., and
 * Helmut Grubmüller. "Full correlation analysis of conformational
 * protein dynamics." Proteins: Structure, Function, and
 * Bioinformatics 70.4 (2008): 1294-1312..
 */
#ifndef FCA_H
#define FCA_H

#include "eigenvec.h"
#include "entropy2D.h"
#include "entropy_basic_histo.h"
#include "entropy_planes.h"
#include "fca_minimizing.h"
#include "frame_set.h"
#include "task_2d.h"
#include "utils.h"

#endif
