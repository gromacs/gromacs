/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016, by the GROMACS development team, led by
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
#include "gmxpre.h"

#include "gromacs/utility/keyvaluetreetransform.h"

#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/utility/keyvaluetree.h"
#include "gromacs/utility/keyvaluetreebuilder.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/strconvert.h"

#include "testutils/refdata.h"

namespace
{

TEST(TreeValueTransformTest, AssignsFromTree)
{
    gmx::KeyValueTreeBuilder     builder;
    builder.rootObject().addValue<std::string>("a", "1 2");
    gmx::KeyValueTreeObject      from = builder.build();

    gmx::KeyValueTreeTransformer transform;
    transform.rules()->addRule()
        .from<std::string>("/a").toObject("/foo").transformWith(
            [] (gmx::KeyValueTreeObjectBuilder *builder, const std::string &value)
            {
                std::vector<std::string> values = gmx::splitString(value);
                builder->addValue<int>("a", gmx::fromString<int>(values[0]));
                builder->addValue<int>("b", gmx::fromString<int>(values[1]));
            });

    gmx::KeyValueTreeObject         result = transform.transform(from);

    gmx::test::TestReferenceData    data;
    gmx::test::TestReferenceChecker checker(data.rootChecker());
    checker.checkKeyValueTreeObject(from, "Input");
    checker.checkKeyValueTreeObject(result, "Tree");
}

} // namespace
