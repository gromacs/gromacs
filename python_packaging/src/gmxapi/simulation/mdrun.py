#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2019, by the GROMACS development team, led by
# Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
# and including many others, as listed in the AUTHORS file in the
# top-level source directory and at http://www.gromacs.org.
#
# GROMACS is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License
# as published by the Free Software Foundation; either version 2.1
# of the License, or (at your option) any later version.
#
# GROMACS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with GROMACS; if not, see
# http://www.gnu.org/licenses, or write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
#
# If you want to redistribute modifications to GROMACS, please
# consider that scientific software is very special. Version
# control is crucial - bugs must be traceable. We will be happy to
# consider code for inclusion in the official distribution, but
# derived work must not be called official GROMACS. Details are found
# in the README & COPYING files - if they are missing, get the
# official version at http://www.gromacs.org.
#
# To help us fund GROMACS development, we humbly ask that you cite
# the research papers on the package. Check out http://www.gromacs.org.

"""mdrun operation module

The mdrun operation (in its first draft) conforms to the user-level API, but does
not use the Python Context resource manager. It uses either the legacy 0.0.7
Context or its own Context, also implemented in this module.
"""

import typing

from gmxapi import exceptions
from gmxapi.operation import AbstractOperation
from gmxapi.operation import function_wrapper
# The following imports are not marked as public API.
# TODO: Resolve public API and exposure.
from . import workflow


def mdrun_dispatcher(context, *, input, label: str = None, **kwargs) -> typing.Type['AbstractOperation']:
    """Dispatch to an appropriate director based on the context and input.

    Runs appropriate director code to set up an operation, returning a handle to
    a simulation operation node.

    Arguments:
        context: Execution context in which to attempt to instantiate an operation from input
        input: operation input. Argument type is context dependent.
        label: human readable tag for the work node to be created.

    Additional key word arguments are passed on to the dispatched factory.

    """
    # TODO: This is a legacy import that should be updated.
    from .context import get_context
    if label is not None:
        raise exceptions.NotImplementedError('sorry... no labels yet')
    try:
        legacy_context = get_context(work=input)
    except Exception:
        legacy_context = None

    def run_session():
        with legacy_context as session:
            session.run()
        return True

    if context is not None and context is legacy_context:
        helper = function_wrapper(output={'data': bool})(run_session)
        return helper(**kwargs)
    else:
        raise exceptions.ValueError('Could not dispatch MD input {} with context {}'.format(input, legacy_context))


def mdrun(input=None) -> typing.Type['AbstractOperation']:
    """MD simulation operation.

    Arguments:
        input : valid simulation input

    Returns:
        runnable operation to perform the specified simulation

    The returned object has a `run()` method to launch the simulation.
    Otherwise, this operation does not yet support the gmxapi data flow model.

    `input` may be a TPR file name.

    Note:
        New function names will be appearing to handle tasks that are separate

        "simulate" is plausibly a dispatcher or base class for various tasks
        dispatched by mdrun. Specific work factories are likely "minimize,"
        "test_particle_insertion," "legacy_simulation" (do_md), or "simulation"
        composition (which may be leap-frog, vv, and other algorithms)
    """
    # If input looks like classic filename input, use the gmxapi 0.0.7
    # implementation.
    from .context import get_context
    try:
        work = workflow.from_tpr(input)
        assert work is not None
        work_context = get_context(work)
        assert work_context is not None
    except exceptions.UsageError:
        # Input was not gmxapi 0.0.7 input.
        work = None
        work_context = None

    # TODO: Inspect input to determine context.
    kwargs = {}
    handle = mdrun_dispatcher(context=work_context, input=work, label=None, **kwargs)

    return handle
