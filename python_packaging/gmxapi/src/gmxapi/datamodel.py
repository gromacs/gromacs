#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2019- The GROMACS Authors
# and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
# Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
# https://www.gnu.org/licenses, or write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
#
# If you want to redistribute modifications to GROMACS, please
# consider that scientific software is very special. Version
# control is crucial - bugs must be traceable. We will be happy to
# consider code for inclusion in the official distribution, but
# derived work must not be called official GROMACS. Details are found
# in the README & COPYING files - if they are missing, get the
# official version at https://www.gromacs.org.
#
# To help us fund GROMACS development, we humbly ask that you cite
# the research papers on the package. Check out https://www.gromacs.org.

"""gmxapi data types and interfaces.

TBD: https://gitlab.com/gromacs/gromacs/-/issues/2993
"""

__all__ = ["ndarray", "NDArray"]

import collections.abc
import typing

import gmxapi.abc
from gmxapi import exceptions


_T = typing.TypeVar("_T")


class NDArray(gmxapi.abc.NDArray, typing.Generic[_T]):
    """N-Dimensional array type."""

    dtype: typing.Type[_T] = None
    shape: typing.Tuple[int] = ()

    def __init__(self, data: typing.Union[_T, typing.Sequence[_T]] = None):
        if data is None:
            self._values = []
            self.dtype = None
            self.shape = (0,)
        else:
            if hasattr(data, "result") or (
                isinstance(data, collections.abc.Iterable)
                and any([hasattr(item, "result") for item in data])
            ):
                raise exceptions.ValueError(
                    "Make a Future of type NDArray instead of NDArray of type Future, or call result() first."
                )
            if isinstance(data, (str, bytes)):
                data = [data]
                length = 1
            else:
                try:
                    length = len(data)
                except TypeError:
                    # data is a scalar
                    length = 1
                    data = [data]
            self._values = data
            if length > 0:
                self.dtype: typing.Type[_T] = type(data[0])
                self.shape = (length,)
        if len(self._values) > 0 or len(self.shape) > 0:
            assert self.shape[0] == len(self._values)

    def to_list(self):
        return self._values

    def __repr__(self):
        return f"<gmxapi.NDArray: dtype={self.dtype}, shape={self.shape}>"

    def __getitem__(self, i: int) -> _T:
        return self._values[i]

    def __len__(self) -> int:
        return len(self._values)


class ArrayFuture(gmxapi.abc.Future, typing.Generic[_T]):
    """Annotation type for gmxapi array results.

    gmxapi interfaces for structured data need to be updated. Until array data
    has a more concrete scheme, it is useful to have an abstract class for use
    in static type hints to distinguish Futures for arrays of certain element
    types from Futures of scalars of certain types.
    """

    _dtype: typing.Type[_T]

    @property
    def dtype(self) -> typing.Type[_T]:
        return self._dtype

    def result(self) -> NDArray[_T]:
        raise NotImplementedError


def ndarray(data=None, shape=None, dtype=None):
    """Create an NDArray object from the provided iterable.

    Arguments:
        data: object supporting sequence, buffer, or Array Interface protocol

    ..  versionadded:: 0.1
        *shape* and *dtype* parameters

    If ``data`` is provided, ``shape`` and ``dtype`` are optional. If ``data`` is not
    provided, both ``shape`` and ``dtype`` are required.

    If ``data`` is provided and shape is provided, ``data`` must be compatible with
    or convertible to ``shape``. See Broadcast Rules in `datamodel` documentation.

    If ``data`` is provided and ``dtype`` is not provided, data type is inferred
    as the narrowest scalar type necessary to hold any element in ``data``.
    ``dtype``, whether inferred or explicit, must be compatible with all elements
    of ``data``.

    The returned object implements the gmxapi N-dimensional Array Interface.
    """
    if data is None:
        array = NDArray()
    else:
        if isinstance(data, NDArray):
            return data
        # data is not None, but may still be an empty sequence.
        length = 0
        try:
            length = len(data)
        except TypeError:
            # data is a scalar
            length = 1
            data = [data]
        array = NDArray(data)
    return array
