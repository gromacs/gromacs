#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2026- The GROMACS Authors
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

"""
gmxstyles.py - Sphinx extension for GROMACS documentation styles.

Writes a small CSS file with global style fixes for the documentation.
"""

from pathlib import Path
from sphinx.application import Sphinx

_STYLES_CSS = """\
/* Note: this filter inverts all article images. Colored images (e.g. plots
   in file-formats.rst) may look incorrect in dark mode. A per-image opt-out
   class can be added in the future if needed. */
[data-theme="dark"] article img {
    filter: invert(1) hue-rotate(180deg);
}

@media (prefers-color-scheme: dark) {
    [data-theme="auto"] article img {
        filter: invert(1) hue-rotate(180deg);
    }
}
"""


def on_build_finished(app: Sphinx, exception) -> None:
    if exception:
        return
    outdir = Path(app.outdir)
    (outdir / "_static").mkdir(parents=True, exist_ok=True)
    (outdir / "_static" / "gmxstyles.css").write_text(_STYLES_CSS, encoding="utf-8")


def setup(app: Sphinx):
    app.connect("build-finished", on_build_finished)
    app.add_css_file("gmxstyles.css")

    return {
        "version": "1.0",
        "parallel_read_safe": True,
        "parallel_write_safe": True,
    }
