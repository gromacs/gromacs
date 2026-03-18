#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2024- The GROMACS Authors
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
gmxcitations.py - Sphinx extension for inline citation hover tooltips.

Parses references.rst at build time, then post-processes every generated HTML
file to inject data-gmx-citation attributes on links that point to citation
anchors. A small CSS file turns those attributes into hover tooltip bubbles.

No new package dependencies. Works on any static hosting (incl. manual.gromacs.org).
The suggested sphinx-hoverxref extension was considered but requires Read the Docs
hosting; since GROMACS docs are self-hosted, a static HTML post-processing approach is used instead.

Setup
-----
1. Place this file in the same directory as gmxsphinx.py.
2. Add "gmxcitations" to the extensions list in conf.cmakein.py.
3. Rebuild the docs.
"""

import logging
import re
from pathlib import Path
from typing import Dict

from sphinx.application import Sphinx
from sphinx.builders.html import StandaloneHTMLBuilder

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Step 1: Parse references.rst -> {label: plain_text}
# ---------------------------------------------------------------------------

_LABEL_RE = re.compile(r"^\.\. _(?P<label>ref\w+):\s*$")
_SUP_RE = re.compile(r":sup:`[^`]+`\s*")
_RST_ROLE_RE = re.compile(r":\w+:`([^`]*)`")
_RAW_HTML_RE = re.compile(r"\.\. raw:: html.*?(?=\n\S|\Z)", re.DOTALL)


def _strip_rst(text: str) -> str:
    """Convert RST inline markup to plain text suitable for an HTML attribute."""
    text = _SUP_RE.sub("", text)
    text = _RST_ROLE_RE.sub(r"\1", text)
    text = re.sub(r"\*{1,2}([^*]+)\*{1,2}", r"\1", text)
    text = re.sub(r"``([^`]+)``", r"\1", text)
    text = re.sub(r"\\ ", " ", text)
    return " ".join(text.split()).strip()


def _parse_references(rst_path: Path) -> Dict[str, str]:
    """Return {label: plain_text_citation} parsed from references.rst."""
    if not rst_path.exists():
        logger.warning(f"[gmxcitations] references.rst not found at {rst_path}")
        return {}

    content = rst_path.read_text(encoding="utf-8")
    content = _RAW_HTML_RE.sub("", content)

    citation_map: Dict[str, str] = {}
    lines = content.splitlines()
    i = 0
    while i < len(lines):
        m = _LABEL_RE.match(lines[i])
        if m:
            label = m.group("label")
            i += 1
            while i < len(lines) and not lines[i].strip():
                i += 1
            text_lines = []
            while i < len(lines):
                stripped = lines[i].strip()
                if not stripped or stripped.startswith(".."):
                    break
                text_lines.append(stripped)
                i += 1
            if text_lines:
                citation_map[label] = _strip_rst(" ".join(text_lines))
        else:
            i += 1

    return citation_map


# ---------------------------------------------------------------------------
# Step 2: Store citation map on the Sphinx environment
# ---------------------------------------------------------------------------


def on_env_get_outdated(app, env, added, changed, removed):
    refs_rst = Path(app.srcdir) / "reference-manual" / "references.rst"
    env.gmx_citation_map = _parse_references(refs_rst)
    print(f"[gmxcitations] loaded {len(env.gmx_citation_map)} citations")
    return []


# ---------------------------------------------------------------------------
# Step 3: Post-process HTML files to inject data-gmx-citation attributes
# ---------------------------------------------------------------------------


def _escape_attr(text: str) -> str:
    return (
        text.replace("&", "&amp;")
        .replace('"', "&quot;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
    )


_HREF_FRAG_RE = re.compile(
    r'(<a\b[^>]*\bhref=["\'][^"\']*#)(ref\w+)(["\'][^>]*>)',
    re.IGNORECASE,
)


def _patch_html(content: str, citation_map: Dict[str, str]) -> str:
    """Add data-gmx-citation attribute to citation anchor links."""
    lower_map = {k.lower(): v for k, v in citation_map.items()}

    def replacer(m):
        fragment = m.group(2).lower()
        if fragment in lower_map:
            escaped = _escape_attr(lower_map[fragment])
            tag = m.group(0)
            return tag[:-1] + f' data-gmx-citation="{escaped}">'
        return m.group(0)

    return _HREF_FRAG_RE.sub(replacer, content)


def on_build_finished(app, exception):
    if exception or not isinstance(app.builder, StandaloneHTMLBuilder):
        return

    citation_map: Dict[str, str] = getattr(app.env, "gmx_citation_map", {})
    if not citation_map:
        logger.warning("[gmxcitations] no citations loaded, skipping HTML patching")
        return

    outdir = Path(app.outdir)

    (outdir / "_static" / "gmxcitations.css").write_text(_TOOLTIP_CSS, encoding="utf-8")

    patched = 0
    for html_file in outdir.rglob("*.html"):
        content = html_file.read_text(encoding="utf-8")
        new_content = _patch_html(content, citation_map)
        if new_content != content:
            html_file.write_text(new_content, encoding="utf-8")
            patched += 1

    print(f"[gmxcitations] patched {patched} HTML files with citation tooltips")


# ---------------------------------------------------------------------------
# Tooltip CSS
# ---------------------------------------------------------------------------

_TOOLTIP_CSS = """\
/* gmxcitations: citation hover tooltips */
a[data-gmx-citation] {
    position: relative;
    cursor: help;
}
a[data-gmx-citation]::after {
    content: attr(data-gmx-citation);
    position: absolute;
    bottom: calc(100% + 6px);
    left: 50%;
    transform: translateX(-50%);
    width: max-content;
    max-width: min(480px, 90vw);
    padding: 6px 10px;
    border-radius: 4px;
    font-size: 0.78rem;
    font-style: normal;
    line-height: 1.45;
    white-space: normal;
    word-break: break-word;
    pointer-events: none;
    opacity: 0;
    transition: opacity 0.15s ease;
    z-index: 1000;
    background: var(--color-tooltip-background, #2b2b2b);
    color: var(--color-tooltip-foreground, #f0f0f0);
    box-shadow: 0 2px 8px rgba(0,0,0,0.35);
}
a[data-gmx-citation]:hover::after,
a[data-gmx-citation]:focus::after {
    opacity: 1;
}
"""


# ---------------------------------------------------------------------------
# Extension entry point
# ---------------------------------------------------------------------------


def setup(app: Sphinx):
    app.connect("env-get-outdated", on_env_get_outdated)
    app.connect("build-finished", on_build_finished)
    app.add_css_file("gmxcitations.css")

    return {
        "version": "1.0",
        "parallel_read_safe": True,
        "parallel_write_safe": True,
    }
