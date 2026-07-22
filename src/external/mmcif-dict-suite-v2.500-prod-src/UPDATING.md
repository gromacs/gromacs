# Vendored RCSB mmCIF Dictionary Suite (trimmed for GROMACS)

- **Upstream:** https://github.com/rcsb/cpp-dict-pack
- **Pinned version:** mmcif-dict-suite **v2.500** (prod-src distribution)
- **License:** see `LICENSE` (upstream RCSB license).

## What is vendored, and what is not

GROMACS uses only the CIF/PDBx **file** reader and writer (via
`src/gromacs/fileio/cifio.cpp`). To keep the GROMACS build tool-free and
offline, this is a **trimmed** copy:

- Kept: the source of the object modules that upstream combines into
  `mmciflib-all` **except** the dictionary-object-file module — i.e.
  `cpp-common`, `cpp-tables`, `cc-regex`, `cpp-cif-file`, `cpp-cif-file-util`,
  and `cpp-cif-parser` (`modules/`).
- Kept: `generated/` — the flex/bison output (`CifParser.c`, `CifScanner.c`,
  `DICParser.c`, `DICScanner.c` and their headers), committed so the GROMACS
  build needs **no flex, bison, or python**.
- Kept: the flex/bison grammar sources (`modules/cpp-cif-parser/src/*.l`,
  `*.y`) as the source of truth for `generated/`.
- Removed: the upstream CMake/dictionary-download/validation machinery, the
  compiled dictionaries, the command-line tools, tests, and `build/`.
- `CMakeLists.txt` here is **GROMACS-authored** (not upstream) and builds a
  single static `mmcif` library.

## Updating to a new upstream release

1. Download the new `mmcif-dict-suite-vX.XXX-prod-src` distribution.
2. Replace `modules/<kept-module>/{include,src}`, `LICENSE`, `README.md` with
   the new versions (keep the module set listed above).
3. Regenerate the committed parser sources:

   ```bash
   ./regenerate-parsers.sh
   ```

   (requires `flex` and `bison` — needed only for this step, never for a normal
   GROMACS build).
4. Update the pinned version above, and reconcile `CMakeLists.txt` source lists
   if upstream added/renamed files (compare against upstream
   `cmake/commonparser.cmake`).
5. Rebuild GROMACS and run the fileio tests.
