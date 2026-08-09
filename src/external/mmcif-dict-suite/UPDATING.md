# Vendored RCSB mmCIF Dictionary Suite (trimmed for GROMACS)

- **Upstream:** https://github.com/rcsb/cpp-dict-pack (umbrella; its `modules/*`
  are git submodules of the per-module repos below).
- **License:** see `LICENSE` (upstream RCSB license).

**Pinned upstream versions** (each module's `master` HEAD, which is also the commit
`cpp-dict-pack` currently pins as a submodule). Upstream distributes numbered
releases such as `v2.500`, but this copy is vendored per module, so the commits
below — not a release number — identify what is here:

| module              | commit    | date       |
| ------------------- | --------- | ---------- |
| `cpp-common`        | `0048354` | 2023-08-16 |
| `cpp-tables`        | `d9d336c` | 2021-09-01 |
| `cc-regex`          | `28b88be` | 2011-09-28 |
| `cpp-cif-file`      | `d01a5c7` | 2024-06-20 |
| `cpp-cif-file-util` | `607d488` | 2022-05-17 |
| `cpp-cif-parser`    | `874d19b` | 2018-01-10 |

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

## Local patches

These sources are **not pristine upstream**. Both changes below are submitted
upstream, so a refresh that already contains them must drop the local copy
rather than re-apply it.

| change | where | upstream |
| ------ | ----- | -------- |
| Replaced `std::bind2nd` / `std::not1` and the `std::unary_function` / `std::binary_function` bases, removed from the standard library in C++17, with type aliases and lambdas. Needed because `cifio.cpp` includes `GenString.h` transitively and is compiled as C++17+. | `cpp-common`: `GenString.h`, `GenString.C`, `GenCont.C` | rcsb/cpp-common#1 |
| Annotated 17 harmless `deadcode.DeadStores` with `[[clang::suppress]]` for the Clang-analyzer CI job (see `CMakeLists.txt` for why a build flag cannot do this). | `git grep '\[\[clang::suppress\]\]' modules/` | rcsb/cpp-tables#1, rcsb/cpp-cif-file#2 |

## Updating to newer upstream versions

Each module is vendored at a specific upstream commit (see the table above). To
refresh:

1. For each module `<m>`, get the desired upstream commit from
   `https://github.com/rcsb/<m>` (or the commit `cpp-dict-pack` pins for `<m>`),
   and overwrite the vendored `modules/<m>/{include,src}` files with it (keep only
   the files already vendored — the trimmed set — unless a new one is required to
   compile). Update the commit table above.
2. Re-check the local patches listed above: build with libc++ at C++17/20 and run
   `scan-build`; drop whatever the new sources no longer need, and (re)add only
   what they still require.
3. Regenerate the committed parser sources **only if** the flex/bison grammar
   (`modules/cpp-cif-parser/src/*.l`, `*.y`) changed:

   ```bash
   ./regenerate-parsers.sh
   ```

   (requires `flex` and `bison` — needed only for this step, never for a normal
   GROMACS build).
4. Update the pinned version above, and reconcile `CMakeLists.txt` source lists
   if upstream added/renamed files (compare against upstream
   `cmake/commonparser.cmake`).
5. Rebuild GROMACS and run the fileio tests.
