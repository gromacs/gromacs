# Vendored RCSB mmCIF Dictionary Suite (trimmed for GROMACS)

- **Upstream:** https://github.com/rcsb/cpp-dict-pack (umbrella; its `modules/*`
  are git submodules of the per-module repos below).
- **License:** Apache License 2.0, see `LICENSE` (upstream RCSB license).

**Pinned upstream versions** (each module's `master` HEAD, which is also the commit
`cpp-dict-pack` currently pins as a submodule). Upstream distributes numbered
releases such as `v2.500`, but this copy is vendored per module, so the commits
below — not a release number — identify what is here:

| module              | commit    | date       |
| ------------------- | --------- | ---------- |
| `cpp-common`        | `6d05a1f` | 2026-08-28 |
| `cpp-tables`        | `2882714` | 2026-08-28 |
| `cc-regex`          | `28b88be` | 2011-09-28 |
| `cpp-cif-file`      | `a160c96` | 2026-08-28 |
| `cpp-cif-file-util` | `607d488` | 2022-05-17 |
| `cpp-cif-parser`    | `874d19b` | 2018-01-10 |

Those six commits are exactly what `cpp-dict-pack` pins as submodules at its
`b227176`, which is suite version **2.50.11**, so this copy corresponds to that
upstream release. The per-module commits above remain the authoritative record,
since only a subset of the suite is vendored here.

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

Only `CMakeLists.txt`, `UPDATING.md`, `regenerate-parsers.sh`, `generated/` and
`msvc/` are GROMACS-made. `msvc/` holds MSVC-only stand-ins: `unistd.h` for the
POSIX header, which MSVC lacks, and two headers that give the variables shared
between the C parsers and the parser classes C linkage, which upstream omits. Everything else is copied unchanged from upstream: the files
under `modules/<m>` from module `<m>` at the commit in the table above, and
`LICENSE` and `README.md` from `cpp-dict-pack` at `b227176`. If a defect is found
in those files, fix it upstream and refresh, rather than patching in place.

## Updating to newer upstream versions

Each module is vendored at a specific upstream commit (see the table above). To
refresh:

1. For each module `<m>`, get the desired upstream commit from
   `https://github.com/rcsb/<m>` (or the commit `cpp-dict-pack` pins for `<m>`),
   and overwrite the vendored `modules/<m>/{include,src}` files with it (keep only
   the files already vendored — the trimmed set — unless a new one is required to
   compile). Update the commit table above, and copy `LICENSE` and `README.md`
   from the matching `cpp-dict-pack` commit.
2. Verify the new sources need nothing added on top: build with libc++ at
   C++17/20 and run `scan-build`. If something does turn out to be needed, fix it
   upstream first and vendor the result.
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
