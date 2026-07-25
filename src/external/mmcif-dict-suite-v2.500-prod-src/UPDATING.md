# Vendored RCSB mmCIF Dictionary Suite (trimmed for GROMACS)

- **Upstream:** https://github.com/rcsb/cpp-dict-pack (umbrella; its `modules/*`
  are git submodules of the per-module repos below).
- **License:** see `LICENSE` (upstream RCSB license).

**Pinned upstream versions** (each module's `master` HEAD, which is also the commit
`cpp-dict-pack` currently pins as a submodule). The directory keeps its historical
`…-v2.500-prod-src` name, but the contents are the per-module snapshots below:

| module | commit | date |
| --- | --- | --- |
| `cpp-common` | `0048354` | 2023-08-16 |
| `cpp-tables` | `d9d336c` | 2021-09-01 |
| `cc-regex` | `28b88be` | 2011-09-28 |
| `cpp-cif-file` | `d01a5c7` | 2024-06-20 |
| `cpp-cif-file-util` | `607d488` | 2022-05-17 |
| `cpp-cif-parser` | `874d19b` | 2018-01-10 |

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

## Local workarounds (to be upstreamed to RCSB)

`CMakeLists.txt` carries a small block of **temporary workarounds** for defects
that live in the upstream RCSB sources, not in GROMACS:

1. **Removed C++ binders.** `GenString.C`, `GenString.h` and `GenCont.C` use
   `std::bind2nd` / `std::not1`, deprecated in C++11 and removed from the
   standard library in C++17. libc++ (unlike libstdc++) drops them in C++17+
   mode, so a C++20 build fails to compile the library. We build only the
   `mmcif` target as C++14 as a bypass.
2. **Static-analyzer dead stores.** The Clang static analyzer reports **17**
   harmless `deadcode.DeadStores` in the vendored code (in `cpp-common`
   `Serializer.C`, `cpp-tables` `ISTable.C`/`ITTable.C`/`TableFile.C`, and
   `cpp-cif-file` `DicFile.C`/`CifParentChild.C`), which fails the analyzer CI
   job. Each flagged statement is annotated in place with the
   `[[clang::suppress]]` attribute (`git grep '\[\[clang::suppress\]\]'` under
   `modules/`). This is non-behavioral: the attribute only silences the analyzer
   and is ignored by non-Clang compilers. A build-flag bypass is not possible
   because the analyzer build compiles with g++ (which rejects Clang-only
   analyzer flags) while only the analysis pass is Clang.

Upstream fix PRs are open: rcsb/cpp-common#1 (the C++17 binders — the priority,
since it is what forces the C++14 build), rcsb/cpp-tables#1 and rcsb/cpp-cif-file#2
(dead stores); tracked by the umbrella issue rcsb/cpp-dict-pack#3. These are
stopgaps: **when** a fix merges and this copy is refreshed to include it, delete
the corresponding workaround (the `CMakeLists.txt` C++14 setting and/or the
relevant `[[clang::suppress]]` annotations) and update this section.

Note: an 18th dead store previously suppressed in `cpp-cif-parser`
`CifParserBase.C` is **gone** with this refresh — upstream `master` already reads
the parse result there — so only 17 remain.

## `Parse(FILE*)` overload — intentionally not carried

Our reader (`ParseCifSimple` in `cpp-cif-file-util`) parses via
`CifParser::Parse(fileName, …)`, which reads and parses the whole file into an
in-memory `CifFile`. It does **not** need the `CifParser::Parse(FILE*, string&)`
overload (rcsb/cpp-cif-parser#1, still open upstream), so this refresh drops it.
Keep that PR on the radar: `Parse(FILE*)` is the entry point we would need if we
later add **sequential / streaming reading of very large CIF files** (parsing from
an open handle instead of loading the whole structure). Re-vendor it at that point.

## Updating to newer upstream versions

Each module is vendored at a specific upstream commit (see the table above). To
refresh:

1. For each module `<m>`, get the desired upstream commit from
   `https://github.com/rcsb/<m>` (or the commit `cpp-dict-pack` pins for `<m>`),
   and overwrite the vendored `modules/<m>/{include,src}` files with it (keep only
   the files already vendored — the trimmed set — unless a new one is required to
   compile). Update the commit table above.
2. Re-check the local workarounds: rebuild the `mmcif` target as C++17/20 and run
   `scan-build`; drop any C++14 setting / `[[clang::suppress]]` that the new code
   no longer needs, and (re)add only what it still flags.
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
