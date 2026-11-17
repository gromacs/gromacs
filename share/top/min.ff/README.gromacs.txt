Origin and versions of this min.ff directory
============================================

This directory is a copy of the MINFF force field directory maintained at

    https://github.com/mholmboe/minff   (directory min.ff/)

Upstream is the reference. When either side changes, compare against the
versions recorded below to see what the other side still matches.

Versions checked and believed to be equivalent
----------------------------------------------

Checked on 2026-11-03:

    GROMACS    branch x2top-minff (merge request !6364)
    MINFF      mholmboe/minff, branch main, commit
               2bed23e3b24fe8f6a8a1b0b331fe64a9ac448810 (2026-11-02,
               "Label every ion-pair block with its source reference"),
               the latest commit that touches min.ff/ there

The 36 files of this directory were compared with min.ff/ of that commit
(.DS_Store files excluded):

  * 34 files are identical byte for byte.
  * ffnonbonded.itp and ffnonbonded_gminff.itp differ only in line endings:
    upstream has CRLF, this copy has LF. With line endings ignored they are
    identical.
  * No other differences: no GROMACS-specific edits have been made. Please
    keep it that way, and make changes upstream and copy them here.

To repeat the check from a clone of the MINFF repository:

    diff -r --strip-trailing-cr -x .DS_Store <minff>/min.ff <gromacs>/share/top/min.ff

This file, README.gromacs.txt, exists only in the GROMACS copy.

Test data derived from this directory
-------------------------------------

The unit tests of gmx x2top run in
src/gromacs/gmxpreprocess/tests/x2top_data/, which holds a small min.ff
with two files taken from this directory:

  * forcefield.itp: all lines other than comments and blank lines are
    identical to forcefield.itp here.
  * atomname2type.n2t: reduced to the six kaolinite lines (H, Al, Si and
    the Ob, Oh and Op oxygens). Each is identical to a line of
    atomname2type.n2t here, as of the commit above.

If the atom types, charges, masses or bond lengths of these lines change in
this directory, the reduced copy and the reference data of the x2top tests
(tests/refdata/X2topTest_*.xml) need to be updated to match.

Related
-------

The modified x2top that the MINFF repository provides for older GROMACS
releases is in x2top_hack/ in the same repository. The features were
submitted to GROMACS as merge request !6364.
