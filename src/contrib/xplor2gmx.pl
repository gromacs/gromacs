#!/usr/bin/perl -w

#
# This script reads an XPLOR input file with distance restraint data
# as sometimes is found in the pdb database (http://www.pdb.org).
# From this input file dihedral restrints should be removed, such that
# only distance restraints are left. The script can handle ambiguous restraints.
# It converts the distance restraints to GROMACS format.
# A typical command line would be
# ./xplor2gmx.pl 33 conf.pdb < restraints.dat > disre.itp
# You can turn off debugging info, but please do verify your output is correct!
#
# David van der Spoel (spoel@gromacs.org), July 2002
#

$debug = 1;

$res0  = shift || die "I need the residue offset\n";
$pdb   = shift || die "I need the name of the pdb file with correct atom numbers\n";

printf "[ distance_restraints ]\n";
printf "; Read an xplor distance restraint file, and output GROMACS distance restraints.\n";
printf "; This also needs a pdb file with correct GROMACS atom numbering.\n";
printf "; This also needs the offset in residues.\n";

$natom = 0;
open(PDB,$pdb) || die "Can not open file $pdb";
while ($line = <PDB>) {
    if (index($line,"ATOM") >= 0) {
	@tmp = split(' ',$line);
	$aname[$natom] = $tmp[2];
	$resnr[$natom] = $tmp[4];
	$natom++;
    }
}
close PDB;
printf "; I found $natom atoms in the pdb file $pdb\n";

@templates = (
 [ "HH#",	"HH11",	"HH12",	"HH21",	"HH22" ],
 [ "HE#",	"HE1",	"HE2" ],
 [ "HD2#",	"HD21",	"HD22" ],
 [ "HD1#", "HD11",	"HD12" ],
 [ "HD#",  "HD1",  "HD2",	"HD11",	"HD12",	"HD13",	"HD21",	"HD22",	"HD23" ],
 [ "HG2#", "HG21",	"HG22",	"HG23" ],
 [ "HG1#",	"HG11",	"HG12",	"HG13" ],
 [ "HG#",  "HG1",	"HG2" ],
 [ "HB#",	"HB1",	"HB2" ],
 [ "HZ#",  "HZ1", "HZ2", "HZ3" ],
 [ "HN",	"H" ],
);

sub expand_template {
    my $atom  = shift(@_);
    my $bdone = 0;
    my @atoms;
    my $jj = 0;
    
    for (my $tt=0; (($tt <= $#templates) && ($bdone == 0)); $tt++) {
	$templ = $templates[$tt];
	if ($atom eq $templ->[0]) {
	    for ($jj = 1; ($jj <= $#{$templ}); $jj++) {
		push @atoms, $templ->[$jj];
	    }
	    $bdone  = 1;
	}
    }
    if ($bdone == 0) {
	push @atoms, $atom;
    }
    if ($debug == 1) {
	my $natom = $#atoms+1;
	printf "; Found $natom elements for atom $atom:";
	for $aa ( @atoms ) {
	    printf " $aa";
	}
	printf "\n";
    }
    return @atoms;
}

if ($debug == 1) {
    printf "; There are $#templates template entries\n";
    for ($tt=0; ($tt <= $#templates); $tt++) {
	$templ  = $templates[$tt];
	$ntempl = $#{$templ};
	printf "; Item $tt ($templates[$tt][0]) has $ntempl entries\n";
    }
}

$myindex = 0;
$linec = 0;
while ($line = <STDIN>) {
    if ((index($line,"assign") >= 0) && (index($line,"!assign") < 0)) {
	@tmp = split('\(',$line);
	# Find first argument
	if (($rhaak  = index($tmp[1],')')) < 0) {
	    printf "No ) in '$tmp[1]'\n";
	}
	$atres1 = substr($tmp[1],0,$rhaak);
	@at1 = split('or',$atres1);
	
	# Find second argument
	if (($rhaak  = index($tmp[2],')')) < 0) {
	    printf "No ) in '$tmp[2]'\n";
	}
	$atres2 = substr($tmp[2],0,$rhaak);
	@at2 = split('or',$atres2);

	@expdata = split('\)',$tmp[2]);
	@dist    = split(' ',$expdata[1]);
		
	foreach $a1 ( @at1 ) {
	    @info1  = split(' ',$a1);
	    $r1     = $info1[1] - $res0;
	    @atoms1 = expand_template($info1[4]);

	    foreach $a2 ( @at2 ) {
		@info2  = split(' ',$a2);
		$r2     = $info2[1] - $res0;
		@atoms2 = expand_template($info2[4]);

		for ($i = 0; ($i < $natom) && ($resnr[$i] < $r1); $i++) { ; }
		for ( ; ($i < $natom) && ($resnr[$i] == $r1); $i++) { 
		    foreach $ii ( @atoms1 ) {
			if ($ii eq $aname[$i]) {
			    for ($j = 0; ($j < $natom) && ($resnr[$j] < $r2); $j++) { ; }
			    for ( ; ($j < $natom) && ($resnr[$j] == $r2); $j++) { 
				foreach $jj ( @atoms2 ) {
				    if ($jj eq $aname[$j]) {
					$dd = 0.1*$dist[0];
					printf("%5d %5d %5d %5d %5d %7.3f %7.3f 0.0 1.0; res $r1 $ii - res $r2 $jj\n",$i,$j,1,$myindex,1,$dd,$dd+0.1);
				    }
				}
			    }
			}
		    }
		}
	    }
	    $myindex++;
	}
    }
    $linec++;
}

print "; A total of $myindex distance restraints were generated from $linec input lines\n";
