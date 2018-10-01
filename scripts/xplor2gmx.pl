#!/usr/bin/perl -w

use strict;
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

# Turn debugging off (0), or on ( > 0).
my $debug = 1;
# Turn atom name translation on and off
my $trans = 1;

my $pdb   = shift || die "I need the name of the pdb file with correct atom numbers\n";
my $core  = shift || "core.ndx";
my $tbl   = "$ENV{GMXDATA}/top/atom_nom.tbl";

printf "[ distance_restraints ]\n";
printf "; Read an xplor distance restraint file, and output GROMACS distance restraints.\n";
printf "; This also needs a pdb file with correct GROMACS atom numbering.\n";
printf "; I used $pdb for the atom numbers\n";
printf "; This also needs the offset in residues.\n";

# Read the pdb file
# If things go wrong, check whether your pdb file is read correctly.
my $natom = 0;
my $nres  = 0;
my @resname;
my @aname;
my @resnr;
open(PDB,$pdb) || die "Can not open file $pdb\n";
while (my $line = <PDB>) {
    if (index($line,"ATOM") >= 0) {
	my @tmp = split(' ',$line);
	$aname[$natom] = $tmp[2];
	$resnr[$natom] = $tmp[4];
	if (!defined $resname[$tmp[4]]) {
	    $resname[$tmp[4]] = $tmp[3];
	    $nres++;
	}
	$natom++;
    }
}
close PDB;
printf "; I found $natom atoms in the pdb file $pdb\n";
printf "; I found $nres residues in the pdb file $pdb\n";
if ($debug > 1) {
    for (my $i = 0; ($i < $natom); $i ++) {
	printf("; %5d  %5s  %5s  %5d\n",$i+1,$aname[$i],
	       $resname[$resnr[$i]],$resnr[$i]);
    }
}

#
# Read the name translation table.
# Source for this comes from: http://www.bmrb.wisc.edu/ref_info/atom_nom.tbl
# It was adapted slightly for GROMACS use, but only as regards formatting,
# not for content
#
open(TBL,$tbl) || die "Can not open atom-name translation table $tbl\n";
my $ntbl=0;
my @tblxplor;
my @tblgmx;
my @tblres;
while (my $line = <TBL>) {
    my @ttt = split('#',$line);
    my @tmp = split(' ',$ttt[0]);
    if ($#tmp == 3) {
	# New table entry
	$tblres[$ntbl] = $tmp[0];
	$tblxplor[$ntbl] = $tmp[1];
	$tblgmx[$ntbl] = $tmp[3];
	$ntbl++;
    }
}
close TBL;
printf "; Read $ntbl entries from $tbl\n";

my $default = "XXX";

my @templates = (
 [ $default, "HA#", "HA1", "HA2" ],
 [ $default, "HA*", "HA1", "HA2" ],
 [ $default, "HB#",  "HB",         "HB1",	"HB2"	],
 [ $default, "HB*",  "HB",         "HB1",	"HB2"	],
 [ $default, "HG#",  "HG",         "HG1",	"HG2",	"HG11",	"HG12",	"HG13",	"HG21",	"HG22",	"HG23"  ],
 [ $default, "HG*",  "HG",         "HG1",	"HG2",	"HG11",	"HG12",	"HG13",	"HG21",	"HG22",	"HG23"  ],
 [ $default, "HG1#", "HG11",	"HG12",	"HG13"	],
 [ $default, "HG1*", "HG11",	"HG12",	"HG13"	],
 [ $default, "HG2#", "HG21",	"HG22",	"HG23"	],
 [ $default, "HG2*", "HG21",	"HG22",	"HG23"	],
 [ $default, "HD#",  "HD1",	"HD2",	"HD11",	"HD12",	"HD13",	"HD21",	"HD22",	"HD23" ],
 [ $default, "HD*",  "HD1",	"HD2",	"HD11",	"HD12",	"HD13",	"HD21",	"HD22",	"HD23" ],
 [ $default, "HD1#", "HD11",	"HD12"	],
 [ "ILE",    "HD1*", "HD1",	"HD2",   "HD3"	],
 [ $default, "HD1*", "HD11",	"HD12"	],
 [ $default, "HD2#", "HD21",	"HD22"	],
 [ $default, "HD2*", "HD21",	"HD22"	],
 [ $default, "HE#",  "HE",      "HE1",	"HE2"	],
 [ "GLN",    "HE*",  "HE21",	"HE22"	],
 [ $default, "HE*",  "HE",        "HE1",	"HE2"	],
 [ $default, "HE2*", "HE2",       "HE21",	"HE22"	],
 [ $default, "HH#",  "HH11",	"HH12",	"HH21",	"HH22" ],
 [ $default, "HH*",  "HH11",	"HH12",	"HH21",	"HH22" ],
 [ $default, "HZ#",  "HZ",         "HZ1",	"HZ2",	"HZ3"	],
 [ $default, "HZ*",  "HZ",         "HZ1",	"HZ2",	"HZ3"	],
 [ $default, "HN",   "H" ]
);

my $ntranslated = 0;
my $nuntransltd = 0;
sub transl_aname {
    my $resnm  = shift;
    my $atom   = shift;

    if ( $trans == 1 ) {    
	for(my $i = 0; ($i < $ntbl); $i ++) {
	    if ($tblres[$i] eq $resnm) {
		if ($tblxplor[$i] eq $atom) {
		    $ntranslated++;
		    return $tblgmx[$i];
		}
	    }
	}
    }
    $nuntransltd++;
    if ($debug > 1) {
	printf "; No name change for $resnm $atom\n";
    }
    return $atom;
}

sub expand_template {
    my $atom  = shift(@_);
    my $rnum  = shift(@_);
    my $bdone = 0;
    my @atoms;
    my $jj = 0;
   
    die("No residue name for residue $rnum") if (!defined ($resname[$rnum]));
    for (my $tt=0; (($tt <= $#templates) && ($bdone == 0)); $tt++) {
	my $templ = $templates[$tt];
        if (($resname[$rnum] eq $templ->[0] ||
             $default eq $templ->[0]) &&
            ($atom eq $templ->[1])) {
	    for ($jj = 2; ($jj <= $#{$templ}); $jj++) {
		push @atoms, transl_aname($resname[$rnum],$templ->[$jj]);
	    }
	    $bdone  = 1;
	}
    }
    if ($bdone == 0) {
	push @atoms, transl_aname($resname[$rnum],$atom);
    }
    if ($debug > 0) {
	my $natom = $#atoms+1;
	printf("; Found $natom elements for atom $resname[$rnum] %d $atom:",
	       $rnum);
	for my $aa ( @atoms ) {
	    printf " $aa";
	}
	printf "\n";
    }
    return @atoms;
}

if ($debug > 1) {
    printf "; There are $#templates template entries\n";
    for (my $tt=0; ($tt <= $#templates); $tt++) {
	my $templ  = $templates[$tt];
        my $ntempl = $#{$templ};
	printf "; Item $tt ($templates[$tt][0] $templates[$tt][1]) has $ntempl entries\n";
    }
}

# This index file holds numbers of restraints involving core atoms
my @protcore = ( "H", "HN", "HA", "HA1", "HA2", "HB", "HB1", "HB2", "HB3", "HG", "HG1", "HG2", "HG3", "N", "O"  );
open(CORE,">$core") || die "Can not open $core\n";
print CORE "[ core-restraints ]\n";
my $ncore = 0;

my $myindex     = 0;
my $linec       = 0;
my $npair       = 0;
my $nunresolved = 0;
while (my $line = <STDIN>) {
    my @ttt = split('!',$line);
    if (index($ttt[0], "dihedral") >= 0) {
        last;
    }
    elsif ((index($ttt[0],"assign") >= 0) && (index($ttt[0],"!assign") < 0)) {
	my @tmp = split('\(',$ttt[0]);
	# Find first argument
        my $rhaak = undef;
	if (($rhaak  = index($tmp[1],')')) < 0) {
	    printf "No ) in '$tmp[1]'\n";
	}
	my $atres1 = substr($tmp[1],0,$rhaak);
	my @at1 = split('or',$atres1);
	
	# Find second argument
	if (($rhaak  = index($tmp[2],')')) < 0) {
	    printf "No ) in '$tmp[2]'\n";
	}
	my $atres2 = substr($tmp[2],0,$rhaak);
	my @at2 = split('or',$atres2);

	my @expdata = split('\)',$tmp[2]);
	my @dist    = split(' ',$expdata[1]);

	my $bOK   = 0;
	my $bCore = 1;
	foreach my $a1 ( @at1 ) {
	    my @info1  = split(' ',$a1);
	    my $r1     = $info1[1];
	    my @atoms1 = expand_template($info1[4],$r1);

	    foreach my $a2 ( @at2 ) {
		my @info2  = split(' ',$a2);
		my $r2     = $info2[1];
		my @atoms2 = expand_template($info2[4],$r2);

                my $i = undef;
		for ($i = 0; ($i < $natom) && ($resnr[$i] < $r1); $i++) { ; }
		for ( ; ($i < $natom) && ($resnr[$i] == $r1); $i++) { 
		    foreach my $ii ( @atoms1 ) {
			if ($ii eq $aname[$i]) {
			    my $bCoreI = 0;
			    for my $pp ( @protcore ) {
				if ($ii eq $pp) {
				    $bCoreI = 1;
				}
			    }
                            my $j = undef;
			    for ($j = 0; ($j < $natom) && ($resnr[$j] < $r2); $j++) { ; }
			    for ( ; ($j < $natom) && ($resnr[$j] == $r2); $j++) { 
				foreach my $jj ( @atoms2 ) {
				    if ($jj eq $aname[$j]) {
					my $dd     = 0.1*$dist[0];
					my $dminus = 0.1*$dist[1];
					my $dplus  = 0.1*$dist[2];
					my $low    = $dd-$dminus;
					my $up1    = $dd+$dplus;
					my $up2    = $up1+1;
					printf("%5d %5d %5d %5d %5d %7.3f %7.3f %7.3f 1.0; res $r1 $ii - res $r2 $jj\n",$i+1,$j+1,1,$myindex,1,$low,$up1,$up2);
					# Do some checks
					$bOK    = 1;
					my $bCoreJ = 0;
					for my $pp ( @protcore ) {
					    if ($jj eq $pp) {
						$bCoreJ = 1;
					    }
					}
					if (($bCoreI == 0) || ($bCoreJ == 0)) {
					    if ($debug > 0) {
						printf "; No core $ii ($bCoreI) $jj ($bCoreJ)\n";
					    }
					    $bCore = 0;
					}
					$npair++;
				    }
				}
			    }
			}
		    }
		}
	    }
	}
	if ($bCore == 1) {
	    printf CORE "$myindex\n";
	    $ncore++;
	}
	if ($bOK == 0) {
	    printf "; UNRESOLVED: $ttt[0]\n";
	    $nunresolved++;
	}
	else {
	    $myindex++;
	}
    }
    $linec++;
}
close CORE;

printf "; A total of $myindex lines with distance restraints were read from $linec input lines\n";
printf "; From this, $npair actual restraints were generated.\n";
printf "; A total of $nunresolved lines could not be interpreted\n";
printf "; In the process $ntranslated atoms had a name change\n";
printf "; However for $nuntransltd no name change could be found\n";
printf "; usually these are either HN->H renamings or not-existing atoms\n";
printf "; generated by the template expansion (e.g. HG#)\n";
printf "; A total of $ncore restraints were classified as core restraints\n";
