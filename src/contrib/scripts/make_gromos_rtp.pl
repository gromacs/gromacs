#!/usr/bin/perl

# usage: make_gromos_rtp.pl mtb43a1.dat > gromos.rtp
# this script tries to make a residue topology database for GROMACS from
# the GROMOS version of this file. It needs ffG96A.atp to fill in the
# atom types for the atom integer codes. It converts until it finds the string
# "#RNMES", which indicates the start of the solvent part
# of mtb43(a,b)1.dat. Solvents are treated differently in GROMACS

# author: Peter Tieleman, 12 March 1999

# read in the atomtypes 
open (TYPES, "gromos.atp") || die "Can't open gromos.atp: $!";
while (<TYPES>) {
    $type++;
    ($gromostype[$type],$rest) = split(' ',$_);
}

$/ = "# RNME"; 
# split input on residue name. The last line of each block of residue data
# has RNME in it. The first line has the name of the residue. 

# write header
print "[ bondedtypes ]\n";
print "; bonds  angles  dihedrals  impropers\n";
print "    2       2          1          2\n\n";

$first = 1;

# loop over the actual residues
while (<>) {
    if ($first) {$first = 0; next;} #skip over the first nonsense residue
    # initialise for a new residue:
    $read_atoms = 0; $j = 0;
    $read_bonds = 0; $k = 0;
    $read_angles = 0; $l = 0;
    $read_dihedrals = 0; $m = 0;
    $read_impropers = 0; $n = 0;
    
    # now $_ has all data for a single residue in it. 
    # split it up in lines
    @residue_data = split('\n',$_);
    $residue_name = @residue_data[1]; 

    # loop over the lines of a residue
    for ($i=0;$i<@residue_data;$i++) {
	# do we have to skip a line?
	if ($skip == 1) {$skip = 0; next;}
	
	# look for the headers of each section, and move to the start
	# of the relevant data.
	
	# if you find ATOM ANM, start reading atoms 
	if ($residue_data[$i] =~ /\#ATOM ANM/) { 
	    $read_atoms = 1;
	}
	# if you find NB, read bonds
	if ($residue_data[$i] =~ /\#\s+NB/) {
	    $read_atoms = 0; $read_bonds = 1; $skip = 1; #skip next line
	}
	# if you find NBA, read angles
	if ($residue_data[$i] =~ /NBA/) {
	    $read_bonds = 0; $read_angles = 1; $skip = 1; #skip next line
	}
	# if you find NIDA, read impropers
	if ($residue_data[$i] =~ /NIDA/) {
	    $read_angles =0; $read_impropers = 1; $skip = 1; #skip next line
	}
	# if you find NDA, read dihedrals
	if ($residue_data[$i] =~ /NDA/) {
	    $read_impropers= 0; $read_dihedrals=1; $skip = 1; #skip next line
	}

	# stop parsing alltogether if we find #RNMES
	if ($residue_data[$i] =~ /\#RNMES/) { last;}
        # skip lines that start with #
	if ( $residue_data[$i] =~ /^\#/) { next; }
        # also skip lines with END on it, and with MTBUILDBLSOLUTE
	if ($residue_data[$i] =~ /END/) { next; }
        if ($residue_data[$i] =~ /BUILD/) {next;}
	
        if (!($read_atoms || $read_bonds || $read_angles || $read_dihedrals
	      || $read_impropers)) { next; } # we're not reading anything yet
	
	if ($read_atoms) {
	    ($atomnr[$j],$atomname[$j],$atomtype[$j],$mass,$charge[$j],
	     $chargegroup[$j],$dummy) = split(' ',$residue_data[$i]);
	    # some lines only have exclusions, check if there is a name
	    if ($atomname[$j] !~ /[A-Z]+/) {$j--;}
	    $j++; 
	}
	if ($read_bonds) {
	    ($bondi[$k],$bondj[$k],$bondtype[$k]) =
		split(' ',$residue_data[$i]);
	    $k++;
	}
	if ($read_angles) {
	    ($anglei[$l],$anglej[$l],$anglek[$l],$angletype[$l]) = 
		split(' ',$residue_data[$i]);
	    $l++;
	}
	
	if ($read_impropers) {
	    ($imp_i[$m],$imp_j[$m],$imp_k[$m],$imp_l[$m],$imp_type[$m]) = 
		split(' ',$residue_data[$i]);
	    $m++;
	}
	
        if ($read_dihedrals) {
	    ($dih_i[$n],$dih_j[$n],$dih_k[$n],$dih_l[$n],$dih_type[$n]) = 
		split(' ',$residue_data[$i]);
	    $n++;
	}
    } 
    $natoms = $j;

    #  print out this residue to the GROMACS file and go to the next residue
    print "[ $residue_name ]\n";
    print " [ atoms ]\n";
    $chargegroup = 0;
    for ($t=0;$t<$j;$t++) {
	print sprintf("%5s %5s %8.3f %5s\n", $atomname[$t], 
		      $gromostype[$atomtype[$t]], $charge[$t], $chargegroup);
	$chargegroup += $chargegroup[$t];
    }
    
    print " [ bonds ]\n";
    for ($t=0;$t<$k;$t++) {
	$bond ="gb_$bondtype[$t]";
	# convert the nrs from GROMOS to atomnames
	if ($bondi[$t] < 0 || $bondj[$t] < 0) { 
	    print STDERR "Skipping specbond $bondi[$t] $bondj[$t]\n";
	    next; 
	}
	if ($bondi[$t] == $natoms+1) {$ati = "+N";} 
	else {$ati = $atomname[$bondi[$t]-1];}
	if ($bondj[$t] == $natoms+1) {$atj = "+N";}
	else {$atj = $atomname[$bondj[$t]-1];}
	# write them out.
	print sprintf("%5s %5s    %-5s\n", $ati, $atj, $bond);
    }
    
    print " [ angles ]\n";
    print ";   ai    aj    ak  gromos type\n";
    for ($t=0;$t<$l;$t++) {

	# convert numbers to names. 
	# i may be in the previous residue
	if ($anglei[$t] == -1) {$ati = "-C";} 
	else {$ati = $atomname[$anglei[$t]-1];}
	# j is always in this residue
	$atj = $atomname[$anglej[$t]-1];
	# k may be in the next residue
	if ($anglek[$t] == $natoms+1) {$atk = "+N";}
	else {$atk = $atomname[$anglek[$t]-1];}

	# print them out
	$angle = "ga_$angletype[$t]";
	print sprintf("%5s %5s %5s    %-5s\n", $ati, $atj, $atk, $angle);
    }
    
    print " [ impropers ]\n";
    print ";   ai    aj    ak    al  gromos type\n";
    for ($t=0;$t<$m;$t++) {

	# convert numbers to names. 
	#  i may be in the next residue
	if ($imp_i[$t] == $natoms+1) {$ati = "+N";} 
	else {$ati = $atomname[$imp_i[$t]-1];}
	# j may be in the next or in the previous residue
	if ($imp_j[$t] == -1) {$atj = "-C";} 
	elsif ($imp_j[$t] == $natoms+1) {$atj = "+N";}
	else {$atj = $atomname[$imp_j[$t]-1];}
	# k may be in the next residue or in the previous residue
	if ($imp_k[$t] == -1) {$atk = "-C";} 
	elsif ($imp_k[$t] == $natoms+1) {$atk = "+N";}
	else {$atk = $atomname[$imp_k[$t]-1];}
	# l may be in the next residue
	if ($imp_l[$t] == $natoms+1) {$atl = "+N";} 
	else {$atl = $atomname[$imp_l[$t]-1];}

	$improper = "gi_$imp_type[$t]";
	print sprintf("%5s %5s %5s %5s    %-5s\n", $ati, $atj, $atk, $atl, 
		      $improper);
    }

    print " [ dihedrals ]\n";
    print ";   ai    aj    ak    al  gromos type\n";
    for ($t=0;$t<$n;$t++) {
	# convert numbers to names. 
	#  i may be in the previous residue
	if ($dih_i[$t] == -1) {$ati = "-C";} 
	elsif ($dih_i[$t] == -2) {$ati = "-CA";}
	else {$ati = $atomname[$dih_i[$t]-1];}
	# j may be in the previous residue
	if ($dih_j[$t] == -1) {$atj = "-C";} 
	else {$atj = $atomname[$dih_j[$t]-1];}
	# k must be in this residue
	$atk = $atomname[$dih_k[$t]-1];
	# l may be in the next residue
	if ($dih_l[$t] == $natoms+1) {$atl = "+N";} 
	else {$atl = $atomname[$dih_l[$t]-1];}

	$dihedral = "gd_$dih_type[$t]";
	print sprintf("%5s %5s %5s %5s    %-5s\n", $ati, $atj, $atk, $atl,
		      $dihedral);
    }
    if ($residue_name eq $last) {last;}
}








