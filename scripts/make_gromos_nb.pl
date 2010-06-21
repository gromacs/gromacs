#!/usr/bin/perl

# usage: make_gromos_nb.pl gromos_atoms > gromos_nb.itp
# this script generates a GROMOS96 nonbonded forcefield file for GROMACS,
# with as input the file gromos_atoms, which contains records of the form
#:        1       O       0.04756 0.8611E-3       1.125E-3        0.0    
#: #CS6 CS12 parameters LJ14PAIR 
#:       0.04756 0.8611E-3 
#:  1   1   2   2   2   2   2   2   2   2   1   1   1   1   1   1   1   1   1   1
#:  2   2   2   2   2   2   2   1   1   1   1   1   2   2   1   1   1   1   1   1
#:  1   1   1
#: #---
# taken directly from ifp43a1.dat. For a description of what the numbers
# mean, see the GROMOS96 manual and the fie ifp43a1.dat

# set the input separator to #--- so we read one atom entry at a time
# make sure the records are actually separated by #---\n and not by say #--\n!!
# don't put a separator at the end of the file, only between records
$/= "#---\n"; 

# start arrays with 1 instead of 0
$[=1; 

$i=1;

# read in all entries (43 in ifp43a1.dat). This will BREAK for other input!

while(<>) {
    @lines = split('\n',$_);
    ($number[$i],$name[$i],$c6[$i],$c12_1[$i],$c12_2[$i],$c12_3[$i]) = 
	split(' ',$lines[1]);
    ($c6_14[$i],$c12_14[$i]) = split(' ',$lines[3]);
    $combination[$i] = $lines[4] . $lines[5] . $lines[6];

    # one type is called P,SI, the same LJ parameters for both P and SI
    # treat P,SI different: create a 44th type for Si, just a copy of Si,P
    # and rename P,SI to P
    if ($name[$i] =~ /P,SI/) {
	$number[44] = 44; $name[44] = "SI"; $c6[44] = $c6[$i];
	$c12_1[44] = $c12_1[$i]; $c12_2[44] = $c12_2[$i]; 
	$c12_3[44] = $c12_3[$i]; $combination[44] = $combination[$i];
	$name[$i] = "P";
	$P = $i; 
    }
    $i++;
}

# now $number[$i] has the atom nr, $name[$i] the name, $c6 the C6 LJ parameter
# $c12_1[$i], $c12_2[$i], $c12_3[$i] the three C12 parameters and 
# $combination[$i] the matrix with 43 elements that tells which C12 parameter
# you need in combination with each of the 44 atomtypes. $i runs from 1 to 44,
# one for each atom type. This goes nicely wrong because of the Stupid SI,P
# entry so we have to give an extra element 44 with the same value as that of
# P

# start printing to the output file: header for gromos-nb.itp
print "[ atomtypes ]\n";
print ";name        mass      charge   ptype            c6           c12\n";

# print out the atomtypes, plus the C6 and C12 for interactions with themselves
# the masses and charges are set to 0, the masses should come from gromos.atp

for ($j=1;$j<=44;$j++) {
# lookup the C12 with itself
    @c12types = split(' ',$combination[$j]);
    $c12types[44] = $c12types[$P]; # SI has the same as P
    if ($c12types[$j] == 1) 
    {$c12 = $c12_1[$j];}
    elsif ($c12types[$j] == 2)
    {$c12 = $c12_2[$j];}
    elsif ($c12types[$j] == 3)
    {$c12 = $c12_3[$j];}
    else {die "Error [atomtypes]: c12 type is not 1,2,3:j=$j,c12=$c12\n";}

    print sprintf("%5s     0.000      0.000     A  %10.8g  %10.8g\n", 
		  $name[$j],$c6[$j]*$c6[$j],$c12*$c12);
}

# Now make the LJ matrix. It's ok to do some double work in shell scripts, 
# trust me. 

print "\n";
print "[ nonbond_params ]\n";
print "; i    j func          c6           c12\n";

for ($j=1;$j<=44;$j++) {
    for ($k=1;$k<$j;$k++) {
# lookup the C12 of j for k
	@c12types_j = split(' ',$combination[$j]);
        $c12types_j[44] = $c12types_j[$P]; # SI has the same as P
	if ($c12types_j[$k] == 1) 
	{$c12_j = $c12_1[$j];}
	elsif ($c12types_j[$k] == 2)
	{$c12_j = $c12_2[$j];}
	elsif ($c12types_j[$k] == 3)
	{$c12_j = $c12_3[$j];}
	else {
	    die "Error [nonbond-params] j=$j,k=$k: c12 type is not one of 1,2,3\n";
	}
# lookup the C12 of k for j
	@c12types_k = split(' ',$combination[$k]);
        $c12types_k[44] = $c12types_k[$P]; # SI has the same as P
	if ($c12types_k[$j] == 1) 
	{$c12_k = $c12_1[$k];}
	elsif ($c12types_k[$j] == 2)
	{$c12_k = $c12_2[$k];}
	elsif ($c12types_k[$j] == 3)
	{$c12_k = $c12_3[$k];}
	else {
	    die "Error [nonbond-params] j=$j,k=$k: c12 type is not one of 1,2,3\n";
	}
    
	print sprintf("%8s %8s  1  %10.8g  %10.8g\n", $name[$j], $name[$k], 
		      $c6[$j]*$c6[$k], $c12_j*$c12_k);
    }
}
    
# Now do the same for the 1-4 interactions
print "\n";
print "[ pairtypes ]\n";
print "; i    j func          c6           c12\n";

for ($j=1;$j<=44;$j++) {
    for ($k=1;$k<=$j;$k++) {
	print sprintf("%8s %8s  1  %10.8g  %10.8g\n", $name[$j], $name[$k], 
		      $c6_14[$j]*$c6_14[$k], $c12_14[$j]*$c12_14[$k]);

    }
}







