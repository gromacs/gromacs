#!/usr/bin/perl

# usage: make_gromos_bon.pl ifp41a1.dat > gromos-bon.itp
# this generates gromos-bon.itp, a list of bonds, angles, impropers
# and dihedrals from ifp41a1.dat. 

while (<>) {
    if (/BONDTYPE/)         {$prefix = "#define gb_"; $bonds = 1; $go = 1;}
    if (/BONDANGLE/)        {$prefix = "#define ga_"; $bonds = 0; $angles = 1;}
    if (/IMPDIH/)           {$prefix = "#define gi_"; $angles= 0; $imp = 1;}
    if (/DIHEDRALTYPECODE/) {$prefix = "#define gd_"; $imp=0; $dih=1;}
    if (/SINGLE/)           {$go = 0;}

    if (/^#/ && $go) {tr/#/;/; print $_;}
    if (/^\d+/ && $go) {
	# need to switch the order of terms for bonds and angles for use
	# GROMACS
	if ($bonds || $angles ) {
	    ($nr,$k,$dist) = split(' ',$_);
	    print sprintf("$prefix%-3d  %10s  %10s\n", $nr, $dist, $k);
	}
	# impropers have the wrong units for k, convert to degrees
	if ($imp) {
	    ($nr,$k,$dist) = split(' ',$_);
	    $k = $k*180*180/(3.141593*3.141593);
	    print sprintf("$prefix%-3d  %10s  %10.5f\n", $nr, $dist, $k);
	}
	# same for dihedrals, also convert phase from cos(phi) to phi
	if ($dih) {
	    ($nr, $k, $phase, $mult) = split(' ',$_);
	    if ($phase > 0) {$phase = 0.0;} 
	    else {$phase = 180.0;}
	    print sprintf("$prefix%-3d %8.3f %10s %10s\n",$nr, $phase, 
			  $k, $mult);
	}
	if (/^[a-zA-Z]+/ && $go) {print ";" . $_;}
    }
}

# add stuff for the termina for now. Can be removed later if we have
# a working termini database

print "\n[ bondtypes ]\n";
print "NL     H       2    gb_2\n";
print "C      OM      2    gb_5\n";
print "OA     H       2    gb_1\n";
print "C      OA      2    gb_12\n";
print "C      O       2    gb_4\n";
print "S      S       2    gb_33\n";
print "NR     FE      2    gb_32\n";

print "\n[ angletypes ]\n";
print "H      NL     H     2   ga_9\n";
print "H      NL     CH1   2   ga_10\n";
print "CH1    C      OM    2   ga_21\n";
print "OM     C      OM    2   ga_37\n";
print "O      C      OA    2   ga_32\n";
print "C      OA     H     2   ga_11\n";
print "CH1    C      O     2   ga_29\n";
print "CH1    CH2    S     2   ga_15\n";
print "CH2    S      S     2   ga_5\n";
print "CH2    C     OM     2   ga_21\n";
print "CR1    NR    FE     2	ga_33\n";
print "NR     FE    NR     2   ga_16\n";

print "\n[ dihedraltypes ]\n";
print "S      S      1   gd_10\n";
print "NR     FE     1   gd_18\n";
print "CH2    S      1   gd_13\n";

