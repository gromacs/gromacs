#!/usr/bin/perl -w

@oldatp = ( "br","c","c1", "c2","c3","ca","cc","ce","cg",
            "cl","cp","cu","cv","cx","cy","cz","f","h1","h2","h3","h4","h5","ha","hc","hn",
            "ho","hp","hs","hw","hx","i","n","n1","n2","n3","n4","na","nb","nc","ne","nh",
            "no","o","oh","on","os","p2","p3","p4","p5","pb","pc","pe","px","py","s","s2",
            "s4","s6","sh","ss","sx","sy");
            
$atp{"br"} = { elem => "Br", qq => "2", bondtype => "br", poltype => "br" };
$atp{"c"}  = { elem => "C",  qq => "2", bondtype => "c", poltype => "c" };
$atp{"c1"} = { elem => "C",  qq => "2", bondtype => "c", poltype => "c" };
$atp{"c2"} = { elem => "C",  qq => "2", bondtype => "c", poltype => "c" };
$atp{"c3"} = { elem => "C",  qq => "2", bondtype => "c", poltype => "c" };
$atp{"ca"} = { elem => "C",  qq => "2", bondtype => "c", poltype => "c" };
$atp{"cc"} = { elem => "C",  qq => "2", bondtype => "c", poltype => "c" };
$atp{"ce"} = { elem => "C",  qq => "2", bondtype => "c", poltype => "c" };
$atp{"cg"} = { elem => "C",  qq => "2", bondtype => "c", poltype => "c" };
$atp{"cl"} = { elem => "Cl", qq => "2", bondtype => "c", poltype => "c" };
$atp{"cp"} = { elem => "C",  qq => "2", bondtype => "c", poltype => "c" };
$atp{"cu"} = { elem => "C",  qq => "2", bondtype => "c", poltype => "c" };
$atp{"cv"} = { elem => "C",  qq => "2", bondtype => "c", poltype => "c" };
$atp{"cx"} = { elem => "C",  qq => "2", bondtype => "c", poltype => "c" };
$atp{"cy"} = { elem => "C",  qq => "2", bondtype => "c", poltype => "c" };
$atp{"cz"} = { elem => "C",  qq => "2", bondtype => "c", poltype => "c" };
$atp{"f"}  = { elem => "F",  qq => "2", bondtype => "f", poltype => "f" };
$atp{"h1"} = { elem => "H",  qq => "2", bondtype => "ha", poltype => "ha" };
$atp{"h2"} = { elem => "H",  qq => "2", bondtype => "ha", poltype => "ha" };
$atp{"h3"} = { elem => "H",  qq => "2", bondtype => "ha", poltype => "ha" };
$atp{"h4"} = { elem => "H",  qq => "2", bondtype => "ha", poltype => "ha" };
$atp{"h5"} = { elem => "H",  qq => "2", bondtype => "ha", poltype => "ha" };
$atp{"ha"} = { elem => "H",  qq => "2", bondtype => "ha", poltype => "ha" };
$atp{"hc"} = { elem => "H",  qq => "2", bondtype => "ha", poltype => "ha" };
$atp{"hn"} = { elem => "H",  qq => "2", bondtype => "hp", poltype => "hp" };
$atp{"ho"} = { elem => "H",  qq => "2", bondtype => "hp", poltype => "hp" };
$atp{"hp"} = { elem => "H",  qq => "2", bondtype => "hp", poltype => "hp" };
$atp{"hs"} = { elem => "H",  qq => "2", bondtype => "hp", poltype => "hp" };
$atp{"hw"} = { elem => "H",  qq => "2", bondtype => "hp", poltype => "hp" };
$atp{"hx"} = { elem => "H",  qq => "2", bondtype => "hp", poltype => "hp" };
$atp{"i"}  = { elem => "I",  qq => "2", bondtype => "i", poltype => "i" };
$atp{"n"}  = { elem => "N",  qq => "2", bondtype => "n", poltype => "n" };
$atp{"n1"} = { elem => "N",  qq => "2", bondtype => "n", poltype => "n" };
$atp{"n2"} = { elem => "N",  qq => "2", bondtype => "n", poltype => "n" };
$atp{"n3"} = { elem => "N",  qq => "2", bondtype => "n", poltype => "n" };
$atp{"n4"} = { elem => "N",  qq => "2", bondtype => "n", poltype => "n" };
$atp{"na"} = { elem => "N",  qq => "2", bondtype => "n", poltype => "n" };
$atp{"nb"} = { elem => "N",  qq => "2", bondtype => "n", poltype => "n" };
$atp{"nc"} = { elem => "N",  qq => "2", bondtype => "n", poltype => "n" };
$atp{"ne"} = { elem => "N",  qq => "2", bondtype => "n", poltype => "n" };
$atp{"nh"} = { elem => "N",  qq => "2", bondtype => "n", poltype => "n" };
$atp{"no"} = { elem => "N",  qq => "2", bondtype => "n", poltype => "n" };
$atp{"o"}  = { elem => "O",  qq => "2", bondtype => "o", poltype => "o" };
$atp{"oh"} = { elem => "O",  qq => "2", bondtype => "o", poltype => "o" };
$atp{"on"} = { elem => "O",  qq => "2", bondtype => "o", poltype => "o" };
$atp{"os"} = { elem => "O",  qq => "2", bondtype => "o", poltype => "o" };
$atp{"p2"} = { elem => "P",  qq => "2", bondtype => "p", poltype => "p" };
$atp{"p3"} = { elem => "P",  qq => "2", bondtype => "p", poltype => "p" };
$atp{"p4"} = { elem => "P",  qq => "2", bondtype => "p", poltype => "p" };
$atp{"p5"} = { elem => "P",  qq => "2", bondtype => "p", poltype => "p" };
$atp{"pb"} = { elem => "P",  qq => "2", bondtype => "p", poltype => "p" };
$atp{"pc"} = { elem => "P",  qq => "2", bondtype => "p", poltype => "p" };
$atp{"pe"} = { elem => "P",  qq => "2", bondtype => "p", poltype => "p" };
$atp{"px"} = { elem => "P",  qq => "2", bondtype => "p", poltype => "p" };
$atp{"py"} = { elem => "P",  qq => "2", bondtype => "p", poltype => "p" };
$atp{"s"}  = { elem => "S",  qq => "2", bondtype => "p", poltype => "s" };
$atp{"s2"} = { elem => "S",  qq => "2", bondtype => "s", poltype => "s" };
$atp{"s4"} = { elem => "S",  qq => "2", bondtype => "s", poltype => "s" };
$atp{"s6"} = { elem => "S",  qq => "2", bondtype => "s", poltype => "s" };
$atp{"sh"} = { elem => "S",  qq => "2", bondtype => "s", poltype => "s" };
$atp{"ss"} = { elem => "S",  qq => "2", bondtype => "s", poltype => "s" };
$atp{"sx"} = { elem => "S",  qq => "2", bondtype => "s", poltype => "s" };
$atp{"sy"} = { elem => "S",  qq => "2", bondtype => "s", poltype => "s" }; 

$mil{"hp"} = "H";
$mil{"ha"} = "H";
$mil{"c"} = "CTE";
$mil{"n"} = "NTE";
$mil{"f"} = "F";
$mil{"o"} = "OTE";
$mil{"p"} = "PTE";
$mil{"s"} = "STE";
$mil{"cl"} = "Cl";
$mil{"br"} = "Br";
$mil{"i"} = "I";

if ( 0 ) {
  foreach $tp ( @oldatp ) {
    $elem = substr($tp,0,1);
    $elem =~ tr/a-z/A-Z/;
    printf("\$atp{\"$tp\"} = { elem => \"%s\", qq => \"2\" };\n",$elem);
    
  }
  exit 1;
}

%ptp = ();
foreach $tp ( keys %atp ) {
  my $pt = $atp{$tp}->{'poltype'};
  $ptp{$pt} = $atp{$tp}->{'elem'};
}

%btp = ();
foreach $tp ( keys %atp ) {
  my $pt = $atp{$tp}->{'bondtype'};
  $btp{$pt} = 1;
}

$gt = "gentop.dat";
open(GT,">$gt") || die("Opening $gt");

printf(GT "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n");
printf(GT "<!DOCTYPE gentop.dtd PUBLIC \"gentop.dtd\" \"gentop.dtd\">\n");
printf(GT "<gentop>\n");

# Atom types
printf(GT "  <atomtypes forcefield=\"GAFF\" function=\"LJ_SR\" combination_rule=\"Arithmetic\" nexclusions=\"3\" fudgeQQ=\"1\" fudgeLJ=\"1\">\n");

foreach $tp ( keys %atp ) {
  printf(GT "    <atomtype elem=\"%s\" description=\"$tp\" atype=\"$tp\" ptype=\"%s\" btype=\"%s\" vdwparams=\"1 1\"/>\n",$atp{$tp}->{"elem"},$atp{$tp}->{"poltype"},$atp{$tp}->{"bondtype"});
}
printf(GT "  </atomtypes>\n");

# Pol types
printf(GT "  <poltypes polarizability_unit=\"A3\" reference=\"Maaren2013a\">\n");

foreach $tp ( keys %ptp ) {
  $mm = (defined $mil{$tp}) ? $mil{$tp} : $tp;
  printf(GT "    <poltype ptype=\"%s\" miller=\"%s\" bosque=\"%s\" polarizability=\"0\" sigma_pol=\"0\"/>\n",$tp,$mm,$ptp{$tp});
}
printf(GT "  </poltypes>\n");

printf(GT "  <bonding_rules>\n");
printf(GT "  </bonding_rules>\n");

printf(GT "  <gt_bonds length_unit=\"pm\" function=\"morse\">\n");
printf(GT "  </gt_bonds>\n");
printf(GT "  <gt_angles angle_unit=\"degree\" function=\"angles\">\n");
printf(GT "  </gt_angles>\n");
printf(GT "  <gt_dihedrals angle_unit=\"degree\" function=\"pdihs\">\n");
printf(GT "  </gt_dihedrals>\n");
printf(GT "  <gt_impropers angle_unit=\"degree\" function=\"idihs\">\n");
printf(GT "  </gt_impropers>\n");

my $bm = "bosque-miller.dat";
open(BM,"$bm") || die("Reading $bm");
@bbb = <BM>;
close BM;
printf(GT "@bbb");

@syms = ( ( "C", "H", "3" ),
          ( "N", "H", "3" ),
          ( "N", "H", "2" ),
          ( "N", "O", "2" ),
          ( "O", "H", "2" ),
          ( "C", "F", "3" ),
          ( "C", "Cl", "3" ),
          ( "C", "Br", "3" ),
          ( "C", "I", "3" ),
          ( "Si", "H", "3" ) );
printf(GT "  <symmetric_charges>\n");
for($k=0; ($k<=$#syms); $k+= 3) {
  printf(GT "    <sym_charge central=\"%s\" attached=\"%s\" numattach=\"%s\"/>\n",$syms[$k],$syms[$k+1],$syms[$k+2]);
}
printf(GT "  </symmetric_charges>\n");

printf(GT "  <eemprops>\n");
foreach $m ( "AXp", "AXg", "AXs" ) {
  foreach $tp ( keys %atp ) {
    printf(GT "    <eemprop model=\"$m\" name=\"$tp\" jaa=\"24.9083\" chi=\"7.49112\" zeta=\"0 40.1 14.2\" charge=\"%s -2 0\" row=\"1 1 1\"/>\n",$atp{$tp}->{"qq"});
  }
}
printf(GT "  </eemprops>\n");
printf(GT "</gentop>\n");
close GT;

printf("Generated $gt\n");
