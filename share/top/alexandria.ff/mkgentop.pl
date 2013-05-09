#!/usr/bin/perl -w

@oldatp = ( "br","c","c1", "c2","c3","ca","cc","ce","cg",
            "cl","cp","cu","cv","cx","cy","cz","f","h1","h2","h3","h4","h5","ha","hc","hn",
            "ho","hp","hs","hw","hx","i","n","n1","n2","n3","n4","na","nb","nc","ne","nh",
            "no","o","oh","on","os","p2","p3","p4","p5","pb","pc","pe","px","py","s","s2",
            "s4","s6","sh","ss","sx","sy");
            
$atp{"br"} = { elem => "Br", qq => "2" };
$atp{"c"}  = { elem => "C", qq => "2" };
$atp{"c1"} = { elem => "C", qq => "2" };
$atp{"c2"} = { elem => "C", qq => "2" };
$atp{"c3"} = { elem => "C", qq => "2" };
$atp{"ca"} = { elem => "C", qq => "2" };
$atp{"cc"} = { elem => "C", qq => "2" };
$atp{"ce"} = { elem => "C", qq => "2" };
$atp{"cg"} = { elem => "C", qq => "2" };
$atp{"cl"} = { elem => "Cl", qq => "2" };
$atp{"cp"} = { elem => "C", qq => "2" };
$atp{"cu"} = { elem => "C", qq => "2" };
$atp{"cv"} = { elem => "C", qq => "2" };
$atp{"cx"} = { elem => "C", qq => "2" };
$atp{"cy"} = { elem => "C", qq => "2" };
$atp{"cz"} = { elem => "C", qq => "2" };
$atp{"f"}  = { elem => "F", qq => "2" };
$atp{"h1"} = { elem => "H", qq => "2" };
$atp{"h2"} = { elem => "H", qq => "2" };
$atp{"h3"} = { elem => "H", qq => "2" };
$atp{"h4"} = { elem => "H", qq => "2" };
$atp{"h5"} = { elem => "H", qq => "2" };
$atp{"ha"} = { elem => "H", qq => "2" };
$atp{"hc"} = { elem => "H", qq => "2" };
$atp{"hn"} = { elem => "H", qq => "2" };
$atp{"ho"} = { elem => "H", qq => "2" };
$atp{"hp"} = { elem => "H", qq => "2" };
$atp{"hs"} = { elem => "H", qq => "2" };
$atp{"hw"} = { elem => "H", qq => "2" };
$atp{"hx"} = { elem => "H", qq => "2" };
$atp{"i"} = { elem => "I", qq => "2" };
$atp{"n"} = { elem => "N", qq => "2" };
$atp{"n1"} = { elem => "N", qq => "2" };
$atp{"n2"} = { elem => "N", qq => "2" };
$atp{"n3"} = { elem => "N", qq => "2" };
$atp{"n4"} = { elem => "N", qq => "2" };
$atp{"na"} = { elem => "N", qq => "2" };
$atp{"nb"} = { elem => "N", qq => "2" };
$atp{"nc"} = { elem => "N", qq => "2" };
$atp{"ne"} = { elem => "N", qq => "2" };
$atp{"nh"} = { elem => "N", qq => "2" };
$atp{"no"} = { elem => "N", qq => "2" };
$atp{"o"} = { elem => "O", qq => "2" };
$atp{"oh"} = { elem => "O", qq => "2" };
$atp{"on"} = { elem => "O", qq => "2" };
$atp{"os"} = { elem => "O", qq => "2" };
$atp{"p2"} = { elem => "P", qq => "2" };
$atp{"p3"} = { elem => "P", qq => "2" };
$atp{"p4"} = { elem => "P", qq => "2" };
$atp{"p5"} = { elem => "P", qq => "2" };
$atp{"pb"} = { elem => "P", qq => "2" };
$atp{"pc"} = { elem => "P", qq => "2" };
$atp{"pe"} = { elem => "P", qq => "2" };
$atp{"px"} = { elem => "P", qq => "2" };
$atp{"py"} = { elem => "P", qq => "2" };
$atp{"s"} = { elem => "S", qq => "2" };
$atp{"s2"} = { elem => "S", qq => "2" };
$atp{"s4"} = { elem => "S", qq => "2" };
$atp{"s6"} = { elem => "S", qq => "2" };
$atp{"sh"} = { elem => "S", qq => "2" };
$atp{"ss"} = { elem => "S", qq => "2" };
$atp{"sx"} = { elem => "S", qq => "2" };
$atp{"sy"} = { elem => "S", qq => "2" };

if ( 0 ) {
  foreach $tp ( @oldatp ) {
    $elem = substr($tp,0,1);
    $elem =~ tr/a-z/A-Z/;
    printf("\$atp{\"$tp\"} = { elem => \"%s\", qq => \"2\" };\n",$elem);
    
  }
  exit 1;
}

$gt = "gentop.dat";
open(GT,">$gt") || die("Opening $gt");

printf(GT "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n");
printf(GT "<!DOCTYPE gentop.dtd PUBLIC \"gentop.dtd\" \"gentop.dtd\">\n");
printf(GT "<gentop>\n");
printf(GT "  <atomtypes forcefield=\"GAFF\" polarizability_unit=\"A3\" function=\"LJ_SR\" combination_rule=\"Arithmetic\" nexclusions=\"3\" fudgeQQ=\"0\" fudgeLJ=\"0\">\n");

foreach $tp ( keys %atp ) {
  printf(GT "    <atype elem=\"%s\" description=\"$tp\" gt_type=\"$tp\" charge=\"0\" polarizability=\"0\" sigma_pol=\"0\" vdwparams=\"1 1\"/>\n",$atp{$tp}->{"elem"});
}
printf(GT "  </atomtypes>\n");
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
