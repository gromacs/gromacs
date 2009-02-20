shopt -s extglob
_anadock_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -ox -od -of -g -h -nice -noxvgr -free -norms -cutoff' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.pdb*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ox) COMPREPLY=( $(compgen -X '!*.pdb*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-od) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-of) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-g) COMPREPLY=( $(compgen -X '!*.log*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _anadock_compl anadock
shopt -s extglob
_do_dssp_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -ssdump -map -o -sc -a -ta -aa -h -nice -b -e -dt -tu -w -noxvgr -sss' -- $c)); return 0; fi
case "$p" in
-tu) COMPREPLY=( $(compgen -W ' ps fs ns us ms s ' -- $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa|gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ssdump) COMPREPLY=( $(compgen -X '!*.dat*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-map) COMPREPLY=( $(compgen -X '!*.map*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.xpm*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-sc) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-a) COMPREPLY=( $(compgen -X '!*.xpm*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ta) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-aa) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _do_dssp_compl do_dssp
shopt -s extglob
_editconf_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -n -o -mead -bf -h -nice -w -ndef -bt -box -angles -d -c -center -translate -rotate -princ -scale -density -pbc -grasp -rvdw -sig56 -vdwread -atom -legend -label' -- $c)); return 0; fi
case "$p" in
-bt) COMPREPLY=( $(compgen -W ' triclinic cubic dodecahedron octahedron ' -- $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|esp|tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|esp)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-mead) COMPREPLY=( $(compgen -X '!*.pqr*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-bf) COMPREPLY=( $(compgen -X '!*.dat*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _editconf_compl editconf
shopt -s extglob
_eneconv_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -o -h -nice -b -e -dt -offset -settime -nosort -scalefac -noerror' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.+(edr|ene)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.+(edr|ene)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _eneconv_compl eneconv
shopt -s extglob
_g_anaeig_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -v -v2 -f -s -n -eig -eig2 -comp -rmsf -proj -2d -3d -filt -extr -over -inpr -h -nice -b -e -dt -tu -w -noxvgr -first -last -skip -max -nframes -split -entropy -temp -nevskip' -- $c)); return 0; fi
case "$p" in
-tu) COMPREPLY=( $(compgen -W ' ps fs ns us ms s ' -- $c ));;
-v) COMPREPLY=( $(compgen -X '!*.+(trr|cpt|trj)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-v2) COMPREPLY=( $(compgen -X '!*.+(trr|cpt|trj)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa|gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-eig) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-eig2) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-comp) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-rmsf) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-proj) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-2d) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-3d) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|esp)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-filt) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-extr) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-over) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-inpr) COMPREPLY=( $(compgen -X '!*.xpm*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_anaeig_compl g_anaeig
shopt -s extglob
_g_analyze_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -ac -msd -cc -dist -av -ee -g -h -nice -w -noxvgr -notime -b -e -n -d -bw -errbar -integrate -aver_start -xydy -regression -luzar -temp -fitstart -smooth -filter -power -nosubav -oneacf -acflen -nonormalize -P -fitfn -ncskip -beginfit -endfit' -- $c)); return 0; fi
case "$p" in
-errbar) COMPREPLY=( $(compgen -W ' none stddev error 90 ' -- $c ));;
-P) COMPREPLY=( $(compgen -W ' 0 1 2 3 ' -- $c ));;
-fitfn) COMPREPLY=( $(compgen -W ' none exp aexp exp_exp vac exp5 exp7 exp9 ' -- $c ));;
-f) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ac) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-msd) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-cc) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-dist) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-av) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ee) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-g) COMPREPLY=( $(compgen -X '!*.log*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_analyze_compl g_analyze
shopt -s extglob
_g_angle_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -n -od -ov -of -ot -oh -oc -or -h -nice -b -e -dt -w -noxvgr -type -all -binwidth -noperiodic -chandler -avercorr -acflen -nonormalize -P -fitfn -ncskip -beginfit -endfit' -- $c)); return 0; fi
case "$p" in
-type) COMPREPLY=( $(compgen -W ' angle dihedral improper ryckaert-bellemans ' -- $c ));;
-P) COMPREPLY=( $(compgen -W ' 0 1 2 3 ' -- $c ));;
-fitfn) COMPREPLY=( $(compgen -W ' none exp aexp exp_exp vac exp5 exp7 exp9 ' -- $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-od) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ov) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-of) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ot) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-oh) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-oc) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-or) COMPREPLY=( $(compgen -X '!*.trr*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_angle_compl g_angle
shopt -s extglob
_g_bond_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -n -s -o -l -d -h -nice -b -e -dt -w -noxvgr -blen -tol -noaver -noaverdist' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa|gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-l) COMPREPLY=( $(compgen -X '!*.log*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-d) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_bond_compl g_bond
shopt -s extglob
_g_bundle_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -ol -od -oz -ot -otr -otl -ok -okr -okl -oa -h -nice -b -e -dt -tu -noxvgr -na -z' -- $c)); return 0; fi
case "$p" in
-tu) COMPREPLY=( $(compgen -W ' ps fs ns us ms s ' -- $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa|gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ol) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-od) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-oz) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ot) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-otr) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-otl) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ok) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-okr) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-okl) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-oa) COMPREPLY=( $(compgen -X '!*.pdb*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_bundle_compl g_bundle
shopt -s extglob
_g_chi_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -s -f -o -p -ss -jc -corr -g -ot -oh -rt -cp -h -nice -b -e -dt -w -noxvgr -r0 -phi -psi -omega -rama -viol -noperiodic -all -rad -shift -binwidth -core_rotamer -maxchi -nonormhisto -ramomega -bfact -chi_prod -HChi -bmax -acflen -nonormalize -P -fitfn -ncskip -beginfit -endfit' -- $c)); return 0; fi
case "$p" in
-maxchi) COMPREPLY=( $(compgen -W ' 0 1 2 3 4 5 6 ' -- $c ));;
-P) COMPREPLY=( $(compgen -W ' 0 1 2 3 ' -- $c ));;
-fitfn) COMPREPLY=( $(compgen -W ' none exp aexp exp_exp vac exp5 exp7 exp9 ' -- $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|esp|tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-p) COMPREPLY=( $(compgen -X '!*.pdb*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ss) COMPREPLY=( $(compgen -X '!*.dat*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-jc) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-corr) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-g) COMPREPLY=( $(compgen -X '!*.log*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ot) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-oh) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-rt) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-cp) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_chi_compl g_chi
shopt -s extglob
_g_cluster_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -dm -o -g -dist -ev -sz -tr -ntr -clid -cl -h -nice -b -e -dt -tu -w -noxvgr -dista -nlevels -cutoff -nofit -max -skip -av -wcl -nst -rmsmin -method -minstruct -binary -M -P -seed -niter -kT' -- $c)); return 0; fi
case "$p" in
-tu) COMPREPLY=( $(compgen -W ' ps fs ns us ms s ' -- $c ));;
-method) COMPREPLY=( $(compgen -W ' linkage jarvis-patrick monte-carlo diagonalization gromos ' -- $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa|gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-dm) COMPREPLY=( $(compgen -X '!*.xpm*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.xpm*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-g) COMPREPLY=( $(compgen -X '!*.log*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-dist) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ev) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-sz) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-tr) COMPREPLY=( $(compgen -X '!*.xpm*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ntr) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-clid) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-cl) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_cluster_compl g_cluster
shopt -s extglob
_g_clustsize_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -o -ow -nc -mc -ac -hc -temp -mcn -h -nice -b -e -dt -tu -w -noxvgr -cut -mol -nopbc -nskip -nlevels -ndf -rgblo -rgbhi' -- $c)); return 0; fi
case "$p" in
-tu) COMPREPLY=( $(compgen -W ' ps fs ns us ms s ' -- $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.tpr*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.xpm*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ow) COMPREPLY=( $(compgen -X '!*.xpm*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-nc) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-mc) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ac) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-hc) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-temp) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-mcn) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_clustsize_compl g_clustsize
shopt -s extglob
_g_confrms_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f1 -f2 -o -n1 -n2 -no -h -nice -w -one -nomw -pbc -nofit -name -label -bfac' -- $c)); return 0; fi
case "$p" in
-f1) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa|gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-f2) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|esp|tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|esp)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n1) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n2) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-no) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_confrms_compl g_confrms
shopt -s extglob
_g_covar_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -o -v -av -l -ascii -xpm -xpma -h -nice -b -e -dt -tu -noxvgr -nofit -ref -mwa -last -nopbc' -- $c)); return 0; fi
case "$p" in
-tu) COMPREPLY=( $(compgen -W ' ps fs ns us ms s ' -- $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa|gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-v) COMPREPLY=( $(compgen -X '!*.+(trr|cpt|trj)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-av) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|esp)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-l) COMPREPLY=( $(compgen -X '!*.log*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ascii) COMPREPLY=( $(compgen -X '!*.dat*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-xpm) COMPREPLY=( $(compgen -X '!*.xpm*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-xpma) COMPREPLY=( $(compgen -X '!*.xpm*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_covar_compl g_covar
shopt -s extglob
_g_current_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -s -n -f -o -caf -dsp -md -mj -mc -h -nice -b -e -dt -w -noxvgr -sh -nonojump -eps -bfit -efit -bvit -evit -tr -temp' -- $c)); return 0; fi
case "$p" in
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa|gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-caf) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-dsp) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-md) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-mj) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-mc) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_current_compl g_current
shopt -s extglob
_g_density_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -n -s -ei -o -h -nice -b -e -dt -w -noxvgr -d -sl -dens -ng -symm -center' -- $c)); return 0; fi
case "$p" in
-dens) COMPREPLY=( $(compgen -W ' mass number charge electron ' -- $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ei) COMPREPLY=( $(compgen -X '!*.dat*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_density_compl g_density
shopt -s extglob
_g_densmap_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -o -h -nice -b -e -dt -w -bin -aver -xmin -xmax -n1 -n2 -amax -rmax -mirror -unit -dmin -dmax' -- $c)); return 0; fi
case "$p" in
-aver) COMPREPLY=( $(compgen -W ' z y x ' -- $c ));;
-unit) COMPREPLY=( $(compgen -W ' nm-3 nm-2 count ' -- $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa|gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.xpm*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_densmap_compl g_densmap
shopt -s extglob
_g_dielectric_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -d -o -c -h -nice -b -e -dt -w -noxvgr -fft -nox1 -eint -bfit -efit -tail -A -tau1 -tau2 -eps0 -epsRF -fix -ffn -nsmooth' -- $c)); return 0; fi
case "$p" in
-ffn) COMPREPLY=( $(compgen -W ' none exp aexp exp_exp vac exp5 exp7 exp9 ' -- $c ));;
-f) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-d) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-c) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_dielectric_compl g_dielectric
shopt -s extglob
_g_dih_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -o -h -nice -b -e -dt -w -sa -mult' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.out*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_dih_compl g_dih
shopt -s extglob
_g_dipoles_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -enx -f -s -n -o -eps -a -d -c -g -adip -dip3d -cos -cmap -q -slab -h -nice -b -e -dt -w -noxvgr -mu -mumax -epsilonRF -skip -temp -corr -nopairs -ncos -axis -sl -gkratom -gkratom2 -rcmax -phi -nlevels -ndegrees -acflen -nonormalize -P -fitfn -ncskip -beginfit -endfit' -- $c)); return 0; fi
case "$p" in
-corr) COMPREPLY=( $(compgen -W ' none mol molsep total ' -- $c ));;
-P) COMPREPLY=( $(compgen -W ' 0 1 2 3 ' -- $c ));;
-fitfn) COMPREPLY=( $(compgen -W ' none exp aexp exp_exp vac exp5 exp7 exp9 ' -- $c ));;
-enx) COMPREPLY=( $(compgen -X '!*.+(edr|ene)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-eps) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-a) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-d) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-c) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-g) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-adip) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-dip3d) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-cos) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-cmap) COMPREPLY=( $(compgen -X '!*.xpm*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-q) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-slab) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_dipoles_compl g_dipoles
shopt -s extglob
_g_disre_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -s -f -ds -da -dn -dm -dr -l -n -q -c -x -h -nice -b -e -dt -w -noxvgr -ntop -maxdr -nlevels -nothird' -- $c)); return 0; fi
case "$p" in
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ds) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-da) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-dn) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-dm) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-dr) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-l) COMPREPLY=( $(compgen -X '!*.log*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-q) COMPREPLY=( $(compgen -X '!*.pdb*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-c) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-x) COMPREPLY=( $(compgen -X '!*.xpm*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_disre_compl g_disre
shopt -s extglob
_g_dist_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -o -lt -h -nice -b -e -dt -noxvgr -dist' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-lt) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_dist_compl g_dist
shopt -s extglob
_g_dyndom_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -o -n -h -nice -firstangle -lastangle -nframe -maxangle -trans -head -tail' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.pdb*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_dyndom_compl g_dyndom
shopt -s extglob
_genbox_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -cp -cs -ci -o -p -h -nice -box -nmol -try -seed -vdwd -shell -maxsol -vel' -- $c)); return 0; fi
case "$p" in
-cp) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|esp|tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-cs) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|esp|tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ci) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|esp|tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|esp)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-p) COMPREPLY=( $(compgen -X '!*.top*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _genbox_compl genbox
shopt -s extglob
_genconf_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -o -trj -h -nice -nbox -dist -seed -rot -shuffle -sort -block -nmolat -maxrot -norenumber' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|esp|tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|esp)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-trj) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _genconf_compl genconf
shopt -s extglob
_g_enemat_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -groups -eref -emat -etot -h -nice -b -e -dt -w -noxvgr -sum -skip -nomean -nlevels -max -min -nocoul -coulr -coul14 -nolj -lj -lj14 -bhamsr -bhamlr -nofree -temp' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.+(edr|ene)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-groups) COMPREPLY=( $(compgen -X '!*.dat*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-eref) COMPREPLY=( $(compgen -X '!*.dat*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-emat) COMPREPLY=( $(compgen -X '!*.xpm*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-etot) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_enemat_compl g_enemat
shopt -s extglob
_g_energy_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -f2 -s -o -viol -pairs -ora -ort -oda -odr -odt -oten -corr -vis -ravg -h -nice -b -e -w -noxvgr -fee -fetemp -zero -sum -dp -mutot -nouni -skip -aver -nmol -ndf -fluc -orinst -ovec -acflen -nonormalize -P -fitfn -ncskip -beginfit -endfit' -- $c)); return 0; fi
case "$p" in
-P) COMPREPLY=( $(compgen -W ' 0 1 2 3 ' -- $c ));;
-fitfn) COMPREPLY=( $(compgen -W ' none exp aexp exp_exp vac exp5 exp7 exp9 ' -- $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(edr|ene)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-f2) COMPREPLY=( $(compgen -X '!*.+(edr|ene)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-viol) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-pairs) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ora) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ort) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-oda) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-odr) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-odt) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-oten) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-corr) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-vis) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ravg) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_energy_compl g_energy
shopt -s extglob
_genion_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -s -table -n -o -g -pot -p -h -nice -noxvgr -np -pname -pq -nn -nname -nq -rmin -norandom -seed -scale -conc -neutral' -- $c)); return 0; fi
case "$p" in
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-table) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|esp)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-g) COMPREPLY=( $(compgen -X '!*.log*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-pot) COMPREPLY=( $(compgen -X '!*.pdb*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-p) COMPREPLY=( $(compgen -X '!*.top*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _genion_compl genion
shopt -s extglob
_genrestr_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -n -o -of -h -nice -fc -freeze -disre -disre_dist -disre_frac -disre_up2 -constr' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|esp|tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.itp*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-of) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _genrestr_compl genrestr
shopt -s extglob
_g_filter_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -ol -oh -h -nice -b -e -dt -w -nf -all -nonojump -fit' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa|gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ol) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-oh) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_filter_compl g_filter
shopt -s extglob
_g_gyrate_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -o -acf -h -nice -b -e -dt -w -noxvgr -nmol -q -p -moi -nz -acflen -nonormalize -P -fitfn -ncskip -beginfit -endfit' -- $c)); return 0; fi
case "$p" in
-P) COMPREPLY=( $(compgen -W ' 0 1 2 3 ' -- $c ));;
-fitfn) COMPREPLY=( $(compgen -W ' none exp aexp exp_exp vac exp5 exp7 exp9 ' -- $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa|gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-acf) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_gyrate_compl g_gyrate
shopt -s extglob
_g_h2order_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -n -nm -s -o -h -nice -b -e -dt -w -noxvgr -d -sl' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-nm) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_h2order_compl g_h2order
shopt -s extglob
_g_hbond_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -num -g -ac -dist -ang -hx -hbn -hbm -don -dan -life -nhbdist -h -nice -b -e -dt -noxvgr -ins -a -r -noda -r2 -abin -rbin -nonitacc -contact -shell -fitstart -temp -smooth -dump -max_hb -nomerge -acflen -nonormalize -P -fitfn -ncskip -beginfit -endfit' -- $c)); return 0; fi
case "$p" in
-P) COMPREPLY=( $(compgen -W ' 0 1 2 3 ' -- $c ));;
-fitfn) COMPREPLY=( $(compgen -W ' none exp aexp exp_exp vac exp5 exp7 exp9 ' -- $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-num) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-g) COMPREPLY=( $(compgen -X '!*.log*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ac) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-dist) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ang) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-hx) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-hbn) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-hbm) COMPREPLY=( $(compgen -X '!*.xpm*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-don) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-dan) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-life) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-nhbdist) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_hbond_compl g_hbond
shopt -s extglob
_g_helix_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -s -n -f -to -cz -co -h -nice -b -e -dt -w -r0 -q -noF -db -prop -ev -ahxstart -ahxend' -- $c)); return 0; fi
case "$p" in
-prop) COMPREPLY=( $(compgen -W ' RAD TWIST RISE LEN NHX DIP RMS CPHI RMSA PHI PSI HB3 HB4 HB5 CD222 ' -- $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-to) COMPREPLY=( $(compgen -X '!*.g87*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-cz) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|esp)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-co) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|esp)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_helix_compl g_helix
shopt -s extglob
_g_helixorient_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -s -f -n -oaxis -ocenter -orise -oradius -otwist -obending -otilt -orot -h -nice -b -e -dt -noxvgr -sidechain -incremental' -- $c)); return 0; fi
case "$p" in
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-oaxis) COMPREPLY=( $(compgen -X '!*.dat*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ocenter) COMPREPLY=( $(compgen -X '!*.dat*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-orise) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-oradius) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-otwist) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-obending) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-otilt) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-orot) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_helixorient_compl g_helixorient
shopt -s extglob
_g_lie_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -o -h -nice -b -e -dt -w -noxvgr -Elj -Eqq -Clj -Cqq -ligand' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.+(edr|ene)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_lie_compl g_lie
shopt -s extglob
_g_mdmat_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -mean -frames -no -h -nice -b -e -dt -noxvgr -t -nlevels' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa|gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-mean) COMPREPLY=( $(compgen -X '!*.xpm*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-frames) COMPREPLY=( $(compgen -X '!*.xpm*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-no) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_mdmat_compl g_mdmat
shopt -s extglob
_g_mindist_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -od -on -o -ox -or -h -nice -b -e -dt -tu -w -noxvgr -matrix -max -d -pi -split -ng -nopbc' -- $c)); return 0; fi
case "$p" in
-tu) COMPREPLY=( $(compgen -W ' ps fs ns us ms s ' -- $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa|gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-od) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-on) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.out*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ox) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-or) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_mindist_compl g_mindist
shopt -s extglob
_g_morph_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f1 -f2 -o -or -n -h -nice -w -noxvgr -ninterm -first -last -nofit' -- $c)); return 0; fi
case "$p" in
-f1) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|esp|tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-f2) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|esp|tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-or) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_morph_compl g_morph
shopt -s extglob
_g_msd_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -o -mol -pdb -h -nice -b -e -dt -tu -w -noxvgr -type -lateral -ten -ngroup -nomw -rmcomm -tpdb -trestart -beginfit -endfit' -- $c)); return 0; fi
case "$p" in
-tu) COMPREPLY=( $(compgen -W ' ps fs ns us ms s ' -- $c ));;
-type) COMPREPLY=( $(compgen -W ' no x y z ' -- $c ));;
-lateral) COMPREPLY=( $(compgen -W ' no x y z ' -- $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa|gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-mol) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-pdb) COMPREPLY=( $(compgen -X '!*.pdb*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_msd_compl g_msd
shopt -s extglob
_gmxcheck_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -f2 -s1 -s2 -c -e -e2 -n -m -h -nice -vdwfac -bonlo -bonhi -tol -ab -lastener' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-f2) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s1) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s2) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-c) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa|gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-e) COMPREPLY=( $(compgen -X '!*.+(edr|ene)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-e2) COMPREPLY=( $(compgen -X '!*.+(edr|ene)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-m) COMPREPLY=( $(compgen -X '!*.tex*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _gmxcheck_compl gmxcheck
shopt -s extglob
_gmxdump_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -s -f -e -cp -om -h -nice -nonr -sys' -- $c)); return 0; fi
case "$p" in
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-e) COMPREPLY=( $(compgen -X '!*.+(edr|ene)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-cp) COMPREPLY=( $(compgen -X '!*.cpt*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-om) COMPREPLY=( $(compgen -X '!*.mdp*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _gmxdump_compl gmxdump
shopt -s extglob
_g_nmeig_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -of -ol -v -h -nice -noxvgr -nom -first -last' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.mtx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa|gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-of) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ol) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-v) COMPREPLY=( $(compgen -X '!*.+(trr|cpt|trj)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_nmeig_compl g_nmeig
shopt -s extglob
_g_nmens_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -v -e -s -n -o -h -nice -noxvgr -temp -seed -num -first -last' -- $c)); return 0; fi
case "$p" in
-v) COMPREPLY=( $(compgen -X '!*.+(trr|cpt|trj)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-e) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa|gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_nmens_compl g_nmens
shopt -s extglob
_g_nmtraj_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -s -v -o -h -nice -eignr -phases -temp -amplitude -nframes' -- $c)); return 0; fi
case "$p" in
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa|gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-v) COMPREPLY=( $(compgen -X '!*.+(trr|cpt|trj)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_nmtraj_compl g_nmtraj
shopt -s extglob
_g_order_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -n -s -o -od -os -Sg -Sk -h -nice -b -e -dt -w -noxvgr -d -sl -szonly -unsat' -- $c)); return 0; fi
case "$p" in
-d) COMPREPLY=( $(compgen -W ' z x y ' -- $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-od) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-os) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-Sg) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-Sk) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_order_compl g_order
shopt -s extglob
_g_polystat_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -s -f -n -o -v -p -h -nice -b -e -dt -tu -w -noxvgr -nomw -pc' -- $c)); return 0; fi
case "$p" in
-tu) COMPREPLY=( $(compgen -W ' ps fs ns us ms s ' -- $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-v) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-p) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_polystat_compl g_polystat
shopt -s extglob
_g_potential_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -n -s -o -oc -of -h -nice -b -e -dt -w -noxvgr -d -sl -cb -ce -tz -spherical -ng -correct' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-oc) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-of) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_potential_compl g_potential
shopt -s extglob
_g_principal_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -a1 -a2 -a3 -om -h -nice -b -e -dt -tu -w -foo' -- $c)); return 0; fi
case "$p" in
-tu) COMPREPLY=( $(compgen -W ' ps fs ns us ms s ' -- $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa|gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-a1) COMPREPLY=( $(compgen -X '!*.dat*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-a2) COMPREPLY=( $(compgen -X '!*.dat*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-a3) COMPREPLY=( $(compgen -X '!*.dat*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-om) COMPREPLY=( $(compgen -X '!*.dat*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_principal_compl g_principal
shopt -s extglob
_g_rama_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -o -h -nice -b -e -dt -w -noxvgr' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_rama_compl g_rama
shopt -s extglob
_g_rdf_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -o -sq -cn -hq -h -nice -b -e -dt -w -noxvgr -bin -com -rdf -nopbc -nonorm -xy -cut -ng -fade -nlevel -startq -endq -energy' -- $c)); return 0; fi
case "$p" in
-rdf) COMPREPLY=( $(compgen -W ' atom mol_com mol_cog res_com res_cog ' -- $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa|gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-sq) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-cn) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-hq) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_rdf_compl g_rdf
shopt -s extglob
_g_rms_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -s -f -f2 -n -o -mir -a -dist -m -bin -bm -h -nice -b -e -dt -tu -w -noxvgr -what -nopbc -fit -prev -split -skip -skip2 -max -min -bmax -bmin -nomw -nlevels -ng' -- $c)); return 0; fi
case "$p" in
-tu) COMPREPLY=( $(compgen -W ' ps fs ns us ms s ' -- $c ));;
-what) COMPREPLY=( $(compgen -W ' rmsd rho rhosc ' -- $c ));;
-fit) COMPREPLY=( $(compgen -W ' rot+trans translation none ' -- $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa|gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-f2) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-mir) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-a) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-dist) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-m) COMPREPLY=( $(compgen -X '!*.xpm*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-bin) COMPREPLY=( $(compgen -X '!*.dat*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-bm) COMPREPLY=( $(compgen -X '!*.xpm*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_rms_compl g_rms
shopt -s extglob
_g_rmsdist_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -equiv -o -rms -scl -mean -nmr3 -nmr6 -noe -h -nice -b -e -dt -w -noxvgr -nlevels -max -nosumh' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa|gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-equiv) COMPREPLY=( $(compgen -X '!*.dat*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-rms) COMPREPLY=( $(compgen -X '!*.xpm*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-scl) COMPREPLY=( $(compgen -X '!*.xpm*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-mean) COMPREPLY=( $(compgen -X '!*.xpm*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-nmr3) COMPREPLY=( $(compgen -X '!*.xpm*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-nmr6) COMPREPLY=( $(compgen -X '!*.xpm*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-noe) COMPREPLY=( $(compgen -X '!*.dat*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_rmsdist_compl g_rmsdist
shopt -s extglob
_g_rmsf_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -q -oq -ox -o -od -oc -dir -h -nice -b -e -dt -w -noxvgr -res -aniso -nofit' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa|gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-q) COMPREPLY=( $(compgen -X '!*.pdb*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-oq) COMPREPLY=( $(compgen -X '!*.pdb*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ox) COMPREPLY=( $(compgen -X '!*.pdb*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-od) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-oc) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-dir) COMPREPLY=( $(compgen -X '!*.log*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_rmsf_compl g_rmsf
shopt -s extglob
_grompp_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -po -c -r -rb -n -p -pp -o -t -e -h -nice -nov -time -normvsbds -maxwarn -zero -norenum' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.mdp*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-po) COMPREPLY=( $(compgen -X '!*.mdp*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-c) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|esp|tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-r) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|esp|tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-rb) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|esp|tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-p) COMPREPLY=( $(compgen -X '!*.top*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-pp) COMPREPLY=( $(compgen -X '!*.top*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-t) COMPREPLY=( $(compgen -X '!*.+(trr|cpt|trj)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-e) COMPREPLY=( $(compgen -X '!*.+(edr|ene)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _grompp_compl grompp
shopt -s extglob
_g_rotacf_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -o -h -nice -b -e -dt -w -noxvgr -d -noaver -acflen -nonormalize -P -fitfn -ncskip -beginfit -endfit' -- $c)); return 0; fi
case "$p" in
-P) COMPREPLY=( $(compgen -W ' 0 1 2 3 ' -- $c ));;
-fitfn) COMPREPLY=( $(compgen -W ' none exp aexp exp_exp vac exp5 exp7 exp9 ' -- $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_rotacf_compl g_rotacf
shopt -s extglob
_g_saltbr_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -h -nice -b -e -dt -t -sep' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_saltbr_compl g_saltbr
shopt -s extglob
_g_sas_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -o -or -oa -tv -q -n -i -h -nice -b -e -dt -w -noxvgr -probe -ndots -qmax -f_index -minarea -nopbc -noprot -dgs' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-or) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-oa) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-tv) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-q) COMPREPLY=( $(compgen -X '!*.pdb*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-i) COMPREPLY=( $(compgen -X '!*.itp*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_sas_compl g_sas
shopt -s extglob
_g_sdf_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -n -s -o -r -h -nice -b -e -dt -mode -triangle -dtri -bin -grid' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa|gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.dat*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-r) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|esp)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_sdf_compl g_sdf
shopt -s extglob
_g_sgangle_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -n -s -oa -od -od1 -od2 -h -nice -b -e -dt -w -noxvgr -one -z' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-oa) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-od) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-od1) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-od2) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_sgangle_compl g_sgangle
shopt -s extglob
_g_sham_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -ge -ene -dist -histo -bin -lp -ls -lsh -lss -map -ls3 -mdata -g -h -nice -w -noxvgr -notime -b -e -ttol -n -d -bw -nosham -tsham -pmin -dim -ngrid -xmin -xmax -pmax -gmax -emin -emax -nlevels -mname' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ge) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ene) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-dist) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-histo) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-bin) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-lp) COMPREPLY=( $(compgen -X '!*.xpm*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ls) COMPREPLY=( $(compgen -X '!*.xpm*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-lsh) COMPREPLY=( $(compgen -X '!*.xpm*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-lss) COMPREPLY=( $(compgen -X '!*.xpm*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-map) COMPREPLY=( $(compgen -X '!*.xpm*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ls3) COMPREPLY=( $(compgen -X '!*.pdb*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-mdata) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-g) COMPREPLY=( $(compgen -X '!*.log*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_sham_compl g_sham
shopt -s extglob
_g_sorient_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -o -no -ro -co -rc -h -nice -b -e -dt -w -noxvgr -com -v23 -rmin -rmax -cbin -rbin -pbc' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa|gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-no) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ro) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-co) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-rc) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_sorient_compl g_sorient
shopt -s extglob
_g_spatial_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -dm -o -g -dist -ev -sz -tr -ntr -clid -cl -h -nice -b -e -dt -tu -w -noxvgr -dista -nlevels -cutoff -nofit -max -skip -av -wcl -nst -rmsmin -method -minstruct -binary -M -P -seed -niter -kT' -- $c)); return 0; fi
case "$p" in
-tu) COMPREPLY=( $(compgen -W ' ps fs ns us ms s ' -- $c ));;
-method) COMPREPLY=( $(compgen -W ' linkage jarvis-patrick monte-carlo diagonalization gromos ' -- $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa|gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-dm) COMPREPLY=( $(compgen -X '!*.xpm*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.xpm*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-g) COMPREPLY=( $(compgen -X '!*.log*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-dist) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ev) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-sz) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-tr) COMPREPLY=( $(compgen -X '!*.xpm*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ntr) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-clid) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-cl) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_spatial_compl g_spatial
shopt -s extglob
_g_spol_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -o -h -nice -b -e -dt -w -noxvgr -com -refat -rmin -rmax -dip -bw' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_spol_compl g_spol
shopt -s extglob
_g_tcaf_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -ot -oa -o -of -oc -ov -h -nice -b -e -dt -w -noxvgr -mol -k34 -wt -acflen -nonormalize -P -fitfn -ncskip -beginfit -endfit' -- $c)); return 0; fi
case "$p" in
-P) COMPREPLY=( $(compgen -W ' 0 1 2 3 ' -- $c ));;
-fitfn) COMPREPLY=( $(compgen -W ' none exp aexp exp_exp vac exp5 exp7 exp9 ' -- $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(trr|cpt|trj)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa|gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ot) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-oa) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-of) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-oc) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ov) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_tcaf_compl g_tcaf
shopt -s extglob
_g_traj_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -ox -oxt -ov -of -ob -ot -ekt -ekr -vd -cv -cf -av -af -h -nice -b -e -dt -tu -w -noxvgr -com -mol -nojump -nox -noy -noz -ng -len -bin -scale' -- $c)); return 0; fi
case "$p" in
-tu) COMPREPLY=( $(compgen -W ' ps fs ns us ms s ' -- $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa|gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ox) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-oxt) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ov) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-of) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ob) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ot) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ekt) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ekr) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-vd) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-cv) COMPREPLY=( $(compgen -X '!*.pdb*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-cf) COMPREPLY=( $(compgen -X '!*.pdb*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-av) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-af) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_traj_compl g_traj
shopt -s extglob
_g_tune_pme_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -p -so -s -o -x -cpi -cpo -c -e -g -dgdl -field -table -tablep -tableb -rerun -tpi -tpid -ei -eo -j -jo -ffout -devout -runav -px -pf -mtx -dn -h -nice -noxvgr -np -r -max -min -fac -ntpr -four -steps -simsteps -launch -deffnm -ddorder -noddcheck -rdd -rcon -dlb -dds -nosum -v -nocompact -seppot -pforce -reprod -cpt -append -maxh -multi -replex -reseed -glas -ionize' -- $c)); return 0; fi
case "$p" in
-ddorder) COMPREPLY=( $(compgen -W ' interleave pp_pme cartesian ' -- $c ));;
-dlb) COMPREPLY=( $(compgen -W ' auto no yes ' -- $c ));;
-p) COMPREPLY=( $(compgen -X '!*.out*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-so) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.+(trr|cpt|trj)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-x) COMPREPLY=( $(compgen -X '!*.xtc*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-cpi) COMPREPLY=( $(compgen -X '!*.cpt*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-cpo) COMPREPLY=( $(compgen -X '!*.cpt*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-c) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|esp|xyz)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-e) COMPREPLY=( $(compgen -X '!*.edr*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-g) COMPREPLY=( $(compgen -X '!*.log*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-dgdl) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-field) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-table) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-tablep) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-tableb) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-rerun) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-tpi) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-tpid) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ei) COMPREPLY=( $(compgen -X '!*.edi*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-eo) COMPREPLY=( $(compgen -X '!*.edo*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-j) COMPREPLY=( $(compgen -X '!*.gct*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-jo) COMPREPLY=( $(compgen -X '!*.gct*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ffout) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-devout) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-runav) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-px) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-pf) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-mtx) COMPREPLY=( $(compgen -X '!*.mtx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-dn) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_tune_pme_compl g_tune_pme
shopt -s extglob
_g_vanhove_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -om -or -ot -h -nice -b -e -dt -w -noxvgr -sqrt -fm -rmax -rbin -mmax -nlevels -nr -fr -rt -ft' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa|gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-om) COMPREPLY=( $(compgen -X '!*.xpm*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-or) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ot) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_vanhove_compl g_vanhove
shopt -s extglob
_g_velacc_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -o -h -nice -b -e -dt -w -noxvgr -m -mol -acflen -nonormalize -P -fitfn -ncskip -beginfit -endfit' -- $c)); return 0; fi
case "$p" in
-P) COMPREPLY=( $(compgen -W ' 0 1 2 3 ' -- $c ));;
-fitfn) COMPREPLY=( $(compgen -W ' none exp aexp exp_exp vac exp5 exp7 exp9 ' -- $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(trr|cpt|trj)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa|gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_velacc_compl g_velacc
shopt -s extglob
_g_wham_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -ix -if -it -ip -o -hist -bsres -bsprof -tab -wcorr -h -nice -noxvgr -min -max -noauto -bins -temp -tol -v -b -e -dt -histonly -boundsonly -nolog -unit -zprof0 -cycl -alpha -flip -hist-eq -nBootstrap -bs-dt -bs-seed -nohistbs -histbs-block -vbs' -- $c)); return 0; fi
case "$p" in
-unit) COMPREPLY=( $(compgen -W ' kJ kCal kT ' -- $c ));;
-cycl) COMPREPLY=( $(compgen -W ' no yes weighted ' -- $c ));;
-ix) COMPREPLY=( $(compgen -X '!*.dat*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-if) COMPREPLY=( $(compgen -X '!*.dat*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-it) COMPREPLY=( $(compgen -X '!*.dat*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ip) COMPREPLY=( $(compgen -X '!*.dat*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-hist) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-bsres) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-bsprof) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-tab) COMPREPLY=( $(compgen -X '!*.dat*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-wcorr) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_wham_compl g_wham
shopt -s extglob
_highway_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -h -nice' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.dat*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _highway_compl highway
shopt -s extglob
_make_edi_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -eig -s -n -tar -ori -o -h -nice -noxvgr -mon -linfix -linacc -flood -radfix -radacc -radcon -outfrq -slope -maxedsteps -deltaF0 -deltaF -tau -eqsteps -Eflnull -T -alpha -linstep -accdir -radstep -restrain -hessian -harmonic' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.+(trr|cpt|trj)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-eig) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa|gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-tar) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|esp|tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ori) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|esp|tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.edi*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _make_edi_compl make_edi
shopt -s extglob
_make_ndx_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -n -o -h -nice -natoms' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|esp|tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _make_ndx_compl make_ndx
shopt -s extglob
_mdrun_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -s -o -x -cpi -cpo -c -e -g -dgdl -field -table -tablep -tableb -rerun -tpi -tpid -ei -eo -j -jo -ffout -devout -runav -px -pf -mtx -dn -h -nice -deffnm -noxvgr -pd -dd -npme -ddorder -noddcheck -rdd -rcon -dlb -dds -nosum -v -nocompact -seppot -pforce -reprod -cpt -append -maxh -multi -replex -reseed -glas -ionize' -- $c)); return 0; fi
case "$p" in
-ddorder) COMPREPLY=( $(compgen -W ' interleave pp_pme cartesian ' -- $c ));;
-dlb) COMPREPLY=( $(compgen -W ' auto no yes ' -- $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.+(trr|cpt|trj)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-x) COMPREPLY=( $(compgen -X '!*.xtc*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-cpi) COMPREPLY=( $(compgen -X '!*.cpt*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-cpo) COMPREPLY=( $(compgen -X '!*.cpt*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-c) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|esp)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-e) COMPREPLY=( $(compgen -X '!*.+(edr|ene)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-g) COMPREPLY=( $(compgen -X '!*.log*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-dgdl) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-field) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-table) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-tablep) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-tableb) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-rerun) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-tpi) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-tpid) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ei) COMPREPLY=( $(compgen -X '!*.edi*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-eo) COMPREPLY=( $(compgen -X '!*.edo*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-j) COMPREPLY=( $(compgen -X '!*.gct*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-jo) COMPREPLY=( $(compgen -X '!*.gct*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ffout) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-devout) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-runav) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-px) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-pf) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-mtx) COMPREPLY=( $(compgen -X '!*.mtx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-dn) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _mdrun_compl mdrun
shopt -s extglob
_mk_angndx_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -s -n -h -nice -type -nohyd' -- $c)); return 0; fi
case "$p" in
-type) COMPREPLY=( $(compgen -W ' angle dihedral improper ryckaert-bellemans ' -- $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _mk_angndx_compl mk_angndx
shopt -s extglob
_ngmx_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -h -nice -b -e -dt' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _ngmx_compl ngmx
shopt -s extglob
_pdb2gmx_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -o -p -i -n -q -h -nice -merge -ff -water -inter -ss -ter -lys -arg -asp -glu -gln -his -angle -dist -una -ignh -missing -v -posrefc -vsite -heavyh -deuterate' -- $c)); return 0; fi
case "$p" in
-water) COMPREPLY=( $(compgen -W ' spc spce tip3p tip4p tip5p f3c ' -- $c ));;
-vsite) COMPREPLY=( $(compgen -W ' none hydrogens aromatics ' -- $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|esp|tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|esp)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-p) COMPREPLY=( $(compgen -X '!*.top*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-i) COMPREPLY=( $(compgen -X '!*.itp*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-q) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|esp)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _pdb2gmx_compl pdb2gmx
shopt -s extglob
_protonate_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -s -f -n -o -h -nice -b -e -dt' -- $c)); return 0; fi
case "$p" in
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa|gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _protonate_compl protonate
shopt -s extglob
_sigeps_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -o -h -nice -w -noxvgr -c6 -cn -pow -sig -eps -A -B -C -qi -qj -sigfac' -- $c)); return 0; fi
case "$p" in
-o) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _sigeps_compl sigeps
shopt -s extglob
_tpbconv_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -s -f -e -n -o -h -nice -nsteps -runtime -time -extend -until -zeroq -nocont' -- $c)); return 0; fi
case "$p" in
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(trr|cpt|trj)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-e) COMPREPLY=( $(compgen -X '!*.+(edr|ene)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _tpbconv_compl tpbconv
shopt -s extglob
_trjcat_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -o -n -demux -h -nice -tu -noxvgr -b -e -dt -prec -novel -settime -nosort -keeplast -cat' -- $c)); return 0; fi
case "$p" in
-tu) COMPREPLY=( $(compgen -W ' ps fs ns us ms s ' -- $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-demux) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _trjcat_compl trjcat
shopt -s extglob
_trjconv_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -o -s -n -fr -sub -drop -h -nice -b -e -tu -w -noxvgr -skip -dt -dump -t0 -timestep -pbc -ur -center -boxcenter -box -trans -shift -fit -ndec -novel -force -trunc -exec -app -split -sep -nzero -ter -dropunder -dropover' -- $c)); return 0; fi
case "$p" in
-tu) COMPREPLY=( $(compgen -W ' ps fs ns us ms s ' -- $c ));;
-pbc) COMPREPLY=( $(compgen -W ' none mol res atom nojump cluster whole ' -- $c ));;
-ur) COMPREPLY=( $(compgen -W ' rect tric compact ' -- $c ));;
-boxcenter) COMPREPLY=( $(compgen -W ' tric rect zero ' -- $c ));;
-fit) COMPREPLY=( $(compgen -W ' none rot+trans rotxy+transxy translation transxy progressive ' -- $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa|gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-fr) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-sub) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-drop) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _trjconv_compl trjconv
shopt -s extglob
_trjorder_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -o -nshell -h -nice -b -e -dt -noxvgr -na -da -com -r' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa|gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-nshell) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _trjorder_compl trjorder
shopt -s extglob
_wheel_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -o -h -nice -r0 -rot0 -T -nonn' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.dat*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.eps*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _wheel_compl wheel
shopt -s extglob
_x2top_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -o -r -h -nice -ff -v -nexcl -noH14 -alldih -remdih -nopairs -name -nopbc -pdbq -noparam -noround -kb -kt -kp' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|esp|tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.top*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-r) COMPREPLY=( $(compgen -X '!*.rtp*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _x2top_compl x2top
shopt -s extglob
_xpm2ps_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -f2 -di -do -o -xpm -h -nice -w -noframe -title -yonce -legend -diag -size -bx -by -rainbow -gradient -skip -zeroline -legoffset -combine -cmin -cmax' -- $c)); return 0; fi
case "$p" in
-title) COMPREPLY=( $(compgen -W ' top once ylabel none ' -- $c ));;
-legend) COMPREPLY=( $(compgen -W ' both first second none ' -- $c ));;
-diag) COMPREPLY=( $(compgen -W ' first second none ' -- $c ));;
-rainbow) COMPREPLY=( $(compgen -W ' no blue red ' -- $c ));;
-combine) COMPREPLY=( $(compgen -W ' halves add sub mult div ' -- $c ));;
-f) COMPREPLY=( $(compgen -X '!*.xpm*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-f2) COMPREPLY=( $(compgen -X '!*.xpm*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-di) COMPREPLY=( $(compgen -X '!*.m2p*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-do) COMPREPLY=( $(compgen -X '!*.m2p*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.eps*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-xpm) COMPREPLY=( $(compgen -X '!*.xpm*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _xpm2ps_compl xpm2ps
shopt -s extglob
_xrama_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -h -nice -b -e -dt' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _xrama_compl xrama
