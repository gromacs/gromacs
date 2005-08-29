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
_cdist_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -s -g -q -d -o -n -dom -h -nice -noengh -bm -am -pm -rr -ar -er -vm -lm -il -dm -im -nm -hm -hb -nobon -nonb -measure -maxdist -add -vir -sm' -- $c)); return 0; fi
case "$p" in
-sm) COMPREPLY=( $(compgen -W ' none tri tetra ' -- $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa|gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-g) COMPREPLY=( $(compgen -X '!*.log*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-q) COMPREPLY=( $(compgen -X '!*.pdb*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-d) COMPREPLY=( $(compgen -X '!*.dat*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.dat*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-dom) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _cdist_compl cdist
shopt -s extglob
_disco_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -g -f -d -do -c -center -n -o -keep -viol -h -nice -nf -nit -nov -nochiral -nopep -lower -weighted -dump -cubic -explicit -fit -nbcheck -nstprint -ranlist -noranlistfirst -lowdev -seed -box -grow' -- $c)); return 0; fi
case "$p" in
-g) COMPREPLY=( $(compgen -X '!*.log*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|esp|tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-d) COMPREPLY=( $(compgen -X '!*.dat*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-do) COMPREPLY=( $(compgen -X '!*.dat*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-c) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|esp)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-center) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|esp)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-keep) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-viol) COMPREPLY=( $(compgen -X '!*.pdb*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _disco_compl disco
shopt -s extglob
_do_dssp_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -ssdump -map -o -sc -a -ta -aa -h -nice -b -e -dt -tu -w -noxvgr -sss' -- $c)); return 0; fi
case "$p" in
-tu) COMPREPLY=( $(compgen -W ' ps fs ns us ms s m h ' -- $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
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
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -n -o -mead -bf -h -nice -w -ndef -bt -box -angles -d -c -center -translate -rotate -princ -scale -density -novol -pbc -grasp -rvdw -sig56 -vdwread -atom -legend -label' -- $c)); return 0; fi
case "$p" in
-bt) COMPREPLY=( $(compgen -W ' tric cubic dodecahedron octahedron ' -- $c ));;
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
_ffscan_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -s -g -table -parm -ga -c -e -o -h -nice -noxvgr -tol -fmax -nocomb -npow -logeps -v -epot -fepot -pres -fpres -fmsf -molsize -nmol' -- $c)); return 0; fi
case "$p" in
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-g) COMPREPLY=( $(compgen -X '!*.log*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-table) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-parm) COMPREPLY=( $(compgen -X '!*.dat*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ga) COMPREPLY=( $(compgen -X '!*.dat*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-c) COMPREPLY=( $(compgen -X '!*.gro*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-e) COMPREPLY=( $(compgen -X '!*.+(edr|ene)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.+(trr|trj)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _ffscan_compl ffscan
shopt -s extglob
_g_anaeig_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -v -v2 -f -s -n -eig -eig2 -comp -rmsf -proj -2d -3d -filt -extr -over -inpr -h -nice -b -e -dt -tu -w -noxvgr -first -last -skip -max -nframes -split' -- $c)); return 0; fi
case "$p" in
-tu) COMPREPLY=( $(compgen -W ' ps fs ns us ms s m h ' -- $c ));;
-v) COMPREPLY=( $(compgen -X '!*.+(trr|trj)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-v2) COMPREPLY=( $(compgen -X '!*.+(trr|trj)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa|gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-eig) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-eig2) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-comp) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-rmsf) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-proj) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-2d) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-3d) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|esp)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-filt) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-extr) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-over) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-inpr) COMPREPLY=( $(compgen -X '!*.xpm*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_anaeig_compl g_anaeig
shopt -s extglob
_g_analyze_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -ac -msd -cc -dist -av -ee -g -h -nice -w -noxvgr -notime -b -e -n -d -bw -errbar -integrate -aver_start -xydy -filter -power -nosubav -oneacf -acflen -nonormalize -P -fitfn -ncskip -beginfit -endfit' -- $c)); return 0; fi
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
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -od -ov -of -ot -oh -oc -or -h -nice -b -e -dt -w -noxvgr -type -all -binwidth -chandler -avercorr -acflen -nonormalize -P -fitfn -ncskip -beginfit -endfit' -- $c)); return 0; fi
case "$p" in
-type) COMPREPLY=( $(compgen -W ' angle dihedral improper ryckaert-bellemans ' -- $c ));;
-P) COMPREPLY=( $(compgen -W ' 0 1 2 3 ' -- $c ));;
-fitfn) COMPREPLY=( $(compgen -W ' none exp aexp exp_exp vac exp5 exp7 exp9 ' -- $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
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
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
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
-tu) COMPREPLY=( $(compgen -W ' ps fs ns us ms s m h ' -- $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
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
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -s -f -o -p -ss -jc -corr -g -ot -oh -rt -cp -h -nice -b -e -dt -w -noxvgr -r0 -phi -psi -omega -rama -viol -all -rad -shift -binwidth -core_rotamer -maxchi -nonormhisto -ramomega -bfact -chi_prod -HChi -bmax -acflen -nonormalize -P -fitfn -ncskip -beginfit -endfit' -- $c)); return 0; fi
case "$p" in
-maxchi) COMPREPLY=( $(compgen -W ' 0 1 2 3 4 5 6 ' -- $c ));;
-P) COMPREPLY=( $(compgen -W ' 0 1 2 3 ' -- $c ));;
-fitfn) COMPREPLY=( $(compgen -W ' none exp aexp exp_exp vac exp5 exp7 exp9 ' -- $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|esp|tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
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
-tu) COMPREPLY=( $(compgen -W ' ps fs ns us ms s m h ' -- $c ));;
-method) COMPREPLY=( $(compgen -W ' linkage jarvis-patrick monte-carlo diagonalization gromos ' -- $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
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
-cl) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_cluster_compl g_cluster
shopt -s extglob
_g_clustsize_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -o -ow -nc -mc -ac -hc -temp -mcn -h -nice -b -e -dt -tu -w -noxvgr -cut -mol -nopbc -nskip -nlevels -ndf -rgblo -rgbhi' -- $c)); return 0; fi
case "$p" in
-tu) COMPREPLY=( $(compgen -W ' ps fs ns us ms s m h ' -- $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
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
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f1 -f2 -o -n1 -n2 -no -h -nice -w -one -nomw -pbc -nofit -name -bfac' -- $c)); return 0; fi
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
-tu) COMPREPLY=( $(compgen -W ' ps fs ns us ms s m h ' -- $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa|gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-v) COMPREPLY=( $(compgen -X '!*.+(trr|trj)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-av) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|esp)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-l) COMPREPLY=( $(compgen -X '!*.log*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ascii) COMPREPLY=( $(compgen -X '!*.dat*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-xpm) COMPREPLY=( $(compgen -X '!*.xpm*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-xpma) COMPREPLY=( $(compgen -X '!*.xpm*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_covar_compl g_covar
shopt -s extglob
_g_density_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -n -s -ei -o -h -nice -b -e -dt -w -noxvgr -d -sl -number -ed -count -ng -symm -center' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
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
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -o -h -nice -b -e -dt -w -bin -nx -nz -amax -rmax -mirror -dmax' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
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
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.out*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_dih_compl g_dih
shopt -s extglob
_g_dipoles_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -enx -f -s -n -o -eps -a -d -c -g -adip -dip3d -cos -q -slab -h -nice -b -e -dt -w -noxvgr -mu -mumax -epsilonRF -skip -temp -avercorr -nopairs -axis -sl -gkratom -acflen -nonormalize -P -fitfn -ncskip -beginfit -endfit' -- $c)); return 0; fi
case "$p" in
-P) COMPREPLY=( $(compgen -W ' 0 1 2 3 ' -- $c ));;
-fitfn) COMPREPLY=( $(compgen -W ' none exp aexp exp_exp vac exp5 exp7 exp9 ' -- $c ));;
-enx) COMPREPLY=( $(compgen -X '!*.+(edr|ene)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
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
-q) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-slab) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_dipoles_compl g_dipoles
shopt -s extglob
_g_disre_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -s -f -ds -da -dn -dm -dr -l -n -q -c -h -nice -b -e -dt -w -noxvgr -ntop' -- $c)); return 0; fi
case "$p" in
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ds) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-da) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-dn) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-dm) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-dr) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-l) COMPREPLY=( $(compgen -X '!*.log*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-q) COMPREPLY=( $(compgen -X '!*.pdb*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-c) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_disre_compl g_disre
shopt -s extglob
_g_dist_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -o -h -nice -b -e -dt -noxvgr -dist' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
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
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -f2 -s -o -viol -pairs -ora -ort -oda -odr -odt -oten -corr -vis -ravg -h -nice -b -e -w -noxvgr -fee -fetemp -zero -sum -dp -mutot -skip -aver -nmol -ndf -fluc -orinst -ovec -acflen -nonormalize -P -fitfn -ncskip -beginfit -endfit' -- $c)); return 0; fi
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
_g_filter_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -ol -oh -h -nice -b -e -dt -w -nf -all -nonojump -fit' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
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
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -o -acf -h -nice -b -e -dt -w -noxvgr -nmol -q -p -moi -acflen -nonormalize -P -fitfn -ncskip -beginfit -endfit' -- $c)); return 0; fi
case "$p" in
-P) COMPREPLY=( $(compgen -W ' 0 1 2 3 ' -- $c ));;
-fitfn) COMPREPLY=( $(compgen -W ' none exp aexp exp_exp vac exp5 exp7 exp9 ' -- $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
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
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
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
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -g -sel -num -ac -dist -ang -hx -hbn -hbm -don -dan -life -h -nice -b -e -dt -noxvgr -ins -a -r -noda -abin -rbin -nonitacc -contact -shell -fitstart -temp -dump -max_hb -nomerge -acflen -nonormalize -P -fitfn -ncskip -beginfit -endfit' -- $c)); return 0; fi
case "$p" in
-P) COMPREPLY=( $(compgen -W ' 0 1 2 3 ' -- $c ));;
-fitfn) COMPREPLY=( $(compgen -W ' none exp aexp exp_exp vac exp5 exp7 exp9 ' -- $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-g) COMPREPLY=( $(compgen -X '!*.log*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-sel) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-num) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ac) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-dist) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ang) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-hx) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-hbn) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-hbm) COMPREPLY=( $(compgen -X '!*.xpm*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-don) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-dan) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-life) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
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
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-to) COMPREPLY=( $(compgen -X '!*.g87*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-cz) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|esp)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-co) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|esp)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_helix_compl g_helix
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
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
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
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -od -on -o -ox -or -h -nice -b -e -dt -tu -w -noxvgr -matrix -max -d -pi -split -ng' -- $c)); return 0; fi
case "$p" in
-tu) COMPREPLY=( $(compgen -W ' ps fs ns us ms s m h ' -- $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
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
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -o -mol -h -nice -b -e -dt -tu -w -noxvgr -type -lateral -ngroup -nomw -trestart -beginfit -endfit' -- $c)); return 0; fi
case "$p" in
-tu) COMPREPLY=( $(compgen -W ' ps fs ns us ms s m h ' -- $c ));;
-type) COMPREPLY=( $(compgen -W ' no x y z ' -- $c ));;
-lateral) COMPREPLY=( $(compgen -W ' no x y z ' -- $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa|gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-mol) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_msd_compl g_msd
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
-v) COMPREPLY=( $(compgen -X '!*.+(trr|trj)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_nmeig_compl g_nmeig
shopt -s extglob
_g_nmens_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -v -e -s -n -o -h -nice -noxvgr -temp -seed -num -first -last' -- $c)); return 0; fi
case "$p" in
-v) COMPREPLY=( $(compgen -X '!*.+(trr|trj)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
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
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -s -v -o -h -nice -eignr -temp -amplitude -nframes' -- $c)); return 0; fi
case "$p" in
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa|gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-v) COMPREPLY=( $(compgen -X '!*.+(trr|trj)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_nmtraj_compl g_nmtraj
shopt -s extglob
_g_order_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -n -s -o -od -os -h -nice -b -e -dt -w -noxvgr -d -sl -szonly -unsat' -- $c)); return 0; fi
case "$p" in
-d) COMPREPLY=( $(compgen -W ' z x y ' -- $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-od) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-os) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_order_compl g_order
shopt -s extglob
_g_potential_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -n -s -o -oc -of -h -nice -b -e -dt -w -noxvgr -d -sl -cb -ce -tz -spherical -ng' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-oc) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-of) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_potential_compl g_potential
shopt -s extglob
_g_rama_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -o -h -nice -b -e -dt -w -noxvgr' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_rama_compl g_rama
shopt -s extglob
_g_rdf_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -o -sq -cn -hq -h -nice -b -e -dt -w -noxvgr -bin -com -nopbc -xy -cut -ng -fade -nlevel -startq -endq -energy' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
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
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -s -f -f2 -n -o -mir -a -dist -m -bin -bm -h -nice -b -e -dt -tu -w -noxvgr -what -nopbc -fit -prev -split -skip -skip2 -max -min -bmax -bmin -nlevels -ng' -- $c)); return 0; fi
case "$p" in
-tu) COMPREPLY=( $(compgen -W ' ps fs ns us ms s m h ' -- $c ));;
-what) COMPREPLY=( $(compgen -W ' rmsd rho rhosc ' -- $c ));;
-fit) COMPREPLY=( $(compgen -W ' rot+trans translation none ' -- $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa|gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-f2) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
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
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
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
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -q -oq -ox -o -od -oc -dir -h -nice -b -e -dt -w -noxvgr -res -aniso' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
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
_g_rotacf_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -o -h -nice -b -e -dt -w -noxvgr -d -noaver -acflen -nonormalize -P -fitfn -ncskip -beginfit -endfit' -- $c)); return 0; fi
case "$p" in
-P) COMPREPLY=( $(compgen -W ' 0 1 2 3 ' -- $c ));;
-fitfn) COMPREPLY=( $(compgen -W ' none exp aexp exp_exp vac exp5 exp7 exp9 ' -- $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
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
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_saltbr_compl g_saltbr
shopt -s extglob
_g_sas_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -o -or -oa -q -n -i -h -nice -b -e -dt -w -noxvgr -solsize -ndots -qmax -f_index -minarea -nopbc -noprot -dgs' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-or) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-oa) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-q) COMPREPLY=( $(compgen -X '!*.pdb*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-i) COMPREPLY=( $(compgen -X '!*.itp*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_sas_compl g_sas
shopt -s extglob
_g_sgangle_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -n -s -oa -od -od1 -od2 -h -nice -b -e -dt -w -noxvgr -pbc -one -z' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
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
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -ge -ene -dist -histo -bin -ls -lsh -lss -map -ls3 -mdata -g -h -nice -w -noxvgr -notime -b -e -ttol -n -d -bw -nosham -tsham -pmin -dim -ngrid -xmin -xmax -gmax -nlevels -mname' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ge) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ene) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-dist) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-histo) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-bin) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
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
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -o -no -ro -co -h -nice -b -e -dt -w -noxvgr -com -rmin -rmax -bin -pbc' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa|gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-no) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ro) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-co) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_sorient_compl g_sorient
shopt -s extglob
_g_tcaf_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -ot -oa -o -of -oc -ov -h -nice -b -e -dt -w -noxvgr -mol -k34 -wt -acflen -nonormalize -P -fitfn -ncskip -beginfit -endfit' -- $c)); return 0; fi
case "$p" in
-P) COMPREPLY=( $(compgen -W ' 0 1 2 3 ' -- $c ));;
-fitfn) COMPREPLY=( $(compgen -W ' none exp aexp exp_exp vac exp5 exp7 exp9 ' -- $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(trr|trj)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
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
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -ox -ov -of -ob -ot -ekt -ekr -vd -cv -cf -h -nice -b -e -dt -tu -w -noxvgr -com -mol -nojump -nox -noy -noz -ng -len -bin -scale' -- $c)); return 0; fi
case "$p" in
-tu) COMPREPLY=( $(compgen -W ' ps fs ns us ms s m h ' -- $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa|gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ox) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ov) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-of) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ob) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ot) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ekt) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ekr) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-vd) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-cv) COMPREPLY=( $(compgen -X '!*.pdb*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-cf) COMPREPLY=( $(compgen -X '!*.pdb*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_traj_compl g_traj
shopt -s extglob
_g_velacc_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -o -h -nice -b -e -dt -w -noxvgr -mol -acflen -nonormalize -P -fitfn -ncskip -beginfit -endfit' -- $c)); return 0; fi
case "$p" in
-P) COMPREPLY=( $(compgen -W ' 0 1 2 3 ' -- $c ));;
-fitfn) COMPREPLY=( $(compgen -W ' none exp aexp exp_exp vac exp5 exp7 exp9 ' -- $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(trr|trj)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa|gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_velacc_compl g_velacc
shopt -s extglob
_g_wham_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -o -hist -h -nice -w -noxvgr -min -max -bins -noprof -temp -flip -tol' -- $c)); return 0; fi
case "$p" in
-o) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-hist) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_wham_compl g_wham
shopt -s extglob
_genbox_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -cp -cs -ci -o -p -h -nice -box -nmol -try -seed -vdwd -shell' -- $c)); return 0; fi
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
-trj) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _genconf_compl genconf
shopt -s extglob
_genion_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -s -table -n -o -g -pot -h -nice -noxvgr -np -pname -pq -nn -nname -nq -rmin -random -seed -scale' -- $c)); return 0; fi
case "$p" in
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-table) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|esp)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-g) COMPREPLY=( $(compgen -X '!*.log*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-pot) COMPREPLY=( $(compgen -X '!*.pdb*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _genion_compl genion
shopt -s extglob
_genpr_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -n -o -of -h -nice -fc -freeze' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|esp|tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.itp*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-of) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _genpr_compl genpr
shopt -s extglob
_gmxcheck_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -f2 -s1 -s2 -c -e -e2 -n -h -nice -vdwfac -bonlo -bonhi -tol -lastener' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-f2) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s1) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s2) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-c) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa|gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-e) COMPREPLY=( $(compgen -X '!*.+(edr|ene)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-e2) COMPREPLY=( $(compgen -X '!*.+(edr|ene)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _gmxcheck_compl gmxcheck
shopt -s extglob
_gmxdump_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -s -f -e -h -nice -nonr' -- $c)); return 0; fi
case "$p" in
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-e) COMPREPLY=( $(compgen -X '!*.+(edr|ene)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _gmxdump_compl gmxdump
shopt -s extglob
_grompp_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -po -c -r -rb -n -deshuf -p -pp -o -t -e -h -nice -nov -time -np -shuffle -sort -normvsbds -load -maxwarn -check14 -norenum' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.mdp*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-po) COMPREPLY=( $(compgen -X '!*.mdp*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-c) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|esp|tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-r) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|esp|tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-rb) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|esp|tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-deshuf) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-p) COMPREPLY=( $(compgen -X '!*.top*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-pp) COMPREPLY=( $(compgen -X '!*.top*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-t) COMPREPLY=( $(compgen -X '!*.+(trr|trj)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-e) COMPREPLY=( $(compgen -X '!*.+(edr|ene)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _grompp_compl grompp
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
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -eig -s -n -tar -ori -o -h -nice -noxvgr -mon -linfix -linacc -radfix -radacc -radcon -flood -outfrq -logfrq -slope -maxedsteps -deltaF0 -deltaF -tau -eqsteps -Eflnull -T -alpha -linstep -accdir -radstep -restrain -hesse -harmonic' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.+(trr|trj)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
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
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -s -o -x -c -e -g -dgdl -field -table -tablep -rerun -tpi -ei -eo -j -jo -ffout -devout -runav -pi -po -pd -pn -mtx -dn -h -nice -deffnm -noxvgr -np -nt -v -nocompact -sepdvdl -multi -replex -reseed -glas -ionize' -- $c)); return 0; fi
case "$p" in
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.+(trr|trj)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-x) COMPREPLY=( $(compgen -X '!*.xtc*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-c) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|esp)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-e) COMPREPLY=( $(compgen -X '!*.+(edr|ene)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-g) COMPREPLY=( $(compgen -X '!*.log*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-dgdl) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-field) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-table) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-tablep) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-rerun) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-tpi) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ei) COMPREPLY=( $(compgen -X '!*.edi*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-eo) COMPREPLY=( $(compgen -X '!*.edo*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-j) COMPREPLY=( $(compgen -X '!*.gct*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-jo) COMPREPLY=( $(compgen -X '!*.gct*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ffout) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-devout) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-runav) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-pi) COMPREPLY=( $(compgen -X '!*.ppa*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-po) COMPREPLY=( $(compgen -X '!*.ppa*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-pd) COMPREPLY=( $(compgen -X '!*.pdo*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-pn) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-mtx) COMPREPLY=( $(compgen -X '!*.mtx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-dn) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _mdrun_compl mdrun
shopt -s extglob
_mk_angndx_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -s -n -h -nice -type' -- $c)); return 0; fi
case "$p" in
-type) COMPREPLY=( $(compgen -W ' angle g96-angle dihedral improper ryckaert-bellemans phi-psi ' -- $c ));;
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
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _ngmx_compl ngmx
shopt -s extglob
_pdb2gmx_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -o -p -i -n -q -h -nice -merge -ff -water -inter -ss -ter -lys -asp -glu -his -angle -dist -una -ignh -missing -v -posrefc -vsite -heavyh -deuterate' -- $c)); return 0; fi
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
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _protonate_compl protonate
shopt -s extglob
_tpbconv_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -s -f -e -n -o -h -nice -time -extend -until -zeroq -nounconstrained' -- $c)); return 0; fi
case "$p" in
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(trr|trj)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
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
-tu) COMPREPLY=( $(compgen -W ' ps fs ns us ms s m h ' -- $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-demux) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _trjcat_compl trjcat
shopt -s extglob
_trjconv_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -o -s -n -fr -sub -drop -h -nice -b -e -tu -w -noxvgr -skip -dt -dump -t0 -timestep -pbc -ur -center -box -shift -fit -ndec -novel -force -trunc -exec -app -split -sep -ter -dropunder -dropover' -- $c)); return 0; fi
case "$p" in
-tu) COMPREPLY=( $(compgen -W ' ps fs ns us ms s m h ' -- $c ));;
-pbc) COMPREPLY=( $(compgen -W ' none whole inbox nojump cluster com ' -- $c ));;
-ur) COMPREPLY=( $(compgen -W ' rect tric compact ' -- $c ));;
-center) COMPREPLY=( $(compgen -W ' no tric rect zero ' -- $c ));;
-fit) COMPREPLY=( $(compgen -W ' none rot+trans translation progressive ' -- $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
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
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
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
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -o -r -h -nice -scale -ff -nexcl -noH14 -alldih -remdih -nopairs -name -nopbc -param -noround -kb -kt -kp' -- $c)); return 0; fi
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
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -f2 -di -do -o -xpm -h -nice -w -noframe -title -yonce -legend -diag -combine -size -bx -by -rainbow -gradient -skip -zeroline -legoffset' -- $c)); return 0; fi
case "$p" in
-title) COMPREPLY=( $(compgen -W ' top once ylabel none ' -- $c ));;
-legend) COMPREPLY=( $(compgen -W ' both first second none ' -- $c ));;
-diag) COMPREPLY=( $(compgen -W ' first second none ' -- $c ));;
-combine) COMPREPLY=( $(compgen -W ' halves add sub mult div ' -- $c ));;
-rainbow) COMPREPLY=( $(compgen -W ' no blue red ' -- $c ));;
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
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _xrama_compl xrama
