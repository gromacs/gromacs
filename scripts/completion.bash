shopt -s extglob
_do_dssp_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -ssdump -map -o -sc -a -ta -aa -h -version -nice -b -e -dt -tu -w -xvg -sss -ver' -- $c)); return 0; fi
case "$p" in
-tu) COMPREPLY=( $(compgen -W ' fs ps ns us ms s ' -- $c ));;
-xvg) COMPREPLY=( $(compgen -W ' xmgrace xmgr none ' -- $c ));;
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
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -n -o -mead -bf -h -version -nice -w -ndef -bt -box -angles -d -c -center -aligncenter -align -translate -rotate -princ -scale -density -pbc -resnr -grasp -rvdw -sig56 -vdwread -atom -legend -label -conect' -- $c)); return 0; fi
case "$p" in
-bt) COMPREPLY=( $(compgen -W ' triclinic cubic dodecahedron octahedron ' -- $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|esp|xyz|tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|esp|xyz)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-mead) COMPREPLY=( $(compgen -X '!*.pqr*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-bf) COMPREPLY=( $(compgen -X '!*.dat*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _editconf_compl editconf
shopt -s extglob
_eneconv_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -o -h -version -nice -b -e -dt -offset -settime -nosort -rmdh -scalefac -noerror' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.edr*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.edr*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _eneconv_compl eneconv
shopt -s extglob
_g_anadock_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -ox -od -of -g -h -version -nice -xvg -free -norms -cutoff' -- $c)); return 0; fi
case "$p" in
-xvg) COMPREPLY=( $(compgen -W ' xmgrace xmgr none ' -- $c ));;
-f) COMPREPLY=( $(compgen -X '!*.pdb*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ox) COMPREPLY=( $(compgen -X '!*.pdb*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-od) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-of) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-g) COMPREPLY=( $(compgen -X '!*.log*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_anadock_compl g_anadock
shopt -s extglob
_g_anaeig_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -v -v2 -f -s -n -eig -eig2 -comp -rmsf -proj -2d -3d -filt -extr -over -inpr -h -version -nice -b -e -dt -tu -w -xvg -first -last -skip -max -nframes -split -entropy -temp -nevskip' -- $c)); return 0; fi
case "$p" in
-tu) COMPREPLY=( $(compgen -W ' fs ps ns us ms s ' -- $c ));;
-xvg) COMPREPLY=( $(compgen -W ' xmgrace xmgr none ' -- $c ));;
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
-3d) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|esp|xyz)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
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
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -ac -msd -cc -dist -av -ee -bal -g -h -version -nice -w -xvg -notime -b -e -n -d -bw -errbar -integrate -aver_start -xydy -regression -luzar -temp -fitstart -fitend -smooth -filter -power -nosubav -oneacf -acflen -nonormalize -P -fitfn -ncskip -beginfit -endfit' -- $c)); return 0; fi
case "$p" in
-xvg) COMPREPLY=( $(compgen -W ' xmgrace xmgr none ' -- $c ));;
-errbar) COMPREPLY=( $(compgen -W ' none stddev error 90 ' -- $c ));;
-P) COMPREPLY=( $(compgen -W ' 0 1 2 3 ' -- $c ));;
-fitfn) COMPREPLY=( $(compgen -W ' none exp aexp exp_exp vac exp5 exp7 exp9 erffit ' -- $c ));;
-f) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ac) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-msd) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-cc) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-dist) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-av) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ee) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-bal) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-g) COMPREPLY=( $(compgen -X '!*.log*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_analyze_compl g_analyze
shopt -s extglob
_g_angle_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -n -od -ov -of -ot -oh -oc -or -h -version -nice -b -e -dt -w -xvg -type -all -binwidth -noperiodic -chandler -avercorr -acflen -nonormalize -P -fitfn -ncskip -beginfit -endfit' -- $c)); return 0; fi
case "$p" in
-xvg) COMPREPLY=( $(compgen -W ' xmgrace xmgr none ' -- $c ));;
-type) COMPREPLY=( $(compgen -W ' angle dihedral improper ryckaert-bellemans ' -- $c ));;
-P) COMPREPLY=( $(compgen -W ' 0 1 2 3 ' -- $c ));;
-fitfn) COMPREPLY=( $(compgen -W ' none exp aexp exp_exp vac exp5 exp7 exp9 erffit ' -- $c ));;
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
_g_bar_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -g -o -oi -oh -h -version -nice -w -xvg -b -e -temp -prec -nbmin -nbmax -nbin -extp' -- $c)); return 0; fi
case "$p" in
-xvg) COMPREPLY=( $(compgen -W ' xmgrace xmgr none ' -- $c ));;
-f) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-g) COMPREPLY=( $(compgen -X '!*.edr*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-oi) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-oh) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_bar_compl g_bar
shopt -s extglob
_g_bond_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -n -s -o -l -d -h -version -nice -b -e -dt -w -xvg -blen -tol -noaver -noaverdist' -- $c)); return 0; fi
case "$p" in
-xvg) COMPREPLY=( $(compgen -W ' xmgrace xmgr none ' -- $c ));;
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
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -ol -od -oz -ot -otr -otl -ok -okr -okl -oa -h -version -nice -b -e -dt -tu -xvg -na -z' -- $c)); return 0; fi
case "$p" in
-tu) COMPREPLY=( $(compgen -W ' fs ps ns us ms s ' -- $c ));;
-xvg) COMPREPLY=( $(compgen -W ' xmgrace xmgr none ' -- $c ));;
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
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -s -f -o -p -ss -jc -corr -g -ot -oh -rt -cp -h -version -nice -b -e -dt -w -xvg -r0 -phi -psi -omega -rama -viol -noperiodic -all -rad -shift -binwidth -core_rotamer -maxchi -nonormhisto -ramomega -bfact -chi_prod -HChi -bmax -acflen -nonormalize -P -fitfn -ncskip -beginfit -endfit' -- $c)); return 0; fi
case "$p" in
-xvg) COMPREPLY=( $(compgen -W ' xmgrace xmgr none ' -- $c ));;
-maxchi) COMPREPLY=( $(compgen -W ' 0 1 2 3 4 5 6 ' -- $c ));;
-P) COMPREPLY=( $(compgen -W ' 0 1 2 3 ' -- $c ));;
-fitfn) COMPREPLY=( $(compgen -W ' none exp aexp exp_exp vac exp5 exp7 exp9 erffit ' -- $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|esp|xyz|tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
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
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -dm -o -g -dist -ev -sz -tr -ntr -clid -cl -h -version -nice -b -e -dt -tu -w -xvg -dista -nlevels -cutoff -nofit -max -skip -av -wcl -nst -rmsmin -method -minstruct -binary -M -P -seed -niter -kT -nopbc' -- $c)); return 0; fi
case "$p" in
-tu) COMPREPLY=( $(compgen -W ' fs ps ns us ms s ' -- $c ));;
-xvg) COMPREPLY=( $(compgen -W ' xmgrace xmgr none ' -- $c ));;
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
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -o -ow -nc -mc -ac -hc -temp -mcn -h -version -nice -b -e -dt -tu -w -xvg -cut -mol -nopbc -nskip -nlevels -ndf -rgblo -rgbhi' -- $c)); return 0; fi
case "$p" in
-tu) COMPREPLY=( $(compgen -W ' fs ps ns us ms s ' -- $c ));;
-xvg) COMPREPLY=( $(compgen -W ' xmgrace xmgr none ' -- $c ));;
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
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f1 -f2 -o -n1 -n2 -no -h -version -nice -w -one -nomw -pbc -nofit -name -label -bfac' -- $c)); return 0; fi
case "$p" in
-f1) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa|gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-f2) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|esp|xyz|tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|esp|xyz)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n1) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n2) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-no) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_confrms_compl g_confrms
shopt -s extglob
_g_covar_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -o -v -av -l -ascii -xpm -xpma -h -version -nice -b -e -dt -tu -xvg -nofit -ref -mwa -last -nopbc' -- $c)); return 0; fi
case "$p" in
-tu) COMPREPLY=( $(compgen -W ' fs ps ns us ms s ' -- $c ));;
-xvg) COMPREPLY=( $(compgen -W ' xmgrace xmgr none ' -- $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa|gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-v) COMPREPLY=( $(compgen -X '!*.+(trr|cpt|trj)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-av) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|esp|xyz)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
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
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -s -n -f -o -caf -dsp -md -mj -mc -h -version -nice -b -e -dt -w -xvg -sh -nonojump -eps -bfit -efit -bvit -evit -tr -temp' -- $c)); return 0; fi
case "$p" in
-xvg) COMPREPLY=( $(compgen -W ' xmgrace xmgr none ' -- $c ));;
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
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -n -s -ei -o -h -version -nice -b -e -dt -w -xvg -d -sl -dens -ng -symm -center' -- $c)); return 0; fi
case "$p" in
-xvg) COMPREPLY=( $(compgen -W ' xmgrace xmgr none ' -- $c ));;
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
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -od -o -h -version -nice -b -e -dt -w -bin -aver -xmin -xmax -n1 -n2 -amax -rmax -mirror -sums -unit -dmin -dmax' -- $c)); return 0; fi
case "$p" in
-aver) COMPREPLY=( $(compgen -W ' z y x ' -- $c ));;
-unit) COMPREPLY=( $(compgen -W ' nm-3 nm-2 count ' -- $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa|gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-od) COMPREPLY=( $(compgen -X '!*.dat*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.xpm*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_densmap_compl g_densmap
shopt -s extglob
_g_densorder_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -s -f -n -o -or -og -Spect -h -version -nice -b -e -dt -w -1d -bw -bwn -order -axis -method -d1 -d2 -tblock -nlevel' -- $c)); return 0; fi
case "$p" in
-method) COMPREPLY=( $(compgen -W ' bisect functional ' -- $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.dat*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-or) COMPREPLY=( $(compgen -X '!*.out*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-og) COMPREPLY=( $(compgen -X '!*.xpm*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-Spect) COMPREPLY=( $(compgen -X '!*.out*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_densorder_compl g_densorder
shopt -s extglob
_g_dielectric_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -d -o -c -h -version -nice -b -e -dt -w -xvg -fft -nox1 -eint -bfit -efit -tail -A -tau1 -tau2 -eps0 -epsRF -fix -ffn -nsmooth' -- $c)); return 0; fi
case "$p" in
-xvg) COMPREPLY=( $(compgen -W ' xmgrace xmgr none ' -- $c ));;
-ffn) COMPREPLY=( $(compgen -W ' none exp aexp exp_exp vac exp5 exp7 exp9 erffit ' -- $c ));;
-f) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-d) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-c) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_dielectric_compl g_dielectric
shopt -s extglob
_g_dipoles_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -en -f -s -n -o -eps -a -d -c -g -adip -dip3d -cos -cmap -q -slab -h -version -nice -b -e -dt -w -xvg -mu -mumax -epsilonRF -skip -temp -corr -nopairs -ncos -axis -sl -gkratom -gkratom2 -rcmax -phi -nlevels -ndegrees -acflen -nonormalize -P -fitfn -ncskip -beginfit -endfit' -- $c)); return 0; fi
case "$p" in
-xvg) COMPREPLY=( $(compgen -W ' xmgrace xmgr none ' -- $c ));;
-corr) COMPREPLY=( $(compgen -W ' none mol molsep total ' -- $c ));;
-P) COMPREPLY=( $(compgen -W ' 0 1 2 3 ' -- $c ));;
-fitfn) COMPREPLY=( $(compgen -W ' none exp aexp exp_exp vac exp5 exp7 exp9 erffit ' -- $c ));;
-en) COMPREPLY=( $(compgen -X '!*.edr*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
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
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -s -f -ds -da -dn -dm -dr -l -n -q -c -x -h -version -nice -b -e -dt -w -xvg -ntop -maxdr -nlevels -nothird' -- $c)); return 0; fi
case "$p" in
-xvg) COMPREPLY=( $(compgen -W ' xmgrace xmgr none ' -- $c ));;
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
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -o -lt -h -version -nice -b -e -dt -xvg -intra -dist' -- $c)); return 0; fi
case "$p" in
-xvg) COMPREPLY=( $(compgen -W ' xmgrace xmgr none ' -- $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-lt) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_dist_compl g_dist
shopt -s extglob
_g_dos_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -vacf -mvacf -dos -g -h -version -nice -b -e -dt -w -xvg -nov -recip -abs -normdos -T -acflen -nonormalize -P -fitfn -ncskip -beginfit -endfit' -- $c)); return 0; fi
case "$p" in
-xvg) COMPREPLY=( $(compgen -W ' xmgrace xmgr none ' -- $c ));;
-P) COMPREPLY=( $(compgen -W ' 0 1 2 3 ' -- $c ));;
-fitfn) COMPREPLY=( $(compgen -W ' none exp aexp exp_exp vac exp5 exp7 exp9 erffit ' -- $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(trr|cpt|trj)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-vacf) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-mvacf) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-dos) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-g) COMPREPLY=( $(compgen -X '!*.log*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_dos_compl g_dos
shopt -s extglob
_g_dyecoupl_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -n -ot -oe -o -rhist -khist -h -version -nice -b -e -tu -w -xvg -pbcdist -norm -bins -R0' -- $c)); return 0; fi
case "$p" in
-tu) COMPREPLY=( $(compgen -W ' fs ps ns us ms s ' -- $c ));;
-xvg) COMPREPLY=( $(compgen -W ' xmgrace xmgr none ' -- $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ot) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-oe) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.dat*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-rhist) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-khist) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_dyecoupl_compl g_dyecoupl
shopt -s extglob
_g_dyndom_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -o -n -h -version -nice -firstangle -lastangle -nframe -maxangle -trans -head -tail' -- $c)); return 0; fi
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
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -cp -cs -ci -o -p -h -version -nice -box -nmol -try -seed -vdwd -shell -maxsol -vel' -- $c)); return 0; fi
case "$p" in
-cp) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|esp|xyz|tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-cs) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|esp|xyz|tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ci) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|esp|xyz|tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|esp|xyz)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-p) COMPREPLY=( $(compgen -X '!*.top*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _genbox_compl genbox
shopt -s extglob
_genconf_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -o -trj -h -version -nice -nbox -dist -seed -rot -shuffle -sort -block -nmolat -maxrot -norenumber' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|esp|xyz|tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|esp|xyz)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-trj) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _genconf_compl genconf
shopt -s extglob
_g_enemat_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -groups -eref -emat -etot -h -version -nice -b -e -dt -w -xvg -sum -skip -nomean -nlevels -max -min -nocoulsr -coullr -coul14 -noljsr -ljlr -lj14 -bhamsr -bhamlr -nofree -temp' -- $c)); return 0; fi
case "$p" in
-xvg) COMPREPLY=( $(compgen -W ' xmgrace xmgr none ' -- $c ));;
-f) COMPREPLY=( $(compgen -X '!*.edr*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
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
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -f2 -s -o -viol -pairs -ora -ort -oda -odr -odt -oten -corr -vis -ravg -odh -h -version -nice -b -e -w -xvg -fee -fetemp -zero -sum -dp -nbmin -nbmax -mutot -skip -aver -nmol -fluct_props -driftcorr -fluc -orinst -ovec -acflen -nonormalize -P -fitfn -ncskip -beginfit -endfit' -- $c)); return 0; fi
case "$p" in
-xvg) COMPREPLY=( $(compgen -W ' xmgrace xmgr none ' -- $c ));;
-P) COMPREPLY=( $(compgen -W ' 0 1 2 3 ' -- $c ));;
-fitfn) COMPREPLY=( $(compgen -W ' none exp aexp exp_exp vac exp5 exp7 exp9 erffit ' -- $c ));;
-f) COMPREPLY=( $(compgen -X '!*.edr*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-f2) COMPREPLY=( $(compgen -X '!*.edr*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
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
-odh) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_energy_compl g_energy
shopt -s extglob
_genion_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -s -table -n -o -g -pot -p -h -version -nice -xvg -np -pname -pq -nn -nname -nq -rmin -norandom -seed -scale -conc -neutral' -- $c)); return 0; fi
case "$p" in
-xvg) COMPREPLY=( $(compgen -W ' xmgrace xmgr none ' -- $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-table) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|esp|xyz)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-g) COMPREPLY=( $(compgen -X '!*.log*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-pot) COMPREPLY=( $(compgen -X '!*.pdb*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-p) COMPREPLY=( $(compgen -X '!*.top*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _genion_compl genion
shopt -s extglob
_genrestr_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -n -o -of -h -version -nice -fc -freeze -disre -disre_dist -disre_frac -disre_up2 -cutoff -constr' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|esp|xyz|tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.itp*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-of) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _genrestr_compl genrestr
shopt -s extglob
_g_filter_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -ol -oh -h -version -nice -b -e -dt -w -nf -all -nonojump -fit' -- $c)); return 0; fi
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
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -o -acf -h -version -nice -b -e -dt -w -xvg -nmol -q -p -moi -nz -acflen -nonormalize -P -fitfn -ncskip -beginfit -endfit' -- $c)); return 0; fi
case "$p" in
-xvg) COMPREPLY=( $(compgen -W ' xmgrace xmgr none ' -- $c ));;
-P) COMPREPLY=( $(compgen -W ' 0 1 2 3 ' -- $c ));;
-fitfn) COMPREPLY=( $(compgen -W ' none exp aexp exp_exp vac exp5 exp7 exp9 erffit ' -- $c ));;
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
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -n -nm -s -o -h -version -nice -b -e -dt -w -xvg -d -sl' -- $c)); return 0; fi
case "$p" in
-xvg) COMPREPLY=( $(compgen -W ' xmgrace xmgr none ' -- $c ));;
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
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -num -g -ac -dist -ang -hx -hbn -hbm -don -dan -life -nhbdist -h -version -nice -b -e -dt -tu -xvg -a -r -noda -r2 -abin -rbin -nonitacc -contact -shell -fitstart -fitstart -temp -smooth -dump -max_hb -nomerge -geminate -diff -nthreads -acflen -nonormalize -P -fitfn -ncskip -beginfit -endfit' -- $c)); return 0; fi
case "$p" in
-tu) COMPREPLY=( $(compgen -W ' fs ps ns us ms s ' -- $c ));;
-xvg) COMPREPLY=( $(compgen -W ' xmgrace xmgr none ' -- $c ));;
-geminate) COMPREPLY=( $(compgen -W ' none dd ad aa a4 ' -- $c ));;
-P) COMPREPLY=( $(compgen -W ' 0 1 2 3 ' -- $c ));;
-fitfn) COMPREPLY=( $(compgen -W ' none exp aexp exp_exp vac exp5 exp7 exp9 erffit ' -- $c ));;
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
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -s -n -f -to -cz -co -h -version -nice -b -e -dt -w -r0 -q -noF -db -prop -ev -ahxstart -ahxend' -- $c)); return 0; fi
case "$p" in
-prop) COMPREPLY=( $(compgen -W ' RAD TWIST RISE LEN NHX DIP RMS CPHI RMSA PHI PSI HB3 HB4 HB5 CD222 ' -- $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-to) COMPREPLY=( $(compgen -X '!*.g87*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-cz) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|esp|xyz)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-co) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|esp|xyz)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_helix_compl g_helix
shopt -s extglob
_g_helixorient_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -s -f -n -oaxis -ocenter -orise -oradius -otwist -obending -otilt -orot -h -version -nice -b -e -dt -xvg -sidechain -incremental' -- $c)); return 0; fi
case "$p" in
-xvg) COMPREPLY=( $(compgen -W ' xmgrace xmgr none ' -- $c ));;
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
_g_hydorder_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -n -s -o -or -Spect -h -version -nice -b -e -dt -w -d -bw -sgang1 -sgang2 -tblock -nlevel' -- $c)); return 0; fi
case "$p" in
-d) COMPREPLY=( $(compgen -W ' z x y ' -- $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.xpm*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-or) COMPREPLY=( $(compgen -X '!*.out*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-Spect) COMPREPLY=( $(compgen -X '!*.out*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_hydorder_compl g_hydorder
shopt -s extglob
_g_kinetics_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -d -d2 -o -o2 -o3 -ee -g -m -h -version -nice -tu -w -xvg -notime -b -e -bfit -efit -T -n -cut -ucut -euf -efu -ei -maxiter -noback -tol -skip -nosplit -nosum -nodiscrete -mult' -- $c)); return 0; fi
case "$p" in
-tu) COMPREPLY=( $(compgen -W ' fs ps ns us ms s ' -- $c ));;
-xvg) COMPREPLY=( $(compgen -W ' xmgrace xmgr none ' -- $c ));;
-f) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-d) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-d2) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o2) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o3) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ee) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-g) COMPREPLY=( $(compgen -X '!*.log*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-m) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_kinetics_compl g_kinetics
shopt -s extglob
_g_lie_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -o -h -version -nice -b -e -dt -w -xvg -Elj -Eqq -Clj -Cqq -ligand' -- $c)); return 0; fi
case "$p" in
-xvg) COMPREPLY=( $(compgen -W ' xmgrace xmgr none ' -- $c ));;
-f) COMPREPLY=( $(compgen -X '!*.edr*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_lie_compl g_lie
shopt -s extglob
_g_mdmat_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -mean -frames -no -h -version -nice -b -e -dt -xvg -t -nlevels' -- $c)); return 0; fi
case "$p" in
-xvg) COMPREPLY=( $(compgen -W ' xmgrace xmgr none ' -- $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa|gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-mean) COMPREPLY=( $(compgen -X '!*.xpm*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-frames) COMPREPLY=( $(compgen -X '!*.xpm*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-no) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_mdmat_compl g_mdmat
shopt -s extglob
_g_membed_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -n -p -o -x -c -e -dat -h -version -nice -xyinit -xyend -zinit -zend -nxy -nz -rad -pieces -asymmetry -ndiff -maxwarn -start -v -mdrun_path' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-p) COMPREPLY=( $(compgen -X '!*.top*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.+(trr|cpt|trj)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-x) COMPREPLY=( $(compgen -X '!*.xtc*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-c) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|esp|xyz)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-e) COMPREPLY=( $(compgen -X '!*.edr*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-dat) COMPREPLY=( $(compgen -X '!*.dat*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_membed_compl g_membed
shopt -s extglob
_g_mindist_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -od -on -o -ox -or -h -version -nice -b -e -dt -tu -w -xvg -matrix -max -d -group -pi -split -ng -nopbc -respertime -printresname' -- $c)); return 0; fi
case "$p" in
-tu) COMPREPLY=( $(compgen -W ' fs ps ns us ms s ' -- $c ));;
-xvg) COMPREPLY=( $(compgen -W ' xmgrace xmgr none ' -- $c ));;
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
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f1 -f2 -o -or -n -h -version -nice -w -xvg -ninterm -first -last -nofit' -- $c)); return 0; fi
case "$p" in
-xvg) COMPREPLY=( $(compgen -W ' xmgrace xmgr none ' -- $c ));;
-f1) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|esp|xyz|tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-f2) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|esp|xyz|tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-or) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_morph_compl g_morph
shopt -s extglob
_g_msd_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -o -mol -pdb -h -version -nice -b -e -tu -w -xvg -type -lateral -ten -ngroup -nomw -rmcomm -tpdb -trestart -beginfit -endfit' -- $c)); return 0; fi
case "$p" in
-tu) COMPREPLY=( $(compgen -W ' fs ps ns us ms s ' -- $c ));;
-xvg) COMPREPLY=( $(compgen -W ' xmgrace xmgr none ' -- $c ));;
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
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -f2 -s1 -s2 -c -e -e2 -n -m -h -version -nice -vdwfac -bonlo -bonhi -rmsd -tol -abstol -ab -lastener' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-f2) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s1) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s2) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-c) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa|gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-e) COMPREPLY=( $(compgen -X '!*.edr*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-e2) COMPREPLY=( $(compgen -X '!*.edr*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-m) COMPREPLY=( $(compgen -X '!*.tex*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _gmxcheck_compl gmxcheck
shopt -s extglob
_gmxdump_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -s -f -e -cp -p -mtx -om -h -version -nice -nonr -sys' -- $c)); return 0; fi
case "$p" in
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-e) COMPREPLY=( $(compgen -X '!*.edr*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-cp) COMPREPLY=( $(compgen -X '!*.cpt*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-p) COMPREPLY=( $(compgen -X '!*.top*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-mtx) COMPREPLY=( $(compgen -X '!*.mtx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-om) COMPREPLY=( $(compgen -X '!*.mdp*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _gmxdump_compl gmxdump
shopt -s extglob
_g_nmeig_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -of -ol -os -qc -v -h -version -nice -xvg -nom -first -last -maxspec -T -constr -width' -- $c)); return 0; fi
case "$p" in
-xvg) COMPREPLY=( $(compgen -W ' xmgrace xmgr none ' -- $c ));;
-f) COMPREPLY=( $(compgen -X '!*.mtx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-of) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ol) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-os) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-qc) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-v) COMPREPLY=( $(compgen -X '!*.+(trr|cpt|trj)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_nmeig_compl g_nmeig
shopt -s extglob
_g_nmens_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -v -e -s -n -o -h -version -nice -xvg -temp -seed -num -first -last' -- $c)); return 0; fi
case "$p" in
-xvg) COMPREPLY=( $(compgen -W ' xmgrace xmgr none ' -- $c ));;
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
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -s -v -o -h -version -nice -eignr -phases -temp -amplitude -nframes' -- $c)); return 0; fi
case "$p" in
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa|gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-v) COMPREPLY=( $(compgen -X '!*.+(trr|cpt|trj)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_nmtraj_compl g_nmtraj
shopt -s extglob
_g_options_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -h -version -nice' -- $c)); return 0; fi
case "$p" in
esac }
complete -F _g_options_compl g_options
shopt -s extglob
_g_order_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -n -nr -s -o -od -ob -os -Sg -Sk -Sgsl -Sksl -h -version -nice -b -e -dt -w -xvg -d -sl -szonly -unsat -permolecule -radial -calcdist' -- $c)); return 0; fi
case "$p" in
-xvg) COMPREPLY=( $(compgen -W ' xmgrace xmgr none ' -- $c ));;
-d) COMPREPLY=( $(compgen -W ' z x y ' -- $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-nr) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-od) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ob) COMPREPLY=( $(compgen -X '!*.pdb*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-os) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-Sg) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-Sk) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-Sgsl) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-Sksl) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_order_compl g_order
shopt -s extglob
_g_pme_error_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -s -o -so -h -version -nice -beta -tune -self -seed -v' -- $c)); return 0; fi
case "$p" in
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.out*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-so) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_pme_error_compl g_pme_error
shopt -s extglob
_g_polystat_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -s -f -n -o -v -p -i -h -version -nice -b -e -dt -tu -w -xvg -nomw -pc' -- $c)); return 0; fi
case "$p" in
-tu) COMPREPLY=( $(compgen -W ' fs ps ns us ms s ' -- $c ));;
-xvg) COMPREPLY=( $(compgen -W ' xmgrace xmgr none ' -- $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-v) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-p) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-i) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_polystat_compl g_polystat
shopt -s extglob
_g_potential_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -n -s -o -oc -of -h -version -nice -b -e -dt -w -xvg -d -sl -cb -ce -tz -spherical -ng -correct' -- $c)); return 0; fi
case "$p" in
-xvg) COMPREPLY=( $(compgen -W ' xmgrace xmgr none ' -- $c ));;
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
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -a1 -a2 -a3 -om -h -version -nice -b -e -dt -tu -w -foo' -- $c)); return 0; fi
case "$p" in
-tu) COMPREPLY=( $(compgen -W ' fs ps ns us ms s ' -- $c ));;
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
_g_protonate_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -s -f -n -o -h -version -nice -b -e -dt' -- $c)); return 0; fi
case "$p" in
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa|gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_protonate_compl g_protonate
shopt -s extglob
_g_rama_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -o -h -version -nice -b -e -dt -w -xvg' -- $c)); return 0; fi
case "$p" in
-xvg) COMPREPLY=( $(compgen -W ' xmgrace xmgr none ' -- $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_rama_compl g_rama
shopt -s extglob
_g_rdf_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -d -o -sq -cn -hq -h -version -nice -b -e -dt -w -xvg -bin -com -surf -rdf -nopbc -nonorm -xy -cut -ng -fade -nlevel -startq -endq -energy' -- $c)); return 0; fi
case "$p" in
-xvg) COMPREPLY=( $(compgen -W ' xmgrace xmgr none ' -- $c ));;
-surf) COMPREPLY=( $(compgen -W ' no mol res ' -- $c ));;
-rdf) COMPREPLY=( $(compgen -W ' atom mol_com mol_cog res_com res_cog ' -- $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa|gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-d) COMPREPLY=( $(compgen -X '!*.dat*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
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
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -s -f -f2 -n -o -mir -a -dist -m -bin -bm -h -version -nice -b -e -dt -tu -w -xvg -what -nopbc -fit -prev -split -skip -skip2 -max -min -bmax -bmin -nomw -nlevels -ng' -- $c)); return 0; fi
case "$p" in
-tu) COMPREPLY=( $(compgen -W ' fs ps ns us ms s ' -- $c ));;
-xvg) COMPREPLY=( $(compgen -W ' xmgrace xmgr none ' -- $c ));;
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
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -equiv -o -rms -scl -mean -nmr3 -nmr6 -noe -h -version -nice -b -e -dt -w -xvg -nlevels -max -nosumh -nopbc' -- $c)); return 0; fi
case "$p" in
-xvg) COMPREPLY=( $(compgen -W ' xmgrace xmgr none ' -- $c ));;
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
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -q -oq -ox -o -od -oc -dir -h -version -nice -b -e -dt -w -xvg -res -aniso -nofit' -- $c)); return 0; fi
case "$p" in
-xvg) COMPREPLY=( $(compgen -W ' xmgrace xmgr none ' -- $c ));;
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
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -po -c -r -rb -n -p -pp -o -t -e -ref -h -version -nice -v -time -normvsbds -maxwarn -zero -norenum' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.mdp*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-po) COMPREPLY=( $(compgen -X '!*.mdp*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-c) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|esp|xyz|tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-r) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|esp|xyz|tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-rb) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|esp|xyz|tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-p) COMPREPLY=( $(compgen -X '!*.top*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-pp) COMPREPLY=( $(compgen -X '!*.top*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-t) COMPREPLY=( $(compgen -X '!*.+(trr|cpt|trj)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-e) COMPREPLY=( $(compgen -X '!*.edr*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ref) COMPREPLY=( $(compgen -X '!*.+(trr|cpt|trj)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _grompp_compl grompp
shopt -s extglob
_g_rotacf_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -o -h -version -nice -b -e -dt -w -xvg -d -noaver -acflen -nonormalize -P -fitfn -ncskip -beginfit -endfit' -- $c)); return 0; fi
case "$p" in
-xvg) COMPREPLY=( $(compgen -W ' xmgrace xmgr none ' -- $c ));;
-P) COMPREPLY=( $(compgen -W ' 0 1 2 3 ' -- $c ));;
-fitfn) COMPREPLY=( $(compgen -W ' none exp aexp exp_exp vac exp5 exp7 exp9 erffit ' -- $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_rotacf_compl g_rotacf
shopt -s extglob
_g_rotmat_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -o -h -version -nice -b -e -dt -w -xvg -ref -skip -fitxy -nomw' -- $c)); return 0; fi
case "$p" in
-xvg) COMPREPLY=( $(compgen -W ' xmgrace xmgr none ' -- $c ));;
-ref) COMPREPLY=( $(compgen -W ' none xyz xy ' -- $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa|gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_rotmat_compl g_rotmat
shopt -s extglob
_g_saltbr_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -h -version -nice -b -e -dt -t -sep' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_saltbr_compl g_saltbr
shopt -s extglob
_g_sans_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -s -f -n -d -pr -sq -prframe -sqframe -h -version -nice -b -e -dt -tu -xvg -mode -mcover -nopbc -startq -endq -qstep -seed -nt' -- $c)); return 0; fi
case "$p" in
-tu) COMPREPLY=( $(compgen -W ' fs ps ns us ms s ' -- $c ));;
-xvg) COMPREPLY=( $(compgen -W ' xmgrace xmgr none ' -- $c ));;
-mode) COMPREPLY=( $(compgen -W ' direct mc ' -- $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-d) COMPREPLY=( $(compgen -X '!*.dat*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-pr) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-sq) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-prframe) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-sqframe) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_sans_compl g_sans
shopt -s extglob
_g_sas_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -o -or -oa -tv -q -n -i -h -version -nice -b -e -dt -w -xvg -probe -ndots -qmax -f_index -minarea -nopbc -noprot -dgs' -- $c)); return 0; fi
case "$p" in
-xvg) COMPREPLY=( $(compgen -W ' xmgrace xmgr none ' -- $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa|gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
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
_g_select_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -sf -n -os -oc -oi -om -on -h -version -nice -b -e -dt -xvg -normpbc -nopbc -select -selrpos -seltype -dump -norm -cfnorm -resnr' -- $c)); return 0; fi
case "$p" in
-xvg) COMPREPLY=( $(compgen -W ' xmgrace xmgr none ' -- $c ));;
-selrpos) COMPREPLY=( $(compgen -W ' atom res_com res_cog mol_com mol_cog whole_res_com whole_res_cog whole_mol_com whole_mol_cog part_res_com part_res_cog part_mol_com part_mol_cog dyn_res_com dyn_res_cog dyn_mol_com dyn_mol_cog ' -- $c ));;
-seltype) COMPREPLY=( $(compgen -W ' atom res_com res_cog mol_com mol_cog whole_res_com whole_res_cog whole_mol_com whole_mol_cog part_res_com part_res_cog part_mol_com part_mol_cog dyn_res_com dyn_res_cog dyn_mol_com dyn_mol_cog ' -- $c ));;
-resnr) COMPREPLY=( $(compgen -W ' number index ' -- $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa|gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-sf) COMPREPLY=( $(compgen -X '!*.dat*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-os) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-oc) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-oi) COMPREPLY=( $(compgen -X '!*.dat*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-om) COMPREPLY=( $(compgen -X '!*.dat*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-on) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_select_compl g_select
shopt -s extglob
_g_sgangle_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -n -s -oa -od -od1 -od2 -h -version -nice -b -e -dt -w -xvg -one -z' -- $c)); return 0; fi
case "$p" in
-xvg) COMPREPLY=( $(compgen -W ' xmgrace xmgr none ' -- $c ));;
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
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -ge -ene -dist -histo -bin -lp -ls -lsh -lss -map -ls3 -mdata -g -h -version -nice -w -xvg -notime -b -e -ttol -n -d -bw -nosham -tsham -pmin -dim -ngrid -xmin -xmax -pmax -gmax -emin -emax -nlevels -mname' -- $c)); return 0; fi
case "$p" in
-xvg) COMPREPLY=( $(compgen -W ' xmgrace xmgr none ' -- $c ));;
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
_g_sigeps_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -o -h -version -nice -w -xvg -c6 -cn -pow -sig -eps -A -B -C -qi -qj -sigfac' -- $c)); return 0; fi
case "$p" in
-xvg) COMPREPLY=( $(compgen -W ' xmgrace xmgr none ' -- $c ));;
-o) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_sigeps_compl g_sigeps
shopt -s extglob
_g_sorient_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -o -no -ro -co -rc -h -version -nice -b -e -dt -w -xvg -com -v23 -rmin -rmax -cbin -rbin -pbc' -- $c)); return 0; fi
case "$p" in
-xvg) COMPREPLY=( $(compgen -W ' xmgrace xmgr none ' -- $c ));;
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
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -s -f -n -h -version -nice -b -e -dt -w -pbc -nodiv -ign -bin -nab' -- $c)); return 0; fi
case "$p" in
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa|gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_spatial_compl g_spatial
shopt -s extglob
_g_spol_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -o -h -version -nice -b -e -dt -w -xvg -com -refat -rmin -rmax -dip -bw' -- $c)); return 0; fi
case "$p" in
-xvg) COMPREPLY=( $(compgen -W ' xmgrace xmgr none ' -- $c ));;
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
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -ot -oa -o -of -oc -ov -h -version -nice -b -e -dt -w -xvg -mol -k34 -wt -acflen -nonormalize -P -fitfn -ncskip -beginfit -endfit' -- $c)); return 0; fi
case "$p" in
-xvg) COMPREPLY=( $(compgen -W ' xmgrace xmgr none ' -- $c ));;
-P) COMPREPLY=( $(compgen -W ' 0 1 2 3 ' -- $c ));;
-fitfn) COMPREPLY=( $(compgen -W ' none exp aexp exp_exp vac exp5 exp7 exp9 erffit ' -- $c ));;
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
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -ox -oxt -ov -of -ob -ot -ekt -ekr -vd -cv -cf -av -af -h -version -nice -b -e -dt -tu -w -xvg -com -nopbc -mol -nojump -nox -noy -noz -ng -len -fp -bin -ctime -scale' -- $c)); return 0; fi
case "$p" in
-tu) COMPREPLY=( $(compgen -W ' fs ps ns us ms s ' -- $c ));;
-xvg) COMPREPLY=( $(compgen -W ' xmgrace xmgr none ' -- $c ));;
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
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -p -err -so -s -o -x -cpi -cpo -c -e -g -dhdl -field -table -tabletf -tablep -tableb -rerun -tpi -tpid -ei -eo -j -jo -ffout -devout -runav -px -pf -ro -ra -rs -rt -mtx -dn -bo -bx -bcpo -bc -be -bg -beo -bdhdl -bfield -btpi -btpid -bjo -bffout -bdevout -brunav -bpx -bpf -bro -bra -brs -brt -bmtx -bdn -h -version -nice -xvg -np -npstring -ntmpi -r -max -min -npme -fix -rmax -rmin -noscalevdw -ntpr -steps -resetstep -simsteps -launch -nobench -noappend -cpnum' -- $c)); return 0; fi
case "$p" in
-xvg) COMPREPLY=( $(compgen -W ' xmgrace xmgr none ' -- $c ));;
-npstring) COMPREPLY=( $(compgen -W ' -np -n none ' -- $c ));;
-npme) COMPREPLY=( $(compgen -W ' auto all subset ' -- $c ));;
-p) COMPREPLY=( $(compgen -X '!*.out*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-err) COMPREPLY=( $(compgen -X '!*.log*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-so) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.+(trr|cpt|trj)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-x) COMPREPLY=( $(compgen -X '!*.xtc*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-cpi) COMPREPLY=( $(compgen -X '!*.cpt*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-cpo) COMPREPLY=( $(compgen -X '!*.cpt*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-c) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|esp|xyz)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-e) COMPREPLY=( $(compgen -X '!*.edr*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-g) COMPREPLY=( $(compgen -X '!*.log*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-dhdl) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-field) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-table) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-tabletf) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-tablep) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-tableb) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-rerun) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-tpi) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-tpid) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ei) COMPREPLY=( $(compgen -X '!*.edi*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-eo) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-j) COMPREPLY=( $(compgen -X '!*.gct*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-jo) COMPREPLY=( $(compgen -X '!*.gct*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ffout) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-devout) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-runav) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-px) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-pf) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ro) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ra) COMPREPLY=( $(compgen -X '!*.log*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-rs) COMPREPLY=( $(compgen -X '!*.log*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-rt) COMPREPLY=( $(compgen -X '!*.log*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-mtx) COMPREPLY=( $(compgen -X '!*.mtx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-dn) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-bo) COMPREPLY=( $(compgen -X '!*.+(trr|cpt|trj)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-bx) COMPREPLY=( $(compgen -X '!*.xtc*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-bcpo) COMPREPLY=( $(compgen -X '!*.cpt*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-bc) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|esp|xyz)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-be) COMPREPLY=( $(compgen -X '!*.edr*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-bg) COMPREPLY=( $(compgen -X '!*.log*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-beo) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-bdhdl) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-bfield) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-btpi) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-btpid) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-bjo) COMPREPLY=( $(compgen -X '!*.gct*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-bffout) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-bdevout) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-brunav) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-bpx) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-bpf) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-bro) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-bra) COMPREPLY=( $(compgen -X '!*.log*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-brs) COMPREPLY=( $(compgen -X '!*.log*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-brt) COMPREPLY=( $(compgen -X '!*.log*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-bmtx) COMPREPLY=( $(compgen -X '!*.mtx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-bdn) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_tune_pme_compl g_tune_pme
shopt -s extglob
_g_vanhove_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -om -or -ot -h -version -nice -b -e -dt -w -xvg -sqrt -fm -rmax -rbin -mmax -nlevels -nr -fr -rt -ft' -- $c)); return 0; fi
case "$p" in
-xvg) COMPREPLY=( $(compgen -W ' xmgrace xmgr none ' -- $c ));;
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
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -o -os -h -version -nice -b -e -dt -w -xvg -m -norecip -mol -acflen -nonormalize -P -fitfn -ncskip -beginfit -endfit' -- $c)); return 0; fi
case "$p" in
-xvg) COMPREPLY=( $(compgen -W ' xmgrace xmgr none ' -- $c ));;
-P) COMPREPLY=( $(compgen -W ' 0 1 2 3 ' -- $c ));;
-fitfn) COMPREPLY=( $(compgen -W ' none exp aexp exp_exp vac exp5 exp7 exp9 erffit ' -- $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(trr|cpt|trj)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa|gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-os) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_velacc_compl g_velacc
shopt -s extglob
_g_wham_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -ix -if -it -ip -o -hist -oiact -iiact -bsres -bsprof -tab -h -version -nice -xvg -min -max -noauto -bins -temp -tol -v -b -e -dt -histonly -boundsonly -nolog -unit -zprof0 -cycl -sym -ac -acsig -ac-trestart -nBootstrap -bs-method -bs-tau -bs-seed -histbs-block -vbs' -- $c)); return 0; fi
case "$p" in
-xvg) COMPREPLY=( $(compgen -W ' xmgrace xmgr none ' -- $c ));;
-unit) COMPREPLY=( $(compgen -W ' kJ kCal kT ' -- $c ));;
-bs-method) COMPREPLY=( $(compgen -W ' b-hist hist traj traj-gauss ' -- $c ));;
-ix) COMPREPLY=( $(compgen -X '!*.dat*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-if) COMPREPLY=( $(compgen -X '!*.dat*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-it) COMPREPLY=( $(compgen -X '!*.dat*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ip) COMPREPLY=( $(compgen -X '!*.dat*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-hist) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-oiact) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-iiact) COMPREPLY=( $(compgen -X '!*.dat*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-bsres) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-bsprof) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-tab) COMPREPLY=( $(compgen -X '!*.dat*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_wham_compl g_wham
shopt -s extglob
_g_wheel_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -o -h -version -nice -r0 -rot0 -T -nonn' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.dat*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.eps*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_wheel_compl g_wheel
shopt -s extglob
_g_x2top_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -o -r -h -version -nice -ff -v -nexcl -noH14 -alldih -remdih -nopairs -name -nopbc -pdbq -noparam -noround -kb -kt -kp' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|esp|xyz|tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.top*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-r) COMPREPLY=( $(compgen -X '!*.rtp*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_x2top_compl g_x2top
shopt -s extglob
_g_xrama_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -h -version -nice -b -e -dt' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _g_xrama_compl g_xrama
shopt -s extglob
_make_edi_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -eig -s -n -tar -ori -o -h -version -nice -xvg -mon -linfix -linacc -radfix -radacc -radcon -flood -outfrq -slope -linstep -accdir -radstep -maxedsteps -eqsteps -deltaF0 -deltaF -tau -Eflnull -T -alpha -restrain -hessian -harmonic -constF' -- $c)); return 0; fi
case "$p" in
-xvg) COMPREPLY=( $(compgen -W ' xmgrace xmgr none ' -- $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(trr|cpt|trj)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-eig) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa|gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-tar) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|esp|xyz|tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ori) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|esp|xyz|tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.edi*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _make_edi_compl make_edi
shopt -s extglob
_make_ndx_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -n -o -h -version -nice -natoms' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|esp|xyz|tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _make_ndx_compl make_ndx
shopt -s extglob
_mdrun_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -s -o -x -cpi -cpo -c -e -g -dhdl -field -table -tabletf -tablep -tableb -rerun -tpi -tpid -ei -eo -j -jo -ffout -devout -runav -px -pf -ro -ra -rs -rt -mtx -dn -multidir -membed -mp -mn -h -version -nice -deffnm -xvg -pd -dd -ddorder -npme -nt -ntmpi -ntomp -ntomp_pme -pin -pinoffset -pinstride -gpu_id -noddcheck -rdd -rcon -dlb -dds -gcom -nb -notunepme -testverlet -v -nocompact -seppot -pforce -reprod -cpt -cpnum -noappend -nsteps -maxh -multi -replex -nex -reseed -ionize' -- $c)); return 0; fi
case "$p" in
-xvg) COMPREPLY=( $(compgen -W ' xmgrace xmgr none ' -- $c ));;
-ddorder) COMPREPLY=( $(compgen -W ' interleave pp_pme cartesian ' -- $c ));;
-pin) COMPREPLY=( $(compgen -W ' auto on off ' -- $c ));;
-dlb) COMPREPLY=( $(compgen -W ' auto no yes ' -- $c ));;
-nb) COMPREPLY=( $(compgen -W ' auto cpu gpu gpu_cpu ' -- $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.+(trr|cpt|trj)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-x) COMPREPLY=( $(compgen -X '!*.xtc*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-cpi) COMPREPLY=( $(compgen -X '!*.cpt*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-cpo) COMPREPLY=( $(compgen -X '!*.cpt*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-c) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|esp|xyz)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-e) COMPREPLY=( $(compgen -X '!*.edr*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-g) COMPREPLY=( $(compgen -X '!*.log*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-dhdl) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-field) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-table) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-tabletf) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-tablep) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-tableb) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-rerun) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-tpi) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-tpid) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ei) COMPREPLY=( $(compgen -X '!*.edi*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-eo) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-j) COMPREPLY=( $(compgen -X '!*.gct*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-jo) COMPREPLY=( $(compgen -X '!*.gct*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ffout) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-devout) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-runav) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-px) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-pf) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ro) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ra) COMPREPLY=( $(compgen -X '!*.log*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-rs) COMPREPLY=( $(compgen -X '!*.log*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-rt) COMPREPLY=( $(compgen -X '!*.log*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-mtx) COMPREPLY=( $(compgen -X '!*.mtx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-dn) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-multidir) COMPREPLY=( $(compgen -X '!*.*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-membed) COMPREPLY=( $(compgen -X '!*.dat*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-mp) COMPREPLY=( $(compgen -X '!*.top*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-mn) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _mdrun_compl mdrun
shopt -s extglob
_mdrun_mpi_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -s -o -x -cpi -cpo -c -e -g -dhdl -field -table -tabletf -tablep -tableb -rerun -tpi -tpid -ei -eo -j -jo -ffout -devout -runav -px -pf -ro -ra -rs -rt -mtx -dn -multidir -membed -mp -mn -h -version -nice -deffnm -xvg -pd -dd -ddorder -npme -nt -ntmpi -ntomp -ntomp_pme -pin -pinoffset -pinstride -gpu_id -noddcheck -rdd -rcon -dlb -dds -gcom -nb -notunepme -testverlet -v -nocompact -seppot -pforce -reprod -cpt -cpnum -noappend -nsteps -maxh -multi -replex -nex -reseed -ionize' -- $c)); return 0; fi
case "$p" in
-xvg) COMPREPLY=( $(compgen -W ' xmgrace xmgr none ' -- $c ));;
-ddorder) COMPREPLY=( $(compgen -W ' interleave pp_pme cartesian ' -- $c ));;
-pin) COMPREPLY=( $(compgen -W ' auto on off ' -- $c ));;
-dlb) COMPREPLY=( $(compgen -W ' auto no yes ' -- $c ));;
-nb) COMPREPLY=( $(compgen -W ' auto cpu gpu gpu_cpu ' -- $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.+(trr|cpt|trj)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-x) COMPREPLY=( $(compgen -X '!*.xtc*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-cpi) COMPREPLY=( $(compgen -X '!*.cpt*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-cpo) COMPREPLY=( $(compgen -X '!*.cpt*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-c) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|esp|xyz)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-e) COMPREPLY=( $(compgen -X '!*.edr*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-g) COMPREPLY=( $(compgen -X '!*.log*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-dhdl) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-field) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-table) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-tabletf) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-tablep) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-tableb) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-rerun) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-tpi) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-tpid) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ei) COMPREPLY=( $(compgen -X '!*.edi*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-eo) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-j) COMPREPLY=( $(compgen -X '!*.gct*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-jo) COMPREPLY=( $(compgen -X '!*.gct*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ffout) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-devout) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-runav) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-px) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-pf) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ro) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ra) COMPREPLY=( $(compgen -X '!*.log*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-rs) COMPREPLY=( $(compgen -X '!*.log*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-rt) COMPREPLY=( $(compgen -X '!*.log*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-mtx) COMPREPLY=( $(compgen -X '!*.mtx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-dn) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-multidir) COMPREPLY=( $(compgen -X '!*.line_buf*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-membed) COMPREPLY=( $(compgen -X '!*.dat*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-mp) COMPREPLY=( $(compgen -X '!*.top*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-mn) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _mdrun_mpi_compl mdrun_mpi
shopt -s extglob
_mk_angndx_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -s -n -h -version -nice -type -nohyd -hq' -- $c)); return 0; fi
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
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -h -version -nice -b -e -dt' -- $c)); return 0; fi
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
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -o -p -i -n -q -h -version -nice -chainsep -merge -ff -water -inter -ss -ter -lys -arg -asp -glu -gln -his -angle -dist -una -ignh -missing -v -posrefc -vsite -heavyh -deuterate -nochargegrp -nocmap -renum -rtpres' -- $c)); return 0; fi
case "$p" in
-chainsep) COMPREPLY=( $(compgen -W ' id_or_ter id_and_ter ter id interactive ' -- $c ));;
-merge) COMPREPLY=( $(compgen -W ' no all interactive ' -- $c ));;
-water) COMPREPLY=( $(compgen -W ' select none spc spce tip3p tip4p tip5p ' -- $c ));;
-vsite) COMPREPLY=( $(compgen -W ' none hydrogens aromatics ' -- $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|esp|xyz|tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|esp|xyz)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-p) COMPREPLY=( $(compgen -X '!*.top*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-i) COMPREPLY=( $(compgen -X '!*.itp*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-q) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|esp|xyz)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _pdb2gmx_compl pdb2gmx
shopt -s extglob
_tpbconv_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -s -f -e -n -o -h -version -nice -extend -until -nsteps -time -zeroq -novel -nocont -init_fep_state' -- $c)); return 0; fi
case "$p" in
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(trr|cpt|trj)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-e) COMPREPLY=( $(compgen -X '!*.edr*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _tpbconv_compl tpbconv
shopt -s extglob
_trjcat_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -o -n -demux -h -version -nice -tu -xvg -b -e -dt -prec -novel -settime -nosort -keeplast -overwrite -cat' -- $c)); return 0; fi
case "$p" in
-tu) COMPREPLY=( $(compgen -W ' fs ps ns us ms s ' -- $c ));;
-xvg) COMPREPLY=( $(compgen -W ' xmgrace xmgr none ' -- $c ));;
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
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -o -s -n -fr -sub -drop -h -version -nice -b -e -tu -w -xvg -skip -dt -round -dump -t0 -timestep -pbc -ur -center -boxcenter -box -clustercenter -trans -shift -fit -ndec -novel -force -trunc -exec -app -split -sep -nzero -dropunder -dropover -conect' -- $c)); return 0; fi
case "$p" in
-tu) COMPREPLY=( $(compgen -W ' fs ps ns us ms s ' -- $c ));;
-xvg) COMPREPLY=( $(compgen -W ' xmgrace xmgr none ' -- $c ));;
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
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -o -nshell -h -version -nice -b -e -dt -xvg -na -da -com -r -z' -- $c)); return 0; fi
case "$p" in
-xvg) COMPREPLY=( $(compgen -W ' xmgrace xmgr none ' -- $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|cpt|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa|gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-nshell) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac }
complete -F _trjorder_compl trjorder
shopt -s extglob
_xpm2ps_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -f2 -di -do -o -xpm -h -version -nice -w -noframe -title -yonce -legend -diag -size -bx -by -rainbow -gradient -skip -zeroline -legoffset -combine -cmin -cmax' -- $c)); return 0; fi
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
