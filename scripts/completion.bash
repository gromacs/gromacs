shopt -s extglob
_do_dssp_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -ssdump -map -o -sc -a -ta -aa -h -X -nice -b -e -dt -tu -w -sss' -- $c)); return 0; fi
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
esac; }; 
complete -F _do_dssp_compl do_dssp
shopt -s extglob
_editconf_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -n -o -bf -h -X -nice -w -ndef -bt -box -angles -d -c -center -rotate -princ -scale -density -pbc -mead -grasp -rvdw -atom -legend -label' -- $c)); return 0; fi
case "$p" in
-bt) COMPREPLY=( $(compgen -W ' tric cubic dodecahedron octahedron ' -- $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-bf) COMPREPLY=( $(compgen -X '!*.dat*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac; }; 
complete -F _editconf_compl editconf
shopt -s extglob
_eneconv_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -o -h -X -nice -b -e -dt -offset -settime -nosort -scalefac -noerror' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.+(edr|ene)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.+(edr|ene)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac; }; 
complete -F _eneconv_compl eneconv
shopt -s extglob
_ffscan_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -s -g -table -parm -ga -h -X -nice -tol -fmax -epot -nocomb -npow -ratio -logeps -nov' -- $c)); return 0; fi
case "$p" in
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-g) COMPREPLY=( $(compgen -X '!*.log*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-table) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-parm) COMPREPLY=( $(compgen -X '!*.dat*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ga) COMPREPLY=( $(compgen -X '!*.dat*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac; }; 
complete -F _ffscan_compl ffscan
shopt -s extglob
_g_anaeig_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -v -v2 -f -s -n -eig1 -eig2 -disp -proj -2d -3d -filt -extr -over -inpr -h -X -nice -b -e -dt -tu -w -first -last -skip -max -nframes -split' -- $c)); return 0; fi
case "$p" in
-tu) COMPREPLY=( $(compgen -W ' ps fs ns us ms s m h ' -- $c ));;
-v) COMPREPLY=( $(compgen -X '!*.+(trr|trj)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-v2) COMPREPLY=( $(compgen -X '!*.+(trr|trj)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa|gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-eig1) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-eig2) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-disp) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-proj) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-2d) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-3d) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-filt) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-extr) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-over) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-inpr) COMPREPLY=( $(compgen -X '!*.xpm*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac; }; 
complete -F _g_anaeig_compl g_anaeig
shopt -s extglob
_g_analyze_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -ac -msd -cc -dist -av -ee -g -h -X -nice -w -notime -b -e -n -d -bw -errbar -power -nosubav -oneacf -acflen -nonormalize -P -fitfn -ncskip -beginfit -endfit' -- $c)); return 0; fi
case "$p" in
-errbar) COMPREPLY=( $(compgen -W ' none stddev error 90 ' -- $c ));;
-P) COMPREPLY=( $(compgen -W ' 0 1 2 3 ' -- $c ));;
-fitfn) COMPREPLY=( $(compgen -W ' none exp aexp exp_exp vac exp5 exp7 ' -- $c ));;
-f) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ac) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-msd) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-cc) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-dist) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-av) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ee) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-g) COMPREPLY=( $(compgen -X '!*.log*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac; }; 
complete -F _g_analyze_compl g_analyze
shopt -s extglob
_g_angle_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -od -ov -of -ot -oh -oc -h -X -nice -b -e -dt -w -type -all -binwidth -chandler -avercorr -acflen -nonormalize -P -fitfn -ncskip -beginfit -endfit' -- $c)); return 0; fi
case "$p" in
-type) COMPREPLY=( $(compgen -W ' angle dihedral improper ryckaert-bellemans ' -- $c ));;
-P) COMPREPLY=( $(compgen -W ' 0 1 2 3 ' -- $c ));;
-fitfn) COMPREPLY=( $(compgen -W ' none exp aexp exp_exp vac exp5 exp7 ' -- $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-od) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ov) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-of) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ot) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-oh) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-oc) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac; }; 
complete -F _g_angle_compl g_angle
shopt -s extglob
_g_bond_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -n -s -o -l -d -h -X -nice -b -e -dt -w -blen -tol -noaver' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa|gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-l) COMPREPLY=( $(compgen -X '!*.log*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-d) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac; }; 
complete -F _g_bond_compl g_bond
shopt -s extglob
_g_bundle_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -ol -od -oz -ot -otr -otl -ok -okr -okl -oa -h -X -nice -b -e -dt -tu -na -z' -- $c)); return 0; fi
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
esac; }; 
complete -F _g_bundle_compl g_bundle
shopt -s extglob
_g_chi_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -c -f -o -p -ss -jc -corr -g -h -X -nice -b -e -dt -w -r0 -phi -psi -omega -rama -viol -all -shift -run -maxchi -nonormhisto -ramomega -bfact -bmax -acflen -nonormalize -P -fitfn -ncskip -beginfit -endfit' -- $c)); return 0; fi
case "$p" in
-maxchi) COMPREPLY=( $(compgen -W ' 0 1 2 3 4 5 6 ' -- $c ));;
-P) COMPREPLY=( $(compgen -W ' 0 1 2 3 ' -- $c ));;
-fitfn) COMPREPLY=( $(compgen -W ' none exp aexp exp_exp vac exp5 exp7 ' -- $c ));;
-c) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-p) COMPREPLY=( $(compgen -X '!*.pdb*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ss) COMPREPLY=( $(compgen -X '!*.dat*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-jc) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-corr) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-g) COMPREPLY=( $(compgen -X '!*.log*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac; }; 
complete -F _g_chi_compl g_chi
shopt -s extglob
_g_cluster_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -dm -o -g -dist -ev -sz -tr -ntr -clid -cl -h -X -nice -b -e -dt -tu -w -dista -nlevels -keepfree -cutoff -nofit -max -skip -av -wcl -nst -rmsmin -method -binary -M -P -seed -niter -kT' -- $c)); return 0; fi
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
esac; }; 
complete -F _g_cluster_compl g_cluster
shopt -s extglob
_g_clustsize_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -n -o -nc -mc -ac -h -X -nice -b -e -dt -w -cut -nskip -nlevels -rgblo -rgbhi' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.xpm*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-nc) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-mc) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ac) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac; }; 
complete -F _g_clustsize_compl g_clustsize
shopt -s extglob
_g_confrms_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f1 -f2 -o -n1 -n2 -h -X -nice -one -pbc' -- $c)); return 0; fi
case "$p" in
-f1) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa|gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-f2) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n1) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n2) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac; }; 
complete -F _g_confrms_compl g_confrms
shopt -s extglob
_g_covar_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -o -v -av -l -xpm -xpma -h -X -nice -b -e -dt -tu -nofit -ref -mwa -last' -- $c)); return 0; fi
case "$p" in
-tu) COMPREPLY=( $(compgen -W ' ps fs ns us ms s m h ' -- $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa|gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-v) COMPREPLY=( $(compgen -X '!*.+(trr|trj)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-av) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-l) COMPREPLY=( $(compgen -X '!*.log*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-xpm) COMPREPLY=( $(compgen -X '!*.xpm*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-xpma) COMPREPLY=( $(compgen -X '!*.xpm*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac; }; 
complete -F _g_covar_compl g_covar
shopt -s extglob
_g_density_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -n -s -ei -o -h -X -nice -b -e -dt -w -d -sl -number -ed -count' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ei) COMPREPLY=( $(compgen -X '!*.dat*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac; }; 
complete -F _g_density_compl g_density
shopt -s extglob
_g_dielectric_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -d -o -c -h -X -nice -b -e -dt -w -fft -nox1 -eint -bfit -efit -tail -A -tau1 -tau2 -eps0 -epsRF -fix -ffn -nsmooth' -- $c)); return 0; fi
case "$p" in
-ffn) COMPREPLY=( $(compgen -W ' none exp aexp exp_exp vac exp5 exp7 ' -- $c ));;
-f) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-d) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-c) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac; }; 
complete -F _g_dielectric_compl g_dielectric
shopt -s extglob
_g_dih_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -o -h -X -nice -b -e -dt -w -sa -mult' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.out*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac; }; 
complete -F _g_dih_compl g_dih
shopt -s extglob
_g_dipoles_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -enx -f -s -n -o -e -a -d -c -g -q -h -X -nice -b -e -dt -w -mu -mumax -epsilonRF -skip -temp -avercorr -gkratom -acflen -nonormalize -P -fitfn -ncskip -beginfit -endfit' -- $c)); return 0; fi
case "$p" in
-P) COMPREPLY=( $(compgen -W ' 0 1 2 3 ' -- $c ));;
-fitfn) COMPREPLY=( $(compgen -W ' none exp aexp exp_exp vac exp5 exp7 ' -- $c ));;
-enx) COMPREPLY=( $(compgen -X '!*.+(edr|ene)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-e) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-a) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-d) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-c) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-g) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-q) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac; }; 
complete -F _g_dipoles_compl g_dipoles
shopt -s extglob
_g_disre_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -s -f -ds -da -dn -dm -dr -l -n -h -X -nice -b -e -dt -w -ntop' -- $c)); return 0; fi
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
esac; }; 
complete -F _g_disre_compl g_disre
shopt -s extglob
_g_dist_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -o -h -X -nice -b -e -dt -dist' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac; }; 
complete -F _g_dist_compl g_dist
shopt -s extglob
_g_dyndom_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -o -n -h -X -nice -firstangle -lastangle -nframe -maxangle -trans -head -tail' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.pdb*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac; }; 
complete -F _g_dyndom_compl g_dyndom
shopt -s extglob
_genbox_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -cp -cs -ci -o -p -h -X -nice -box -nmol -try -seed -vdwd -shell' -- $c)); return 0; fi
case "$p" in
-cp) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-cs) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ci) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-p) COMPREPLY=( $(compgen -X '!*.top*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac; }; 
complete -F _genbox_compl genbox
shopt -s extglob
_genconf_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -o -trj -h -X -nice -nbox -dist -seed -rot -shuffle -sort -block -nmolat -maxrot -renumber' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-trj) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac; }; 
complete -F _genconf_compl genconf
shopt -s extglob
_g_enemat_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -groups -eref -emat -etot -h -X -nice -b -e -dt -w -sum -skip -nomean -nlevels -max -min -nocoul -coulr -coul14 -nolj -lj14 -bham -nofree -temp' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.+(edr|ene)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-groups) COMPREPLY=( $(compgen -X '!*.dat*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-eref) COMPREPLY=( $(compgen -X '!*.dat*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-emat) COMPREPLY=( $(compgen -X '!*.xpm*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-etot) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac; }; 
complete -F _g_enemat_compl g_enemat
shopt -s extglob
_g_energy_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -f2 -s -o -viol -pairs -ora -ort -oda -odr -odt -corr -vis -ravg -h -X -nice -b -e -w -fee -fetemp -zero -sum -dp -mutot -skip -aver -nmol -ndf -fluc -orinst -acflen -nonormalize -P -fitfn -ncskip -beginfit -endfit' -- $c)); return 0; fi
case "$p" in
-P) COMPREPLY=( $(compgen -W ' 0 1 2 3 ' -- $c ));;
-fitfn) COMPREPLY=( $(compgen -W ' none exp aexp exp_exp vac exp5 exp7 ' -- $c ));;
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
-corr) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-vis) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ravg) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac; }; 
complete -F _g_energy_compl g_energy
shopt -s extglob
_genion_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -s -n -o -g -pot -h -X -nice -np -pname -pq -nn -nname -nq -rmin -random -seed' -- $c)); return 0; fi
case "$p" in
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-g) COMPREPLY=( $(compgen -X '!*.log*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-pot) COMPREPLY=( $(compgen -X '!*.pdb*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac; }; 
complete -F _genion_compl genion
shopt -s extglob
_genpr_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -n -o -h -X -nice -fc' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.itp*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac; }; 
complete -F _genpr_compl genpr
shopt -s extglob
_g_gyrate_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -o -n -h -X -nice -b -e -dt -w -q -p' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa|gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac; }; 
complete -F _g_gyrate_compl g_gyrate
shopt -s extglob
_g_h2order_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -n -nm -s -o -h -X -nice -b -e -dt -w -d -sl' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-nm) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac; }; 
complete -F _g_h2order_compl g_h2order
shopt -s extglob
_g_hbond_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -g -sel -num -ac -dist -ang -hx -hbn -hbm -da -h -X -nice -b -e -dt -ins -a -r -abin -rbin -nonitacc -shell' -- $c)); return 0; fi
case "$p" in
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
-da) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac; }; 
complete -F _g_hbond_compl g_hbond
shopt -s extglob
_g_helix_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -s -n -f -to -cz -co -h -X -nice -b -e -dt -w -r0 -q -noF -db -ev -ahxstart -ahxend' -- $c)); return 0; fi
case "$p" in
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-to) COMPREPLY=( $(compgen -X '!*.g87*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-cz) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-co) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac; }; 
complete -F _g_helix_compl g_helix
shopt -s extglob
_g_lie_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -o -h -X -nice -b -e -dt -w -Elj -Eqq -Clj -Cqq -ligand' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.+(edr|ene)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac; }; 
complete -F _g_lie_compl g_lie
shopt -s extglob
_g_mdmat_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -mean -frames -no -h -X -nice -b -e -dt -t -nlevels' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa|gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-mean) COMPREPLY=( $(compgen -X '!*.xpm*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-frames) COMPREPLY=( $(compgen -X '!*.xpm*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-no) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac; }; 
complete -F _g_mdmat_compl g_mdmat
shopt -s extglob
_g_mindist_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -od -on -o -h -X -nice -b -e -dt -w -matrix -d -pi' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa|gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-od) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-on) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.out*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac; }; 
complete -F _g_mindist_compl g_mindist
shopt -s extglob
_g_morph_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f1 -f2 -o -or -n -h -X -nice -w -ninterm -first -last -nofit' -- $c)); return 0; fi
case "$p" in
-f1) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-f2) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-or) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac; }; 
complete -F _g_morph_compl g_morph
shopt -s extglob
_g_msd_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -o -mol -h -X -nice -b -e -dt -tu -w -type -lateral -ngroup -nomw -trestart -beginfit -endfit' -- $c)); return 0; fi
case "$p" in
-tu) COMPREPLY=( $(compgen -W ' ps fs ns us ms s m h ' -- $c ));;
-type) COMPREPLY=( $(compgen -W ' no x y z ' -- $c ));;
-lateral) COMPREPLY=( $(compgen -W ' no x y z ' -- $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa|gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-mol) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac; }; 
complete -F _g_msd_compl g_msd
shopt -s extglob
_gmxcheck_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -f2 -s1 -s2 -c -e -e2 -h -X -nice -vdwfac -bonlo -bonhi -tol -lastener' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-f2) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s1) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s2) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-c) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa|gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-e) COMPREPLY=( $(compgen -X '!*.+(edr|ene)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-e2) COMPREPLY=( $(compgen -X '!*.+(edr|ene)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac; }; 
complete -F _gmxcheck_compl gmxcheck
shopt -s extglob
_gmxdump_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -s -f -e -h -X -nice -nonr' -- $c)); return 0; fi
case "$p" in
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-e) COMPREPLY=( $(compgen -X '!*.+(edr|ene)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac; }; 
complete -F _gmxdump_compl gmxdump
shopt -s extglob
_g_nmeig_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -o -v -h -X -nice -nom -first -last' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.mtx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa|gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-v) COMPREPLY=( $(compgen -X '!*.+(trr|trj)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac; }; 
complete -F _g_nmeig_compl g_nmeig
shopt -s extglob
_g_nmens_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -v -e -s -n -o -h -X -nice -temp -seed -num -first -last' -- $c)); return 0; fi
case "$p" in
-v) COMPREPLY=( $(compgen -X '!*.+(trr|trj)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-e) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa|gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac; }; 
complete -F _g_nmens_compl g_nmens
shopt -s extglob
_g_order_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -n -s -o -od -os -h -X -nice -b -e -dt -w -d -sl -szonly -unsat' -- $c)); return 0; fi
case "$p" in
-d) COMPREPLY=( $(compgen -W ' z x y ' -- $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-od) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-os) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac; }; 
complete -F _g_order_compl g_order
shopt -s extglob
_g_potential_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -n -s -o -oc -of -h -X -nice -b -e -dt -w -d -sl -cb -ce -tz -spherical' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-oc) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-of) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac; }; 
complete -F _g_potential_compl g_potential
shopt -s extglob
_g_rama_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -o -h -X -nice -b -e -dt -w' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac; }; 
complete -F _g_rama_compl g_rama
shopt -s extglob
_g_rdf_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -o -sq -cn -hq -image -h -X -nice -b -e -dt -w -bin -com -cut -fade -grid -nlevel -wave' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa|gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-sq) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-cn) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-hq) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-image) COMPREPLY=( $(compgen -X '!*.xpm*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac; }; 
complete -F _g_rdf_compl g_rdf
shopt -s extglob
_g_rms_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -s -f -f2 -n -o -mir -a -dist -m -bin -bm -h -X -nice -b -e -dt -tu -w -what -nopbc -nofit -prev -split -skip -skip2 -max -min -bmax -bmin -nlevels' -- $c)); return 0; fi
case "$p" in
-tu) COMPREPLY=( $(compgen -W ' ps fs ns us ms s m h ' -- $c ));;
-what) COMPREPLY=( $(compgen -W ' rmsd rho rhosc ' -- $c ));;
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
esac; }; 
complete -F _g_rms_compl g_rms
shopt -s extglob
_g_rmsdist_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -equiv -o -rms -scl -mean -nmr3 -nmr6 -noe -h -X -nice -b -e -dt -w -nlevels -max -nosumh' -- $c)); return 0; fi
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
esac; }; 
complete -F _g_rmsdist_compl g_rmsdist
shopt -s extglob
_g_rmsf_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -q -oq -ox -o -od -oc -dir -h -X -nice -b -e -dt -w -res -aniso' -- $c)); return 0; fi
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
esac; }; 
complete -F _g_rmsf_compl g_rmsf
shopt -s extglob
_grompp_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -po -c -r -n -deshuf -p -pp -o -t -h -X -nice -nov -time -np -shuffle -sort -normdumbds -load -maxwarn -check14' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.mdp*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-po) COMPREPLY=( $(compgen -X '!*.mdp*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-c) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-r) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-deshuf) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-p) COMPREPLY=( $(compgen -X '!*.top*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-pp) COMPREPLY=( $(compgen -X '!*.top*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-t) COMPREPLY=( $(compgen -X '!*.+(trr|trj)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac; }; 
complete -F _grompp_compl grompp
shopt -s extglob
_g_rotacf_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -o -h -X -nice -b -e -dt -w -d -noaver -acflen -nonormalize -P -fitfn -ncskip -beginfit -endfit' -- $c)); return 0; fi
case "$p" in
-P) COMPREPLY=( $(compgen -W ' 0 1 2 3 ' -- $c ));;
-fitfn) COMPREPLY=( $(compgen -W ' none exp aexp exp_exp vac exp5 exp7 ' -- $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac; }; 
complete -F _g_rotacf_compl g_rotacf
shopt -s extglob
_g_saltbr_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -h -X -nice -b -e -dt -t -sep' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac; }; 
complete -F _g_saltbr_compl g_saltbr
shopt -s extglob
_g_sas_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -o -r -q -ao -n -i -h -X -nice -b -e -dt -w -solsize -ndots -qmax -minarea -skip -nopbc -noprot -dgs' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-r) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-q) COMPREPLY=( $(compgen -X '!*.pdb*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ao) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-i) COMPREPLY=( $(compgen -X '!*.itp*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac; }; 
complete -F _g_sas_compl g_sas
shopt -s extglob
_g_sgangle_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -n -s -oa -od -od1 -od2 -h -X -nice -b -e -dt -w' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-oa) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-od) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-od1) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-od2) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac; }; 
complete -F _g_sgangle_compl g_sgangle
shopt -s extglob
_g_sorient_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -o -no -ro -co -h -X -nice -b -e -dt -w -com -rmin -rmax -nbin' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa|gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-no) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ro) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-co) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac; }; 
complete -F _g_sorient_compl g_sorient
shopt -s extglob
_g_tcaf_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -ot -oa -o -of -oc -ov -h -X -nice -b -e -dt -w -mol -k34 -wt' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.+(trr|trj)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa|gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ot) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-oa) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-of) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-oc) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-ov) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac; }; 
complete -F _g_tcaf_compl g_tcaf
shopt -s extglob
_g_traj_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -ox -ov -of -ob -ot -ekt -ekr -cv -cf -h -X -nice -b -e -dt -tu -w -com -mol -nojump -nox -noy -noz -len' -- $c)); return 0; fi
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
-cv) COMPREPLY=( $(compgen -X '!*.pdb*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-cf) COMPREPLY=( $(compgen -X '!*.pdb*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac; }; 
complete -F _g_traj_compl g_traj
shopt -s extglob
_g_velacc_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -o -h -X -nice -b -e -dt -w -mol -acflen -nonormalize -P -fitfn -ncskip -beginfit -endfit' -- $c)); return 0; fi
case "$p" in
-P) COMPREPLY=( $(compgen -W ' 0 1 2 3 ' -- $c ));;
-fitfn) COMPREPLY=( $(compgen -W ' none exp aexp exp_exp vac exp5 exp7 ' -- $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(trr|trj)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa|gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac; }; 
complete -F _g_velacc_compl g_velacc
shopt -s extglob
_highway_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -a -h -X -nice -b -e -dt' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.dat*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-a) COMPREPLY=( $(compgen -X '!*.dat*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac; }; 
complete -F _highway_compl highway
shopt -s extglob
_make_ndx_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -n -o -h -X -nice' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac; }; 
complete -F _make_ndx_compl make_ndx
shopt -s extglob
_mdrun_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -s -o -x -c -e -g -dgdl -table -rerun -ei -eo -j -jo -ffout -devout -runav -pi -po -pd -pn -mtx -h -X -nice -deffnm -np -v -nocompact -multi -glas -ionize' -- $c)); return 0; fi
case "$p" in
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.+(trr|trj)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-x) COMPREPLY=( $(compgen -X '!*.xtc*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-c) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-e) COMPREPLY=( $(compgen -X '!*.+(edr|ene)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-g) COMPREPLY=( $(compgen -X '!*.log*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-dgdl) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-table) COMPREPLY=( $(compgen -X '!*.xvg*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-rerun) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
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
esac; }; 
complete -F _mdrun_compl mdrun
shopt -s extglob
_mk_angndx_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -s -n -h -X -nice -type' -- $c)); return 0; fi
case "$p" in
-type) COMPREPLY=( $(compgen -W ' angle g96-angle dihedral improper ryckaert-bellemans phi-psi ' -- $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac; }; 
complete -F _mk_angndx_compl mk_angndx
shopt -s extglob
_ngmx_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -h -X -nice -b -e -dt' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac; }; 
complete -F _ngmx_compl ngmx
shopt -s extglob
_pdb2gmx_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -o -p -i -n -q -h -X -nice -merge -inter -ss -ter -lys -asp -glu -his -angle -dist -posrefc -una -nosort -H14 -ignh -alldih -dummy -heavyh -deuterate' -- $c)); return 0; fi
case "$p" in
-dummy) COMPREPLY=( $(compgen -W ' none hydrogens aromatics ' -- $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-p) COMPREPLY=( $(compgen -X '!*.top*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-i) COMPREPLY=( $(compgen -X '!*.itp*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-q) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac; }; 
complete -F _pdb2gmx_compl pdb2gmx
shopt -s extglob
_protonate_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -s -f -n -o -h -X -nice -b -e -dt' -- $c)); return 0; fi
case "$p" in
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa|gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac; }; 
complete -F _protonate_compl protonate
shopt -s extglob
_tpbconv_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -s -f -n -o -h -X -nice -time -extend -until -zeroq -nounconstrained' -- $c)); return 0; fi
case "$p" in
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(trr|trj)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac; }; 
complete -F _tpbconv_compl tpbconv
shopt -s extglob
_trjcat_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -o -n -h -X -nice -tu -b -e -dt -prec -novel -settime -nosort -cat' -- $c)); return 0; fi
case "$p" in
-tu) COMPREPLY=( $(compgen -W ' ps fs ns us ms s m h ' -- $c ));;
-o) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac; }; 
complete -F _trjcat_compl trjcat
shopt -s extglob
_trjconv_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -o -s -n -fr -h -X -nice -b -e -tu -w -skip -dt -dump -t0 -timestep -pbc -ur -center -box -shift -fit -pfit -ndec -novel -force -trunc -exec -app -split -sep' -- $c)); return 0; fi
case "$p" in
-tu) COMPREPLY=( $(compgen -W ' ps fs ns us ms s m h ' -- $c ));;
-pbc) COMPREPLY=( $(compgen -W ' none whole inbox nojump ' -- $c ));;
-ur) COMPREPLY=( $(compgen -W ' rect tric compact ' -- $c ));;
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa|gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-fr) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac; }; 
complete -F _trjconv_compl trjconv
shopt -s extglob
_trjorder_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -n -o -h -X -nice -b -e -dt -na -da' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa|gro|g96|pdb|brk|ent)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-n) COMPREPLY=( $(compgen -X '!*.ndx*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac; }; 
complete -F _trjorder_compl trjorder
shopt -s extglob
_wheel_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -o -h -X -nice -r0 -rot0 -T -nonn' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.dat*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.eps*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac; }; 
complete -F _wheel_compl wheel
shopt -s extglob
_x2top_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -o -r -h -X -nice -scale -nexcl -H14 -alldih -nopairs -name -nopbc -param -noround -kb -kt -kp' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.+(gro|g96|pdb|brk|ent|tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-o) COMPREPLY=( $(compgen -X '!*.top*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-r) COMPREPLY=( $(compgen -X '!*.rtp*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac; }; 
complete -F _x2top_compl x2top
shopt -s extglob
_xpm2ps_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -f2 -di -do -o -xpm -h -X -nice -w -noframe -title -yonce -legend -diag -combine -bx -by -rainbow -gradient -skip -zeroline' -- $c)); return 0; fi
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
esac; }; 
complete -F _xpm2ps_compl xpm2ps
shopt -s extglob
_xrama_compl() {
local p c
COMPREPLY=() c=${COMP_WORDS[COMP_CWORD]} p=${COMP_WORDS[COMP_CWORD-1]}
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen  -W ' -f -s -h -X -nice -b -e -dt' -- $c)); return 0; fi
case "$p" in
-f) COMPREPLY=( $(compgen -X '!*.+(xtc|trr|trj|gro|g96|pdb|g87)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
-s) COMPREPLY=( $(compgen -X '!*.+(tpr|tpb|tpa)*(.gz|.Z)' -f $c ; compgen -S '/' -X '.*' -d $c ));;
esac; }; 
complete -F _xrama_compl xrama
