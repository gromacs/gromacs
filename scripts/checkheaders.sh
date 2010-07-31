#! /bin/bash

incldir="../include"
rm -f gromacs
ln -s $incldir gromacs

if [ -z "$1" ]; then
  files=$(cd $incldir; find -name "*.h" | sed 's/^\./gromacs/')
else
  files="$@"
fi

for i in $files; do
  echo $i
  cat << EOF > t.c
#include <$i>
int main(){
  return 0;
}
EOF
  gcc -I. -c t.c -D bool=int || echo "Failed"
done
rm -f gromacs t.[co]
