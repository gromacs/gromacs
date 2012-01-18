mkdir releng
cd releng
git init
git fetch ssh://jenkins@gerrit.gromacs.org:29418/releng refs/heads/1
git checkout -q -f FETCH_HEAD
cd ..
python -u releng/GerritBuild.py
