mkdir releng
cd releng
git init && git fetch git://git.gromacs.org/releng.git refs/heads/4.6.0 && git checkout -q -f FETCH_HEAD && cd .. && python -u releng/GerritBuild.py
