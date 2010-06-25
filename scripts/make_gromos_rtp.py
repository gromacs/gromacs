#!/usr/freeware/bin/python

# usage:  make_gromos_rtp.py > ffG43a1.rtp
# this script tries to make a residue topology database for GROMACS from
# the GROMOS version of this file. It needs ffG43a1.atp to fill in the
# atom types for the atom integer codes. It converts until it finds the string
# "#RNMES", which indicates the start of the solvent part
# of mtb43a1.dat. Solvents are treated differently in GROMACS
# The vaiables GROMOS_FILE and ATOMTYPE_FILE has to be specified as required.

# The output file requires some manual modifications:
#   add ACE and NH2
#   add 4 dihedrals with parameters 0 0 2 in HEME
#   remove I in front of all atoms names in SO42-
#   make a copy of H2O called HOH

# author: Klaas Dijkstra, Univ. Groningen
# Oct. 2000
# minor modifications by Berk Hess


from string import atoi, split

GROMOS_FILE   = 'mtb43a1.dat'
ATOMTYPE_FILE = 'ffG43a1.atp'

START  = '# RNME\012'
STOP   = '#RNMES\012'

NATOM  = '# NMAT,NLIN\012'
EXCLUS = '#ATOM                               MAE MSAE\012'
ATOM1  = '#ATOM ANM  IACM MASS        CGMICGM MAE MSAE\012'
ATOM2  = '#ATOM ANM  IACM MASS        CGM ICGM  MAE MSAE\012'
ATOMtr = '#ATOM ANM  IACM MASS        CGMICGM\012'
NB     = '#  NB\012'
NBA    = '# NBA\012'
NIDA   = '# NIDA\012'
NDA    = '# NDA\012'


#### functions
def translate(t):               # translate number into atom code
    if    t ==  0        : t = '-O'
    elif  t == -1        : t = '-C'
    elif  t == -2        : t = '-CA'
    elif  t == f.natom+1 : t = '+N'
    else                 : t = f.atoms[t-1][1]
    return t


def findbonds(atomnumber):
    atom = eval(f.atoms[atomnumber][0])

    onebond = []
    twobond = []
    ext     = [] 
    bond    = []
    ind     = []
    excl    = []


    for i in f.bonds:
        a, b = i[0:2]
        bond.append(eval(a),eval(b))
        bond.append(eval(b),eval(a))

    for i in bond:
        if i[0] == atom: onebond.append(i[1])

    for i in bond:
        for j in onebond:
            if i[0] == j and atom < i[1]: twobond.append(i[1])

    ext = onebond
    ext.extend(twobond)

    for i in ext:
        if i<atom:
           ind.append(i)
    for i in ind:
        ext.remove(i)

    ext.sort()

    for i in ext:
        excl.append(`i`)

    return excl


#### general class to unravel a gromos topology file
class Cin:
   def __init__(self, filename):
       f=open(filename)
       self.G = open(filename).readlines()
       self.nlines = len(self.G)

   def index(self):
       " Find all indexes of lines with string '# RNME\012'and '#RNMES\012' "
       self.ind = []
       for i in range(self.nlines):
           if self.G[i] == STOP:
              self.ind.append(i)
              return
           if self.G[i] == START:
              self.ind.append(i)

   def mkres(self):
       " Returns list of residues "
       self.all = []
       length = len(self.ind)
       for i in range(length):
           res = []
           start = self.ind[i] + 1
           try:    stop  = self.ind[i+1] - 2
           except: return
           for j in range(start,stop):
               res.append(self.G[j])
           self.all.append(res)

   def gets(self,list,string):
       " Returns linenumber of string in list"
       for i in range(len(list)):
                if list[i] == string:
                   return i
       print >> sys.stderr "Could not find string",string,"in list of length",len(list)
       return -1

#--------------------------#
# unravel gromos list

   def getRESname(self, res):
       " Get the name of the current residue "
       self.residue = split(res[0])

   def getNATOM(self, res):
       " Get number of atoms, number of preceding excl. "
       ind = self.gets(res, NATOM)
       self.natom, self.nlin = split(res[ind+1])
       self.natom = atoi(self.natom)
       self.nlin   = atoi(self.nlin)

   def getEXCLUS(self, res):
       " Get preceding exclusions "
       self.exclus = []
       ind = self.gets(res, EXCLUS)
       for i in range(self.nlin):
           self.exclus.append(split(res[ind+i+1]))

   def getATOM(self, res):
       " Get atoms"
       self.atoms = []
       try:     ind = self.gets(res, ATOM1)
       except:  ind = self.gets(res, ATOM2)
       i = 0
       cntr = 0
       while cntr < self.natom - self.nlin:
           i = i + 1
           line = split(res[ind+i])		# get next line
           noflo = (atoi(line[6])-1)/8		# if MAE > 8 cont. on next line
           if noflo < 0: noflo = 0
           for j in range(noflo):
               i = i + 1
               line1 = split(res[ind+i])	# overflow line
               "print line1"
               line = line + line1
           self.atoms.append(line)
           cntr = cntr + 1

   def getATOMtr(self, res):
       " Get trailing atoms"
       self.atomtr = []
       try:    ind = self.gets(res, ATOMtr)
       except: return
       for i in range(self.nlin):
           self.atomtr.append(split(res[ind+i+1]))
           self.atoms.append(split(res[ind+i+1]))

   def getBOND(self, res):
       " Get bonds"
       self.bonds = []
       ind = self.gets(res, NB)
       self.nb = atoi(split(res[ind+1])[0])
       j = 0
       for i in range(self.nb):
           line = split(res[ind+i+j+3])
           if line[0] == '#':
              line = split(res[ind+i+j+4])
              j = j+1
           self.bonds.append(line)
       
   def getNBA(self, res):
       " Get bond angles"
       self.ba = []
       ind = self.gets(res, NBA)
       self.nba = atoi(split(res[ind+1])[0])
       j = 0
       for i in range(self.nba):
           line = split(res[ind+i+j+3])
           if line[0] == '#':
              line = split(res[ind+i+j+4])
              j = j + 1
           self.ba.append(line)
       
   def getNIDA(self, res):
       " Get improper dihedrals"
       self.ida = []
       ind = self.gets(res, NIDA)
       self.nida = atoi(split(res[ind+1])[0])
       j = 0
       for i in range(self.nida):
           line = split(res[ind+i+j+3])
           if line[0] == '#':
              line = split(res[ind+i+j+4])
              j = j + 1
           self.ida.append(line)

   def getNDA(self, res):
       " Get dihedrals"
       self.da = []
       ind = self.gets(res, NDA)
       j = 0
       self.nda = atoi(split(res[ind+1])[0])
       for i in range(self.nda):
           line = split(res[ind+i+j+3])
           if line[0] == '#':
              line = split(res[ind+i+j+4])
              j = j + 1
           self.da.append(line)


#-----------------------------#
# main program

typ       = open(ATOMTYPE_FILE)		# translate numbers to atoms
typelines = typ.readlines()
for i in range(len(typelines)):
    typelines[i]=split(typelines[i])

f=Cin(GROMOS_FILE)			# bind class instance
f.index()				# mark all residues (f.ind)
f.mkres()				# put all residues into list (f.all)

start = 0; stop = 92

" The rtp header "
print "[ bondedtypes ]"
print "; bonds  angles  dihedrals  impropers"
print "    2       2          1          2"


for resnum in range(start,stop):	# loop through all residues
    f.getRESname(f.all[resnum])		# residue name
    f.getNATOM  (f.all[resnum])		# number of atoms

    if f.nlin != 0:			# 0 for a separate molecule
          f.getEXCLUS(f.all[resnum])	# number of exclusions

    f.getATOM   (f.all[resnum])		# atoms			=> f.atoms
    f.getATOMtr (f.all[resnum])		# trailing atoms	=> f.atomtr
    f.getBOND   (f.all[resnum])		# bonds			=> f.bonds
    f.getNBA    (f.all[resnum])		# bond angles		=> f.ba
    f.getNIDA   (f.all[resnum])		# improper dihedrals	=> f.ida
    f.getNDA    (f.all[resnum])		# dihedrals		=> f.da
    


    # output to Gromacs format
    #-------------------------#

#####
    # atoms
    print ""
    print "[",f.residue[0],"]"
    print " [ atoms ]"
    chargegroup = 0
    for j in range(f.natom - f.nlin):
      try:
        atomtype = atoi(f.atoms [j][2]) - 1
        atomfield = typelines[atomtype][0]
        print "%5s %5s %11s %5s" % \
                  (f.atoms [j][1],atomfield,f.atoms [j][4],chargegroup)
        chargegroup = chargegroup + atoi(f.atoms[j][5])
      except:
        print j

#####
    # trailing atoms
    for j in range(f.nlin):
        atomtype = atoi(f.atomtr [j][2]) - 1
        atomfield = typelines[atomtype][0]
        print "%5s %5s %11s %5s" % \
                  (f.atomtr[j][1],atomfield,f.atomtr[j][4][:-2],chargegroup)
        chargegroup = chargegroup + atoi(f.atomtr[j][5])
    
#####
    # bonds
    print " [ bonds ]"
    for j in range(f.nb):
        t1 = atoi(f.bonds [j][0])
        t2 = atoi(f.bonds [j][1])
        " Only special bonds go to 0 or less "
	if t1 >= 1  and t2 >= 1:
            t1 = translate(t1)
            t2 = translate(t2)
            print "%5s %5s    gb_%-5s" % \
                  (t1, t2, f.bonds[j][2])

#####
    # exclusions
    ne = 0
    for j in range(f.natom - f.nlin):
        aaa = findbonds(j)
        bbb = f.atoms[j][7:]
        for i in aaa:
            if i in bbb: bbb.remove(i)
        for i in bbb:
            " Ignore special exclusions "
            t1 = atoi(i)
            if t1 >= 0:
                t1 = translate(t1)
                t2 = atoi(f.atoms[j][0])
                t2 = translate(t2)
                if ne == 0:  print " [ exclusions ]\n;  ai    aj"
                print "%5s %5s" % (t2,t1)
                ne = ne + 1

#####
    # angles
    print " [ angles ]"
    print ";  ai    aj    ak   gromos type"
    for j in range(f.nba):
        t1 = atoi(f.ba [j][0])
        t2 = atoi(f.ba [j][1])
        t3 = atoi(f.ba [j][2])
        if t1 >= -2 and t2 >= -2 and t3 >= -2:
            t1 = translate(t1)
            t2 = translate(t2)
            t3 = translate(t3)
            print "%5s %5s %5s     ga_%-5s" % \
                  (t1,t2,t3,f.ba[j][3])

#####
    # improper dihedrals
    print " [ impropers ]"
    print ";  ai    aj    ak    al   gromos type"
    for j in range(f.nida):
        t1 = atoi(f.ida [j][0])
        t2 = atoi(f.ida [j][1])
        t3 = atoi(f.ida [j][2])
        t4 = atoi(f.ida [j][3])
        if t1 >= -2 and t2 >= -2 and t3 >= -2 and t4 >= -2:
            t1 = translate(t1)
            t2 = translate(t2)
            t3 = translate(t3)
            t4 = translate(t4)
            print "%5s %5s %5s %5s     gi_%-5s" % \
                  (t1,t2,t3,t4,f.ida[j][4])

#####
    # dihedrals
    print " [ dihedrals ]"
    print ";  ai    aj    ak    al   gromos type"
    for j in range(f.nda):
        t1 = atoi(f.da [j][0])
        t2 = atoi(f.da [j][1])
        t3 = atoi(f.da [j][2])
        t4 = atoi(f.da [j][3])
        if t1 >= -2 and t2 >= -2 and t3 >= -2 and t4 >= -2:
            t1 = translate(t1)
            t2 = translate(t2)
            t3 = translate(t3)
            t4 = translate(t4)
            print "%5s %5s %5s %5s     gd_%-5s" % \
                  (t1,t2,t3,t4,f.da[j][4])
