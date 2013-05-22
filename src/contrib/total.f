C Chemical shift calculation, to read in multiple NMR structures
C (with protons) and calculate sd for each
C Will also read Xray structures and add protons
C Current limit is 5000 heavy atoms 3000 protons and 50 rings
C Author - Mike Williamson Jan 94
C see M P Williamson and T Asakura, J Magn Reson Ser B 101 63-71 1993
	DIMENSION SHIFT(3000)	!Shift values for each proton
	DIMENSION ANISCO(3000),ANISCN(3000),SHIFTE(3000)
	DIMENSION sum(3000),exptln(3000)
	DIMENSION EXPTL(3000)	!Exptl shifts (external file)
	DIMENSION IRES(5000),RES(5000),ANAME(5000),FINAL(3000)
	DIMENSION X(5000), Y(5000), Z(5000) !Heavy atom input
	DIMENSION IRESH(3000),RESH(3000),ANAMEH(3000)
	DIMENSION XH(3000), YH(3000), ZH(3000) !H atom input
	DIMENSION vx(3),vy(3),vz(3),CEN(3),rh(3),pregam(3),yy(3)
C (Anisotropy tensor plus anisotropy centre)
	DIMENSION vca(3),vc(3),vo(3),vn(3),calctemp(3000)
	DIMENSION datres(179),datnam(179),datshift(179)
	DIMENSION RCNET(3000)
	DIMENSION IACODE(50),IRATPT(9,50),IRATS(50)
	DIMENSION shiftt(3000),exptlt(3000),VCHN(3)
	dimension rn(3,50),cent(3,50),signr(50)
        dimension iexp(2000),atexp(2000),pexp(2000),eexp(2000)
        CHARACTER ANAME*4,RES*3,FN1*60,FN3*40,junk*80
	CHARACTER ans*1,FN4*40,pexp*3,resh*3,anameh*4,prott*3
	CHARACTER datres*3,datnam*4,YN*1,atexp*4,FN9*60,RFILE*60
	integer AT1,AT2,AT3,INDX(3000)
	COMMON x,y,z,xh,yh,zh
C Set values for anisotropies, atomic charges and multiplication
C factor for electric field: shift=-Ez*sigmaE
	DATA XCO1/-13.0/,XCO2/-4.0/,XCN1/-11.0/,XCN2/1.40/
	DATA sigmaE/0.60/,Qc/1.60/,Qn/-1.70/,Qo/-2.30/
	DATA Qhn/0.70/
C NB The atomic charges are in esu
C
C******************FILE INPUT************************************
	iwarning=1
	type *,'Input file?'
	read(5,999)FN1
997	format (I)
98	FORMAT (I5)
999	format(A)
        type *,'Output file?'
        read(5,999)FN3
	type *,'Calculate for HA[1], HN[2] or all protons[3]?'
	read(5,997)icalculate 
	print *,'icalculate=',icalculate
C NB Actually calculates all protons, just prints differently
99	FORMAT(A80)
	type *,' Are you comparing calc. to experimental shifts? [N]'
	read(5,970)YN
970	format(A)
        OPEN(1,file=FN1,STATUS='OLD',readonly)
	OPEN(3,file=FN3)
	type *,' Random file ?'
	read(5,999)RFILE
	open(4,file=RFILE,status='old',readonly)
	do 25 I=1,179		!Random coil shifts
	read(4,998) datres(I),datnam(I),datshift(I)
998	format(A3,1X,A4,F5.2)
25	continue
        if(yn.eq.'Y'.or.yn.eq.'y') then
	type *,' File containing experimental shifts?'
	read(5,999)FN4
	OPEN(UNIT=8,FILE=FN4,STATUS='OLD')
C This file should contain 2-line header, no. of protons (I5) plus list
        type *,'Name of protein in this file?'
	read(5,992) prott
992     format(A3)
	READ (8,99) JUNK
	READ (8,99) JUNK
	READ (8,98) Iprot
	do 900 II=1,Iprot
	read (8,996) iexp(ii),atexp(ii),pexp(ii),eexp(ii)
996	format(I3,7X,A4,5X,A3,27X,F9.5)
900     continue
        endif
        imodel=0
776     iheavyatom=0
        iproton=0
C PDB file should end with TER or END or preferably both
        do 10 I=1,100000
	READ (1,99,end=777) JUNK
        if (junk(1:5).eq.'MODEL') then
         read(junk(6:14),'(I9)') nmodel
         imodel=1
        endif
        if (junk(1:6).eq.'ENDMDL') goto 11 
	if (junk(1:4).eq.'END') then
	  if (imodel.eq.1) goto 777
          goto 11
        endif
C Replace D by H for neutron structures
	if (junk(14:14).eq.'D') junk(14:14)='H'
        if (junk(1:4).eq.'TER ') goto 11
        if (junk(1:4).eq.'ATOM') then
	 if(junk(27:27).ne.' ') print 601,junk(23:26)
601       format('WARNING: substituted residue present at ',A4)
	 if(junk(22:22).ne.' '.and.iwarning.eq.1) then
	  iwarning=2
	  type *,'WARNING: more than one chain present'
	 endif
	 if(junk(17:17).ne.' ') print 602,junk(23:26)
602       format('WARNING: alternate conformations at residue ',A4)
C Next line reads HN in as heavy atom, for Electric field calc.
         if((junk(14:14).ne.'H').or.(junk(13:15).eq.' H ')) then
         iheavyatom=iheavyatom+1
         read(junk(13:16),'(A4)') aname(iheavyatom)
         read(junk(18:20),'(A3)') res(iheavyatom)
         read(junk(23:26),'(I4)') ires(iheavyatom)
         read(junk(31:54),'(3F8.0)') x(iheavyatom),
     .     y(iheavyatom),z(iheavyatom)
	 endif
         if(junk(14:14).eq.'H') then
         if((icalculate.eq.1).and.(junk(14:15).ne.'HA'))goto 10 
         if((icalculate.eq.2).and.(junk(14:15).ne.'H '))goto 10 
         iproton=iproton+1 
         read(junk(13:16),'(A4)') anameh(iproton)
         read(junk(18:20),'(A3)') resh(iproton)
         read(junk(23:26),'(I4)') iresh(iproton)
         read(junk(31:38),'(F8.3)') xh(iproton)
         read(junk(39:46),'(F8.3)') yh(iproton)
         read(junk(47:54),'(F8.3)') zh(iproton)
         endif
        endif
10	 CONTINUE
11       continue
C To avoid calculating for water molecules
         if(res(1).eq.'WAT'.or.res(2).eq.'WAT') goto 777
C Now see if protons are present and add them if not
C
C       Next Line hacked (DvdS, 12/94)
	if ((iproton.eq.0) .or. (icalculate.eq.3)) then
	type *,'Adding protons'
	type *,'Print out file with protons?'
	read(5,607)ANS
	call addprot(x,y,z,xh,yh,zh,Iheavyatom,Iproton,
     &   aname,anameh,res,resh,ires,iresh)
	if(ans.eq.'Y'.or.ans.eq.'y') then
         do 670 I=1,60
	  if(FN1(I:I).eq.' ')goto 671
670      continue
671      ilength=I-1
	  if(ilength.lt.40) then
           FN9=FN1(1:ilength-4)//'_protonated.pdb'
	  else
	   FN9=FN1(ilength-8:ilength-4)//'_protonated.pdn'
	  endif
	  print 669,FN9
669       format('Output going to ',A)
	  open(9,file=FN9)
	  write(9,660)FN1 
660       format('REMARK  Generated from ',60A)
	  iresstart=ires(1)
	  iresend=ires(iheavyatom)
	  iline=0
	  do 668 icount=iresstart,iresend
	  do 661 I=1,iheavyatom
	  if(aname(I).eq.' H  ') goto 661
	  if(ires(I).eq.icount) then 
	  iline=iline+1
          write(9,662)iline,aname(I),res(I),ires(I),x(I),y(I),z(I)
	  endif
662       format('ATOM',I7,1X,A4,1X,A3,2X,I4,4X,3F8.3)
661       continue
          do 663 I=1,iproton
          if(iresh(I).eq.icount) then
	  iline=iline+1
          write(9,662)iline,anameh(I),resh(I),iresh(I),xh(I),yh(I),zh(I)
	  endif
663       continue
668       continue
          write(9,665)
          write(9,666)
665       format('TER')
666       format('END')
        end if
       end if
         do 20 I=1,iproton
	 shift(I)=0.0		!initialise
	 anisco(I)=0.0
	 aniscn(I)=0.0
	 shiftE(I)=0.0
20	 continue
607	FORMAT(A1)
	type *,' All atoms read...initialising'
        if(imodel.eq.1) type 774, nmodel
774     format(' Calculating shifts for model',I9)
C***********Calculate number of aromatic residues, rearrange order of
C ring atoms, and set pointers to line numbers of aromatic atoms
C (NB each Trp counts as two rings)
	narom=0
	do 105 I=1,iheavyatom
	if((res(I).eq.'TRP'.or.res(I).eq.'TYR'.or.res(I).eq.'PHE'.
     1   or.res(I).eq.'HIS').and.(aname(I).eq.' CB ')) THEN
	 narom = narom + 1
         do 102 Iring=1,6
          IRATPT(Iring,narom)=0
102      continue
	 IACODE(narom)=IRES(I)
	  IF ((res(I).eq.'PHE').or.(res(I).eq.'TYR')) THEN
	    IRATS(narom)=6
	    do 101 K=1,20
	    IF(ires(I+K).eq.ires(I).and.
     &       aname(I+K).eq.' CG ') IRATPT(1,narom)=I+K
	    IF(ires(I+K).eq.ires(I).and.
     &       aname(I+K).eq.' CD1') IRATPT(2,narom)=I+K
	    IF(ires(I+K).eq.ires(I).and.
     &       aname(I+K).eq.' CE1') IRATPT(3,narom)=I+K
	    IF(ires(I+K).eq.ires(I).and.
     &       aname(I+K).eq.' CZ ') IRATPT(4,narom)=I+K
	    IF(ires(I+K).eq.ires(I).and.
     &       aname(I+K).eq.' CE2') IRATPT(5,narom)=I+K
	    IF(ires(I+K).eq.ires(I).and.
     &       aname(I+K).eq.' CD2') IRATPT(6,narom)=I+K
101	    continue
	  ELSE IF (res(I).eq.'HIS') THEN
	    IRATS(narom)=5
	    do 103 K=1,20
	    IF(ires(I+K).eq.ires(I).and.
     &       aname(I+K).eq.' CG ') IRATPT(1,narom)=I+K
	    IF(ires(I+K).eq.ires(I).and.
     &       aname(I+K).eq.' CD2') IRATPT(2,narom)=I+K
	    IF(ires(I+K).eq.ires(I).and.
     &       aname(I+K).eq.' NE2') IRATPT(3,narom)=I+K
	    IF(ires(I+K).eq.ires(I).and.
     &       aname(I+K).eq.' CE1') IRATPT(4,narom)=I+K
	    IF(ires(I+K).eq.ires(I).and.
     &       aname(I+K).eq.' ND1') IRATPT(5,narom)=I+K
103	    continue
	  ELSE
C Trp counts as two aromatic rings (5 then 6)
	    IRATS(narom)=5
	    do 104 K=1,20        
	    IF(ires(I+K).eq.ires(I).and.
     &       aname(I+K).eq.' CG ') IRATPT(1,narom)=I+K
	    IF(ires(I+K).eq.ires(I).and.
     &       aname(I+K).eq.' CD1') IRATPT(2,narom)=I+K
	    IF(ires(I+K).eq.ires(I).and.
     &       aname(I+K).eq.' NE1') IRATPT(3,narom)=I+K
	    IF(ires(I+K).eq.ires(I).and.
     &       aname(I+K).eq.' CE2') IRATPT(4,narom)=I+K
	    IF(ires(I+K).eq.ires(I).and.
     &       aname(I+K).eq.' CD2') IRATPT(5,narom)=I+K
104	    continue
            narom = narom + 1
	    IACODE(narom)=IRES(I)
	    IRATS(narom)=6
	    do 106 K=1,25
	    IF(ires(I+K).eq.ires(I).and.
     &       aname(I+K).eq.' CE2') IRATPT(1,narom)=I+K
	    IF(ires(I+K).eq.ires(I).and.
     &       aname(I+K).eq.' CD2') IRATPT(2,narom)=I+K
	    IF(ires(I+K).eq.ires(I).and.
     &       aname(I+K).eq.' CE3') IRATPT(3,narom)=I+K
	    IF(ires(I+K).eq.ires(I).and.
     &       aname(I+K).eq.' CZ3') IRATPT(4,narom)=I+K
	    IF(ires(I+K).eq.ires(I).and.
     &       aname(I+K).eq.' CH2') IRATPT(5,narom)=I+K
	    IF(ires(I+K).eq.ires(I).and.
     &       aname(I+K).eq.' CZ2') IRATPT(6,narom)=I+K
106	    continue
	  ENDIF
	endif	!End of loop for each aromatic residue
105	continue
        do 107 J=1,narom
	if(iratpt(1,j).eq.0.or.iratpt(2,j).eq.0.or.iratpt(3,j).eq.0.
     &   or.iratpt(4,j).eq.0.or.iratpt(5,j).eq.0) print 960,
     &   ires(iacode(j)),res(ires(iacode(j)))
107     continue
960	format(1X,'Ring atom missing for ',I4,1X,A3)
	DO 115 L=1,iproton	!Initialise all ring current shifts to zero
	  RCNET(L)=0.
115	CONTINUE
C ******************************************************************
	type *,' Calculating ring current shifts....'
C
      DO 116 J=1,NAROM		!FOR EACH AROMATIC RESIDUE
C Now set up rn (ring normal) and cent for each ring.
C Determined using atoms 1, 3 and 5 only of ring.
C Much of this is lifted from a version of AMBER provided by Dave Case
        call plane(IRATPT(1,j),IRATPT(3,j),IRATPT(5,j),x,y,z,rn(1,j),
     .               cent(1,j))
c  From pts 1, 3 and 5 calculates rn(ring normal), drn (deriv. of
c  ring normal) and centre of ring
c
c  --   check on signr of normal vector
c
            signr(j) = 1.0
            d1cx = x(IRATPT(1,j)) - cent(1,j)
            d1cy = y(IRATPT(1,j)) - cent(2,j)
            d1cz = z(IRATPT(1,j)) - cent(3,j)
            d2cx = x(IRATPT(3,j)) - cent(1,j)
            d2cy = y(IRATPT(3,j)) - cent(2,j)
            d2cz = z(IRATPT(3,j)) - cent(3,j)
            vp1 = d1cy*d2cz - d1cz*d2cy
            vp2 = d1cz*d2cx - d1cx*d2cz
            vp3 = d1cx*d2cy - d1cy*d2cx
            if ((vp3*rn(3,j)+vp2*rn(2,j)+vp1*rn(1,j))
     .               .gt.0.0) signr(j) = -1.0
  119     DO 120 L=1,iproton	!For each proton
	ip=L
C Next line includes self aromatic shift for HA only
         IF(IRESH(L).EQ.IRES(IRATPT(1,J)).AND.(ANAMEH(L)(2:3).NE.
     &     'HA'.AND.ANAMEH(L)(2:3).NE.'H ')) GOTO 120
c
c  --     skip rings whose centre is more than 15A away
c
            relc = (xh(ip)-cent(1,j))**2 + (yh(ip)-cent(2,j))**2
     .            +(zh(ip)-cent(3,j))**2
            if (relc.gt.225.) go to 120
c
c  --   loop over pairs of bonded atoms in the ring
c
          do 80 k=1,irats(j)
            kp1 = k + 1
            if (kp1.gt.irats(j)) kp1 = 1
c
            call hm(ip,iratpt(k,j),iratpt(kp1,j),x,y,z,xh,yh,zh,
     . rn(1,j),shifthm)
            rcnet(ip) = rcnet(ip)+ signr(j)*5.4548*shifthm
   80     continue
  100   continue
c
120	continue !PROTON LOOP
  116 CONTINUE
C *************BEGIN ANISOTROPY CALCULATION**********************
C ****Locate all backbone C=O and calculate anisotropic shift****
	type *,' Calculating CO anisotropy...'
	do 30 I=1,iheavyatom
	  if (aname(I).ne.' CA ') goto 30
	  DO 35 J=1,15
	    IF ((ANAME(I+J).EQ.' C  ').AND.(IRES(I).EQ.IRES(I+J)))
     1       then
	     Inext=I+J
	     GOTO 38
	    ENDIF
35	  CONTINUE
	  print 95,IRES(I)
95	  format (' C atom not found for residue ',I5)
	goto 30
38	    IF (ANAME(Inext+1).NE.' O  ') then
 	     print 94,IRES(I)
94	     format (' O atom not found for residue ',I5)
	     goto 30
	    ENDIF
C Atoms CA, C and O now located at I, Inext and (Inext+1)
	  at3=I	  
	  at1=Inext
	  at2=Inext+1
	CALL VEC(at1,at2,at3,vx,vy,vz,r)
C Returns vectors vx,vy,vz, and distance r between at1 and at2
	CALL CENTRE(vz,x(at1),y(at1),z(at1),1.1,cen)
C
C Now loop through all H, calculating shift due to this C=O
	do 40 K=1,iproton
	  IF(IRESH(K).EQ.IRES(AT1)) GOTO 40	!SELF SHIFT
	  CALL RHAT(xh(k),yh(k),zh(k),cen,rh,rlen)
C Returns vector rh, length rlen, between HA and CEN
	  IF(RLEN.GT.12.0) GOTO 40	!distance cutoff
	  CALL VPROD(rh,vy,pregam,stheta,tempab)
C Returns sine of angle theta between rh and vy
	  CALL VPROD(pregam,vx,yy,sgamma,tempab)
	  calc1=XCO1*((3.0*stheta*stheta)-2.0)
	  calc2=XCO2*(1.0-(3.0*stheta*stheta*sgamma*sgamma))
	  calc3=(calc1+calc2)/(3.0*rlen*rlen*rlen)
	  shift(k)=shift(k)-calc3	  
	  anisco(k)=anisco(k)-calc3
40	continue
30	continue
C	if(XCO1.ne.9999.9) goto 9999 !To skip sidechain CO calculation
C *********************************************************
C *******************Sidechain CO from Asn, Gln************
	do 530 I=1,iheavyatom
	if(res(I).ne.'ASN'.and.res(I).ne.'GLN') goto 530
	  if (aname(I).ne.' CB ') goto 530
	  if(res(I).eq.'ASN') then
	   if (aname(I+1).ne.' CG '.or.aname(I+2).ne.' OD1') then
		print 595, IRES(I)
	        goto 530
595	  format (' Missing atoms for ASN ',I5)
	   endif
	  at1=I+1
	  at2=I+2
	  at3=I
	  else if(res(I).eq.'GLN') then
	   if (aname(I+1).ne.' CG '.or.aname(I+2).ne.' CD '.or.
     1      aname(I+3).ne.' OE1')then
		print 596, IRES(I)
	        goto 530
596	  format (' Missing atoms for GLN ',I5)
	   endif
	  at1=I+2
	  at2=I+3
	  at3=I+1
	 endif
	CALL VEC(at1,at2,at3,vx,vy,vz,r)
	CALL CENTRE(vz,x(at1),y(at1),z(at1),1.1,cen)
C Now loop through all HA, calculating shift due to this C=O
	do 540 K=1,iproton
	  CALL RHAT(xh(k),yh(k),zh(k),cen,rh,rlen)
	  IF(RLEN.GT.12.0) GOTO 540
	  CALL VPROD(rh,vy,pregam,stheta,tempab)
	  CALL VPROD(pregam,vx,yy,sgamma,tempab)
	  calc1=XCO1*((3.0*stheta*stheta)-2.0)
	  calc2=XCO2*(1.0-(3.0*stheta*stheta*sgamma*sgamma))
	  calc3=(calc1+calc2)/(3.0*rlen*rlen*rlen)
	  shift(k)=shift(k)-calc3	  
	  anisco(k)=anisco(k)-calc3
540	continue
530	continue
C *****************************************************
C *****************Sidechain CO from Asp, Glu**********
	do 630 I=1,iheavyatom
	if(res(I).ne.'ASP'.and.res(I).ne.'GLU') goto 630
	  if (aname(I).ne.' CB ') goto 630
	  if(res(I).eq.'ASP') then
	   if (aname(I+1).ne.' CG '.or.aname(I+2).ne.' OD1'.or.
     1       aname(I+3).ne.' OD2') then
		print 695, IRES(I)
	        goto 630
695	  format (' Missing atoms for ASP ',I5)
	   endif
	  at1=I+1
	  at2=I+2
	  at3=I
	  else if(res(I).eq.'GLU') then
	   if (aname(I+1).ne.' CG '.or.aname(I+2).ne.' CD '.or.
     1      aname(I+3).ne.' OE1'.or.aname(I+4).ne.' OE2')then
		print 696, IRES(I)
	        goto 630
696	  format (' Missing atoms for GLU ',I5)
	   endif
	  at1=I+2
	  at2=I+3
	  at3=I+1
	 endif
	do 650 Itmp=1,iproton
	calctemp(Itmp)=0.0
650	continue
	do 667 Icalc=0,1		!loop for two C-O
	at2=at2+Icalc
	CALL VEC(at1,at2,at3,vx,vy,vz,r)
	CALL CENTRE(vz,x(at1),y(at1),z(at1),1.1,cen)
	do 640 K=1,iproton
	  CALL RHAT(xh(k),yh(k),zh(k),cen,rh,rlen)
	  IF(RLEN.GT.12.0) GOTO 640
	  CALL VPROD(rh,vy,pregam,stheta,tempab)
	  CALL VPROD(pregam,vx,yy,sgamma,tempab)
	  calc1=XCO1*((3.0*stheta*stheta)-2.0)
	  calc2=XCO2*(1.0-(3.0*stheta*stheta*sgamma*sgamma))
	  calc3=(calc1+calc2)/(3.0*rlen*rlen*rlen)
	calctemp(K)=calctemp(K)+calc3
640	continue
667	continue
	  do 664 kk=1,iproton
	  shift(kk)=shift(kk)-(calctemp(kk)/2.0)
	  anisco(kk)=anisco(kk)-(calctemp(kk)/2.0)
664	  continue	  
630	continue
9999	continue
C *******NOW DO (O=)C-N TOO *****************************
	type *,' Calculating CN anisotropy...'
	do 130 I=1,iheavyatom
	  if (aname(I).ne.' C  ') goto 130
	    IF (ANAME(I+1).NE.' O  ') then
 	     print 94,IRES(I)
	     goto 130
	    ENDIF
          do 131 j=1,23
            if((aname(I+j).eq.' N  ').and.(ires(I).eq.(ires(i+j)-1)))
     1       then
	     Inext=I+J
	     GOTO 132
	    ENDIF
131	  CONTINUE
	     if(ires(I).ne.ires(iheavyatom)) then
 	     print 194,IRES(I)
194	     format(' N atom not found after residue ',I5)
	     endif
	     goto 130
132	  at1=I	  
	  at2=Inext
	  at3=I+1
	CALL VEC(at1,at2,at3,vx,vy,vz,r)
	rbond=r*0.85	!i.e. 85% along C-N bond
	CALL CENTRE(vz,x(at1),y(at1),z(at1),rbond,cen)
	do 140 K=1,iproton
	  CALL RHAT(xh(k),yh(k),zh(k),cen,rh,rlen)
	  IF(RLEN.GT.12.0) GOTO 140
        if(anameh(K).eq.' H  '.and.(iresh(k).eq.ires(Inext)))
     $    goto 140
	  CALL VPROD(rh,vy,pregam,stheta,tempab)
	  CALL VPROD(pregam,vx,yy,sgamma,tempab)
	  calc1=XCN1*((3.0*stheta*stheta)-2.0)
	  calc2=XCN2*(1.0-(3.0*stheta*stheta*sgamma*sgamma))
	  calc3=(calc1+calc2)/(3.0*rlen*rlen*rlen)
	  shift(k)=shift(k)-calc3	  
	  aniscn(k)=aniscn(k)-calc3
140	continue
130	continue
C *************************************************
C *******Calculate electric field component********
C First find C(alpha)-H(alpha) pair, then work through
C  all C, O and N finding field from them
	type *,' Calculating electric field shift...'
	do 2000 ie=1,iproton
C skip for NH proton
        if(anameh(ie).eq.' H  ')goto 2000
	iha=ie
	Ez=0.0
	do 2010 je=1,iheavyatom	
C Find attached heavy atom
	  if((ires(je).eq.iresh(ie)).and.(aname(je)(3:3).eq.
     1      anameh(ie)(3:3))) then
	  ica=je
	  goto 2020
	  endif
2010	continue
2020	continue
	call VEC2(iha,ica,vca,rca)
C Returns vector vca and length rca between H(ie) and CA(je)
C Now go through each C and add shift due to this
	do 2030 je=1,iheavyatom
	if(ires(je).eq.iresh(iha)) goto 2030
	if(aname(je).eq.' C  ') then
	  call VEC2(iha,je,vc,rcc)
	   if(rcc.gt.6.0) goto 2030	!Ignore for distance>6A
	  call VSCAL(vca,vc,sc,cthet)
	  Efact=Qc/(rcc*rcc)
	  Ez=Ez+(cthet*Efact)
	ELSE if(aname(je).eq.' N  ') then
	  call VEC2(iha,je,vn,rcn)
	   if(rcn.gt.6.0) goto 2030
	  call VSCAL(vca,vn,sc,cthet)
	  Efact=Qn/(rcn*rcn)
	  Ez=Ez+(cthet*Efact)
	ELSE if(aname(je).eq.' O  ') then
	  call VEC2(iha,je,vo,rco)
	   if(rco.gt.6.0) goto 2030
	  call VSCAL(vca,vo,sc,cthet)
	  Efact=Qo/(rco*rco)
	  Ez=Ez+(cthet*Efact)
	ELSE if(aname(je).eq.' H  ') then
	  call VEC2(iha,je,vchn,rchn)
	   if(rchn.gt.6.0) goto 2030
	  call VSCAL(vca,vchn,sc,cthet)
	  Efact=Qhn/(rchn*rchn)
	  Ez=Ez+(cthet*Efact)
	endif
2030	continue
	shift(ie)=shift(ie)+(Ez*sigmaE)
	shiftE(ie)=shiftE(ie)+(Ez*sigmaE)
2000	continue
C ****************END OF CALCULATIONS PROPER********************
C *Correction of Gly HA shift by 0.22 ppm for random coil effect,
C subtract 0.20 from HN shift, subtract 0.65 from HA shift,
C convert to single precision, add random coil shift
	do 2040 kkk=1,iproton
	shift(kkk)=shift(kkk) + rcnet(kkk)
	if(anameh(kkk)(2:3).eq.'HA') shift(kkk)=shift(kkk)-0.65
	if(anameh(kkk)(2:3).eq.'H ') shift(kkk)=shift(kkk)-0.20
	sum(kkk)=rcnet(kkk)+anisco(kkk)+aniscn(kkk)+shiftE(kkk)
	if((resh(kkk).eq.'GLY').and.(anameh(kkk).eq.'1HA '.or.
     1   anameh(kkk).eq.'2HA '.or.anameh(kkk).eq.' HA1'.or.
     2   anameh(kkk).eq.' HA2')) shift(kkk)=shift(kkk)+0.22
	  do 2050 L=1,179
	   IF ((datres(L).eq.resh(kkk)).and.(datnam(L)(2:4).eq.
     1     anameh(kkk)(2:4))) THEN
C This next bit replaces HA shifts for aromatics by 4.45 ppm
C  (gives roughly best fit to data)
	    if((resh(kkk).eq.'TYR'.or.resh(kkk).eq.'PHE'.or.
     1      resh(kkk).eq.'TRP'.or.resh(kkk).eq.'HIS').and.anameh(kkk).
     2      eq.' HA ') then
	      final(kkk)=shift(kkk)+4.45 
	      shift(kkk)=shift(kkk)+4.45-datshift(L)
	    ELSE
	      final(kkk)=shift(kkk) + datshift(L)
	    ENDIF
	  ENDIF
2050	  continue
2040	continue
9990    continue

C ***************** Output SHIFT to channel 3**********************
	if(yn.ne.'Y'.and.yn.ne.'y') goto 2100
        do 909 ii=1,iprot
	DO 910 KK=1,iproton
	IF(IEXP(ii).EQ.IRESH(KK).AND.pexp(ii).eq.prott.and. 
     1   atexp(ii)(1:2).eq.anameh(kk)(2:3)) THEN
	  exptln(kk)=eexp(ii)
C exptln is the experimental shift: now subtract random coil 
	  do 2051 L=1,179
	   IF ((datres(L).eq.resh(kk)).and.(datnam(L)(1:3).eq.
     1     anameh(kk)(1:3))) then
	   exptln(kk)=exptln(kk) - datshift(L)
	   ENDIF
2051	  continue
	  exptl(kk)=exptln(kk)
	  goto 909
	ENDIF
910	continue
909	continue	
        if(imodel.eq.1) goto 501
	write (3,988)
988	format(' Proton name, followed by calc. and observed shift and
     . difference')
	write (3,2052)
2052	format(' (as shift - random coil shift)')
	do 500 iprnt=1,iproton
	  if(exptln(iprnt).ne.0.00)
     1     write(3,991) iprnt,iresh(iprnt),resh(iprnt),anameh(iprnt),
     1     shift(iprnt),exptl(iprnt),(shift(iprnt)-exptl(iprnt))
991	format(I5,I5,1X,A3,2X,A4,3F10.6)
500	continue
501     continue
        if(imodel.eq.1) then
	write(3,775) nmodel 
775     format(' Result for model number',I8)
        endif
	CALL STATS(shift,EXPTL,iproton,anameh)
	type *,'Do you want output sorted? [N]'
	read(5,970)ans
	if(ans.ne.'Y'.and.ans.ne.'y') goto 1000
	itemp=0
	do 505 I=1,iproton
	if(exptl(I).eq.0.0) goto 505
	itemp=itemp+1
	shiftt(itemp)=shift(I)
	exptlt(itemp)=exptl(I)
505	continue
	CALL INDEXX(itemp,shiftt,indx)
	write (3,990) 
C990	format('    Calc       Exptl       Diff')
990	format('    Calc       Diff')
	do 506 I=1,itemp
C	write (3,989) shiftt(indx(I)),exptlt(indx(I)),
	write (3,983) shiftt(indx(I)),
     1   (shiftt(indx(I))-exptlt(indx(I)))
989	format(3F10.4)
983	format(2F10.4)
506	continue
	goto 1000
C **************Normal output starts here***************************
2100	continue
C	write (3,987) FN1
987	format(' Shift calculated for protons from ',A)
C	write (3,986)
986	format(' Added 0.22 for Gly, subtracted 0.65 (HA only), added 
     1random coil shift')
C	write (3,985)
985	format('           Proton       ring     anisCO    anisCN     
     1sigmaE   sum       final')
	do 510 ip=1,iproton
	if(icalculate.eq.1.and.anameh(ip)(2:3).ne.'HA') goto 510
	if(icalculate.eq.2.and.anameh(ip)(2:3).ne.'H ') goto 510
C Now do averaging for Phe,Tyr,methyls
	if((resh(ip).eq.'VAL'.and.(anameh(ip).eq.'1HG1'.or.anameh(ip).
     1eq.'1HG2')).or.(resh(ip).eq.'LEU'.and.(anameh(ip).eq.'1HD1'.or.
     2anameh(ip).eq.'1HD2')).or.(resh(ip).eq.'ILE'.and.(anameh(ip).eq.
     3'1HG2'.or.anameh(ip).eq.'1HD1')).or.(resh(ip).eq.'ALA'.and.
     4anameh(ip).eq.'1HB ').or.(resh(ip).eq.'THR'.and.anameh(ip).eq.
     5'1HG2')) THEN
	  final(ip)=final(ip)-sum(ip)
	  sum(ip)=(sum(ip)+sum(ip+1)+sum(ip+2))/3.0
	  anisco(ip)=(anisco(ip)+anisco(ip+1)+anisco(ip+2))/3.0
	  aniscn(ip)=(aniscn(ip)+aniscn(ip+1)+aniscn(ip+2))/3.0
	  shiftE(ip)=(shiftE(ip)+shiftE(ip+1)+shiftE(ip+2))/3.0
	  shift(ip)=(shift(ip)+shift(ip+1)+shift(ip+2))/3.0
	  final(ip)=final(ip)+sum(ip)
	  final(ip+1)=0.0
	  final(ip+2)=0.0
	endif
	if(resh(ip).eq.'TYR'.and.anameh(ip).eq.' HD1') then
	  if(anameh(ip+1).eq.' HD2') then
	   final(ip)=final(ip)-sum(ip)
	   sum(ip)=(sum(ip)+sum(ip+1))/2.0
	   anisco(ip)=(anisco(ip)+anisco(ip+1))/2.0
	   aniscn(ip)=(aniscn(ip)+aniscn(ip+1))/2.0
	   shiftE(ip)=(shiftE(ip)+shiftE(ip+1))/2.0
	   shift(ip)=(shift(ip)+shift(ip+1))/2.0
	   final(ip)=final(ip)+sum(ip)
	   final(ip+1)=0.0
	  else
	   type 940,resh(ip),iresh(ip)
940	  format(' Ring protons in unexpected order for ',A3,I4)
	  endif
	endif
	if(resh(ip).eq.'TYR'.and.anameh(ip).eq.' HE1') then
	  if(anameh(ip+1).eq.' HE2') then
	   final(ip)=final(ip)-sum(ip)
	   sum(ip)=(sum(ip)+sum(ip+1))/2.0
	   anisco(ip)=(anisco(ip)+anisco(ip+1))/2.0
	   aniscn(ip)=(aniscn(ip)+aniscn(ip+1))/2.0
	   shiftE(ip)=(shiftE(ip)+shiftE(ip+1))/2.0
	   shift(ip)=(shift(ip)+shift(ip+1))/2.0
	   final(ip)=final(ip)+sum(ip)
	   final(ip+1)=0.0
	  else
	   type 940,resh(ip),iresh(ip)
	  endif
	endif
	if(resh(ip).eq.'PHE'.and.anameh(ip).eq.' HD1') then
	  if(anameh(ip+1).eq.' HD2') then
	   final(ip)=final(ip)-sum(ip)
	   sum(ip)=(sum(ip)+sum(ip+1))/2.0
	   anisco(ip)=(anisco(ip)+anisco(ip+1))/2.0
	   aniscn(ip)=(aniscn(ip)+aniscn(ip+1))/2.0
	   shiftE(ip)=(shiftE(ip)+shiftE(ip+1))/2.0
	   shift(ip)=(shift(ip)+shift(ip+1))/2.0
	   final(ip)=final(ip)+sum(ip)
	   final(ip+1)=0.0
	  else
	   type 940,resh(ip),iresh(ip)
	  endif
	endif
	if(resh(ip).eq.'PHE'.and.anameh(ip).eq.' HE1') then
	  if(anameh(ip+1).eq.' HE2') then
	   final(ip)=final(ip)-sum(ip)
	   sum(ip)=(sum(ip)+sum(ip+1))/2.0
	   anisco(ip)=(anisco(ip)+anisco(ip+1))/2.0
	   aniscn(ip)=(aniscn(ip)+aniscn(ip+1))/2.0
	   shiftE(ip)=(shiftE(ip)+shiftE(ip+1))/2.0
	   shift(ip)=(shift(ip)+shift(ip+1))/2.0
	   final(ip)=final(ip)+sum(ip)
	   final(ip+1)=0.0
	  else
	   type 940,resh(ip),iresh(ip)
	  endif
	endif
C ******Write out results*********************
	if(final(ip).ne.0.0) write(3,984)ip,iresh(ip),resh(ip),
     1   anameh(ip),rcnet(ip),anisco(ip),aniscn(ip),shiftE(ip),
     2   shift(ip),final(ip)
984	format(I5,I5,1X,A3,2X,A4,6F10.5)
510	CONTINUE
1000	continue
        imodel=1 !to make it work for single structure
        goto 776
C       Hacked DvdS 10/98
 777	continue
	print *,'Closing unit 3'
	close(3)
C
	STOP
	END
C *******************************************************
	SUBROUTINE VEC(at1,at2,at3,vx,vy,vz,r)
	COMMON X,Y,Z
	DIMENSION X(5000),Y(5000),Z(5000),vx(3),vy(3),vz(3)
	DIMENSION vvx(3),vvy(3),vvz(3)
	INTEGER at1,at2,at3
	x1=x(at1)
	y1=y(at1)
	z1=z(at1)
	x2=x(at2)
	y2=y(at2)
	z2=z(at2)
	x3=x(at3)
	y3=y(at3)
	z3=z(at3)
	vvz(1)=x2-x1
	vvz(2)=y2-y1
	vvz(3)=z2-z1
	CALL UNITV(vvz,vz,r)
	vvz(1)=x3-x1
	vvz(2)=y3-y1
	vvz(3)=z3-z1
	CALL VPROD(vz,vvz,vvy,sthet,tempab)	
	vvy(1)=tempab
	CALL UNITV(vvy,vy,D)
	CALL VPROD(vz,vy,vvx,sthet,tempab)
	CALL UNITV(vvx,vx,D)
	return
	end
C ******************************************************
	SUBROUTINE UNITV(U,U1,D)
	DIMENSION U(3),U1(3)
	D=0.0
	DO 1 JJ=1,3
	D=D+(U(JJ)*U(JJ))
1	CONTINUE
	D=SQRT(D)
	If (D.EQ.0.0) type *,' D is zero in UNITV'
	DO 2 JJ=1,3
	U1(JJ)=U(JJ)/D
2	CONTINUE
	RETURN
	END
C **********************************************************
	SUBROUTINE VPROD(A,B,AB,STHET,tempab)
	DIMENSION A(3),B(3),AB(3)
	AB(2)=A(3)*B(1) - A(1)*B(3)
	AB(1)=A(2)*B(3) - A(3)*B(2)
	AB(3)=A(1)*B(2) - A(2)*B(1)
	tempab=AB(1)
c This is a stupid thing to do, but it's the only way I can get it
c  to behave!
	R1=SQRT(AB(1)*AB(1) + AB(2)*AB(2) + AB(3)*AB(3))
	RA=SQRT(A(1)*A(1) + A(2)*A(2) + A(3)*A(3))
	RB=SQRT(B(1)*B(1) + B(2)*B(2) + B(3)*B(3))
	  IF(RA.EQ.0.0.OR.RB.EQ.0.0) THEN
	    type *,' VPROD Zero divide...ignore',RA,RB
	    STHET=0.0
	  ELSE
	    STHET=R1/(RA*RB)
	  ENDIF
	RETURN
	END
C *******************************************************
	SUBROUTINE CENTRE(vz,a1,a2,a3,dist,cen)
	DIMENSION vz(3),cen(3)
	cen(1)=a1+(vz(1)*dist)
	cen(2)=a2+(vz(2)*dist)
	cen(3)=a3+(vz(3)*dist)
	return
	end
C ********************************************************
	SUBROUTINE RHAT(ax,ay,az,cen,rh,rlen)
	DIMENSION cen(3),rh(3)
	rh(1)=ax-cen(1)
	rh(2)=ay-cen(2)
	rh(3)=az-cen(3)
	rlen=sqrt(rh(1)*rh(1)+rh(2)*rh(2)+rh(3)*rh(3))
	return
	end
C ********************************************************
	SUBROUTINE STATS(SSHIFT,EXPTL,iproton,anameh)
	dimension sshift(3000),exptl(3000),anameh(3000)
	dimension zshift(3000),zexptl(3000)
C Calculates mean and sd of iproton values of shift and exptl
C  also mean and sd of (shift-exptl) and regression coeff
C Values of exptl of 0.00 presumably not found in protha.out
C and are excluded
	itemp=0
	do 5 I=1,iproton
	if(exptl(I).eq.0.0) goto 5
	itemp=itemp+1
	if(anameh(I).eq.'1HA '.and.anameh(I+1).eq.'2HA ') then
	  sshift(I)=(sshift(I)+sshift(I+1))/2.0
	endif
	zshift(itemp)=sshift(I)
	zexptl(itemp)=exptl(I)
5	continue
	sm1=0.0
	sm2=0.0
	sms1=0.0
	sms2=0.0
	sm3=0.0
	sms3=0.0
	sm4=0.0
C sm1,sm2 are accumulated sums of shift and exptl
C sms1,sms2 are accumulated sums of shift**2 and exptl**2
	do 600 L=1,itemp
	sm1=sm1+zshift(L)
	sm2=sm2+zexptl(L)
	sm3=sm3+(zshift(L)-zexptl(L))
600	continue
	sm1=sm1/float(itemp)
	sm2=sm2/float(itemp)
	sm3=sm3/float(itemp)
	do 605 L=1,itemp
	sms1=sms1+((zshift(L)-sm1)*(zshift(L)-sm1))
	sms2=sms2+((zexptl(L)-sm2)*(zexptl(L)-sm2))
	sms3=sms3+((zshift(L)-zexptl(L)-sm3)*(zshift(L)-zexptl(L)-sm3))
605	continue
	sms1=sqrt(sms1/float(itemp))
	sms2=sqrt(sms2/float(itemp))
	sms3=sqrt(sms3/float(itemp))
	do 610 L=1,itemp
	sm4=sm4+((zshift(L)-sm1)*(zexptl(L)-sm2))
610	continue
	sm4=sm4/(float(itemp)*sms1*sms2)
	write (3,1000) sm1,sm2
1000	format(' Means of calc and exptl shifts: ',F6.4,
     1   1X,F6.4)
	write (3,1001) sms1,sms2
1001	format(' SD of calc and exptl shifts: ',F6.4,
     1   1X,F6.4)
	write (3,1002) sm3,sms3
1002	format(' Mean and SD of (calc-exptl shift): ',F6.4,
     1   1X,F6.4)
	write (3,1003) sm4,itemp
1003	format (' Correlation coefficient calc vs exptl: ',F7.5,
     1   ' (',I4,' protons )')
	return
	end
C ************************************************
	SUBROUTINE VEC2(iha,ica,vca,rca)
	DIMENSION vca(3),x(5000),y(5000),z(5000)
	DIMENSION xh(3000),yh(3000),zh(3000)
	COMMON x,y,z,xh,yh,zh
	vca(1)=xh(iha)-x(ica)
	vca(2)=yh(iha)-y(ica)
	vca(3)=zh(iha)-z(ica)
	rca=sqrt((vca(1)*vca(1))+(vca(2)*vca(2))+(vca(3)*vca(3)))
	return
	end
C *******************************************
       SUBROUTINE VSCAL(A,B,SC,CTHET)
       DIMENSION A(3),B(3)
       SC=A(1)*B(1)+A(2)*B(2)+A(3)*B(3)
       RA=SQRT(A(1)*A(1)+A(2)*A(2)+A(3)*A(3))
       RB=SQRT(B(1)*B(1)+B(2)*B(2)+B(3)*B(3))
       CTHET=SC/(RA*RB)
       RETURN
       END
C **********************************************
	SUBROUTINE INDEXX(N,ARRIN,INDX)
C Indexes an array ARRIN of length N, i.e. outputs array INDX
C such that ARRIN(INDX(J)) is in ascending order.
C Taken from 'Numerical Recipes', W.H.Press et al, CUP 1989
	DIMENSION ARRIN(3000),INDX(3000)
	DO 11 J=1,N
		INDX(J)=J
11		CONTINUE
	IF(N.EQ.1) RETURN
	L=N/2+1
	IR=N
10	CONTINUE
	   IF(L.GT.1) THEN
		L=L-1
		INDXT=INDX(L)
		Q=ARRIN(INDXT)
	   ELSE
		INDXT=INDX(IR)
		Q=ARRIN(INDXT)
		INDX(IR)=INDX(1)
		IR=IR-1
		IF(IR.EQ.1) THEN
		  INDX(1)=INDXT
		  RETURN
		ENDIF
	   ENDIF
	   I=L
	   J=L+L
20	   IF(J.LE.IR) THEN
		IF(J.LT.IR) THEN
		  IF(ARRIN(INDX(J)).LT.ARRIN(INDX(J+1)))J=J+1
		ENDIF
		IF(Q.LT.ARRIN(INDX(J))) THEN
		  INDX(I)=INDX(J)
		  I=J
		  J=J+J
		ELSE
		  J=IR+1
		ENDIF
	   GO TO 20
	   ENDIF
	   INDX(I)=INDXT
	GOTO 10
	END
c**********************************************************
      subroutine hm(ip,ir1,ir2,x,y,z,xh,yh,zh,rn,shhm)
c
c  Subroutine Haigh-Mallion:
c  -- given proton ip, and ring atoms ir1 and ir2, with 
c        coordinates x and ring normal vector rn;
c  -- return Haigh-Maillion shift contribution shhm
c
	dimension x(5000),y(5000),z(5000),xh(3000),yh(3000)
      dimension rn(3),r1(3),r2(3),zh(3000)
c
c --- extract coordinates from array x,y,z
c
	r1(1)=x(ir1) - xh(ip)
	r1(2)=y(ir1) - yh(ip)
	r1(3)=z(ir1) - zh(ip)
	r2(1)=x(ir2) - xh(ip)
	r2(2)=y(ir2) - yh(ip)
	r2(3)=z(ir2) - zh(ip)
c
c  -- compute triple scalar product: later versions could
c       hard-code this for efficiency, but this form is 
c       easier to read. (Triple scalar product=vol. of 
c       parallelepiped, =area of desired triangle)
c
      s12 = r1(1)*(r2(2)*rn(3) - r2(3)*rn(2))
     .    + r1(2)*(r2(3)*rn(1) - r2(1)*rn(3))
     .    + r1(3)*(r2(1)*rn(2) - r2(2)*rn(1))
c
c  -- get radial factors
c
      r1sq = r1(1)*r1(1) + r1(2)*r1(2) + r1(3)*r1(3)
      r2sq = r2(1)*r2(1) + r2(2)*r2(2) + r2(3)*r2(3)
      if (r1sq.eq.0.0 .or. r2sq.eq.0.0) then
        write(6,*) 'Geometry error in hm: ',ip,r1,r2
        stop
      end if
      temp = (1./r1sq**1.5 + 1./r2sq**1.5)*0.5	!distance bit
      shhm = s12*temp
      return
      end
C****************************************************************
      subroutine plane (i1, i2, i3, x,y,z, rn, cent)
c**************************************************************
c  --- given three atoms i1,i2 and i3, and coordinates x,y,z
c    - returns rn(i) [i=1,2,3] components of the normalized
c      vector normal to the plane containing the three
c      points and cent(1,2,3), the centre of the three atom.
c
      dimension x(5000),y(5000),z(5000) 
      dimension rn(3),cent(3)
c
        x1 = x(i1)
        y1 = y(i1)
        z1 = z(i1)
        x2 = x(i2)
        y2 = y(i2)
        z2 = z(i2)
        x3 = x(i3)
        y3 = y(i3)
        z3 = z(i3)
c
c       ----- coefficients of the equation for the plane of atoms 1-3
c
        ax = y1*z2 - y2*z1 + y3*z1 - y1*z3 + y2*z3 - y3*z2
        ay = -(x1*z2 - x2*z1 + x3*z1 - x1*z3 + x2*z3 - x3*z2)
        az = x1*y2 - x2*y1 + x3*y1 - x1*y3 + x2*y3 - x3*y2
        anorm = 1./sqrt(ax*ax + ay*ay + az*az)
c
c       ----- normalize to standard form for plane equation (i.e. such
c       ----- that length of the vector "a" is unity
c
        rn(1) = ax*anorm	!ring normal unit vector
        rn(2) = ay*anorm	!=(2-1)x(3-1)
        rn(3) = az*anorm	!divided by its magnitude
c
        cent(1) = (x1+x2+x3)/3.
        cent(2) = (y1+y2+y3)/3.
        cent(3) = (z1+z2+z3)/3.
c
        return
        end
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	subroutine addprot(x,y,z,xh,yh,zh,Iheavyatom,Iproton,
     &   aname,anameh,res,resh,ires,iresh)
C Does not add most exchangeable protons: N-terminal NH3,
C  Glu,Gln,Asp,Asn sidechain, Arg sidechain,Tyr HH, His HD1,
C  Trp HE1, Ser & Thr OH, Cys SH, Lys NH3
C Goes through heavy atom list and adds in simple geometric way
C NB Looks in neighbourhood of attached heavy atom for bonding
C atoms. If the atom ordering is odd, this could cause problems.
C Usually this would mean that the proton does not get added;
C very occasionally it could lead to a proton being put on wrong.
	DIMENSION x(5000),y(5000),z(5000),xh(3000),yh(3000),zh(3000)
	DIMENSION aname(5000),res(5000),ires(5000)
	DIMENSION anameh(3000),resh(3000),iresh(3000)
	CHARACTER aname*4,anameh*4,res*3,resh*3
	iproton=0
C NB This effectively deletes any protons that already were in the file
	do 10 I=1,iheavyatom
	iat2=0
	iat3=0
	iat4=0
	if(aname(I).eq.' N  '.and.ires(I).gt.1.and.res(I).
     &    ne.'PRO') then
	 iat1=I
	 do 11 J=1,13
	 if(aname(I-J).eq.' C  '.and.ires(I-J).eq.ires(I)-1) then
	 iat2=I-J
	 goto 12
         endif
11       continue 
12       continue
         do 13 J=1,15
	 if(aname(I+J).eq.' CA '.and.ires(I+J).eq.ires(I)) then
	 iat3=I+J
	 goto 14
	 endif
13       continue
14       continue
	 if(iat2.gt.0.and.iat3.gt.0) then
	 iproton=iproton+1
	 iheavyatom=iheavyatom+1
	 call addsp2(x,y,z,xh,yh,zh,iproton,iat1,iat2,iat3,1)
         anameh(iproton)=' H  '
         resh(iproton)=res(iat1)
         iresh(iproton)=ires(iat1)
	 x(iheavyatom)=xh(iproton)
	 y(iheavyatom)=yh(iproton)
	 z(iheavyatom)=zh(iproton)
	 aname(iheavyatom)=' H  '
	 res(iheavyatom)=res(iat1)
	 ires(iheavyatom)=ires(iat1)
C The final 1 signals that it is an NH ie use NH bond distance
	 endif
        else if (aname(I).eq.' CA '.and.res(I).ne.'GLY') then
	 iat1=I
         do 21 J=-5,3
	 if(aname(I+J).eq.' N  '.and.ires(I+J).eq.ires(I)) then
	 iat2=I+J
	 goto 22
	 endif
21	 continue
22	 continue
	 do 23 J=1,15
	 if(aname(I+J).eq.' C  '.and.ires(I+J).eq.ires(I)) then
	 iat3=I+J
	 goto 24
	 endif
23	 continue
24	 continue
	 do 25 J=1,15
	 if(aname(I+J).eq.' CB '.and.ires(I+J).eq.ires(I)) then
	 iat4=I+J
	 goto 26
	 endif
25	 continue
26	 continue
         if(iat2.gt.0.and.iat3.gt.0.and.iat4.gt.0) then
	 iproton=iproton+1
	 call addonesp3(x,y,z,xh,yh,zh,iproton,iat1,iat2,iat3,iat4)
	 anameh(iproton)=' HA '
	 resh(iproton)=res(iat1)
	 iresh(iproton)=ires(iat1)
	 endif
	else if((aname(I).eq.' CB '.and.(res(I).eq.'THR'.or.
     &   res(I).eq.'VAL'.or.res(I).eq.'ILE')).or.(aname(I).
     &   eq.' CG '.and.res(I).eq.'LEU')) then
	 iat1=I
	 do 31 J=1,6
	 if(aname(I-J).eq.' CA '.and.ires(I-J).eq.ires(I).
     &    and.res(I).ne.'LEU') then
	 iat2=I-J
	 goto 32
	 endif
31       continue
32       continue
         do 33 J=1,6
         if(res(I-J).eq.'LEU'.and.aname(I-J).eq.' CB ') then
         iat2=I-J
         goto 34
         endif
33       continue
34       continue
         do 35 J=1,3
         if(res(I+J).eq.'THR'.and.aname(I+J)(2:3).eq.'OG') 
     $    iat3=I+J
         if(res(I+J).eq.'THR'.and.aname(I+J)(2:3).eq.'CG')
     $    iat4=I+J
         if(res(I+J).eq.'VAL'.and.aname(I+J).eq.' CG1') iat3=I+J
         if(res(I+J).eq.'VAL'.and.aname(I+J).eq.' CG2') iat4=I+J
         if(res(I+J).eq.'ILE'.and.aname(I+J).eq.' CG1') iat3=I+J
         if(res(I+J).eq.'ILE'.and.aname(I+J).eq.' CG2') iat4=I+J
         if(res(I+J).eq.'LEU'.and.aname(I+J).eq.' CD1') iat3=I+J
         if(res(I+J).eq.'LEU'.and.aname(I+J).eq.' CD2') iat4=I+J
35       continue
         if(iat2.gt.0.and.iat3.gt.0.and.iat4.gt.0) then
	 iproton=iproton+1
	 call addonesp3(x,y,z,xh,yh,zh,iproton,iat1,iat2,iat3,iat4)
         if(res(I).eq.'LEU') then
          anameh(iproton)=' HG '
         else
          anameh(iproton)=' HB '
         endif
	 resh(iproton)=res(iat1)
	 iresh(iproton)=ires(iat1)
         endif
        else if(aname(I).eq.' CA '.and.res(I).eq.'GLY') then
         iat1=I
         do 41 J=-2,2
         if(ires(I+J).eq.ires(I).and.aname(I+J).eq.' N  ')iat3=I+J
         if(ires(I+J).eq.ires(I).and.aname(I+J).eq.' C  ')iat2=I+J
41       continue
         if(iat2.gt.0.and.iat3.gt.0) then
	 iproton=iproton+2
         call addtwosp3(x,y,z,xh,yh,zh,iproton,iat1,iat2,iat3)
	 anameh(iproton-1)='1HA '
	 anameh(iproton)='2HA '
	 resh(iproton-1)=res(iat1)
	 resh(iproton)=res(iat1)
	 iresh(iproton-1)=ires(iat1)
	 iresh(iproton)=ires(iat1)
         endif
        else if(aname(I).eq.' CB '.and.(res(I).eq.'CYS'.or.res(I).
     1   eq.'ASP'.or.res(I).eq.'GLU'.or.res(I).eq.'PHE'.or.res(I).
     2   eq.'HIS'.or.res(I).eq.'LYS'.or.res(I).eq.'LEU'.or.res(I).
     3   eq.'MET'.or.res(I).eq.'ASN'.or.res(I).eq.'PRO'.or.res(I).
     4   eq.'GLN'.or.res(I).eq.'ARG'.or.res(I).eq.'SER'.or.res(I).
     5   eq.'TRP'.or.res(I).eq.'TYR')) then
         iat1=I
         do 42 J=-3,3
         if(ires(I+J).eq.ires(I).and.aname(I+J).eq.' CA ')iat2=I+J
         if(ires(I+J).eq.ires(I).and.aname(I+J)(3:3).eq.'G')
     &    iat3=I+J
42      continue
         if(iat2.gt.0.and.iat3.gt.0) then
	 iproton=iproton+2
         call addtwosp3(x,y,z,xh,yh,zh,iproton,iat1,iat2,iat3)
	 anameh(iproton-1)='1HB '
	 anameh(iproton)='2HB '
	 resh(iproton-1)=res(iat1)
	 resh(iproton)=res(iat1)
	 iresh(iproton-1)=ires(iat1)
	 iresh(iproton)=ires(iat1)
         endif
        else if((aname(I).eq.' CG '.and.(res(I).eq.'GLU'.or.res(I).
     1   eq.'LYS'.or.res(I).eq.'MET'.or.res(I).eq.'PRO'.or.res(I).
     2   eq.'GLN'.or.res(I).eq.'ARG')).or.(aname(I).eq.' CG1'.and.
     3   res(I).eq.'ILE')) then
         iat1=I
         do 43 J=-3,3
         if(ires(I+J).eq.ires(I).and.aname(I+J).eq.' CB ')iat2=I+J
         if(ires(I+J).eq.ires(I).and.aname(I+J)(3:3).eq.'D')
     &    iat3=I+J
43      continue
         if(iat2.gt.0.and.iat3.gt.0) then
	 iproton=iproton+2
         call addtwosp3(x,y,z,xh,yh,zh,iproton,iat1,iat2,iat3)
	 anameh(iproton-1)='1HG '
	 anameh(iproton)='2HG '
         if(res(I).eq.'ILE')anameh(iproton-1)='1HG1'
         if(res(I).eq.'ILE')anameh(iproton)='2HG1'
	 resh(iproton-1)=res(iat1)
	 resh(iproton)=res(iat1)
	 iresh(iproton-1)=ires(iat1)
	 iresh(iproton)=ires(iat1)
         endif
        else if(aname(I).eq.' CD '.and.(res(I).eq.'LYS'.or.res(I).
     1   eq.'ARG'.or.res(I).eq.'PRO')) then
         iat1=I
         do 44 J=-6,3
	 if(ires(I+J).eq.ires(I).and.aname(I+J).eq.' CG ')iat2=I+J
	 if(ires(I+J).eq.ires(I).and.aname(I+J)(3:3).eq.'E')
     &   iat3=I+J 
         if(ires(I+J).eq.ires(I).and.res(I).eq.'PRO'.and.
     &    aname(I+J).eq.' N  ')iat3=I+J
44       continue
         if(iat2.gt.0.and.iat3.gt.0) then
	 iproton=iproton+2
         call addtwosp3(x,y,z,xh,yh,zh,iproton,iat1,iat2,iat3)
	 anameh(iproton-1)='1HD '
	 anameh(iproton)='2HD '
	 resh(iproton-1)=res(iat1)
	 resh(iproton)=res(iat1)
	 iresh(iproton-1)=ires(iat1)
	 iresh(iproton)=ires(iat1)
         endif
        else if(aname(I).eq.' CE '.and.res(I).eq.'LYS') then
         iat1=I
         do 45 J=-3,3
	 if(ires(I+J).eq.ires(I).and.aname(I+J).eq.' CD ')iat2=I+J
	 if(ires(I+J).eq.ires(I).and.aname(I+J)(3:3).eq.'Z')
     &   iat3=I+J 
45       continue
         if(iat2.gt.0.and.iat3.gt.0) then
	 iproton=iproton+2
         call addtwosp3(x,y,z,xh,yh,zh,iproton,iat1,iat2,iat3)
	 anameh(iproton-1)='1HE '
	 anameh(iproton)='2HE '
	 resh(iproton-1)=res(iat1)
	 resh(iproton)=res(iat1)
	 iresh(iproton-1)=ires(iat1)
	 iresh(iproton)=ires(iat1)
         endif
        else if(aname(I).eq.' CB '.and.res(I).eq.'ALA') then
         iat1=I
         do 51 J=1,5
         if(ires(I-J).eq.ires(I).and.aname(I-J).eq.' CA ')iat2=I-J
         if(ires(I-J).eq.ires(I).and.aname(I-J).eq.' N  ')iat3=I-J
51       continue
         if(iat2.gt.0.and.iat3.gt.0)goto 52
       else if((aname(I)(2:3).eq.'CG'.and.(res(I).eq.'VAL'.or.res(I).
     1  eq.'THR')).or.(aname(I).eq.' CG2'.and.res(I).eq.'ILE'))then
        iat1=I
        do 53 J=1,5
        if(ires(I-J).eq.ires(I).and.aname(I-J).eq.' CB ')iat2=I-J
        if(ires(I-J).eq.ires(I).and.aname(I-J).eq.' CA ')iat3=I-J
53      continue
        if(iat2.gt.0.and.iat3.gt.0)goto 52
       else if(aname(I)(2:3).eq.'CD'.and.(res(I).eq.'LEU'.or.res(I).
     1   eq.'ILE'))then
        iat1=I
        do 54 J=1,5
        if(ires(I-J).eq.ires(I).and.aname(I-J)(2:3).eq.'CG')iat2=I-J
        if(ires(I-J).eq.ires(I).and.aname(I-J).eq.' CB ')iat3=I-J
54      continue
        if(iat2.gt.0.and.iat3.gt.0)goto 52
       else if(aname(I)(2:3).eq.'CE'.and.res(I).eq.'MET')then
        iat1=I
        do 55 J=1,5
        if(ires(I-J).eq.ires(I).and.aname(I-J).eq.' SD ')iat2=I-J
        if(ires(I-J).eq.ires(I).and.aname(I-J).eq.' CG ')iat3=I-J
55      continue
        if(iat2.gt.0.and.iat3.gt.0) then
52      iproton=iproton+3
        call addthreesp3(x,y,z,xh,yh,zh,iproton,iat1,iat2,iat3)
        anameh(iproton-2)='1H  '
        anameh(iproton-1)='2H  '
        anameh(iproton)='3H  '
        anameh(iproton-2)(3:4)=aname(iat1)(3:4)
        anameh(iproton-1)(3:4)=aname(iat1)(3:4)
        anameh(iproton)(3:4)=aname(iat1)(3:4)
        resh(iproton-2)=res(iat1)
        resh(iproton-1)=res(iat1)
        resh(iproton)=res(iat1)
	iresh(iproton-2)=ires(iat1)
	iresh(iproton-1)=ires(iat1)
	iresh(iproton)=ires(iat1)
        endif
       else if((res(I).eq.'TYR'.or.res(I).eq.'PHE'.or.res(I).eq.
     1  'TRP').and.aname(I).eq.' CD1')then
        iat1=I
        do 61 J=-5,5
        if(ires(I+J).eq.ires(I).and.aname(I+J).eq.' CG ')iat2=I+J
	if(ires(I+J).eq.ires(I).and.aname(I+J)(3:4).eq.'E1')iat3=I+J
61      continue
        if(iat2.gt.0.and.iat3.gt.0) goto 62
       else if((res(I).eq.'TYR'.or.res(I).eq.'PHE'.or.res(I).eq.
     1  'HIS').and.aname(I).eq.' CD2') then
        iat1=I
        do 63 J=-5,5
        if(ires(I+J).eq.ires(I).and.aname(I+J).eq.' CG ')iat2=I+J
	if(ires(I+J).eq.ires(I).and.aname(I+J)(3:4).eq.'E2')iat3=I+J
63      continue
        if(iat2.gt.0.and.iat3.gt.0) goto 62
       else if((res(I).eq.'TYR'.or.res(I).eq.'PHE').and.aname(I).
     1  eq.' CE1') then
        iat1=I
        do 64 J=-5,5
        if(ires(I+J).eq.ires(I).and.aname(I+J).eq.' CD1')iat2=I+J
	if(ires(I+J).eq.ires(I).and.aname(I+J).eq.' CZ ')iat3=I+J
64      continue
        if(iat2.gt.0.and.iat3.gt.0) goto 62
       else if((res(I).eq.'TYR'.or.res(I).eq.'PHE').and.aname(I).
     1  eq.' CE2') then
        iat1=I
        do 65 J=-5,5
        if(ires(I+J).eq.ires(I).and.aname(I+J).eq.' CD2')iat2=I+J
	if(ires(I+J).eq.ires(I).and.aname(I+J).eq.' CZ ')iat3=I+J
65      continue
        if(iat2.gt.0.and.iat3.gt.0) goto 62
       else if((res(I).eq.'PHE').and.aname(I).eq.' CZ ') then
        iat1=I
        do 66 J=-5,5
        if(ires(I+J).eq.ires(I).and.aname(I+J).eq.' CE1')iat2=I+J
	if(ires(I+J).eq.ires(I).and.aname(I+J).eq.' CE2')iat3=I+J
66      continue
        if(iat2.gt.0.and.iat3.gt.0) goto 62
       else if(res(I).eq.'HIS'.and.aname(I).eq.' CE1') then
        iat1=I
        do 67 J=-5,5
        if(ires(I+J).eq.ires(I).and.aname(I+J).eq.' ND1')iat2=I+J
	if(ires(I+J).eq.ires(I).and.aname(I+J).eq.' NE2')iat3=I+J
67      continue
        if(iat2.gt.0.and.iat3.gt.0) goto 62
       else if(res(I).eq.'TRP'.and.aname(I).eq.' CE3') then
        iat1=I
        do 68 J=-8,8
        if(ires(I+J).eq.ires(I).and.aname(I+J).eq.' CD2')iat2=I+J
	if(ires(I+J).eq.ires(I).and.aname(I+J).eq.' CZ3')iat3=I+J
68      continue
        if(iat2.gt.0.and.iat3.gt.0) goto 62
       else if(res(I).eq.'TRP'.and.aname(I).eq.' CZ3') then
        iat1=I
        do 69 J=-8,8
        if(ires(I+J).eq.ires(I).and.aname(I+J).eq.' CE3')iat2=I+J
	if(ires(I+J).eq.ires(I).and.aname(I+J).eq.' CH2')iat3=I+J
69      continue
        if(iat2.gt.0.and.iat3.gt.0) goto 62
       else if(res(I).eq.'TRP'.and.aname(I).eq.' CH2') then
        iat1=I
        do 70 J=-8,8
        if(ires(I+J).eq.ires(I).and.aname(I+J).eq.' CZ3')iat2=I+J
	if(ires(I+J).eq.ires(I).and.aname(I+J).eq.' CZ2')iat3=I+J
70      continue
        if(iat2.gt.0.and.iat3.gt.0) goto 62
       else if(res(I).eq.'TRP'.and.aname(I).eq.' CZ2') then
        iat1=I
        do 71 J=-8,8
        if(ires(I+J).eq.ires(I).and.aname(I+J).eq.' CH2')iat2=I+J
	if(ires(I+J).eq.ires(I).and.aname(I+J).eq.' CE2')iat3=I+J
71      continue
        if(iat2.gt.0.and.iat3.gt.0) then
62      iproton=iproton+1
        call addsp2(x,y,z,xh,yh,zh,iproton,iat1,iat2,iat3,2)
        anameh(iproton)=' H  '
        anameh(iproton)(3:4)=aname(iat1)(3:4)
        resh(iproton)=res(iat1)
	iresh(iproton)=ires(iat1)
        endif
       else
C presumably a heavy atom that doesn't have protons, or one
C with exchangeable protons that I can't be bothered to do
       endif
10     continue
       return
       end
C===================================================================
	subroutine addsp2(x,y,z,xh,yh,zh,Ih,Iat1,Iat2,Iat3,Ibond)
	DIMENSION x(5000),y(5000),z(5000)
	DIMENSION xh(3000),yh(3000),zh(3000)
	DIMENSION r12(3),u12(3),r13(3),u13(3),r14(3),u14(3)
C If Ibond=1, it's a NH bond, if Ibond=2 it's a CH (aromatic) bond
	if(Ibond.eq.1)blength=0.98
	if(Ibond.eq.2)blength=1.08
	r12(1)=x(Iat1)-x(Iat2)
	r12(2)=y(Iat1)-y(Iat2)
	r12(3)=z(Iat1)-z(Iat2)
	call unitv(r12,u12,d)
	r13(1)=x(Iat1)-x(Iat3)
	r13(2)=y(Iat1)-y(Iat3)
	r13(3)=z(Iat1)-z(Iat3)
	call unitv(r13,u13,d)
	r14(1)=u12(1)+u13(1)
	r14(2)=u12(2)+u13(2)
	r14(3)=u12(3)+u13(3)
	call unitv(r14,u14,d)
	xh(Ih)=x(Iat1)+blength*u14(1)
	yh(Ih)=y(Iat1)+blength*u14(2)
	zh(Ih)=z(Iat1)+blength*u14(3)
	return
	end
C------------------------------------------------------------------
	subroutine addonesp3(x,y,z,xh,yh,zh,Ih,Iat1,Iat2,Iat3,Iat4)
	DIMENSION x(5000),y(5000),z(5000)
	DIMENSION xh(3000),yh(3000),zh(3000)
	DIMENSION r12(3),u12(3),r13(3),u13(3),r14(3),u14(3)
	DIMENSION r15(3),u15(3)
	blength=1.08
	r12(1)=x(Iat1)-x(Iat2)
	r12(2)=y(Iat1)-y(Iat2)
	r12(3)=z(Iat1)-z(Iat2)
	call unitv(r12,u12,d)
	r13(1)=x(Iat1)-x(Iat3)
	r13(2)=y(Iat1)-y(Iat3)
	r13(3)=z(Iat1)-z(Iat3)
	call unitv(r13,u13,d)
	r14(1)=x(Iat1)-x(Iat4)
	r14(2)=y(Iat1)-y(Iat4)
	r14(3)=z(Iat1)-z(Iat4)
	call unitv(r14,u14,d)
	r15(1)=u12(1)+u13(1)+u14(1)
	r15(2)=u12(2)+u13(2)+u14(2)
	r15(3)=u12(3)+u13(3)+u14(3)
	call unitv(r15,u15,d)
	xh(Ih)=x(Iat1)+blength*u15(1)
	yh(Ih)=y(Iat1)+blength*u15(2)
	zh(Ih)=z(Iat1)+blength*u15(3)
	return
	end
C------------------------------------------------------------------
	subroutine addtwosp3(x,y,z,xh,yh,zh,Ih,Iat1,Iat2,Iat3)
C This is based on the GEN routine in VNMR (J. Hoch)
	DIMENSION x(5000),y(5000),z(5000)
	DIMENSION xh(3000),yh(3000),zh(3000)
	DIMENSION r12(3),r13(3),side(3),us(3)
	DIMENSION add(3),uadd(3),sub(3),usub(3)
	blength=1.08
	alpha=54.75*3.14159/180.   !H-C-H=109.5
	r12(1)=x(Iat1)-x(Iat2)
	r12(2)=y(Iat1)-y(Iat2)
	r12(3)=z(Iat1)-z(Iat2)
	r13(1)=x(Iat1)-x(Iat3)
	r13(2)=y(Iat1)-y(Iat3)
	r13(3)=z(Iat1)-z(Iat3)
	sub(1)=r12(1)-r13(1)
	sub(2)=r12(2)-r13(2)
	sub(3)=r12(3)-r13(3)
	call unitv(sub,usub,d)
	add(1)=r12(1)+r13(1)
	add(2)=r12(2)+r13(2)
	add(3)=r12(3)+r13(3)
	call unitv(add,uadd,d)
	call vprod(usub,uadd,side,sintemp,tmp)
	call unitv(side,us,d)
	xh(Ih-1)=x(Iat1)+blength*(uadd(1)*cos(alpha)-us(1)*sin(alpha))
	yh(Ih-1)=y(Iat1)+blength*(uadd(2)*cos(alpha)-us(2)*sin(alpha))
	zh(Ih-1)=z(Iat1)+blength*(uadd(3)*cos(alpha)-us(3)*sin(alpha))
	xh(Ih)=x(Iat1)+blength*(uadd(1)*cos(alpha)+us(1)*sin(alpha))
	yh(Ih)=y(Iat1)+blength*(uadd(2)*cos(alpha)+us(2)*sin(alpha))
	zh(Ih)=z(Iat1)+blength*(uadd(3)*cos(alpha)+us(3)*sin(alpha))
	return
	end
C------------------------------------------------------------------
	subroutine addthreesp3(x,y,z,xh,yh,zh,Ih,Iat1,Iat2,Iat3)
C This is based on the GEN routine in VNMR (J. Hoch)
	DIMENSION x(5000),y(5000),z(5000)
	DIMENSION xh(3000),yh(3000),zh(3000),twist(3),utw(3)
	DIMENSION r12(3),u12(3),r23(3),perp(3),side(3),uside(3)
	REAL perp
	blength=1.08
	beta=70.5*3.14159/180.  !180-109.5
	cosbeta=cos(beta)
	sinbeta=sin(beta)
	cosgam=0.866     !cos(30)
c The methyl protons should be staggered to Iat3
	r12(1)=x(Iat1)-x(Iat2)
	r12(2)=y(Iat1)-y(Iat2)
	r12(3)=z(Iat1)-z(Iat2)
	call unitv(r12,u12,d)
	r23(1)=x(Iat2)-x(Iat3)
	r23(2)=y(Iat2)-y(Iat3)
	r23(3)=z(Iat2)-z(Iat3)
	vert=(r12(1)*r23(1)+r12(2)*r23(2)+r12(3)*r23(3))/d
	perp(1)=u12(1)*vert
	perp(2)=u12(2)*vert
	perp(3)=u12(3)*vert
	side(1)=r23(1)-perp(1)
	side(2)=r23(2)-perp(2)
	side(3)=r23(3)-perp(3)
	call unitv(side,uside,d)
	call vprod(u12,uside,twist,sintemp,tmp)
	call unitv(twist,utw,d)
	u12(1)=u12(1)*cosbeta
	u12(2)=u12(2)*cosbeta
	u12(3)=u12(3)*cosbeta
	xh(Ih-2)=x(Iat1)+blength*(u12(1)+uside(1)*sinbeta)
	yh(Ih-2)=y(Iat1)+blength*(u12(2)+uside(2)*sinbeta)
	zh(Ih-2)=z(Iat1)+blength*(u12(3)+uside(3)*sinbeta)
	xh(Ih-1)=x(Iat1)+blength*(u12(1)-sinbeta*(
     &   uside(1)/2.0-utw(1)*cosgam))
	yh(Ih-1)=y(Iat1)+blength*(u12(2)-sinbeta*(
     &   uside(2)/2.0-utw(2)*cosgam))
	zh(Ih-1)=z(Iat1)+blength*(u12(3)-sinbeta*(
     &   uside(3)/2.0-utw(3)*cosgam))
C NB The /2.0 is actually *sin(30)
	xh(Ih)=x(Iat1)+blength*(u12(1)-sinbeta*(
     &   uside(1)/2.0+utw(1)*cosgam))
	yh(Ih)=y(Iat1)+blength*(u12(2)-sinbeta*(
     &   uside(2)/2.0+utw(2)*cosgam))
	zh(Ih)=z(Iat1)+blength*(u12(3)-sinbeta*(
     &   uside(3)/2.0+utw(3)*cosgam))
	return
	end
