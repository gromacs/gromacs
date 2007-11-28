##
## $Id$
##
## Gromacs 4.0                         Copyright (c) 1991-2003 
## David van der Spoel, Erik Lindahl
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License
## as published by the Free Software Foundation; either version 2
## of the License, or (at your option) any later version.
##
## To help us fund GROMACS development, we humbly ask that you cite
## the research papers on the package. Check out http://www.gromacs.org
## 
## And Hey:
## Gnomes, ROck Monsters And Chili Sauce
##




.globl nb_kernel302_ia32_3dnow
.globl _nb_kernel302_ia32_3dnow
nb_kernel302_ia32_3dnow:        
_nb_kernel302_ia32_3dnow:       
.set nb302_p_nri, 8
.set nb302_iinr, 12
.set nb302_jindex, 16
.set nb302_jjnr, 20
.set nb302_shift, 24
.set nb302_shiftvec, 28
.set nb302_fshift, 32
.set nb302_gid, 36
.set nb302_pos, 40
.set nb302_faction, 44
.set nb302_charge, 48
.set nb302_p_facel, 52
.set nb302_p_krf, 56
.set nb302_p_crf, 60
.set nb302_Vc, 64
.set nb302_type, 68
.set nb302_p_ntype, 72
.set nb302_vdwparam, 76
.set nb302_Vvdw, 80
.set nb302_p_tabscale, 84
.set nb302_VFtab, 88
.set nb302_invsqrta, 92
.set nb302_dvda, 96
.set nb302_p_gbtabscale, 100
.set nb302_GBtab, 104
.set nb302_p_nthreads, 108
.set nb302_count, 112
.set nb302_mtx, 116
.set nb302_outeriter, 120
.set nb302_inneriter, 124
.set nb302_work, 128
                        ## stack offsets for local variables 
.set nb302_is3, 0
.set nb302_ii3, 4
.set nb302_ixO, 8
.set nb302_iyO, 12
.set nb302_izO, 16
.set nb302_ixH, 20
.set nb302_iyH, 28
.set nb302_izH, 36
.set nb302_qqOO, 44
.set nb302_qqOH, 52
.set nb302_qqHH, 60
.set nb302_two, 68
.set nb302_n1, 76
.set nb302_tsc, 84
.set nb302_vctot, 92
.set nb302_innerjjnr, 100
.set nb302_innerk, 104
.set nb302_fixO, 108
.set nb302_fiyO, 112
.set nb302_fizO, 116
.set nb302_fixH, 120
.set nb302_fiyH, 128
.set nb302_fizH, 136
.set nb302_dxO, 144
.set nb302_dyO, 148
.set nb302_dzO, 152
.set nb302_dxH, 156
.set nb302_dyH, 164
.set nb302_dzH, 172
.set nb302_tmprsqH, 180
.set nb302_n, 188                           ## idx for outer loop
.set nb302_nn1, 192                         ## number of outer iterations
.set nb302_nri, 196
.set nb302_nouter, 200
.set nb302_ninner, 204
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $208,%esp          ## local stack space 
        femms
        movl nb302_p_nri(%ebp),%ecx
        movl nb302_p_facel(%ebp),%esi
        movl nb302_p_tabscale(%ebp),%edi
        movl (%ecx),%ecx
        movl %ecx,nb302_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb302_nouter(%esp)
        movl %eax,nb302_ninner(%esp)

        ## assume we have at least one i particle - start directly      

        movl  nb302_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx),%ebx           ## ebx=ii 

        movl  nb302_charge(%ebp),%edx
        movd  (%esi),%mm1       ## mm1=facel 
        movd  (%edx,%ebx,4),%mm2    ## mm2=charge[ii0] (O) 
        movd  4(%edx,%ebx,4),%mm3       ## mm2=charge[ii0+1] (H)  
        movq  %mm2,%mm4
        pfmul %mm1,%mm4
        movq  %mm3,%mm6
        pfmul %mm1,%mm6
        movq  %mm4,%mm5
        pfmul %mm2,%mm4                 ## mm4=qqOO*facel 
        pfmul %mm3,%mm5                 ## mm5=qqOH*facel 
        pfmul %mm3,%mm6                 ## mm6=qqHH*facel 
        punpckldq %mm5,%mm5         ## spread to both halves 
        punpckldq %mm6,%mm6         ## spread to both halves 
        movq  %mm4,nb302_qqOO(%esp)
        movq  %mm5,nb302_qqOH(%esp)
        movq  %mm6,nb302_qqHH(%esp)
        movd  (%edi),%mm3
        punpckldq %mm3,%mm3
        movq  %mm3,nb302_tsc(%esp)
        movl $0x40000000,%eax
        movl %eax,nb302_two(%esp)
        movl %eax,nb302_two+4(%esp)

_nb_kernel302_ia32_3dnow.nb302_threadloop: 
        movl  nb302_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel302_ia32_3dnow.nb302_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel302_ia32_3dnow.nb302_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb302_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb302_n(%esp)
        movl %ebx,nb302_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel302_ia32_3dnow.nb302_outerstart
        jmp _nb_kernel302_ia32_3dnow.nb302_end

_nb_kernel302_ia32_3dnow.nb302_outerstart: 
        ## ebx contains number of outer iterations
        addl nb302_nouter(%esp),%ebx
        movl %ebx,nb302_nouter(%esp)

_nb_kernel302_ia32_3dnow.nb302_outer: 
        movl  nb302_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb302_is3(%esp)      ## store is3 

        movl  nb302_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movq  (%eax,%ebx,4),%mm5        ## move shX/shY to mm5 and shZ to mm6. 
        movd  8(%eax,%ebx,4),%mm6
        movq  %mm5,%mm0
        movq  %mm5,%mm1
        movq  %mm6,%mm2
        punpckldq %mm0,%mm0         ## also expand shX,Y,Z in mm0--mm2. 
        punpckhdq %mm1,%mm1
        punpckldq %mm2,%mm2

        movl  nb302_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx,%esi,4),%ebx    ## ebx=ii 

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb302_pos(%ebp),%eax      ## eax = base of pos[] 

        pfadd (%eax,%ebx,4),%mm5    ## ix = shX + posX (and iy too) 
        movd  8(%eax,%ebx,4),%mm7       ## cant use direct memory add for 4 bytes (iz) 
        movl  %ebx,nb302_ii3(%esp)          ## (use mm7 as temp. storage for iz.) 
        pfadd %mm7,%mm6
        movq  %mm5,nb302_ixO(%esp)
        movq  %mm6,nb302_izO(%esp)

        movd  12(%eax,%ebx,4),%mm3
        movd  16(%eax,%ebx,4),%mm4
        movd  20(%eax,%ebx,4),%mm5
        punpckldq  24(%eax,%ebx,4),%mm3
        punpckldq  28(%eax,%ebx,4),%mm4
        punpckldq  32(%eax,%ebx,4),%mm5    ## coords of H1 in low mm3-mm5, H2 in high 

        pfadd %mm3,%mm0
        pfadd %mm4,%mm1
        pfadd %mm5,%mm2
        movq %mm0,nb302_ixH(%esp)
        movq %mm1,nb302_iyH(%esp)
        movq %mm2,nb302_izH(%esp)

        ## clear vctot and i forces 
        pxor  %mm7,%mm7
        movq  %mm7,nb302_vctot(%esp)
        movq  %mm7,nb302_fixO(%esp)
        movq  %mm7,nb302_fizO(%esp)
        movq  %mm7,nb302_fixH(%esp)
        movq  %mm7,nb302_fiyH(%esp)
        movq  %mm7,nb302_fizH(%esp)

        movl  nb302_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 
        movl  %edx,nb302_innerk(%esp)      ## number of innerloop atoms 
        addl  nb302_ninner(%esp),%edx
        movl  %edx,nb302_ninner(%esp)

        movl  nb302_pos(%ebp),%esi
        movl  nb302_faction(%ebp),%edi
        movl  nb302_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb302_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
_nb_kernel302_ia32_3dnow.nb302_inner_loop: 
        ## a single j particle iteration here - compare with the unrolled code for comments. 
        movl  nb302_innerjjnr(%esp),%eax
        movl  (%eax),%eax       ## eax=jnr offset 
        addl $4,nb302_innerjjnr(%esp)             ## advance pointer 

        leal  (%eax,%eax,2),%eax

        movq  (%esi,%eax,4),%mm0
        movd  8(%esi,%eax,4),%mm1
        ## copy & expand to mm2-mm4 for the H interactions 
        movq  %mm0,%mm2
        movq  %mm0,%mm3
        movq  %mm1,%mm4
        punpckldq %mm2,%mm2
        punpckhdq %mm3,%mm3
        punpckldq %mm4,%mm4

        pfsubr nb302_ixO(%esp),%mm0
        pfsubr nb302_izO(%esp),%mm1

        movq  %mm0,nb302_dxO(%esp)
        pfmul %mm0,%mm0
        movd  %mm1,nb302_dzO(%esp)
        pfmul %mm1,%mm1
        pfacc %mm0,%mm0
        pfadd %mm1,%mm0         ## mm0=rsqO 

        punpckldq %mm2,%mm2
        punpckldq %mm3,%mm3
        punpckldq %mm4,%mm4 ## mm2-mm4 is jx-jz 
        pfsubr nb302_ixH(%esp),%mm2
        pfsubr nb302_iyH(%esp),%mm3
        pfsubr nb302_izH(%esp),%mm4   ## mm2-mm4 is dxH-dzH 

        movq %mm2,nb302_dxH(%esp)
        movq %mm3,nb302_dyH(%esp)
        movq %mm4,nb302_dzH(%esp)
        pfmul %mm2,%mm2
        pfmul %mm3,%mm3
        pfmul %mm4,%mm4

        pfadd %mm2,%mm3
        pfadd %mm4,%mm3         ## mm3=rsqH 
        movq %mm3,nb302_tmprsqH(%esp)

        pfrsqrt %mm0,%mm1

        movq %mm1,%mm2
        pfmul %mm1,%mm1
        pfrsqit1 %mm0,%mm1
        pfrcpit2 %mm2,%mm1      ## mm1=invsqrt  
        pfmul %mm1,%mm0         ## mm0=rsq  

        pfmul nb302_tsc(%esp),%mm0
        pf2iw %mm0,%mm4
        movd %mm4,nb302_n1(%esp)
        pi2fd %mm4,%mm4
        pfsub %mm4,%mm0              ## now mm0 is eps and mm4 n0 
        movq  %mm0,%mm2
        pfmul %mm2,%mm2         ## mm0 is eps, mm2 eps2 

        ## coulomb table 
        movl nb302_VFtab(%ebp),%edx
        movl nb302_n1(%esp),%ecx
        shll $2,%ecx

        ## load all values we need 
        movd (%edx,%ecx,4),%mm4
        movd 4(%edx,%ecx,4),%mm5
        movd 8(%edx,%ecx,4),%mm6
        movd 12(%edx,%ecx,4),%mm7

        pfmul %mm0,%mm6 ## mm6 = Geps           
        pfmul %mm2,%mm7 ## mm7 = Heps2 

        pfadd %mm6,%mm5
        pfadd %mm7,%mm5 ## mm5 = Fp 

        pfmul nb302_two(%esp),%mm7      ## two*Heps2 
        pfadd %mm6,%mm7
        pfadd %mm5,%mm7 ## mm7=FF 

        pfmul %mm0,%mm5 ## mm5=eps*Fp 
        pfadd %mm4,%mm5 ##  mm5= VV 

        pfmul nb302_qqOO(%esp),%mm5     ## vcoul=qq*VV 
        pfmul nb302_qqOO(%esp),%mm7     ## fijC=qq*FF 

        ## update vctot directly, use mm3 for fscal sum. 
        pfadd nb302_vctot(%esp),%mm5
        movq %mm5,nb302_vctot(%esp)
        movq %mm7,%mm3

        ## change sign of fscal and multiply with rinv  
        pxor %mm0,%mm0
        pfsubr %mm0,%mm3
        pfmul nb302_tsc(%esp),%mm3
        pfmul %mm1,%mm3   ## mm3 is total fscal (for the oxygen) now 

        ## Ready with the oxygen - potential is updated, fscal is in mm3. 
        ## time for hydrogens! 

        movq nb302_tmprsqH(%esp),%mm0

        pfrsqrt %mm0,%mm1
        pswapd %mm0,%mm0
        pfrsqrt %mm0,%mm2
        pswapd %mm0,%mm0
        punpckldq %mm2,%mm1     ## seeds are in mm1 now, and rsq in mm0. 

        movq %mm1,%mm2
        pfmul %mm1,%mm1
        pfrsqit1 %mm0,%mm1
        pfrcpit2 %mm2,%mm1      ## mm1=invsqrt 

        pfmul %mm1,%mm0         ## mm0=r 
        pfmul nb302_tsc(%esp),%mm0
        pf2iw %mm0,%mm4
        movq %mm4,nb302_n1(%esp)
        pi2fd %mm4,%mm4
        pfsub %mm4,%mm0              ## now mm0 is eps and mm4 n0 
        movq  %mm0,%mm2
        pfmul %mm2,%mm2         ## mm0 is eps, mm2 eps2 

        ## coulomb table 
        movl nb302_VFtab(%ebp),%edx
        movl nb302_n1(%esp),%ecx
        shll $2,%ecx
        ## load all values we need 
        movd (%edx,%ecx,4),%mm4
        movd 4(%edx,%ecx,4),%mm5
        movd 8(%edx,%ecx,4),%mm6
        movd 12(%edx,%ecx,4),%mm7
        movl nb302_n1+4(%esp),%ecx
        shll $2,%ecx
        punpckldq (%edx,%ecx,4),%mm4
        punpckldq 4(%edx,%ecx,4),%mm5
        punpckldq 8(%edx,%ecx,4),%mm6
        punpckldq 12(%edx,%ecx,4),%mm7

        pfmul %mm0,%mm6 ## mm6 = Geps           
        pfmul %mm2,%mm7 ## mm7 = Heps2 

        pfadd %mm6,%mm5
        pfadd %mm7,%mm5 ## mm5 = Fp 

        pfmul nb302_two(%esp),%mm7      ## two*Heps2 
        pfadd %mm6,%mm7
        pfadd %mm5,%mm7 ## mm7=FF 

        pfmul %mm0,%mm5 ## mm5=eps*Fp 
        pfadd %mm4,%mm5 ##  mm5= VV 

        pfmul nb302_qqOH(%esp),%mm5     ## vcoul=qq*VV 
        pfmul nb302_qqOH(%esp),%mm7     ## fijC=qq*FF 
        ## update vctot 
        pfadd nb302_vctot(%esp),%mm5
        movq %mm5,nb302_vctot(%esp)

        ## change sign of fijC and multiply by rinv 
        pxor %mm4,%mm4
        pfsub %mm7,%mm4
        pfmul nb302_tsc(%esp),%mm4
        pfmul %mm1,%mm4   ## mm4 is total fscal (for the hydrogens) now         

        ## spread oxygen fscalar to both positions 
        punpckldq %mm3,%mm3
        ## calc vectorial force for O 
        movq nb302_dxO(%esp),%mm0
        movd nb302_dzO(%esp),%mm1
        pfmul %mm3,%mm0
        pfmul %mm3,%mm1

        ## calc vectorial force for H's 
        movq nb302_dxH(%esp),%mm5
        movq nb302_dyH(%esp),%mm6
        movq nb302_dzH(%esp),%mm7
        pfmul %mm4,%mm5
        pfmul %mm4,%mm6
        pfmul %mm4,%mm7

        ## update iO particle force 
        movq nb302_fixO(%esp),%mm2
        movd nb302_fizO(%esp),%mm3
        pfadd %mm0,%mm2
        pfadd %mm1,%mm3
        movq %mm2,nb302_fixO(%esp)
        movd %mm3,nb302_fizO(%esp)

        ## update iH forces 
        movq nb302_fixH(%esp),%mm2
        movq nb302_fiyH(%esp),%mm3
        movq nb302_fizH(%esp),%mm4
        pfadd %mm5,%mm2
        pfadd %mm6,%mm3
        pfadd %mm7,%mm4
        movq %mm2,nb302_fixH(%esp)
        movq %mm3,nb302_fiyH(%esp)
        movq %mm4,nb302_fizH(%esp)

        ## pack j forces from H in the same form as the oxygen force. 
        pfacc %mm6,%mm5         ## mm5(l)=fjx(H1+ h2) mm5(h)=fjy(H1+ h2) 
        pfacc %mm7,%mm7         ## mm7(l)=fjz(H1+ h2) 

        pfadd %mm5,%mm0         ## add up total force on j particle.  
        pfadd %mm7,%mm1

        ## update j particle force 
        movq (%edi,%eax,4),%mm2
        movd 8(%edi,%eax,4),%mm3
        pfsub %mm0,%mm2
        pfsub %mm1,%mm3
        movq %mm2,(%edi,%eax,4)
        movd %mm3,8(%edi,%eax,4)

        ## interactions with j H1 

        movq  12(%esi,%eax,4),%mm0
        movd  20(%esi,%eax,4),%mm1
        ## copy & expand to mm2-mm4 for the H interactions 
        movq  %mm0,%mm2
        movq  %mm0,%mm3
        movq  %mm1,%mm4
        punpckldq %mm2,%mm2
        punpckhdq %mm3,%mm3
        punpckldq %mm4,%mm4

        pfsubr nb302_ixO(%esp),%mm0
        pfsubr nb302_izO(%esp),%mm1

        movq  %mm0,nb302_dxO(%esp)
        pfmul %mm0,%mm0
        movd  %mm1,nb302_dzO(%esp)
        pfmul %mm1,%mm1
        pfacc %mm1,%mm0
        pfadd %mm1,%mm0         ## mm0=rsqO 

        punpckldq %mm2,%mm2
        punpckldq %mm3,%mm3
        punpckldq %mm4,%mm4 ## mm2-mm4 is jx-jz 
        pfsubr nb302_ixH(%esp),%mm2
        pfsubr nb302_iyH(%esp),%mm3
        pfsubr nb302_izH(%esp),%mm4   ## mm2-mm4 is dxH-dzH 

        movq %mm2,nb302_dxH(%esp)
        movq %mm3,nb302_dyH(%esp)
        movq %mm4,nb302_dzH(%esp)
        pfmul %mm2,%mm2
        pfmul %mm3,%mm3
        pfmul %mm4,%mm4

        pfadd %mm2,%mm3
        pfadd %mm4,%mm3         ## mm3=rsqH 
        movq %mm3,nb302_tmprsqH(%esp)

        pfrsqrt %mm0,%mm1

        movq %mm1,%mm2
        pfmul %mm1,%mm1
        pfrsqit1 %mm0,%mm1
        pfrcpit2 %mm2,%mm1      ## mm1=invsqrt 
        pfmul %mm1,%mm0         ## mm0=rsq  

        pfmul nb302_tsc(%esp),%mm0
        pf2iw %mm0,%mm4
        movd %mm4,nb302_n1(%esp)
        pi2fd %mm4,%mm4
        pfsub %mm4,%mm0              ## now mm0 is eps and mm4 n0 
        movq  %mm0,%mm2
        pfmul %mm2,%mm2         ## mm0 is eps, mm2 eps2 

        ## coulomb table 
        movl nb302_VFtab(%ebp),%edx
        movl nb302_n1(%esp),%ecx
        shll $2,%ecx

        ## load all values we need 
        movd (%edx,%ecx,4),%mm4
        movd 4(%edx,%ecx,4),%mm5
        movd 8(%edx,%ecx,4),%mm6
        movd 12(%edx,%ecx,4),%mm7

        pfmul %mm0,%mm6 ## mm6 = Geps           
        pfmul %mm2,%mm7 ## mm7 = Heps2 

        pfadd %mm6,%mm5
        pfadd %mm7,%mm5 ## mm5 = Fp 

        pfmul nb302_two(%esp),%mm7      ## two*Heps2 
        pfadd %mm6,%mm7
        pfadd %mm5,%mm7 ## mm7=FF 

        pfmul %mm0,%mm5 ## mm5=eps*Fp 
        pfadd %mm4,%mm5 ##  mm5= VV 

        pfmul nb302_qqOH(%esp),%mm5     ## vcoul=qq*VV 
        pfmul nb302_qqOH(%esp),%mm7     ## fijC=qq*FF 

        ## update vctot  directly, force is moved to mm3 
        pfadd nb302_vctot(%esp),%mm5
        movq %mm5,nb302_vctot(%esp)
        pxor %mm3,%mm3
        pfsub %mm7,%mm3
        pfmul nb302_tsc(%esp),%mm3
        pfmul %mm1,%mm3   ## mm3 is total fscal (for the oxygen) now 

        movq nb302_tmprsqH(%esp),%mm0

        pfrsqrt %mm0,%mm1
        pswapd %mm0,%mm0
        pfrsqrt %mm0,%mm2
        pswapd %mm0,%mm0
        punpckldq %mm2,%mm1     ## seeds are in mm1 now, and rsq in mm0. 

        movq %mm1,%mm2
        pfmul %mm1,%mm1
        pfrsqit1 %mm0,%mm1
        pfrcpit2 %mm2,%mm1      ## mm1=invsqrt 

        pfmul %mm1,%mm0         ## mm0=r 
        pfmul nb302_tsc(%esp),%mm0
        pf2iw %mm0,%mm4
        movq %mm4,nb302_n1(%esp)
        pi2fd %mm4,%mm4
        pfsub %mm4,%mm0              ## now mm0 is eps and mm4 n0 
        movq  %mm0,%mm2
        pfmul %mm2,%mm2         ## mm0 is eps, mm2 eps2 

        ## coulomb table 
        movl nb302_VFtab(%ebp),%edx
        movl nb302_n1(%esp),%ecx
        shll $2,%ecx
        ## load all values we need 
        movd (%edx,%ecx,4),%mm4
        movd 4(%edx,%ecx,4),%mm5
        movd 8(%edx,%ecx,4),%mm6
        movd 12(%edx,%ecx,4),%mm7
        movl nb302_n1+4(%esp),%ecx
        shll $2,%ecx
        punpckldq (%edx,%ecx,4),%mm4
        punpckldq 4(%edx,%ecx,4),%mm5
        punpckldq 8(%edx,%ecx,4),%mm6
        punpckldq 12(%edx,%ecx,4),%mm7


        pfmul %mm0,%mm6 ## mm6 = Geps           
        pfmul %mm2,%mm7 ## mm7 = Heps2 

        pfadd %mm6,%mm5
        pfadd %mm7,%mm5 ## mm5 = Fp 

        pfmul nb302_two(%esp),%mm7      ## two*Heps2 
        pfadd %mm6,%mm7
        pfadd %mm5,%mm7 ## mm7=FF 

        pfmul %mm0,%mm5 ## mm5=eps*Fp 
        pfadd %mm4,%mm5 ##  mm5= VV 

        pfmul nb302_qqHH(%esp),%mm5     ## vcoul=qq*VV 
        pfmul nb302_qqHH(%esp),%mm7     ## fijC=qq*FF 
        ## update vctot 
        pfadd nb302_vctot(%esp),%mm5
        movq %mm5,nb302_vctot(%esp)

        ## change sign of fijC and multiply by rinv 
        pxor %mm4,%mm4
        pfsub %mm7,%mm4
        pfmul nb302_tsc(%esp),%mm4
        pfmul %mm1,%mm4   ## mm4 is total fscal (for the hydrogens) now                 

        ## spread oxygen fscalar to both positions 
        punpckldq %mm3,%mm3
        ## calc vectorial force for O 
        movq nb302_dxO(%esp),%mm0
        movd nb302_dzO(%esp),%mm1
        pfmul %mm3,%mm0
        pfmul %mm3,%mm1

        ## calc vectorial force for H's 
        movq nb302_dxH(%esp),%mm5
        movq nb302_dyH(%esp),%mm6
        movq nb302_dzH(%esp),%mm7
        pfmul %mm4,%mm5
        pfmul %mm4,%mm6
        pfmul %mm4,%mm7

        ## update iO particle force 
        movq nb302_fixO(%esp),%mm2
        movd nb302_fizO(%esp),%mm3
        pfadd %mm0,%mm2
        pfadd %mm1,%mm3
        movq %mm2,nb302_fixO(%esp)
        movd %mm3,nb302_fizO(%esp)

        ## update iH forces 
        movq nb302_fixH(%esp),%mm2
        movq nb302_fiyH(%esp),%mm3
        movq nb302_fizH(%esp),%mm4
        pfadd %mm5,%mm2
        pfadd %mm6,%mm3
        pfadd %mm7,%mm4
        movq %mm2,nb302_fixH(%esp)
        movq %mm3,nb302_fiyH(%esp)
        movq %mm4,nb302_fizH(%esp)

        ## pack j forces from H in the same form as the oxygen force. 
        pfacc %mm6,%mm5         ## mm5(l)=fjx(H1+ h2) mm5(h)=fjy(H1+ h2) 
        pfacc %mm7,%mm7         ## mm7(l)=fjz(H1+ h2) 

        pfadd %mm5,%mm0         ## add up total force on j particle.  
        pfadd %mm7,%mm1

        ## update j particle force 
        movq 12(%edi,%eax,4),%mm2
        movd 20(%edi,%eax,4),%mm3
        pfsub %mm0,%mm2
        pfsub %mm1,%mm3
        movq %mm2,12(%edi,%eax,4)
        movd %mm3,20(%edi,%eax,4)

        ## interactions with j H2 
        movq  24(%esi,%eax,4),%mm0
        movd  32(%esi,%eax,4),%mm1
        ## copy & expand to mm2-mm4 for the H interactions 
        movq  %mm0,%mm2
        movq  %mm0,%mm3
        movq  %mm1,%mm4
        punpckldq %mm2,%mm2
        punpckhdq %mm3,%mm3
        punpckldq %mm4,%mm4

        pfsubr nb302_ixO(%esp),%mm0
        pfsubr nb302_izO(%esp),%mm1

        movq  %mm0,nb302_dxO(%esp)
        pfmul %mm0,%mm0
        movd  %mm1,nb302_dzO(%esp)
        pfmul %mm1,%mm1
        pfacc %mm1,%mm0
        pfadd %mm1,%mm0         ## mm0=rsqO 

        punpckldq %mm2,%mm2
        punpckldq %mm3,%mm3
        punpckldq %mm4,%mm4 ## mm2-mm4 is jx-jz 
        pfsubr nb302_ixH(%esp),%mm2
        pfsubr nb302_iyH(%esp),%mm3
        pfsubr nb302_izH(%esp),%mm4   ## mm2-mm4 is dxH-dzH 

        movq %mm2,nb302_dxH(%esp)
        movq %mm3,nb302_dyH(%esp)
        movq %mm4,nb302_dzH(%esp)
        pfmul %mm2,%mm2
        pfmul %mm3,%mm3
        pfmul %mm4,%mm4

        pfadd %mm2,%mm3
        pfadd %mm4,%mm3         ## mm3=rsqH 
        movq %mm3,nb302_tmprsqH(%esp)

        pfrsqrt %mm0,%mm1

        movq %mm1,%mm2
        pfmul %mm1,%mm1
        pfrsqit1 %mm0,%mm1
        pfrcpit2 %mm2,%mm1      ## mm1=invsqrt 
        pfmul %mm1,%mm0

        pfmul nb302_tsc(%esp),%mm0
        pf2iw %mm0,%mm4
        movd %mm4,nb302_n1(%esp)
        pi2fd %mm4,%mm4
        pfsub %mm4,%mm0              ## now mm0 is eps and mm4 n0 
        movq  %mm0,%mm2
        pfmul %mm2,%mm2         ## mm0 is eps, mm2 eps2 

        ## coulomb table 
        movl nb302_VFtab(%ebp),%edx
        movl nb302_n1(%esp),%ecx
        shll $2,%ecx

        ## load all values we need 
        movd (%edx,%ecx,4),%mm4
        movd 4(%edx,%ecx,4),%mm5
        movd 8(%edx,%ecx,4),%mm6
        movd 12(%edx,%ecx,4),%mm7

        pfmul %mm0,%mm6 ## mm6 = Geps           
        pfmul %mm2,%mm7 ## mm7 = Heps2 

        pfadd %mm6,%mm5
        pfadd %mm7,%mm5 ## mm5 = Fp 

        pfmul nb302_two(%esp),%mm7      ## two*Heps2 
        pfadd %mm6,%mm7
        pfadd %mm5,%mm7 ## mm7=FF 

        pfmul %mm0,%mm5 ## mm5=eps*Fp 
        pfadd %mm4,%mm5 ##  mm5= VV 

        pfmul nb302_qqOH(%esp),%mm5     ## vcoul=qq*VV 
        pfmul nb302_qqOH(%esp),%mm7     ## fijC=qq*FF 

        ## update vctot directly, use mm3 for fscal sum. 
        pfadd nb302_vctot(%esp),%mm5
        movq %mm5,nb302_vctot(%esp)
        pxor %mm3,%mm3
        pfsub %mm7,%mm3
        pfmul nb302_tsc(%esp),%mm3
        pfmul %mm1,%mm3   ## mm3 is total fscal (for the oxygen) now 

        movq nb302_tmprsqH(%esp),%mm0

        pfrsqrt %mm0,%mm1
        pswapd %mm0,%mm0
        pfrsqrt %mm0,%mm2
        pswapd %mm0,%mm0
        punpckldq %mm2,%mm1     ## seeds are in mm1 now, and rsq in mm0. 

        movq %mm1,%mm2
        pfmul %mm1,%mm1
        pfrsqit1 %mm0,%mm1
        pfrcpit2 %mm2,%mm1      ## mm1=invsqrt 

        pfmul %mm1,%mm0         ## mm0=r 
        pfmul nb302_tsc(%esp),%mm0
        pf2iw %mm0,%mm4
        movq %mm4,nb302_n1(%esp)
        pi2fd %mm4,%mm4
        pfsub %mm4,%mm0              ## now mm0 is eps and mm4 n0 
        movq  %mm0,%mm2
        pfmul %mm2,%mm2         ## mm0 is eps, mm2 eps2 

        ## coulomb table 
        movl nb302_VFtab(%ebp),%edx
        movl nb302_n1(%esp),%ecx
        shll $2,%ecx
        ## load all values we need 
        movd (%edx,%ecx,4),%mm4
        movd 4(%edx,%ecx,4),%mm5
        movd 8(%edx,%ecx,4),%mm6
        movd 12(%edx,%ecx,4),%mm7
        movl nb302_n1+4(%esp),%ecx
        shll $2,%ecx
        punpckldq (%edx,%ecx,4),%mm4
        punpckldq 4(%edx,%ecx,4),%mm5
        punpckldq 8(%edx,%ecx,4),%mm6
        punpckldq 12(%edx,%ecx,4),%mm7


        pfmul %mm0,%mm6 ## mm6 = Geps           
        pfmul %mm2,%mm7 ## mm7 = Heps2 

        pfadd %mm6,%mm5
        pfadd %mm7,%mm5 ## mm5 = Fp 

        pfmul nb302_two(%esp),%mm7      ## two*Heps2 
        pfadd %mm6,%mm7
        pfadd %mm5,%mm7 ## mm7=FF 

        pfmul %mm0,%mm5 ## mm5=eps*Fp 
        pfadd %mm4,%mm5 ##  mm5= VV 

        pfmul nb302_qqHH(%esp),%mm5     ## vcoul=qq*VV 
        pfmul nb302_qqHH(%esp),%mm7     ## fijC=qq*FF 
        ## update vctot 
        pfadd nb302_vctot(%esp),%mm5
        movq %mm5,nb302_vctot(%esp)

        ## change sign of fijC and multiply by rinv 
        pxor %mm4,%mm4
        pfsub %mm7,%mm4
        pfmul nb302_tsc(%esp),%mm4
        pfmul %mm1,%mm4   ## mm4 is total fscal (for the hydrogens) now         

        ## spread oxygen fscalar to both positions 
        punpckldq %mm3,%mm3
        ## calc vectorial force for O 
        movq nb302_dxO(%esp),%mm0
        movd nb302_dzO(%esp),%mm1
        pfmul %mm3,%mm0
        pfmul %mm3,%mm1

        ## calc vectorial force for H's 
        movq nb302_dxH(%esp),%mm5
        movq nb302_dyH(%esp),%mm6
        movq nb302_dzH(%esp),%mm7
        pfmul %mm4,%mm5
        pfmul %mm4,%mm6
        pfmul %mm4,%mm7

        ## update iO particle force 
        movq nb302_fixO(%esp),%mm2
        movd nb302_fizO(%esp),%mm3
        pfadd %mm0,%mm2
        pfadd %mm1,%mm3
        movq %mm2,nb302_fixO(%esp)
        movd %mm3,nb302_fizO(%esp)

        ## update iH forces 
        movq nb302_fixH(%esp),%mm2
        movq nb302_fiyH(%esp),%mm3
        movq nb302_fizH(%esp),%mm4
        pfadd %mm5,%mm2
        pfadd %mm6,%mm3
        pfadd %mm7,%mm4
        movq %mm2,nb302_fixH(%esp)
        movq %mm3,nb302_fiyH(%esp)
        movq %mm4,nb302_fizH(%esp)

        ## pack j forces from H in the same form as the oxygen force. 
        pfacc %mm6,%mm5         ## mm5(l)=fjx(H1+ h2) mm5(h)=fjy(H1+ h2) 
        pfacc %mm7,%mm7         ## mm7(l)=fjz(H1+ h2) 

        pfadd %mm5,%mm0         ## add up total force on j particle.  
        pfadd %mm7,%mm1

        ## update j particle force 
        movq 24(%edi,%eax,4),%mm2
        movd 32(%edi,%eax,4),%mm3
        pfsub %mm0,%mm2
        pfsub %mm1,%mm3
        movq %mm2,24(%edi,%eax,4)
        movd %mm3,32(%edi,%eax,4)

        ##  done  - one more? 
        decl nb302_innerk(%esp)
        jz  _nb_kernel302_ia32_3dnow.nb302_updateouterdata
        jmp _nb_kernel302_ia32_3dnow.nb302_inner_loop
_nb_kernel302_ia32_3dnow.nb302_updateouterdata: 
        movl  nb302_ii3(%esp),%ecx

        movq  (%edi,%ecx,4),%mm6       ## increment iO force  
        movd  8(%edi,%ecx,4),%mm7
        pfadd nb302_fixO(%esp),%mm6
        pfadd nb302_fizO(%esp),%mm7
        movq  %mm6,(%edi,%ecx,4)
        movd  %mm7,8(%edi,%ecx,4)

        movq  nb302_fixH(%esp),%mm0
        movq  nb302_fiyH(%esp),%mm3
        movq  nb302_fizH(%esp),%mm1
        movq  %mm0,%mm2
        punpckldq %mm3,%mm0     ## mm0(l)=fxH1, mm0(h)=fyH1 
        punpckhdq %mm3,%mm2     ## mm2(l)=fxH2, mm2(h)=fyH2 
        movq %mm1,%mm3
        pswapd %mm3,%mm3
        ## mm1 is fzH1 
        ## mm3 is fzH2 

        movq  12(%edi,%ecx,4),%mm6          ## increment iH1 force  
        movd  20(%edi,%ecx,4),%mm7
        pfadd %mm0,%mm6
        pfadd %mm1,%mm7
        movq  %mm6,12(%edi,%ecx,4)
        movd  %mm7,20(%edi,%ecx,4)

        movq  24(%edi,%ecx,4),%mm6          ## increment iH2 force 
        movd  32(%edi,%ecx,4),%mm7
        pfadd %mm2,%mm6
        pfadd %mm3,%mm7
        movq  %mm6,24(%edi,%ecx,4)
        movd  %mm7,32(%edi,%ecx,4)


        movl  nb302_fshift(%ebp),%ebx      ## increment fshift force 
        movl  nb302_is3(%esp),%edx

        movq  (%ebx,%edx,4),%mm6
        movd  8(%ebx,%edx,4),%mm7
        pfadd nb302_fixO(%esp),%mm6
        pfadd nb302_fizO(%esp),%mm7
        pfadd %mm0,%mm6
        pfadd %mm1,%mm7
        pfadd %mm2,%mm6
        pfadd %mm3,%mm7
        movq  %mm6,(%ebx,%edx,4)
        movd  %mm7,8(%ebx,%edx,4)

        ## get n from stack
        movl nb302_n(%esp),%esi
        ## get group index for i particle 
        movl  nb302_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        movq  nb302_vctot(%esp),%mm7
        pfacc %mm7,%mm7           ## get and sum the two parts of total potential 

        movl  nb302_Vc(%ebp),%eax
        movd  (%eax,%edx,4),%mm6
        pfadd %mm7,%mm6
        movd  %mm6,(%eax,%edx,4)          ## increment vc[gid] 

        ## finish if last 
        movl nb302_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel302_ia32_3dnow.nb302_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb302_n(%esp)
        jmp _nb_kernel302_ia32_3dnow.nb302_outer
_nb_kernel302_ia32_3dnow.nb302_outerend: 
        ## check if more outer neighborlists remain
        movl  nb302_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel302_ia32_3dnow.nb302_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel302_ia32_3dnow.nb302_threadloop
_nb_kernel302_ia32_3dnow.nb302_end: 
        femms
        movl nb302_nouter(%esp),%eax
        movl nb302_ninner(%esp),%ebx
        movl nb302_outeriter(%ebp),%ecx
        movl nb302_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        addl $208,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret






.globl nb_kernel302nf_ia32_3dnow
.globl _nb_kernel302nf_ia32_3dnow
nb_kernel302nf_ia32_3dnow:      
_nb_kernel302nf_ia32_3dnow:     
.set nb302nf_p_nri, 8
.set nb302nf_iinr, 12
.set nb302nf_jindex, 16
.set nb302nf_jjnr, 20
.set nb302nf_shift, 24
.set nb302nf_shiftvec, 28
.set nb302nf_fshift, 32
.set nb302nf_gid, 36
.set nb302nf_pos, 40
.set nb302nf_faction, 44
.set nb302nf_charge, 48
.set nb302nf_p_facel, 52
.set nb302nf_p_krf, 56
.set nb302nf_p_crf, 60
.set nb302nf_Vc, 64
.set nb302nf_type, 68
.set nb302nf_p_ntype, 72
.set nb302nf_vdwparam, 76
.set nb302nf_Vvdw, 80
.set nb302nf_p_tabscale, 84
.set nb302nf_VFtab, 88
.set nb302nf_invsqrta, 92
.set nb302nf_dvda, 96
.set nb302nf_p_gbtabscale, 100
.set nb302nf_GBtab, 104
.set nb302nf_p_nthreads, 108
.set nb302nf_count, 112
.set nb302nf_mtx, 116
.set nb302nf_outeriter, 120
.set nb302nf_inneriter, 124
.set nb302nf_work, 128
                        ## stack offsets for local variables 
.set nb302nf_is3, 0
.set nb302nf_ii3, 4
.set nb302nf_ixO, 8
.set nb302nf_iyO, 12
.set nb302nf_izO, 16
.set nb302nf_ixH, 20
.set nb302nf_iyH, 28
.set nb302nf_izH, 36
.set nb302nf_qqOO, 44
.set nb302nf_qqOH, 52
.set nb302nf_qqHH, 60
.set nb302nf_n1, 68
.set nb302nf_tsc, 76
.set nb302nf_vctot, 84
.set nb302nf_innerjjnr, 92
.set nb302nf_innerk, 96
.set nb302nf_tmprsqH, 100
.set nb302nf_n, 108                         ## idx for outer loop
.set nb302nf_nn1, 112                       ## number of outer iterations
.set nb302nf_nri, 116
.set nb302nf_nouter, 120
.set nb302nf_ninner, 124
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $128,%esp          ## local stack space 
        femms
        movl nb302nf_p_nri(%ebp),%ecx
        movl nb302nf_p_facel(%ebp),%esi
        movl nb302nf_p_tabscale(%ebp),%edi
        movl (%ecx),%ecx
        movl %ecx,nb302nf_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb302nf_nouter(%esp)
        movl %eax,nb302nf_ninner(%esp)

        ## assume we have at least one i particle - start directly      

        movl  nb302nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx),%ebx           ## ebx=ii 

        movl  nb302nf_charge(%ebp),%edx
        movd  (%esi),%mm1       ## mm1=facel 
        movd  (%edx,%ebx,4),%mm2    ## mm2=charge[ii0] (O) 
        movd  4(%edx,%ebx,4),%mm3       ## mm2=charge[ii0+1] (H)  
        movq  %mm2,%mm4
        pfmul %mm1,%mm4
        movq  %mm3,%mm6
        pfmul %mm1,%mm6
        movq  %mm4,%mm5
        pfmul %mm2,%mm4                 ## mm4=qqOO*facel 
        pfmul %mm3,%mm5                 ## mm5=qqOH*facel 
        pfmul %mm3,%mm6                 ## mm6=qqHH*facel 
        punpckldq %mm5,%mm5         ## spread to both halves 
        punpckldq %mm6,%mm6         ## spread to both halves 
        movq  %mm4,nb302nf_qqOO(%esp)
        movq  %mm5,nb302nf_qqOH(%esp)
        movq  %mm6,nb302nf_qqHH(%esp)
        movd  (%edi),%mm3
        punpckldq %mm3,%mm3
        movq  %mm3,nb302nf_tsc(%esp)
_nb_kernel302nf_ia32_3dnow.nb302nf_threadloop: 
        movl  nb302nf_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel302nf_ia32_3dnow.nb302nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel302nf_ia32_3dnow.nb302nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb302nf_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb302nf_n(%esp)
        movl %ebx,nb302nf_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel302nf_ia32_3dnow.nb302nf_outerstart
        jmp _nb_kernel302nf_ia32_3dnow.nb302nf_end

_nb_kernel302nf_ia32_3dnow.nb302nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb302nf_nouter(%esp),%ebx
        movl %ebx,nb302nf_nouter(%esp)

_nb_kernel302nf_ia32_3dnow.nb302nf_outer: 
        movl  nb302nf_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb302nf_is3(%esp)            ## store is3 

        movl  nb302nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movq  (%eax,%ebx,4),%mm5        ## move shX/shY to mm5 and shZ to mm6. 
        movd  8(%eax,%ebx,4),%mm6
        movq  %mm5,%mm0
        movq  %mm5,%mm1
        movq  %mm6,%mm2
        punpckldq %mm0,%mm0         ## also expand shX,Y,Z in mm0--mm2. 
        punpckhdq %mm1,%mm1
        punpckldq %mm2,%mm2

        movl  nb302nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx,%esi,4),%ebx    ## ebx=ii 

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb302nf_pos(%ebp),%eax      ## eax = base of pos[] 

        pfadd (%eax,%ebx,4),%mm5    ## ix = shX + posX (and iy too) 
        movd  8(%eax,%ebx,4),%mm7       ## cant use direct memory add for 4 bytes (iz) 
        movl  %ebx,nb302nf_ii3(%esp)        ## (use mm7 as temp. storage for iz.) 
        pfadd %mm7,%mm6
        movq  %mm5,nb302nf_ixO(%esp)
        movq  %mm6,nb302nf_izO(%esp)

        movd  12(%eax,%ebx,4),%mm3
        movd  16(%eax,%ebx,4),%mm4
        movd  20(%eax,%ebx,4),%mm5
        punpckldq  24(%eax,%ebx,4),%mm3
        punpckldq  28(%eax,%ebx,4),%mm4
        punpckldq  32(%eax,%ebx,4),%mm5    ## coords of H1 in low mm3-mm5, H2 in high 

        pfadd %mm3,%mm0
        pfadd %mm4,%mm1
        pfadd %mm5,%mm2
        movq %mm0,nb302nf_ixH(%esp)
        movq %mm1,nb302nf_iyH(%esp)
        movq %mm2,nb302nf_izH(%esp)

        ## clear vctot and i forces 
        pxor  %mm7,%mm7
        movq  %mm7,nb302nf_vctot(%esp)

        movl  nb302nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 
        movl  %edx,nb302nf_innerk(%esp)      ## number of innerloop atoms 
        addl  nb302nf_ninner(%esp),%edx
        movl  %edx,nb302nf_ninner(%esp)

        movl  nb302nf_pos(%ebp),%esi
        movl  nb302nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb302nf_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
_nb_kernel302nf_ia32_3dnow.nb302nf_inner_loop: 
        ## a single j particle iteration here - compare with the unrolled code for comments. 
        movl  nb302nf_innerjjnr(%esp),%eax
        movl  (%eax),%eax       ## eax=jnr offset 
        addl $4,nb302nf_innerjjnr(%esp)             ## advance pointer 

        leal  (%eax,%eax,2),%eax

        movq  (%esi,%eax,4),%mm0
        movd  8(%esi,%eax,4),%mm1
        ## copy & expand to mm2-mm4 for the H interactions 
        movq  %mm0,%mm2
        movq  %mm0,%mm3
        movq  %mm1,%mm4
        punpckldq %mm2,%mm2
        punpckhdq %mm3,%mm3
        punpckldq %mm4,%mm4

        pfsubr nb302nf_ixO(%esp),%mm0
        pfsubr nb302nf_izO(%esp),%mm1

        pfmul %mm0,%mm0
        pfmul %mm1,%mm1
        pfacc %mm0,%mm0
        pfadd %mm1,%mm0         ## mm0=rsqO 

        punpckldq %mm2,%mm2
        punpckldq %mm3,%mm3
        punpckldq %mm4,%mm4 ## mm2-mm4 is jx-jz 
        pfsubr nb302nf_ixH(%esp),%mm2
        pfsubr nb302nf_iyH(%esp),%mm3
        pfsubr nb302nf_izH(%esp),%mm4   ## mm2-mm4 is dxH-dzH 

        pfmul %mm2,%mm2
        pfmul %mm3,%mm3
        pfmul %mm4,%mm4

        pfadd %mm2,%mm3
        pfadd %mm4,%mm3         ## mm3=rsqH 
        movq %mm3,nb302nf_tmprsqH(%esp)

        pfrsqrt %mm0,%mm1

        movq %mm1,%mm2
        pfmul %mm1,%mm1
        pfrsqit1 %mm0,%mm1
        pfrcpit2 %mm2,%mm1      ## mm1=invsqrt  
        pfmul %mm1,%mm0         ## mm0=rsq  

        pfmul nb302nf_tsc(%esp),%mm0
        pf2iw %mm0,%mm4
        movd %mm4,nb302nf_n1(%esp)
        pi2fd %mm4,%mm4
        pfsub %mm4,%mm0              ## now mm0 is eps and mm4 n0 
        movq  %mm0,%mm2
        pfmul %mm2,%mm2         ## mm0 is eps, mm2 eps2 

        ## coulomb table 
        movl nb302nf_VFtab(%ebp),%edx
        movl nb302nf_n1(%esp),%ecx
        shll $2,%ecx

        ## load all values we need 
        movd (%edx,%ecx,4),%mm4
        movd 4(%edx,%ecx,4),%mm5
        movd 8(%edx,%ecx,4),%mm6
        movd 12(%edx,%ecx,4),%mm7

        pfmul %mm0,%mm6 ## mm6 = Geps           
        pfmul %mm2,%mm7 ## mm7 = Heps2 

        pfadd %mm6,%mm5
        pfadd %mm7,%mm5 ## mm5 = Fp 

        pfmul %mm0,%mm5 ## mm5=eps*Fp 
        pfadd %mm4,%mm5 ##  mm5= VV 

        pfmul nb302nf_qqOO(%esp),%mm5   ## vcoul=qq*VV 
        ## update vctot directly, use mm3 for fscal sum. 
        pfadd nb302nf_vctot(%esp),%mm5
        movq %mm5,nb302nf_vctot(%esp)

        ## time for hydrogens! 

        movq nb302nf_tmprsqH(%esp),%mm0

        pfrsqrt %mm0,%mm1
        pswapd %mm0,%mm0
        pfrsqrt %mm0,%mm2
        pswapd %mm0,%mm0
        punpckldq %mm2,%mm1     ## seeds are in mm1 now, and rsq in mm0. 

        movq %mm1,%mm2
        pfmul %mm1,%mm1
        pfrsqit1 %mm0,%mm1
        pfrcpit2 %mm2,%mm1      ## mm1=invsqrt 

        pfmul %mm1,%mm0         ## mm0=r 
        pfmul nb302nf_tsc(%esp),%mm0
        pf2iw %mm0,%mm4
        movq %mm4,nb302nf_n1(%esp)
        pi2fd %mm4,%mm4
        pfsub %mm4,%mm0              ## now mm0 is eps and mm4 n0 
        movq  %mm0,%mm2
        pfmul %mm2,%mm2         ## mm0 is eps, mm2 eps2 

        ## coulomb table 
        movl nb302nf_VFtab(%ebp),%edx
        movl nb302nf_n1(%esp),%ecx
        shll $2,%ecx
        ## load all values we need 
        movd (%edx,%ecx,4),%mm4
        movd 4(%edx,%ecx,4),%mm5
        movd 8(%edx,%ecx,4),%mm6
        movd 12(%edx,%ecx,4),%mm7
        movl nb302nf_n1+4(%esp),%ecx
        shll $2,%ecx
        punpckldq (%edx,%ecx,4),%mm4
        punpckldq 4(%edx,%ecx,4),%mm5
        punpckldq 8(%edx,%ecx,4),%mm6
        punpckldq 12(%edx,%ecx,4),%mm7

        pfmul %mm0,%mm6 ## mm6 = Geps           
        pfmul %mm2,%mm7 ## mm7 = Heps2 

        pfadd %mm6,%mm5
        pfadd %mm7,%mm5 ## mm5 = Fp 

        pfmul %mm0,%mm5 ## mm5=eps*Fp 
        pfadd %mm4,%mm5 ##  mm5= VV 

        pfmul nb302nf_qqOH(%esp),%mm5   ## vcoul=qq*VV 
        ## update vctot 
        pfadd nb302nf_vctot(%esp),%mm5
        movq %mm5,nb302nf_vctot(%esp)

        ## interactions with j H1 

        movq  12(%esi,%eax,4),%mm0
        movd  20(%esi,%eax,4),%mm1
        ## copy & expand to mm2-mm4 for the H interactions 
        movq  %mm0,%mm2
        movq  %mm0,%mm3
        movq  %mm1,%mm4
        punpckldq %mm2,%mm2
        punpckhdq %mm3,%mm3
        punpckldq %mm4,%mm4

        pfsubr nb302nf_ixO(%esp),%mm0
        pfsubr nb302nf_izO(%esp),%mm1

        pfmul %mm0,%mm0
        pfmul %mm1,%mm1
        pfacc %mm1,%mm0
        pfadd %mm1,%mm0         ## mm0=rsqO 

        punpckldq %mm2,%mm2
        punpckldq %mm3,%mm3
        punpckldq %mm4,%mm4 ## mm2-mm4 is jx-jz 
        pfsubr nb302nf_ixH(%esp),%mm2
        pfsubr nb302nf_iyH(%esp),%mm3
        pfsubr nb302nf_izH(%esp),%mm4   ## mm2-mm4 is dxH-dzH 

        pfmul %mm2,%mm2
        pfmul %mm3,%mm3
        pfmul %mm4,%mm4

        pfadd %mm2,%mm3
        pfadd %mm4,%mm3         ## mm3=rsqH 
        movq %mm3,nb302nf_tmprsqH(%esp)

        pfrsqrt %mm0,%mm1

        movq %mm1,%mm2
        pfmul %mm1,%mm1
        pfrsqit1 %mm0,%mm1
        pfrcpit2 %mm2,%mm1      ## mm1=invsqrt 
        pfmul %mm1,%mm0         ## mm0=rsq  

        pfmul nb302nf_tsc(%esp),%mm0
        pf2iw %mm0,%mm4
        movd %mm4,nb302nf_n1(%esp)
        pi2fd %mm4,%mm4
        pfsub %mm4,%mm0              ## now mm0 is eps and mm4 n0 
        movq  %mm0,%mm2
        pfmul %mm2,%mm2         ## mm0 is eps, mm2 eps2 

        ## coulomb table 
        movl nb302nf_VFtab(%ebp),%edx
        movl nb302nf_n1(%esp),%ecx
        shll $2,%ecx

        ## load all values we need 
        movd (%edx,%ecx,4),%mm4
        movd 4(%edx,%ecx,4),%mm5
        movd 8(%edx,%ecx,4),%mm6
        movd 12(%edx,%ecx,4),%mm7

        pfmul %mm0,%mm6 ## mm6 = Geps           
        pfmul %mm2,%mm7 ## mm7 = Heps2 

        pfadd %mm6,%mm5
        pfadd %mm7,%mm5 ## mm5 = Fp 

        pfmul %mm0,%mm5 ## mm5=eps*Fp 
        pfadd %mm4,%mm5 ##  mm5= VV 

        pfmul nb302nf_qqOH(%esp),%mm5   ## vcoul=qq*VV 

        ## update vctot  directly, force is moved to mm3 
        pfadd nb302nf_vctot(%esp),%mm5
        movq %mm5,nb302nf_vctot(%esp)

        movq nb302nf_tmprsqH(%esp),%mm0

        pfrsqrt %mm0,%mm1
        pswapd %mm0,%mm0
        pfrsqrt %mm0,%mm2
        pswapd %mm0,%mm0
        punpckldq %mm2,%mm1     ## seeds are in mm1 now, and rsq in mm0. 

        movq %mm1,%mm2
        pfmul %mm1,%mm1
        pfrsqit1 %mm0,%mm1
        pfrcpit2 %mm2,%mm1      ## mm1=invsqrt 

        pfmul %mm1,%mm0         ## mm0=r 
        pfmul nb302nf_tsc(%esp),%mm0
        pf2iw %mm0,%mm4
        movq %mm4,nb302nf_n1(%esp)
        pi2fd %mm4,%mm4
        pfsub %mm4,%mm0              ## now mm0 is eps and mm4 n0 
        movq  %mm0,%mm2
        pfmul %mm2,%mm2         ## mm0 is eps, mm2 eps2 

        ## coulomb table 
        movl nb302nf_VFtab(%ebp),%edx
        movl nb302nf_n1(%esp),%ecx
        shll $2,%ecx
        ## load all values we need 
        movd (%edx,%ecx,4),%mm4
        movd 4(%edx,%ecx,4),%mm5
        movd 8(%edx,%ecx,4),%mm6
        movd 12(%edx,%ecx,4),%mm7
        movl nb302nf_n1+4(%esp),%ecx
        shll $2,%ecx
        punpckldq (%edx,%ecx,4),%mm4
        punpckldq 4(%edx,%ecx,4),%mm5
        punpckldq 8(%edx,%ecx,4),%mm6
        punpckldq 12(%edx,%ecx,4),%mm7


        pfmul %mm0,%mm6 ## mm6 = Geps           
        pfmul %mm2,%mm7 ## mm7 = Heps2 

        pfadd %mm6,%mm5
        pfadd %mm7,%mm5 ## mm5 = Fp 

        pfmul %mm0,%mm5 ## mm5=eps*Fp 
        pfadd %mm4,%mm5 ##  mm5= VV 

        pfmul nb302nf_qqHH(%esp),%mm5   ## vcoul=qq*VV 
        ## update vctot 
        pfadd nb302nf_vctot(%esp),%mm5
        movq %mm5,nb302nf_vctot(%esp)

        ## interactions with j H2 
        movq  24(%esi,%eax,4),%mm0
        movd  32(%esi,%eax,4),%mm1
        ## copy & expand to mm2-mm4 for the H interactions 
        movq  %mm0,%mm2
        movq  %mm0,%mm3
        movq  %mm1,%mm4
        punpckldq %mm2,%mm2
        punpckhdq %mm3,%mm3
        punpckldq %mm4,%mm4

        pfsubr nb302nf_ixO(%esp),%mm0
        pfsubr nb302nf_izO(%esp),%mm1

        pfmul %mm0,%mm0
        pfmul %mm1,%mm1
        pfacc %mm1,%mm0
        pfadd %mm1,%mm0         ## mm0=rsqO 

        punpckldq %mm2,%mm2
        punpckldq %mm3,%mm3
        punpckldq %mm4,%mm4 ## mm2-mm4 is jx-jz 
        pfsubr nb302nf_ixH(%esp),%mm2
        pfsubr nb302nf_iyH(%esp),%mm3
        pfsubr nb302nf_izH(%esp),%mm4   ## mm2-mm4 is dxH-dzH 

        pfmul %mm2,%mm2
        pfmul %mm3,%mm3
        pfmul %mm4,%mm4

        pfadd %mm2,%mm3
        pfadd %mm4,%mm3         ## mm3=rsqH 
        movq %mm3,nb302nf_tmprsqH(%esp)

        pfrsqrt %mm0,%mm1

        movq %mm1,%mm2
        pfmul %mm1,%mm1
        pfrsqit1 %mm0,%mm1
        pfrcpit2 %mm2,%mm1      ## mm1=invsqrt 
        pfmul %mm1,%mm0

        pfmul nb302nf_tsc(%esp),%mm0
        pf2iw %mm0,%mm4
        movd %mm4,nb302nf_n1(%esp)
        pi2fd %mm4,%mm4
        pfsub %mm4,%mm0              ## now mm0 is eps and mm4 n0 
        movq  %mm0,%mm2
        pfmul %mm2,%mm2         ## mm0 is eps, mm2 eps2 

        ## coulomb table 
        movl nb302nf_VFtab(%ebp),%edx
        movl nb302nf_n1(%esp),%ecx
        shll $2,%ecx

        ## load all values we need 
        movd (%edx,%ecx,4),%mm4
        movd 4(%edx,%ecx,4),%mm5
        movd 8(%edx,%ecx,4),%mm6
        movd 12(%edx,%ecx,4),%mm7

        pfmul %mm0,%mm6 ## mm6 = Geps           
        pfmul %mm2,%mm7 ## mm7 = Heps2 

        pfadd %mm6,%mm5
        pfadd %mm7,%mm5 ## mm5 = Fp 

        pfmul %mm0,%mm5 ## mm5=eps*Fp 
        pfadd %mm4,%mm5 ##  mm5= VV 

        pfmul nb302nf_qqOH(%esp),%mm5   ## vcoul=qq*VV 

        ## update vctot directly, use mm3 for fscal sum. 
        pfadd nb302nf_vctot(%esp),%mm5
        movq %mm5,nb302nf_vctot(%esp)

        movq nb302nf_tmprsqH(%esp),%mm0

        pfrsqrt %mm0,%mm1
        pswapd %mm0,%mm0
        pfrsqrt %mm0,%mm2
        pswapd %mm0,%mm0
        punpckldq %mm2,%mm1     ## seeds are in mm1 now, and rsq in mm0. 

        movq %mm1,%mm2
        pfmul %mm1,%mm1
        pfrsqit1 %mm0,%mm1
        pfrcpit2 %mm2,%mm1      ## mm1=invsqrt 

        pfmul %mm1,%mm0         ## mm0=r 
        pfmul nb302nf_tsc(%esp),%mm0
        pf2iw %mm0,%mm4
        movq %mm4,nb302nf_n1(%esp)
        pi2fd %mm4,%mm4
        pfsub %mm4,%mm0              ## now mm0 is eps and mm4 n0 
        movq  %mm0,%mm2
        pfmul %mm2,%mm2         ## mm0 is eps, mm2 eps2 

        ## coulomb table 
        movl nb302nf_VFtab(%ebp),%edx
        movl nb302nf_n1(%esp),%ecx
        shll $2,%ecx
        ## load all values we need 
        movd (%edx,%ecx,4),%mm4
        movd 4(%edx,%ecx,4),%mm5
        movd 8(%edx,%ecx,4),%mm6
        movd 12(%edx,%ecx,4),%mm7
        movl nb302nf_n1+4(%esp),%ecx
        shll $2,%ecx
        punpckldq (%edx,%ecx,4),%mm4
        punpckldq 4(%edx,%ecx,4),%mm5
        punpckldq 8(%edx,%ecx,4),%mm6
        punpckldq 12(%edx,%ecx,4),%mm7


        pfmul %mm0,%mm6 ## mm6 = Geps           
        pfmul %mm2,%mm7 ## mm7 = Heps2 

        pfadd %mm6,%mm5
        pfadd %mm7,%mm5 ## mm5 = Fp 

        pfmul %mm0,%mm5 ## mm5=eps*Fp 
        pfadd %mm4,%mm5 ##  mm5= VV 

        pfmul nb302nf_qqHH(%esp),%mm5   ## vcoul=qq*VV 
        ## update vctot 
        pfadd nb302nf_vctot(%esp),%mm5
        movq %mm5,nb302nf_vctot(%esp)

        ##  done  - one more? 
        decl nb302nf_innerk(%esp)
        jz  _nb_kernel302nf_ia32_3dnow.nb302nf_updateouterdata
        jmp _nb_kernel302nf_ia32_3dnow.nb302nf_inner_loop
_nb_kernel302nf_ia32_3dnow.nb302nf_updateouterdata: 
        ## get n from stack
        movl nb302nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb302nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        movq  nb302nf_vctot(%esp),%mm7
        pfacc %mm7,%mm7           ## get and sum the two parts of total potential 

        movl  nb302nf_Vc(%ebp),%eax
        movd  (%eax,%edx,4),%mm6
        pfadd %mm7,%mm6
        movd  %mm6,(%eax,%edx,4)          ## increment vc[gid] 

        ## finish if last 
        movl nb302nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel302nf_ia32_3dnow.nb302nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb302nf_n(%esp)
        jmp _nb_kernel302nf_ia32_3dnow.nb302nf_outer
_nb_kernel302nf_ia32_3dnow.nb302nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb302nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel302nf_ia32_3dnow.nb302nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel302nf_ia32_3dnow.nb302nf_threadloop
_nb_kernel302nf_ia32_3dnow.nb302nf_end: 
        femms
        movl nb302nf_nouter(%esp),%eax
        movl nb302nf_ninner(%esp),%ebx
        movl nb302nf_outeriter(%ebp),%ecx
        movl nb302nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        addl $128,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret




