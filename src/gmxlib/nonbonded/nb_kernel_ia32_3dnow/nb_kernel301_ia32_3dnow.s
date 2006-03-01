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




.globl nb_kernel301_ia32_3dnow
.globl _nb_kernel301_ia32_3dnow
nb_kernel301_ia32_3dnow:        
_nb_kernel301_ia32_3dnow:       
.set nb301_p_nri, 8
.set nb301_iinr, 12
.set nb301_jindex, 16
.set nb301_jjnr, 20
.set nb301_shift, 24
.set nb301_shiftvec, 28
.set nb301_fshift, 32
.set nb301_gid, 36
.set nb301_pos, 40
.set nb301_faction, 44
.set nb301_charge, 48
.set nb301_p_facel, 52
.set nb301_p_krf, 56
.set nb301_p_crf, 60
.set nb301_Vc, 64
.set nb301_type, 68
.set nb301_p_ntype, 72
.set nb301_vdwparam, 76
.set nb301_Vvdw, 80
.set nb301_p_tabscale, 84
.set nb301_VFtab, 88
.set nb301_invsqrta, 92
.set nb301_dvda, 96
.set nb301_p_gbtabscale, 100
.set nb301_GBtab, 104
.set nb301_p_nthreads, 108
.set nb301_count, 112
.set nb301_mtx, 116
.set nb301_outeriter, 120
.set nb301_inneriter, 124
.set nb301_work, 128
                        ## stack offsets for local variables 
.set nb301_is3, 0
.set nb301_ii3, 4
.set nb301_ixO, 8
.set nb301_iyO, 12
.set nb301_izO, 16
.set nb301_ixH, 20
.set nb301_iyH, 28
.set nb301_izH, 36
.set nb301_iqO, 44
.set nb301_iqH, 52
.set nb301_qqO, 60
.set nb301_qqH, 68
.set nb301_vctot, 76
.set nb301_two, 84
.set nb301_n1, 92
.set nb301_tsc, 100
.set nb301_innerjjnr, 108
.set nb301_innerk, 112
.set nb301_fixO, 116
.set nb301_fiyO, 120
.set nb301_fizO, 124
.set nb301_fixH, 128
.set nb301_fiyH, 136
.set nb301_fizH, 144
.set nb301_dxO, 152
.set nb301_dyO, 156
.set nb301_dzO, 160
.set nb301_dxH, 164
.set nb301_dyH, 172
.set nb301_dzH, 180
.set nb301_tmprsqH, 188
.set nb301_n, 196                           ## idx for outer loop
.set nb301_nn1, 200                         ## number of outer iterations
.set nb301_nri, 204
.set nb301_nouter, 208
.set nb301_ninner, 212
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $216,%esp          ## local stack space 
        femms

        movl nb301_p_nri(%ebp),%ecx
        movl nb301_p_facel(%ebp),%esi
        movl nb301_p_tabscale(%ebp),%edi
        movl (%ecx),%ecx
        movl %ecx,nb301_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb301_nouter(%esp)
        movl %eax,nb301_ninner(%esp)

        movl  nb301_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx),%ebx           ## ebx=ii 

        movl  nb301_charge(%ebp),%edx
        movd  (%esi),%mm1
        movd  (%edx,%ebx,4),%mm2    ## mm2=charge[ii0] 
        pfmul %mm1,%mm2
        movq  %mm2,nb301_iqO(%esp)          ## iqO = facel*charge[ii] 

        movd  4(%edx,%ebx,4),%mm2       ## mm2=charge[ii0+1] 
        pfmul %mm1,%mm2
        punpckldq %mm2,%mm2         ## spread to both halves 
        movq  %mm2,nb301_iqH(%esp)          ## iqH = facel*charge[ii0+1] 

        movd  (%edi),%mm4
        punpckldq %mm4,%mm4         ## spread to both halves 
        movq  %mm4,nb301_tsc(%esp)
        movl $0x40000000,%eax
        movl %eax,nb301_two(%esp)
        movl %eax,nb301_two+4(%esp)
_nb_kernel301_ia32_3dnow.nb301_threadloop: 
        movl  nb301_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel301_ia32_3dnow.nb301_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel301_ia32_3dnow.nb301_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb301_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb301_n(%esp)
        movl %ebx,nb301_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel301_ia32_3dnow.nb301_outerstart
        jmp _nb_kernel301_ia32_3dnow.nb301_end

_nb_kernel301_ia32_3dnow.nb301_outerstart: 
        ## ebx contains number of outer iterations
        addl nb301_nouter(%esp),%ebx
        movl %ebx,nb301_nouter(%esp)

_nb_kernel301_ia32_3dnow.nb301_outer: 
        movl  nb301_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb301_is3(%esp)      ## store is3 

        movl  nb301_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movq  (%eax,%ebx,4),%mm5        ## move shX/shY to mm5 and shZ to mm6. 
        movd  8(%eax,%ebx,4),%mm6
        movq  %mm5,%mm0
        movq  %mm5,%mm1
        movq  %mm6,%mm2
        punpckldq %mm0,%mm0         ## also expand shX,Y,Z in mm0--mm2. 
        punpckhdq %mm1,%mm1
        punpckldq %mm2,%mm2

        movl  nb301_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx,%esi,4),%ebx    ## ebx=ii 

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb301_pos(%ebp),%eax      ## eax = base of pos[] 

        pfadd (%eax,%ebx,4),%mm5    ## ix = shX + posX (and iy too) 
        movd  8(%eax,%ebx,4),%mm7       ## cant use direct memory add for 4 bytes (iz) 
        movl  %ebx,nb301_ii3(%esp)          ## (use mm7 as temp. storage for iz.) 
        pfadd %mm7,%mm6
        movq  %mm5,nb301_ixO(%esp)
        movq  %mm6,nb301_izO(%esp)

        movd  12(%eax,%ebx,4),%mm3
        movd  16(%eax,%ebx,4),%mm4
        movd  20(%eax,%ebx,4),%mm5
        punpckldq  24(%eax,%ebx,4),%mm3
        punpckldq  28(%eax,%ebx,4),%mm4
        punpckldq  32(%eax,%ebx,4),%mm5    ## coords of H1 in low mm3-mm5, H2 in high 

        pfadd %mm3,%mm0
        pfadd %mm4,%mm1
        pfadd %mm5,%mm2
        movq %mm0,nb301_ixH(%esp)
        movq %mm1,nb301_iyH(%esp)
        movq %mm2,nb301_izH(%esp)

        ## clear vctot and i forces 
        pxor  %mm7,%mm7
        movq  %mm7,nb301_vctot(%esp)
        movq  %mm7,nb301_fixO(%esp)
        movd  %mm7,nb301_fizO(%esp)
        movq  %mm7,nb301_fixH(%esp)
        movq  %mm7,nb301_fiyH(%esp)
        movq  %mm7,nb301_fizH(%esp)

        movl  nb301_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 
        movl  %edx,nb301_innerk(%esp)
        addl  nb301_ninner(%esp),%edx
        movl  %edx,nb301_ninner(%esp)

        movl  nb301_pos(%ebp),%esi
        movl  nb301_faction(%ebp),%edi
        movl  nb301_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb301_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
_nb_kernel301_ia32_3dnow.nb301_inner_loop: 
        ## a single j particle iteration 
        movl  nb301_innerjjnr(%esp),%eax
        movl  (%eax),%eax        ## eax=jnr offset 
        addl $4,nb301_innerjjnr(%esp)             ## advance pointer 
        prefetch 16(%ecx)          ## prefetch data - trial and error says 16 is best 

        movl nb301_charge(%ebp),%ecx
        movd (%ecx,%eax,4),%mm7
        punpckldq %mm7,%mm7
        movq %mm7,%mm6
        pfmul nb301_iqO(%esp),%mm6
        pfmul nb301_iqH(%esp),%mm7       ## mm6=qqO, mm7=qqH 
        movd %mm6,nb301_qqO(%esp)
        movq %mm7,nb301_qqH(%esp)

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

        pfsubr nb301_ixO(%esp),%mm0
        pfsubr nb301_izO(%esp),%mm1

        movq  %mm0,nb301_dxO(%esp)
        pfmul %mm0,%mm0
        movd  %mm1,nb301_dzO(%esp)
        pfmul %mm1,%mm1
        pfacc %mm1,%mm0
        pfadd %mm1,%mm0         ## mm0=rsqO 

        punpckldq %mm2,%mm2
        punpckldq %mm3,%mm3
        punpckldq %mm4,%mm4 ## mm2-mm4 is jx-jz 
        pfsubr nb301_ixH(%esp),%mm2
        pfsubr nb301_iyH(%esp),%mm3
        pfsubr nb301_izH(%esp),%mm4   ## mm2-mm4 is dxH-dzH 

        movq %mm2,nb301_dxH(%esp)
        movq %mm3,nb301_dyH(%esp)
        movq %mm4,nb301_dzH(%esp)
        pfmul %mm2,%mm2
        pfmul %mm3,%mm3
        pfmul %mm4,%mm4

        pfadd %mm2,%mm3
        pfadd %mm4,%mm3         ## mm3=rsqH 
        movq %mm3,nb301_tmprsqH(%esp)

        pfrsqrt %mm0,%mm1

        movq %mm1,%mm2
        pfmul %mm1,%mm1
        pfrsqit1 %mm0,%mm1
        pfrcpit2 %mm2,%mm1      ## mm1=invsqrt 

        pfmul %mm1,%mm0         ## mm0=r 

        pfmul nb301_tsc(%esp),%mm0
        pf2iw %mm0,%mm4
        movd %mm4,nb301_n1(%esp)
        pi2fd %mm4,%mm4
        pfsub %mm4,%mm0              ## now mm0 is eps and mm4 n0 
        movq  %mm0,%mm2
        pfmul %mm2,%mm2         ## mm0 is eps, mm2 eps2 

        ## coulomb table 
        movl nb301_VFtab(%ebp),%edx
        movl nb301_n1(%esp),%ecx
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

        pfmul nb301_two(%esp),%mm7      ## two*Heps2 
        pfadd %mm6,%mm7
        pfadd %mm5,%mm7 ## mm7=FF 

        pfmul %mm0,%mm5 ## mm5=eps*Fp 
        pfadd %mm4,%mm5 ##  mm5= VV 

        pfmul nb301_qqO(%esp),%mm5      ## vcoul=qq*VV 
        pfmul nb301_qqO(%esp),%mm7      ## fijC=qq*FF 
        ## update vctot directly, use mm3 for fscal sum. 
        pfadd nb301_vctot(%esp),%mm5
        movq %mm5,nb301_vctot(%esp)
        movq %mm7,%mm3

        ## change sign of fscal and multiply with rinv  
        pxor %mm0,%mm0
        pfsubr %mm0,%mm3
        pfmul nb301_tsc(%esp),%mm3
        pfmul %mm1,%mm3   ## mm3 is total fscal (for the oxygen) now    

        ## Ready with the oxygen - potential is updated, fscal is in mm3. 
        ## now do the two hydrogens. 

        movq nb301_tmprsqH(%esp),%mm0   ## mm0=rsqH 

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
        pfmul nb301_tsc(%esp),%mm0
        pf2iw %mm0,%mm4
        movq %mm4,nb301_n1(%esp)
        pi2fd %mm4,%mm4
        pfsub %mm4,%mm0              ## now mm0 is eps and mm4 n0 
        movq  %mm0,%mm2
        pfmul %mm2,%mm2         ## mm0 is eps, mm2 eps2 

        ## coulomb table 
        movl nb301_VFtab(%ebp),%edx
        movl nb301_n1(%esp),%ecx
        shll $2,%ecx
        ## load all values we need 
        movd (%edx,%ecx,4),%mm4
        movd 4(%edx,%ecx,4),%mm5
        movd 8(%edx,%ecx,4),%mm6
        movd 12(%edx,%ecx,4),%mm7
        movl nb301_n1+4(%esp),%ecx
        shll $2,%ecx
        punpckldq (%edx,%ecx,4),%mm4
        punpckldq 4(%edx,%ecx,4),%mm5
        punpckldq 8(%edx,%ecx,4),%mm6
        punpckldq 12(%edx,%ecx,4),%mm7

        pfmul %mm0,%mm6 ## mm6 = Geps           
        pfmul %mm2,%mm7 ## mm7 = Heps2 

        pfadd %mm6,%mm5
        pfadd %mm7,%mm5 ## mm5 = Fp 

        pfmul nb301_two(%esp),%mm7      ## two*Heps2 
        pfadd %mm6,%mm7
        pfadd %mm5,%mm7 ## mm7=FF 

        pfmul %mm0,%mm5 ## mm5=eps*Fp 
        pfadd %mm4,%mm5 ##  mm5= VV 

        pfmul nb301_qqH(%esp),%mm5      ## vcoul=qq*VV 
        pfmul nb301_qqH(%esp),%mm7      ## fijC=qq*FF 

        ## update vctot 
        pfadd nb301_vctot(%esp),%mm5
        movq %mm5,nb301_vctot(%esp)

        ## change sign of fijC and multiply by rinv 
        pxor %mm4,%mm4
        pfsub %mm7,%mm4
        pfmul nb301_tsc(%esp),%mm4
        pfmul %mm1,%mm4   ## mm4 is total fscal (for the hydrogens) now         

        ## spread oxygen fscalar to both positions 
        punpckldq %mm3,%mm3
        ## calc vectorial force for O 
        prefetchw (%edi,%eax,4) ## prefetch faction to cache  
        movq nb301_dxO(%esp),%mm0
        movd nb301_dzO(%esp),%mm1
        pfmul %mm3,%mm0
        pfmul %mm3,%mm1

        ## calc vectorial force for H's 
        movq nb301_dxH(%esp),%mm5
        movq nb301_dyH(%esp),%mm6
        movq nb301_dzH(%esp),%mm7
        pfmul %mm4,%mm5
        pfmul %mm4,%mm6
        pfmul %mm4,%mm7

        ## update iO particle force 
        movq nb301_fixO(%esp),%mm2
        movd nb301_fizO(%esp),%mm3
        pfadd %mm0,%mm2
        pfadd %mm1,%mm3
        movq %mm2,nb301_fixO(%esp)
        movd %mm3,nb301_fizO(%esp)

        ## update iH forces 
        movq nb301_fixH(%esp),%mm2
        movq nb301_fiyH(%esp),%mm3
        movq nb301_fizH(%esp),%mm4
        pfadd %mm5,%mm2
        pfadd %mm6,%mm3
        pfadd %mm7,%mm4
        movq %mm2,nb301_fixH(%esp)
        movq %mm3,nb301_fiyH(%esp)
        movq %mm4,nb301_fizH(%esp)

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

        ##  done  - one more? 
        decl nb301_innerk(%esp)
        jz  _nb_kernel301_ia32_3dnow.nb301_updateouterdata
        jmp _nb_kernel301_ia32_3dnow.nb301_inner_loop
_nb_kernel301_ia32_3dnow.nb301_updateouterdata: 
        movl  nb301_ii3(%esp),%ecx

        movq  (%edi,%ecx,4),%mm6       ## increment iO force  
        movd  8(%edi,%ecx,4),%mm7
        pfadd nb301_fixO(%esp),%mm6
        pfadd nb301_fizO(%esp),%mm7
        movq  %mm6,(%edi,%ecx,4)
        movd  %mm7,8(%edi,%ecx,4)

        movq  nb301_fixH(%esp),%mm0
        movq  nb301_fiyH(%esp),%mm3
        movq  nb301_fizH(%esp),%mm1
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


        movl  nb301_fshift(%ebp),%ebx      ## increment fshift force 
        movl  nb301_is3(%esp),%edx

        movq  (%ebx,%edx,4),%mm6
        movd  8(%ebx,%edx,4),%mm7
        pfadd nb301_fixO(%esp),%mm6
        pfadd nb301_fizO(%esp),%mm7
        pfadd %mm0,%mm6
        pfadd %mm1,%mm7
        pfadd %mm2,%mm6
        pfadd %mm3,%mm7
        movq  %mm6,(%ebx,%edx,4)
        movd  %mm7,8(%ebx,%edx,4)

        ## get n from stack
        movl nb301_n(%esp),%esi
        ## get group index for i particle 
        movl  nb301_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        movq  nb301_vctot(%esp),%mm7
        pfacc %mm7,%mm7           ## get and sum the two parts of total potential 

        movl  nb301_Vc(%ebp),%eax
        movd  (%eax,%edx,4),%mm6
        pfadd %mm7,%mm6
        movd  %mm6,(%eax,%edx,4)          ## increment vc[gid] 

        ## finish if last 
        movl nb301_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jecxz _nb_kernel301_ia32_3dnow.nb301_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb301_n(%esp)
        jmp _nb_kernel301_ia32_3dnow.nb301_outer
_nb_kernel301_ia32_3dnow.nb301_outerend: 
        ## check if more outer neighborlists remain
        movl  nb301_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jecxz _nb_kernel301_ia32_3dnow.nb301_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel301_ia32_3dnow.nb301_threadloop
_nb_kernel301_ia32_3dnow.nb301_end: 
        femms

        movl nb301_nouter(%esp),%eax
        movl nb301_ninner(%esp),%ebx
        movl nb301_outeriter(%ebp),%ecx
        movl nb301_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        addl $216,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret




.globl nb_kernel301nf_ia32_3dnow
.globl _nb_kernel301nf_ia32_3dnow
nb_kernel301nf_ia32_3dnow:      
_nb_kernel301nf_ia32_3dnow:     
.set nb301nf_p_nri, 8
.set nb301nf_iinr, 12
.set nb301nf_jindex, 16
.set nb301nf_jjnr, 20
.set nb301nf_shift, 24
.set nb301nf_shiftvec, 28
.set nb301nf_fshift, 32
.set nb301nf_gid, 36
.set nb301nf_pos, 40
.set nb301nf_faction, 44
.set nb301nf_charge, 48
.set nb301nf_p_facel, 52
.set nb301nf_p_krf, 56
.set nb301nf_p_crf, 60
.set nb301nf_Vc, 64
.set nb301nf_type, 68
.set nb301nf_p_ntype, 72
.set nb301nf_vdwparam, 76
.set nb301nf_Vvdw, 80
.set nb301nf_p_tabscale, 84
.set nb301nf_VFtab, 88
.set nb301nf_invsqrta, 92
.set nb301nf_dvda, 96
.set nb301nf_p_gbtabscale, 100
.set nb301nf_GBtab, 104
.set nb301nf_p_nthreads, 108
.set nb301nf_count, 112
.set nb301nf_mtx, 116
.set nb301nf_outeriter, 120
.set nb301nf_inneriter, 124
.set nb301nf_work, 128
                        ## stack offsets for local variables 
.set nb301nf_is3, 0
.set nb301nf_ii3, 4
.set nb301nf_ixO, 8
.set nb301nf_iyO, 12
.set nb301nf_izO, 16
.set nb301nf_ixH, 20
.set nb301nf_iyH, 28
.set nb301nf_izH, 36
.set nb301nf_iqO, 44
.set nb301nf_iqH, 52
.set nb301nf_qqO, 60
.set nb301nf_qqH, 68
.set nb301nf_vctot, 76
.set nb301nf_n1, 84
.set nb301nf_tsc, 92
.set nb301nf_innerjjnr, 100
.set nb301nf_innerk, 104
.set nb301nf_tmprsqH, 108
.set nb301nf_n, 116                         ## idx for outer loop
.set nb301nf_nn1, 120                       ## number of outer iterations
.set nb301nf_nri, 124
.set nb301nf_nouter, 128
.set nb301nf_ninner, 132
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $136,%esp          ## local stack space 
        femms

        movl nb301nf_p_nri(%ebp),%ecx
        movl nb301nf_p_facel(%ebp),%esi
        movl nb301nf_p_tabscale(%ebp),%edi
        movl (%ecx),%ecx
        movl %ecx,nb301nf_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb301nf_nouter(%esp)
        movl %eax,nb301nf_ninner(%esp)

        movl  nb301nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx),%ebx           ## ebx=ii 

        movl  nb301nf_charge(%ebp),%edx
        movd  (%esi),%mm1
        movd  (%edx,%ebx,4),%mm2    ## mm2=charge[ii0] 
        pfmul %mm1,%mm2
        movq  %mm2,nb301nf_iqO(%esp)        ## iqO = facel*charge[ii] 

        movd  4(%edx,%ebx,4),%mm2       ## mm2=charge[ii0+1] 
        pfmul %mm1,%mm2
        punpckldq %mm2,%mm2         ## spread to both halves 
        movq  %mm2,nb301nf_iqH(%esp)        ## iqH = facel*charge[ii0+1] 

        movd  (%edi),%mm4
        punpckldq %mm4,%mm4         ## spread to both halves 
        movq  %mm4,nb301nf_tsc(%esp)
        ## assume we have at least one i particle - start directly       
_nb_kernel301nf_ia32_3dnow.nb301nf_threadloop: 
        movl  nb301nf_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel301nf_ia32_3dnow.nb301nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel301nf_ia32_3dnow.nb301nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb301nf_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb301nf_n(%esp)
        movl %ebx,nb301nf_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel301nf_ia32_3dnow.nb301nf_outerstart
        jmp _nb_kernel301nf_ia32_3dnow.nb301nf_end

_nb_kernel301nf_ia32_3dnow.nb301nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb301nf_nouter(%esp),%ebx
        movl %ebx,nb301nf_nouter(%esp)

_nb_kernel301nf_ia32_3dnow.nb301nf_outer: 
        movl  nb301nf_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb301nf_is3(%esp)            ## store is3 

        movl  nb301nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movq  (%eax,%ebx,4),%mm5        ## move shX/shY to mm5 and shZ to mm6. 
        movd  8(%eax,%ebx,4),%mm6
        movq  %mm5,%mm0
        movq  %mm5,%mm1
        movq  %mm6,%mm2
        punpckldq %mm0,%mm0         ## also expand shX,Y,Z in mm0--mm2. 
        punpckhdq %mm1,%mm1
        punpckldq %mm2,%mm2

        movl  nb301nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx,%esi,4),%ebx    ## ebx=ii 

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb301nf_pos(%ebp),%eax      ## eax = base of pos[] 

        pfadd (%eax,%ebx,4),%mm5    ## ix = shX + posX (and iy too) 
        movd  8(%eax,%ebx,4),%mm7       ## cant use direct memory add for 4 bytes (iz) 
        movl  %ebx,nb301nf_ii3(%esp)        ## (use mm7 as temp. storage for iz.) 
        pfadd %mm7,%mm6
        movq  %mm5,nb301nf_ixO(%esp)
        movq  %mm6,nb301nf_izO(%esp)

        movd  12(%eax,%ebx,4),%mm3
        movd  16(%eax,%ebx,4),%mm4
        movd  20(%eax,%ebx,4),%mm5
        punpckldq  24(%eax,%ebx,4),%mm3
        punpckldq  28(%eax,%ebx,4),%mm4
        punpckldq  32(%eax,%ebx,4),%mm5    ## coords of H1 in low mm3-mm5, H2 in high 

        pfadd %mm3,%mm0
        pfadd %mm4,%mm1
        pfadd %mm5,%mm2
        movq %mm0,nb301nf_ixH(%esp)
        movq %mm1,nb301nf_iyH(%esp)
        movq %mm2,nb301nf_izH(%esp)

        ## clear vctot and i forces 
        pxor  %mm7,%mm7
        movq  %mm7,nb301nf_vctot(%esp)

        movl  nb301nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 
        movl  %edx,nb301nf_innerk(%esp)
        addl  nb301nf_ninner(%esp),%edx
        movl  %edx,nb301nf_ninner(%esp)

        movl  nb301nf_pos(%ebp),%esi
        movl  nb301nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb301nf_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
_nb_kernel301nf_ia32_3dnow.nb301nf_inner_loop: 
        ## a single j particle iteration 
        movl  nb301nf_innerjjnr(%esp),%eax
        movl  (%eax),%eax        ## eax=jnr offset 
        addl $4,nb301nf_innerjjnr(%esp)             ## advance pointer 
        prefetch 16(%ecx)          ## prefetch data - trial and error says 16 is best 

        movl nb301nf_charge(%ebp),%ecx
        movd (%ecx,%eax,4),%mm7
        punpckldq %mm7,%mm7
        movq %mm7,%mm6
        pfmul nb301nf_iqO(%esp),%mm6
        pfmul nb301nf_iqH(%esp),%mm7     ## mm6=qqO, mm7=qqH 
        movd %mm6,nb301nf_qqO(%esp)
        movq %mm7,nb301nf_qqH(%esp)

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

        pfsubr nb301nf_ixO(%esp),%mm0
        pfsubr nb301nf_izO(%esp),%mm1

        pfmul %mm0,%mm0
        pfmul %mm1,%mm1
        pfacc %mm1,%mm0
        pfadd %mm1,%mm0         ## mm0=rsqO 

        punpckldq %mm2,%mm2
        punpckldq %mm3,%mm3
        punpckldq %mm4,%mm4 ## mm2-mm4 is jx-jz 
        pfsubr nb301nf_ixH(%esp),%mm2
        pfsubr nb301nf_iyH(%esp),%mm3
        pfsubr nb301nf_izH(%esp),%mm4   ## mm2-mm4 is dxH-dzH 

        pfmul %mm2,%mm2
        pfmul %mm3,%mm3
        pfmul %mm4,%mm4

        pfadd %mm2,%mm3
        pfadd %mm4,%mm3         ## mm3=rsqH 
        movq %mm3,nb301nf_tmprsqH(%esp)

        pfrsqrt %mm0,%mm1

        movq %mm1,%mm2
        pfmul %mm1,%mm1
        pfrsqit1 %mm0,%mm1
        pfrcpit2 %mm2,%mm1      ## mm1=invsqrt 

        pfmul %mm1,%mm0         ## mm0=r 

        pfmul nb301nf_tsc(%esp),%mm0
        pf2iw %mm0,%mm4
        movd %mm4,nb301nf_n1(%esp)
        pi2fd %mm4,%mm4
        pfsub %mm4,%mm0              ## now mm0 is eps and mm4 n0 
        movq  %mm0,%mm2
        pfmul %mm2,%mm2         ## mm0 is eps, mm2 eps2 

        ## coulomb table 
        movl nb301nf_VFtab(%ebp),%edx
        movl nb301nf_n1(%esp),%ecx
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

        pfmul nb301nf_qqO(%esp),%mm5    ## vcoul=qq*VV 
        ## update vctot directly 
        pfadd nb301nf_vctot(%esp),%mm5
        movq %mm5,nb301nf_vctot(%esp)

        ## now do the two hydrogens. 
        movq nb301nf_tmprsqH(%esp),%mm0   ## mm0=rsqH 

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
        pfmul nb301nf_tsc(%esp),%mm0
        pf2iw %mm0,%mm4
        movq %mm4,nb301nf_n1(%esp)
        pi2fd %mm4,%mm4
        pfsub %mm4,%mm0              ## now mm0 is eps and mm4 n0 
        movq  %mm0,%mm2
        pfmul %mm2,%mm2         ## mm0 is eps, mm2 eps2 

        ## coulomb table 
        movl nb301nf_VFtab(%ebp),%edx
        movl nb301nf_n1(%esp),%ecx
        shll $2,%ecx
        ## load all values we need 
        movd (%edx,%ecx,4),%mm4
        movd 4(%edx,%ecx,4),%mm5
        movd 8(%edx,%ecx,4),%mm6
        movd 12(%edx,%ecx,4),%mm7
        movl nb301nf_n1+4(%esp),%ecx
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

        pfmul nb301nf_qqH(%esp),%mm5    ## vcoul=qq*VV 

        ## update vctot 
        pfadd nb301nf_vctot(%esp),%mm5
        movq %mm5,nb301nf_vctot(%esp)

        ##  done  - one more? 
        decl nb301nf_innerk(%esp)
        jz  _nb_kernel301nf_ia32_3dnow.nb301nf_updateouterdata
        jmp _nb_kernel301nf_ia32_3dnow.nb301nf_inner_loop
_nb_kernel301nf_ia32_3dnow.nb301nf_updateouterdata: 
        ## get n from stack
        movl nb301nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb301nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        movq  nb301nf_vctot(%esp),%mm7
        pfacc %mm7,%mm7           ## get and sum the two parts of total potential 

        movl  nb301nf_Vc(%ebp),%eax
        movd  (%eax,%edx,4),%mm6
        pfadd %mm7,%mm6
        movd  %mm6,(%eax,%edx,4)          ## increment vc[gid] 

        ## finish if last 
        movl nb301nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jecxz _nb_kernel301nf_ia32_3dnow.nb301nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb301nf_n(%esp)
        jmp _nb_kernel301nf_ia32_3dnow.nb301nf_outer
_nb_kernel301nf_ia32_3dnow.nb301nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb301nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jecxz _nb_kernel301nf_ia32_3dnow.nb301nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel301nf_ia32_3dnow.nb301nf_threadloop
_nb_kernel301nf_ia32_3dnow.nb301nf_end: 
        femms

        movl nb301nf_nouter(%esp),%eax
        movl nb301nf_ninner(%esp),%ebx
        movl nb301nf_outeriter(%ebp),%ecx
        movl nb301nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        addl $136,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret



