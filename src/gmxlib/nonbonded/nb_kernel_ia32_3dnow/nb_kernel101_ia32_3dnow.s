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


.globl nb_kernel101_ia32_3dnow
.globl _nb_kernel101_ia32_3dnow
nb_kernel101_ia32_3dnow:        
_nb_kernel101_ia32_3dnow:       
.set nb101_p_nri, 8
.set nb101_iinr, 12
.set nb101_jindex, 16
.set nb101_jjnr, 20
.set nb101_shift, 24
.set nb101_shiftvec, 28
.set nb101_fshift, 32
.set nb101_gid, 36
.set nb101_pos, 40
.set nb101_faction, 44
.set nb101_charge, 48
.set nb101_p_facel, 52
.set nb101_p_krf, 56
.set nb101_p_crf, 60
.set nb101_Vc, 64
.set nb101_type, 68
.set nb101_p_ntype, 72
.set nb101_vdwparam, 76
.set nb101_Vvdw, 80
.set nb101_p_tabscale, 84
.set nb101_VFtab, 88
.set nb101_invsqrta, 92
.set nb101_dvda, 96
.set nb101_p_gbtabscale, 100
.set nb101_GBtab, 104
.set nb101_p_nthreads, 108
.set nb101_count, 112
.set nb101_mtx, 116
.set nb101_outeriter, 120
.set nb101_inneriter, 124
.set nb101_work, 128
                        ## stack offsets for local variables 
.set nb101_is3, 0
.set nb101_ii3, 4
.set nb101_ixO, 8
.set nb101_iyO, 12
.set nb101_izO, 16
.set nb101_ixH, 20
.set nb101_iyH, 28
.set nb101_izH, 36
.set nb101_iqO, 44
.set nb101_iqH, 52
.set nb101_vctot, 60
.set nb101_innerjjnr, 68
.set nb101_innerk, 72
.set nb101_fixO, 76
.set nb101_fiyO, 80
.set nb101_fizO, 84
.set nb101_fixH, 88
.set nb101_fiyH, 96
.set nb101_fizH, 104
.set nb101_dxO, 112
.set nb101_dyO, 116
.set nb101_dzO, 120
.set nb101_dxH, 124
.set nb101_dyH, 132
.set nb101_dzH, 140
.set nb101_n, 148                           ## idx for outer loop
.set nb101_nn1, 152                         ## number of outer iterations
.set nb101_nri, 156
.set nb101_nouter, 160
.set nb101_ninner, 164
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $168,%esp          ## local stack space 
        femms

        movl nb101_p_nri(%ebp),%ecx
        movl nb101_p_facel(%ebp),%esi
        movl (%ecx),%ecx
        movl %ecx,nb101_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb101_nouter(%esp)
        movl %eax,nb101_ninner(%esp)

        ## assume we have at least one i particle - start directly      

        movl  nb101_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx),%ebx           ## ebx=ii 

        movl  nb101_charge(%ebp),%edx
        movd  (%esi),%mm1
        movd  (%edx,%ebx,4),%mm2    ## mm2=charge[ii0] 
        pfmul %mm1,%mm2
        movq  %mm2,nb101_iqO(%esp)          ## iqO = facel*charge[ii] 

        movd  4(%edx,%ebx,4),%mm2       ## mm2=charge[ii0+1] 
        pfmul %mm1,%mm2
        punpckldq %mm2,%mm2         ## spread to both halves 
        movq  %mm2,nb101_iqH(%esp)          ## iqH = facel*charge[ii0+1] 
_nb_kernel101_ia32_3dnow.nb101_threadloop: 
        movl  nb101_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel101_ia32_3dnow.nb101_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel101_ia32_3dnow.nb101_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb101_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb101_n(%esp)
        movl %ebx,nb101_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel101_ia32_3dnow.nb101_outerstart
        jmp _nb_kernel101_ia32_3dnow.nb101_end

_nb_kernel101_ia32_3dnow.nb101_outerstart: 
        ## ebx contains number of outer iterations
        addl nb101_nouter(%esp),%ebx
        movl %ebx,nb101_nouter(%esp)

_nb_kernel101_ia32_3dnow.nb101_outer: 
        movl  nb101_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb101_is3(%esp)      ## store is3 

        movl  nb101_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movq  (%eax,%ebx,4),%mm5        ## move shX/shY to mm5 and shZ to mm6 
        movd  8(%eax,%ebx,4),%mm6
        movq  %mm5,%mm0
        movq  %mm5,%mm1
        movq  %mm6,%mm2
        punpckldq %mm0,%mm0         ## also expand shX,Y,Z in mm0--mm2 
        punpckhdq %mm1,%mm1
        punpckldq %mm2,%mm2

        movl  nb101_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx,%esi,4),%ebx            ## ebx=ii 

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb101_pos(%ebp),%eax      ## eax = base of pos[] 

        pfadd (%eax,%ebx,4),%mm5    ## ix = shX + posX (and iy too) 
        movd  8(%eax,%ebx,4),%mm7       ## cant use direct memory add for 4 bytes (iz) 
        movl  %ebx,nb101_ii3(%esp)          ## (use mm7 as temp storage for iz) 
        pfadd %mm7,%mm6
        movq  %mm5,nb101_ixO(%esp)
        movq  %mm6,nb101_izO(%esp)

        movd  12(%eax,%ebx,4),%mm3
        movd  16(%eax,%ebx,4),%mm4
        movd  20(%eax,%ebx,4),%mm5
        punpckldq  24(%eax,%ebx,4),%mm3
        punpckldq  28(%eax,%ebx,4),%mm4
        punpckldq  32(%eax,%ebx,4),%mm5    ## coords of H1 in low mm3-mm5, H2 in high 

        pfadd %mm3,%mm0
        pfadd %mm4,%mm1
        pfadd %mm5,%mm2
        movq %mm0,nb101_ixH(%esp)
        movq %mm1,nb101_iyH(%esp)
        movq %mm2,nb101_izH(%esp)

        ## clear vctot and i forces 
        pxor  %mm7,%mm7
        movq  %mm7,nb101_vctot(%esp)
        movq  %mm7,nb101_fixO(%esp)
        movd  %mm7,nb101_fizO(%esp)
        movq  %mm7,nb101_fixH(%esp)
        movq  %mm7,nb101_fiyH(%esp)
        movq  %mm7,nb101_fizH(%esp)

        movl  nb101_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 
        movl  %edx,nb101_innerk(%esp)      ## number of innerloop atoms 
        addl  nb101_ninner(%esp),%edx
        movl  %edx,nb101_ninner(%esp)

        movl  nb101_pos(%ebp),%esi
        movl  nb101_faction(%ebp),%edi
        movl  nb101_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb101_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
_nb_kernel101_ia32_3dnow.nb101_inner_loop: 
        ## a single j particle iteration here - compare with the unrolled code for comments 
        movl  nb101_innerjjnr(%esp),%eax
        movl  (%eax),%eax       ## eax=jnr offset 
        addl $4,nb101_innerjjnr(%esp)             ## advance pointer 
        prefetch 16(%ecx)          ## prefetch data - trial and error says 16 is best 

        movl nb101_charge(%ebp),%ecx
        movd (%ecx,%eax,4),%mm7
        punpckldq %mm7,%mm7
        movq %mm7,%mm6
        pfmul nb101_iqO(%esp),%mm6
        pfmul nb101_iqH(%esp),%mm7      ## mm6=qqO, mm7=qqH 

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

        pfsubr nb101_ixO(%esp),%mm0
        pfsubr nb101_izO(%esp),%mm1

        movq  %mm0,nb101_dxO(%esp)
        pfmul %mm0,%mm0
        movd  %mm1,nb101_dzO(%esp)
        pfmul %mm1,%mm1
        pfacc %mm1,%mm0
        pfadd %mm1,%mm0         ## mm0=rsqO 

        punpckldq %mm2,%mm2
        punpckldq %mm3,%mm3
        punpckldq %mm4,%mm4 ## mm2-mm4 is jx-jz 
        pfsubr nb101_ixH(%esp),%mm2
        pfsubr nb101_iyH(%esp),%mm3
        pfsubr nb101_izH(%esp),%mm4   ## mm2-mm4 is dxH-dzH 

        movq %mm2,nb101_dxH(%esp)
        movq %mm3,nb101_dyH(%esp)
        movq %mm4,nb101_dzH(%esp)
        pfmul %mm2,%mm2
        pfmul %mm3,%mm3
        pfmul %mm4,%mm4

        pfadd %mm2,%mm3
        pfadd %mm4,%mm3         ## mm3=rsqH 

        pfrsqrt %mm0,%mm1

        movq %mm1,%mm2
        pfmul %mm1,%mm1
        pfrsqit1 %mm0,%mm1
        pfrcpit2 %mm2,%mm1      ## mm1=invsqrt 
        movq  %mm1,%mm4
        pfmul %mm4,%mm4         ## mm4=invsq 
        ## calculate potential and scalar force 
        pfmul %mm1,%mm6         ## mm6=vcoul 
        pfmul %mm6,%mm4         ## mm4=fscalar  

        pfrsqrt %mm3,%mm5
        pswapd %mm3,%mm3
        pfrsqrt %mm3,%mm2
        pswapd %mm3,%mm3
        punpckldq %mm2,%mm5     ## seeds are in mm5 now, and rsq in mm3 

        movq %mm5,%mm2
        pfmul %mm5,%mm5
        pfrsqit1 %mm3,%mm5
        pfrcpit2 %mm2,%mm5      ## mm5=invsqrt 
        movq %mm5,%mm3
        pfmul %mm3,%mm3         ## mm3=invsq 
        pfmul %mm5,%mm7         ## mm7=vcoul 
        pfmul %mm7,%mm3         ## mm3=fscal for the two H's 

        ## update vctot 
        pfadd %mm6,%mm7
        pfadd nb101_vctot(%esp),%mm7
        movq %mm7,nb101_vctot(%esp)

        ## spread oxygen fscalar to both positions 
        punpckldq %mm4,%mm4
        ## calc vectorial force for O 
        prefetchw (%edi,%eax,4) ## prefetch faction to cache  
        movq nb101_dxO(%esp),%mm0
        movd nb101_dzO(%esp),%mm1
        pfmul %mm4,%mm0
        pfmul %mm4,%mm1

        ## calc vectorial force for H's 
        movq nb101_dxH(%esp),%mm5
        movq nb101_dyH(%esp),%mm6
        movq nb101_dzH(%esp),%mm7
        pfmul %mm3,%mm5
        pfmul %mm3,%mm6
        pfmul %mm3,%mm7

        ## update iO particle force 
        movq nb101_fixO(%esp),%mm2
        movd nb101_fizO(%esp),%mm3
        pfadd %mm0,%mm2
        pfadd %mm1,%mm3
        movq %mm2,nb101_fixO(%esp)
        movd %mm3,nb101_fizO(%esp)

        ## update iH forces 
        movq nb101_fixH(%esp),%mm2
        movq nb101_fiyH(%esp),%mm3
        movq nb101_fizH(%esp),%mm4
        pfadd %mm5,%mm2
        pfadd %mm6,%mm3
        pfadd %mm7,%mm4
        movq %mm2,nb101_fixH(%esp)
        movq %mm3,nb101_fiyH(%esp)
        movq %mm4,nb101_fizH(%esp)

        ## pack j forces from H in the same form as the oxygen force 
        pfacc %mm6,%mm5         ## mm5(l)=fjx(H1+ h2) mm5(h)=fjy(H1+ h2) 
        pfacc %mm7,%mm7         ## mm7(l)=fjz(H1+ h2) 

        pfadd %mm5,%mm0         ## add up total force on j particle  
        pfadd %mm7,%mm1

        ## update j particle force 
        movq (%edi,%eax,4),%mm2
        movd 8(%edi,%eax,4),%mm3
        pfsub %mm0,%mm2
        pfsub %mm1,%mm3
        movq %mm2,(%edi,%eax,4)
        movd %mm3,8(%edi,%eax,4)

        ##  done  - one more? 
        decl nb101_innerk(%esp)
        jz  _nb_kernel101_ia32_3dnow.nb101_updateouterdata
        jmp _nb_kernel101_ia32_3dnow.nb101_inner_loop
_nb_kernel101_ia32_3dnow.nb101_updateouterdata: 
        movl  nb101_ii3(%esp),%ecx

        movq  (%edi,%ecx,4),%mm6       ## increment iO force  
        movd  8(%edi,%ecx,4),%mm7
        pfadd nb101_fixO(%esp),%mm6
        pfadd nb101_fizO(%esp),%mm7
        movq  %mm6,(%edi,%ecx,4)
        movd  %mm7,8(%edi,%ecx,4)

        movq  nb101_fixH(%esp),%mm0
        movq  nb101_fiyH(%esp),%mm3
        movq  nb101_fizH(%esp),%mm1
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


        movl  nb101_fshift(%ebp),%ebx      ## increment fshift force 
        movl  nb101_is3(%esp),%edx

        movq  (%ebx,%edx,4),%mm6
        movd  8(%ebx,%edx,4),%mm7
        pfadd nb101_fixO(%esp),%mm6
        pfadd nb101_fizO(%esp),%mm7
        pfadd %mm0,%mm6
        pfadd %mm1,%mm7
        pfadd %mm2,%mm6
        pfadd %mm3,%mm7
        movq  %mm6,(%ebx,%edx,4)
        movd  %mm7,8(%ebx,%edx,4)

        ## get n from stack
        movl nb101_n(%esp),%esi
        ## get group index for i particle 
        movl  nb101_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        movq  nb101_vctot(%esp),%mm7
        pfacc %mm7,%mm7           ## get and sum the two parts of total potential 

        movl  nb101_Vc(%ebp),%eax
        movd  (%eax,%edx,4),%mm6
        pfadd %mm7,%mm6
        movd  %mm6,(%eax,%edx,4)          ## increment vc[gid] 

       ## finish if last 
        movl nb101_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jecxz _nb_kernel101_ia32_3dnow.nb101_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb101_n(%esp)
        jmp _nb_kernel101_ia32_3dnow.nb101_outer
_nb_kernel101_ia32_3dnow.nb101_outerend: 
        ## check if more outer neighborlists remain
        movl  nb101_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jecxz _nb_kernel101_ia32_3dnow.nb101_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel101_ia32_3dnow.nb101_threadloop
_nb_kernel101_ia32_3dnow.nb101_end: 
        femms
        movl nb101_nouter(%esp),%eax
        movl nb101_ninner(%esp),%ebx
        movl nb101_outeriter(%ebp),%ecx
        movl nb101_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        addl $168,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret





.globl nb_kernel101nf_ia32_3dnow
.globl _nb_kernel101nf_ia32_3dnow
nb_kernel101nf_ia32_3dnow:      
_nb_kernel101nf_ia32_3dnow:     
.set nb101nf_p_nri, 8
.set nb101nf_iinr, 12
.set nb101nf_jindex, 16
.set nb101nf_jjnr, 20
.set nb101nf_shift, 24
.set nb101nf_shiftvec, 28
.set nb101nf_fshift, 32
.set nb101nf_gid, 36
.set nb101nf_pos, 40
.set nb101nf_faction, 44
.set nb101nf_charge, 48
.set nb101nf_p_facel, 52
.set nb101nf_p_krf, 56
.set nb101nf_p_crf, 60
.set nb101nf_Vc, 64
.set nb101nf_type, 68
.set nb101nf_p_ntype, 72
.set nb101nf_vdwparam, 76
.set nb101nf_Vvdw, 80
.set nb101nf_p_tabscale, 84
.set nb101nf_VFtab, 88
.set nb101nf_invsqrta, 92
.set nb101nf_dvda, 96
.set nb101nf_p_gbtabscale, 100
.set nb101nf_GBtab, 104
.set nb101nf_p_nthreads, 108
.set nb101nf_count, 112
.set nb101nf_mtx, 116
.set nb101nf_outeriter, 120
.set nb101nf_inneriter, 124
.set nb101nf_work, 128
                        ## stack offsets for local variables 
.set nb101nf_is3, 0
.set nb101nf_ii3, 4
.set nb101nf_ixO, 8
.set nb101nf_iyO, 12
.set nb101nf_izO, 16
.set nb101nf_ixH, 20
.set nb101nf_iyH, 28
.set nb101nf_izH, 36
.set nb101nf_iqO, 44
.set nb101nf_iqH, 52
.set nb101nf_vctot, 60
.set nb101nf_innerjjnr, 68
.set nb101nf_innerk, 72
.set nb101nf_n, 76                         ## idx for outer loop
.set nb101nf_nn1, 80                       ## number of outer iterations
.set nb101nf_nri, 84
.set nb101nf_nouter, 88
.set nb101nf_ninner, 92
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $96,%esp           ## local stack space 
        femms

        movl nb101nf_p_nri(%ebp),%ecx
        movl nb101nf_p_facel(%ebp),%esi
        movl (%ecx),%ecx
        movl %ecx,nb101nf_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb101nf_nouter(%esp)
        movl %eax,nb101nf_ninner(%esp)

        ## assume we have at least one i particle - start directly      

        movl  nb101nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx),%ebx           ## ebx=ii 

        movl  nb101nf_charge(%ebp),%edx
        movd  (%esi),%mm1
        movd  (%edx,%ebx,4),%mm2    ## mm2=charge[ii0] 
        pfmul %mm1,%mm2
        movq  %mm2,nb101nf_iqO(%esp)        ## iqO = facel*charge[ii] 

        movd  4(%edx,%ebx,4),%mm2       ## mm2=charge[ii0+1] 
        pfmul %mm1,%mm2
        punpckldq %mm2,%mm2         ## spread to both halves 
        movq  %mm2,nb101nf_iqH(%esp)        ## iqH = facel*charge[ii0+1] 
_nb_kernel101nf_ia32_3dnow.nb101nf_threadloop: 
        movl  nb101nf_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel101nf_ia32_3dnow.nb101nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel101nf_ia32_3dnow.nb101nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb101nf_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb101nf_n(%esp)
        movl %ebx,nb101nf_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel101nf_ia32_3dnow.nb101nf_outerstart
        jmp _nb_kernel101nf_ia32_3dnow.nb101nf_end

_nb_kernel101nf_ia32_3dnow.nb101nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb101nf_nouter(%esp),%ebx
        movl %ebx,nb101nf_nouter(%esp)

_nb_kernel101nf_ia32_3dnow.nb101nf_outer: 
        movl  nb101nf_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb101nf_is3(%esp)            ## store is3 

        movl  nb101nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movq  (%eax,%ebx,4),%mm5        ## move shX/shY to mm5 and shZ to mm6 
        movd  8(%eax,%ebx,4),%mm6
        movq  %mm5,%mm0
        movq  %mm5,%mm1
        movq  %mm6,%mm2
        punpckldq %mm0,%mm0         ## also expand shX,Y,Z in mm0--mm2 
        punpckhdq %mm1,%mm1
        punpckldq %mm2,%mm2

        movl  nb101nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx,%esi,4),%ebx            ## ebx=ii 

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb101nf_pos(%ebp),%eax      ## eax = base of pos[] 

        pfadd (%eax,%ebx,4),%mm5    ## ix = shX + posX (and iy too) 
        movd  8(%eax,%ebx,4),%mm7       ## cant use direct memory add for 4 bytes (iz) 
        movl  %ebx,nb101nf_ii3(%esp)        ## (use mm7 as temp storage for iz) 
        pfadd %mm7,%mm6
        movq  %mm5,nb101nf_ixO(%esp)
        movq  %mm6,nb101nf_izO(%esp)

        movd  12(%eax,%ebx,4),%mm3
        movd  16(%eax,%ebx,4),%mm4
        movd  20(%eax,%ebx,4),%mm5
        punpckldq  24(%eax,%ebx,4),%mm3
        punpckldq  28(%eax,%ebx,4),%mm4
        punpckldq  32(%eax,%ebx,4),%mm5
        ## coords of H1 in low mm3-mm5, H2 in high 

        pfadd %mm3,%mm0
        pfadd %mm4,%mm1
        pfadd %mm5,%mm2
        movq %mm0,nb101nf_ixH(%esp)
        movq %mm1,nb101nf_iyH(%esp)
        movq %mm2,nb101nf_izH(%esp)

        ## clear vctot 
        pxor  %mm7,%mm7
        movq  %mm7,nb101nf_vctot(%esp)

        movl  nb101nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 
        movl  %edx,nb101nf_innerk(%esp)      ## number of innerloop atoms 
        addl  nb101nf_ninner(%esp),%edx
        movl  %edx,nb101nf_ninner(%esp)

        movl  nb101nf_pos(%ebp),%esi
        movl  nb101nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb101nf_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
_nb_kernel101nf_ia32_3dnow.nb101nf_inner_loop: 
        ## a single j particle iteration here 
        movl  nb101nf_innerjjnr(%esp),%eax
        movl  (%eax),%eax       ## eax=jnr offset 
        addl $4,nb101nf_innerjjnr(%esp)             ## advance pointer 
        ## prefetch data - trial and error says 16 is best 
        prefetch 16(%ecx)

        movl nb101nf_charge(%ebp),%ecx
        movd (%ecx,%eax,4),%mm7
        punpckldq %mm7,%mm7
        movq %mm7,%mm6
        pfmul nb101nf_iqO(%esp),%mm6
        pfmul nb101nf_iqH(%esp),%mm7    ## mm6=qqO, mm7=qqH 

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

        pfsubr nb101nf_ixO(%esp),%mm0
        pfsubr nb101nf_izO(%esp),%mm1

        pfmul %mm0,%mm0
        pfmul %mm1,%mm1
        pfacc %mm1,%mm0
        pfadd %mm1,%mm0         ## mm0=rsqO 

        punpckldq %mm2,%mm2
        punpckldq %mm3,%mm3
        punpckldq %mm4,%mm4 ## mm2-mm4 is jx-jz 
        pfsubr nb101nf_ixH(%esp),%mm2
        pfsubr nb101nf_iyH(%esp),%mm3
        pfsubr nb101nf_izH(%esp),%mm4   ## mm2-mm4 is dxH-dzH 

        pfmul %mm2,%mm2
        pfmul %mm3,%mm3
        pfmul %mm4,%mm4

        pfadd %mm2,%mm3
        pfadd %mm4,%mm3         ## mm3=rsqH 

        pfrsqrt %mm0,%mm1

        movq %mm1,%mm2
        pfmul %mm1,%mm1
        pfrsqit1 %mm0,%mm1
        pfrcpit2 %mm2,%mm1      ## mm1=invsqrt 
        ## calculate potential and scalar force 
        pfmul %mm1,%mm6         ## mm6=vcoul 

        pfrsqrt %mm3,%mm5
        pswapd %mm3,%mm3
        pfrsqrt %mm3,%mm2
        pswapd %mm3,%mm3
        punpckldq %mm2,%mm5     ## seeds are in mm5 now, and rsq in mm3 

        movq %mm5,%mm2
        pfmul %mm5,%mm5
        pfrsqit1 %mm3,%mm5
        pfrcpit2 %mm2,%mm5      ## mm5=invsqrt 
        pfmul %mm5,%mm7         ## mm7=vcoul 

        ## update vctot 
        pfadd %mm6,%mm7
        pfadd nb101nf_vctot(%esp),%mm7
        movq %mm7,nb101nf_vctot(%esp)

        ##  done  - one more? 
        decl nb101nf_innerk(%esp)
        jz  _nb_kernel101nf_ia32_3dnow.nb101nf_updateouterdata
        jmp _nb_kernel101nf_ia32_3dnow.nb101nf_inner_loop
_nb_kernel101nf_ia32_3dnow.nb101nf_updateouterdata: 
        ## get n from stack
        movl nb101nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb101nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        movq  nb101nf_vctot(%esp),%mm7
        pfacc %mm7,%mm7           ## get and sum the two parts of total potential 

        movl  nb101nf_Vc(%ebp),%eax
        movd  (%eax,%edx,4),%mm6
        pfadd %mm7,%mm6
        movd  %mm6,(%eax,%edx,4)          ## increment vc[gid] 

       ## finish if last 
        movl nb101nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jecxz _nb_kernel101nf_ia32_3dnow.nb101nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb101nf_n(%esp)
        jmp _nb_kernel101nf_ia32_3dnow.nb101nf_outer
_nb_kernel101nf_ia32_3dnow.nb101nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb101nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jecxz _nb_kernel101nf_ia32_3dnow.nb101nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel101nf_ia32_3dnow.nb101nf_threadloop
_nb_kernel101nf_ia32_3dnow.nb101nf_end: 
        femms
        movl nb101nf_nouter(%esp),%eax
        movl nb101nf_ninner(%esp),%ebx
        movl nb101nf_outeriter(%ebp),%ecx
        movl nb101nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        addl $96,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret




