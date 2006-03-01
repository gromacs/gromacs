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




.globl nb_kernel102_ia32_3dnow
.globl _nb_kernel102_ia32_3dnow
nb_kernel102_ia32_3dnow:        
_nb_kernel102_ia32_3dnow:       
.set nb102_p_nri, 8
.set nb102_iinr, 12
.set nb102_jindex, 16
.set nb102_jjnr, 20
.set nb102_shift, 24
.set nb102_shiftvec, 28
.set nb102_fshift, 32
.set nb102_gid, 36
.set nb102_pos, 40
.set nb102_faction, 44
.set nb102_charge, 48
.set nb102_p_facel, 52
.set nb102_p_krf, 56
.set nb102_p_crf, 60
.set nb102_Vc, 64
.set nb102_type, 68
.set nb102_p_ntype, 72
.set nb102_vdwparam, 76
.set nb102_Vvdw, 80
.set nb102_p_tabscale, 84
.set nb102_VFtab, 88
.set nb102_invsqrta, 92
.set nb102_dvda, 96
.set nb102_p_gbtabscale, 100
.set nb102_GBtab, 104
.set nb102_p_nthreads, 108
.set nb102_count, 112
.set nb102_mtx, 116
.set nb102_outeriter, 120
.set nb102_inneriter, 124
.set nb102_work, 128
                        ## stack offsets for local variables 
.set nb102_is3, 0
.set nb102_ii3, 4
.set nb102_ixO, 8
.set nb102_iyO, 12
.set nb102_izO, 16
.set nb102_ixH, 20
.set nb102_iyH, 28
.set nb102_izH, 36
.set nb102_qqOO, 44
.set nb102_qqOH, 52
.set nb102_qqHH, 60
.set nb102_vctot, 68
.set nb102_innerjjnr, 76
.set nb102_innerk, 80
.set nb102_fixO, 84
.set nb102_fiyO, 88
.set nb102_fizO, 92
.set nb102_fixH, 96
.set nb102_fiyH, 104
.set nb102_fizH, 112
.set nb102_dxO, 120
.set nb102_dyO, 124
.set nb102_dzO, 128
.set nb102_dxH, 132
.set nb102_dyH, 140
.set nb102_dzH, 148
.set nb102_n, 156                           ## idx for outer loop
.set nb102_nn1, 160                         ## number of outer iterations
.set nb102_nri, 164
.set nb102_nouter, 168
.set nb102_ninner, 172
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $176,%esp          ## local stack space 
        femms
        movl nb102_p_nri(%ebp),%ecx
        movl nb102_p_facel(%ebp),%esi
        movl (%ecx),%ecx
        movl %ecx,nb102_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb102_nouter(%esp)
        movl %eax,nb102_ninner(%esp)

        ## assume we have at least one i particle - start directly      

        movl  nb102_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx),%ebx           ## ebx=ii 

        movl  nb102_charge(%ebp),%edx
        movd  (%esi),%mm1              ## mm1=facel 
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
        movq  %mm4,nb102_qqOO(%esp)
        movq  %mm5,nb102_qqOH(%esp)
        movq  %mm6,nb102_qqHH(%esp)
_nb_kernel102_ia32_3dnow.nb102_threadloop: 
        movl  nb102_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel102_ia32_3dnow.nb102_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel102_ia32_3dnow.nb102_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb102_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb102_n(%esp)
        movl %ebx,nb102_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel102_ia32_3dnow.nb102_outerstart
        jmp _nb_kernel102_ia32_3dnow.nb102_end

_nb_kernel102_ia32_3dnow.nb102_outerstart: 
        ## ebx contains number of outer iterations
        addl nb102_nouter(%esp),%ebx
        movl %ebx,nb102_nouter(%esp)

_nb_kernel102_ia32_3dnow.nb102_outer: 
        movl  nb102_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb102_is3(%esp)      ## store is3 

        movl  nb102_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movq  (%eax,%ebx,4),%mm5        ## move shX/shY to mm5 and shZ to mm6 
        movd  8(%eax,%ebx,4),%mm6
        movq  %mm5,%mm0
        movq  %mm5,%mm1
        movq  %mm6,%mm2
        punpckldq %mm0,%mm0         ## also expand shX,Y,Z in mm0--mm2 
        punpckhdq %mm1,%mm1
        punpckldq %mm2,%mm2

        movl  nb102_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx,%esi,4),%ebx            ## ebx=ii 

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb102_pos(%ebp),%eax      ## eax = base of pos[] 

        pfadd (%eax,%ebx,4),%mm5    ## ix = shX + posX (and iy too) 
        movd  8(%eax,%ebx,4),%mm7       ## cant use direct memory add for 4 bytes (iz) 
        movl  %ebx,nb102_ii3(%esp)          ## (use mm7 as temp storage for iz) 
        pfadd %mm7,%mm6
        movq  %mm5,nb102_ixO(%esp)
        movq  %mm6,nb102_izO(%esp)

        movd  12(%eax,%ebx,4),%mm3
        movd  16(%eax,%ebx,4),%mm4
        movd  20(%eax,%ebx,4),%mm5
        punpckldq  24(%eax,%ebx,4),%mm3
        punpckldq  28(%eax,%ebx,4),%mm4
        punpckldq  32(%eax,%ebx,4),%mm5    ## coords of H1 in low mm3-mm5, H2 in high 

        pfadd %mm3,%mm0
        pfadd %mm4,%mm1
        pfadd %mm5,%mm2
        movq %mm0,nb102_ixH(%esp)
        movq %mm1,nb102_iyH(%esp)
        movq %mm2,nb102_izH(%esp)

        ## clear vctot and i forces 
        pxor  %mm7,%mm7
        movq  %mm7,nb102_vctot(%esp)
        movq  %mm7,nb102_fixO(%esp)
        movq  %mm7,nb102_fizO(%esp)
        movq  %mm7,nb102_fixH(%esp)
        movq  %mm7,nb102_fiyH(%esp)
        movq  %mm7,nb102_fizH(%esp)

        movl  nb102_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 
        movl  %edx,nb102_innerk(%esp)      ## number of innerloop atoms 
        addl  nb102_ninner(%esp),%edx
        movl  %edx,nb102_ninner(%esp)

        movl  nb102_pos(%ebp),%esi
        movl  nb102_faction(%ebp),%edi
        movl  nb102_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb102_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
_nb_kernel102_ia32_3dnow.nb102_inner_loop: 
        ## a single j particle iteration here - compare with the unrolled code for comments 
        movl  nb102_innerjjnr(%esp),%eax
        movl  (%eax),%eax       ## eax=jnr offset 
    addl $4,nb102_innerjjnr(%esp)             ## advance pointer 

        movd  nb102_qqOO(%esp),%mm6
        movq  nb102_qqOH(%esp),%mm7

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

        pfsubr nb102_ixO(%esp),%mm0
        pfsubr nb102_izO(%esp),%mm1

        movq  %mm0,nb102_dxO(%esp)
        pfmul %mm0,%mm0
        movd  %mm1,nb102_dzO(%esp)
        pfmul %mm1,%mm1
        pfacc %mm0,%mm0
        pfadd %mm1,%mm0         ## mm0=rsqO 

        punpckldq %mm2,%mm2
        punpckldq %mm3,%mm3
        punpckldq %mm4,%mm4 ## mm2-mm4 is jx-jz 
        pfsubr nb102_ixH(%esp),%mm2
        pfsubr nb102_iyH(%esp),%mm3
        pfsubr nb102_izH(%esp),%mm4   ## mm2-mm4 is dxH-dzH 

        movq %mm2,nb102_dxH(%esp)
        movq %mm3,nb102_dyH(%esp)
        movq %mm4,nb102_dzH(%esp)
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
        pfadd nb102_vctot(%esp),%mm7
        movq %mm7,nb102_vctot(%esp)

        ## spread oxygen fscalar to both positions 
        punpckldq %mm4,%mm4
        ## calc vectorial force for O 
        movq nb102_dxO(%esp),%mm0
        movd nb102_dzO(%esp),%mm1
        pfmul %mm4,%mm0
        pfmul %mm4,%mm1

        ## calc vectorial force for H's 
        movq nb102_dxH(%esp),%mm5
        movq nb102_dyH(%esp),%mm6
        movq nb102_dzH(%esp),%mm7
        pfmul %mm3,%mm5
        pfmul %mm3,%mm6
        pfmul %mm3,%mm7

        ## update iO particle force 
        movq nb102_fixO(%esp),%mm2
        movd nb102_fizO(%esp),%mm3
        pfadd %mm0,%mm2
        pfadd %mm1,%mm3
        movq %mm2,nb102_fixO(%esp)
        movd %mm3,nb102_fizO(%esp)

        ## update iH forces 
        movq nb102_fixH(%esp),%mm2
        movq nb102_fiyH(%esp),%mm3
        movq nb102_fizH(%esp),%mm4
        pfadd %mm5,%mm2
        pfadd %mm6,%mm3
        pfadd %mm7,%mm4
        movq %mm2,nb102_fixH(%esp)
        movq %mm3,nb102_fiyH(%esp)
        movq %mm4,nb102_fizH(%esp)

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

        movd nb102_qqOH(%esp),%mm6
        movq nb102_qqHH(%esp),%mm7

        pfsubr nb102_ixO(%esp),%mm0
        pfsubr nb102_izO(%esp),%mm1

        movq  %mm0,nb102_dxO(%esp)
        pfmul %mm0,%mm0
        movd  %mm1,nb102_dzO(%esp)
        pfmul %mm1,%mm1
        pfacc %mm1,%mm0
        pfadd %mm1,%mm0         ## mm0=rsqO 

        punpckldq %mm2,%mm2
        punpckldq %mm3,%mm3
        punpckldq %mm4,%mm4 ## mm2-mm4 is jx-jz 
        pfsubr nb102_ixH(%esp),%mm2
        pfsubr nb102_iyH(%esp),%mm3
        pfsubr nb102_izH(%esp),%mm4   ## mm2-mm4 is dxH-dzH 

        movq %mm2,nb102_dxH(%esp)
        movq %mm3,nb102_dyH(%esp)
        movq %mm4,nb102_dzH(%esp)
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
        pfadd nb102_vctot(%esp),%mm7
        movq %mm7,nb102_vctot(%esp)

        ## spread oxygen fscalar to both positions 
        punpckldq %mm4,%mm4
        ## calc vectorial force for O 
        movq nb102_dxO(%esp),%mm0
        movd nb102_dzO(%esp),%mm1
        pfmul %mm4,%mm0
        pfmul %mm4,%mm1

        ## calc vectorial force for H's 
        movq nb102_dxH(%esp),%mm5
        movq nb102_dyH(%esp),%mm6
        movq nb102_dzH(%esp),%mm7
        pfmul %mm3,%mm5
        pfmul %mm3,%mm6
        pfmul %mm3,%mm7

        ## update iO particle force 
        movq nb102_fixO(%esp),%mm2
        movd nb102_fizO(%esp),%mm3
        pfadd %mm0,%mm2
        pfadd %mm1,%mm3
        movq %mm2,nb102_fixO(%esp)
        movd %mm3,nb102_fizO(%esp)

        ## update iH forces 
        movq nb102_fixH(%esp),%mm2
        movq nb102_fiyH(%esp),%mm3
        movq nb102_fizH(%esp),%mm4
        pfadd %mm5,%mm2
        pfadd %mm6,%mm3
        pfadd %mm7,%mm4
        movq %mm2,nb102_fixH(%esp)
        movq %mm3,nb102_fiyH(%esp)
        movq %mm4,nb102_fizH(%esp)

        ## pack j forces from H in the same form as the oxygen force 
        pfacc %mm6,%mm5         ## mm5(l)=fjx(H1+ h2) mm5(h)=fjy(H1+ h2) 
        pfacc %mm7,%mm7         ## mm7(l)=fjz(H1+ h2) 

        pfadd %mm5,%mm0         ## add up total force on j particle  
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

        movd nb102_qqOH(%esp),%mm6
        movq nb102_qqHH(%esp),%mm7

        pfsubr nb102_ixO(%esp),%mm0
        pfsubr nb102_izO(%esp),%mm1

        movq  %mm0,nb102_dxO(%esp)
        pfmul %mm0,%mm0
        movd  %mm1,nb102_dzO(%esp)
        pfmul %mm1,%mm1
        pfacc %mm1,%mm0
        pfadd %mm1,%mm0         ## mm0=rsqO 

        punpckldq %mm2,%mm2
        punpckldq %mm3,%mm3
        punpckldq %mm4,%mm4 ## mm2-mm4 is jx-jz 
        pfsubr nb102_ixH(%esp),%mm2
        pfsubr nb102_iyH(%esp),%mm3
        pfsubr nb102_izH(%esp),%mm4   ## mm2-mm4 is dxH-dzH 

        movq %mm2,nb102_dxH(%esp)
        movq %mm3,nb102_dyH(%esp)
        movq %mm4,nb102_dzH(%esp)
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
        pfadd nb102_vctot(%esp),%mm7
        movq %mm7,nb102_vctot(%esp)

        ## spread oxygen fscalar to both positions 
        punpckldq %mm4,%mm4
        ## calc vectorial force for O 
        movq nb102_dxO(%esp),%mm0
        movd nb102_dzO(%esp),%mm1
        pfmul %mm4,%mm0
        pfmul %mm4,%mm1

        ## calc vectorial force for H's 
        movq nb102_dxH(%esp),%mm5
        movq nb102_dyH(%esp),%mm6
        movq nb102_dzH(%esp),%mm7
        pfmul %mm3,%mm5
        pfmul %mm3,%mm6
        pfmul %mm3,%mm7

        ## update iO particle force 
        movq nb102_fixO(%esp),%mm2
        movd nb102_fizO(%esp),%mm3
        pfadd %mm0,%mm2
        pfadd %mm1,%mm3
        movq %mm2,nb102_fixO(%esp)
        movd %mm3,nb102_fizO(%esp)

        ## update iH forces 
        movq nb102_fixH(%esp),%mm2
        movq nb102_fiyH(%esp),%mm3
        movq nb102_fizH(%esp),%mm4
        pfadd %mm5,%mm2
        pfadd %mm6,%mm3
        pfadd %mm7,%mm4
        movq %mm2,nb102_fixH(%esp)
        movq %mm3,nb102_fiyH(%esp)
        movq %mm4,nb102_fizH(%esp)

        ## pack j forces from H in the same form as the oxygen force 
        pfacc %mm6,%mm5         ## mm5(l)=fjx(H1+ h2) mm5(h)=fjy(H1+ h2) 
        pfacc %mm7,%mm7         ## mm7(l)=fjz(H1+ h2) 

        pfadd %mm5,%mm0         ## add up total force on j particle  
        pfadd %mm7,%mm1

        ## update j particle force 
        movq 24(%edi,%eax,4),%mm2
        movd 32(%edi,%eax,4),%mm3
        pfsub %mm0,%mm2
        pfsub %mm1,%mm3
        movq %mm2,24(%edi,%eax,4)
        movd %mm3,32(%edi,%eax,4)

        ##  done  - one more? 
        decl nb102_innerk(%esp)
        jz  _nb_kernel102_ia32_3dnow.nb102_updateouterdata
        jmp _nb_kernel102_ia32_3dnow.nb102_inner_loop
_nb_kernel102_ia32_3dnow.nb102_updateouterdata: 
        movl  nb102_ii3(%esp),%ecx

        movq  (%edi,%ecx,4),%mm6       ## increment iO force  
        movd  8(%edi,%ecx,4),%mm7
        pfadd nb102_fixO(%esp),%mm6
        pfadd nb102_fizO(%esp),%mm7
        movq  %mm6,(%edi,%ecx,4)
        movd  %mm7,8(%edi,%ecx,4)

        movq  nb102_fixH(%esp),%mm0
        movq  nb102_fiyH(%esp),%mm3
        movq  nb102_fizH(%esp),%mm1
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


        movl  nb102_fshift(%ebp),%ebx      ## increment fshift force 
        movl  nb102_is3(%esp),%edx

        movq  (%ebx,%edx,4),%mm6
        movd  8(%ebx,%edx,4),%mm7
        pfadd nb102_fixO(%esp),%mm6
        pfadd nb102_fizO(%esp),%mm7
        pfadd %mm0,%mm6
        pfadd %mm1,%mm7
        pfadd %mm2,%mm6
        pfadd %mm3,%mm7
        movq  %mm6,(%ebx,%edx,4)
        movd  %mm7,8(%ebx,%edx,4)

        ## get n from stack
        movl nb102_n(%esp),%esi
        ## get group index for i particle 
        movl  nb102_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        movq  nb102_vctot(%esp),%mm7
        pfacc %mm7,%mm7           ## get and sum the two parts of total potential 

        movl  nb102_Vc(%ebp),%eax
        movd  (%eax,%edx,4),%mm6
        pfadd %mm7,%mm6
        movd  %mm6,(%eax,%edx,4)          ## increment vc[gid] 
       ## finish if last 
        movl nb102_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jecxz _nb_kernel102_ia32_3dnow.nb102_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb102_n(%esp)
        jmp _nb_kernel102_ia32_3dnow.nb102_outer
_nb_kernel102_ia32_3dnow.nb102_outerend: 
        ## check if more outer neighborlists remain
        movl  nb102_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jecxz _nb_kernel102_ia32_3dnow.nb102_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel102_ia32_3dnow.nb102_threadloop
_nb_kernel102_ia32_3dnow.nb102_end: 
        femms
        movl nb102_nouter(%esp),%eax
        movl nb102_ninner(%esp),%ebx
        movl nb102_outeriter(%ebp),%ecx
        movl nb102_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)
        addl $176,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret






.globl nb_kernel102nf_ia32_3dnow
.globl _nb_kernel102nf_ia32_3dnow
nb_kernel102nf_ia32_3dnow:      
_nb_kernel102nf_ia32_3dnow:     
.set nb102nf_p_nri, 8
.set nb102nf_iinr, 12
.set nb102nf_jindex, 16
.set nb102nf_jjnr, 20
.set nb102nf_shift, 24
.set nb102nf_shiftvec, 28
.set nb102nf_fshift, 32
.set nb102nf_gid, 36
.set nb102nf_pos, 40
.set nb102nf_faction, 44
.set nb102nf_charge, 48
.set nb102nf_p_facel, 52
.set nb102nf_p_krf, 56
.set nb102nf_p_crf, 60
.set nb102nf_Vc, 64
.set nb102nf_type, 68
.set nb102nf_p_ntype, 72
.set nb102nf_vdwparam, 76
.set nb102nf_Vvdw, 80
.set nb102nf_p_tabscale, 84
.set nb102nf_VFtab, 88
.set nb102nf_invsqrta, 92
.set nb102nf_dvda, 96
.set nb102nf_p_gbtabscale, 100
.set nb102nf_GBtab, 104
.set nb102nf_p_nthreads, 108
.set nb102nf_count, 112
.set nb102nf_mtx, 116
.set nb102nf_outeriter, 120
.set nb102nf_inneriter, 124
.set nb102nf_work, 128
                        ## stack offsets for local variables 
.set nb102nf_is3, 0
.set nb102nf_ii3, 4
.set nb102nf_ixO, 8
.set nb102nf_iyO, 12
.set nb102nf_izO, 16
.set nb102nf_ixH, 20
.set nb102nf_iyH, 28
.set nb102nf_izH, 36
.set nb102nf_qqOO, 44
.set nb102nf_qqOH, 52
.set nb102nf_qqHH, 60
.set nb102nf_vctot, 68
.set nb102nf_innerjjnr, 76
.set nb102nf_innerk, 80
.set nb102nf_n, 84                         ## idx for outer loop
.set nb102nf_nn1, 88                       ## number of outer iterations
.set nb102nf_nri, 92
.set nb102nf_nouter, 96
.set nb102nf_ninner, 100
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $104,%esp          ## local stack space 
        femms
        movl nb102nf_p_nri(%ebp),%ecx
        movl nb102nf_p_facel(%ebp),%esi
        movl (%ecx),%ecx
        movl %ecx,nb102nf_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb102nf_nouter(%esp)
        movl %eax,nb102nf_ninner(%esp)
        ## assume we have at least one i particle - start directly      

        movl  nb102nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx),%ebx           ## ebx=ii 

        movl  nb102nf_charge(%ebp),%edx
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
        movq  %mm4,nb102nf_qqOO(%esp)
        movq  %mm5,nb102nf_qqOH(%esp)
        movq  %mm6,nb102nf_qqHH(%esp)
_nb_kernel102nf_ia32_3dnow.nb102nf_threadloop: 
        movl  nb102nf_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel102nf_ia32_3dnow.nb102nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel102nf_ia32_3dnow.nb102nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb102nf_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb102nf_n(%esp)
        movl %ebx,nb102nf_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel102nf_ia32_3dnow.nb102nf_outerstart
        jmp _nb_kernel102nf_ia32_3dnow.nb102nf_end

_nb_kernel102nf_ia32_3dnow.nb102nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb102nf_nouter(%esp),%ebx
        movl %ebx,nb102nf_nouter(%esp)

_nb_kernel102nf_ia32_3dnow.nb102nf_outer: 
        movl  nb102nf_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb102nf_is3(%esp)            ## store is3 

        movl  nb102nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movq  (%eax,%ebx,4),%mm5        ## move shX/shY to mm5 and shZ to mm6 
        movd  8(%eax,%ebx,4),%mm6
        movq  %mm5,%mm0
        movq  %mm5,%mm1
        movq  %mm6,%mm2
        punpckldq %mm0,%mm0         ## also expand shX,Y,Z in mm0--mm2 
        punpckhdq %mm1,%mm1
        punpckldq %mm2,%mm2

        movl  nb102nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx,%esi,4),%ebx            ## ebx=ii 

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb102nf_pos(%ebp),%eax      ## eax = base of pos[] 

        pfadd (%eax,%ebx,4),%mm5    ## ix = shX + posX (and iy too) 
        movd  8(%eax,%ebx,4),%mm7       ## cant use direct memory add for 4 bytes (iz) 
        movl  %ebx,nb102nf_ii3(%esp)        ## (use mm7 as temp storage for iz) 
        pfadd %mm7,%mm6
        movq  %mm5,nb102nf_ixO(%esp)
        movq  %mm6,nb102nf_izO(%esp)

        movd  12(%eax,%ebx,4),%mm3
        movd  16(%eax,%ebx,4),%mm4
        movd  20(%eax,%ebx,4),%mm5
        punpckldq  24(%eax,%ebx,4),%mm3
        punpckldq  28(%eax,%ebx,4),%mm4
        punpckldq  32(%eax,%ebx,4),%mm5    ## coords of H1 in low mm3-mm5, H2 in high 

        pfadd %mm3,%mm0
        pfadd %mm4,%mm1
        pfadd %mm5,%mm2
        movq %mm0,nb102nf_ixH(%esp)
        movq %mm1,nb102nf_iyH(%esp)
        movq %mm2,nb102nf_izH(%esp)

        ## clear vctot
        pxor  %mm7,%mm7
        movq  %mm7,nb102nf_vctot(%esp)

        movl  nb102nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 
        movl  %edx,nb102nf_innerk(%esp)      ## number of innerloop atoms 
        addl  nb102nf_ninner(%esp),%edx
        movl  %edx,nb102nf_ninner(%esp)

        movl  nb102nf_pos(%ebp),%esi
        movl  nb102nf_faction(%ebp),%edi
        movl  nb102nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb102nf_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
_nb_kernel102nf_ia32_3dnow.nb102nf_inner_loop: 
        ## a single j particle iteration here - compare with the unrolled code for comments 
        movl  nb102nf_innerjjnr(%esp),%eax
        movl  (%eax),%eax       ## eax=jnr offset 
        addl $4,nb102nf_innerjjnr(%esp)             ## advance pointer 

        movd  nb102nf_qqOO(%esp),%mm6
        movq  nb102nf_qqOH(%esp),%mm7

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

        pfsubr nb102nf_ixO(%esp),%mm0
        pfsubr nb102nf_izO(%esp),%mm1

        pfmul %mm0,%mm0
        pfmul %mm1,%mm1
        pfacc %mm0,%mm0
        pfadd %mm1,%mm0         ## mm0=rsqO 

        punpckldq %mm2,%mm2
        punpckldq %mm3,%mm3
        punpckldq %mm4,%mm4 ## mm2-mm4 is jx-jz 
        pfsubr nb102nf_ixH(%esp),%mm2
        pfsubr nb102nf_iyH(%esp),%mm3
        pfsubr nb102nf_izH(%esp),%mm4   ## mm2-mm4 is dxH-dzH 

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
        ## calculate potential 
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
        pfadd nb102nf_vctot(%esp),%mm7
        movq %mm7,nb102nf_vctot(%esp)

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

        movd nb102nf_qqOH(%esp),%mm6
        movq nb102nf_qqHH(%esp),%mm7

        pfsubr nb102nf_ixO(%esp),%mm0
        pfsubr nb102nf_izO(%esp),%mm1

        pfmul %mm0,%mm0
        pfmul %mm1,%mm1
        pfacc %mm1,%mm0
        pfadd %mm1,%mm0         ## mm0=rsqO 

        punpckldq %mm2,%mm2
        punpckldq %mm3,%mm3
        punpckldq %mm4,%mm4 ## mm2-mm4 is jx-jz 
        pfsubr nb102nf_ixH(%esp),%mm2
        pfsubr nb102nf_iyH(%esp),%mm3
        pfsubr nb102nf_izH(%esp),%mm4   ## mm2-mm4 is dxH-dzH 

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
        ## calculate potential 
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
        pfadd nb102nf_vctot(%esp),%mm7
        movq %mm7,nb102nf_vctot(%esp)

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

        movd nb102nf_qqOH(%esp),%mm6
        movq nb102nf_qqHH(%esp),%mm7

        pfsubr nb102nf_ixO(%esp),%mm0
        pfsubr nb102nf_izO(%esp),%mm1

        pfmul %mm0,%mm0
        pfmul %mm1,%mm1
        pfacc %mm1,%mm0
        pfadd %mm1,%mm0         ## mm0=rsqO 

        punpckldq %mm2,%mm2
        punpckldq %mm3,%mm3
        punpckldq %mm4,%mm4 ## mm2-mm4 is jx-jz 
        pfsubr nb102nf_ixH(%esp),%mm2
        pfsubr nb102nf_iyH(%esp),%mm3
        pfsubr nb102nf_izH(%esp),%mm4   ## mm2-mm4 is dxH-dzH 

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
        ## calculate potential 
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
        pfadd nb102nf_vctot(%esp),%mm7
        movq %mm7,nb102nf_vctot(%esp)

        ##  done  - one more? 
        decl nb102nf_innerk(%esp)
        jz  _nb_kernel102nf_ia32_3dnow.nb102nf_updateouterdata
        jmp _nb_kernel102nf_ia32_3dnow.nb102nf_inner_loop
_nb_kernel102nf_ia32_3dnow.nb102nf_updateouterdata: 
        ## get n from stack
        movl nb102nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb102nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        movq  nb102nf_vctot(%esp),%mm7
        pfacc %mm7,%mm7           ## get and sum the two parts of total potential 

        movl  nb102nf_Vc(%ebp),%eax
        movd  (%eax,%edx,4),%mm6
        pfadd %mm7,%mm6
        movd  %mm6,(%eax,%edx,4)          ## increment vc[gid] 
       ## finish if last 
        movl nb102nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jecxz _nb_kernel102nf_ia32_3dnow.nb102nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb102nf_n(%esp)
        jmp _nb_kernel102nf_ia32_3dnow.nb102nf_outer
_nb_kernel102nf_ia32_3dnow.nb102nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb102nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jecxz _nb_kernel102nf_ia32_3dnow.nb102nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel102nf_ia32_3dnow.nb102nf_threadloop
_nb_kernel102nf_ia32_3dnow.nb102nf_end: 
        femms

        movl nb102nf_nouter(%esp),%eax
        movl nb102nf_ninner(%esp),%ebx
        movl nb102nf_outeriter(%ebp),%ecx
        movl nb102nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        addl $104,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret

