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



.globl nb_kernel112_ia32_3dnow
.globl _nb_kernel112_ia32_3dnow
nb_kernel112_ia32_3dnow:        
_nb_kernel112_ia32_3dnow:       
.set nb112_p_nri, 8
.set nb112_iinr, 12
.set nb112_jindex, 16
.set nb112_jjnr, 20
.set nb112_shift, 24
.set nb112_shiftvec, 28
.set nb112_fshift, 32
.set nb112_gid, 36
.set nb112_pos, 40
.set nb112_faction, 44
.set nb112_charge, 48
.set nb112_p_facel, 52
.set nb112_p_krf, 56
.set nb112_p_crf, 60
.set nb112_Vc, 64
.set nb112_type, 68
.set nb112_p_ntype, 72
.set nb112_vdwparam, 76
.set nb112_Vvdw, 80
.set nb112_p_tabscale, 84
.set nb112_VFtab, 88
.set nb112_invsqrta, 92
.set nb112_dvda, 96
.set nb112_p_gbtabscale, 100
.set nb112_GBtab, 104
.set nb112_p_nthreads, 108
.set nb112_count, 112
.set nb112_mtx, 116
.set nb112_outeriter, 120
.set nb112_inneriter, 124
.set nb112_work, 128
                        ## stack offsets for local variables 
.set nb112_is3, 0
.set nb112_ii3, 4
.set nb112_ixO, 8
.set nb112_iyO, 12
.set nb112_izO, 16
.set nb112_ixH, 20
.set nb112_iyH, 28
.set nb112_izH, 36
.set nb112_qqOO, 44
.set nb112_qqOH, 52
.set nb112_qqHH, 60
.set nb112_c6, 68
.set nb112_c12, 76
.set nb112_six, 84
.set nb112_twelve, 92
.set nb112_vctot, 100
.set nb112_Vvdwtot, 108
.set nb112_innerjjnr, 116
.set nb112_innerk, 120
.set nb112_fixO, 124
.set nb112_fiyO, 128
.set nb112_fizO, 132
.set nb112_fixH, 136
.set nb112_fiyH, 144
.set nb112_fizH, 152
.set nb112_dxO, 160
.set nb112_dyO, 164
.set nb112_dzO, 168
.set nb112_dxH, 172
.set nb112_dyH, 180
.set nb112_dzH, 188
.set nb112_n, 196                           ## idx for outer loop
.set nb112_nn1, 200                         ## number of outer iterations
.set nb112_nri, 204
.set nb112_ntype, 208
.set nb112_nouter, 212
.set nb112_ninner, 216
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $220,%esp          ## local stack space 
        femms

        movl nb112_p_nri(%ebp),%ecx
        movl nb112_p_ntype(%ebp),%edx
        movl nb112_p_facel(%ebp),%esi
        movl (%ecx),%ecx
        movl (%edx),%edx
        movl %ecx,nb112_nri(%esp)
        movl %edx,nb112_ntype(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb112_nouter(%esp)
        movl %eax,nb112_ninner(%esp)

        ## assume we have at least one i particle - start directly      

        movl  nb112_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx),%ebx                   ## ebx=ii 

        movl  nb112_charge(%ebp),%edx
        movd  (%esi),%mm1                   ## mm1=facel 
        movd  (%edx,%ebx,4),%mm2            ## mm2=charge[ii0] (O) 
        movd  4(%edx,%ebx,4),%mm3           ## mm2=charge[ii0+1] (H)  
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
        movq  %mm4,nb112_qqOO(%esp)
        movq  %mm5,nb112_qqOH(%esp)
        movq  %mm6,nb112_qqHH(%esp)
        movl  nb112_type(%ebp),%edx
        movl  (%edx,%ebx,4),%ecx
        shll  %ecx
        movl  %ecx,%edx
        imull nb112_ntype(%esp),%ecx
        addl  %ecx,%edx
        movl  nb112_vdwparam(%ebp),%eax
        movd  (%eax,%edx,4),%mm0
        movd  4(%eax,%edx,4),%mm1
        movq  %mm0,nb112_c6(%esp)
        movq  %mm1,nb112_c12(%esp)
        ## move data to local stack  
        movl $0x40c00000,%eax ## fp 6.0
        movl $0x41400000,%ebx ## fp 12.0

        movl %eax,nb112_six(%esp)
        movl %eax,nb112_six+4(%esp)
        movl %ebx,nb112_twelve(%esp)
        movl %ebx,nb112_twelve+4(%esp)

_nb_kernel112_ia32_3dnow.nb112_threadloop: 
        movl  nb112_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel112_ia32_3dnow.nb112_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel112_ia32_3dnow.nb112_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb112_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb112_n(%esp)
        movl %ebx,nb112_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel112_ia32_3dnow.nb112_outerstart
        jmp _nb_kernel112_ia32_3dnow.nb112_end

_nb_kernel112_ia32_3dnow.nb112_outerstart: 
        ## ebx contains number of outer iterations
        addl nb112_nouter(%esp),%ebx
        movl %ebx,nb112_nouter(%esp)

_nb_kernel112_ia32_3dnow.nb112_outer: 
        movl  nb112_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb112_is3(%esp)      ## store is3 

        movl  nb112_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movq  (%eax,%ebx,4),%mm5        ## move shX/shY to mm5 and shZ to mm6. 
        movd  8(%eax,%ebx,4),%mm6
        movq  %mm5,%mm0
        movq  %mm5,%mm1
        movq  %mm6,%mm2
        punpckldq %mm0,%mm0         ## also expand shX,Y,Z in mm0--mm2. 
        punpckhdq %mm1,%mm1
        punpckldq %mm2,%mm2

        movl  nb112_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx,%esi,4),%ebx            ## ebx=ii 

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb112_pos(%ebp),%eax      ## eax = base of pos[] 

        pfadd (%eax,%ebx,4),%mm5    ## ix = shX + posX (and iy too) 
        movd  8(%eax,%ebx,4),%mm7       ## cant use direct memory add for 4 bytes (iz) 
        movl  %ebx,nb112_ii3(%esp)          ## (use mm7 as temp. storage for iz.) 
        pfadd %mm7,%mm6
        movq  %mm5,nb112_ixO(%esp)
        movq  %mm6,nb112_izO(%esp)

        movd  12(%eax,%ebx,4),%mm3
        movd  16(%eax,%ebx,4),%mm4
        movd  20(%eax,%ebx,4),%mm5
        punpckldq  24(%eax,%ebx,4),%mm3
        punpckldq  28(%eax,%ebx,4),%mm4
        punpckldq  32(%eax,%ebx,4),%mm5    ## coords of H1 in low mm3-mm5, H2 in high 

        pfadd %mm3,%mm0
        pfadd %mm4,%mm1
        pfadd %mm5,%mm2
        movq %mm0,nb112_ixH(%esp)
        movq %mm1,nb112_iyH(%esp)
        movq %mm2,nb112_izH(%esp)

        ## clear vctot and i forces 
        pxor  %mm7,%mm7
        movq  %mm7,nb112_vctot(%esp)
        movq  %mm7,nb112_Vvdwtot(%esp)
        movq  %mm7,nb112_fixO(%esp)
        movq  %mm7,nb112_fizO(%esp)
        movq  %mm7,nb112_fixH(%esp)
        movq  %mm7,nb112_fiyH(%esp)
        movq  %mm7,nb112_fizH(%esp)

        movl  nb112_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 
        movl  %edx,nb112_innerk(%esp)      ## number of innerloop atoms 
        addl  nb112_ninner(%esp),%edx
        movl  %edx,nb112_ninner(%esp)

        movl  nb112_pos(%ebp),%esi
        movl  nb112_faction(%ebp),%edi
        movl  nb112_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb112_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
_nb_kernel112_ia32_3dnow.nb112_inner_loop: 
        ## a single j particle iteration here - compare with the unrolled code for comments. 
        movl  nb112_innerjjnr(%esp),%eax
        movl  (%eax),%eax       ## eax=jnr offset 
        addl $4,nb112_innerjjnr(%esp)             ## advance pointer 

        movd  nb112_qqOO(%esp),%mm6
        movq  nb112_qqOH(%esp),%mm7

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

        pfsubr nb112_ixO(%esp),%mm0
        pfsubr nb112_izO(%esp),%mm1

        movq  %mm0,nb112_dxO(%esp)
        pfmul %mm0,%mm0
        movd  %mm1,nb112_dzO(%esp)
        pfmul %mm1,%mm1
        pfacc %mm0,%mm0
        pfadd %mm1,%mm0         ## mm0=rsqO 

        punpckldq %mm2,%mm2
        punpckldq %mm3,%mm3
        punpckldq %mm4,%mm4 ## mm2-mm4 is jx-jz 
        pfsubr nb112_ixH(%esp),%mm2
        pfsubr nb112_iyH(%esp),%mm3
        pfsubr nb112_izH(%esp),%mm4   ## mm2-mm4 is dxH-dzH 

        movq %mm2,nb112_dxH(%esp)
        movq %mm3,nb112_dyH(%esp)
        movq %mm4,nb112_dzH(%esp)
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

        movq %mm4,%mm2
        pfmul %mm4,%mm2
        pfmul %mm4,%mm2
        movq %mm2,%mm0
        pfmul %mm0,%mm0
        pfmul nb112_c6(%esp),%mm2
        pfmul nb112_c12(%esp),%mm0
        movq %mm0,%mm5
        pfsub %mm2,%mm5         ## Vvdw 

        pfmul nb112_six(%esp),%mm2
        pfmul nb112_twelve(%esp),%mm0

        pfsub %mm2,%mm0

        ## calculate potential and scalar force 
        pfmul %mm1,%mm6         ## mm6=vcoul 
        pfadd %mm6,%mm0
        pfmul %mm0,%mm4         ## mm4=fscalar  

        ## update nb potential 
        pfadd nb112_Vvdwtot(%esp),%mm5
        movq %mm5,nb112_Vvdwtot(%esp)

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
        pfmul %mm7,%mm3         ## mm3=fscal for the two H's. 

        ## update vctot 
        pfadd %mm6,%mm7
        pfadd nb112_vctot(%esp),%mm7
        movq %mm7,nb112_vctot(%esp)

        ## spread oxygen fscalar to both positions 
        punpckldq %mm4,%mm4
        ## calc vectorial force for O 
        movq nb112_dxO(%esp),%mm0
        movd nb112_dzO(%esp),%mm1
        pfmul %mm4,%mm0
        pfmul %mm4,%mm1

        ## calc vectorial force for H's 
        movq nb112_dxH(%esp),%mm5
        movq nb112_dyH(%esp),%mm6
        movq nb112_dzH(%esp),%mm7
        pfmul %mm3,%mm5
        pfmul %mm3,%mm6
        pfmul %mm3,%mm7

        ## update iO particle force 
        movq nb112_fixO(%esp),%mm2
        movd nb112_fizO(%esp),%mm3
        pfadd %mm0,%mm2
        pfadd %mm1,%mm3
        movq %mm2,nb112_fixO(%esp)
        movd %mm3,nb112_fizO(%esp)

        ## update iH forces 
        movq nb112_fixH(%esp),%mm2
        movq nb112_fiyH(%esp),%mm3
        movq nb112_fizH(%esp),%mm4
        pfadd %mm5,%mm2
        pfadd %mm6,%mm3
        pfadd %mm7,%mm4
        movq %mm2,nb112_fixH(%esp)
        movq %mm3,nb112_fiyH(%esp)
        movq %mm4,nb112_fizH(%esp)

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

        movd nb112_qqOH(%esp),%mm6
        movq nb112_qqHH(%esp),%mm7

        pfsubr nb112_ixO(%esp),%mm0
        pfsubr nb112_izO(%esp),%mm1

        movq  %mm0,nb112_dxO(%esp)
        pfmul %mm0,%mm0
        movd  %mm1,nb112_dzO(%esp)
        pfmul %mm1,%mm1
        pfacc %mm1,%mm0
        pfadd %mm1,%mm0         ## mm0=rsqO 

        punpckldq %mm2,%mm2
        punpckldq %mm3,%mm3
        punpckldq %mm4,%mm4 ## mm2-mm4 is jx-jz 
        pfsubr nb112_ixH(%esp),%mm2
        pfsubr nb112_iyH(%esp),%mm3
        pfsubr nb112_izH(%esp),%mm4   ## mm2-mm4 is dxH-dzH 

        movq %mm2,nb112_dxH(%esp)
        movq %mm3,nb112_dyH(%esp)
        movq %mm4,nb112_dzH(%esp)
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
        pfmul %mm7,%mm3         ## mm3=fscal for the two H's. 

        ## update vctot 
        pfadd %mm6,%mm7
        pfadd nb112_vctot(%esp),%mm7
        movq %mm7,nb112_vctot(%esp)

        ## spread oxygen fscalar to both positions 
        punpckldq %mm4,%mm4
        ## calc vectorial force for O 
        movq nb112_dxO(%esp),%mm0
        movd nb112_dzO(%esp),%mm1
        pfmul %mm4,%mm0
        pfmul %mm4,%mm1

        ## calc vectorial force for H's 
        movq nb112_dxH(%esp),%mm5
        movq nb112_dyH(%esp),%mm6
        movq nb112_dzH(%esp),%mm7
        pfmul %mm3,%mm5
        pfmul %mm3,%mm6
        pfmul %mm3,%mm7

        ## update iO particle force 
        movq nb112_fixO(%esp),%mm2
        movd nb112_fizO(%esp),%mm3
        pfadd %mm0,%mm2
        pfadd %mm1,%mm3
        movq %mm2,nb112_fixO(%esp)
        movd %mm3,nb112_fizO(%esp)

        ## update iH forces 
        movq nb112_fixH(%esp),%mm2
        movq nb112_fiyH(%esp),%mm3
        movq nb112_fizH(%esp),%mm4
        pfadd %mm5,%mm2
        pfadd %mm6,%mm3
        pfadd %mm7,%mm4
        movq %mm2,nb112_fixH(%esp)
        movq %mm3,nb112_fiyH(%esp)
        movq %mm4,nb112_fizH(%esp)

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

        movd nb112_qqOH(%esp),%mm6
        movq nb112_qqHH(%esp),%mm7

        pfsubr nb112_ixO(%esp),%mm0
        pfsubr nb112_izO(%esp),%mm1

        movq  %mm0,nb112_dxO(%esp)
        pfmul %mm0,%mm0
        movd  %mm1,nb112_dzO(%esp)
        pfmul %mm1,%mm1
        pfacc %mm1,%mm0
        pfadd %mm1,%mm0         ## mm0=rsqO 

        punpckldq %mm2,%mm2
        punpckldq %mm3,%mm3
        punpckldq %mm4,%mm4 ## mm2-mm4 is jx-jz 
        pfsubr nb112_ixH(%esp),%mm2
        pfsubr nb112_iyH(%esp),%mm3
        pfsubr nb112_izH(%esp),%mm4   ## mm2-mm4 is dxH-dzH 

        movq %mm2,nb112_dxH(%esp)
        movq %mm3,nb112_dyH(%esp)
        movq %mm4,nb112_dzH(%esp)
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
        punpckldq %mm2,%mm5     ## seeds are in mm5 now, and rsq in mm3. 

        movq %mm5,%mm2
        pfmul %mm5,%mm5
        pfrsqit1 %mm3,%mm5
        pfrcpit2 %mm2,%mm5      ## mm5=invsqrt 
        movq %mm5,%mm3
        pfmul %mm3,%mm3         ## mm3=invsq 
        pfmul %mm5,%mm7         ## mm7=vcoul 
        pfmul %mm7,%mm3         ## mm3=fscal for the two H's. 

        ## update vctot 
        pfadd %mm6,%mm7
        pfadd nb112_vctot(%esp),%mm7
        movq %mm7,nb112_vctot(%esp)

        ## spread oxygen fscalar to both positions 
        punpckldq %mm4,%mm4
        ## calc vectorial force for O 
        movq nb112_dxO(%esp),%mm0
        movd nb112_dzO(%esp),%mm1
        pfmul %mm4,%mm0
        pfmul %mm4,%mm1

        ## calc vectorial force for H's 
        movq nb112_dxH(%esp),%mm5
        movq nb112_dyH(%esp),%mm6
        movq nb112_dzH(%esp),%mm7
        pfmul %mm3,%mm5
        pfmul %mm3,%mm6
        pfmul %mm3,%mm7

        ## update iO particle force 
        movq nb112_fixO(%esp),%mm2
        movd nb112_fizO(%esp),%mm3
        pfadd %mm0,%mm2
        pfadd %mm1,%mm3
        movq %mm2,nb112_fixO(%esp)
        movd %mm3,nb112_fizO(%esp)

        ## update iH forces 
        movq nb112_fixH(%esp),%mm2
        movq nb112_fiyH(%esp),%mm3
        movq nb112_fizH(%esp),%mm4
        pfadd %mm5,%mm2
        pfadd %mm6,%mm3
        pfadd %mm7,%mm4
        movq %mm2,nb112_fixH(%esp)
        movq %mm3,nb112_fiyH(%esp)
        movq %mm4,nb112_fizH(%esp)

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
        decl nb112_innerk(%esp)
        jz  _nb_kernel112_ia32_3dnow.nb112_updateouterdata
        jmp _nb_kernel112_ia32_3dnow.nb112_inner_loop
_nb_kernel112_ia32_3dnow.nb112_updateouterdata: 
        movl  nb112_ii3(%esp),%ecx

        movq  (%edi,%ecx,4),%mm6       ## increment iO force  
        movd  8(%edi,%ecx,4),%mm7
        pfadd nb112_fixO(%esp),%mm6
        pfadd nb112_fizO(%esp),%mm7
        movq  %mm6,(%edi,%ecx,4)
        movd  %mm7,8(%edi,%ecx,4)

        movq  nb112_fixH(%esp),%mm0
        movq  nb112_fiyH(%esp),%mm3
        movq  nb112_fizH(%esp),%mm1
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


        movl  nb112_fshift(%ebp),%ebx      ## increment fshift force 
        movl  nb112_is3(%esp),%edx

        movq  (%ebx,%edx,4),%mm6
        movd  8(%ebx,%edx,4),%mm7
        pfadd nb112_fixO(%esp),%mm6
        pfadd nb112_fizO(%esp),%mm7
        pfadd %mm0,%mm6
        pfadd %mm1,%mm7
        pfadd %mm2,%mm6
        pfadd %mm3,%mm7
        movq  %mm6,(%ebx,%edx,4)
        movd  %mm7,8(%ebx,%edx,4)

        ## get n from stack
        movl nb112_n(%esp),%esi
        ## get group index for i particle 
        movl  nb112_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        movq  nb112_vctot(%esp),%mm7
        pfacc %mm7,%mm7           ## get and sum the two parts of total potential 

        movl  nb112_Vc(%ebp),%eax
        movd  (%eax,%edx,4),%mm6
        pfadd %mm7,%mm6
        movd  %mm6,(%eax,%edx,4)          ## increment vc[gid] 

        movq  nb112_Vvdwtot(%esp),%mm7
        pfacc %mm7,%mm7           ## get and sum the two parts of total potential 

        movl  nb112_Vvdw(%ebp),%eax
        movd  (%eax,%edx,4),%mm6
        pfadd %mm7,%mm6
        movd  %mm6,(%eax,%edx,4)          ## increment Vvdwtot[gid] 
       ## finish if last 
        movl nb112_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jecxz _nb_kernel112_ia32_3dnow.nb112_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb112_n(%esp)
        jmp _nb_kernel112_ia32_3dnow.nb112_outer
_nb_kernel112_ia32_3dnow.nb112_outerend: 
        ## check if more outer neighborlists remain
        movl  nb112_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jecxz _nb_kernel112_ia32_3dnow.nb112_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel112_ia32_3dnow.nb112_threadloop
_nb_kernel112_ia32_3dnow.nb112_end: 
        femms

        movl nb112_nouter(%esp),%eax
        movl nb112_ninner(%esp),%ebx
        movl nb112_outeriter(%ebp),%ecx
        movl nb112_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        addl $220,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret


.globl nb_kernel112nf_ia32_3dnow
.globl _nb_kernel112nf_ia32_3dnow
nb_kernel112nf_ia32_3dnow:      
_nb_kernel112nf_ia32_3dnow:     
.set nb112nf_p_nri, 8
.set nb112nf_iinr, 12
.set nb112nf_jindex, 16
.set nb112nf_jjnr, 20
.set nb112nf_shift, 24
.set nb112nf_shiftvec, 28
.set nb112nf_fshift, 32
.set nb112nf_gid, 36
.set nb112nf_pos, 40
.set nb112nf_faction, 44
.set nb112nf_charge, 48
.set nb112nf_p_facel, 52
.set nb112nf_p_krf, 56
.set nb112nf_p_crf, 60
.set nb112nf_Vc, 64
.set nb112nf_type, 68
.set nb112nf_p_ntype, 72
.set nb112nf_vdwparam, 76
.set nb112nf_Vvdw, 80
.set nb112nf_p_tabscale, 84
.set nb112nf_VFtab, 88
.set nb112nf_invsqrta, 92
.set nb112nf_dvda, 96
.set nb112nf_p_gbtabscale, 100
.set nb112nf_GBtab, 104
.set nb112nf_p_nthreads, 108
.set nb112nf_count, 112
.set nb112nf_mtx, 116
.set nb112nf_outeriter, 120
.set nb112nf_inneriter, 124
.set nb112nf_work, 128
                        ## stack offsets for local variables 
.set nb112nf_is3, 0
.set nb112nf_ii3, 4
.set nb112nf_ixO, 8
.set nb112nf_iyO, 12
.set nb112nf_izO, 16
.set nb112nf_ixH, 20
.set nb112nf_iyH, 28
.set nb112nf_izH, 36
.set nb112nf_qqOO, 44
.set nb112nf_qqOH, 52
.set nb112nf_qqHH, 60
.set nb112nf_c6, 68
.set nb112nf_c12, 76
.set nb112nf_vctot, 84
.set nb112nf_Vvdwtot, 92
.set nb112nf_innerjjnr, 100
.set nb112nf_innerk, 104
.set nb112nf_n, 108                         ## idx for outer loop
.set nb112nf_nn1, 112                       ## number of outer iterations
.set nb112nf_nri, 116
.set nb112nf_ntype, 120
.set nb112nf_nouter, 124
.set nb112nf_ninner, 128
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $132,%esp          ## local stack space 
        femms
        movl nb112nf_p_nri(%ebp),%ecx
        movl nb112nf_p_ntype(%ebp),%edx
        movl nb112nf_p_facel(%ebp),%esi
        movl (%ecx),%ecx
        movl (%edx),%edx
        movl %ecx,nb112nf_nri(%esp)
        movl %edx,nb112nf_ntype(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb112nf_nouter(%esp)
        movl %eax,nb112nf_ninner(%esp)

        ## assume we have at least one i particle - start directly      

        movl  nb112nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx),%ebx           ## ebx=ii 

        movl  nb112nf_charge(%ebp),%edx
        movd  (%esi),%mm1               ## mm1=facel 
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
        movq  %mm4,nb112nf_qqOO(%esp)
        movq  %mm5,nb112nf_qqOH(%esp)
        movq  %mm6,nb112nf_qqHH(%esp)
        movl  nb112nf_type(%ebp),%edx
        movl  (%edx,%ebx,4),%ecx
        shll  %ecx
        movl  %ecx,%edx
        imull nb112nf_ntype(%esp),%ecx
        addl  %ecx,%edx
        movl  nb112nf_vdwparam(%ebp),%eax
        movd  (%eax,%edx,4),%mm0
        movd  4(%eax,%edx,4),%mm1
        movq  %mm0,nb112nf_c6(%esp)
        movq  %mm1,nb112nf_c12(%esp)

_nb_kernel112nf_ia32_3dnow.nb112nf_threadloop: 
        movl  nb112nf_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel112nf_ia32_3dnow.nb112nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel112nf_ia32_3dnow.nb112nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb112nf_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb112nf_n(%esp)
        movl %ebx,nb112nf_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel112nf_ia32_3dnow.nb112nf_outerstart
        jmp _nb_kernel112nf_ia32_3dnow.nb112nf_end

_nb_kernel112nf_ia32_3dnow.nb112nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb112nf_nouter(%esp),%ebx
        movl %ebx,nb112nf_nouter(%esp)

_nb_kernel112nf_ia32_3dnow.nb112nf_outer: 
        movl  nb112nf_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb112nf_is3(%esp)            ## store is3 

        movl  nb112nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movq  (%eax,%ebx,4),%mm5        ## move shX/shY to mm5 and shZ to mm6. 
        movd  8(%eax,%ebx,4),%mm6
        movq  %mm5,%mm0
        movq  %mm5,%mm1
        movq  %mm6,%mm2
        punpckldq %mm0,%mm0         ## also expand shX,Y,Z in mm0--mm2. 
        punpckhdq %mm1,%mm1
        punpckldq %mm2,%mm2

        movl  nb112nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx,%esi,4),%ebx            ## ebx=ii 

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb112nf_pos(%ebp),%eax      ## eax = base of pos[] 

        pfadd (%eax,%ebx,4),%mm5    ## ix = shX + posX (and iy too) 
        movd  8(%eax,%ebx,4),%mm7       ## cant use direct memory add for 4 bytes (iz) 
        movl  %ebx,nb112nf_ii3(%esp)        ## (use mm7 as temp. storage for iz.) 
        pfadd %mm7,%mm6
        movq  %mm5,nb112nf_ixO(%esp)
        movq  %mm6,nb112nf_izO(%esp)

        movd  12(%eax,%ebx,4),%mm3
        movd  16(%eax,%ebx,4),%mm4
        movd  20(%eax,%ebx,4),%mm5
        punpckldq  24(%eax,%ebx,4),%mm3
        punpckldq  28(%eax,%ebx,4),%mm4
        punpckldq  32(%eax,%ebx,4),%mm5    ## coords of H1 in low mm3-mm5, H2 in high 

        pfadd %mm3,%mm0
        pfadd %mm4,%mm1
        pfadd %mm5,%mm2
        movq %mm0,nb112nf_ixH(%esp)
        movq %mm1,nb112nf_iyH(%esp)
        movq %mm2,nb112nf_izH(%esp)

        ## clear vctot and i forces 
        pxor  %mm7,%mm7
        movq  %mm7,nb112nf_vctot(%esp)
        movq  %mm7,nb112nf_Vvdwtot(%esp)

        movl  nb112nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 
        movl  %edx,nb112nf_innerk(%esp)      ## number of innerloop atoms 
        addl  nb112nf_ninner(%esp),%edx
        movl  %edx,nb112nf_ninner(%esp)

        movl  nb112nf_pos(%ebp),%esi
        movl  nb112nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb112nf_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
_nb_kernel112nf_ia32_3dnow.nb112nf_inner_loop: 
        ## a single j particle iteration here - compare with the unrolled code for comments. 
        movl  nb112nf_innerjjnr(%esp),%eax
        movl  (%eax),%eax       ## eax=jnr offset 
        addl $4,nb112nf_innerjjnr(%esp)             ## advance pointer 

        movd  nb112nf_qqOO(%esp),%mm6
        movq  nb112nf_qqOH(%esp),%mm7

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

        pfsubr nb112nf_ixO(%esp),%mm0
        pfsubr nb112nf_izO(%esp),%mm1

        pfmul %mm0,%mm0
        pfmul %mm1,%mm1
        pfacc %mm0,%mm0
        pfadd %mm1,%mm0         ## mm0=rsqO 

        punpckldq %mm2,%mm2
        punpckldq %mm3,%mm3
        punpckldq %mm4,%mm4 ## mm2-mm4 is jx-jz 
        pfsubr nb112nf_ixH(%esp),%mm2
        pfsubr nb112nf_iyH(%esp),%mm3
        pfsubr nb112nf_izH(%esp),%mm4   ## mm2-mm4 is dxH-dzH 

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

        movq %mm4,%mm2
        pfmul %mm4,%mm2
        pfmul %mm4,%mm2
        movq %mm2,%mm0
        pfmul %mm0,%mm0
        pfmul nb112nf_c6(%esp),%mm2
        pfmul nb112nf_c12(%esp),%mm0
        movq %mm0,%mm5
        pfsub %mm2,%mm5         ## Vvdw 

        ## calculate potential and scalar force 
        pfmul %mm1,%mm6         ## mm6=vcoul 
        ## update nb potential 
        pfadd nb112nf_Vvdwtot(%esp),%mm5
        movq %mm5,nb112nf_Vvdwtot(%esp)

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
        pfadd nb112nf_vctot(%esp),%mm7
        movq %mm7,nb112nf_vctot(%esp)

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

        movd nb112nf_qqOH(%esp),%mm6
        movq nb112nf_qqHH(%esp),%mm7

        pfsubr nb112nf_ixO(%esp),%mm0
        pfsubr nb112nf_izO(%esp),%mm1

        pfmul %mm0,%mm0
        pfmul %mm1,%mm1
        pfacc %mm1,%mm0
        pfadd %mm1,%mm0         ## mm0=rsqO 

        punpckldq %mm2,%mm2
        punpckldq %mm3,%mm3
        punpckldq %mm4,%mm4 ## mm2-mm4 is jx-jz 
        pfsubr nb112nf_ixH(%esp),%mm2
        pfsubr nb112nf_iyH(%esp),%mm3
        pfsubr nb112nf_izH(%esp),%mm4   ## mm2-mm4 is dxH-dzH 

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
        pfadd nb112nf_vctot(%esp),%mm7
        movq %mm7,nb112nf_vctot(%esp)

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

        movd nb112nf_qqOH(%esp),%mm6
        movq nb112nf_qqHH(%esp),%mm7

        pfsubr nb112nf_ixO(%esp),%mm0
        pfsubr nb112nf_izO(%esp),%mm1

        pfmul %mm0,%mm0
        pfmul %mm1,%mm1
        pfacc %mm1,%mm0
        pfadd %mm1,%mm0         ## mm0=rsqO 

        punpckldq %mm2,%mm2
        punpckldq %mm3,%mm3
        punpckldq %mm4,%mm4 ## mm2-mm4 is jx-jz 
        pfsubr nb112nf_ixH(%esp),%mm2
        pfsubr nb112nf_iyH(%esp),%mm3
        pfsubr nb112nf_izH(%esp),%mm4   ## mm2-mm4 is dxH-dzH 

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
        punpckldq %mm2,%mm5     ## seeds are in mm5 now, and rsq in mm3. 

        movq %mm5,%mm2
        pfmul %mm5,%mm5
        pfrsqit1 %mm3,%mm5
        pfrcpit2 %mm2,%mm5      ## mm5=invsqrt 
        pfmul %mm5,%mm7         ## mm7=vcoul 

        ## update vctot 
        pfadd %mm6,%mm7
        pfadd nb112nf_vctot(%esp),%mm7
        movq %mm7,nb112nf_vctot(%esp)

        ##  done  - one more? 
        decl nb112nf_innerk(%esp)
        jz  _nb_kernel112nf_ia32_3dnow.nb112nf_updateouterdata
        jmp _nb_kernel112nf_ia32_3dnow.nb112nf_inner_loop
_nb_kernel112nf_ia32_3dnow.nb112nf_updateouterdata: 
        ## get n from stack
        movl nb112nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb112nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        movq  nb112nf_vctot(%esp),%mm7
        pfacc %mm7,%mm7           ## get and sum the two parts of total potential 

        movl  nb112nf_Vc(%ebp),%eax
        movd  (%eax,%edx,4),%mm6
        pfadd %mm7,%mm6
        movd  %mm6,(%eax,%edx,4)          ## increment vc[gid] 

        movq  nb112nf_Vvdwtot(%esp),%mm7
        pfacc %mm7,%mm7           ## get and sum the two parts of total potential 

        movl  nb112nf_Vvdw(%ebp),%eax
        movd  (%eax,%edx,4),%mm6
        pfadd %mm7,%mm6
        movd  %mm6,(%eax,%edx,4)          ## increment Vvdwtot[gid] 
        ## finish if last 
        movl nb112nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jecxz _nb_kernel112nf_ia32_3dnow.nb112nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb112nf_n(%esp)
        jmp _nb_kernel112nf_ia32_3dnow.nb112nf_outer
_nb_kernel112nf_ia32_3dnow.nb112nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb112nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jecxz _nb_kernel112nf_ia32_3dnow.nb112nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel112nf_ia32_3dnow.nb112nf_threadloop
_nb_kernel112nf_ia32_3dnow.nb112nf_end: 
        femms

        movl nb112nf_nouter(%esp),%eax
        movl nb112nf_ninner(%esp),%ebx
        movl nb112nf_outeriter(%ebp),%ecx
        movl nb112nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        addl $132,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret


