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



.globl nb_kernel332_ia32_3dnow
.globl _nb_kernel332_ia32_3dnow
nb_kernel332_ia32_3dnow:        
_nb_kernel332_ia32_3dnow:       
.set nb332_p_nri, 8
.set nb332_iinr, 12
.set nb332_jindex, 16
.set nb332_jjnr, 20
.set nb332_shift, 24
.set nb332_shiftvec, 28
.set nb332_fshift, 32
.set nb332_gid, 36
.set nb332_pos, 40
.set nb332_faction, 44
.set nb332_charge, 48
.set nb332_p_facel, 52
.set nb332_p_krf, 56
.set nb332_p_crf, 60
.set nb332_Vc, 64
.set nb332_type, 68
.set nb332_p_ntype, 72
.set nb332_vdwparam, 76
.set nb332_Vvdw, 80
.set nb332_p_tabscale, 84
.set nb332_VFtab, 88
.set nb332_invsqrta, 92
.set nb332_dvda, 96
.set nb332_p_gbtabscale, 100
.set nb332_GBtab, 104
.set nb332_p_nthreads, 108
.set nb332_count, 112
.set nb332_mtx, 116
.set nb332_outeriter, 120
.set nb332_inneriter, 124
.set nb332_work, 128
                        ## stack offsets for local variables 
.set nb332_is3, 0
.set nb332_ii3, 4
.set nb332_ixO, 8
.set nb332_iyO, 12
.set nb332_izO, 16
.set nb332_ixH, 20
.set nb332_iyH, 28
.set nb332_izH, 36
.set nb332_qqOO, 44
.set nb332_qqOH, 52
.set nb332_qqHH, 60
.set nb332_c6, 68
.set nb332_c12, 76
.set nb332_two, 84
.set nb332_n1, 92
.set nb332_tsc, 100
.set nb332_vctot, 108
.set nb332_Vvdwtot, 116
.set nb332_innerjjnr, 124
.set nb332_innerk, 128
.set nb332_fixO, 132
.set nb332_fiyO, 136
.set nb332_fizO, 140
.set nb332_fixH, 144
.set nb332_fiyH, 152
.set nb332_fizH, 160
.set nb332_dxO, 168
.set nb332_dyO, 172
.set nb332_dzO, 176
.set nb332_dxH, 180
.set nb332_dyH, 188
.set nb332_dzH, 196
.set nb332_tmprsqH, 204
.set nb332_n, 212                           ## idx for outer loop
.set nb332_nn1, 216                         ## number of outer iterations
.set nb332_nri, 220
.set nb332_nouter, 224
.set nb332_ninner, 228
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $232,%esp          ## local stack space 
        femms
        movl nb332_p_nri(%ebp),%ecx
        movl nb332_p_facel(%ebp),%esi
        movl nb332_p_tabscale(%ebp),%edi
        movl (%ecx),%ecx
        movl %ecx,nb332_nri(%esp)
        movd  (%esi),%mm1       ## facel
        movl nb332_p_ntype(%ebp),%esi

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb332_nouter(%esp)
        movl %eax,nb332_ninner(%esp)

        ## assume we have at least one i particle - start directly      

        movl  nb332_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx),%ebx           ## ebx=ii 

        movl  nb332_charge(%ebp),%edx
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
        movq  %mm4,nb332_qqOO(%esp)
        movq  %mm5,nb332_qqOH(%esp)
        movq  %mm6,nb332_qqHH(%esp)
        movl  nb332_type(%ebp),%edx
        movl  (%edx,%ebx,4),%ecx
        shll  %ecx
        movl  %ecx,%edx
        imull (%esi),%ecx
        addl  %ecx,%edx
        movl  nb332_vdwparam(%ebp),%eax
        movd  (%eax,%edx,4),%mm0
        movd  4(%eax,%edx,4),%mm1
        movq  %mm0,nb332_c6(%esp)
        movq  %mm1,nb332_c12(%esp)
        movd  (%edi),%mm3
        punpckldq %mm3,%mm3
        movq  %mm3,nb332_tsc(%esp)
        movl $0x40000000,%eax
        movl %eax,nb332_two(%esp)
        movl %eax,nb332_two+4(%esp)
_nb_kernel332_ia32_3dnow.nb332_threadloop: 
        movl  nb332_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel332_ia32_3dnow.nb332_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel332_ia32_3dnow.nb332_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb332_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb332_n(%esp)
        movl %ebx,nb332_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel332_ia32_3dnow.nb332_outerstart
        jmp _nb_kernel332_ia32_3dnow.nb332_end

_nb_kernel332_ia32_3dnow.nb332_outerstart: 
        ## ebx contains number of outer iterations
        addl nb332_nouter(%esp),%ebx
        movl %ebx,nb332_nouter(%esp)

_nb_kernel332_ia32_3dnow.nb332_outer: 
        movl  nb332_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb332_is3(%esp)      ## store is3 

        movl  nb332_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movq  (%eax,%ebx,4),%mm5        ## move shX/shY to mm5 and shZ to mm6. 
        movd  8(%eax,%ebx,4),%mm6
        movq  %mm5,%mm0
        movq  %mm5,%mm1
        movq  %mm6,%mm2
        punpckldq %mm0,%mm0         ## also expand shX,Y,Z in mm0--mm2. 
        punpckhdq %mm1,%mm1
        punpckldq %mm2,%mm2

        movl  nb332_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx,%esi,4),%ebx    ## ebx=ii 

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb332_pos(%ebp),%eax      ## eax = base of pos[] 

        pfadd (%eax,%ebx,4),%mm5    ## ix = shX + posX (and iy too) 
        movd  8(%eax,%ebx,4),%mm7       ## cant use direct memory add for 4 bytes (iz) 
        movl  %ebx,nb332_ii3(%esp)          ## (use mm7 as temp. storage for iz.) 
        pfadd %mm7,%mm6
        movq  %mm5,nb332_ixO(%esp)
        movq  %mm6,nb332_izO(%esp)

        movd  12(%eax,%ebx,4),%mm3
        movd  16(%eax,%ebx,4),%mm4
        movd  20(%eax,%ebx,4),%mm5
        punpckldq  24(%eax,%ebx,4),%mm3
        punpckldq  28(%eax,%ebx,4),%mm4
        punpckldq  32(%eax,%ebx,4),%mm5    ## coords of H1 in low mm3-mm5, H2 in high 

        pfadd %mm3,%mm0
        pfadd %mm4,%mm1
        pfadd %mm5,%mm2
        movq %mm0,nb332_ixH(%esp)
        movq %mm1,nb332_iyH(%esp)
        movq %mm2,nb332_izH(%esp)

        ## clear vctot and i forces 
        pxor  %mm7,%mm7
        movq  %mm7,nb332_vctot(%esp)
        movq  %mm7,nb332_Vvdwtot(%esp)
        movq  %mm7,nb332_fixO(%esp)
        movq  %mm7,nb332_fizO(%esp)
        movq  %mm7,nb332_fixH(%esp)
        movq  %mm7,nb332_fiyH(%esp)
        movq  %mm7,nb332_fizH(%esp)

        movl  nb332_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx     ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 
        movl  %edx,nb332_innerk(%esp)
        addl  nb332_ninner(%esp),%edx
        movl  %edx,nb332_ninner(%esp)

        movl  nb332_pos(%ebp),%esi
        movl  nb332_faction(%ebp),%edi
        movl  nb332_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb332_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
_nb_kernel332_ia32_3dnow.nb332_inner_loop: 
        ## a single j particle iteration here - compare with the unrolled code for comments. 
        movl  nb332_innerjjnr(%esp),%eax
        movl  (%eax),%eax       ## eax=jnr offset 
    addl $4,nb332_innerjjnr(%esp)             ## advance pointer 

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

        pfsubr nb332_ixO(%esp),%mm0
        pfsubr nb332_izO(%esp),%mm1

        movq  %mm0,nb332_dxO(%esp)
        pfmul %mm0,%mm0
        movd  %mm1,nb332_dzO(%esp)
        pfmul %mm1,%mm1
        pfacc %mm0,%mm0
        pfadd %mm1,%mm0         ## mm0=rsqO 

        punpckldq %mm2,%mm2
        punpckldq %mm3,%mm3
        punpckldq %mm4,%mm4 ## mm2-mm4 is jx-jz 
        pfsubr nb332_ixH(%esp),%mm2
        pfsubr nb332_iyH(%esp),%mm3
        pfsubr nb332_izH(%esp),%mm4   ## mm2-mm4 is dxH-dzH 

        movq %mm2,nb332_dxH(%esp)
        movq %mm3,nb332_dyH(%esp)
        movq %mm4,nb332_dzH(%esp)
        pfmul %mm2,%mm2
        pfmul %mm3,%mm3
        pfmul %mm4,%mm4

        pfadd %mm2,%mm3
        pfadd %mm4,%mm3         ## mm3=rsqH 
        movq %mm3,nb332_tmprsqH(%esp)

    pfrsqrt %mm0,%mm1

    movq %mm1,%mm2
    pfmul %mm1,%mm1
    pfrsqit1 %mm0,%mm1
    pfrcpit2 %mm2,%mm1  ## mm1=invsqrt  
        pfmul %mm1,%mm0         ## mm0=rsq  

        pfmul nb332_tsc(%esp),%mm0
        pf2iw %mm0,%mm4
        movd %mm4,nb332_n1(%esp)
        pi2fd %mm4,%mm4
        pfsub %mm4,%mm0              ## now mm0 is eps and mm4 n0 
        movq  %mm0,%mm2
        pfmul %mm2,%mm2         ## mm0 is eps, mm2 eps2 

        ## coulomb table 
        movl nb332_VFtab(%ebp),%edx
        movl nb332_n1(%esp),%ecx
        leal (%ecx,%ecx,2),%ecx
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

        pfmul nb332_two(%esp),%mm7      ## two*Heps2 
        pfadd %mm6,%mm7
        pfadd %mm5,%mm7 ## mm7=FF 

        pfmul %mm0,%mm5 ## mm5=eps*Fp 
        pfadd %mm4,%mm5 ##  mm5= VV 

        pfmul nb332_qqOO(%esp),%mm5     ## vcoul=qq*VV 
        pfmul nb332_qqOO(%esp),%mm7     ## fijC=qq*FF 

        ## update vctot directly, use mm3 for fscal sum. 
        pfadd nb332_vctot(%esp),%mm5
        movq %mm5,nb332_vctot(%esp)
        movq %mm7,%mm3

        ## dispersion table 
        ## load all the table values we need 
        movd 16(%edx,%ecx,4),%mm4
        movd 20(%edx,%ecx,4),%mm5
        movd 24(%edx,%ecx,4),%mm6
        movd 28(%edx,%ecx,4),%mm7
        pfmul %mm0,%mm6 ## mm6 = Geps           
        pfmul %mm2,%mm7 ## mm7 = Heps2 
        pfadd %mm6,%mm5
        pfadd %mm7,%mm5 ## mm5 = Fp 
        pfmul nb332_two(%esp),%mm7      ## two*Heps2 
        pfadd %mm6,%mm7
        pfadd %mm5,%mm7 ## mm7=FF 
        pfmul %mm0,%mm5 ## mm5=eps*Fp 
        pfadd %mm4,%mm5 ##  mm5= VV     

        movq nb332_c6(%esp),%mm4
        pfmul %mm4,%mm7 ## fijD 
        pfmul %mm4,%mm5 ## Vvdw6            
        pfadd %mm7,%mm3 ## add to fscal  

        ## update Vvdwtot to release mm5! 
        pfadd nb332_Vvdwtot(%esp),%mm5        ## add the earlier value 
        movq %mm5,nb332_Vvdwtot(%esp)         ## store the sum       

        ## repulsion table 
        ## load all the table values we need 
        movd 32(%edx,%ecx,4),%mm4
        movd 36(%edx,%ecx,4),%mm5
        movd 40(%edx,%ecx,4),%mm6
        movd 44(%edx,%ecx,4),%mm7

        pfmul %mm0,%mm6 ## mm6 = Geps           
        pfmul %mm2,%mm7 ## mm7 = Heps2 
        pfadd %mm6,%mm5
        pfadd %mm7,%mm5 ## mm5 = Fp 
        pfmul nb332_two(%esp),%mm7      ## two*Heps2 
        pfadd %mm6,%mm7
        pfadd %mm5,%mm7 ## mm7=FF 
        pfmul %mm0,%mm5 ## mm5=eps*Fp 
        pfadd %mm4,%mm5 ##  mm5= VV 

        movq nb332_c12(%esp),%mm6
        pfmul %mm6,%mm7 ## fijR 
        pfmul %mm6,%mm5 ## Vvdw12 
        pfadd %mm7,%mm3 ## total fscal fijC+ fijD+ fijR 

        ## change sign of fscal and multiply with rinv  
    pxor %mm0,%mm0
        pfsubr %mm0,%mm3
        pfmul nb332_tsc(%esp),%mm3
        pfmul %mm1,%mm3   ## mm3 is total fscal (for the oxygen) now 

        ## update Vvdwtot  
        pfadd nb332_Vvdwtot(%esp),%mm5        ## add the earlier value 
        movq %mm5,nb332_Vvdwtot(%esp)         ## store the sum       

        ## Ready with the oxygen - potential is updated, fscal is in mm3.
         ## time for hydrogens!


        movq nb332_tmprsqH(%esp),%mm0

        pfrsqrt %mm0,%mm1
        pswapd %mm0,%mm0
        pfrsqrt %mm0,%mm2
        pswapd %mm0,%mm0
        punpckldq %mm2,%mm1     ## seeds are in mm1 now, and rsq in mm0. 

        movq %mm1,%mm2
        pfmul %mm1,%mm1
    pfrsqit1 %mm0,%mm1
    pfrcpit2 %mm2,%mm1  ## mm1=invsqrt 

        pfmul %mm1,%mm0         ## mm0=r 
        pfmul nb332_tsc(%esp),%mm0
        pf2iw %mm0,%mm4
        movq %mm4,nb332_n1(%esp)
        pi2fd %mm4,%mm4
        pfsub %mm4,%mm0              ## now mm0 is eps and mm4 n0 
        movq  %mm0,%mm2
        pfmul %mm2,%mm2         ## mm0 is eps, mm2 eps2 

        ## coulomb table 
        movl nb332_VFtab(%ebp),%edx
        movl nb332_n1(%esp),%ecx
        leal (%ecx,%ecx,2),%ecx
        shll $2,%ecx
        ## load all values we need 
        movd (%edx,%ecx,4),%mm4
        movd 4(%edx,%ecx,4),%mm5
        movd 8(%edx,%ecx,4),%mm6
        movd 12(%edx,%ecx,4),%mm7
        movl nb332_n1+4(%esp),%ecx
        leal (%ecx,%ecx,2),%ecx
        shll $2,%ecx
        punpckldq (%edx,%ecx,4),%mm4
        punpckldq 4(%edx,%ecx,4),%mm5
        punpckldq 8(%edx,%ecx,4),%mm6
        punpckldq 12(%edx,%ecx,4),%mm7

        pfmul %mm0,%mm6 ## mm6 = Geps           
        pfmul %mm2,%mm7 ## mm7 = Heps2 

        pfadd %mm6,%mm5
        pfadd %mm7,%mm5 ## mm5 = Fp 

        pfmul nb332_two(%esp),%mm7      ## two*Heps2 
        pfadd %mm6,%mm7
        pfadd %mm5,%mm7 ## mm7=FF 

        pfmul %mm0,%mm5 ## mm5=eps*Fp 
        pfadd %mm4,%mm5 ##  mm5= VV 

        pfmul nb332_qqOH(%esp),%mm5     ## vcoul=qq*VV 
        pfmul nb332_qqOH(%esp),%mm7     ## fijC=qq*FF 
        ## update vctot 
        pfadd nb332_vctot(%esp),%mm5
        movq %mm5,nb332_vctot(%esp)

        ## change sign of fijC and multiply by rinv 
    pxor %mm4,%mm4
        pfsub %mm7,%mm4
        pfmul nb332_tsc(%esp),%mm4
        pfmul %mm1,%mm4   ## mm4 is total fscal (for the hydrogens) now         

        ## spread oxygen fscalar to both positions 
        punpckldq %mm3,%mm3
        ## calc vectorial force for O 
        movq nb332_dxO(%esp),%mm0
        movd nb332_dzO(%esp),%mm1
        pfmul %mm3,%mm0
        pfmul %mm3,%mm1

        ## calc vectorial force for H's 
        movq nb332_dxH(%esp),%mm5
        movq nb332_dyH(%esp),%mm6
        movq nb332_dzH(%esp),%mm7
        pfmul %mm4,%mm5
        pfmul %mm4,%mm6
        pfmul %mm4,%mm7

        ## update iO particle force 
        movq nb332_fixO(%esp),%mm2
        movd nb332_fizO(%esp),%mm3
        pfadd %mm0,%mm2
        pfadd %mm1,%mm3
        movq %mm2,nb332_fixO(%esp)
        movd %mm3,nb332_fizO(%esp)

        ## update iH forces 
        movq nb332_fixH(%esp),%mm2
        movq nb332_fiyH(%esp),%mm3
        movq nb332_fizH(%esp),%mm4
        pfadd %mm5,%mm2
        pfadd %mm6,%mm3
        pfadd %mm7,%mm4
        movq %mm2,nb332_fixH(%esp)
        movq %mm3,nb332_fiyH(%esp)
        movq %mm4,nb332_fizH(%esp)

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

        pfsubr nb332_ixO(%esp),%mm0
        pfsubr nb332_izO(%esp),%mm1

        movq  %mm0,nb332_dxO(%esp)
        pfmul %mm0,%mm0
        movd  %mm1,nb332_dzO(%esp)
        pfmul %mm1,%mm1
        pfacc %mm1,%mm0
        pfadd %mm1,%mm0         ## mm0=rsqO 

        punpckldq %mm2,%mm2
        punpckldq %mm3,%mm3
        punpckldq %mm4,%mm4 ## mm2-mm4 is jx-jz 
        pfsubr nb332_ixH(%esp),%mm2
        pfsubr nb332_iyH(%esp),%mm3
        pfsubr nb332_izH(%esp),%mm4   ## mm2-mm4 is dxH-dzH 

        movq %mm2,nb332_dxH(%esp)
        movq %mm3,nb332_dyH(%esp)
        movq %mm4,nb332_dzH(%esp)
        pfmul %mm2,%mm2
        pfmul %mm3,%mm3
        pfmul %mm4,%mm4

        pfadd %mm2,%mm3
        pfadd %mm4,%mm3         ## mm3=rsqH 
        movq %mm3,nb332_tmprsqH(%esp)

    pfrsqrt %mm0,%mm1

    movq %mm1,%mm2
    pfmul %mm1,%mm1
    pfrsqit1 %mm0,%mm1
    pfrcpit2 %mm2,%mm1  ## mm1=invsqrt 
        pfmul %mm1,%mm0         ## mm0=rsq  

        pfmul nb332_tsc(%esp),%mm0
        pf2iw %mm0,%mm4
        movd %mm4,nb332_n1(%esp)
        pi2fd %mm4,%mm4
        pfsub %mm4,%mm0              ## now mm0 is eps and mm4 n0 
        movq  %mm0,%mm2
        pfmul %mm2,%mm2         ## mm0 is eps, mm2 eps2 

        ## coulomb table 
        movl nb332_VFtab(%ebp),%edx
        movl nb332_n1(%esp),%ecx
        leal (%ecx,%ecx,2),%ecx
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

        pfmul nb332_two(%esp),%mm7      ## two*Heps2 
        pfadd %mm6,%mm7
        pfadd %mm5,%mm7 ## mm7=FF 

        pfmul %mm0,%mm5 ## mm5=eps*Fp 
        pfadd %mm4,%mm5 ##  mm5= VV 

        pfmul nb332_qqOH(%esp),%mm5     ## vcoul=qq*VV 
        pfmul nb332_qqOH(%esp),%mm7     ## fijC=qq*FF 

        ## update vctot directly, force is moved to mm3. 
        pfadd nb332_vctot(%esp),%mm5
        movq %mm5,nb332_vctot(%esp)
        pxor %mm3,%mm3
        pfsub %mm7,%mm3
        pfmul nb332_tsc(%esp),%mm3
        pfmul %mm1,%mm3   ## mm3 is total fscal (for the oxygen) now 

        movq nb332_tmprsqH(%esp),%mm0

        pfrsqrt %mm0,%mm1
        pswapd %mm0,%mm0
        pfrsqrt %mm0,%mm2
        pswapd %mm0,%mm0
        punpckldq %mm2,%mm1     ## seeds are in mm1 now, and rsq in mm0. 

        movq %mm1,%mm2
        pfmul %mm1,%mm1
    pfrsqit1 %mm0,%mm1
    pfrcpit2 %mm2,%mm1  ## mm1=invsqrt 

        pfmul %mm1,%mm0         ## mm0=r 
        pfmul nb332_tsc(%esp),%mm0
        pf2iw %mm0,%mm4
        movq %mm4,nb332_n1(%esp)
        pi2fd %mm4,%mm4
        pfsub %mm4,%mm0              ## now mm0 is eps and mm4 n0 
        movq  %mm0,%mm2
        pfmul %mm2,%mm2         ## mm0 is eps, mm2 eps2 

        ## coulomb table 
        movl nb332_VFtab(%ebp),%edx
        movl nb332_n1(%esp),%ecx
        leal (%ecx,%ecx,2),%ecx
        shll $2,%ecx
        ## load all values we need 
        movd (%edx,%ecx,4),%mm4
        movd 4(%edx,%ecx,4),%mm5
        movd 8(%edx,%ecx,4),%mm6
        movd 12(%edx,%ecx,4),%mm7
        movl nb332_n1+4(%esp),%ecx
        leal (%ecx,%ecx,2),%ecx
        shll $2,%ecx
        punpckldq (%edx,%ecx,4),%mm4
        punpckldq 4(%edx,%ecx,4),%mm5
        punpckldq 8(%edx,%ecx,4),%mm6
        punpckldq 12(%edx,%ecx,4),%mm7


        pfmul %mm0,%mm6 ## mm6 = Geps           
        pfmul %mm2,%mm7 ## mm7 = Heps2 

        pfadd %mm6,%mm5
        pfadd %mm7,%mm5 ## mm5 = Fp 

        pfmul nb332_two(%esp),%mm7      ## two*Heps2 
        pfadd %mm6,%mm7
        pfadd %mm5,%mm7 ## mm7=FF 

        pfmul %mm0,%mm5 ## mm5=eps*Fp 
        pfadd %mm4,%mm5 ##  mm5= VV 

        pfmul nb332_qqHH(%esp),%mm5     ## vcoul=qq*VV 
        pfmul nb332_qqHH(%esp),%mm7     ## fijC=qq*FF 
        ## update vctot 
        pfadd nb332_vctot(%esp),%mm5
        movq %mm5,nb332_vctot(%esp)

        ## change sign of fijC and multiply by rinv 
    pxor %mm4,%mm4
        pfsub %mm7,%mm4
        pfmul nb332_tsc(%esp),%mm4
        pfmul %mm1,%mm4   ## mm4 is total fscal (for the hydrogens) now                 

        ## spread oxygen fscalar to both positions 
        punpckldq %mm3,%mm3
        ## calc vectorial force for O 
        movq nb332_dxO(%esp),%mm0
        movd nb332_dzO(%esp),%mm1
        pfmul %mm3,%mm0
        pfmul %mm3,%mm1

        ## calc vectorial force for H's 
        movq nb332_dxH(%esp),%mm5
        movq nb332_dyH(%esp),%mm6
        movq nb332_dzH(%esp),%mm7
        pfmul %mm4,%mm5
        pfmul %mm4,%mm6
        pfmul %mm4,%mm7

        ## update iO particle force 
        movq nb332_fixO(%esp),%mm2
        movd nb332_fizO(%esp),%mm3
        pfadd %mm0,%mm2
        pfadd %mm1,%mm3
        movq %mm2,nb332_fixO(%esp)
        movd %mm3,nb332_fizO(%esp)

        ## update iH forces 
        movq nb332_fixH(%esp),%mm2
        movq nb332_fiyH(%esp),%mm3
        movq nb332_fizH(%esp),%mm4
        pfadd %mm5,%mm2
        pfadd %mm6,%mm3
        pfadd %mm7,%mm4
        movq %mm2,nb332_fixH(%esp)
        movq %mm3,nb332_fiyH(%esp)
        movq %mm4,nb332_fizH(%esp)

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

        pfsubr nb332_ixO(%esp),%mm0
        pfsubr nb332_izO(%esp),%mm1

        movq  %mm0,nb332_dxO(%esp)
        pfmul %mm0,%mm0
        movd  %mm1,nb332_dzO(%esp)
        pfmul %mm1,%mm1
        pfacc %mm1,%mm0
        pfadd %mm1,%mm0         ## mm0=rsqO 

        punpckldq %mm2,%mm2
        punpckldq %mm3,%mm3
        punpckldq %mm4,%mm4 ## mm2-mm4 is jx-jz 
        pfsubr nb332_ixH(%esp),%mm2
        pfsubr nb332_iyH(%esp),%mm3
        pfsubr nb332_izH(%esp),%mm4   ## mm2-mm4 is dxH-dzH 

        movq %mm2,nb332_dxH(%esp)
        movq %mm3,nb332_dyH(%esp)
        movq %mm4,nb332_dzH(%esp)
        pfmul %mm2,%mm2
        pfmul %mm3,%mm3
        pfmul %mm4,%mm4

        pfadd %mm2,%mm3
        pfadd %mm4,%mm3         ## mm3=rsqH 
        movq %mm3,nb332_tmprsqH(%esp)

    pfrsqrt %mm0,%mm1

    movq %mm1,%mm2
    pfmul %mm1,%mm1
    pfrsqit1 %mm0,%mm1
    pfrcpit2 %mm2,%mm1  ## mm1=invsqrt 
        pfmul %mm1,%mm0

        pfmul nb332_tsc(%esp),%mm0
        pf2iw %mm0,%mm4
        movd %mm4,nb332_n1(%esp)
        pi2fd %mm4,%mm4
        pfsub %mm4,%mm0              ## now mm0 is eps and mm4 n0 
        movq  %mm0,%mm2
        pfmul %mm2,%mm2         ## mm0 is eps, mm2 eps2 

        ## coulomb table 
        movl nb332_VFtab(%ebp),%edx
        movl nb332_n1(%esp),%ecx
        leal (%ecx,%ecx,2),%ecx
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

        pfmul nb332_two(%esp),%mm7      ## two*Heps2 
        pfadd %mm6,%mm7
        pfadd %mm5,%mm7 ## mm7=FF 

        pfmul %mm0,%mm5 ## mm5=eps*Fp 
        pfadd %mm4,%mm5 ##  mm5= VV 

        pfmul nb332_qqOH(%esp),%mm5     ## vcoul=qq*VV 
        pfmul nb332_qqOH(%esp),%mm7     ## fijC=qq*FF 

        ## update vctot directly, use mm3 for fscal sum 
        pfadd nb332_vctot(%esp),%mm5
        movq %mm5,nb332_vctot(%esp)
        pxor %mm3,%mm3
        pfsub %mm7,%mm3
        pfmul nb332_tsc(%esp),%mm3
        pfmul %mm1,%mm3    ## mm3 is total fscal (for the oxygen) now   

        movq nb332_tmprsqH(%esp),%mm0

        pfrsqrt %mm0,%mm1
        pswapd %mm0,%mm0
        pfrsqrt %mm0,%mm2
        pswapd %mm0,%mm0
        punpckldq %mm2,%mm1     ## seeds are in mm1 now, and rsq in mm0. 

        movq %mm1,%mm2
        pfmul %mm1,%mm1
    pfrsqit1 %mm0,%mm1
    pfrcpit2 %mm2,%mm1  ## mm1=invsqrt 

        pfmul %mm1,%mm0         ## mm0=r 
        pfmul nb332_tsc(%esp),%mm0
        pf2iw %mm0,%mm4
        movq %mm4,nb332_n1(%esp)
        pi2fd %mm4,%mm4
        pfsub %mm4,%mm0              ## now mm0 is eps and mm4 n0 
        movq  %mm0,%mm2
        pfmul %mm2,%mm2         ## mm0 is eps, mm2 eps2 

        ## coulomb table 
        movl nb332_VFtab(%ebp),%edx
        movl nb332_n1(%esp),%ecx
        leal (%ecx,%ecx,2),%ecx
        shll $2,%ecx
        ## load all values we need 
        movd (%edx,%ecx,4),%mm4
        movd 4(%edx,%ecx,4),%mm5
        movd 8(%edx,%ecx,4),%mm6
        movd 12(%edx,%ecx,4),%mm7
        movl nb332_n1+4(%esp),%ecx ## mm5 = Fp 
        leal (%ecx,%ecx,2),%ecx
        shll $2,%ecx
        punpckldq (%edx,%ecx,4),%mm4
        punpckldq 4(%edx,%ecx,4),%mm5
        punpckldq 8(%edx,%ecx,4),%mm6
        punpckldq 12(%edx,%ecx,4),%mm7


        pfmul %mm0,%mm6 ## mm6 = Geps           
        pfmul %mm2,%mm7 ## mm7 = Heps2 

        pfadd %mm6,%mm5
        pfadd %mm7,%mm5 ## mm5 = Fp 

        pfmul nb332_two(%esp),%mm7      ## two*Heps2 
        pfadd %mm6,%mm7
        pfadd %mm5,%mm7 ## mm7=FF 

        pfmul %mm0,%mm5 ## mm5=eps*Fp 
        pfadd %mm4,%mm5 ##  mm5= VV 

        pfmul nb332_qqHH(%esp),%mm5     ## vcoul=qq*VV 
        pfmul nb332_qqHH(%esp),%mm7     ## fijC=qq*FF 
        ## update vctot 
        pfadd nb332_vctot(%esp),%mm5
        movq %mm5,nb332_vctot(%esp)

        ## change sign of fijC and multiply by rinv 
    pxor %mm4,%mm4
        pfsub %mm7,%mm4
        pfmul nb332_tsc(%esp),%mm4
        pfmul %mm1,%mm4   ## mm4 is total fscal (for the hydrogens) now         

        ## spread oxygen fscalar to both positions 
        punpckldq %mm3,%mm3
        ## calc vectorial force for O 
        movq nb332_dxO(%esp),%mm0
        movd nb332_dzO(%esp),%mm1
        pfmul %mm3,%mm0
        pfmul %mm3,%mm1

        ## calc vectorial force for H's 
        movq nb332_dxH(%esp),%mm5
        movq nb332_dyH(%esp),%mm6
        movq nb332_dzH(%esp),%mm7
        pfmul %mm4,%mm5
        pfmul %mm4,%mm6
        pfmul %mm4,%mm7

        ## update iO particle force 
        movq nb332_fixO(%esp),%mm2
        movd nb332_fizO(%esp),%mm3
        pfadd %mm0,%mm2
        pfadd %mm1,%mm3
        movq %mm2,nb332_fixO(%esp)
        movd %mm3,nb332_fizO(%esp)

        ## update iH forces 
        movq nb332_fixH(%esp),%mm2
        movq nb332_fiyH(%esp),%mm3
        movq nb332_fizH(%esp),%mm4
        pfadd %mm5,%mm2
        pfadd %mm6,%mm3
        pfadd %mm7,%mm4
        movq %mm2,nb332_fixH(%esp)
        movq %mm3,nb332_fiyH(%esp)
        movq %mm4,nb332_fizH(%esp)

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
        decl nb332_innerk(%esp)
        jz  _nb_kernel332_ia32_3dnow.nb332_updateouterdata
        jmp _nb_kernel332_ia32_3dnow.nb332_inner_loop
_nb_kernel332_ia32_3dnow.nb332_updateouterdata: 
        movl  nb332_ii3(%esp),%ecx

        movq  (%edi,%ecx,4),%mm6       ## increment iO force  
        movd  8(%edi,%ecx,4),%mm7
        pfadd nb332_fixO(%esp),%mm6
        pfadd nb332_fizO(%esp),%mm7
        movq  %mm6,(%edi,%ecx,4)
        movd  %mm7,8(%edi,%ecx,4)

        movq  nb332_fixH(%esp),%mm0
        movq  nb332_fiyH(%esp),%mm3
        movq  nb332_fizH(%esp),%mm1
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


        movl  nb332_fshift(%ebp),%ebx      ## increment fshift force 
        movl  nb332_is3(%esp),%edx

        movq  (%ebx,%edx,4),%mm6
        movd  8(%ebx,%edx,4),%mm7
        pfadd nb332_fixO(%esp),%mm6
        pfadd nb332_fizO(%esp),%mm7
        pfadd %mm0,%mm6
        pfadd %mm1,%mm7
        pfadd %mm2,%mm6
        pfadd %mm3,%mm7
        movq  %mm6,(%ebx,%edx,4)
        movd  %mm7,8(%ebx,%edx,4)

        ## get n from stack
        movl nb332_n(%esp),%esi
        ## get group index for i particle 
        movl  nb332_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        movq  nb332_vctot(%esp),%mm7
        pfacc %mm7,%mm7           ## get and sum the two parts of total potential 

        movl  nb332_Vc(%ebp),%eax
        movd  (%eax,%edx,4),%mm6
        pfadd %mm7,%mm6
        movd  %mm6,(%eax,%edx,4)          ## increment vc[gid] 

        movq  nb332_Vvdwtot(%esp),%mm7
        pfacc %mm7,%mm7           ## get and sum the two parts of total potential 

        movl  nb332_Vvdw(%ebp),%eax
        movd  (%eax,%edx,4),%mm6
        pfadd %mm7,%mm6
        movd  %mm6,(%eax,%edx,4)          ## increment Vvdwtot[gid] 
        ## finish if last 
        movl nb332_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jecxz _nb_kernel332_ia32_3dnow.nb332_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb332_n(%esp)
        jmp _nb_kernel332_ia32_3dnow.nb332_outer
_nb_kernel332_ia32_3dnow.nb332_outerend: 
        ## check if more outer neighborlists remain
        movl  nb332_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jecxz _nb_kernel332_ia32_3dnow.nb332_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel332_ia32_3dnow.nb332_threadloop
_nb_kernel332_ia32_3dnow.nb332_end: 
        femms
        movl nb332_nouter(%esp),%eax
        movl nb332_ninner(%esp),%ebx
        movl nb332_outeriter(%ebp),%ecx
        movl nb332_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        addl $232,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret



.globl nb_kernel332nf_ia32_3dnow
.globl _nb_kernel332nf_ia32_3dnow
nb_kernel332nf_ia32_3dnow:      
_nb_kernel332nf_ia32_3dnow:     
.set nb332nf_p_nri, 8
.set nb332nf_iinr, 12
.set nb332nf_jindex, 16
.set nb332nf_jjnr, 20
.set nb332nf_shift, 24
.set nb332nf_shiftvec, 28
.set nb332nf_fshift, 32
.set nb332nf_gid, 36
.set nb332nf_pos, 40
.set nb332nf_faction, 44
.set nb332nf_charge, 48
.set nb332nf_p_facel, 52
.set nb332nf_p_krf, 56
.set nb332nf_p_crf, 60
.set nb332nf_Vc, 64
.set nb332nf_type, 68
.set nb332nf_p_ntype, 72
.set nb332nf_vdwparam, 76
.set nb332nf_Vvdw, 80
.set nb332nf_p_tabscale, 84
.set nb332nf_VFtab, 88
.set nb332nf_invsqrta, 92
.set nb332nf_dvda, 96
.set nb332nf_p_gbtabscale, 100
.set nb332nf_GBtab, 104
.set nb332nf_p_nthreads, 108
.set nb332nf_count, 112
.set nb332nf_mtx, 116
.set nb332nf_outeriter, 120
.set nb332nf_inneriter, 124
.set nb332nf_work, 128
                        ## stack offsets for local variables 
.set nb332nf_is3, 0
.set nb332nf_ii3, 4
.set nb332nf_ixO, 8
.set nb332nf_iyO, 12
.set nb332nf_izO, 16
.set nb332nf_ixH, 20
.set nb332nf_iyH, 28
.set nb332nf_izH, 36
.set nb332nf_qqOO, 44
.set nb332nf_qqOH, 52
.set nb332nf_qqHH, 60
.set nb332nf_c6, 68
.set nb332nf_c12, 76
.set nb332nf_n1, 84
.set nb332nf_tsc, 92
.set nb332nf_vctot, 100
.set nb332nf_Vvdwtot, 108
.set nb332nf_innerjjnr, 116
.set nb332nf_innerk, 120
.set nb332nf_tmprsqH, 124
.set nb332nf_n, 132                         ## idx for outer loop
.set nb332nf_nn1, 136                       ## number of outer iterations
.set nb332nf_nri, 140
.set nb332nf_nouter, 144
.set nb332nf_ninner, 148
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $152,%esp          ## local stack space 
        femms
        movl nb332nf_p_nri(%ebp),%ecx
        movl nb332nf_p_facel(%ebp),%esi
        movl nb332nf_p_tabscale(%ebp),%edi
        movl (%ecx),%ecx
        movl %ecx,nb332nf_nri(%esp)
        movd  (%esi),%mm1       ## facel
        movl nb332nf_p_ntype(%ebp),%esi

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb332nf_nouter(%esp)
        movl %eax,nb332nf_ninner(%esp)

        ## assume we have at least one i particle - start directly      

        movl  nb332nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx),%ebx           ## ebx=ii 

        movl  nb332nf_charge(%ebp),%edx
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
        movq  %mm4,nb332nf_qqOO(%esp)
        movq  %mm5,nb332nf_qqOH(%esp)
        movq  %mm6,nb332nf_qqHH(%esp)
        movl  nb332nf_type(%ebp),%edx
        movl  (%edx,%ebx,4),%ecx
        shll  %ecx
        movl  %ecx,%edx
        imull (%esi),%ecx
        addl  %ecx,%edx
        movl  nb332nf_vdwparam(%ebp),%eax
        movd  (%eax,%edx,4),%mm0
        movd  4(%eax,%edx,4),%mm1
        movq  %mm0,nb332nf_c6(%esp)
        movq  %mm1,nb332nf_c12(%esp)
        movd  (%edi),%mm3
        punpckldq %mm3,%mm3
        movq  %mm3,nb332nf_tsc(%esp)
_nb_kernel332nf_ia32_3dnow.nb332nf_threadloop: 
        movl  nb332nf_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel332nf_ia32_3dnow.nb332nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel332nf_ia32_3dnow.nb332nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb332nf_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb332nf_n(%esp)
        movl %ebx,nb332nf_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel332nf_ia32_3dnow.nb332nf_outerstart
        jmp _nb_kernel332nf_ia32_3dnow.nb332nf_end

_nb_kernel332nf_ia32_3dnow.nb332nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb332nf_nouter(%esp),%ebx
        movl %ebx,nb332nf_nouter(%esp)

_nb_kernel332nf_ia32_3dnow.nb332nf_outer: 
        movl  nb332nf_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb332nf_is3(%esp)            ## store is3 

        movl  nb332nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movq  (%eax,%ebx,4),%mm5        ## move shX/shY to mm5 and shZ to mm6. 
        movd  8(%eax,%ebx,4),%mm6
        movq  %mm5,%mm0
        movq  %mm5,%mm1
        movq  %mm6,%mm2
        punpckldq %mm0,%mm0         ## also expand shX,Y,Z in mm0--mm2. 
        punpckhdq %mm1,%mm1
        punpckldq %mm2,%mm2

        movl  nb332nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx,%esi,4),%ebx    ## ebx=ii 

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb332nf_pos(%ebp),%eax      ## eax = base of pos[] 

        pfadd (%eax,%ebx,4),%mm5    ## ix = shX + posX (and iy too) 
        movd  8(%eax,%ebx,4),%mm7       ## cant use direct memory add for 4 bytes (iz) 
        movl  %ebx,nb332nf_ii3(%esp)        ## (use mm7 as temp. storage for iz.) 
        pfadd %mm7,%mm6
        movq  %mm5,nb332nf_ixO(%esp)
        movq  %mm6,nb332nf_izO(%esp)

        movd  12(%eax,%ebx,4),%mm3
        movd  16(%eax,%ebx,4),%mm4
        movd  20(%eax,%ebx,4),%mm5
        punpckldq  24(%eax,%ebx,4),%mm3
        punpckldq  28(%eax,%ebx,4),%mm4
        punpckldq  32(%eax,%ebx,4),%mm5    ## coords of H1 in low mm3-mm5, H2 in high 

        pfadd %mm3,%mm0
        pfadd %mm4,%mm1
        pfadd %mm5,%mm2
        movq %mm0,nb332nf_ixH(%esp)
        movq %mm1,nb332nf_iyH(%esp)
        movq %mm2,nb332nf_izH(%esp)

        ## clear vctot and i forces 
        pxor  %mm7,%mm7
        movq  %mm7,nb332nf_vctot(%esp)
        movq  %mm7,nb332nf_Vvdwtot(%esp)

        movl  nb332nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx     ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 
        movl  %edx,nb332nf_innerk(%esp)
        addl  nb332nf_ninner(%esp),%edx
        movl  %edx,nb332nf_ninner(%esp)

        movl  nb332nf_pos(%ebp),%esi
        movl  nb332nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb332nf_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
_nb_kernel332nf_ia32_3dnow.nb332nf_inner_loop: 
        ## a single j particle iteration here - compare with the unrolled code for comments. 
        movl  nb332nf_innerjjnr(%esp),%eax
        movl  (%eax),%eax       ## eax=jnr offset 
    addl $4,nb332nf_innerjjnr(%esp)             ## advance pointer 

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

        pfsubr nb332nf_ixO(%esp),%mm0
        pfsubr nb332nf_izO(%esp),%mm1

        pfmul %mm0,%mm0
        pfmul %mm1,%mm1
        pfacc %mm0,%mm0
        pfadd %mm1,%mm0         ## mm0=rsqO 

        punpckldq %mm2,%mm2
        punpckldq %mm3,%mm3
        punpckldq %mm4,%mm4 ## mm2-mm4 is jx-jz 
        pfsubr nb332nf_ixH(%esp),%mm2
        pfsubr nb332nf_iyH(%esp),%mm3
        pfsubr nb332nf_izH(%esp),%mm4   ## mm2-mm4 is dxH-dzH 

        pfmul %mm2,%mm2
        pfmul %mm3,%mm3
        pfmul %mm4,%mm4

        pfadd %mm2,%mm3
        pfadd %mm4,%mm3         ## mm3=rsqH 
        movq %mm3,nb332nf_tmprsqH(%esp)

    pfrsqrt %mm0,%mm1

    movq %mm1,%mm2
    pfmul %mm1,%mm1
    pfrsqit1 %mm0,%mm1
    pfrcpit2 %mm2,%mm1  ## mm1=invsqrt  
        pfmul %mm1,%mm0         ## mm0=rsq  

        pfmul nb332nf_tsc(%esp),%mm0
        pf2iw %mm0,%mm4
        movd %mm4,nb332nf_n1(%esp)
        pi2fd %mm4,%mm4
        pfsub %mm4,%mm0              ## now mm0 is eps and mm4 n0 
        movq  %mm0,%mm2
        pfmul %mm2,%mm2         ## mm0 is eps, mm2 eps2 

        ## coulomb table 
        movl nb332nf_VFtab(%ebp),%edx
        movl nb332nf_n1(%esp),%ecx
        leal (%ecx,%ecx,2),%ecx
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

        pfmul nb332nf_qqOO(%esp),%mm5   ## vcoul=qq*VV 
        ## update vctot directly, use mm3 for fscal sum. 
        pfadd nb332nf_vctot(%esp),%mm5
        movq %mm5,nb332nf_vctot(%esp)

        ## dispersion table 
        ## load all the table values we need 
        movd 16(%edx,%ecx,4),%mm4
        movd 20(%edx,%ecx,4),%mm5
        movd 24(%edx,%ecx,4),%mm6
        movd 28(%edx,%ecx,4),%mm7
        pfmul %mm0,%mm6 ## mm6 = Geps           
        pfmul %mm2,%mm7 ## mm7 = Heps2 
        pfadd %mm6,%mm5
        pfadd %mm7,%mm5 ## mm5 = Fp 
        pfmul %mm0,%mm5 ## mm5=eps*Fp 
        pfadd %mm4,%mm5 ##  mm5= VV     

        movq nb332nf_c6(%esp),%mm4
        pfmul %mm4,%mm5 ## Vvdw6            
        ## update Vvdwtot to release mm5! 
        pfadd nb332nf_Vvdwtot(%esp),%mm5        ## add the earlier value 
        movq %mm5,nb332nf_Vvdwtot(%esp)         ## store the sum       

        ## repulsion table 
        ## load all the table values we need 
        movd 32(%edx,%ecx,4),%mm4
        movd 36(%edx,%ecx,4),%mm5
        movd 40(%edx,%ecx,4),%mm6
        movd 44(%edx,%ecx,4),%mm7

        pfmul %mm0,%mm6 ## mm6 = Geps           
        pfmul %mm2,%mm7 ## mm7 = Heps2 
        pfadd %mm6,%mm5
        pfadd %mm7,%mm5 ## mm5 = Fp 
        pfmul %mm0,%mm5 ## mm5=eps*Fp 
        pfadd %mm4,%mm5 ##  mm5= VV 

        movq nb332nf_c12(%esp),%mm6
        pfmul %mm6,%mm5 ## Vvdw12 
        ## change sign of fscal and multiply with rinv  
        ## update Vvdwtot  
        pfadd nb332nf_Vvdwtot(%esp),%mm5        ## add the earlier value 
        movq %mm5,nb332nf_Vvdwtot(%esp)         ## store the sum       

        ## Ready with the oxygen - time for hydrogens 

        movq nb332nf_tmprsqH(%esp),%mm0

        pfrsqrt %mm0,%mm1
        pswapd %mm0,%mm0
        pfrsqrt %mm0,%mm2
        pswapd %mm0,%mm0
        punpckldq %mm2,%mm1     ## seeds are in mm1 now, and rsq in mm0. 

        movq %mm1,%mm2
        pfmul %mm1,%mm1
    pfrsqit1 %mm0,%mm1
    pfrcpit2 %mm2,%mm1  ## mm1=invsqrt 

        pfmul %mm1,%mm0         ## mm0=r 
        pfmul nb332nf_tsc(%esp),%mm0
        pf2iw %mm0,%mm4
        movq %mm4,nb332nf_n1(%esp)
        pi2fd %mm4,%mm4
        pfsub %mm4,%mm0              ## now mm0 is eps and mm4 n0 
        movq  %mm0,%mm2
        pfmul %mm2,%mm2         ## mm0 is eps, mm2 eps2 

        ## coulomb table 
        movl nb332nf_VFtab(%ebp),%edx
        movl nb332nf_n1(%esp),%ecx
        leal (%ecx,%ecx,2),%ecx
        shll $2,%ecx
        ## load all values we need 
        movd (%edx,%ecx,4),%mm4
        movd 4(%edx,%ecx,4),%mm5
        movd 8(%edx,%ecx,4),%mm6
        movd 12(%edx,%ecx,4),%mm7
        movl nb332nf_n1+4(%esp),%ecx
        leal (%ecx,%ecx,2),%ecx
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

        pfmul nb332nf_qqOH(%esp),%mm5   ## vcoul=qq*VV 
        ## update vctot 
        pfadd nb332nf_vctot(%esp),%mm5
        movq %mm5,nb332nf_vctot(%esp)

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

        pfsubr nb332nf_ixO(%esp),%mm0
        pfsubr nb332nf_izO(%esp),%mm1

        pfmul %mm0,%mm0
        pfmul %mm1,%mm1
        pfacc %mm1,%mm0
        pfadd %mm1,%mm0         ## mm0=rsqO 

        punpckldq %mm2,%mm2
        punpckldq %mm3,%mm3
        punpckldq %mm4,%mm4 ## mm2-mm4 is jx-jz 
        pfsubr nb332nf_ixH(%esp),%mm2
        pfsubr nb332nf_iyH(%esp),%mm3
        pfsubr nb332nf_izH(%esp),%mm4   ## mm2-mm4 is dxH-dzH 

        pfmul %mm2,%mm2
        pfmul %mm3,%mm3
        pfmul %mm4,%mm4

        pfadd %mm2,%mm3
        pfadd %mm4,%mm3         ## mm3=rsqH 
        movq %mm3,nb332nf_tmprsqH(%esp)

    pfrsqrt %mm0,%mm1

    movq %mm1,%mm2
    pfmul %mm1,%mm1
    pfrsqit1 %mm0,%mm1
    pfrcpit2 %mm2,%mm1  ## mm1=invsqrt 
        pfmul %mm1,%mm0         ## mm0=rsq  

        pfmul nb332nf_tsc(%esp),%mm0
        pf2iw %mm0,%mm4
        movd %mm4,nb332nf_n1(%esp)
        pi2fd %mm4,%mm4
        pfsub %mm4,%mm0              ## now mm0 is eps and mm4 n0 
        movq  %mm0,%mm2
        pfmul %mm2,%mm2         ## mm0 is eps, mm2 eps2 

        ## coulomb table 
        movl nb332nf_VFtab(%ebp),%edx
        movl nb332nf_n1(%esp),%ecx
        leal (%ecx,%ecx,2),%ecx
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

        pfmul nb332nf_qqOH(%esp),%mm5   ## vcoul=qq*VV 
        ## update vctot directly 
        pfadd nb332nf_vctot(%esp),%mm5
        movq %mm5,nb332nf_vctot(%esp)

        movq nb332nf_tmprsqH(%esp),%mm0
        pfrsqrt %mm0,%mm1
        pswapd %mm0,%mm0
        pfrsqrt %mm0,%mm2
        pswapd %mm0,%mm0
        punpckldq %mm2,%mm1     ## seeds are in mm1 now, and rsq in mm0. 

        movq %mm1,%mm2
        pfmul %mm1,%mm1
    pfrsqit1 %mm0,%mm1
    pfrcpit2 %mm2,%mm1  ## mm1=invsqrt 

        pfmul %mm1,%mm0         ## mm0=r 
        pfmul nb332nf_tsc(%esp),%mm0
        pf2iw %mm0,%mm4
        movq %mm4,nb332nf_n1(%esp)
        pi2fd %mm4,%mm4
        pfsub %mm4,%mm0              ## now mm0 is eps and mm4 n0 
        movq  %mm0,%mm2
        pfmul %mm2,%mm2         ## mm0 is eps, mm2 eps2 

        ## coulomb table 
        movl nb332nf_VFtab(%ebp),%edx
        movl nb332nf_n1(%esp),%ecx
        leal (%ecx,%ecx,2),%ecx
        shll $2,%ecx
        ## load all values we need 
        movd (%edx,%ecx,4),%mm4
        movd 4(%edx,%ecx,4),%mm5
        movd 8(%edx,%ecx,4),%mm6
        movd 12(%edx,%ecx,4),%mm7
        movl nb332nf_n1+4(%esp),%ecx
        leal (%ecx,%ecx,2),%ecx
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

        pfmul nb332nf_qqHH(%esp),%mm5   ## vcoul=qq*VV 
        ## update vctot 
        pfadd nb332nf_vctot(%esp),%mm5
        movq %mm5,nb332nf_vctot(%esp)

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

        pfsubr nb332nf_ixO(%esp),%mm0
        pfsubr nb332nf_izO(%esp),%mm1

        pfmul %mm0,%mm0
        pfmul %mm1,%mm1
        pfacc %mm1,%mm0
        pfadd %mm1,%mm0         ## mm0=rsqO 

        punpckldq %mm2,%mm2
        punpckldq %mm3,%mm3
        punpckldq %mm4,%mm4 ## mm2-mm4 is jx-jz 
        pfsubr nb332nf_ixH(%esp),%mm2
        pfsubr nb332nf_iyH(%esp),%mm3
        pfsubr nb332nf_izH(%esp),%mm4   ## mm2-mm4 is dxH-dzH 

        pfmul %mm2,%mm2
        pfmul %mm3,%mm3
        pfmul %mm4,%mm4

        pfadd %mm2,%mm3
        pfadd %mm4,%mm3         ## mm3=rsqH 
        movq %mm3,nb332nf_tmprsqH(%esp)

    pfrsqrt %mm0,%mm1

    movq %mm1,%mm2
    pfmul %mm1,%mm1
    pfrsqit1 %mm0,%mm1
    pfrcpit2 %mm2,%mm1  ## mm1=invsqrt 
        pfmul %mm1,%mm0

        pfmul nb332nf_tsc(%esp),%mm0
        pf2iw %mm0,%mm4
        movd %mm4,nb332nf_n1(%esp)
        pi2fd %mm4,%mm4
        pfsub %mm4,%mm0              ## now mm0 is eps and mm4 n0 
        movq  %mm0,%mm2
        pfmul %mm2,%mm2         ## mm0 is eps, mm2 eps2 

        ## coulomb table 
        movl nb332nf_VFtab(%ebp),%edx
        movl nb332nf_n1(%esp),%ecx
        leal (%ecx,%ecx,2),%ecx
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

        pfmul nb332nf_qqOH(%esp),%mm5   ## vcoul=qq*VV 
        ## update vctot directly 
        pfadd nb332nf_vctot(%esp),%mm5
        movq %mm5,nb332nf_vctot(%esp)

        movq nb332nf_tmprsqH(%esp),%mm0
        pfrsqrt %mm0,%mm1
        pswapd %mm0,%mm0
        pfrsqrt %mm0,%mm2
        pswapd %mm0,%mm0
        punpckldq %mm2,%mm1     ## seeds are in mm1 now, and rsq in mm0. 

        movq %mm1,%mm2
        pfmul %mm1,%mm1
    pfrsqit1 %mm0,%mm1
    pfrcpit2 %mm2,%mm1  ## mm1=invsqrt 

        pfmul %mm1,%mm0         ## mm0=r 
        pfmul nb332nf_tsc(%esp),%mm0
        pf2iw %mm0,%mm4
        movq %mm4,nb332nf_n1(%esp)
        pi2fd %mm4,%mm4
        pfsub %mm4,%mm0              ## now mm0 is eps and mm4 n0 
        movq  %mm0,%mm2
        pfmul %mm2,%mm2         ## mm0 is eps, mm2 eps2 

        ## coulomb table 
        movl nb332nf_VFtab(%ebp),%edx
        movl nb332nf_n1(%esp),%ecx
        leal (%ecx,%ecx,2),%ecx
        shll $2,%ecx
        ## load all values we need 
        movd (%edx,%ecx,4),%mm4
        movd 4(%edx,%ecx,4),%mm5
        movd 8(%edx,%ecx,4),%mm6
        movd 12(%edx,%ecx,4),%mm7
        movl nb332nf_n1+4(%esp),%ecx  ## mm5 = Fp 
        leal (%ecx,%ecx,2),%ecx
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

        pfmul nb332nf_qqHH(%esp),%mm5   ## vcoul=qq*VV 
        ## update vctot 
        pfadd nb332nf_vctot(%esp),%mm5
        movq %mm5,nb332nf_vctot(%esp)

        ##  done  - one more? 
        decl nb332nf_innerk(%esp)
        jz  _nb_kernel332nf_ia32_3dnow.nb332nf_updateouterdata
        jmp _nb_kernel332nf_ia32_3dnow.nb332nf_inner_loop
_nb_kernel332nf_ia32_3dnow.nb332nf_updateouterdata: 
        ## get n from stack
        movl nb332nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb332nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        movq  nb332nf_vctot(%esp),%mm7
        pfacc %mm7,%mm7           ## get and sum the two parts of total potential 

        movl  nb332nf_Vc(%ebp),%eax
        movd  (%eax,%edx,4),%mm6
        pfadd %mm7,%mm6
        movd  %mm6,(%eax,%edx,4)          ## increment vc[gid] 

        movq  nb332nf_Vvdwtot(%esp),%mm7
        pfacc %mm7,%mm7           ## get and sum the two parts of total potential 

        movl  nb332nf_Vvdw(%ebp),%eax
        movd  (%eax,%edx,4),%mm6
        pfadd %mm7,%mm6
        movd  %mm6,(%eax,%edx,4)          ## increment Vvdwtot[gid] 
        ## finish if last 
        movl nb332nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jecxz _nb_kernel332nf_ia32_3dnow.nb332nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb332nf_n(%esp)
        jmp _nb_kernel332nf_ia32_3dnow.nb332nf_outer
_nb_kernel332nf_ia32_3dnow.nb332nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb332nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jecxz _nb_kernel332nf_ia32_3dnow.nb332nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel332nf_ia32_3dnow.nb332nf_threadloop
_nb_kernel332nf_ia32_3dnow.nb332nf_end: 
        femms
        movl nb332nf_nouter(%esp),%eax
        movl nb332nf_ninner(%esp),%ebx
        movl nb332nf_outeriter(%ebp),%ecx
        movl nb332nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        addl $152,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret



