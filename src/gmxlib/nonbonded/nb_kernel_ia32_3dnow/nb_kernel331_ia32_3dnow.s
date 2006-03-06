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




.globl nb_kernel331_ia32_3dnow
.globl _nb_kernel331_ia32_3dnow
nb_kernel331_ia32_3dnow:        
_nb_kernel331_ia32_3dnow:       
.set nb331_p_nri, 8
.set nb331_iinr, 12
.set nb331_jindex, 16
.set nb331_jjnr, 20
.set nb331_shift, 24
.set nb331_shiftvec, 28
.set nb331_fshift, 32
.set nb331_gid, 36
.set nb331_pos, 40
.set nb331_faction, 44
.set nb331_charge, 48
.set nb331_p_facel, 52
.set nb331_p_krf, 56
.set nb331_p_crf, 60
.set nb331_Vc, 64
.set nb331_type, 68
.set nb331_p_ntype, 72
.set nb331_vdwparam, 76
.set nb331_Vvdw, 80
.set nb331_p_tabscale, 84
.set nb331_VFtab, 88
.set nb331_invsqrta, 92
.set nb331_dvda, 96
.set nb331_p_gbtabscale, 100
.set nb331_GBtab, 104
.set nb331_p_nthreads, 108
.set nb331_count, 112
.set nb331_mtx, 116
.set nb331_outeriter, 120
.set nb331_inneriter, 124
.set nb331_work, 128
                        ## stack offsets for local variables 
.set nb331_is3, 0
.set nb331_ii3, 4
.set nb331_ixO, 8
.set nb331_iyO, 12
.set nb331_izO, 16
.set nb331_ixH, 20
.set nb331_iyH, 28
.set nb331_izH, 36
.set nb331_iqO, 44
.set nb331_iqH, 52
.set nb331_qqO, 60
.set nb331_qqH, 68
.set nb331_vctot, 76
.set nb331_Vvdwtot, 84
.set nb331_c6, 92
.set nb331_c12, 100
.set nb331_two, 108
.set nb331_n1, 116
.set nb331_tsc, 124
.set nb331_ntia, 132
.set nb331_innerjjnr, 140
.set nb331_innerk, 144
.set nb331_fixO, 148
.set nb331_fiyO, 152
.set nb331_fizO, 156
.set nb331_fixH, 160
.set nb331_fiyH, 168
.set nb331_fizH, 176
.set nb331_dxO, 184
.set nb331_dyO, 188
.set nb331_dzO, 192
.set nb331_dxH, 196
.set nb331_dyH, 204
.set nb331_dzH, 212
.set nb331_tmprsqH, 220
.set nb331_n, 228                           ## idx for outer loop
.set nb331_nn1, 232                         ## number of outer iterations
.set nb331_nri, 236
.set nb331_nouter, 240
.set nb331_ninner, 244
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $248,%esp          ## local stack space 
        femms

        movl nb331_p_nri(%ebp),%ecx
        movl nb331_p_facel(%ebp),%esi
        movl nb331_p_tabscale(%ebp),%edi
        movl (%ecx),%ecx
        movl %ecx,nb331_nri(%esp)
        movd  (%esi),%mm1       ## facel
        movl nb331_p_ntype(%ebp),%esi

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb331_nouter(%esp)
        movl %eax,nb331_ninner(%esp)

        movl  nb331_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx),%ebx           ## ebx=ii 

        movl  nb331_charge(%ebp),%edx
        movd  (%edx,%ebx,4),%mm2    ## mm2=charge[ii0] 
        pfmul %mm1,%mm2
        movq  %mm2,nb331_iqO(%esp)          ## iqO = facel*charge[ii] 

        movd  4(%edx,%ebx,4),%mm2       ## mm2=charge[ii0+1] 
        pfmul %mm1,%mm2
        punpckldq %mm2,%mm2         ## spread to both halves 
        movq  %mm2,nb331_iqH(%esp)          ## iqH = facel*charge[ii0+1] 

        movl  nb331_type(%ebp),%edx
        movl  (%edx,%ebx,4),%ecx
        shll  %ecx
        imull (%esi),%ecx     ## ecx = ntia = 2*ntype*type[ii0]  
        movl  %ecx,nb331_ntia(%esp)

        movq  (%edi),%mm4
        punpckldq %mm4,%mm4         ## spread to both halves 
        movq  %mm4,nb331_tsc(%esp)
        movl $0x40000000,%eax
        movl %eax,nb331_two(%esp)
        movl %eax,nb331_two+4(%esp)

_nb_kernel331_ia32_3dnow.nb331_threadloop: 
        movl  nb331_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel331_ia32_3dnow.nb331_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel331_ia32_3dnow.nb331_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb331_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb331_n(%esp)
        movl %ebx,nb331_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel331_ia32_3dnow.nb331_outerstart
        jmp _nb_kernel331_ia32_3dnow.nb331_end

_nb_kernel331_ia32_3dnow.nb331_outerstart: 
        ## ebx contains number of outer iterations
        addl nb331_nouter(%esp),%ebx
        movl %ebx,nb331_nouter(%esp)

_nb_kernel331_ia32_3dnow.nb331_outer: 
        movl  nb331_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb331_is3(%esp)      ## store is3 

        movl  nb331_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movq  (%eax,%ebx,4),%mm5        ## move shX/shY to mm5 and shZ to mm6. 
        movd  8(%eax,%ebx,4),%mm6
        movq  %mm5,%mm0
        movq  %mm5,%mm1
        movq  %mm6,%mm2
        punpckldq %mm0,%mm0         ## also expand shX,Y,Z in mm0--mm2. 
        punpckhdq %mm1,%mm1
        punpckldq %mm2,%mm2

        movl  nb331_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx,%esi,4),%ebx    ## ebx=ii 

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb331_pos(%ebp),%eax      ## eax = base of pos[] 

        pfadd (%eax,%ebx,4),%mm5    ## ix = shX + posX (and iy too) 
        movd  8(%eax,%ebx,4),%mm7       ## cant use direct memory add for 4 bytes (iz) 
        movl  %ebx,nb331_ii3(%esp)          ## (use mm7 as temp. storage for iz.) 
        pfadd %mm7,%mm6
        movq  %mm5,nb331_ixO(%esp)
        movq  %mm6,nb331_izO(%esp)

        movd  12(%eax,%ebx,4),%mm3
        movd  16(%eax,%ebx,4),%mm4
        movd  20(%eax,%ebx,4),%mm5
        punpckldq  24(%eax,%ebx,4),%mm3
        punpckldq  28(%eax,%ebx,4),%mm4
        punpckldq  32(%eax,%ebx,4),%mm5    ## coords of H1 in low mm3-mm5, H2 in high 

        pfadd %mm3,%mm0
        pfadd %mm4,%mm1
        pfadd %mm5,%mm2
        movq %mm0,nb331_ixH(%esp)
        movq %mm1,nb331_iyH(%esp)
        movq %mm2,nb331_izH(%esp)

        ## clear vctot and i forces 
        pxor  %mm7,%mm7
        movq  %mm7,nb331_vctot(%esp)
        movq  %mm7,nb331_Vvdwtot(%esp)
        movq  %mm7,nb331_fixO(%esp)
        movd  %mm7,nb331_fizO(%esp)
        movq  %mm7,nb331_fixH(%esp)
        movq  %mm7,nb331_fiyH(%esp)
        movq  %mm7,nb331_fizH(%esp)

        movl  nb331_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx     ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 
        movl  %edx,nb331_innerk(%esp)
        addl  nb331_ninner(%esp),%edx
        movl  %edx,nb331_ninner(%esp)

        movl  nb331_pos(%ebp),%esi
        movl  nb331_faction(%ebp),%edi
        movl  nb331_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb331_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
_nb_kernel331_ia32_3dnow.nb331_inner_loop: 
        ## a single j particle iteration here - compare with the unrolled code for comments. 
        movl  nb331_innerjjnr(%esp),%eax
        movl  (%eax),%eax       ## eax=jnr offset 
    addl $4,nb331_innerjjnr(%esp)             ## advance pointer 
        prefetch 16(%ecx)          ## prefetch data - trial and error says 16 is best 

        movl nb331_charge(%ebp),%ecx
        movd (%ecx,%eax,4),%mm7
        punpckldq %mm7,%mm7
        movq %mm7,%mm6
        pfmul nb331_iqO(%esp),%mm6
        pfmul nb331_iqH(%esp),%mm7      ## mm6=qqO, mm7=qqH 
        movd %mm6,nb331_qqO(%esp)
        movq %mm7,nb331_qqH(%esp)

        movl nb331_type(%ebp),%ecx
        movl (%ecx,%eax,4),%edx          ## type [jnr] 
        movl nb331_vdwparam(%ebp),%ecx
        shll %edx
        addl nb331_ntia(%esp),%edx           ## tja = ntia + 2*type 
        movd (%ecx,%edx,4),%mm5         ## mm5 = 1st c6                 
        movq %mm5,nb331_c6(%esp)
        movd 4(%ecx,%edx,4),%mm5        ## mm5 = 1st c12                
        movq %mm5,nb331_c12(%esp)

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

        pfsubr nb331_ixO(%esp),%mm0
        pfsubr nb331_izO(%esp),%mm1

        movq  %mm0,nb331_dxO(%esp)
        pfmul %mm0,%mm0
        movd  %mm1,nb331_dzO(%esp)
        pfmul %mm1,%mm1
        pfacc %mm1,%mm0
        pfadd %mm1,%mm0         ## mm0=rsqO 

        punpckldq %mm2,%mm2
        punpckldq %mm3,%mm3
        punpckldq %mm4,%mm4 ## mm2-mm4 is jx-jz 
        pfsubr nb331_ixH(%esp),%mm2
        pfsubr nb331_iyH(%esp),%mm3
        pfsubr nb331_izH(%esp),%mm4   ## mm2-mm4 is dxH-dzH 

        movq %mm2,nb331_dxH(%esp)
        movq %mm3,nb331_dyH(%esp)
        movq %mm4,nb331_dzH(%esp)
        pfmul %mm2,%mm2
        pfmul %mm3,%mm3
        pfmul %mm4,%mm4

        pfadd %mm2,%mm3
        pfadd %mm4,%mm3         ## mm3=rsqH 
        movq %mm3,nb331_tmprsqH(%esp)

    pfrsqrt %mm0,%mm1

    movq %mm1,%mm2
    pfmul %mm1,%mm1
    pfrsqit1 %mm0,%mm1
    pfrcpit2 %mm2,%mm1  ## mm1=invsqrt 

        pfmul %mm1,%mm0         ## mm0=r 

        pfmul nb331_tsc(%esp),%mm0
        pf2iw %mm0,%mm4
        movd %mm4,nb331_n1(%esp)
        pi2fd %mm4,%mm4
        pfsub %mm4,%mm0              ## now mm0 is eps and mm4 n0 
        movq  %mm0,%mm2
        pfmul %mm2,%mm2         ## mm0 is eps, mm2 eps2 

        ## coulomb table 
        movl nb331_VFtab(%ebp),%edx
        movl nb331_n1(%esp),%ecx
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

        pfmul nb331_two(%esp),%mm7      ## two*Heps2 
        pfadd %mm6,%mm7
        pfadd %mm5,%mm7 ## mm7=FF 

        pfmul %mm0,%mm5 ## mm5=eps*Fp 
        pfadd %mm4,%mm5 ##  mm5= VV 

        pfmul nb331_qqO(%esp),%mm5      ## vcoul=qq*VV 
        pfmul nb331_qqO(%esp),%mm7      ## fijC=qq*FF 

        ## update vctot directly, use mm3 for fscal sum. 
        pfadd nb331_vctot(%esp),%mm5
        movq %mm5,nb331_vctot(%esp)
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
        pfmul nb331_two(%esp),%mm7      ## two*Heps2 
        pfadd %mm6,%mm7
        pfadd %mm5,%mm7 ## mm7=FF 
        pfmul %mm0,%mm5 ## mm5=eps*Fp 
        pfadd %mm4,%mm5 ##  mm5= VV     

        movq nb331_c6(%esp),%mm4
        pfmul %mm4,%mm7 ## fijD 
        pfmul %mm4,%mm5 ## Vvdw6            
        pfadd %mm7,%mm3 ## add to fscal  

        ## update Vvdwtot to release mm5! 
        pfadd nb331_Vvdwtot(%esp),%mm5        ## add the earlier value 
        movq %mm5,nb331_Vvdwtot(%esp)         ## store the sum       

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
        pfmul nb331_two(%esp),%mm7      ## two*Heps2 
        pfadd %mm6,%mm7
        pfadd %mm5,%mm7 ## mm7=FF 
        pfmul %mm0,%mm5 ## mm5=eps*Fp 
        pfadd %mm4,%mm5 ##  mm5= VV 

        movq nb331_c12(%esp),%mm6
        pfmul %mm6,%mm7 ## fijR 
        pfmul %mm6,%mm5 ## Vvdw12 
        pfadd %mm7,%mm3 ## total fscal fijC+ fijD+ fijR 

        ## change sign of fscal and multiply with rinv  
    pxor %mm0,%mm0
        pfsubr %mm0,%mm3
        pfmul nb331_tsc(%esp),%mm3
        pfmul %mm1,%mm3   ## mm3 is total fscal (for the oxygen) now 

        ## update Vvdwtot  
        pfadd nb331_Vvdwtot(%esp),%mm5        ## add the earlier value 
        movq %mm5,nb331_Vvdwtot(%esp)         ## store the sum       

        ## Ready with the oxygen - potential is updated, fscal is in mm3. 
        ## now do the two hydrogens. 
        movq nb331_tmprsqH(%esp),%mm0   ## mm0=rsqH 

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
        pfmul nb331_tsc(%esp),%mm0
        pf2iw %mm0,%mm4
        movq %mm4,nb331_n1(%esp)
        pi2fd %mm4,%mm4
        pfsub %mm4,%mm0              ## now mm0 is eps and mm4 n0 
        movq  %mm0,%mm2
        pfmul %mm2,%mm2         ## mm0 is eps, mm2 eps2 

        ## coulomb table 
        movl nb331_VFtab(%ebp),%edx
        movl nb331_n1(%esp),%ecx
        leal (%ecx,%ecx,2),%ecx
        shll $2,%ecx
        ## load all values we need 
        movd (%edx,%ecx,4),%mm4
        movd 4(%edx,%ecx,4),%mm5
        movd 8(%edx,%ecx,4),%mm6
        movd 12(%edx,%ecx,4),%mm7
        movl nb331_n1+4(%esp),%ecx
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

        pfmul nb331_two(%esp),%mm7      ## two*Heps2 
        pfadd %mm6,%mm7
        pfadd %mm5,%mm7 ## mm7=FF 

        pfmul %mm0,%mm5 ## mm5=eps*Fp 
        pfadd %mm4,%mm5 ##  mm5= VV 

        pfmul nb331_qqH(%esp),%mm5      ## vcoul=qq*VV 
        pfmul nb331_qqH(%esp),%mm7      ## fijC=qq*FF 
        ## update vctot 
        pfadd nb331_vctot(%esp),%mm5
        movq %mm5,nb331_vctot(%esp)

        ## change sign of fijC and multiply by rinv 
    pxor %mm4,%mm4
        pfsub %mm7,%mm4
        pfmul nb331_tsc(%esp),%mm4
        pfmul %mm1,%mm4   ## mm4 is total fscal (for the hydrogens) now         

        ## spread oxygen fscalar to both positions 
        punpckldq %mm3,%mm3
        ## calc vectorial force for O 
        prefetchw (%edi,%eax,4) ## prefetch faction to cache  
        movq nb331_dxO(%esp),%mm0
        movd nb331_dzO(%esp),%mm1
        pfmul %mm3,%mm0
        pfmul %mm3,%mm1

        ## calc vectorial force for H's 
        movq nb331_dxH(%esp),%mm5
        movq nb331_dyH(%esp),%mm6
        movq nb331_dzH(%esp),%mm7
        pfmul %mm4,%mm5
        pfmul %mm4,%mm6
        pfmul %mm4,%mm7

        ## update iO particle force 
        movq nb331_fixO(%esp),%mm2
        movd nb331_fizO(%esp),%mm3
        pfadd %mm0,%mm2
        pfadd %mm1,%mm3
        movq %mm2,nb331_fixO(%esp)
        movd %mm3,nb331_fizO(%esp)

        ## update iH forces 
        movq nb331_fixH(%esp),%mm2
        movq nb331_fiyH(%esp),%mm3
        movq nb331_fizH(%esp),%mm4
        pfadd %mm5,%mm2
        pfadd %mm6,%mm3
        pfadd %mm7,%mm4
        movq %mm2,nb331_fixH(%esp)
        movq %mm3,nb331_fiyH(%esp)
        movq %mm4,nb331_fizH(%esp)

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
        decl nb331_innerk(%esp)
        jz  _nb_kernel331_ia32_3dnow.nb331_updateouterdata
        jmp _nb_kernel331_ia32_3dnow.nb331_inner_loop
_nb_kernel331_ia32_3dnow.nb331_updateouterdata: 
        movl  nb331_ii3(%esp),%ecx

        movq  (%edi,%ecx,4),%mm6       ## increment iO force  
        movd  8(%edi,%ecx,4),%mm7
        pfadd nb331_fixO(%esp),%mm6
        pfadd nb331_fizO(%esp),%mm7
        movq  %mm6,(%edi,%ecx,4)
        movd  %mm7,8(%edi,%ecx,4)

        movq  nb331_fixH(%esp),%mm0
        movq  nb331_fiyH(%esp),%mm3
        movq  nb331_fizH(%esp),%mm1
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


        movl  nb331_fshift(%ebp),%ebx      ## increment fshift force 
        movl  nb331_is3(%esp),%edx

        movq  (%ebx,%edx,4),%mm6
        movd  8(%ebx,%edx,4),%mm7
        pfadd nb331_fixO(%esp),%mm6
        pfadd nb331_fizO(%esp),%mm7
        pfadd %mm0,%mm6
        pfadd %mm1,%mm7
        pfadd %mm2,%mm6
        pfadd %mm3,%mm7
        movq  %mm6,(%ebx,%edx,4)
        movd  %mm7,8(%ebx,%edx,4)

        ## get n from stack
        movl nb331_n(%esp),%esi
        ## get group index for i particle 
        movl  nb331_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        movq  nb331_vctot(%esp),%mm7
        pfacc %mm7,%mm7           ## get and sum the two parts of total potential 

        movl  nb331_Vc(%ebp),%eax
        movd  (%eax,%edx,4),%mm6
        pfadd %mm7,%mm6
        movd  %mm6,(%eax,%edx,4)          ## increment vc[gid] 

        movq  nb331_Vvdwtot(%esp),%mm7
        pfacc %mm7,%mm7           ## same for Vvdw 

        movl  nb331_Vvdw(%ebp),%eax
        movd  (%eax,%edx,4),%mm6
        pfadd %mm7,%mm6
        movd  %mm6,(%eax,%edx,4)          ## increment Vvdw[gid] 
        ## finish if last 
        movl nb331_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jecxz _nb_kernel331_ia32_3dnow.nb331_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb331_n(%esp)
        jmp _nb_kernel331_ia32_3dnow.nb331_outer
_nb_kernel331_ia32_3dnow.nb331_outerend: 
        ## check if more outer neighborlists remain
        movl  nb331_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jecxz _nb_kernel331_ia32_3dnow.nb331_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel331_ia32_3dnow.nb331_threadloop
_nb_kernel331_ia32_3dnow.nb331_end: 
        femms
        movl nb331_nouter(%esp),%eax
        movl nb331_ninner(%esp),%ebx
        movl nb331_outeriter(%ebp),%ecx
        movl nb331_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        addl $248,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret




.globl nb_kernel331nf_ia32_3dnow
.globl _nb_kernel331nf_ia32_3dnow
nb_kernel331nf_ia32_3dnow:      
_nb_kernel331nf_ia32_3dnow:     
.set nb331nf_p_nri, 8
.set nb331nf_iinr, 12
.set nb331nf_jindex, 16
.set nb331nf_jjnr, 20
.set nb331nf_shift, 24
.set nb331nf_shiftvec, 28
.set nb331nf_fshift, 32
.set nb331nf_gid, 36
.set nb331nf_pos, 40
.set nb331nf_faction, 44
.set nb331nf_charge, 48
.set nb331nf_p_facel, 52
.set nb331nf_p_krf, 56
.set nb331nf_p_crf, 60
.set nb331nf_Vc, 64
.set nb331nf_type, 68
.set nb331nf_p_ntype, 72
.set nb331nf_vdwparam, 76
.set nb331nf_Vvdw, 80
.set nb331nf_p_tabscale, 84
.set nb331nf_VFtab, 88
.set nb331nf_invsqrta, 92
.set nb331nf_dvda, 96
.set nb331nf_p_gbtabscale, 100
.set nb331nf_GBtab, 104
.set nb331nf_p_nthreads, 108
.set nb331nf_count, 112
.set nb331nf_mtx, 116
.set nb331nf_outeriter, 120
.set nb331nf_inneriter, 124
.set nb331nf_work, 128
                        ## stack offsets for local variables 
.set nb331nf_is3, 0
.set nb331nf_ii3, 4
.set nb331nf_ixO, 8
.set nb331nf_iyO, 12
.set nb331nf_izO, 16
.set nb331nf_ixH, 20
.set nb331nf_iyH, 28
.set nb331nf_izH, 36
.set nb331nf_iqO, 44
.set nb331nf_iqH, 52
.set nb331nf_qqO, 60
.set nb331nf_qqH, 68
.set nb331nf_vctot, 76
.set nb331nf_Vvdwtot, 84
.set nb331nf_c6, 92
.set nb331nf_c12, 100
.set nb331nf_n1, 108
.set nb331nf_tsc, 116
.set nb331nf_ntia, 124
.set nb331nf_innerjjnr, 128
.set nb331nf_innerk, 132
.set nb331nf_tmprsqH, 136
.set nb331nf_n, 144                         ## idx for outer loop
.set nb331nf_nn1, 148                       ## number of outer iterations
.set nb331nf_nri, 152
.set nb331nf_nouter, 156
.set nb331nf_ninner, 160
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $164,%esp          ## local stack space 
        femms

        movl nb331nf_p_nri(%ebp),%ecx
        movl nb331nf_p_facel(%ebp),%esi
        movl nb331nf_p_tabscale(%ebp),%edi
        movl (%ecx),%ecx
        movl %ecx,nb331nf_nri(%esp)
        movd  (%esi),%mm1       ## facel
        movl nb331nf_p_ntype(%ebp),%esi

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb331nf_nouter(%esp)
        movl %eax,nb331nf_ninner(%esp)

        movl  nb331nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx),%ebx           ## ebx=ii 

        movl  nb331nf_charge(%ebp),%edx

        movd  (%edx,%ebx,4),%mm2    ## mm2=charge[ii0] 
        pfmul %mm1,%mm2
        movq  %mm2,nb331nf_iqO(%esp)        ## iqO = facel*charge[ii] 

        movd  4(%edx,%ebx,4),%mm2       ## mm2=charge[ii0+1] 
        pfmul %mm1,%mm2
        punpckldq %mm2,%mm2         ## spread to both halves 
        movq  %mm2,nb331nf_iqH(%esp)        ## iqH = facel*charge[ii0+1] 

        movl  nb331nf_type(%ebp),%edx
        movl  (%edx,%ebx,4),%ecx
        shll  %ecx
        imull (%esi),%ecx     ## ecx = ntia = 2*ntype*type[ii0]  
        movl  %ecx,nb331nf_ntia(%esp)

        movq  (%edi),%mm4
        punpckldq %mm4,%mm4         ## spread to both halves 
        movq  %mm4,nb331nf_tsc(%esp)
_nb_kernel331nf_ia32_3dnow.nb331nf_threadloop: 
        movl  nb331nf_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel331nf_ia32_3dnow.nb331nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel331nf_ia32_3dnow.nb331nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb331nf_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb331nf_n(%esp)
        movl %ebx,nb331nf_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel331nf_ia32_3dnow.nb331nf_outerstart
        jmp _nb_kernel331nf_ia32_3dnow.nb331nf_end

_nb_kernel331nf_ia32_3dnow.nb331nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb331nf_nouter(%esp),%ebx
        movl %ebx,nb331nf_nouter(%esp)

_nb_kernel331nf_ia32_3dnow.nb331nf_outer: 
        movl  nb331nf_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb331nf_is3(%esp)            ## store is3 

        movl  nb331nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movq  (%eax,%ebx,4),%mm5        ## move shX/shY to mm5 and shZ to mm6. 
        movd  8(%eax,%ebx,4),%mm6
        movq  %mm5,%mm0
        movq  %mm5,%mm1
        movq  %mm6,%mm2
        punpckldq %mm0,%mm0         ## also expand shX,Y,Z in mm0--mm2. 
        punpckhdq %mm1,%mm1
        punpckldq %mm2,%mm2

        movl  nb331nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx,%esi,4),%ebx    ## ebx=ii 

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb331nf_pos(%ebp),%eax      ## eax = base of pos[] 

        pfadd (%eax,%ebx,4),%mm5    ## ix = shX + posX (and iy too) 
        movd  8(%eax,%ebx,4),%mm7       ## cant use direct memory add for 4 bytes (iz) 
        movl  %ebx,nb331nf_ii3(%esp)        ## (use mm7 as temp. storage for iz.) 
        pfadd %mm7,%mm6
        movq  %mm5,nb331nf_ixO(%esp)
        movq  %mm6,nb331nf_izO(%esp)

        movd  12(%eax,%ebx,4),%mm3
        movd  16(%eax,%ebx,4),%mm4
        movd  20(%eax,%ebx,4),%mm5
        punpckldq  24(%eax,%ebx,4),%mm3
        punpckldq  28(%eax,%ebx,4),%mm4
        punpckldq  32(%eax,%ebx,4),%mm5    ## coords of H1 in low mm3-mm5, H2 in high 

        pfadd %mm3,%mm0
        pfadd %mm4,%mm1
        pfadd %mm5,%mm2
        movq %mm0,nb331nf_ixH(%esp)
        movq %mm1,nb331nf_iyH(%esp)
        movq %mm2,nb331nf_izH(%esp)

        ## clear vctot and i forces 
        pxor  %mm7,%mm7
        movq  %mm7,nb331nf_vctot(%esp)
        movq  %mm7,nb331nf_Vvdwtot(%esp)

        movl  nb331nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 
        movl  %edx,nb331nf_innerk(%esp)
        addl  nb331nf_ninner(%esp),%edx
        movl  %edx,nb331nf_ninner(%esp)

        movl  nb331nf_pos(%ebp),%esi
        movl  nb331nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb331nf_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
_nb_kernel331nf_ia32_3dnow.nb331nf_inner_loop: 
        ## a single j particle iteration here - compare with the unrolled code for comments. 
        movl  nb331nf_innerjjnr(%esp),%eax
        movl  (%eax),%eax       ## eax=jnr offset 
    addl $4,nb331nf_innerjjnr(%esp)             ## advance pointer 
        prefetch 16(%ecx)          ## prefetch data - trial and error says 16 is best 

        movl nb331nf_charge(%ebp),%ecx
        movd (%ecx,%eax,4),%mm7
        punpckldq %mm7,%mm7
        movq %mm7,%mm6
        pfmul nb331nf_iqO(%esp),%mm6
        pfmul nb331nf_iqH(%esp),%mm7    ## mm6=qqO, mm7=qqH 
        movd %mm6,nb331nf_qqO(%esp)
        movq %mm7,nb331nf_qqH(%esp)

        movl nb331nf_type(%ebp),%ecx
        movl (%ecx,%eax,4),%edx          ## type [jnr] 
        movl nb331nf_vdwparam(%ebp),%ecx
        shll %edx
        addl nb331nf_ntia(%esp),%edx         ## tja = ntia + 2*type 
        movd (%ecx,%edx,4),%mm5         ## mm5 = 1st c6                 
        movq %mm5,nb331nf_c6(%esp)
        movd 4(%ecx,%edx,4),%mm5        ## mm5 = 1st c12                
        movq %mm5,nb331nf_c12(%esp)

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

        pfsubr nb331nf_ixO(%esp),%mm0
        pfsubr nb331nf_izO(%esp),%mm1

        pfmul %mm0,%mm0
        pfmul %mm1,%mm1
        pfacc %mm1,%mm0
        pfadd %mm1,%mm0         ## mm0=rsqO 

        punpckldq %mm2,%mm2
        punpckldq %mm3,%mm3
        punpckldq %mm4,%mm4 ## mm2-mm4 is jx-jz 
        pfsubr nb331nf_ixH(%esp),%mm2
        pfsubr nb331nf_iyH(%esp),%mm3
        pfsubr nb331nf_izH(%esp),%mm4   ## mm2-mm4 is dxH-dzH 

        pfmul %mm2,%mm2
        pfmul %mm3,%mm3
        pfmul %mm4,%mm4

        pfadd %mm2,%mm3
        pfadd %mm4,%mm3         ## mm3=rsqH 
        movq %mm3,nb331nf_tmprsqH(%esp)

    pfrsqrt %mm0,%mm1

    movq %mm1,%mm2
    pfmul %mm1,%mm1
    pfrsqit1 %mm0,%mm1
    pfrcpit2 %mm2,%mm1  ## mm1=invsqrt 

        pfmul %mm1,%mm0         ## mm0=r 

        pfmul nb331nf_tsc(%esp),%mm0
        pf2iw %mm0,%mm4
        movd %mm4,nb331nf_n1(%esp)
        pi2fd %mm4,%mm4
        pfsub %mm4,%mm0              ## now mm0 is eps and mm4 n0 
        movq  %mm0,%mm2
        pfmul %mm2,%mm2         ## mm0 is eps, mm2 eps2 

        ## coulomb table 
        movl nb331nf_VFtab(%ebp),%edx
        movl nb331nf_n1(%esp),%ecx
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

        pfmul nb331nf_qqO(%esp),%mm5    ## vcoul=qq*VV 

        ## update vctot directly 
        pfadd nb331nf_vctot(%esp),%mm5
        movq %mm5,nb331nf_vctot(%esp)

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

        movq nb331nf_c6(%esp),%mm4
        pfmul %mm4,%mm5 ## Vvdw6            
        ## update Vvdwtot to release mm5! 
        pfadd nb331nf_Vvdwtot(%esp),%mm5        ## add the earlier value 
        movq %mm5,nb331nf_Vvdwtot(%esp)         ## store the sum       

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

        movq nb331nf_c12(%esp),%mm6
        pfmul %mm6,%mm5 ## Vvdw12 
        ## update Vvdwtot  
        pfadd nb331nf_Vvdwtot(%esp),%mm5        ## add the earlier value 
        movq %mm5,nb331nf_Vvdwtot(%esp)         ## store the sum       

        ## now do the two hydrogens. 
        movq nb331nf_tmprsqH(%esp),%mm0   ## mm0=rsqH 

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
        pfmul nb331nf_tsc(%esp),%mm0
        pf2iw %mm0,%mm4
        movq %mm4,nb331nf_n1(%esp)
        pi2fd %mm4,%mm4
        pfsub %mm4,%mm0              ## now mm0 is eps and mm4 n0 
        movq  %mm0,%mm2
        pfmul %mm2,%mm2         ## mm0 is eps, mm2 eps2 

        ## coulomb table 
        movl nb331nf_VFtab(%ebp),%edx
        movl nb331nf_n1(%esp),%ecx
        leal (%ecx,%ecx,2),%ecx
        shll $2,%ecx
        ## load all values we need 
        movd (%edx,%ecx,4),%mm4
        movd 4(%edx,%ecx,4),%mm5
        movd 8(%edx,%ecx,4),%mm6
        movd 12(%edx,%ecx,4),%mm7
        movl nb331nf_n1+4(%esp),%ecx
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

        pfmul nb331nf_qqH(%esp),%mm5    ## vcoul=qq*VV 
        ## update vctot 
        pfadd nb331nf_vctot(%esp),%mm5
        movq %mm5,nb331nf_vctot(%esp)

        ##  done  - one more? 
        decl nb331nf_innerk(%esp)
        jz  _nb_kernel331nf_ia32_3dnow.nb331nf_updateouterdata
        jmp _nb_kernel331nf_ia32_3dnow.nb331nf_inner_loop
_nb_kernel331nf_ia32_3dnow.nb331nf_updateouterdata: 
        ## get n from stack
        movl nb331nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb331nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        movq  nb331nf_vctot(%esp),%mm7
        pfacc %mm7,%mm7           ## get and sum the two parts of total potential 

        movl  nb331nf_Vc(%ebp),%eax
        movd  (%eax,%edx,4),%mm6
        pfadd %mm7,%mm6
        movd  %mm6,(%eax,%edx,4)          ## increment vc[gid] 

        movq  nb331nf_Vvdwtot(%esp),%mm7
        pfacc %mm7,%mm7           ## same for Vvdw 

        movl  nb331nf_Vvdw(%ebp),%eax
        movd  (%eax,%edx,4),%mm6
        pfadd %mm7,%mm6
        movd  %mm6,(%eax,%edx,4)          ## increment Vvdw[gid] 
        ## finish if last 
        movl nb331nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jecxz _nb_kernel331nf_ia32_3dnow.nb331nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb331nf_n(%esp)
        jmp _nb_kernel331nf_ia32_3dnow.nb331nf_outer
_nb_kernel331nf_ia32_3dnow.nb331nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb331nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jecxz _nb_kernel331nf_ia32_3dnow.nb331nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel331nf_ia32_3dnow.nb331nf_threadloop
_nb_kernel331nf_ia32_3dnow.nb331nf_end: 
        femms
        movl nb331nf_nouter(%esp),%eax
        movl nb331nf_ninner(%esp),%ebx
        movl nb331nf_outeriter(%ebp),%ecx
        movl nb331nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        addl $164,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret



