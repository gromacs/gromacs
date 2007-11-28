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


.globl nb_kernel311_ia32_3dnow
.globl _nb_kernel311_ia32_3dnow
nb_kernel311_ia32_3dnow:        
_nb_kernel311_ia32_3dnow:       
.set nb311_p_nri, 8
.set nb311_iinr, 12
.set nb311_jindex, 16
.set nb311_jjnr, 20
.set nb311_shift, 24
.set nb311_shiftvec, 28
.set nb311_fshift, 32
.set nb311_gid, 36
.set nb311_pos, 40
.set nb311_faction, 44
.set nb311_charge, 48
.set nb311_p_facel, 52
.set nb311_p_krf, 56
.set nb311_p_crf, 60
.set nb311_Vc, 64
.set nb311_type, 68
.set nb311_p_ntype, 72
.set nb311_vdwparam, 76
.set nb311_Vvdw, 80
.set nb311_p_tabscale, 84
.set nb311_VFtab, 88
.set nb311_invsqrta, 92
.set nb311_dvda, 96
.set nb311_p_gbtabscale, 100
.set nb311_GBtab, 104
.set nb311_p_nthreads, 108
.set nb311_count, 112
.set nb311_mtx, 116
.set nb311_outeriter, 120
.set nb311_inneriter, 124
.set nb311_work, 128
                        ## stack offsets for local variables 
.set nb311_is3, 0
.set nb311_ii3, 4
.set nb311_ixO, 8
.set nb311_iyO, 12
.set nb311_izO, 16
.set nb311_ixH, 20
.set nb311_iyH, 28
.set nb311_izH, 36
.set nb311_iqO, 44
.set nb311_iqH, 52
.set nb311_qqO, 60
.set nb311_qqH, 68
.set nb311_vctot, 76
.set nb311_Vvdwtot, 84
.set nb311_c6, 92
.set nb311_c12, 100
.set nb311_six, 108
.set nb311_twelve, 116
.set nb311_two, 124
.set nb311_n1, 132
.set nb311_tsc, 140
.set nb311_ntia, 148
.set nb311_innerjjnr, 156
.set nb311_innerk, 160
.set nb311_fixO, 164
.set nb311_fiyO, 168
.set nb311_fizO, 172
.set nb311_fixH, 176
.set nb311_fiyH, 184
.set nb311_fizH, 192
.set nb311_dxO, 200
.set nb311_dyO, 204
.set nb311_dzO, 208
.set nb311_dxH, 212
.set nb311_dyH, 220
.set nb311_dzH, 228
.set nb311_tmprsqH, 236
.set nb311_n, 244                           ## idx for outer loop
.set nb311_nn1, 248                         ## number of outer iterations
.set nb311_nri, 252
.set nb311_nouter, 256
.set nb311_ninner, 260
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $264,%esp          ## local stack space 
        femms

        movl nb311_p_nri(%ebp),%ecx
        movl nb311_p_facel(%ebp),%esi
        movl nb311_p_tabscale(%ebp),%edi
        movl (%ecx),%ecx
        movl %ecx,nb311_nri(%esp)
        movd  (%esi),%mm1       ## facel
        movl nb311_p_ntype(%ebp),%esi

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb311_nouter(%esp)
        movl %eax,nb311_ninner(%esp)

        movl  nb311_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx),%ebx           ## ebx=ii 

        movl  nb311_charge(%ebp),%edx
        movd  (%edx,%ebx,4),%mm2    ## mm2=charge[ii0] 
        pfmul %mm1,%mm2
        movq  %mm2,nb311_iqO(%esp)          ## iqO = facel*charge[ii] 

        movd  4(%edx,%ebx,4),%mm2       ## mm2=charge[ii0+1] 
        pfmul %mm1,%mm2
        punpckldq %mm2,%mm2         ## spread to both halves 
        movq  %mm2,nb311_iqH(%esp)          ## iqH = facel*charge[ii0+1] 

        movl  nb311_type(%ebp),%edx
        movl  (%edx,%ebx,4),%edx
        shll  %edx
        movl  %edx,%ecx
        imull (%esi),%ecx     ## ecx = ntia = 2*ntype*type[ii0]  
        movl  %ecx,nb311_ntia(%esp)

        movq  (%edi),%mm6
        punpckldq %mm6,%mm6         ## spread to both halves 
        movq  %mm6,nb311_tsc(%esp)
        movl $0x40000000,%eax  ## 2.0
        movl %eax,nb311_two(%esp)
        movl %eax,nb311_two+4(%esp)
        movl $0x40c00000,%ebx ## 6.0
        movl %ebx,nb311_six(%esp)
        movl %ebx,nb311_six+4(%esp)
        movl $0x41400000,%ecx  ## 12.0
        movl %ecx,nb311_twelve(%esp)
        movl %ecx,nb311_twelve+4(%esp)

_nb_kernel311_ia32_3dnow.nb311_threadloop: 
        movl  nb311_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel311_ia32_3dnow.nb311_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel311_ia32_3dnow.nb311_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb311_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb311_n(%esp)
        movl %ebx,nb311_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel311_ia32_3dnow.nb311_outerstart
        jmp _nb_kernel311_ia32_3dnow.nb311_end

_nb_kernel311_ia32_3dnow.nb311_outerstart: 
        ## ebx contains number of outer iterations
        addl nb311_nouter(%esp),%ebx
        movl %ebx,nb311_nouter(%esp)

_nb_kernel311_ia32_3dnow.nb311_outer: 
        movl  nb311_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb311_is3(%esp)      ## store is3 

        movl  nb311_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movq  (%eax,%ebx,4),%mm5        ## move shX/shY to mm5 and shZ to mm6. 
        movd  8(%eax,%ebx,4),%mm6
        movq  %mm5,%mm0
        movq  %mm5,%mm1
        movq  %mm6,%mm2
        punpckldq %mm0,%mm0         ## also expand shX,Y,Z in mm0--mm2. 
        punpckhdq %mm1,%mm1
        punpckldq %mm2,%mm2

        movl  nb311_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx,%esi,4),%ebx    ## ebx=ii 

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb311_pos(%ebp),%eax      ## eax = base of pos[] 

        pfadd (%eax,%ebx,4),%mm5    ## ix = shX + posX (and iy too) 
        movd  8(%eax,%ebx,4),%mm7       ## cant use direct memory add for 4 bytes (iz) 
        movl  %ebx,nb311_ii3(%esp)          ## (use mm7 as temp. storage for iz.) 
        pfadd %mm7,%mm6
        movq  %mm5,nb311_ixO(%esp)
        movq  %mm6,nb311_izO(%esp)

        movd  12(%eax,%ebx,4),%mm3
        movd  16(%eax,%ebx,4),%mm4
        movd  20(%eax,%ebx,4),%mm5
        punpckldq  24(%eax,%ebx,4),%mm3
        punpckldq  28(%eax,%ebx,4),%mm4
        punpckldq  32(%eax,%ebx,4),%mm5    ## coords of H1 in low mm3-mm5, H2 in high 

        pfadd %mm3,%mm0
        pfadd %mm4,%mm1
        pfadd %mm5,%mm2
        movq %mm0,nb311_ixH(%esp)
        movq %mm1,nb311_iyH(%esp)
        movq %mm2,nb311_izH(%esp)

        ## clear vctot and i forces 
        pxor  %mm7,%mm7
        movq  %mm7,nb311_vctot(%esp)
        movq  %mm7,nb311_Vvdwtot(%esp)
        movq  %mm7,nb311_fixO(%esp)
        movd  %mm7,nb311_fizO(%esp)
        movq  %mm7,nb311_fixH(%esp)
        movq  %mm7,nb311_fiyH(%esp)
        movq  %mm7,nb311_fizH(%esp)

        movl  nb311_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 
        movl  %edx,nb311_innerk(%esp)
        addl  nb311_ninner(%esp),%edx
        movl  %edx,nb311_ninner(%esp)

        movl  nb311_pos(%ebp),%esi
        movl  nb311_faction(%ebp),%edi
        movl  nb311_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb311_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
_nb_kernel311_ia32_3dnow.nb311_inner_loop: 
        ## a single j particle iteration here - compare with the unrolled code for comments. 
        movl  nb311_innerjjnr(%esp),%eax
        movl  (%eax),%eax       ## eax=jnr offset 
    addl $4,nb311_innerjjnr(%esp)             ## advance pointer 
        prefetch 16(%ecx)          ## prefetch data - trial and error says 16 is best 

        movl nb311_charge(%ebp),%ecx
        movd (%ecx,%eax,4),%mm7
        punpckldq %mm7,%mm7
        movq %mm7,%mm6
        pfmul nb311_iqO(%esp),%mm6
        pfmul nb311_iqH(%esp),%mm7      ## mm6=qqO, mm7=qqH 
        movd %mm6,nb311_qqO(%esp)
        movq %mm7,nb311_qqH(%esp)

        movl nb311_type(%ebp),%ecx
        movl (%ecx,%eax,4),%edx          ## type [jnr] 
        movl nb311_vdwparam(%ebp),%ecx
        shll %edx
        addl nb311_ntia(%esp),%edx           ## tja = ntia + 2*type 
        movd (%ecx,%edx,4),%mm5         ## mm5 = 1st c6                 
        movq %mm5,nb311_c6(%esp)
        movd 4(%ecx,%edx,4),%mm5        ## mm5 = 1st c12                
        movq %mm5,nb311_c12(%esp)

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

        pfsubr nb311_ixO(%esp),%mm0
        pfsubr nb311_izO(%esp),%mm1

        movq  %mm0,nb311_dxO(%esp)
        pfmul %mm0,%mm0
        movd  %mm1,nb311_dzO(%esp)
        pfmul %mm1,%mm1
        pfacc %mm1,%mm0
        pfadd %mm1,%mm0         ## mm0=rsqO 

        punpckldq %mm2,%mm2
        punpckldq %mm3,%mm3
        punpckldq %mm4,%mm4 ## mm2-mm4 is jx-jz 
        pfsubr nb311_ixH(%esp),%mm2
        pfsubr nb311_iyH(%esp),%mm3
        pfsubr nb311_izH(%esp),%mm4   ## mm2-mm4 is dxH-dzH 

        movq %mm2,nb311_dxH(%esp)
        movq %mm3,nb311_dyH(%esp)
        movq %mm4,nb311_dzH(%esp)
        pfmul %mm2,%mm2
        pfmul %mm3,%mm3
        pfmul %mm4,%mm4

        pfadd %mm2,%mm3
        pfadd %mm4,%mm3         ## mm3=rsqH 
        movq %mm3,nb311_tmprsqH(%esp)

    pfrsqrt %mm0,%mm1

    movq %mm1,%mm2
    pfmul %mm1,%mm1
    pfrsqit1 %mm0,%mm1
    pfrcpit2 %mm2,%mm1  ## mm1=invsqrt 

        pfmul %mm1,%mm0         ## mm0=r 

        pfmul nb311_tsc(%esp),%mm0
        pf2iw %mm0,%mm4
        movd %mm4,nb311_n1(%esp)
        pi2fd %mm4,%mm4
        pfsub %mm4,%mm0              ## now mm0 is eps and mm4 n0 
        movq  %mm0,%mm2
        pfmul %mm2,%mm2         ## mm0 is eps, mm2 eps2 

        ## coulomb table 
        movl nb311_VFtab(%ebp),%edx
        movl nb311_n1(%esp),%ecx
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

        pfmul nb311_two(%esp),%mm7      ## two*Heps2 
        pfadd %mm6,%mm7
        pfadd %mm5,%mm7 ## mm7=FF 

        pfmul %mm0,%mm5 ## mm5=eps*Fp 
        pfadd %mm4,%mm5 ##  mm5= VV 

        pfmul nb311_qqO(%esp),%mm5      ## vcoul=qq*VV 
        pfmul nb311_qqO(%esp),%mm7      ## fijC=qq*FF 
        ## update vctot directly, use mm3 for fscal sum. 
        pfadd nb311_vctot(%esp),%mm5
        movq %mm5,nb311_vctot(%esp)

        movq %mm7,%mm3
        pfmul nb311_tsc(%esp),%mm3

        ## nontabulated LJ - mm1 is invsqrt. - keep mm1! 
        movq %mm1,%mm0
        pfmul %mm0,%mm0         ## mm0 is invsq 
        movq %mm0,%mm2
        pfmul %mm0,%mm2
        pfmul %mm0,%mm2         ## mm2 = rinvsix 
        movq %mm2,%mm4
        pfmul %mm4,%mm4         ## mm4=rinvtwelve 

        pfmul nb311_c12(%esp),%mm4
        pfmul nb311_c6(%esp),%mm2
        movq %mm4,%mm5
        pfsub %mm2,%mm5         ## mm5=Vvdw12-Vvdw6 

        pfmul nb311_six(%esp),%mm2
        pfmul nb311_twelve(%esp),%mm4
        pfsub %mm2,%mm4
        pfmul %mm1,%mm4   ## mm4=(12*Vvdw12-6*Vvdw6)*rinv11 

        pfsubr %mm4,%mm3
        pfmul %mm1,%mm3   ## mm3 is total fscal (for the oxygen) now 

        ## update Vvdwtot  
        pfadd nb311_Vvdwtot(%esp),%mm5        ## add the earlier value 
        movq %mm5,nb311_Vvdwtot(%esp)         ## store the sum       

        ## Ready with the oxygen - potential is updated, fscal is in mm3. 
        ## now do the two hydrogens. 
        movq nb311_tmprsqH(%esp),%mm0   ## mm0=rsqH 

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
        pfmul nb311_tsc(%esp),%mm0
        pf2iw %mm0,%mm4
        movq %mm4,nb311_n1(%esp)
        pi2fd %mm4,%mm4
        pfsub %mm4,%mm0              ## now mm0 is eps and mm4 n0 
        movq  %mm0,%mm2
        pfmul %mm2,%mm2         ## mm0 is eps, mm2 eps2 

        ## coulomb table 
        movl nb311_VFtab(%ebp),%edx
        movl nb311_n1(%esp),%ecx
        shll $2,%ecx
        ## load all values we need 
        movd (%edx,%ecx,4),%mm4
        movd 4(%edx,%ecx,4),%mm5
        movd 8(%edx,%ecx,4),%mm6
        movd 12(%edx,%ecx,4),%mm7
        movl nb311_n1+4(%esp),%ecx
        shll $2,%ecx
        punpckldq (%edx,%ecx,4),%mm4
        punpckldq 4(%edx,%ecx,4),%mm5
        punpckldq 8(%edx,%ecx,4),%mm6
        punpckldq 12(%edx,%ecx,4),%mm7

        pfmul %mm0,%mm6 ## mm6 = Geps           
        pfmul %mm2,%mm7 ## mm7 = Heps2 

        pfadd %mm6,%mm5
        pfadd %mm7,%mm5 ## mm5 = Fp 

        pfmul nb311_two(%esp),%mm7      ## two*Heps2 
        pfadd %mm6,%mm7
        pfadd %mm5,%mm7 ## mm7=FF 

        pfmul %mm0,%mm5 ## mm5=eps*Fp 
        pfadd %mm4,%mm5 ##  mm5= VV 

        pfmul nb311_qqH(%esp),%mm5      ## vcoul=qq*VV 
        pfmul nb311_qqH(%esp),%mm7      ## fijC=qq*FF 
        ## update vctot 
        pfadd nb311_vctot(%esp),%mm5
        movq %mm5,nb311_vctot(%esp)

        ## change sign of fijC and multiply by rinv 
    pxor %mm4,%mm4
        pfsub %mm7,%mm4
        pfmul nb311_tsc(%esp),%mm4
        pfmul %mm1,%mm4   ## mm4 is total fscal (for the hydrogens) now         

        ## spread oxygen fscalar to both positions 
        punpckldq %mm3,%mm3
        ## calc vectorial force for O 
        prefetchw (%edi,%eax,4) ## prefetch faction to cache  
        movq nb311_dxO(%esp),%mm0
        movd nb311_dzO(%esp),%mm1
        pfmul %mm3,%mm0
        pfmul %mm3,%mm1

        ## calc vectorial force for H's 
        movq nb311_dxH(%esp),%mm5
        movq nb311_dyH(%esp),%mm6
        movq nb311_dzH(%esp),%mm7
        pfmul %mm4,%mm5
        pfmul %mm4,%mm6
        pfmul %mm4,%mm7

        ## update iO particle force 
        movq nb311_fixO(%esp),%mm2
        movd nb311_fizO(%esp),%mm3
        pfadd %mm0,%mm2
        pfadd %mm1,%mm3
        movq %mm2,nb311_fixO(%esp)
        movd %mm3,nb311_fizO(%esp)

        ## update iH forces 
        movq nb311_fixH(%esp),%mm2
        movq nb311_fiyH(%esp),%mm3
        movq nb311_fizH(%esp),%mm4
        pfadd %mm5,%mm2
        pfadd %mm6,%mm3
        pfadd %mm7,%mm4
        movq %mm2,nb311_fixH(%esp)
        movq %mm3,nb311_fiyH(%esp)
        movq %mm4,nb311_fizH(%esp)

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
        decl nb311_innerk(%esp)
        jz  _nb_kernel311_ia32_3dnow.nb311_updateouterdata
        jmp _nb_kernel311_ia32_3dnow.nb311_inner_loop
_nb_kernel311_ia32_3dnow.nb311_updateouterdata: 
        movl  nb311_ii3(%esp),%ecx

        movq  (%edi,%ecx,4),%mm6       ## increment iO force  
        movd  8(%edi,%ecx,4),%mm7
        pfadd nb311_fixO(%esp),%mm6
        pfadd nb311_fizO(%esp),%mm7
        movq  %mm6,(%edi,%ecx,4)
        movd  %mm7,8(%edi,%ecx,4)

        movq  nb311_fixH(%esp),%mm0
        movq  nb311_fiyH(%esp),%mm3
        movq  nb311_fizH(%esp),%mm1
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


        movl  nb311_fshift(%ebp),%ebx      ## increment fshift force 
        movl  nb311_is3(%esp),%edx

        movq  (%ebx,%edx,4),%mm6
        movd  8(%ebx,%edx,4),%mm7
        pfadd nb311_fixO(%esp),%mm6
        pfadd nb311_fizO(%esp),%mm7
        pfadd %mm0,%mm6
        pfadd %mm1,%mm7
        pfadd %mm2,%mm6
        pfadd %mm3,%mm7
        movq  %mm6,(%ebx,%edx,4)
        movd  %mm7,8(%ebx,%edx,4)

        ## get n from stack
        movl nb311_n(%esp),%esi
        ## get group index for i particle 
        movl  nb311_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        movq  nb311_vctot(%esp),%mm7
        pfacc %mm7,%mm7           ## get and sum the two parts of total potential 

        movl  nb311_Vc(%ebp),%eax
        movd  (%eax,%edx,4),%mm6
        pfadd %mm7,%mm6
        movd  %mm6,(%eax,%edx,4)          ## increment vc[gid] 

        movq  nb311_Vvdwtot(%esp),%mm7
        pfacc %mm7,%mm7           ## same for Vvdw 

        movl  nb311_Vvdw(%ebp),%eax
        movd  (%eax,%edx,4),%mm6
        pfadd %mm7,%mm6
        movd  %mm6,(%eax,%edx,4)          ## increment Vvdw[gid] 
        ## finish if last 
        movl nb311_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel311_ia32_3dnow.nb311_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb311_n(%esp)
        jmp _nb_kernel311_ia32_3dnow.nb311_outer
_nb_kernel311_ia32_3dnow.nb311_outerend: 
        ## check if more outer neighborlists remain
        movl  nb311_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel311_ia32_3dnow.nb311_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel311_ia32_3dnow.nb311_threadloop
_nb_kernel311_ia32_3dnow.nb311_end: 
        femms
        movl nb311_nouter(%esp),%eax
        movl nb311_ninner(%esp),%ebx
        movl nb311_outeriter(%ebp),%ecx
        movl nb311_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        addl $264,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret


.globl nb_kernel311nf_ia32_3dnow
.globl _nb_kernel311nf_ia32_3dnow
nb_kernel311nf_ia32_3dnow:      
_nb_kernel311nf_ia32_3dnow:     
.set nb311nf_p_nri, 8
.set nb311nf_iinr, 12
.set nb311nf_jindex, 16
.set nb311nf_jjnr, 20
.set nb311nf_shift, 24
.set nb311nf_shiftvec, 28
.set nb311nf_fshift, 32
.set nb311nf_gid, 36
.set nb311nf_pos, 40
.set nb311nf_faction, 44
.set nb311nf_charge, 48
.set nb311nf_p_facel, 52
.set nb311nf_p_krf, 56
.set nb311nf_p_crf, 60
.set nb311nf_Vc, 64
.set nb311nf_type, 68
.set nb311nf_p_ntype, 72
.set nb311nf_vdwparam, 76
.set nb311nf_Vvdw, 80
.set nb311nf_p_tabscale, 84
.set nb311nf_VFtab, 88
.set nb311nf_invsqrta, 92
.set nb311nf_dvda, 96
.set nb311nf_p_gbtabscale, 100
.set nb311nf_GBtab, 104
.set nb311nf_p_nthreads, 108
.set nb311nf_count, 112
.set nb311nf_mtx, 116
.set nb311nf_outeriter, 120
.set nb311nf_inneriter, 124
.set nb311nf_work, 128
                        ## stack offsets for local variables 
.set nb311nf_is3, 0
.set nb311nf_ii3, 4
.set nb311nf_ixO, 8
.set nb311nf_iyO, 12
.set nb311nf_izO, 16
.set nb311nf_ixH, 20
.set nb311nf_iyH, 28
.set nb311nf_izH, 36
.set nb311nf_iqO, 44
.set nb311nf_iqH, 52
.set nb311nf_qqO, 60
.set nb311nf_qqH, 68
.set nb311nf_vctot, 76
.set nb311nf_Vvdwtot, 84
.set nb311nf_c6, 92
.set nb311nf_c12, 100
.set nb311nf_n1, 108
.set nb311nf_tsc, 116
.set nb311nf_ntia, 124
.set nb311nf_innerjjnr, 128
.set nb311nf_innerk, 132
.set nb311nf_tmprsqH, 136
.set nb311nf_n, 144                         ## idx for outer loop
.set nb311nf_nn1, 148                       ## number of outer iterations
.set nb311nf_nri, 152
.set nb311nf_nouter, 156
.set nb311nf_ninner, 160
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

        movl nb311nf_p_nri(%ebp),%ecx
        movl nb311nf_p_facel(%ebp),%esi
        movl nb311nf_p_tabscale(%ebp),%edi
        movl (%ecx),%ecx
        movl %ecx,nb311nf_nri(%esp)
        movd  (%esi),%mm1       ## facel
        movl nb311nf_p_ntype(%ebp),%esi

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb311nf_nouter(%esp)
        movl %eax,nb311nf_ninner(%esp)


        movl  nb311nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx),%ebx           ## ebx=ii 

        movl  nb311nf_charge(%ebp),%edx
        movd  (%edx,%ebx,4),%mm2    ## mm2=charge[ii0] 
        pfmul %mm1,%mm2
        movq  %mm2,nb311nf_iqO(%esp)        ## iqO = facel*charge[ii] 

        movd  4(%edx,%ebx,4),%mm2       ## mm2=charge[ii0+1] 
        pfmul %mm1,%mm2
        punpckldq %mm2,%mm2         ## spread to both halves 
        movq  %mm2,nb311nf_iqH(%esp)        ## iqH = facel*charge[ii0+1] 

        movl  nb311nf_type(%ebp),%edx
        movl  (%edx,%ebx,4),%edx
        shll  %edx
        movl  %edx,%ecx
        imull (%esi),%ecx     ## ecx = ntia = 2*ntype*type[ii0]  
        movl  %ecx,nb311nf_ntia(%esp)

        movq  (%edi),%mm6
        punpckldq %mm6,%mm6         ## spread to both halves 
        movq  %mm6,nb311nf_tsc(%esp)
_nb_kernel311nf_ia32_3dnow.nb311nf_threadloop: 
        movl  nb311nf_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel311nf_ia32_3dnow.nb311nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel311nf_ia32_3dnow.nb311nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb311nf_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb311nf_n(%esp)
        movl %ebx,nb311nf_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel311nf_ia32_3dnow.nb311nf_outerstart
        jmp _nb_kernel311nf_ia32_3dnow.nb311nf_end

_nb_kernel311nf_ia32_3dnow.nb311nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb311nf_nouter(%esp),%ebx
        movl %ebx,nb311nf_nouter(%esp)

_nb_kernel311nf_ia32_3dnow.nb311nf_outer: 
        movl  nb311nf_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb311nf_is3(%esp)            ## store is3 

        movl  nb311nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movq  (%eax,%ebx,4),%mm5        ## move shX/shY to mm5 and shZ to mm6. 
        movd  8(%eax,%ebx,4),%mm6
        movq  %mm5,%mm0
        movq  %mm5,%mm1
        movq  %mm6,%mm2
        punpckldq %mm0,%mm0         ## also expand shX,Y,Z in mm0--mm2. 
        punpckhdq %mm1,%mm1
        punpckldq %mm2,%mm2

        movl  nb311nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx,%esi,4),%ebx    ## ebx=ii 

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb311nf_pos(%ebp),%eax      ## eax = base of pos[] 

        pfadd (%eax,%ebx,4),%mm5    ## ix = shX + posX (and iy too) 
        movd  8(%eax,%ebx,4),%mm7       ## cant use direct memory add for 4 bytes (iz) 
        movl  %ebx,nb311nf_ii3(%esp)        ## (use mm7 as temp. storage for iz.) 
        pfadd %mm7,%mm6
        movq  %mm5,nb311nf_ixO(%esp)
        movq  %mm6,nb311nf_izO(%esp)

        movd  12(%eax,%ebx,4),%mm3
        movd  16(%eax,%ebx,4),%mm4
        movd  20(%eax,%ebx,4),%mm5
        punpckldq  24(%eax,%ebx,4),%mm3
        punpckldq  28(%eax,%ebx,4),%mm4
        punpckldq  32(%eax,%ebx,4),%mm5    ## coords of H1 in low mm3-mm5, H2 in high 

        pfadd %mm3,%mm0
        pfadd %mm4,%mm1
        pfadd %mm5,%mm2
        movq %mm0,nb311nf_ixH(%esp)
        movq %mm1,nb311nf_iyH(%esp)
        movq %mm2,nb311nf_izH(%esp)

        ## clear vctot and i forces 
        pxor  %mm7,%mm7
        movq  %mm7,nb311nf_vctot(%esp)
        movq  %mm7,nb311nf_Vvdwtot(%esp)

        movl  nb311nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx     ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 
        movl  %edx,nb311nf_innerk(%esp)
        addl  nb311nf_ninner(%esp),%edx
        movl  %edx,nb311nf_ninner(%esp)

        movl  nb311nf_pos(%ebp),%esi
        movl  nb311nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb311nf_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
_nb_kernel311nf_ia32_3dnow.nb311nf_inner_loop: 
        ## a single j particle iteration here - compare with the unrolled code for comments. 
        movl  nb311nf_innerjjnr(%esp),%eax
        movl  (%eax),%eax       ## eax=jnr offset 
        addl $4,nb311nf_innerjjnr(%esp)             ## advance pointer 
        prefetch 16(%ecx)          ## prefetch data - trial and error says 16 is best 

        movl nb311nf_charge(%ebp),%ecx
        movd (%ecx,%eax,4),%mm7
        punpckldq %mm7,%mm7
        movq %mm7,%mm6
        pfmul nb311nf_iqO(%esp),%mm6
        pfmul nb311nf_iqH(%esp),%mm7    ## mm6=qqO, mm7=qqH 
        movd %mm6,nb311nf_qqO(%esp)
        movq %mm7,nb311nf_qqH(%esp)

        movl nb311nf_type(%ebp),%ecx
        movl (%ecx,%eax,4),%edx          ## type [jnr] 
        movl nb311nf_vdwparam(%ebp),%ecx
        shll %edx
        addl nb311nf_ntia(%esp),%edx         ## tja = ntia + 2*type 
        movd (%ecx,%edx,4),%mm5         ## mm5 = 1st c6                 
        movq %mm5,nb311nf_c6(%esp)
        movd 4(%ecx,%edx,4),%mm5        ## mm5 = 1st c12                
        movq %mm5,nb311nf_c12(%esp)

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

        pfsubr nb311nf_ixO(%esp),%mm0
        pfsubr nb311nf_izO(%esp),%mm1

        pfmul %mm0,%mm0
        pfmul %mm1,%mm1
        pfacc %mm1,%mm0
        pfadd %mm1,%mm0         ## mm0=rsqO 

        punpckldq %mm2,%mm2
        punpckldq %mm3,%mm3
        punpckldq %mm4,%mm4 ## mm2-mm4 is jx-jz 
        pfsubr nb311nf_ixH(%esp),%mm2
        pfsubr nb311nf_iyH(%esp),%mm3
        pfsubr nb311nf_izH(%esp),%mm4   ## mm2-mm4 is dxH-dzH 

        pfmul %mm2,%mm2
        pfmul %mm3,%mm3
        pfmul %mm4,%mm4

        pfadd %mm2,%mm3
        pfadd %mm4,%mm3         ## mm3=rsqH 
        movq %mm3,nb311nf_tmprsqH(%esp)

        pfrsqrt %mm0,%mm1

        movq %mm1,%mm2
        pfmul %mm1,%mm1
        pfrsqit1 %mm0,%mm1
        pfrcpit2 %mm2,%mm1      ## mm1=invsqrt 

        pfmul %mm1,%mm0         ## mm0=r 

        pfmul nb311nf_tsc(%esp),%mm0
        pf2iw %mm0,%mm4
        movd %mm4,nb311nf_n1(%esp)
        pi2fd %mm4,%mm4
        pfsub %mm4,%mm0              ## now mm0 is eps and mm4 n0 
        movq  %mm0,%mm2
        pfmul %mm2,%mm2         ## mm0 is eps, mm2 eps2 

        ## coulomb table 
        movl nb311nf_VFtab(%ebp),%edx
        movl nb311nf_n1(%esp),%ecx
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

        pfmul nb311nf_qqO(%esp),%mm5    ## vcoul=qq*VV 
        ## update vctot directly 
        pfadd nb311nf_vctot(%esp),%mm5
        movq %mm5,nb311nf_vctot(%esp)

        ## nontabulated LJ - mm1 is invsqrt. - keep mm1! 
        movq %mm1,%mm0
        pfmul %mm0,%mm0         ## mm0 is invsq 
        movq %mm0,%mm2
        pfmul %mm0,%mm2
        pfmul %mm0,%mm2         ## mm2 = rinvsix 
        movq %mm2,%mm4
        pfmul %mm4,%mm4         ## mm4=rinvtwelve 

        pfmul nb311nf_c12(%esp),%mm4
        pfmul nb311nf_c6(%esp),%mm2
        pfsub %mm2,%mm4         ## mm4=Vvdw12-Vvdw6 

        ## update Vvdwtot  
        pfadd nb311nf_Vvdwtot(%esp),%mm4        ## add the earlier value 
        movq %mm4,nb311nf_Vvdwtot(%esp)         ## store the sum       

                ## now do the two hydrogens. 
        movq nb311nf_tmprsqH(%esp),%mm0   ## mm0=rsqH 

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
        pfmul nb311nf_tsc(%esp),%mm0
        pf2iw %mm0,%mm4
        movq %mm4,nb311nf_n1(%esp)
        pi2fd %mm4,%mm4
        pfsub %mm4,%mm0              ## now mm0 is eps and mm4 n0 
        movq  %mm0,%mm2
        pfmul %mm2,%mm2         ## mm0 is eps, mm2 eps2 

        ## coulomb table 
        movl nb311nf_VFtab(%ebp),%edx
        movl nb311nf_n1(%esp),%ecx
        shll $2,%ecx
        ## load all values we need 
        movd (%edx,%ecx,4),%mm4
        movd 4(%edx,%ecx,4),%mm5
        movd 8(%edx,%ecx,4),%mm6
        movd 12(%edx,%ecx,4),%mm7
        movl nb311nf_n1+4(%esp),%ecx
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

        pfmul nb311nf_qqH(%esp),%mm5    ## vcoul=qq*VV 
        ## update vctot 
        pfadd nb311nf_vctot(%esp),%mm5
        movq %mm5,nb311nf_vctot(%esp)

        ##  done  - one more? 
        decl nb311nf_innerk(%esp)
        jz  _nb_kernel311nf_ia32_3dnow.nb311nf_updateouterdata
        jmp _nb_kernel311nf_ia32_3dnow.nb311nf_inner_loop
_nb_kernel311nf_ia32_3dnow.nb311nf_updateouterdata: 
        ## get n from stack
        movl nb311nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb311nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        movq  nb311nf_vctot(%esp),%mm7
        pfacc %mm7,%mm7           ## get and sum the two parts of total potential 

        movl  nb311nf_Vc(%ebp),%eax
        movd  (%eax,%edx,4),%mm6
        pfadd %mm7,%mm6
        movd  %mm6,(%eax,%edx,4)          ## increment vc[gid] 

        movq  nb311nf_Vvdwtot(%esp),%mm7
        pfacc %mm7,%mm7           ## same for Vvdw 

        movl  nb311nf_Vvdw(%ebp),%eax
        movd  (%eax,%edx,4),%mm6
        pfadd %mm7,%mm6
        movd  %mm6,(%eax,%edx,4)          ## increment Vvdw[gid] 
        ## finish if last 
        movl nb311nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel311nf_ia32_3dnow.nb311nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb311nf_n(%esp)
        jmp _nb_kernel311nf_ia32_3dnow.nb311nf_outer
_nb_kernel311nf_ia32_3dnow.nb311nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb311nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel311nf_ia32_3dnow.nb311nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel311nf_ia32_3dnow.nb311nf_threadloop
_nb_kernel311nf_ia32_3dnow.nb311nf_end: 
        femms
        movl nb311nf_nouter(%esp),%eax
        movl nb311nf_ninner(%esp),%ebx
        movl nb311nf_outeriter(%ebp),%ecx
        movl nb311nf_inneriter(%ebp),%edx
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

