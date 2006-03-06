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




.globl nb_kernel111_ia32_3dnow
.globl _nb_kernel111_ia32_3dnow
nb_kernel111_ia32_3dnow:        
_nb_kernel111_ia32_3dnow:       
.set nb111_p_nri, 8
.set nb111_iinr, 12
.set nb111_jindex, 16
.set nb111_jjnr, 20
.set nb111_shift, 24
.set nb111_shiftvec, 28
.set nb111_fshift, 32
.set nb111_gid, 36
.set nb111_pos, 40
.set nb111_faction, 44
.set nb111_charge, 48
.set nb111_p_facel, 52
.set nb111_p_krf, 56
.set nb111_p_crf, 60
.set nb111_Vc, 64
.set nb111_type, 68
.set nb111_p_ntype, 72
.set nb111_vdwparam, 76
.set nb111_Vvdw, 80
.set nb111_p_tabscale, 84
.set nb111_VFtab, 88
.set nb111_invsqrta, 92
.set nb111_dvda, 96
.set nb111_p_gbtabscale, 100
.set nb111_GBtab, 104
.set nb111_p_nthreads, 108
.set nb111_count, 112
.set nb111_mtx, 116
.set nb111_outeriter, 120
.set nb111_inneriter, 124
.set nb111_work, 128
                        ## stack offsets for local variables 
.set nb111_is3, 0
.set nb111_ii3, 4
.set nb111_ixO, 8
.set nb111_iyO, 12
.set nb111_izO, 16
.set nb111_ixH, 20
.set nb111_iyH, 28
.set nb111_izH, 36
.set nb111_iqO, 44
.set nb111_iqH, 52
.set nb111_vctot, 60
.set nb111_Vvdwtot, 68
.set nb111_c6, 76
.set nb111_c12, 84
.set nb111_six, 92
.set nb111_twelve, 100
.set nb111_ntia, 108
.set nb111_innerjjnr, 116
.set nb111_innerk, 120
.set nb111_fixO, 124
.set nb111_fiyO, 128
.set nb111_fizO, 132
.set nb111_fixH, 136
.set nb111_fiyH, 144
.set nb111_fizH, 152
.set nb111_dxO, 160
.set nb111_dyO, 164
.set nb111_dzO, 168
.set nb111_dxH, 172
.set nb111_dyH, 180
.set nb111_dzH, 188
.set nb111_n, 196                           ## idx for outer loop
.set nb111_nn1, 200                         ## number of outer iterations
.set nb111_nri, 204
.set nb111_ntype, 208
.set nb111_nouter, 212
.set nb111_ninner, 216
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
        movl nb111_p_nri(%ebp),%ecx
        movl nb111_p_ntype(%ebp),%edx
        movl nb111_p_facel(%ebp),%esi
        movl (%ecx),%ecx
        movl (%edx),%edx
        movl %ecx,nb111_nri(%esp)
        movl %edx,nb111_ntype(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb111_nouter(%esp)
        movl %eax,nb111_ninner(%esp)

        ## assume we have at least one i particle - start directly      

        movl  nb111_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx),%ebx           ## ebx=ii 

        movl  nb111_charge(%ebp),%edx
        movd  (%esi),%mm1
        movd  (%edx,%ebx,4),%mm2    ## mm2=charge[ii0] 
        pfmul %mm1,%mm2
        movq  %mm2,nb111_iqO(%esp)          ## iqO = facel*charge[ii] 

        movd  4(%edx,%ebx,4),%mm2       ## mm2=charge[ii0+1] 
        pfmul %mm1,%mm2
        punpckldq %mm2,%mm2         ## spread to both halves 
        movq  %mm2,nb111_iqH(%esp)          ## iqH = facel*charge[ii0+1] 

        movl  nb111_type(%ebp),%edx
        movl  (%edx,%ebx,4),%ecx
        shll  %ecx
        imull nb111_ntype(%esp),%ecx        ## ecx = ntia = 2*ntype*type[ii0]  
        movl  %ecx,nb111_ntia(%esp)

        ## move data to local stack  
        movl $0x40c00000,%eax ## fp 6.0
        movl $0x41400000,%ebx ## fp 12.0

        movl %eax,nb111_six(%esp)
        movl %eax,nb111_six+4(%esp)
        movl %ebx,nb111_twelve(%esp)
        movl %ebx,nb111_twelve+4(%esp)
_nb_kernel111_ia32_3dnow.nb111_threadloop: 
        movl  nb111_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel111_ia32_3dnow.nb111_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel111_ia32_3dnow.nb111_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb111_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb111_n(%esp)
        movl %ebx,nb111_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel111_ia32_3dnow.nb111_outerstart
        jmp _nb_kernel111_ia32_3dnow.nb111_end

_nb_kernel111_ia32_3dnow.nb111_outerstart: 
        ## ebx contains number of outer iterations
        addl nb111_nouter(%esp),%ebx
        movl %ebx,nb111_nouter(%esp)

_nb_kernel111_ia32_3dnow.nb111_outer: 
        movl  nb111_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb111_is3(%esp)      ## store is3 

        movl  nb111_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movq  (%eax,%ebx,4),%mm5        ## move shX/shY to mm5 and shZ to mm6. 
        movd  8(%eax,%ebx,4),%mm6
        movq  %mm5,%mm0
        movq  %mm5,%mm1
        movq  %mm6,%mm2
        punpckldq %mm0,%mm0         ## also expand shX,Y,Z in mm0--mm2. 
        punpckhdq %mm1,%mm1
        punpckldq %mm2,%mm2

        movl  nb111_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx,%esi,4),%ebx            ## ebx=ii 

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb111_pos(%ebp),%eax      ## eax = base of pos[] 

        pfadd (%eax,%ebx,4),%mm5    ## ix = shX + posX (and iy too) 
        movd  8(%eax,%ebx,4),%mm7       ## cant use direct memory add for 4 bytes (iz) 
        movl  %ebx,nb111_ii3(%esp)          ## (use mm7 as temp. storage for iz.) 
        pfadd %mm7,%mm6
        movq  %mm5,nb111_ixO(%esp)
        movq  %mm6,nb111_izO(%esp)

        movd  12(%eax,%ebx,4),%mm3
        movd  16(%eax,%ebx,4),%mm4
        movd  20(%eax,%ebx,4),%mm5
        punpckldq  24(%eax,%ebx,4),%mm3
        punpckldq  28(%eax,%ebx,4),%mm4
        punpckldq  32(%eax,%ebx,4),%mm5    ## coords of H1 in low mm3-mm5, H2 in high 

        pfadd %mm3,%mm0
        pfadd %mm4,%mm1
        pfadd %mm5,%mm2
        movq %mm0,nb111_ixH(%esp)
        movq %mm1,nb111_iyH(%esp)
        movq %mm2,nb111_izH(%esp)

        ## clear vctot and i forces 
        pxor  %mm7,%mm7
        movq  %mm7,nb111_vctot(%esp)
        movq  %mm7,nb111_Vvdwtot(%esp)
        movq  %mm7,nb111_fixO(%esp)
        movd  %mm7,nb111_fizO(%esp)
        movq  %mm7,nb111_fixH(%esp)
        movq  %mm7,nb111_fiyH(%esp)
        movq  %mm7,nb111_fizH(%esp)

        movl  nb111_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 
        movl  %edx,nb111_innerk(%esp)      ## number of innerloop atoms 
        addl  nb111_ninner(%esp),%edx
        movl  %edx,nb111_ninner(%esp)

        movl  nb111_pos(%ebp),%esi
        movl  nb111_faction(%ebp),%edi
        movl  nb111_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb111_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
_nb_kernel111_ia32_3dnow.nb111_inner_loop: 
        ## a single j particle iteration here - compare with the unrolled code for comments. 
        movl  nb111_innerjjnr(%esp),%eax
        movl  (%eax),%eax       ## eax=jnr offset 
        addl $4,nb111_innerjjnr(%esp)             ## advance pointer 
        prefetch 16(%ecx)          ## prefetch data - trial and error says 16 is best 

        movl nb111_charge(%ebp),%ecx
        movd (%ecx,%eax,4),%mm7
        punpckldq %mm7,%mm7
        movq %mm7,%mm6
        pfmul nb111_iqO(%esp),%mm6
        pfmul nb111_iqH(%esp),%mm7      ## mm6=qqO, mm7=qqH 

        movl nb111_type(%ebp),%ecx
        movl (%ecx,%eax,4),%edx          ## type [jnr] 
        movl nb111_vdwparam(%ebp),%ecx
        shll %edx
        addl nb111_ntia(%esp),%edx           ## tja = ntia + 2*type 
        movd (%ecx,%edx,4),%mm5         ## mm5 = 1st c6                 
        movq %mm5,nb111_c6(%esp)
        movd 4(%ecx,%edx,4),%mm5        ## mm5 = 1st c12                
        movq %mm5,nb111_c12(%esp)

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

        pfsubr nb111_ixO(%esp),%mm0
        pfsubr nb111_izO(%esp),%mm1

        movq  %mm0,nb111_dxO(%esp)
        pfmul %mm0,%mm0
        movd  %mm1,nb111_dzO(%esp)
        pfmul %mm1,%mm1
        pfacc %mm1,%mm0
        pfadd %mm1,%mm0         ## mm0=rsqO 

        punpckldq %mm2,%mm2
        punpckldq %mm3,%mm3
        punpckldq %mm4,%mm4 ## mm2-mm4 is jx-jz 
        pfsubr nb111_ixH(%esp),%mm2
        pfsubr nb111_iyH(%esp),%mm3
        pfsubr nb111_izH(%esp),%mm4   ## mm2-mm4 is dxH-dzH 

        movq %mm2,nb111_dxH(%esp)
        movq %mm3,nb111_dyH(%esp)
        movq %mm4,nb111_dzH(%esp)
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

        movq  %mm4,%mm0
        pfmul %mm4,%mm0
        pfmul %mm4,%mm0         ## mm0=rinvsix 
        movq  %mm0,%mm2
        pfmul %mm2,%mm2         ## mm2=rintwelve 

        ## calculate potential and scalar force 
        pfmul %mm1,%mm6         ## mm6=vcoul 
        movq  %mm6,%mm1         ## use mm1 for fscal sum 

        ## LJ for the oxygen 
        pfmul nb111_c6(%esp),%mm0
        pfmul nb111_c12(%esp),%mm2

        ## calc nb potential 
        movq %mm2,%mm5
        pfsub %mm0,%mm5

        ## calc nb force 
        pfmul nb111_six(%esp),%mm0
        pfmul nb111_twelve(%esp),%mm2

        ## increment scalar force 
        pfsub %mm0,%mm1
        pfadd %mm2,%mm1
        pfmul %mm1,%mm4         ## total scalar force on oxygen. 

        ## update nb potential 
        pfadd nb111_Vvdwtot(%esp),%mm5
        movq %mm5,nb111_Vvdwtot(%esp)

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
        pfadd nb111_vctot(%esp),%mm7
        movq %mm7,nb111_vctot(%esp)

        ## spread oxygen fscalar to both positions 
        punpckldq %mm4,%mm4
        ## calc vectorial force for O 
        prefetchw (%edi,%eax,4) ## prefetch faction to cache  
        movq nb111_dxO(%esp),%mm0
        movd nb111_dzO(%esp),%mm1
        pfmul %mm4,%mm0
        pfmul %mm4,%mm1

        ## calc vectorial force for H's 
        movq nb111_dxH(%esp),%mm5
        movq nb111_dyH(%esp),%mm6
        movq nb111_dzH(%esp),%mm7
        pfmul %mm3,%mm5
        pfmul %mm3,%mm6
        pfmul %mm3,%mm7

        ## update iO particle force 
        movq nb111_fixO(%esp),%mm2
        movd nb111_fizO(%esp),%mm3
        pfadd %mm0,%mm2
        pfadd %mm1,%mm3
        movq %mm2,nb111_fixO(%esp)
        movd %mm3,nb111_fizO(%esp)

        ## update iH forces 
        movq nb111_fixH(%esp),%mm2
        movq nb111_fiyH(%esp),%mm3
        movq nb111_fizH(%esp),%mm4
        pfadd %mm5,%mm2
        pfadd %mm6,%mm3
        pfadd %mm7,%mm4
        movq %mm2,nb111_fixH(%esp)
        movq %mm3,nb111_fiyH(%esp)
        movq %mm4,nb111_fizH(%esp)

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
        decl nb111_innerk(%esp)
        jz  _nb_kernel111_ia32_3dnow.nb111_updateouterdata
        jmp _nb_kernel111_ia32_3dnow.nb111_inner_loop
_nb_kernel111_ia32_3dnow.nb111_updateouterdata: 
        movl  nb111_ii3(%esp),%ecx

        movq  (%edi,%ecx,4),%mm6       ## increment iO force  
        movd  8(%edi,%ecx,4),%mm7
        pfadd nb111_fixO(%esp),%mm6
        pfadd nb111_fizO(%esp),%mm7
        movq  %mm6,(%edi,%ecx,4)
        movd  %mm7,8(%edi,%ecx,4)

        movq  nb111_fixH(%esp),%mm0
        movq  nb111_fiyH(%esp),%mm3
        movq  nb111_fizH(%esp),%mm1
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


        movl  nb111_fshift(%ebp),%ebx      ## increment fshift force 
        movl  nb111_is3(%esp),%edx

        movq  (%ebx,%edx,4),%mm6
        movd  8(%ebx,%edx,4),%mm7
        pfadd nb111_fixO(%esp),%mm6
        pfadd nb111_fizO(%esp),%mm7
        pfadd %mm0,%mm6
        pfadd %mm1,%mm7
        pfadd %mm2,%mm6
        pfadd %mm3,%mm7
        movq  %mm6,(%ebx,%edx,4)
        movd  %mm7,8(%ebx,%edx,4)

        ## get n from stack
        movl nb111_n(%esp),%esi
        ## get group index for i particle 
        movl  nb111_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        movq  nb111_vctot(%esp),%mm7
        pfacc %mm7,%mm7           ## get and sum the two parts of total potential 

        movl  nb111_Vc(%ebp),%eax
        movd  (%eax,%edx,4),%mm6
        pfadd %mm7,%mm6
        movd  %mm6,(%eax,%edx,4)          ## increment vc[gid] 

        movq  nb111_Vvdwtot(%esp),%mm7
        pfacc %mm7,%mm7           ## same for Vvdw 

        movl  nb111_Vvdw(%ebp),%eax
        movd  (%eax,%edx,4),%mm6
        pfadd %mm7,%mm6
        movd  %mm6,(%eax,%edx,4)          ## increment Vvdw[gid] 
       ## finish if last 
        movl nb111_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jecxz _nb_kernel111_ia32_3dnow.nb111_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb111_n(%esp)
        jmp _nb_kernel111_ia32_3dnow.nb111_outer
_nb_kernel111_ia32_3dnow.nb111_outerend: 
        ## check if more outer neighborlists remain
        movl  nb111_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jecxz _nb_kernel111_ia32_3dnow.nb111_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel111_ia32_3dnow.nb111_threadloop
_nb_kernel111_ia32_3dnow.nb111_end: 
        femms

        movl nb111_nouter(%esp),%eax
        movl nb111_ninner(%esp),%ebx
        movl nb111_outeriter(%ebp),%ecx
        movl nb111_inneriter(%ebp),%edx
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




.globl nb_kernel111nf_ia32_3dnow
.globl _nb_kernel111nf_ia32_3dnow
nb_kernel111nf_ia32_3dnow:      
_nb_kernel111nf_ia32_3dnow:     
.set nb111nf_p_nri, 8
.set nb111nf_iinr, 12
.set nb111nf_jindex, 16
.set nb111nf_jjnr, 20
.set nb111nf_shift, 24
.set nb111nf_shiftvec, 28
.set nb111nf_fshift, 32
.set nb111nf_gid, 36
.set nb111nf_pos, 40
.set nb111nf_faction, 44
.set nb111nf_charge, 48
.set nb111nf_p_facel, 52
.set nb111nf_p_krf, 56
.set nb111nf_p_crf, 60
.set nb111nf_Vc, 64
.set nb111nf_type, 68
.set nb111nf_p_ntype, 72
.set nb111nf_vdwparam, 76
.set nb111nf_Vvdw, 80
.set nb111nf_p_tabscale, 84
.set nb111nf_VFtab, 88
.set nb111nf_invsqrta, 92
.set nb111nf_dvda, 96
.set nb111nf_p_gbtabscale, 100
.set nb111nf_GBtab, 104
.set nb111nf_p_nthreads, 108
.set nb111nf_count, 112
.set nb111nf_mtx, 116
.set nb111nf_outeriter, 120
.set nb111nf_inneriter, 124
.set nb111nf_work, 128
                        ## stack offsets for local variables 
.set nb111nf_is3, 0
.set nb111nf_ii3, 4
.set nb111nf_ixO, 8
.set nb111nf_iyO, 12
.set nb111nf_izO, 16
.set nb111nf_ixH, 20
.set nb111nf_iyH, 28
.set nb111nf_izH, 36
.set nb111nf_iqO, 44
.set nb111nf_iqH, 52
.set nb111nf_vctot, 60
.set nb111nf_Vvdwtot, 68
.set nb111nf_c6, 76
.set nb111nf_c12, 84
.set nb111nf_ntia, 92
.set nb111nf_innerjjnr, 96
.set nb111nf_innerk, 100
.set nb111nf_n, 104                         ## idx for outer loop
.set nb111nf_nn1, 108                       ## number of outer iterations
.set nb111nf_nri, 112
.set nb111nf_ntype, 116
.set nb111nf_nouter, 120
.set nb111nf_ninner, 124
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
        ## zero iteration counters
        movl nb111nf_p_nri(%ebp),%ecx
        movl nb111nf_p_ntype(%ebp),%edx
        movl nb111nf_p_facel(%ebp),%esi
        movl (%ecx),%ecx
        movl (%edx),%edx
        movl %ecx,nb111nf_nri(%esp)
        movl %edx,nb111nf_ntype(%esp)

        movl $0,%eax
        movl %eax,nb111nf_nouter(%esp)
        movl %eax,nb111nf_ninner(%esp)

        ## assume we have at least one i particle - start directly      

        movl  nb111nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx),%ebx           ## ebx=ii 

        movl  nb111nf_charge(%ebp),%edx
        movd  (%esi),%mm1
        movd  (%edx,%ebx,4),%mm2    ## mm2=charge[ii0] 
        pfmul %mm1,%mm2
        movq  %mm2,nb111nf_iqO(%esp)        ## iqO = facel*charge[ii] 

        movd  4(%edx,%ebx,4),%mm2       ## mm2=charge[ii0+1] 
        pfmul %mm1,%mm2
        punpckldq %mm2,%mm2         ## spread to both halves 
        movq  %mm2,nb111nf_iqH(%esp)        ## iqH = facel*charge[ii0+1] 

        movl  nb111nf_type(%ebp),%edx
        movl  (%edx,%ebx,4),%ecx
        shll  %ecx
        imull nb111nf_ntype(%esp),%ecx        ## ecx = ntia = 2*ntype*type[ii0]  
        movl  %ecx,nb111nf_ntia(%esp)

_nb_kernel111nf_ia32_3dnow.nb111nf_threadloop: 
        movl  nb111nf_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel111nf_ia32_3dnow.nb111nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel111nf_ia32_3dnow.nb111nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb111nf_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb111nf_n(%esp)
        movl %ebx,nb111nf_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel111nf_ia32_3dnow.nb111nf_outerstart
        jmp _nb_kernel111nf_ia32_3dnow.nb111nf_end

_nb_kernel111nf_ia32_3dnow.nb111nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb111nf_nouter(%esp),%ebx
        movl %ebx,nb111nf_nouter(%esp)

_nb_kernel111nf_ia32_3dnow.nb111nf_outer: 
        movl  nb111nf_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb111nf_is3(%esp)            ## store is3 

        movl  nb111nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movq  (%eax,%ebx,4),%mm5        ## move shX/shY to mm5 and shZ to mm6. 
        movd  8(%eax,%ebx,4),%mm6
        movq  %mm5,%mm0
        movq  %mm5,%mm1
        movq  %mm6,%mm2
        punpckldq %mm0,%mm0         ## also expand shX,Y,Z in mm0--mm2. 
        punpckhdq %mm1,%mm1
        punpckldq %mm2,%mm2

        movl  nb111nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx,%esi,4),%ebx    ## ebx=ii 

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb111nf_pos(%ebp),%eax      ## eax = base of pos[] 

        pfadd (%eax,%ebx,4),%mm5    ## ix = shX + posX (and iy too) 
        movd  8(%eax,%ebx,4),%mm7       ## cant use direct memory add for 4 bytes (iz) 
        movl  %ebx,nb111nf_ii3(%esp)        ## (use mm7 as temp. storage for iz.) 
        pfadd %mm7,%mm6
        movq  %mm5,nb111nf_ixO(%esp)
        movq  %mm6,nb111nf_izO(%esp)

        movd  12(%eax,%ebx,4),%mm3
        movd  16(%eax,%ebx,4),%mm4
        movd  20(%eax,%ebx,4),%mm5
        punpckldq  24(%eax,%ebx,4),%mm3
        punpckldq  28(%eax,%ebx,4),%mm4
        punpckldq  32(%eax,%ebx,4),%mm5    ## coords of H1 in low mm3-mm5, H2 in high 

        pfadd %mm3,%mm0
        pfadd %mm4,%mm1
        pfadd %mm5,%mm2
        movq %mm0,nb111nf_ixH(%esp)
        movq %mm1,nb111nf_iyH(%esp)
        movq %mm2,nb111nf_izH(%esp)

        ## clear vctot and i forces 
        pxor  %mm7,%mm7
        movq  %mm7,nb111nf_vctot(%esp)
        movq  %mm7,nb111nf_Vvdwtot(%esp)

        movl  nb111nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 
        movl  %edx,nb111nf_innerk(%esp)      ## number of innerloop atoms 
        addl  nb111nf_ninner(%esp),%edx
        movl  %edx,nb111nf_ninner(%esp)

        movl  nb111nf_pos(%ebp),%esi
        movl  nb111nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb111nf_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
_nb_kernel111nf_ia32_3dnow.nb111nf_inner_loop: 
        ## a single j particle iteration here - compare with the unrolled code for comments. 
        movl  nb111nf_innerjjnr(%esp),%eax
        movl  (%eax),%eax       ## eax=jnr offset 
        addl $4,nb111nf_innerjjnr(%esp)             ## advance pointer 
        prefetch 16(%ecx)          ## prefetch data - trial and error says 16 is best 

        movl nb111nf_charge(%ebp),%ecx
        movd (%ecx,%eax,4),%mm7
        punpckldq %mm7,%mm7
        movq %mm7,%mm6
        pfmul nb111nf_iqO(%esp),%mm6
        pfmul nb111nf_iqH(%esp),%mm7    ## mm6=qqO, mm7=qqH 

        movl nb111nf_type(%ebp),%ecx
        movl (%ecx,%eax,4),%edx          ## type [jnr] 
        movl nb111nf_vdwparam(%ebp),%ecx
        shll %edx
        addl nb111nf_ntia(%esp),%edx         ## tja = ntia + 2*type 
        movd (%ecx,%edx,4),%mm5         ## mm5 = 1st c6                 
        movq %mm5,nb111nf_c6(%esp)
        movd 4(%ecx,%edx,4),%mm5        ## mm5 = 1st c12                
        movq %mm5,nb111nf_c12(%esp)

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

        pfsubr nb111nf_ixO(%esp),%mm0
        pfsubr nb111nf_izO(%esp),%mm1

        pfmul %mm0,%mm0
        pfmul %mm1,%mm1
        pfacc %mm1,%mm0
        pfadd %mm1,%mm0         ## mm0=rsqO 

        punpckldq %mm2,%mm2
        punpckldq %mm3,%mm3
        punpckldq %mm4,%mm4 ## mm2-mm4 is jx-jz 
        pfsubr nb111nf_ixH(%esp),%mm2
        pfsubr nb111nf_iyH(%esp),%mm3
        pfsubr nb111nf_izH(%esp),%mm4   ## mm2-mm4 is dxH-dzH 

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

        movq  %mm4,%mm0
        pfmul %mm4,%mm0
        pfmul %mm4,%mm0         ## mm0=rinvsix 
        movq  %mm0,%mm2
        pfmul %mm2,%mm2         ## mm2=rintwelve 

        ## calculate potential and scalar force 
        pfmul %mm1,%mm6         ## mm6=vcoul 
        movq  %mm6,%mm1         ## use mm1 for fscal sum 

        ## LJ for the oxygen 
        pfmul nb111nf_c6(%esp),%mm0
        pfmul nb111nf_c12(%esp),%mm2

        ## calc nb potential 
        pfsub %mm0,%mm2
        ## update nb potential 
        pfadd nb111nf_Vvdwtot(%esp),%mm2
        movq %mm2,nb111nf_Vvdwtot(%esp)

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
        pfadd nb111nf_vctot(%esp),%mm7
        movq %mm7,nb111nf_vctot(%esp)

        ##  done  - one more? 
        decl nb111nf_innerk(%esp)
        jz  _nb_kernel111nf_ia32_3dnow.nb111nf_updateouterdata
        jmp _nb_kernel111nf_ia32_3dnow.nb111nf_inner_loop
_nb_kernel111nf_ia32_3dnow.nb111nf_updateouterdata: 
        ## get n from stack
        movl nb111nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb111nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        movq  nb111nf_vctot(%esp),%mm7
        pfacc %mm7,%mm7           ## get and sum the two parts of total potential 

        movl  nb111nf_Vc(%ebp),%eax
        movd  (%eax,%edx,4),%mm6
        pfadd %mm7,%mm6
        movd  %mm6,(%eax,%edx,4)          ## increment vc[gid] 

        movq  nb111nf_Vvdwtot(%esp),%mm7
        pfacc %mm7,%mm7           ## same for Vvdw 

        movl  nb111nf_Vvdw(%ebp),%eax
        movd  (%eax,%edx,4),%mm6
        pfadd %mm7,%mm6
        movd  %mm6,(%eax,%edx,4)          ## increment Vvdw[gid] 
        ## finish if last 
        movl nb111nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jecxz _nb_kernel111nf_ia32_3dnow.nb111nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb111nf_n(%esp)
        jmp _nb_kernel111nf_ia32_3dnow.nb111nf_outer
_nb_kernel111nf_ia32_3dnow.nb111nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb111nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jecxz _nb_kernel111nf_ia32_3dnow.nb111nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel111nf_ia32_3dnow.nb111nf_threadloop
_nb_kernel111nf_ia32_3dnow.nb111nf_end: 
        femms

        movl nb111nf_nouter(%esp),%eax
        movl nb111nf_ninner(%esp),%ebx
        movl nb111nf_outeriter(%ebp),%ecx
        movl nb111nf_inneriter(%ebp),%edx
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



