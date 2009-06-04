##
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




.globl nb_kernel030_ia32_sse2
.globl _nb_kernel030_ia32_sse2
nb_kernel030_ia32_sse2: 
_nb_kernel030_ia32_sse2:        
.set nb030_p_nri, 8
.set nb030_iinr, 12
.set nb030_jindex, 16
.set nb030_jjnr, 20
.set nb030_shift, 24
.set nb030_shiftvec, 28
.set nb030_fshift, 32
.set nb030_gid, 36
.set nb030_pos, 40
.set nb030_faction, 44
.set nb030_charge, 48
.set nb030_p_facel, 52
.set nb030_argkrf, 56
.set nb030_argcrf, 60
.set nb030_Vc, 64
.set nb030_type, 68
.set nb030_p_ntype, 72
.set nb030_vdwparam, 76
.set nb030_Vvdw, 80
.set nb030_p_tabscale, 84
.set nb030_VFtab, 88
.set nb030_invsqrta, 92
.set nb030_dvda, 96
.set nb030_p_gbtabscale, 100
.set nb030_GBtab, 104
.set nb030_p_nthreads, 108
.set nb030_count, 112
.set nb030_mtx, 116
.set nb030_outeriter, 120
.set nb030_inneriter, 124
.set nb030_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse2 use 
.set nb030_ix, 0
.set nb030_iy, 16
.set nb030_iz, 32
.set nb030_dx, 48
.set nb030_dy, 64
.set nb030_dz, 80
.set nb030_two, 96
.set nb030_tsc, 112
.set nb030_c6, 128
.set nb030_c12, 144
.set nb030_fscal, 160
.set nb030_Vvdwtot, 176
.set nb030_fix, 192
.set nb030_fiy, 208
.set nb030_fiz, 224
.set nb030_half, 240
.set nb030_three, 256
.set nb030_is3, 272
.set nb030_ii3, 276
.set nb030_ntia, 280
.set nb030_innerjjnr, 284
.set nb030_innerk, 288
.set nb030_n, 292
.set nb030_nn1, 296
.set nb030_nri, 300
.set nb030_ntype, 304
.set nb030_nouter, 308
.set nb030_ninner, 312
.set nb030_salign, 316
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $320,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb030_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb030_p_nri(%ebp),%ecx
        movl nb030_p_ntype(%ebp),%edi
        movl (%ecx),%ecx
        movl (%edi),%edi
        movl %ecx,nb030_nri(%esp)
        movl %edi,nb030_ntype(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb030_nouter(%esp)
        movl %eax,nb030_ninner(%esp)


        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double 0.5 IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb030_half(%esp)
        movl %ebx,nb030_half+4(%esp)
        movsd nb030_half(%esp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## 1.0
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## 2.0
        addpd  %xmm2,%xmm3      ## 3.0
        movapd %xmm1,nb030_half(%esp)
        movapd %xmm2,nb030_two(%esp)
        movapd %xmm3,nb030_three(%esp)
        movl nb030_p_tabscale(%ebp),%eax
        movsd (%eax),%xmm3
        shufpd $0,%xmm3,%xmm3
        movapd %xmm3,nb030_tsc(%esp)

_nb_kernel030_ia32_sse2.nb030_threadloop: 
        movl  nb030_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel030_ia32_sse2.nb030_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel030_ia32_sse2.nb030_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb030_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb030_n(%esp)
        movl %ebx,nb030_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel030_ia32_sse2.nb030_outerstart
        jmp _nb_kernel030_ia32_sse2.nb030_end

_nb_kernel030_ia32_sse2.nb030_outerstart: 
        ## ebx contains number of outer iterations
        addl nb030_nouter(%esp),%ebx
        movl %ebx,nb030_nouter(%esp)

_nb_kernel030_ia32_sse2.nb030_outer: 
        movl  nb030_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb030_is3(%esp)      ## store is3 

        movl  nb030_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movsd (%eax,%ebx,8),%xmm0
        movsd 8(%eax,%ebx,8),%xmm1
        movsd 16(%eax,%ebx,8),%xmm2

        movl  nb030_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx,%esi,4),%ebx            ## ebx =ii 

        movl  nb030_type(%ebp),%edx
        movl  (%edx,%ebx,4),%edx
        imull nb030_ntype(%esp),%edx
        shll  %edx
        movl  %edx,nb030_ntia(%esp)

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb030_pos(%ebp),%eax      ## eax = base of pos[]  

        addsd (%eax,%ebx,8),%xmm0
        addsd 8(%eax,%ebx,8),%xmm1
        addsd 16(%eax,%ebx,8),%xmm2

        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2

        movapd %xmm0,nb030_ix(%esp)
        movapd %xmm1,nb030_iy(%esp)
        movapd %xmm2,nb030_iz(%esp)

        movl  %ebx,nb030_ii3(%esp)

        ## clear tot potential and i forces 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb030_Vvdwtot(%esp)
        movapd %xmm4,nb030_fix(%esp)
        movapd %xmm4,nb030_fiy(%esp)
        movapd %xmm4,nb030_fiz(%esp)

        movl  nb030_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb030_pos(%ebp),%esi
        movl  nb030_faction(%ebp),%edi
        movl  nb030_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb030_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb030_ninner(%esp),%ecx
        movl  %ecx,nb030_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb030_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel030_ia32_sse2.nb030_unroll_loop
        jmp   _nb_kernel030_ia32_sse2.nb030_checksingle
_nb_kernel030_ia32_sse2.nb030_unroll_loop: 
        ## twice unrolled innerloop here 
        movl  nb030_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx
        addl $8,nb030_innerjjnr(%esp)             ## advance pointer (unrolled 2) 

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1

        movl nb030_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl (%esi,%ebx,4),%ebx
        movl nb030_vdwparam(%ebp),%esi
        shll %eax
        shll %ebx
        movl nb030_ntia(%esp),%edi
        addl %edi,%eax
        addl %edi,%ebx

        movlpd (%esi,%eax,8),%xmm6      ## c6a
        movlpd (%esi,%ebx,8),%xmm7      ## c6b
        movhpd 8(%esi,%eax,8),%xmm6     ## c6a c12a 
        movhpd 8(%esi,%ebx,8),%xmm7     ## c6b c12b 

        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movd  %mm0,%eax
        movd  %mm1,%ebx

        movapd %xmm4,nb030_c6(%esp)
        movapd %xmm6,nb030_c12(%esp)

        movl nb030_pos(%ebp),%esi               ## base of pos[] 
        leal  (%eax,%eax,2),%eax        ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx

        ## move two coordinates to xmm0-xmm2    
        movlpd (%esi,%eax,8),%xmm0
        movlpd 8(%esi,%eax,8),%xmm1
        movlpd 16(%esi,%eax,8),%xmm2
        movhpd (%esi,%ebx,8),%xmm0
        movhpd 8(%esi,%ebx,8),%xmm1
        movhpd 16(%esi,%ebx,8),%xmm2

        ## move nb030_ix-iz to xmm4-xmm6 
        movapd nb030_ix(%esp),%xmm4
        movapd nb030_iy(%esp),%xmm5
        movapd nb030_iz(%esp),%xmm6

        ## calc dr 
        subpd %xmm0,%xmm4
        subpd %xmm1,%xmm5
        subpd %xmm2,%xmm6

        ## store dr 
        movapd %xmm4,nb030_dx(%esp)
        movapd %xmm5,nb030_dy(%esp)
        movapd %xmm6,nb030_dz(%esp)
        ## square it 
        mulpd %xmm4,%xmm4
        mulpd %xmm5,%xmm5
        mulpd %xmm6,%xmm6
        addpd %xmm5,%xmm4
        addpd %xmm6,%xmm4
        ## rsq in xmm4 

        cvtpd2ps %xmm4,%xmm5
        rsqrtps %xmm5,%xmm5
        cvtps2pd %xmm5,%xmm2    ## lu in low xmm2 

        ## lookup seed in xmm2 
        movapd %xmm2,%xmm5      ## copy of lu 
        mulpd %xmm2,%xmm2       ## lu*lu 
        movapd nb030_three(%esp),%xmm1
        mulpd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb030_half(%esp),%xmm0
        subpd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm1
        mulpd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulpd %xmm1,%xmm1       ## lu*lu 
        movapd nb030_three(%esp),%xmm2
        mulpd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb030_half(%esp),%xmm0
        subpd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm2
        mulpd %xmm2,%xmm0       ## xmm0=iter2 of rinv (new lu) 

        mulpd %xmm0,%xmm4       ## xmm4=r 
        mulpd nb030_tsc(%esp),%xmm4

        cvttpd2pi %xmm4,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm5
        subpd %xmm5,%xmm4
        movapd %xmm4,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $3,%mm6           ## idx *= 8 

        movd %eax,%mm0
        movd %ebx,%mm1

        movl nb030_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx

        ## dispersion 
        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movlpd (%esi,%ebx,8),%xmm3      ## Y2 
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        
        movhpd 8(%esi,%ebx,8),%xmm3     ## Y2 F2 
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 16(%esi,%eax,8),%xmm6    ## G1
        movlpd 16(%esi,%ebx,8),%xmm3    ## G2
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        
        movhpd 24(%esi,%ebx,8),%xmm3    ## G2 H2 
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## dispersion table ready, in xmm4-xmm7         
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        mulpd  nb030_two(%esp),%xmm7    ## two*Heps2 
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 

        movapd nb030_c6(%esp),%xmm4
        mulpd  %xmm4,%xmm7       ## fijD 
        mulpd  %xmm4,%xmm5       ## Vvdw6 

        ## put scalar force on stack Update Vvdwtot directly 
        addpd  nb030_Vvdwtot(%esp),%xmm5
        movapd %xmm7,nb030_fscal(%esp)
        movapd %xmm5,nb030_Vvdwtot(%esp)

        ## repulsion 
        movlpd 32(%esi,%eax,8),%xmm4    ## Y1   
        movlpd 32(%esi,%ebx,8),%xmm3    ## Y2 
        movhpd 40(%esi,%eax,8),%xmm4    ## Y1 F1        
        movhpd 40(%esi,%ebx,8),%xmm3    ## Y2 F2 

        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 48(%esi,%eax,8),%xmm6    ## G1
        movlpd 48(%esi,%ebx,8),%xmm3    ## G2
        movhpd 56(%esi,%eax,8),%xmm6    ## G1 H1        
        movhpd 56(%esi,%ebx,8),%xmm3    ## G2 H2 

        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 

        ## table ready, in xmm4-xmm7    
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        mulpd  nb030_two(%esp),%xmm7    ## two*Heps2 
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 

        movapd nb030_c12(%esp),%xmm4
        mulpd  %xmm4,%xmm7
        mulpd  %xmm4,%xmm5
        addpd  nb030_fscal(%esp),%xmm7

        addpd  nb030_Vvdwtot(%esp),%xmm5
        movapd %xmm5,nb030_Vvdwtot(%esp)
        xorpd  %xmm4,%xmm4

        mulpd nb030_tsc(%esp),%xmm7
        mulpd %xmm0,%xmm7
        subpd  %xmm7,%xmm4

        movapd nb030_dx(%esp),%xmm0
        movapd nb030_dy(%esp),%xmm1
        movapd nb030_dz(%esp),%xmm2

        movd %mm0,%eax
        movd %mm1,%ebx

        movl   nb030_faction(%ebp),%edi
        mulpd  %xmm4,%xmm0
        mulpd  %xmm4,%xmm1
        mulpd  %xmm4,%xmm2
        ## xmm0-xmm2 contains tx-tz (partial force) 
        ## now update f_i 
        movapd nb030_fix(%esp),%xmm3
        movapd nb030_fiy(%esp),%xmm4
        movapd nb030_fiz(%esp),%xmm5
        addpd  %xmm0,%xmm3
        addpd  %xmm1,%xmm4
        addpd  %xmm2,%xmm5
        movapd %xmm3,nb030_fix(%esp)
        movapd %xmm4,nb030_fiy(%esp)
        movapd %xmm5,nb030_fiz(%esp)
        ## the fj's - start by accumulating forces from memory 
        movlpd (%edi,%eax,8),%xmm3
        movlpd 8(%edi,%eax,8),%xmm4
        movlpd 16(%edi,%eax,8),%xmm5
        movhpd (%edi,%ebx,8),%xmm3
        movhpd 8(%edi,%ebx,8),%xmm4
        movhpd 16(%edi,%ebx,8),%xmm5
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        movlpd %xmm3,(%edi,%eax,8)
        movlpd %xmm4,8(%edi,%eax,8)
        movlpd %xmm5,16(%edi,%eax,8)
        movhpd %xmm3,(%edi,%ebx,8)
        movhpd %xmm4,8(%edi,%ebx,8)
        movhpd %xmm5,16(%edi,%ebx,8)

        ## should we do one more iteration? 
        subl $2,nb030_innerk(%esp)
        jl    _nb_kernel030_ia32_sse2.nb030_checksingle
        jmp   _nb_kernel030_ia32_sse2.nb030_unroll_loop

_nb_kernel030_ia32_sse2.nb030_checksingle:      
        movl  nb030_innerk(%esp),%edx
        andl  $1,%edx
        jnz    _nb_kernel030_ia32_sse2.nb030_dosingle
        jmp    _nb_kernel030_ia32_sse2.nb030_updateouterdata
_nb_kernel030_ia32_sse2.nb030_dosingle: 
        movl  nb030_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax

        movd  %eax,%mm0         ## use mmx registers as temp storage 

        movl nb030_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl nb030_vdwparam(%ebp),%esi
        shll %eax
        movl nb030_ntia(%esp),%edi
        addl %edi,%eax

        movlpd (%esi,%eax,8),%xmm6      ## c6a
        movhpd 8(%esi,%eax,8),%xmm6     ## c6a c12a 

        xorpd %xmm7,%xmm7
        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movd  %mm0,%eax

        movapd %xmm4,nb030_c6(%esp)
        movapd %xmm6,nb030_c12(%esp)

        movl nb030_pos(%ebp),%esi               ## base of pos[] 
        leal  (%eax,%eax,2),%eax        ## replace jnr with j3 

        ## move coordinates to xmm0-xmm2        
        movlpd (%esi,%eax,8),%xmm0
        movlpd 8(%esi,%eax,8),%xmm1
        movlpd 16(%esi,%eax,8),%xmm2

        ## move nb030_ix-iz to xmm4-xmm6 
        movapd nb030_ix(%esp),%xmm4
        movapd nb030_iy(%esp),%xmm5
        movapd nb030_iz(%esp),%xmm6

        ## calc dr 
        subsd %xmm0,%xmm4
        subsd %xmm1,%xmm5
        subsd %xmm2,%xmm6

        ## store dr 
        movapd %xmm4,nb030_dx(%esp)
        movapd %xmm5,nb030_dy(%esp)
        movapd %xmm6,nb030_dz(%esp)
        ## square it 
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        mulsd %xmm6,%xmm6
        addsd %xmm5,%xmm4
        addsd %xmm6,%xmm4
        ## rsq in xmm4 

        cvtsd2ss %xmm4,%xmm5
        rsqrtss %xmm5,%xmm5
        cvtss2sd %xmm5,%xmm2    ## lu in low xmm2 

        ## lookup seed in xmm2 
        movapd %xmm2,%xmm5      ## copy of lu 
        mulsd %xmm2,%xmm2       ## lu*lu 
        movapd nb030_three(%esp),%xmm1
        mulsd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb030_half(%esp),%xmm0
        subsd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm1
        mulsd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulsd %xmm1,%xmm1       ## lu*lu 
        movapd nb030_three(%esp),%xmm2
        mulsd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb030_half(%esp),%xmm0
        subsd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm2
        mulsd %xmm2,%xmm0       ## xmm0=iter2 of rinv (new lu) 

        mulsd %xmm0,%xmm4       ## xmm4=r 
        mulsd nb030_tsc(%esp),%xmm4

        movd %eax,%mm0

        cvttsd2si %xmm4,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm5
        subsd %xmm5,%xmm4
        movapd %xmm4,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 
        shll $3,%eax

        movl nb030_VFtab(%ebp),%esi

        ## dispersion 
        movlpd (%esi,%eax,8),%xmm4      ## Y1
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1 
        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1  
        unpckhpd %xmm3,%xmm5    ## F1  

        movlpd 16(%esi,%eax,8),%xmm6    ## G1   
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        

        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1  
        unpckhpd %xmm3,%xmm7    ## H1  
        ## dispersion table ready, in xmm4-xmm7         
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb030_two(%esp),%xmm7    ## two*Heps2 
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 

        movapd nb030_c6(%esp),%xmm4
        mulsd  %xmm4,%xmm7       ## fijD 
        mulsd  %xmm4,%xmm5       ## Vvdw6 

        ## put scalar force on stack Update Vvdwtot directly 
        addsd  nb030_Vvdwtot(%esp),%xmm5
        movlpd %xmm7,nb030_fscal(%esp)
        movlpd %xmm5,nb030_Vvdwtot(%esp)

        ## repulsion 
        movlpd 32(%esi,%eax,8),%xmm4    ## Y1
        movhpd 40(%esi,%eax,8),%xmm4    ## Y1 F1 

        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1  
        unpckhpd %xmm3,%xmm5    ## F1  

        movlpd 48(%esi,%eax,8),%xmm6    ## G1   
        movhpd 56(%esi,%eax,8),%xmm6    ## G1 H1        
        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1  
        unpckhpd %xmm3,%xmm7    ## H1  

        ## table ready, in xmm4-xmm7    
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb030_two(%esp),%xmm7    ## two*Heps2 
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 

        movapd nb030_c12(%esp),%xmm4
        mulsd  %xmm4,%xmm7
        mulsd  %xmm4,%xmm5
        addsd  nb030_fscal(%esp),%xmm7

        addsd  nb030_Vvdwtot(%esp),%xmm5
        movlpd %xmm5,nb030_Vvdwtot(%esp)
        xorpd  %xmm4,%xmm4

        mulsd nb030_tsc(%esp),%xmm7
        mulsd %xmm0,%xmm7
        subsd  %xmm7,%xmm4

        movapd nb030_dx(%esp),%xmm0
        movapd nb030_dy(%esp),%xmm1
        movapd nb030_dz(%esp),%xmm2

        movd %mm0,%eax

        movl   nb030_faction(%ebp),%edi
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2
        ## xmm0-xmm2 contains tx-tz (partial force) 
        ## now update f_i 
        movapd nb030_fix(%esp),%xmm3
        movapd nb030_fiy(%esp),%xmm4
        movapd nb030_fiz(%esp),%xmm5
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm5
        movlpd %xmm3,nb030_fix(%esp)
        movlpd %xmm4,nb030_fiy(%esp)
        movlpd %xmm5,nb030_fiz(%esp)
        ## the fj's - start by accumulating forces from memory 
        movlpd (%edi,%eax,8),%xmm3
        movlpd 8(%edi,%eax,8),%xmm4
        movlpd 16(%edi,%eax,8),%xmm5
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        movlpd %xmm3,(%edi,%eax,8)
        movlpd %xmm4,8(%edi,%eax,8)
        movlpd %xmm5,16(%edi,%eax,8)
_nb_kernel030_ia32_sse2.nb030_updateouterdata: 
        movl  nb030_ii3(%esp),%ecx
        movl  nb030_faction(%ebp),%edi
        movl  nb030_fshift(%ebp),%esi
        movl  nb030_is3(%esp),%edx

        ## accumulate i forces in xmm0, xmm1, xmm2 
        movapd nb030_fix(%esp),%xmm0
        movapd nb030_fiy(%esp),%xmm1
        movapd nb030_fiz(%esp),%xmm2

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addsd  %xmm3,%xmm0
        addsd  %xmm4,%xmm1
        addsd  %xmm5,%xmm2 ## sum is in low xmm0-xmm2 

        ## increment i force 
        movsd  (%edi,%ecx,8),%xmm3
        movsd  8(%edi,%ecx,8),%xmm4
        movsd  16(%edi,%ecx,8),%xmm5
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm5
        movsd  %xmm3,(%edi,%ecx,8)
        movsd  %xmm4,8(%edi,%ecx,8)
        movsd  %xmm5,16(%edi,%ecx,8)

        ## increment fshift force  
        movsd  (%esi,%edx,8),%xmm3
        movsd  8(%esi,%edx,8),%xmm4
        movsd  16(%esi,%edx,8),%xmm5
        addsd %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm5
        movsd  %xmm3,(%esi,%edx,8)
        movsd  %xmm4,8(%esi,%edx,8)
        movsd  %xmm5,16(%esi,%edx,8)

        ## get n from stack
        movl nb030_n(%esp),%esi
        ## get group index for i particle 
        movl  nb030_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total lj energy and update it 
        movapd nb030_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb030_Vvdw(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## finish if last 
        movl nb030_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel030_ia32_sse2.nb030_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb030_n(%esp)
        jmp _nb_kernel030_ia32_sse2.nb030_outer
_nb_kernel030_ia32_sse2.nb030_outerend: 
        ## check if more outer neighborlists remain
        movl  nb030_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel030_ia32_sse2.nb030_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel030_ia32_sse2.nb030_threadloop
_nb_kernel030_ia32_sse2.nb030_end: 
        emms

        movl nb030_nouter(%esp),%eax
        movl nb030_ninner(%esp),%ebx
        movl nb030_outeriter(%ebp),%ecx
        movl nb030_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb030_salign(%esp),%eax
        addl %eax,%esp
        addl $320,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret



.globl nb_kernel030nf_ia32_sse2
.globl _nb_kernel030nf_ia32_sse2
nb_kernel030nf_ia32_sse2:       
_nb_kernel030nf_ia32_sse2:      
.set nb030nf_p_nri, 8
.set nb030nf_iinr, 12
.set nb030nf_jindex, 16
.set nb030nf_jjnr, 20
.set nb030nf_shift, 24
.set nb030nf_shiftvec, 28
.set nb030nf_fshift, 32
.set nb030nf_gid, 36
.set nb030nf_pos, 40
.set nb030nf_faction, 44
.set nb030nf_charge, 48
.set nb030nf_p_facel, 52
.set nb030nf_argkrf, 56
.set nb030nf_argcrf, 60
.set nb030nf_Vc, 64
.set nb030nf_type, 68
.set nb030nf_p_ntype, 72
.set nb030nf_vdwparam, 76
.set nb030nf_Vvdw, 80
.set nb030nf_p_tabscale, 84
.set nb030nf_VFtab, 88
.set nb030nf_invsqrta, 92
.set nb030nf_dvda, 96
.set nb030nf_p_gbtabscale, 100
.set nb030nf_GBtab, 104
.set nb030nf_p_nthreads, 108
.set nb030nf_count, 112
.set nb030nf_mtx, 116
.set nb030nf_outeriter, 120
.set nb030nf_inneriter, 124
.set nb030nf_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb030nf_ix, 0
.set nb030nf_iy, 16
.set nb030nf_iz, 32
.set nb030nf_tsc, 48
.set nb030nf_c6, 64
.set nb030nf_c12, 80
.set nb030nf_Vvdwtot, 96
.set nb030nf_half, 112
.set nb030nf_three, 128
.set nb030nf_is3, 144
.set nb030nf_ii3, 148
.set nb030nf_ntia, 152
.set nb030nf_innerjjnr, 156
.set nb030nf_innerk, 160
.set nb030nf_n, 164
.set nb030nf_nn1, 168
.set nb030nf_nri, 172
.set nb030nf_ntype, 176
.set nb030nf_nouter, 180
.set nb030nf_ninner, 184
.set nb030nf_salign, 188
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $192,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb030nf_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb030nf_p_nri(%ebp),%ecx
        movl nb030nf_p_ntype(%ebp),%edi
        movl (%ecx),%ecx
        movl (%edi),%edi
        movl %ecx,nb030nf_nri(%esp)
        movl %edi,nb030nf_ntype(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb030nf_nouter(%esp)
        movl %eax,nb030nf_ninner(%esp)


        ## create constant floating-point factors on stack
        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double 0.5 IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb030nf_half(%esp)
        movl %ebx,nb030nf_half+4(%esp)
        movsd nb030nf_half(%esp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## 1.0
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## 2.0
        addpd  %xmm2,%xmm3      ## 3.0
        movapd %xmm1,nb030nf_half(%esp)
        movapd %xmm3,nb030nf_three(%esp)

        movl nb030nf_p_tabscale(%ebp),%eax
        movsd (%eax),%xmm3
        shufpd $0,%xmm3,%xmm3
        movapd %xmm3,nb030nf_tsc(%esp)

_nb_kernel030nf_ia32_sse2.nb030nf_threadloop: 
        movl  nb030nf_count(%ebp),%esi          ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel030nf_ia32_sse2.nb030nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                           ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel030nf_ia32_sse2.nb030nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb030nf_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb030nf_n(%esp)
        movl %ebx,nb030nf_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel030nf_ia32_sse2.nb030nf_outerstart
        jmp _nb_kernel030nf_ia32_sse2.nb030nf_end

_nb_kernel030nf_ia32_sse2.nb030nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb030nf_nouter(%esp),%ebx
        movl %ebx,nb030nf_nouter(%esp)

_nb_kernel030nf_ia32_sse2.nb030nf_outer: 
        movl  nb030nf_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 

        movl  nb030nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movsd (%eax,%ebx,8),%xmm0
        movsd 8(%eax,%ebx,8),%xmm1
        movsd 16(%eax,%ebx,8),%xmm2

        movl  nb030nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx,%esi,4),%ebx    ## ebx =ii 

        movl  nb030nf_type(%ebp),%edx
        movl  (%edx,%ebx,4),%edx
        imull nb030nf_ntype(%esp),%edx
        shll  %edx
        movl  %edx,nb030nf_ntia(%esp)

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb030nf_pos(%ebp),%eax      ## eax = base of pos[]  

        addsd (%eax,%ebx,8),%xmm0
        addsd 8(%eax,%ebx,8),%xmm1
        addsd 16(%eax,%ebx,8),%xmm2

        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2

        movapd %xmm0,nb030nf_ix(%esp)
        movapd %xmm1,nb030nf_iy(%esp)
        movapd %xmm2,nb030nf_iz(%esp)

        movl  %ebx,nb030nf_ii3(%esp)

        ## clear tot potential 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb030nf_Vvdwtot(%esp)

        movl  nb030nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb030nf_pos(%ebp),%esi
        movl  nb030nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb030nf_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb030nf_ninner(%esp),%ecx
        movl  %ecx,nb030nf_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb030nf_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel030nf_ia32_sse2.nb030nf_unroll_loop
        jmp   _nb_kernel030nf_ia32_sse2.nb030nf_checksingle
_nb_kernel030nf_ia32_sse2.nb030nf_unroll_loop: 
        ## twice unrolled innerloop here 
        movl  nb030nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx
        addl $8,nb030nf_innerjjnr(%esp)             ## advance pointer (unrolled 2) 

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1

        movl nb030nf_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl (%esi,%ebx,4),%ebx
        movl nb030nf_vdwparam(%ebp),%esi
        shll %eax
        shll %ebx
        movl nb030nf_ntia(%esp),%edi
        addl %edi,%eax
        addl %edi,%ebx

        movlpd (%esi,%eax,8),%xmm6      ## c6a
        movlpd (%esi,%ebx,8),%xmm7      ## c6b
        movhpd 8(%esi,%eax,8),%xmm6     ## c6a c12a 
        movhpd 8(%esi,%ebx,8),%xmm7     ## c6b c12b 

        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movd  %mm0,%eax
        movd  %mm1,%ebx

        movapd %xmm4,nb030nf_c6(%esp)
        movapd %xmm6,nb030nf_c12(%esp)

        movl nb030nf_pos(%ebp),%esi             ## base of pos[] 
        leal  (%eax,%eax,2),%eax        ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx

        ## move two coordinates to xmm0-xmm2    
        movlpd (%esi,%eax,8),%xmm0
        movlpd 8(%esi,%eax,8),%xmm1
        movlpd 16(%esi,%eax,8),%xmm2
        movhpd (%esi,%ebx,8),%xmm0
        movhpd 8(%esi,%ebx,8),%xmm1
        movhpd 16(%esi,%ebx,8),%xmm2

        ## move nb030nf_ix-iz to xmm4-xmm6 
        movapd nb030nf_ix(%esp),%xmm4
        movapd nb030nf_iy(%esp),%xmm5
        movapd nb030nf_iz(%esp),%xmm6

        ## calc dr 
        subpd %xmm0,%xmm4
        subpd %xmm1,%xmm5
        subpd %xmm2,%xmm6

        ## square it 
        mulpd %xmm4,%xmm4
        mulpd %xmm5,%xmm5
        mulpd %xmm6,%xmm6
        addpd %xmm5,%xmm4
        addpd %xmm6,%xmm4
        ## rsq in xmm4 

        cvtpd2ps %xmm4,%xmm5
        rsqrtps %xmm5,%xmm5
        cvtps2pd %xmm5,%xmm2    ## lu in low xmm2 

        ## lookup seed in xmm2 
        movapd %xmm2,%xmm5      ## copy of lu 
        mulpd %xmm2,%xmm2       ## lu*lu 
        movapd nb030nf_three(%esp),%xmm1
        mulpd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb030nf_half(%esp),%xmm0
        subpd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm1
        mulpd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulpd %xmm1,%xmm1       ## lu*lu 
        movapd nb030nf_three(%esp),%xmm2
        mulpd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb030nf_half(%esp),%xmm0
        subpd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm2
        mulpd %xmm2,%xmm0       ## xmm0=iter2 of rinv (new lu) 

        mulpd %xmm0,%xmm4       ## xmm4=r 
        mulpd nb030nf_tsc(%esp),%xmm4

        cvttpd2pi %xmm4,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm5
        subpd %xmm5,%xmm4
        movapd %xmm4,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $3,%mm6           ## idx *= 8 

        movd %eax,%mm0
        movd %ebx,%mm1

        movl nb030nf_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx

        ## dispersion 
        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movlpd (%esi,%ebx,8),%xmm3      ## Y2 
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        
        movhpd 8(%esi,%ebx,8),%xmm3     ## Y2 F2 
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 16(%esi,%eax,8),%xmm6    ## G1
        movlpd 16(%esi,%ebx,8),%xmm3    ## G2
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        
        movhpd 24(%esi,%ebx,8),%xmm3    ## G2 H2 

        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## dispersion table ready, in xmm4-xmm7         
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 

        mulpd  nb030nf_c6(%esp),%xmm5   ## Vvdw6 

        ## Update Vvdwtot directly 
        addpd  nb030nf_Vvdwtot(%esp),%xmm5
        movapd %xmm5,nb030nf_Vvdwtot(%esp)

        ## repulsion 
        movlpd 32(%esi,%eax,8),%xmm4    ## Y1   
        movlpd 32(%esi,%ebx,8),%xmm3    ## Y2 
        movhpd 40(%esi,%eax,8),%xmm4    ## Y1 F1        
        movhpd 40(%esi,%ebx,8),%xmm3    ## Y2 F2 
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 48(%esi,%eax,8),%xmm6    ## G1
        movlpd 48(%esi,%ebx,8),%xmm3    ## G2
        movhpd 56(%esi,%eax,8),%xmm6    ## G1 H1        
        movhpd 56(%esi,%ebx,8),%xmm3    ## G2 H2 

        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 

        ## table ready, in xmm4-xmm7    
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 

        mulpd  nb030nf_c12(%esp),%xmm5

        addpd  nb030nf_Vvdwtot(%esp),%xmm5
        movapd %xmm5,nb030nf_Vvdwtot(%esp)

        ## should we do one more iteration? 
        subl $2,nb030nf_innerk(%esp)
        jl    _nb_kernel030nf_ia32_sse2.nb030nf_checksingle
        jmp   _nb_kernel030nf_ia32_sse2.nb030nf_unroll_loop

_nb_kernel030nf_ia32_sse2.nb030nf_checksingle:  
        movl  nb030nf_innerk(%esp),%edx
        andl  $1,%edx
        jnz    _nb_kernel030nf_ia32_sse2.nb030nf_dosingle
        jmp    _nb_kernel030nf_ia32_sse2.nb030nf_updateouterdata
_nb_kernel030nf_ia32_sse2.nb030nf_dosingle: 
        movl  nb030nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax

        movd  %eax,%mm0         ## use mmx registers as temp storage 

        movl nb030nf_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl nb030nf_vdwparam(%ebp),%esi
        shll %eax
        movl nb030nf_ntia(%esp),%edi
        addl %edi,%eax

        movlpd (%esi,%eax,8),%xmm6      ## c6a
        movhpd 8(%esi,%eax,8),%xmm6     ## c6a c12a 

        xorpd %xmm7,%xmm7
        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movd  %mm0,%eax

        movapd %xmm4,nb030nf_c6(%esp)
        movapd %xmm6,nb030nf_c12(%esp)

        movl nb030nf_pos(%ebp),%esi             ## base of pos[] 
        leal  (%eax,%eax,2),%eax        ## replace jnr with j3 

        ## move coordinates to xmm0-xmm2        
        movlpd (%esi,%eax,8),%xmm0
        movlpd 8(%esi,%eax,8),%xmm1
        movlpd 16(%esi,%eax,8),%xmm2

        ## move nb030nf_ix-iz to xmm4-xmm6 
        movapd nb030nf_ix(%esp),%xmm4
        movapd nb030nf_iy(%esp),%xmm5
        movapd nb030nf_iz(%esp),%xmm6

        ## calc dr 
        subsd %xmm0,%xmm4
        subsd %xmm1,%xmm5
        subsd %xmm2,%xmm6

        ## square it 
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        mulsd %xmm6,%xmm6
        addsd %xmm5,%xmm4
        addsd %xmm6,%xmm4
        ## rsq in xmm4 

        cvtsd2ss %xmm4,%xmm5
        rsqrtss %xmm5,%xmm5
        cvtss2sd %xmm5,%xmm2    ## lu in low xmm2 

        ## lookup seed in xmm2 
        movapd %xmm2,%xmm5      ## copy of lu 
        mulsd %xmm2,%xmm2       ## lu*lu 
        movapd nb030nf_three(%esp),%xmm1
        mulsd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb030nf_half(%esp),%xmm0
        subsd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm1
        mulsd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulsd %xmm1,%xmm1       ## lu*lu 
        movapd nb030nf_three(%esp),%xmm2
        mulsd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb030nf_half(%esp),%xmm0
        subsd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm2
        mulsd %xmm2,%xmm0       ## xmm0=iter2 of rinv (new lu) 

        mulsd %xmm0,%xmm4       ## xmm4=r 
        mulsd nb030nf_tsc(%esp),%xmm4

        movd %eax,%mm0

        cvttsd2si %xmm4,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm5
        subsd %xmm5,%xmm4
        movapd %xmm4,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 
        shll $3,%eax

        movl nb030nf_VFtab(%ebp),%esi

        ## dispersion 
        movlpd (%esi,%eax,8),%xmm4      ## Y1
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1 
        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1  
        unpckhpd %xmm3,%xmm5    ## F1  

        movlpd 16(%esi,%eax,8),%xmm6    ## G1
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1 
        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1  
        unpckhpd %xmm3,%xmm7    ## H1  
        ## dispersion table ready, in xmm4-xmm7         
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 

        mulsd  nb030nf_c6(%esp),%xmm5  ## Vvdw6 

        ## Update Vvdwtot directly 
        addsd  nb030nf_Vvdwtot(%esp),%xmm5
        movlpd %xmm5,nb030nf_Vvdwtot(%esp)

        ## repulsion 
        movlpd 32(%esi,%eax,8),%xmm4    ## Y1 
        movhpd 40(%esi,%eax,8),%xmm4    ## Y1 F1 
        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1  
        unpckhpd %xmm3,%xmm5    ## F1  

        movlpd 48(%esi,%eax,8),%xmm6    ## G1
        movhpd 56(%esi,%eax,8),%xmm6    ## G1 H1        
        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1  
        unpckhpd %xmm3,%xmm7    ## H1  

        ## table ready, in xmm4-xmm7    
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 

        mulsd  nb030nf_c12(%esp),%xmm5

        addsd  nb030nf_Vvdwtot(%esp),%xmm5
        movlpd %xmm5,nb030nf_Vvdwtot(%esp)

_nb_kernel030nf_ia32_sse2.nb030nf_updateouterdata: 
        ## get n from stack
        movl nb030nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb030nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total lj energy and update it 
        movapd nb030nf_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb030nf_Vvdw(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## finish if last 
        movl nb030nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel030nf_ia32_sse2.nb030nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb030nf_n(%esp)
        jmp _nb_kernel030nf_ia32_sse2.nb030nf_outer
_nb_kernel030nf_ia32_sse2.nb030nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb030nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel030nf_ia32_sse2.nb030nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel030nf_ia32_sse2.nb030nf_threadloop
_nb_kernel030nf_ia32_sse2.nb030nf_end: 
        emms

        movl nb030nf_nouter(%esp),%eax
        movl nb030nf_ninner(%esp),%ebx
        movl nb030nf_outeriter(%ebp),%ecx
        movl nb030nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb030nf_salign(%esp),%eax
        addl %eax,%esp
        addl $192,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret



