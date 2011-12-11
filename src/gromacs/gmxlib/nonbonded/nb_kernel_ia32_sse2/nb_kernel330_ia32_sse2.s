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


.globl nb_kernel330_ia32_sse2
.globl _nb_kernel330_ia32_sse2
nb_kernel330_ia32_sse2: 
_nb_kernel330_ia32_sse2:        
.set nb330_p_nri, 8
.set nb330_iinr, 12
.set nb330_jindex, 16
.set nb330_jjnr, 20
.set nb330_shift, 24
.set nb330_shiftvec, 28
.set nb330_fshift, 32
.set nb330_gid, 36
.set nb330_pos, 40
.set nb330_faction, 44
.set nb330_charge, 48
.set nb330_p_facel, 52
.set nb330_argkrf, 56
.set nb330_argcrf, 60
.set nb330_Vc, 64
.set nb330_type, 68
.set nb330_p_ntype, 72
.set nb330_vdwparam, 76
.set nb330_Vvdw, 80
.set nb330_p_tabscale, 84
.set nb330_VFtab, 88
.set nb330_invsqrta, 92
.set nb330_dvda, 96
.set nb330_p_gbtabscale, 100
.set nb330_GBtab, 104
.set nb330_p_nthreads, 108
.set nb330_count, 112
.set nb330_mtx, 116
.set nb330_outeriter, 120
.set nb330_inneriter, 124
.set nb330_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse2 use 
.set nb330_ix, 0
.set nb330_iy, 16
.set nb330_iz, 32
.set nb330_iq, 48
.set nb330_dx, 64
.set nb330_dy, 80
.set nb330_dz, 96
.set nb330_two, 112
.set nb330_tsc, 128
.set nb330_qq, 144
.set nb330_c6, 160
.set nb330_c12, 176
.set nb330_fscal, 192
.set nb330_vctot, 208
.set nb330_Vvdwtot, 224
.set nb330_fix, 240
.set nb330_fiy, 256
.set nb330_fiz, 272
.set nb330_half, 288
.set nb330_three, 304
.set nb330_is3, 320
.set nb330_ii3, 324
.set nb330_ntia, 328
.set nb330_innerjjnr, 332
.set nb330_innerk, 336
.set nb330_n, 340
.set nb330_nn1, 344
.set nb330_nri, 348
.set nb330_facel, 352                         ## uses 8 bytes
.set nb330_ntype, 360
.set nb330_nouter, 364
.set nb330_ninner, 368
.set nb330_salign, 372
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $376,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb330_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb330_p_nri(%ebp),%ecx
        movl nb330_p_facel(%ebp),%esi
        movl nb330_p_ntype(%ebp),%edi
        movl (%ecx),%ecx
        movsd (%esi),%xmm7
        movl (%edi),%edi
        movl %ecx,nb330_nri(%esp)
        movsd %xmm7,nb330_facel(%esp)
        movl %edi,nb330_ntype(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb330_nouter(%esp)
        movl %eax,nb330_ninner(%esp)


        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double 0.5 IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb330_half(%esp)
        movl %ebx,nb330_half+4(%esp)
        movsd nb330_half(%esp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## 1.0
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## 2.0
        addpd  %xmm2,%xmm3      ## 3.0
        movapd %xmm1,nb330_half(%esp)
        movapd %xmm2,nb330_two(%esp)
        movapd %xmm3,nb330_three(%esp)
        movl nb330_p_tabscale(%ebp),%eax
        movsd (%eax),%xmm3
        shufpd $0,%xmm3,%xmm3
        movapd %xmm3,nb330_tsc(%esp)

_nb_kernel330_ia32_sse2.nb330_threadloop: 
        movl  nb330_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel330_ia32_sse2.nb330_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel330_ia32_sse2.nb330_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb330_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb330_n(%esp)
        movl %ebx,nb330_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel330_ia32_sse2.nb330_outerstart
        jmp _nb_kernel330_ia32_sse2.nb330_end

_nb_kernel330_ia32_sse2.nb330_outerstart: 
        ## ebx contains number of outer iterations
        addl nb330_nouter(%esp),%ebx
        movl %ebx,nb330_nouter(%esp)

_nb_kernel330_ia32_sse2.nb330_outer: 
        movl  nb330_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb330_is3(%esp)      ## store is3 

        movl  nb330_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movsd (%eax,%ebx,8),%xmm0
        movsd 8(%eax,%ebx,8),%xmm1
        movsd 16(%eax,%ebx,8),%xmm2

        movl  nb330_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx,%esi,4),%ebx    ## ebx =ii 

        movl  nb330_charge(%ebp),%edx
        movsd (%edx,%ebx,8),%xmm3
        mulsd nb330_facel(%esp),%xmm3
        shufpd $0,%xmm3,%xmm3

        movl  nb330_type(%ebp),%edx
        movl  (%edx,%ebx,4),%edx
        imull nb330_ntype(%esp),%edx
        shll  %edx
        movl  %edx,nb330_ntia(%esp)

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb330_pos(%ebp),%eax      ## eax = base of pos[]  

        addsd (%eax,%ebx,8),%xmm0
        addsd 8(%eax,%ebx,8),%xmm1
        addsd 16(%eax,%ebx,8),%xmm2

        movapd %xmm3,nb330_iq(%esp)

        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2

        movapd %xmm0,nb330_ix(%esp)
        movapd %xmm1,nb330_iy(%esp)
        movapd %xmm2,nb330_iz(%esp)

        movl  %ebx,nb330_ii3(%esp)

        ## clear vctot and i forces 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb330_vctot(%esp)
        movapd %xmm4,nb330_Vvdwtot(%esp)
        movapd %xmm4,nb330_fix(%esp)
        movapd %xmm4,nb330_fiy(%esp)
        movapd %xmm4,nb330_fiz(%esp)

        movl  nb330_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb330_pos(%ebp),%esi
        movl  nb330_faction(%ebp),%edi
        movl  nb330_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb330_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb330_ninner(%esp),%ecx
        movl  %ecx,nb330_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb330_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel330_ia32_sse2.nb330_unroll_loop
        jmp   _nb_kernel330_ia32_sse2.nb330_checksingle
_nb_kernel330_ia32_sse2.nb330_unroll_loop: 
        ## twice unrolled innerloop here 
        movl  nb330_innerjjnr(%esp),%edx     ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx
        addl $8,nb330_innerjjnr(%esp)                   ## advance pointer (unrolled 2) 

        movl nb330_charge(%ebp),%esi     ## base of charge[] 

        movlpd (%esi,%eax,8),%xmm3
        movhpd (%esi,%ebx,8),%xmm3

        movapd nb330_iq(%esp),%xmm2
        mulpd  %xmm2,%xmm3
        movapd %xmm3,nb330_qq(%esp)

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1

        movl nb330_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl (%esi,%ebx,4),%ebx
        movl nb330_vdwparam(%ebp),%esi
        shll %eax
        shll %ebx
        movl nb330_ntia(%esp),%edi
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
        movapd %xmm4,nb330_c6(%esp)
        movapd %xmm6,nb330_c12(%esp)

        movl nb330_pos(%ebp),%esi               ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx

        ## move two coordinates to xmm0-xmm2 
        movlpd (%esi,%eax,8),%xmm0
        movlpd 8(%esi,%eax,8),%xmm1
        movlpd 16(%esi,%eax,8),%xmm2
        movhpd (%esi,%ebx,8),%xmm0
        movhpd 8(%esi,%ebx,8),%xmm1
        movhpd 16(%esi,%ebx,8),%xmm2

        movl   nb330_faction(%ebp),%edi

        ## move nb330_ix-iz to xmm4-xmm6 
        movapd nb330_ix(%esp),%xmm4
        movapd nb330_iy(%esp),%xmm5
        movapd nb330_iz(%esp),%xmm6

        ## calc dr 
        subpd %xmm0,%xmm4
        subpd %xmm1,%xmm5
        subpd %xmm2,%xmm6

        ## store dr 
        movapd %xmm4,nb330_dx(%esp)
        movapd %xmm5,nb330_dy(%esp)
        movapd %xmm6,nb330_dz(%esp)
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
        movapd nb330_three(%esp),%xmm1
        mulpd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb330_half(%esp),%xmm0
        subpd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm1
        mulpd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulpd %xmm1,%xmm1       ## lu*lu 
        movapd nb330_three(%esp),%xmm2
        mulpd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb330_half(%esp),%xmm0
        subpd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm2
        mulpd %xmm2,%xmm0       ## xmm0=iter2 of rinv (new lu) 
        mulpd %xmm0,%xmm4       ## xmm4=r 
        mulpd nb330_tsc(%esp),%xmm4

        cvttpd2pi %xmm4,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm5
        subpd %xmm5,%xmm4
        movapd %xmm4,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 

        movd %eax,%mm0
        movd %ebx,%mm1

        movl nb330_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 
        leal  (%eax,%eax,2),%eax        ## idx*=3 (12 total now) 
        leal  (%ebx,%ebx,2),%ebx        ## idx*=3 (12 total now) 

        ## Coulomb 
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
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        mulpd  nb330_two(%esp),%xmm7    ## two*Heps2 
        movapd nb330_qq(%esp),%xmm3
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulpd  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and mm3 fijC 
        ## increment vcoul - then we can get rid of mm5 
        ## update vctot 
        addpd  nb330_vctot(%esp),%xmm5
        movapd %xmm5,nb330_vctot(%esp)

        ## put scalar force on stack temporarily 
        movapd %xmm3,nb330_fscal(%esp)

        ## Dispersion 
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
        ## Dispersion table ready, in xmm4-xmm7                 
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        mulpd  nb330_two(%esp),%xmm7    ## two*Heps2 
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 

        movapd nb330_c6(%esp),%xmm4
        mulpd  %xmm4,%xmm7       ## fijD 
        mulpd  %xmm4,%xmm5       ## Vvdw6 
        addpd  nb330_fscal(%esp),%xmm7   ## add to fscal 

        ## put scalar force back on stack Update Vvdwtot directly 
        addpd  nb330_Vvdwtot(%esp),%xmm5
        movapd %xmm7,nb330_fscal(%esp)
        movapd %xmm5,nb330_Vvdwtot(%esp)

        ## Repulsion 
        movlpd 64(%esi,%eax,8),%xmm4    ## Y1   
        movlpd 64(%esi,%ebx,8),%xmm3    ## Y2
        movhpd 72(%esi,%eax,8),%xmm4    ## Y1 F1        
        movhpd 72(%esi,%ebx,8),%xmm3    ## Y2 F2 

        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 80(%esi,%eax,8),%xmm6    ## G1
        movlpd 80(%esi,%ebx,8),%xmm3    ## G2
        movhpd 88(%esi,%eax,8),%xmm6    ## G1 H1        
        movhpd 88(%esi,%ebx,8),%xmm3    ## G2 H2 

        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## Dispersion table ready, in xmm4-xmm7                 
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        mulpd  nb330_two(%esp),%xmm7    ## two*Heps2 
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 

        movapd nb330_c12(%esp),%xmm4
        mulpd  %xmm4,%xmm7 ## fijR 
        mulpd  %xmm4,%xmm5 ## Vvdw12 
        addpd  nb330_fscal(%esp),%xmm7

        addpd  nb330_Vvdwtot(%esp),%xmm5
        movapd %xmm5,nb330_Vvdwtot(%esp)
        xorpd  %xmm4,%xmm4

        mulpd nb330_tsc(%esp),%xmm7
        mulpd %xmm0,%xmm7
        subpd %xmm7,%xmm4

        movapd nb330_dx(%esp),%xmm0
        movapd nb330_dy(%esp),%xmm1
        movapd nb330_dz(%esp),%xmm2

        movd %mm0,%eax
        movd %mm1,%ebx

        movl   nb330_faction(%ebp),%edi
        mulpd  %xmm4,%xmm0
        mulpd  %xmm4,%xmm1
        mulpd  %xmm4,%xmm2
        ## xmm0-xmm2 contains tx-tz (partial force) 
        ## now update f_i 
        movapd nb330_fix(%esp),%xmm3
        movapd nb330_fiy(%esp),%xmm4
        movapd nb330_fiz(%esp),%xmm5
        addpd  %xmm0,%xmm3
        addpd  %xmm1,%xmm4
        addpd  %xmm2,%xmm5
        movapd %xmm3,nb330_fix(%esp)
        movapd %xmm4,nb330_fiy(%esp)
        movapd %xmm5,nb330_fiz(%esp)
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
        subl $2,nb330_innerk(%esp)
        jl    _nb_kernel330_ia32_sse2.nb330_checksingle
        jmp   _nb_kernel330_ia32_sse2.nb330_unroll_loop
_nb_kernel330_ia32_sse2.nb330_checksingle: 
        movl  nb330_innerk(%esp),%edx
        andl  $1,%edx
        jnz    _nb_kernel330_ia32_sse2.nb330_dosingle
        jmp    _nb_kernel330_ia32_sse2.nb330_updateouterdata
_nb_kernel330_ia32_sse2.nb330_dosingle: 
        movl nb330_charge(%ebp),%esi
        movl nb330_pos(%ebp),%edi
        movl  nb330_innerjjnr(%esp),%ecx
        movl  (%ecx),%eax
        xorpd  %xmm3,%xmm3
        movlpd (%esi,%eax,8),%xmm3      ## xmm6(0) has the charge       
        mulpd  nb330_iq(%esp),%xmm3
        movapd %xmm3,nb330_qq(%esp)

        movd  %eax,%mm0         ## use mmx registers as temp storage 

        movl nb330_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl nb330_vdwparam(%ebp),%esi
        shll %eax
        movl nb330_ntia(%esp),%edi
        addl %edi,%eax

        movlpd (%esi,%eax,8),%xmm6      ## c6a
        movhpd 8(%esi,%eax,8),%xmm6     ## c6a c12a 

        xorpd %xmm7,%xmm7
        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movd  %mm0,%eax
        movd  %mm1,%ebx
        movapd %xmm4,nb330_c6(%esp)
        movapd %xmm6,nb330_c12(%esp)

        movl nb330_pos(%ebp),%esi               ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 

        ## move two coordinates to xmm0-xmm2 
        movlpd (%esi,%eax,8),%xmm0
        movlpd 8(%esi,%eax,8),%xmm1
        movlpd 16(%esi,%eax,8),%xmm2

        movl   nb330_faction(%ebp),%edi

        ## move nb330_ix-iz to xmm4-xmm6 
        movapd nb330_ix(%esp),%xmm4
        movapd nb330_iy(%esp),%xmm5
        movapd nb330_iz(%esp),%xmm6

        ## calc dr 
        subsd %xmm0,%xmm4
        subsd %xmm1,%xmm5
        subsd %xmm2,%xmm6

        ## store dr 
        movapd %xmm4,nb330_dx(%esp)
        movapd %xmm5,nb330_dy(%esp)
        movapd %xmm6,nb330_dz(%esp)
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
        movapd nb330_three(%esp),%xmm1
        mulsd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb330_half(%esp),%xmm0
        subsd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm1
        mulsd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulsd %xmm1,%xmm1       ## lu*lu 
        movapd nb330_three(%esp),%xmm2
        mulsd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb330_half(%esp),%xmm0
        subsd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm2
        mulsd %xmm2,%xmm0       ## xmm0=iter2 of rinv (new lu) 
        mulsd %xmm0,%xmm4       ## xmm4=r 
        mulsd nb330_tsc(%esp),%xmm4

        movd %eax,%mm0
        cvttsd2si %xmm4,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm5
        subsd %xmm5,%xmm4
        movapd %xmm4,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb330_VFtab(%ebp),%esi
        leal  (%eax,%eax,2),%eax        ## idx*=3 (12 total now) 

        ## Coulomb 
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
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb330_two(%esp),%xmm7    ## two*Heps2 
        movapd nb330_qq(%esp),%xmm3
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulsd  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and mm3 fijC 
        ## increment vcoul - then we can get rid of mm5 
        ## update vctot 
        addsd  nb330_vctot(%esp),%xmm5
        movlpd %xmm5,nb330_vctot(%esp)

        ## put scalar force on stack temporarily 
        movapd %xmm3,nb330_fscal(%esp)

        ## Dispersion 
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
        ## Dispersion table ready, in xmm4-xmm7                 
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb330_two(%esp),%xmm7    ## two*Heps2 
        movapd nb330_qq(%esp),%xmm3
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 

        movapd nb330_c6(%esp),%xmm4
        mulsd  %xmm4,%xmm7       ## fijD 
        mulsd  %xmm4,%xmm5       ## Vvdw6 
        addsd  nb330_fscal(%esp),%xmm7   ## add to fscal 

        ## put scalar force back on stack Update Vvdwtot directly 
        addsd  nb330_Vvdwtot(%esp),%xmm5
        movlpd %xmm7,nb330_fscal(%esp)
        movlpd %xmm5,nb330_Vvdwtot(%esp)

        ## Repulsion 
        movlpd 64(%esi,%eax,8),%xmm4    ## Y1   
        movhpd 72(%esi,%eax,8),%xmm4    ## Y1 F1        

        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 
        unpckhpd %xmm3,%xmm5    ## F1 

        movlpd 80(%esi,%eax,8),%xmm6    ## G1   
        movhpd 88(%esi,%eax,8),%xmm6    ## G1 H1        

        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## Dispersion table ready, in xmm4-xmm7                 
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb330_two(%esp),%xmm7    ## two*Heps2 
        movapd nb330_qq(%esp),%xmm3
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 

        movapd nb330_c12(%esp),%xmm4
        mulsd  %xmm4,%xmm7 ## fijR 
        mulsd  %xmm4,%xmm5 ## Vvdw12 
        addsd  nb330_fscal(%esp),%xmm7

        addsd  nb330_Vvdwtot(%esp),%xmm5
        movlpd %xmm5,nb330_Vvdwtot(%esp)
        xorpd  %xmm4,%xmm4

        mulsd nb330_tsc(%esp),%xmm7
        mulsd %xmm0,%xmm7
        subsd %xmm7,%xmm4

        movapd nb330_dx(%esp),%xmm0
        movapd nb330_dy(%esp),%xmm1
        movapd nb330_dz(%esp),%xmm2

        movd %mm0,%eax

        movl   nb330_faction(%ebp),%edi
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2
        ## xmm0-xmm2 contains tx-tz (partial force) 
        ## now update f_i 
        movapd nb330_fix(%esp),%xmm3
        movapd nb330_fiy(%esp),%xmm4
        movapd nb330_fiz(%esp),%xmm5
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm5
        movlpd %xmm3,nb330_fix(%esp)
        movlpd %xmm4,nb330_fiy(%esp)
        movlpd %xmm5,nb330_fiz(%esp)
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
_nb_kernel330_ia32_sse2.nb330_updateouterdata: 
        movl  nb330_ii3(%esp),%ecx
        movl  nb330_faction(%ebp),%edi
        movl  nb330_fshift(%ebp),%esi
        movl  nb330_is3(%esp),%edx

        ## accumulate i forces in xmm0, xmm1, xmm2 
        movapd nb330_fix(%esp),%xmm0
        movapd nb330_fiy(%esp),%xmm1
        movapd nb330_fiz(%esp),%xmm2

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
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm5
        movsd  %xmm3,(%esi,%edx,8)
        movsd  %xmm4,8(%esi,%edx,8)
        movsd  %xmm5,16(%esi,%edx,8)

        ## get n from stack
        movl nb330_n(%esp),%esi
        ## get group index for i particle 
        movl  nb330_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb330_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb330_Vc(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## accumulate total lj energy and update it 
        movapd nb330_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb330_Vvdw(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## finish if last 
        movl nb330_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel330_ia32_sse2.nb330_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb330_n(%esp)
        jmp _nb_kernel330_ia32_sse2.nb330_outer
_nb_kernel330_ia32_sse2.nb330_outerend: 
        ## check if more outer neighborlists remain
        movl  nb330_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel330_ia32_sse2.nb330_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel330_ia32_sse2.nb330_threadloop
_nb_kernel330_ia32_sse2.nb330_end: 
        emms

        movl nb330_nouter(%esp),%eax
        movl nb330_ninner(%esp),%ebx
        movl nb330_outeriter(%ebp),%ecx
        movl nb330_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb330_salign(%esp),%eax
        addl %eax,%esp
        addl $376,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret



.globl nb_kernel330nf_ia32_sse2
.globl _nb_kernel330nf_ia32_sse2
nb_kernel330nf_ia32_sse2:       
_nb_kernel330nf_ia32_sse2:      
.set nb330nf_p_nri, 8
.set nb330nf_iinr, 12
.set nb330nf_jindex, 16
.set nb330nf_jjnr, 20
.set nb330nf_shift, 24
.set nb330nf_shiftvec, 28
.set nb330nf_fshift, 32
.set nb330nf_gid, 36
.set nb330nf_pos, 40
.set nb330nf_faction, 44
.set nb330nf_charge, 48
.set nb330nf_p_facel, 52
.set nb330nf_argkrf, 56
.set nb330nf_argcrf, 60
.set nb330nf_Vc, 64
.set nb330nf_type, 68
.set nb330nf_p_ntype, 72
.set nb330nf_vdwparam, 76
.set nb330nf_Vvdw, 80
.set nb330nf_p_tabscale, 84
.set nb330nf_VFtab, 88
.set nb330nf_invsqrta, 92
.set nb330nf_dvda, 96
.set nb330nf_p_gbtabscale, 100
.set nb330nf_GBtab, 104
.set nb330nf_p_nthreads, 108
.set nb330nf_count, 112
.set nb330nf_mtx, 116
.set nb330nf_outeriter, 120
.set nb330nf_inneriter, 124
.set nb330nf_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb330nf_ix, 0
.set nb330nf_iy, 16
.set nb330nf_iz, 32
.set nb330nf_iq, 48
.set nb330nf_tsc, 64
.set nb330nf_qq, 80
.set nb330nf_c6, 96
.set nb330nf_c12, 112
.set nb330nf_vctot, 128
.set nb330nf_Vvdwtot, 144
.set nb330nf_half, 160
.set nb330nf_three, 176
.set nb330nf_is3, 192
.set nb330nf_ii3, 196
.set nb330nf_ntia, 200
.set nb330nf_innerjjnr, 204
.set nb330nf_innerk, 208
.set nb330nf_n, 212
.set nb330nf_nn1, 216
.set nb330nf_nri, 220
.set nb330nf_facel, 224                       ## uses 8 bytes
.set nb330nf_ntype, 232
.set nb330nf_nouter, 236
.set nb330nf_ninner, 240
.set nb330nf_salign, 244
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $248,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb330nf_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb330nf_p_nri(%ebp),%ecx
        movl nb330nf_p_facel(%ebp),%esi
        movl nb330nf_p_ntype(%ebp),%edi
        movl (%ecx),%ecx
        movsd (%esi),%xmm7
        movl (%edi),%edi
        movl %ecx,nb330nf_nri(%esp)
        movsd %xmm7,nb330nf_facel(%esp)
        movl %edi,nb330nf_ntype(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb330nf_nouter(%esp)
        movl %eax,nb330nf_ninner(%esp)


        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double 0.5 IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb330nf_half(%esp)
        movl %ebx,nb330nf_half+4(%esp)
        movsd nb330nf_half(%esp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## 1.0
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## 2.0
        addpd  %xmm2,%xmm3      ## 3.0
        movapd %xmm1,nb330nf_half(%esp)
        movapd %xmm3,nb330nf_three(%esp)
        movl nb330nf_p_tabscale(%ebp),%eax
        movsd (%eax),%xmm3
        shufpd $0,%xmm3,%xmm3
        movapd %xmm3,nb330nf_tsc(%esp)

_nb_kernel330nf_ia32_sse2.nb330nf_threadloop: 
        movl  nb330nf_count(%ebp),%esi          ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel330nf_ia32_sse2.nb330nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                           ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel330nf_ia32_sse2.nb330nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb330nf_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb330nf_n(%esp)
        movl %ebx,nb330nf_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel330nf_ia32_sse2.nb330nf_outerstart
        jmp _nb_kernel330nf_ia32_sse2.nb330nf_end

_nb_kernel330nf_ia32_sse2.nb330nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb330nf_nouter(%esp),%ebx
        movl %ebx,nb330nf_nouter(%esp)

_nb_kernel330nf_ia32_sse2.nb330nf_outer: 
        movl  nb330nf_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 

        movl  nb330nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movsd (%eax,%ebx,8),%xmm0
        movsd 8(%eax,%ebx,8),%xmm1
        movsd 16(%eax,%ebx,8),%xmm2

        movl  nb330nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx,%esi,4),%ebx    ## ebx =ii 

        movl  nb330nf_charge(%ebp),%edx
        movsd (%edx,%ebx,8),%xmm3
        mulsd nb330nf_facel(%esp),%xmm3
        shufpd $0,%xmm3,%xmm3

        movl  nb330nf_type(%ebp),%edx
        movl  (%edx,%ebx,4),%edx
        imull nb330nf_ntype(%esp),%edx
        shll  %edx
        movl  %edx,nb330nf_ntia(%esp)

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb330nf_pos(%ebp),%eax      ## eax = base of pos[]  

        addsd (%eax,%ebx,8),%xmm0
        addsd 8(%eax,%ebx,8),%xmm1
        addsd 16(%eax,%ebx,8),%xmm2

        movapd %xmm3,nb330nf_iq(%esp)

        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2

        movapd %xmm0,nb330nf_ix(%esp)
        movapd %xmm1,nb330nf_iy(%esp)
        movapd %xmm2,nb330nf_iz(%esp)

        movl  %ebx,nb330nf_ii3(%esp)

        ## clear vctot 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb330nf_vctot(%esp)
        movapd %xmm4,nb330nf_Vvdwtot(%esp)

        movl  nb330nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb330nf_pos(%ebp),%esi
        movl  nb330nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb330nf_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb330nf_ninner(%esp),%ecx
        movl  %ecx,nb330nf_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb330nf_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel330nf_ia32_sse2.nb330nf_unroll_loop
        jmp   _nb_kernel330nf_ia32_sse2.nb330nf_checksingle
_nb_kernel330nf_ia32_sse2.nb330nf_unroll_loop: 
        ## twice unrolled innerloop here 
        movl  nb330nf_innerjjnr(%esp),%edx     ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx
        addl $8,nb330nf_innerjjnr(%esp)                 ## advance pointer (unrolled 2) 

        movl nb330nf_charge(%ebp),%esi     ## base of charge[] 

        movlpd (%esi,%eax,8),%xmm3
        movhpd (%esi,%ebx,8),%xmm3

        movapd nb330nf_iq(%esp),%xmm2
        mulpd  %xmm2,%xmm3
        movapd %xmm3,nb330nf_qq(%esp)

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1

        movl nb330nf_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl (%esi,%ebx,4),%ebx
        movl nb330nf_vdwparam(%ebp),%esi
        shll %eax
        shll %ebx
        movl nb330nf_ntia(%esp),%edi
        addl %edi,%eax
        addl %edi,%ebx

        movlpd (%esi,%eax,8),%xmm6     ## c6a 
        movlpd (%esi,%ebx,8),%xmm7     # # c6b
        movhpd 8(%esi,%eax,8),%xmm6    # # c6a c12a
        movhpd 8(%esi,%ebx,8),%xmm7    # # c6b c12b

        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movd  %mm0,%eax
        movd  %mm1,%ebx
        movapd %xmm4,nb330nf_c6(%esp)
        movapd %xmm6,nb330nf_c12(%esp)

        movl nb330nf_pos(%ebp),%esi             ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx

        ## move two coordinates to xmm0-xmm2 
        movlpd (%esi,%eax,8),%xmm0
        movlpd 8(%esi,%eax,8),%xmm1
        movlpd 16(%esi,%eax,8),%xmm2
        movhpd (%esi,%ebx,8),%xmm0
        movhpd 8(%esi,%ebx,8),%xmm1
        movhpd 16(%esi,%ebx,8),%xmm2

        ## move nb330nf_ix-iz to xmm4-xmm6 
        movapd nb330nf_ix(%esp),%xmm4
        movapd nb330nf_iy(%esp),%xmm5
        movapd nb330nf_iz(%esp),%xmm6

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
        movapd nb330nf_three(%esp),%xmm1
        mulpd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb330nf_half(%esp),%xmm0
        subpd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm1
        mulpd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulpd %xmm1,%xmm1       ## lu*lu 
        movapd nb330nf_three(%esp),%xmm2
        mulpd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb330nf_half(%esp),%xmm0
        subpd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm2
        mulpd %xmm2,%xmm0       ## xmm0=iter2 of rinv (new lu) 
        mulpd %xmm0,%xmm4       ## xmm4=r 
        mulpd nb330nf_tsc(%esp),%xmm4

        cvttpd2pi %xmm4,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm5
        subpd %xmm5,%xmm4
        movapd %xmm4,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 

        movl nb330nf_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 
        leal  (%eax,%eax,2),%eax        ## idx*=3 (12 total now) 
        leal  (%ebx,%ebx,2),%ebx        ## idx*=3 (12 total now) 

        ## Coulomb 
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
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb330nf_qq(%esp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV 
        ## at this point mm5 contains vcoul 
        ## increment vcoul - then we can get rid of mm5 
        ## update vctot 
        addpd  nb330nf_vctot(%esp),%xmm5
        movapd %xmm5,nb330nf_vctot(%esp)

        ## Dispersion 
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
        ## Dispersion table ready, in xmm4-xmm7                 
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 

        mulpd  nb330nf_c6(%esp),%xmm5   ## Vvdw6 

        addpd  nb330nf_Vvdwtot(%esp),%xmm5
        movapd %xmm5,nb330nf_Vvdwtot(%esp)

        ## Repulsion 
        movlpd 64(%esi,%eax,8),%xmm4    ## Y1   
        movlpd 64(%esi,%ebx,8),%xmm3    ## Y2
        movhpd 72(%esi,%eax,8),%xmm4    ## Y1 F1        
        movhpd 72(%esi,%ebx,8),%xmm3    ## Y2 F2 

        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 80(%esi,%eax,8),%xmm6    ## G1
        movlpd 80(%esi,%ebx,8),%xmm3    ## G2
        movhpd 88(%esi,%eax,8),%xmm6    ## G1 H1        
        movhpd 88(%esi,%ebx,8),%xmm3    ## G2 H2 

        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## Dispersion table ready, in xmm4-xmm7                 
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 

        mulpd  nb330nf_c12(%esp),%xmm5   ## Vvdw12 

        addpd  nb330nf_Vvdwtot(%esp),%xmm5
        movapd %xmm5,nb330nf_Vvdwtot(%esp)

        ## should we do one more iteration? 
        subl $2,nb330nf_innerk(%esp)
        jl    _nb_kernel330nf_ia32_sse2.nb330nf_checksingle
        jmp   _nb_kernel330nf_ia32_sse2.nb330nf_unroll_loop
_nb_kernel330nf_ia32_sse2.nb330nf_checksingle: 
        movl  nb330nf_innerk(%esp),%edx
        andl  $1,%edx
        jnz    _nb_kernel330nf_ia32_sse2.nb330nf_dosingle
        jmp    _nb_kernel330nf_ia32_sse2.nb330nf_updateouterdata
_nb_kernel330nf_ia32_sse2.nb330nf_dosingle: 
        movl nb330nf_charge(%ebp),%esi
        movl nb330nf_pos(%ebp),%edi
        movl  nb330nf_innerjjnr(%esp),%ecx
        movl  (%ecx),%eax
        xorpd  %xmm3,%xmm3
        movlpd (%esi,%eax,8),%xmm3      ## xmm6(0) has the charge       
        mulpd  nb330nf_iq(%esp),%xmm3
        movapd %xmm3,nb330nf_qq(%esp)

        movd  %eax,%mm0         ## use mmx registers as temp storage 

        movl nb330nf_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl nb330nf_vdwparam(%ebp),%esi
        shll %eax
        movl nb330nf_ntia(%esp),%edi
        addl %edi,%eax

        movlpd (%esi,%eax,8),%xmm6      ## c6a
        movhpd 8(%esi,%eax,8),%xmm6     ## c6a c12a 

        xorpd %xmm7,%xmm7
        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movd  %mm0,%eax
        movd  %mm1,%ebx
        movapd %xmm4,nb330nf_c6(%esp)
        movapd %xmm6,nb330nf_c12(%esp)

        movl nb330nf_pos(%ebp),%esi             ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 

        ## move two coordinates to xmm0-xmm2 
        movlpd (%esi,%eax,8),%xmm0
        movlpd 8(%esi,%eax,8),%xmm1
        movlpd 16(%esi,%eax,8),%xmm2

        ## move nb330nf_ix-iz to xmm4-xmm6 
        movapd nb330nf_ix(%esp),%xmm4
        movapd nb330nf_iy(%esp),%xmm5
        movapd nb330nf_iz(%esp),%xmm6

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
        movapd nb330nf_three(%esp),%xmm1
        mulsd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb330nf_half(%esp),%xmm0
        subsd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm1
        mulsd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulsd %xmm1,%xmm1       ## lu*lu 
        movapd nb330nf_three(%esp),%xmm2
        mulsd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb330nf_half(%esp),%xmm0
        subsd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm2
        mulsd %xmm2,%xmm0       ## xmm0=iter2 of rinv (new lu) 
        mulsd %xmm0,%xmm4       ## xmm4=r 
        mulsd nb330nf_tsc(%esp),%xmm4

        movd %eax,%mm0
        cvttsd2si %xmm4,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm5
        subsd %xmm5,%xmm4
        movapd %xmm4,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb330nf_VFtab(%ebp),%esi
        leal  (%eax,%eax,2),%eax        ## idx*=3 (12 total now) 

        ## Coulomb 
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
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb330nf_qq(%esp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 
        ## increment vcoul - then we can get rid of mm5 
        ## update vctot 
        addsd  nb330nf_vctot(%esp),%xmm5
        movlpd %xmm5,nb330nf_vctot(%esp)

        ## Dispersion 
        movsd 32(%esi,%eax,8),%xmm4     ## Y1   
        movsd 40(%esi,%eax,8),%xmm5     ## Y1 F1        
        movsd 48(%esi,%eax,8),%xmm6     ## G1   
        movsd 56(%esi,%eax,8),%xmm7     ## G1 H1        
        ## Dispersion table ready, in xmm4-xmm7                 
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb330nf_qq(%esp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 

        mulsd  nb330nf_c6(%esp),%xmm5    ## Vvdw6 

        ## Update Vvdwtot directly 
        addsd  nb330nf_Vvdwtot(%esp),%xmm5
        movlpd %xmm5,nb330nf_Vvdwtot(%esp)

        ## Repulsion 
        movsd 64(%esi,%eax,8),%xmm4     ## Y1   
        movsd 72(%esi,%eax,8),%xmm5     ## Y1 F1        
        movsd 80(%esi,%eax,8),%xmm6     ## G1   
        movsd 88(%esi,%eax,8),%xmm7     ## G1 H1        
        ## Dispersion table ready, in xmm4-xmm7                 
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb330nf_qq(%esp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 

        mulsd  nb330nf_c12(%esp),%xmm5   ## Vvdw12 

        addsd  nb330nf_Vvdwtot(%esp),%xmm5
        movlpd %xmm5,nb330nf_Vvdwtot(%esp)

_nb_kernel330nf_ia32_sse2.nb330nf_updateouterdata: 
        ## get n from stack
        movl nb330nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb330nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb330nf_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb330nf_Vc(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## accumulate total lj energy and update it 
        movapd nb330nf_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb330nf_Vvdw(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## finish if last 
        movl nb330nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel330nf_ia32_sse2.nb330nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb330nf_n(%esp)
        jmp _nb_kernel330nf_ia32_sse2.nb330nf_outer
_nb_kernel330nf_ia32_sse2.nb330nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb330nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel330nf_ia32_sse2.nb330nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel330nf_ia32_sse2.nb330nf_threadloop
_nb_kernel330nf_ia32_sse2.nb330nf_end: 
        emms

        movl nb330nf_nouter(%esp),%eax
        movl nb330nf_ninner(%esp),%ebx
        movl nb330nf_outeriter(%ebp),%ecx
        movl nb330nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb330nf_salign(%esp),%eax
        addl %eax,%esp
        addl $248,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret


