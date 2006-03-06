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




.globl nb_kernel410_ia32_sse2
.globl _nb_kernel410_ia32_sse2
nb_kernel410_ia32_sse2: 
_nb_kernel410_ia32_sse2:        
.set nb410_p_nri, 8
.set nb410_iinr, 12
.set nb410_jindex, 16
.set nb410_jjnr, 20
.set nb410_shift, 24
.set nb410_shiftvec, 28
.set nb410_fshift, 32
.set nb410_gid, 36
.set nb410_pos, 40
.set nb410_faction, 44
.set nb410_charge, 48
.set nb410_p_facel, 52
.set nb410_argkrf, 56
.set nb410_argcrf, 60
.set nb410_Vc, 64
.set nb410_type, 68
.set nb410_p_ntype, 72
.set nb410_vdwparam, 76
.set nb410_Vvdw, 80
.set nb410_p_tabscale, 84
.set nb410_VFtab, 88
.set nb410_invsqrta, 92
.set nb410_dvda, 96
.set nb410_p_gbtabscale, 100
.set nb410_GBtab, 104
.set nb410_p_nthreads, 108
.set nb410_count, 112
.set nb410_mtx, 116
.set nb410_outeriter, 120
.set nb410_inneriter, 124
.set nb410_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse2 use 
.set nb410_ix, 0
.set nb410_iy, 16
.set nb410_iz, 32
.set nb410_iq, 48
.set nb410_dx, 64
.set nb410_dy, 80
.set nb410_dz, 96
.set nb410_two, 112
.set nb410_six, 128
.set nb410_twelve, 144
.set nb410_gbtsc, 160
.set nb410_qq, 176
.set nb410_c6, 192
.set nb410_c12, 208
.set nb410_fscal, 224
.set nb410_vctot, 240
.set nb410_Vvdwtot, 256
.set nb410_fix, 272
.set nb410_fiy, 288
.set nb410_fiz, 304
.set nb410_half, 320
.set nb410_three, 336
.set nb410_r, 352
.set nb410_isai, 368
.set nb410_isaprod, 384
.set nb410_dvdasum, 400
.set nb410_gbscale, 416
.set nb410_ii, 432
.set nb410_is3, 436
.set nb410_ii3, 440
.set nb410_ntia, 444
.set nb410_innerjjnr, 448
.set nb410_innerk, 452
.set nb410_n, 456
.set nb410_nn1, 460
.set nb410_nri, 464
.set nb410_facel, 472                         ## uses 8 bytes
.set nb410_ntype, 480
.set nb410_nouter, 484
.set nb410_ninner, 488
.set nb410_salign, 492
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $496,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb410_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb410_p_nri(%ebp),%ecx
        movl nb410_p_facel(%ebp),%esi
        movl nb410_p_ntype(%ebp),%edi
        movl (%ecx),%ecx
        movsd (%esi),%xmm7
        movl (%edi),%edi
        movl %ecx,nb410_nri(%esp)
        movsd %xmm7,nb410_facel(%esp)
        movl %edi,nb410_ntype(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb410_nouter(%esp)
        movl %eax,nb410_ninner(%esp)


        movl nb410_p_gbtabscale(%ebp),%eax
        movsd (%eax),%xmm5
        shufpd $0,%xmm5,%xmm5
        movapd %xmm5,nb410_gbtsc(%esp)
        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double 0.5 IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb410_half(%esp)
        movl %ebx,nb410_half+4(%esp)
        movsd nb410_half(%esp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## 1.0
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## 2.0
        addpd  %xmm2,%xmm3      ## 3.0
        movapd %xmm3,%xmm4
        addpd  %xmm4,%xmm4      ## 6.0
        movapd %xmm4,%xmm5
        addpd  %xmm5,%xmm5      ## 12.0
        movapd %xmm1,nb410_half(%esp)
        movapd %xmm2,nb410_two(%esp)
        movapd %xmm3,nb410_three(%esp)
        movapd %xmm4,nb410_six(%esp)
        movapd %xmm5,nb410_twelve(%esp)

_nb_kernel410_ia32_sse2.nb410_threadloop: 
        movl  nb410_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel410_ia32_sse2.nb410_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel410_ia32_sse2.nb410_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb410_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb410_n(%esp)
        movl %ebx,nb410_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel410_ia32_sse2.nb410_outerstart
        jmp _nb_kernel410_ia32_sse2.nb410_end

_nb_kernel410_ia32_sse2.nb410_outerstart: 
        ## ebx contains number of outer iterations
        addl nb410_nouter(%esp),%ebx
        movl %ebx,nb410_nouter(%esp)

_nb_kernel410_ia32_sse2.nb410_outer: 
        movl  nb410_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb410_is3(%esp)      ## store is3 

        movl  nb410_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movsd (%eax,%ebx,8),%xmm0
        movsd 8(%eax,%ebx,8),%xmm1
        movsd 16(%eax,%ebx,8),%xmm2

        movl  nb410_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx,%esi,4),%ebx    ## ebx =ii 
        movl  %ebx,nb410_ii(%esp)

        movl  nb410_charge(%ebp),%edx
        movsd (%edx,%ebx,8),%xmm3
        mulsd nb410_facel(%esp),%xmm3
        shufpd $0,%xmm3,%xmm3

        movl  nb410_invsqrta(%ebp),%edx         ## load invsqrta[ii]
        movsd (%edx,%ebx,8),%xmm4
        shufpd $0,%xmm4,%xmm4

        movl  nb410_type(%ebp),%edx
        movl  (%edx,%ebx,4),%edx
        imull nb410_ntype(%esp),%edx
        shll  %edx
        movl  %edx,nb410_ntia(%esp)

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb410_pos(%ebp),%eax      ## eax = base of pos[]  

        addsd (%eax,%ebx,8),%xmm0
        addsd 8(%eax,%ebx,8),%xmm1
        addsd 16(%eax,%ebx,8),%xmm2

        movapd %xmm3,nb410_iq(%esp)
        movapd %xmm4,nb410_isai(%esp)

        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2

        movapd %xmm0,nb410_ix(%esp)
        movapd %xmm1,nb410_iy(%esp)
        movapd %xmm2,nb410_iz(%esp)

        movl  %ebx,nb410_ii3(%esp)

        ## clear vctot and i forces 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb410_vctot(%esp)
        movapd %xmm4,nb410_Vvdwtot(%esp)
        movapd %xmm4,nb410_dvdasum(%esp)
        movapd %xmm4,nb410_fix(%esp)
        movapd %xmm4,nb410_fiy(%esp)
        movapd %xmm4,nb410_fiz(%esp)

        movl  nb410_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb410_pos(%ebp),%esi
        movl  nb410_faction(%ebp),%edi
        movl  nb410_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb410_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb410_ninner(%esp),%ecx
        movl  %ecx,nb410_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb410_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel410_ia32_sse2.nb410_unroll_loop
        jmp   _nb_kernel410_ia32_sse2.nb410_checksingle
_nb_kernel410_ia32_sse2.nb410_unroll_loop: 
        ## twice unrolled innerloop here 
        movl  nb410_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx
        addl $8,nb410_innerjjnr(%esp)             ## advance pointer (unrolled 2) 

        ## load isaj
        movl nb410_invsqrta(%ebp),%esi
        movlpd (%esi,%eax,8),%xmm2
        movhpd (%esi,%ebx,8),%xmm2
        mulpd  nb410_isai(%esp),%xmm2
        movapd %xmm2,nb410_isaprod(%esp)
        movapd %xmm2,%xmm1
        mulpd nb410_gbtsc(%esp),%xmm1
        movapd %xmm1,nb410_gbscale(%esp)

        movl nb410_charge(%ebp),%esi     ## base of charge[] 
        movlpd (%esi,%eax,8),%xmm3
        movhpd (%esi,%ebx,8),%xmm3

        mulpd nb410_iq(%esp),%xmm2
        mulpd  %xmm2,%xmm3
        movapd %xmm3,nb410_qq(%esp)

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1

        movl nb410_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl (%esi,%ebx,4),%ebx
        movl nb410_vdwparam(%ebp),%esi
        shll %eax
        shll %ebx
        movl nb410_ntia(%esp),%edi
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
        movapd %xmm4,nb410_c6(%esp)
        movapd %xmm6,nb410_c12(%esp)

        movl nb410_pos(%ebp),%esi        ## base of pos[] 

        movd  %eax,%mm2
        movd  %ebx,%mm3
        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx

        ## move two coordinates to xmm0-xmm2    
        movlpd (%esi,%eax,8),%xmm0
        movlpd 8(%esi,%eax,8),%xmm1
        movlpd 16(%esi,%eax,8),%xmm2
        movhpd (%esi,%ebx,8),%xmm0
        movhpd 8(%esi,%ebx,8),%xmm1
        movhpd 16(%esi,%ebx,8),%xmm2

        ## move ix-iz to xmm4-xmm6 
        movapd nb410_ix(%esp),%xmm4
        movapd nb410_iy(%esp),%xmm5
        movapd nb410_iz(%esp),%xmm6

        ## calc dr 
        subpd %xmm0,%xmm4
        subpd %xmm1,%xmm5
        subpd %xmm2,%xmm6

        ## store dr 
        movapd %xmm4,nb410_dx(%esp)
        movapd %xmm5,nb410_dy(%esp)
        movapd %xmm6,nb410_dz(%esp)
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
        movapd nb410_three(%esp),%xmm1
        mulpd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb410_half(%esp),%xmm0
        subpd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm1
        mulpd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulpd %xmm1,%xmm1       ## lu*lu 
        movapd nb410_three(%esp),%xmm2
        mulpd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb410_half(%esp),%xmm0
        subpd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm2
        mulpd %xmm2,%xmm0       ## xmm0=rinv 

        mulpd %xmm0,%xmm4       ## xmm4=r 
        movapd %xmm4,nb410_r(%esp)
        mulpd nb410_gbscale(%esp),%xmm4

        cvttpd2pi %xmm4,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm5
        subpd %xmm5,%xmm4
        movapd %xmm4,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 

        movd %eax,%mm0
        movd %ebx,%mm1

        movl nb410_GBtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movapd (%esi,%eax,8),%xmm4      ## Y1 F1        
        movapd (%esi,%ebx,8),%xmm3      ## Y2 F2 
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movapd 16(%esi,%eax,8),%xmm6    ## G1 H1        
        movapd 16(%esi,%ebx,8),%xmm3    ## G2 H2 
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        mulpd  nb410_two(%esp),%xmm7    ## two*Heps2 
        movapd nb410_qq(%esp),%xmm3
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulpd  %xmm7,%xmm3 ## fijC=FF*qq 
        ## get jnr from regs
        movd %mm2,%ecx
        movd %mm3,%edx
        movl nb410_dvda(%ebp),%esi

        ## Calculate dVda
        xorpd %xmm7,%xmm7
        mulpd nb410_gbscale(%esp),%xmm3
        movapd %xmm3,%xmm6
        mulpd  nb410_r(%esp),%xmm6
        addpd  %xmm5,%xmm6
        addpd  nb410_vctot(%esp),%xmm5
        movapd %xmm5,nb410_vctot(%esp)

        ## xmm6=(vcoul+fijC*r)
        subpd  %xmm6,%xmm7
        movapd %xmm7,%xmm6

        ## update dvdasum
        addpd  nb410_dvdasum(%esp),%xmm7
        movapd %xmm7,nb410_dvdasum(%esp)

        ## update j atoms dvdaj
        movhlps %xmm6,%xmm7
        addsd  (%esi,%ecx,8),%xmm6
        addsd  (%esi,%edx,8),%xmm7
        movsd  %xmm6,(%esi,%ecx,8)
        movsd  %xmm7,(%esi,%edx,8)

        ## L-J 
        movapd %xmm0,%xmm4
        mulpd  %xmm0,%xmm4      ## xmm4=rinvsq 

        movapd %xmm4,%xmm6
        mulpd  %xmm4,%xmm6

        mulpd  %xmm4,%xmm6      ## xmm6=rinvsix 
        movapd %xmm6,%xmm4
        mulpd  %xmm4,%xmm4      ## xmm4=rinvtwelve 
        mulpd  nb410_c6(%esp),%xmm6
        mulpd  nb410_c12(%esp),%xmm4
        movapd nb410_Vvdwtot(%esp),%xmm7
        addpd  %xmm4,%xmm7
        mulpd  nb410_twelve(%esp),%xmm4
        subpd  %xmm6,%xmm7
        mulpd  nb410_six(%esp),%xmm6
        movapd %xmm7,nb410_Vvdwtot(%esp)
        subpd  %xmm6,%xmm4
        mulpd  %xmm0,%xmm4
        subpd  %xmm3,%xmm4
        mulpd  %xmm0,%xmm4

        movapd nb410_dx(%esp),%xmm0
        movapd nb410_dy(%esp),%xmm1
        movapd nb410_dz(%esp),%xmm2

        movd %mm0,%eax
        movd %mm1,%ebx

        movl   nb410_faction(%ebp),%edi
        mulpd  %xmm4,%xmm0
        mulpd  %xmm4,%xmm1
        mulpd  %xmm4,%xmm2
        ## xmm0-xmm2 contains tx-tz (partial force) 
        ## now update f_i 
        movapd nb410_fix(%esp),%xmm3
        movapd nb410_fiy(%esp),%xmm4
        movapd nb410_fiz(%esp),%xmm5
        addpd  %xmm0,%xmm3
        addpd  %xmm1,%xmm4
        addpd  %xmm2,%xmm5
        movapd %xmm3,nb410_fix(%esp)
        movapd %xmm4,nb410_fiy(%esp)
        movapd %xmm5,nb410_fiz(%esp)
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
        subl $2,nb410_innerk(%esp)
        jl    _nb_kernel410_ia32_sse2.nb410_checksingle
        jmp   _nb_kernel410_ia32_sse2.nb410_unroll_loop
_nb_kernel410_ia32_sse2.nb410_checksingle: 
        movl  nb410_innerk(%esp),%edx
        andl  $1,%edx
        jnz    _nb_kernel410_ia32_sse2.nb410_dosingle
        jmp    _nb_kernel410_ia32_sse2.nb410_updateouterdata
_nb_kernel410_ia32_sse2.nb410_dosingle: 
        movl nb410_charge(%ebp),%esi
        movl nb410_invsqrta(%ebp),%edx
        movl nb410_pos(%ebp),%edi
        movl  nb410_innerjjnr(%esp),%ecx
        movl  (%ecx),%eax

        xorpd  %xmm6,%xmm6
        movapd %xmm6,%xmm7
        movsd  (%edx,%eax,8),%xmm7
        movlpd (%esi,%eax,8),%xmm6      ## xmm6(0) has the charge
        mulsd  nb410_isai(%esp),%xmm7
        movapd %xmm7,nb410_isaprod(%esp)
        movapd %xmm7,%xmm1
        mulpd nb410_gbtsc(%esp),%xmm1
        movapd %xmm1,nb410_gbscale(%esp)

        mulsd  nb410_iq(%esp),%xmm7
        mulsd  %xmm7,%xmm6
        movapd %xmm6,nb410_qq(%esp)

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movl nb410_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl nb410_vdwparam(%ebp),%esi
        shll %eax
        movl nb410_ntia(%esp),%edi
        addl %edi,%eax

        movlpd (%esi,%eax,8),%xmm6      ## c6a
        movhpd 8(%esi,%eax,8),%xmm6     ## c6a c12a 
        xorpd %xmm7,%xmm7
        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movd  %mm0,%eax
        movapd %xmm4,nb410_c6(%esp)
        movapd %xmm6,nb410_c12(%esp)

        movl nb410_pos(%ebp),%esi        ## base of pos[]

        movd  %eax,%mm2
        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 

        ## move coordinates to xmm0-xmm2        
        movlpd (%esi,%eax,8),%xmm0
        movlpd 8(%esi,%eax,8),%xmm1
        movlpd 16(%esi,%eax,8),%xmm2

        ## move ix-iz to xmm4-xmm6 
        movapd nb410_ix(%esp),%xmm4
        movapd nb410_iy(%esp),%xmm5
        movapd nb410_iz(%esp),%xmm6

        ## calc dr 
        subsd %xmm0,%xmm4
        subsd %xmm1,%xmm5
        subsd %xmm2,%xmm6

        ## store dr 
        movapd %xmm4,nb410_dx(%esp)
        movapd %xmm5,nb410_dy(%esp)
        movapd %xmm6,nb410_dz(%esp)
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
        movapd nb410_three(%esp),%xmm1
        mulsd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb410_half(%esp),%xmm0
        subsd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm1
        mulsd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulsd %xmm1,%xmm1       ## lu*lu 
        movapd nb410_three(%esp),%xmm2
        mulsd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb410_half(%esp),%xmm0
        subsd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm2
        mulsd %xmm2,%xmm0       ## xmm0=rinv 

        mulsd %xmm0,%xmm4       ## xmm4=r 
        movapd %xmm4,nb410_r(%esp)
        mulsd nb410_gbscale(%esp),%xmm4

        movd %eax,%mm0
        cvttsd2si %xmm4,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm5
        subsd %xmm5,%xmm4
        movapd %xmm4,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 

        movl nb410_GBtab(%ebp),%esi

        movapd (%esi,%eax,8),%xmm4      ## Y1 F1        
        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 
        unpckhpd %xmm3,%xmm5    ## F1 

        movapd 16(%esi,%eax,8),%xmm6    ## G1 H1        
        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb410_two(%esp),%xmm7    ## two*Heps2 
        movapd nb410_qq(%esp),%xmm3
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulsd  %xmm7,%xmm3 ## fijC=FF*qq 
        ## get jnr from regs
        movd %mm2,%ebx
        movl nb410_dvda(%ebp),%esi

        ## Calculate dVda
        xorpd %xmm7,%xmm7
        mulsd nb410_gbscale(%esp),%xmm3
        movsd %xmm3,%xmm6
        mulsd  nb410_r(%esp),%xmm6
        addsd  %xmm5,%xmm6
        addsd  nb410_vctot(%esp),%xmm5
        movsd %xmm5,nb410_vctot(%esp)

        ## xmm6=(vcoul+fijC*r)
        subpd %xmm7,%xmm7
        movsd %xmm7,%xmm6

        ## update dvdasum
        addsd  nb410_dvdasum(%esp),%xmm7
        movsd %xmm7,nb410_dvdasum(%esp)

        ## update j atoms dvdaj
        addsd  (%esi,%ebx,8),%xmm6
        movsd  %xmm6,(%esi,%ebx,8)

        ## L-J 
        movapd %xmm0,%xmm4
        mulsd  %xmm0,%xmm4      ## xmm4=rinvsq 


        movapd %xmm4,%xmm6
        mulsd  %xmm4,%xmm6

        mulsd  %xmm4,%xmm6      ## xmm6=rinvsix 
        movapd %xmm6,%xmm4
        mulsd  %xmm4,%xmm4      ## xmm4=rinvtwelve 
        mulsd  nb410_c6(%esp),%xmm6
        mulsd  nb410_c12(%esp),%xmm4
        movapd nb410_Vvdwtot(%esp),%xmm7
        addsd  %xmm4,%xmm7
        mulsd  nb410_twelve(%esp),%xmm4
        subsd  %xmm6,%xmm7
        mulsd  nb410_six(%esp),%xmm6
        movlpd %xmm7,nb410_Vvdwtot(%esp)
        subsd  %xmm6,%xmm4
        mulsd  %xmm0,%xmm4
        subsd  %xmm3,%xmm4
        mulsd  %xmm0,%xmm4

        movapd nb410_dx(%esp),%xmm0
        movapd nb410_dy(%esp),%xmm1
        movapd nb410_dz(%esp),%xmm2

        movd %mm0,%eax

        movl   nb410_faction(%ebp),%edi
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2
        ## xmm0-xmm2 contains tx-tz (partial force) 
        ## now update f_i 
        movapd nb410_fix(%esp),%xmm3
        movapd nb410_fiy(%esp),%xmm4
        movapd nb410_fiz(%esp),%xmm5
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm5
        movlpd %xmm3,nb410_fix(%esp)
        movlpd %xmm4,nb410_fiy(%esp)
        movlpd %xmm5,nb410_fiz(%esp)
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

_nb_kernel410_ia32_sse2.nb410_updateouterdata: 
        movl  nb410_ii3(%esp),%ecx
        movl  nb410_faction(%ebp),%edi
        movl  nb410_fshift(%ebp),%esi
        movl  nb410_is3(%esp),%edx

        ## accumulate i forces in xmm0, xmm1, xmm2 
        movapd nb410_fix(%esp),%xmm0
        movapd nb410_fiy(%esp),%xmm1
        movapd nb410_fiz(%esp),%xmm2

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
        movl nb410_n(%esp),%esi
        ## get group index for i particle 
        movl  nb410_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb410_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb410_Vc(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## accumulate total lj energy and update it 
        movapd nb410_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb410_Vvdw(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## accumulate dVda and update it 
        movapd nb410_dvdasum(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        movl nb410_ii(%esp),%edx
        movl nb410_dvda(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        movsd %xmm7,(%eax,%edx,8)

        ## finish if last 
        movl nb410_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jecxz _nb_kernel410_ia32_sse2.nb410_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb410_n(%esp)
        jmp _nb_kernel410_ia32_sse2.nb410_outer
_nb_kernel410_ia32_sse2.nb410_outerend: 
        ## check if more outer neighborlists remain
        movl  nb410_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jecxz _nb_kernel410_ia32_sse2.nb410_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel410_ia32_sse2.nb410_threadloop
_nb_kernel410_ia32_sse2.nb410_end: 
        emms

        movl nb410_nouter(%esp),%eax
        movl nb410_ninner(%esp),%ebx
        movl nb410_outeriter(%ebp),%ecx
        movl nb410_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb410_salign(%esp),%eax
        addl %eax,%esp
        addl $496,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret







.globl nb_kernel410nf_ia32_sse2
.globl _nb_kernel410nf_ia32_sse2
nb_kernel410nf_ia32_sse2:       
_nb_kernel410nf_ia32_sse2:      
.set nb410nf_p_nri, 8
.set nb410nf_iinr, 12
.set nb410nf_jindex, 16
.set nb410nf_jjnr, 20
.set nb410nf_shift, 24
.set nb410nf_shiftvec, 28
.set nb410nf_fshift, 32
.set nb410nf_gid, 36
.set nb410nf_pos, 40
.set nb410nf_faction, 44
.set nb410nf_charge, 48
.set nb410nf_p_facel, 52
.set nb410nf_argkrf, 56
.set nb410nf_argcrf, 60
.set nb410nf_Vc, 64
.set nb410nf_type, 68
.set nb410nf_p_ntype, 72
.set nb410nf_vdwparam, 76
.set nb410nf_Vvdw, 80
.set nb410nf_p_tabscale, 84
.set nb410nf_VFtab, 88
.set nb410nf_invsqrta, 92
.set nb410nf_dvda, 96
.set nb410nf_p_gbtabscale, 100
.set nb410nf_GBtab, 104
.set nb410nf_p_nthreads, 108
.set nb410nf_count, 112
.set nb410nf_mtx, 116
.set nb410nf_outeriter, 120
.set nb410nf_inneriter, 124
.set nb410nf_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse2 use 
.set nb410nf_ix, 0
.set nb410nf_iy, 16
.set nb410nf_iz, 32
.set nb410nf_iq, 48
.set nb410nf_two, 64
.set nb410nf_gbtsc, 80
.set nb410nf_qq, 96
.set nb410nf_c6, 112
.set nb410nf_c12, 128
.set nb410nf_vctot, 144
.set nb410nf_Vvdwtot, 160
.set nb410nf_half, 176
.set nb410nf_three, 192
.set nb410nf_r, 208
.set nb410nf_isai, 224
.set nb410nf_isaprod, 240
.set nb410nf_gbscale, 256
.set nb410nf_ii, 272
.set nb410nf_is3, 276
.set nb410nf_ii3, 280
.set nb410nf_ntia, 284
.set nb410nf_innerjjnr, 288
.set nb410nf_innerk, 292
.set nb410nf_n, 296
.set nb410nf_nn1, 300
.set nb410nf_nri, 304
.set nb410nf_facel, 312                       ## uses 8 bytes
.set nb410nf_ntype, 320
.set nb410nf_nouter, 324
.set nb410nf_ninner, 328
.set nb410nf_salign, 332
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $336,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb410nf_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb410nf_p_nri(%ebp),%ecx
        movl nb410nf_p_facel(%ebp),%esi
        movl nb410nf_p_ntype(%ebp),%edi
        movl (%ecx),%ecx
        movsd (%esi),%xmm7
        movl (%edi),%edi
        movl %ecx,nb410nf_nri(%esp)
        movsd %xmm7,nb410nf_facel(%esp)
        movl %edi,nb410nf_ntype(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb410nf_nouter(%esp)
        movl %eax,nb410nf_ninner(%esp)


        movl nb410nf_p_gbtabscale(%ebp),%eax
        movsd (%eax),%xmm5
        shufpd $0,%xmm5,%xmm5
        movapd %xmm5,nb410nf_gbtsc(%esp)
        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double 0.5 IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb410nf_half(%esp)
        movl %ebx,nb410nf_half+4(%esp)
        movsd nb410nf_half(%esp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## 1.0
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## 2.0
        addpd  %xmm2,%xmm3      ## 3.0
        movapd %xmm1,nb410nf_half(%esp)
        movapd %xmm2,nb410nf_two(%esp)
        movapd %xmm3,nb410nf_three(%esp)

_nb_kernel410nf_ia32_sse2.nb410nf_threadloop: 
        movl  nb410nf_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel410nf_ia32_sse2.nb410nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel410nf_ia32_sse2.nb410nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb410nf_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb410nf_n(%esp)
        movl %ebx,nb410nf_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
                movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel410nf_ia32_sse2.nb410nf_outerstart
        jmp _nb_kernel410nf_ia32_sse2.nb410nf_end

_nb_kernel410nf_ia32_sse2.nb410nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb410nf_nouter(%esp),%ebx
        movl %ebx,nb410nf_nouter(%esp)

_nb_kernel410nf_ia32_sse2.nb410nf_outer: 
        movl  nb410nf_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb410nf_is3(%esp)            ## store is3 

        movl  nb410nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movsd (%eax,%ebx,8),%xmm0
        movsd 8(%eax,%ebx,8),%xmm1
        movsd 16(%eax,%ebx,8),%xmm2

        movl  nb410nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx,%esi,4),%ebx    ## ebx =ii 
        movl  %ebx,nb410nf_ii(%esp)

        movl  nb410nf_charge(%ebp),%edx
        movsd (%edx,%ebx,8),%xmm3
        mulsd nb410nf_facel(%esp),%xmm3
        shufpd $0,%xmm3,%xmm3

        movl  nb410nf_invsqrta(%ebp),%edx       ## load invsqrta[ii]
        movsd (%edx,%ebx,8),%xmm4
        shufpd $0,%xmm4,%xmm4

        movl  nb410nf_type(%ebp),%edx
        movl  (%edx,%ebx,4),%edx
        imull nb410nf_ntype(%esp),%edx
        shll  %edx
    movl  %edx,nb410nf_ntia(%esp)

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb410nf_pos(%ebp),%eax      ## eax = base of pos[]  

        addsd (%eax,%ebx,8),%xmm0
        addsd 8(%eax,%ebx,8),%xmm1
        addsd 16(%eax,%ebx,8),%xmm2

        movapd %xmm3,nb410nf_iq(%esp)
        movapd %xmm4,nb410nf_isai(%esp)

        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2

        movapd %xmm0,nb410nf_ix(%esp)
        movapd %xmm1,nb410nf_iy(%esp)
        movapd %xmm2,nb410nf_iz(%esp)

        movl  %ebx,nb410nf_ii3(%esp)

        ## clear vctot and Vvdwtot
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb410nf_vctot(%esp)
        movapd %xmm4,nb410nf_Vvdwtot(%esp)

        movl  nb410nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb410nf_pos(%ebp),%esi
        movl  nb410nf_faction(%ebp),%edi
        movl  nb410nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb410nf_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb410nf_ninner(%esp),%ecx
        movl  %ecx,nb410nf_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb410nf_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel410nf_ia32_sse2.nb410nf_unroll_loop
        jmp   _nb_kernel410nf_ia32_sse2.nb410nf_checksingle
_nb_kernel410nf_ia32_sse2.nb410nf_unroll_loop: 
        ## twice unrolled innerloop here 
        movl  nb410nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx
        addl $8,nb410nf_innerjjnr(%esp)             ## advance pointer (unrolled 2) 

        ## load isaj
        movl nb410nf_invsqrta(%ebp),%esi
        movlpd (%esi,%eax,8),%xmm2
        movhpd (%esi,%ebx,8),%xmm2
        mulpd  nb410nf_isai(%esp),%xmm2
        movapd %xmm2,nb410nf_isaprod(%esp)
        movapd %xmm2,%xmm1
        mulpd nb410nf_gbtsc(%esp),%xmm1
        movapd %xmm1,nb410nf_gbscale(%esp)

        movl nb410nf_charge(%ebp),%esi     ## base of charge[] 
        movlpd (%esi,%eax,8),%xmm3
        movhpd (%esi,%ebx,8),%xmm3

        mulpd nb410nf_iq(%esp),%xmm2
        mulpd  %xmm2,%xmm3
        movapd %xmm3,nb410nf_qq(%esp)

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1

        movl nb410nf_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl (%esi,%ebx,4),%ebx
        movl nb410nf_vdwparam(%ebp),%esi
        shll %eax
        shll %ebx
        movl nb410nf_ntia(%esp),%edi
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
        movapd %xmm4,nb410nf_c6(%esp)
        movapd %xmm6,nb410nf_c12(%esp)

        movl nb410nf_pos(%ebp),%esi        ## base of pos[] 

        movd  %eax,%mm2
        movd  %ebx,%mm3
        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx

        ## move two coordinates to xmm0-xmm2    
        movlpd (%esi,%eax,8),%xmm0
        movlpd 8(%esi,%eax,8),%xmm1
        movlpd 16(%esi,%eax,8),%xmm2
        movhpd (%esi,%ebx,8),%xmm0
        movhpd 8(%esi,%ebx,8),%xmm1
        movhpd 16(%esi,%ebx,8),%xmm2

        ## move ix-iz to xmm4-xmm6 
        movapd nb410nf_ix(%esp),%xmm4
        movapd nb410nf_iy(%esp),%xmm5
        movapd nb410nf_iz(%esp),%xmm6

        ## calc dr 
        subpd %xmm0,%xmm4
        subpd %xmm1,%xmm5
        subpd %xmm2,%xmm6

        ## square dr 
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
        movapd nb410nf_three(%esp),%xmm1
        mulpd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb410nf_half(%esp),%xmm0
        subpd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm1
        mulpd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulpd %xmm1,%xmm1       ## lu*lu 
        movapd nb410nf_three(%esp),%xmm2
        mulpd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb410nf_half(%esp),%xmm0
        subpd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm2
        mulpd %xmm2,%xmm0       ## xmm0=rinv 

        mulpd %xmm0,%xmm4       ## xmm4=r 
        movapd %xmm4,nb410nf_r(%esp)
        mulpd nb410nf_gbscale(%esp),%xmm4

        cvttpd2pi %xmm4,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm5
        subpd %xmm5,%xmm4
        movapd %xmm4,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 

        movd %eax,%mm0
        movd %ebx,%mm1

        movl nb410nf_GBtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movapd (%esi,%eax,8),%xmm4      ## Y1 F1        
        movapd (%esi,%ebx,8),%xmm3      ## Y2 F2 
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movapd 16(%esi,%eax,8),%xmm6    ## G1 H1        
        movapd 16(%esi,%ebx,8),%xmm3    ## G2 H2 
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb410nf_qq(%esp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  

        addpd  nb410nf_vctot(%esp),%xmm5
        movapd %xmm5,nb410nf_vctot(%esp)

        ## L-J 
        movapd %xmm0,%xmm4
        mulpd  %xmm0,%xmm4      ## xmm4=rinvsq 

        movapd %xmm4,%xmm6
        mulpd  %xmm4,%xmm6

        mulpd  %xmm4,%xmm6      ## xmm6=rinvsix 
        movapd %xmm6,%xmm4
        mulpd  %xmm4,%xmm4      ## xmm4=rinvtwelve 
        mulpd  nb410nf_c6(%esp),%xmm6
        mulpd  nb410nf_c12(%esp),%xmm4
        movapd nb410nf_Vvdwtot(%esp),%xmm7
        addpd  %xmm4,%xmm7
        subpd  %xmm6,%xmm7
        movapd %xmm7,nb410nf_Vvdwtot(%esp)

        ## should we do one more iteration? 
        subl $2,nb410nf_innerk(%esp)
        jl    _nb_kernel410nf_ia32_sse2.nb410nf_checksingle
        jmp   _nb_kernel410nf_ia32_sse2.nb410nf_unroll_loop
_nb_kernel410nf_ia32_sse2.nb410nf_checksingle: 
        movl  nb410nf_innerk(%esp),%edx
        andl  $1,%edx
        jnz    _nb_kernel410nf_ia32_sse2.nb410nf_dosingle
        jmp    _nb_kernel410nf_ia32_sse2.nb410nf_updateouterdata
_nb_kernel410nf_ia32_sse2.nb410nf_dosingle: 
        movl nb410nf_charge(%ebp),%esi
        movl nb410nf_invsqrta(%ebp),%edx
        movl nb410nf_pos(%ebp),%edi
        movl  nb410nf_innerjjnr(%esp),%ecx
        movl  (%ecx),%eax

        xorpd  %xmm6,%xmm6
        movapd %xmm6,%xmm7
        movsd  (%edx,%eax,8),%xmm7
        movlpd (%esi,%eax,8),%xmm6      ## xmm6(0) has the charge
        mulsd  nb410nf_isai(%esp),%xmm7
        movapd %xmm7,nb410nf_isaprod(%esp)
        movapd %xmm7,%xmm1
        mulpd nb410nf_gbtsc(%esp),%xmm1
        movapd %xmm1,nb410nf_gbscale(%esp)

        mulsd  nb410nf_iq(%esp),%xmm7
        mulsd  %xmm7,%xmm6
        movapd %xmm6,nb410nf_qq(%esp)

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movl nb410nf_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl nb410nf_vdwparam(%ebp),%esi
        shll %eax
        movl nb410nf_ntia(%esp),%edi
        addl %edi,%eax

        movlpd (%esi,%eax,8),%xmm6      ## c6a
        movhpd 8(%esi,%eax,8),%xmm6     ## c6a c12a 

        xorpd %xmm7,%xmm7
        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movd  %mm0,%eax
        movapd %xmm4,nb410nf_c6(%esp)
        movapd %xmm6,nb410nf_c12(%esp)

        movl nb410nf_pos(%ebp),%esi        ## base of pos[]

        movd  %eax,%mm2
        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 

        ## move coordinates to xmm0-xmm2        
        movlpd (%esi,%eax,8),%xmm0
        movlpd 8(%esi,%eax,8),%xmm1
        movlpd 16(%esi,%eax,8),%xmm2

        ## move ix-iz to xmm4-xmm6 
        movapd nb410nf_ix(%esp),%xmm4
        movapd nb410nf_iy(%esp),%xmm5
        movapd nb410nf_iz(%esp),%xmm6

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
        movapd nb410nf_three(%esp),%xmm1
        mulsd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb410nf_half(%esp),%xmm0
        subsd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm1
        mulsd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulsd %xmm1,%xmm1       ## lu*lu 
        movapd nb410nf_three(%esp),%xmm2
        mulsd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb410nf_half(%esp),%xmm0
        subsd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm2
        mulsd %xmm2,%xmm0       ## xmm0=rinv 

        mulsd %xmm0,%xmm4       ## xmm4=r 
        movapd %xmm4,nb410nf_r(%esp)
        mulsd nb410nf_gbscale(%esp),%xmm4

        movd %eax,%mm0
        cvttsd2si %xmm4,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm5
        subsd %xmm5,%xmm4
        movapd %xmm4,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 

        movl nb410nf_GBtab(%ebp),%esi

        movapd (%esi,%eax,8),%xmm4      ## Y1 F1        
        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 
        unpckhpd %xmm3,%xmm5    ## F1 

        movapd 16(%esi,%eax,8),%xmm6    ## G1 H1        
        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb410nf_qq(%esp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  

        addsd  nb410nf_vctot(%esp),%xmm5
        movsd %xmm5,nb410nf_vctot(%esp)

        ## L-J 
        movapd %xmm0,%xmm4
        mulsd  %xmm0,%xmm4      ## xmm4=rinvsq 


        movapd %xmm4,%xmm6
        mulsd  %xmm4,%xmm6

        mulsd  %xmm4,%xmm6      ## xmm6=rinvsix 
        movapd %xmm6,%xmm4
        mulsd  %xmm4,%xmm4      ## xmm4=rinvtwelve 
        mulsd  nb410nf_c6(%esp),%xmm6
        mulsd  nb410nf_c12(%esp),%xmm4
        movapd nb410nf_Vvdwtot(%esp),%xmm7
        addsd  %xmm4,%xmm7
        subsd  %xmm6,%xmm7
        movlpd %xmm7,nb410nf_Vvdwtot(%esp)

_nb_kernel410nf_ia32_sse2.nb410nf_updateouterdata: 
        movl  nb410nf_ii3(%esp),%ecx
        movl  nb410nf_is3(%esp),%edx

        ## get n from stack
        movl nb410nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb410nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb410nf_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb410nf_Vc(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## accumulate total lj energy and update it 
        movapd nb410nf_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb410nf_Vvdw(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## finish if last 
        movl nb410nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jecxz _nb_kernel410nf_ia32_sse2.nb410nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb410nf_n(%esp)
        jmp _nb_kernel410nf_ia32_sse2.nb410nf_outer
_nb_kernel410nf_ia32_sse2.nb410nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb410nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jecxz _nb_kernel410nf_ia32_sse2.nb410nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel410nf_ia32_sse2.nb410nf_threadloop
_nb_kernel410nf_ia32_sse2.nb410nf_end: 
        emms

        movl nb410nf_nouter(%esp),%eax
        movl nb410nf_ninner(%esp),%ebx
        movl nb410nf_outeriter(%ebp),%ecx
        movl nb410nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb410nf_salign(%esp),%eax
        addl %eax,%esp
        addl $336,%esp
        popl %edi
        popl %esi
    popl %edx
    popl %ecx
    popl %ebx
    popl %eax
        leave
        ret



