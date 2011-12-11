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




.globl nb_kernel331_ia32_sse2
.globl _nb_kernel331_ia32_sse2
nb_kernel331_ia32_sse2: 
_nb_kernel331_ia32_sse2:        
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
.set nb331_argkrf, 56
.set nb331_argcrf, 60
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
        ## bottom of stack is cache-aligned for sse2 use 
.set nb331_ixO, 0
.set nb331_iyO, 16
.set nb331_izO, 32
.set nb331_ixH1, 48
.set nb331_iyH1, 64
.set nb331_izH1, 80
.set nb331_ixH2, 96
.set nb331_iyH2, 112
.set nb331_izH2, 128
.set nb331_iqO, 144
.set nb331_iqH, 160
.set nb331_dxO, 176
.set nb331_dyO, 192
.set nb331_dzO, 208
.set nb331_dxH1, 224
.set nb331_dyH1, 240
.set nb331_dzH1, 256
.set nb331_dxH2, 272
.set nb331_dyH2, 288
.set nb331_dzH2, 304
.set nb331_qqO, 320
.set nb331_qqH, 336
.set nb331_rinvO, 352
.set nb331_rinvH1, 368
.set nb331_rinvH2, 384
.set nb331_rO, 400
.set nb331_rH1, 416
.set nb331_rH2, 432
.set nb331_tsc, 448
.set nb331_two, 464
.set nb331_c6, 480
.set nb331_c12, 496
.set nb331_vctot, 512
.set nb331_Vvdwtot, 528
.set nb331_fixO, 544
.set nb331_fiyO, 560
.set nb331_fizO, 576
.set nb331_fixH1, 592
.set nb331_fiyH1, 608
.set nb331_fizH1, 624
.set nb331_fixH2, 640
.set nb331_fiyH2, 656
.set nb331_fizH2, 672
.set nb331_fjx, 688
.set nb331_fjy, 704
.set nb331_fjz, 720
.set nb331_half, 736
.set nb331_three, 752
.set nb331_is3, 768
.set nb331_ii3, 772
.set nb331_ntia, 776
.set nb331_innerjjnr, 780
.set nb331_innerk, 784
.set nb331_n, 788
.set nb331_nn1, 792
.set nb331_nri, 796
.set nb331_nouter, 800
.set nb331_ninner, 804
.set nb331_salign, 808
        pushl %ebp
        movl %esp,%ebp
    pushl %eax
    pushl %ebx
    pushl %ecx
    pushl %edx
        pushl %esi
        pushl %edi
        subl $812,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb331_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb331_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb331_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb331_nouter(%esp)
        movl %eax,nb331_ninner(%esp)


        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double 0.5 IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb331_half(%esp)
        movl %ebx,nb331_half+4(%esp)
        movsd nb331_half(%esp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## 1.0
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## 2.0
        addpd  %xmm2,%xmm3      ## 3.0
        movapd %xmm1,nb331_half(%esp)
        movapd %xmm2,nb331_two(%esp)
        movapd %xmm3,nb331_three(%esp)
        movl nb331_p_tabscale(%ebp),%eax
        movsd (%eax),%xmm5
        shufpd $0,%xmm5,%xmm5
        movapd %xmm5,nb331_tsc(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb331_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx),%ebx           ## ebx =ii 

        movl  nb331_charge(%ebp),%edx
        movsd (%edx,%ebx,8),%xmm3
        movsd 8(%edx,%ebx,8),%xmm4
        movl nb331_p_facel(%ebp),%esi
        movsd (%esi),%xmm5
        mulsd  %xmm5,%xmm3
        mulsd  %xmm5,%xmm4

        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        movapd %xmm3,nb331_iqO(%esp)
        movapd %xmm4,nb331_iqH(%esp)

        movl  nb331_type(%ebp),%edx
        movl  (%edx,%ebx,4),%ecx
        shll  %ecx
        movl nb331_p_ntype(%ebp),%edi
        imull (%edi),%ecx     ## ecx = ntia = 2*ntype*type[ii0] 
        movl  %ecx,nb331_ntia(%esp)
_nb_kernel331_ia32_sse2.nb331_threadloop: 
        movl  nb331_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel331_ia32_sse2.nb331_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel331_ia32_sse2.nb331_spinlock

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
        jg  _nb_kernel331_ia32_sse2.nb331_outerstart
        jmp _nb_kernel331_ia32_sse2.nb331_end

_nb_kernel331_ia32_sse2.nb331_outerstart: 
        ## ebx contains number of outer iterations
        addl nb331_nouter(%esp),%ebx
        movl %ebx,nb331_nouter(%esp)

_nb_kernel331_ia32_sse2.nb331_outer: 
        movl  nb331_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb331_is3(%esp)      ## store is3 

        movl  nb331_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movsd (%eax,%ebx,8),%xmm0
        movsd 8(%eax,%ebx,8),%xmm1
        movsd 16(%eax,%ebx,8),%xmm2

        movl  nb331_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx,%esi,4),%ebx    ## ebx =ii 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb331_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb331_ii3(%esp)

        addsd (%eax,%ebx,8),%xmm3
        addsd 8(%eax,%ebx,8),%xmm4
        addsd 16(%eax,%ebx,8),%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb331_ixO(%esp)
        movapd %xmm4,nb331_iyO(%esp)
        movapd %xmm5,nb331_izO(%esp)

        movsd %xmm0,%xmm3
        movsd %xmm1,%xmm4
        movsd %xmm2,%xmm5
        addsd 24(%eax,%ebx,8),%xmm0
        addsd 32(%eax,%ebx,8),%xmm1
        addsd 40(%eax,%ebx,8),%xmm2
        addsd 48(%eax,%ebx,8),%xmm3
        addsd 56(%eax,%ebx,8),%xmm4
        addsd 64(%eax,%ebx,8),%xmm5

        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm0,nb331_ixH1(%esp)
        movapd %xmm1,nb331_iyH1(%esp)
        movapd %xmm2,nb331_izH1(%esp)
        movapd %xmm3,nb331_ixH2(%esp)
        movapd %xmm4,nb331_iyH2(%esp)
        movapd %xmm5,nb331_izH2(%esp)

        ## clear vctot and i forces 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb331_vctot(%esp)
        movapd %xmm4,nb331_Vvdwtot(%esp)
        movapd %xmm4,nb331_fixO(%esp)
        movapd %xmm4,nb331_fiyO(%esp)
        movapd %xmm4,nb331_fizO(%esp)
        movapd %xmm4,nb331_fixH1(%esp)
        movapd %xmm4,nb331_fiyH1(%esp)
        movapd %xmm4,nb331_fizH1(%esp)
        movapd %xmm4,nb331_fixH2(%esp)
        movapd %xmm4,nb331_fiyH2(%esp)
        movapd %xmm4,nb331_fizH2(%esp)

        movl  nb331_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx     ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb331_pos(%ebp),%esi
        movl  nb331_faction(%ebp),%edi
        movl  nb331_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb331_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb331_ninner(%esp),%ecx
        movl  %ecx,nb331_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb331_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel331_ia32_sse2.nb331_unroll_loop
        jmp   _nb_kernel331_ia32_sse2.nb331_checksingle
_nb_kernel331_ia32_sse2.nb331_unroll_loop: 
        ## twice unrolled innerloop here 
        movl  nb331_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx

        addl $8,nb331_innerjjnr(%esp)             ## advance pointer (unrolled 2) 

        movl nb331_charge(%ebp),%esi     ## base of charge[] 

        movlpd (%esi,%eax,8),%xmm3
        movhpd (%esi,%ebx,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb331_iqO(%esp),%xmm3
        mulpd  nb331_iqH(%esp),%xmm4
        movapd  %xmm3,nb331_qqO(%esp)
        movapd  %xmm4,nb331_qqH(%esp)

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1
        movl nb331_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl (%esi,%ebx,4),%ebx
        movl nb331_vdwparam(%ebp),%esi
        shll %eax
        shll %ebx
        movl nb331_ntia(%esp),%edi
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
        movapd %xmm4,nb331_c6(%esp)
        movapd %xmm6,nb331_c12(%esp)

        movl nb331_pos(%ebp),%esi        ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx

        ## move two coordinates to xmm0-xmm2    
        movlpd (%esi,%eax,8),%xmm0
        movlpd 8(%esi,%eax,8),%xmm1
        movlpd 16(%esi,%eax,8),%xmm2
        movhpd (%esi,%ebx,8),%xmm0
        movhpd 8(%esi,%ebx,8),%xmm1
        movhpd 16(%esi,%ebx,8),%xmm2

        ## move ixO-izO to xmm4-xmm6 
        movapd nb331_ixO(%esp),%xmm4
        movapd nb331_iyO(%esp),%xmm5
        movapd nb331_izO(%esp),%xmm6

        ## calc dr 
        subpd %xmm0,%xmm4
        subpd %xmm1,%xmm5
        subpd %xmm2,%xmm6

        ## store dr 
        movapd %xmm4,nb331_dxO(%esp)
        movapd %xmm5,nb331_dyO(%esp)
        movapd %xmm6,nb331_dzO(%esp)
        ## square it 
        mulpd %xmm4,%xmm4
        mulpd %xmm5,%xmm5
        mulpd %xmm6,%xmm6
        addpd %xmm5,%xmm4
        addpd %xmm6,%xmm4
        movapd %xmm4,%xmm7
        ## rsqO in xmm7 

        ## move ixH1-izH1 to xmm4-xmm6 
        movapd nb331_ixH1(%esp),%xmm4
        movapd nb331_iyH1(%esp),%xmm5
        movapd nb331_izH1(%esp),%xmm6

        ## calc dr 
        subpd %xmm0,%xmm4
        subpd %xmm1,%xmm5
        subpd %xmm2,%xmm6

        ## store dr 
        movapd %xmm4,nb331_dxH1(%esp)
        movapd %xmm5,nb331_dyH1(%esp)
        movapd %xmm6,nb331_dzH1(%esp)
        ## square it 
        mulpd %xmm4,%xmm4
        mulpd %xmm5,%xmm5
        mulpd %xmm6,%xmm6
        addpd %xmm5,%xmm6
        addpd %xmm4,%xmm6
        ## rsqH1 in xmm6 

        ## move ixH2-izH2 to xmm3-xmm5  
        movapd nb331_ixH2(%esp),%xmm3
        movapd nb331_iyH2(%esp),%xmm4
        movapd nb331_izH2(%esp),%xmm5

        ## calc dr 
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5

        ## store dr 
        movapd %xmm3,nb331_dxH2(%esp)
        movapd %xmm4,nb331_dyH2(%esp)
        movapd %xmm5,nb331_dzH2(%esp)
        ## square it 
        mulpd %xmm3,%xmm3
        mulpd %xmm4,%xmm4
        mulpd %xmm5,%xmm5
        addpd %xmm4,%xmm5
        addpd %xmm3,%xmm5
        ## rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

        ## start with rsqO - put seed in xmm2 
        cvtpd2ps %xmm7,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb331_three(%esp),%xmm4
        mulpd   %xmm7,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulpd   nb331_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm7,%xmm2
        movapd %xmm4,%xmm3
        mulpd %xmm4,%xmm4       ## lu*lu 
        mulpd %xmm4,%xmm2       ## rsq*lu*lu 
        movapd nb331_three(%esp),%xmm4
        subpd %xmm2,%xmm4       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulpd nb331_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,nb331_rinvO(%esp)         ## rinvO in xmm4 
        mulpd   %xmm4,%xmm7
        movapd  %xmm7,nb331_rO(%esp)    ## r in xmm7 

        ## rsqH1 - seed in xmm2 
        cvtpd2ps %xmm6,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb331_three(%esp),%xmm4
        mulpd   %xmm6,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulpd   nb331_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm6,%xmm2
        movapd %xmm4,%xmm3
        mulpd %xmm4,%xmm4       ## lu*lu 
        mulpd %xmm4,%xmm2       ## rsq*lu*lu 
        movapd nb331_three(%esp),%xmm4
        subpd %xmm2,%xmm4       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulpd nb331_half(%esp),%xmm4   ## rinv 
        movapd %xmm4,nb331_rinvH1(%esp)         ## rinvH1 
        mulpd  %xmm4,%xmm6
        movapd %xmm6,nb331_rH1(%esp)    ## rH1 

        ## rsqH2 - seed in xmm2 
        cvtpd2ps %xmm5,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb331_three(%esp),%xmm4
        mulpd   %xmm5,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulpd   nb331_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm5,%xmm2
        movapd %xmm4,%xmm3
        mulpd %xmm4,%xmm4       ## lu*lu 
        mulpd %xmm4,%xmm2       ## rsq*lu*lu 
        movapd nb331_three(%esp),%xmm4
        subpd %xmm2,%xmm4       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulpd nb331_half(%esp),%xmm4   ## rinv 
        movapd %xmm4,nb331_rinvH2(%esp)   ## rinv 
        mulpd %xmm4,%xmm5
        movapd %xmm5,nb331_rH2(%esp)   ## r 

        ## do O interactions 
        ## rO is still in xmm7 
        movapd nb331_rinvO(%esp),%xmm0
        mulpd   nb331_tsc(%esp),%xmm7
        cvttpd2pi %xmm7,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movd %eax,%mm0
        movd %ebx,%mm1
        movl nb331_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 
        leal  (%eax,%eax,2),%eax ## idx *= 3 (total *=12 now) 
        leal  (%ebx,%ebx,2),%ebx

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
    mulpd  nb331_two(%esp),%xmm7         ## two*Heps2 
    movapd nb331_qqO(%esp),%xmm0
    addpd  %xmm6,%xmm7
    addpd  %xmm5,%xmm7 ## xmm7=FF 
    mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addpd  %xmm4,%xmm5 ## xmm5=VV 
    mulpd  %xmm0,%xmm5 ## vcoul=qq*VV  
    mulpd  %xmm7,%xmm0 ## fijC=FF*qq 
        ## at this point xmm5 contains vcoul and xmm0 fijC 
        ## increment vcoul - then we can get rid of mm5 
        ## update vctot 
        addpd  nb331_vctot(%esp),%xmm5
        movapd %xmm5,nb331_vctot(%esp)

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
        mulpd  nb331_two(%esp),%xmm7    ## two*Heps2 
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 

        movapd nb331_c6(%esp),%xmm4
        mulpd  %xmm4,%xmm7       ## fijD 
        mulpd  %xmm4,%xmm5       ## Vvdw6 
        addpd  %xmm7,%xmm0      ## fscal summation 

        addpd  nb331_Vvdwtot(%esp),%xmm5
        movapd %xmm5,nb331_Vvdwtot(%esp)

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
        mulpd  nb331_two(%esp),%xmm7    ## two*Heps2 
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 

        movapd nb331_c12(%esp),%xmm4
        mulpd  %xmm4,%xmm7 ## fijR 
        mulpd  %xmm4,%xmm5 ## Vvdw12 
        addpd  %xmm0,%xmm7      ## fscal summation 
        mulpd nb331_tsc(%esp),%xmm7
        mulpd nb331_rinvO(%esp),%xmm7

        addpd  nb331_Vvdwtot(%esp),%xmm5
        movapd %xmm5,nb331_Vvdwtot(%esp)

        xorpd %xmm4,%xmm4
        subpd %xmm7,%xmm4
        movapd nb331_dxO(%esp),%xmm0
        movapd nb331_dyO(%esp),%xmm1
        movapd nb331_dzO(%esp),%xmm2

        movd %mm0,%eax
        movd %mm1,%ebx

        movapd nb331_dxO(%esp),%xmm0
        movapd nb331_dyO(%esp),%xmm1
        movapd nb331_dzO(%esp),%xmm2
        mulpd  %xmm4,%xmm0
        mulpd  %xmm4,%xmm1
        mulpd  %xmm4,%xmm2      ## tx in xmm0-xmm2 

        ## update O forces 
        movapd nb331_fixO(%esp),%xmm3
        movapd nb331_fiyO(%esp),%xmm4
        movapd nb331_fizO(%esp),%xmm7
        addpd  %xmm0,%xmm3
        addpd  %xmm1,%xmm4
        addpd  %xmm2,%xmm7
        movapd %xmm3,nb331_fixO(%esp)
        movapd %xmm4,nb331_fiyO(%esp)
        movapd %xmm7,nb331_fizO(%esp)
        ## update j forces with water O 
        movapd %xmm0,nb331_fjx(%esp)
        movapd %xmm1,nb331_fjy(%esp)
        movapd %xmm2,nb331_fjz(%esp)

        ## Done with O interactions - now H1! 
        movapd nb331_rH1(%esp),%xmm7
        mulpd nb331_tsc(%esp),%xmm7
        cvttpd2pi %xmm7,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb331_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 
        leal  (%eax,%eax,2),%eax ## idx *= 3 (total *=12 now)   
        leal  (%ebx,%ebx,2),%ebx

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
        mulpd  nb331_two(%esp),%xmm7    ## two*Heps2 
        movapd nb331_qqH(%esp),%xmm3
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulpd  %xmm7,%xmm3 ## fijC=FF*qq 
    ## at this point mm5 contains vcoul and xmm3 fijC 
    ## increment vcoul 
        xorpd  %xmm4,%xmm4
    addpd  nb331_vctot(%esp),%xmm5
        mulpd  nb331_rinvH1(%esp),%xmm3
    movapd %xmm5,nb331_vctot(%esp)
        mulpd  nb331_tsc(%esp),%xmm3
        subpd %xmm3,%xmm4

        movapd nb331_dxH1(%esp),%xmm0
        movapd nb331_dyH1(%esp),%xmm1
        movapd nb331_dzH1(%esp),%xmm2
        mulpd  %xmm4,%xmm0
        mulpd  %xmm4,%xmm1
        mulpd  %xmm4,%xmm2

        ## update H1 forces 
        movapd nb331_fixH1(%esp),%xmm3
        movapd nb331_fiyH1(%esp),%xmm4
        movapd nb331_fizH1(%esp),%xmm7
        addpd  %xmm0,%xmm3
        addpd  %xmm1,%xmm4
        addpd  %xmm2,%xmm7
        movapd %xmm3,nb331_fixH1(%esp)
        movapd %xmm4,nb331_fiyH1(%esp)
        movapd %xmm7,nb331_fizH1(%esp)
        ## update j forces with water H1 
        addpd  nb331_fjx(%esp),%xmm0
        addpd  nb331_fjy(%esp),%xmm1
        addpd  nb331_fjz(%esp),%xmm2
        movapd %xmm0,nb331_fjx(%esp)
        movapd %xmm1,nb331_fjy(%esp)
        movapd %xmm2,nb331_fjz(%esp)

        ## Done with H1, finally we do H2 interactions 
        movapd nb331_rH2(%esp),%xmm7
        mulpd   nb331_tsc(%esp),%xmm7
        cvttpd2pi %xmm7,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb331_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 
        leal  (%eax,%eax,2),%eax ## idx *= 3 (total *=12 now)
        leal  (%ebx,%ebx,2),%ebx

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
        mulpd  nb331_two(%esp),%xmm7    ## two*Heps2 
        movapd nb331_qqH(%esp),%xmm3
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulpd  %xmm7,%xmm3 ## fijC=FF*qq 
    ## at this point mm5 contains vcoul and xmm3 fijC 
    ## increment vcoul 
        xorpd  %xmm4,%xmm4
    addpd  nb331_vctot(%esp),%xmm5
        mulpd  nb331_rinvH2(%esp),%xmm3
    movapd %xmm5,nb331_vctot(%esp)
        mulpd  nb331_tsc(%esp),%xmm3
        subpd  %xmm3,%xmm4

        movapd nb331_dxH2(%esp),%xmm0
        movapd nb331_dyH2(%esp),%xmm1
        movapd nb331_dzH2(%esp),%xmm2
        mulpd  %xmm4,%xmm0
        mulpd  %xmm4,%xmm1
        mulpd  %xmm4,%xmm2

    movd %mm0,%eax
    movd %mm1,%ebx

        ## update H2 forces 
        movapd nb331_fixH2(%esp),%xmm3
        movapd nb331_fiyH2(%esp),%xmm4
        movapd nb331_fizH2(%esp),%xmm7
        addpd  %xmm0,%xmm3
        addpd  %xmm1,%xmm4
        addpd  %xmm2,%xmm7
        movapd %xmm3,nb331_fixH2(%esp)
        movapd %xmm4,nb331_fiyH2(%esp)
        movapd %xmm7,nb331_fizH2(%esp)

        movl nb331_faction(%ebp),%edi
        ## update j forces 
        addpd  nb331_fjx(%esp),%xmm0
        addpd  nb331_fjy(%esp),%xmm1
        addpd  nb331_fjz(%esp),%xmm2

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
        subl $2,nb331_innerk(%esp)
        jl    _nb_kernel331_ia32_sse2.nb331_checksingle
        jmp   _nb_kernel331_ia32_sse2.nb331_unroll_loop
_nb_kernel331_ia32_sse2.nb331_checksingle: 
        movl  nb331_innerk(%esp),%edx
        andl  $1,%edx
        jnz   _nb_kernel331_ia32_sse2.nb331_dosingle
        jmp   _nb_kernel331_ia32_sse2.nb331_updateouterdata
_nb_kernel331_ia32_sse2.nb331_dosingle: 
        movl  nb331_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax

        movl nb331_charge(%ebp),%esi     ## base of charge[] 
        xorpd %xmm3,%xmm3
        movlpd (%esi,%eax,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb331_iqO(%esp),%xmm3
        mulpd  nb331_iqH(%esp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movapd  %xmm3,nb331_qqO(%esp)
        movapd  %xmm4,nb331_qqH(%esp)

        movl nb331_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl nb331_vdwparam(%ebp),%esi
        shll %eax
        movl nb331_ntia(%esp),%edi
        addl %edi,%eax

        movlpd (%esi,%eax,8),%xmm6      ## c6a
        movhpd 8(%esi,%eax,8),%xmm6     ## c6a c12a 

        xorpd %xmm7,%xmm7
        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movd  %mm0,%eax
        movapd %xmm4,nb331_c6(%esp)
        movapd %xmm6,nb331_c12(%esp)

        movl nb331_pos(%ebp),%esi        ## base of pos[] 
        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 

        ## move coords to xmm0-xmm2 
        movlpd (%esi,%eax,8),%xmm0
        movlpd 8(%esi,%eax,8),%xmm1
        movlpd 16(%esi,%eax,8),%xmm2

        ## move ixO-izO to xmm4-xmm6 
        movapd nb331_ixO(%esp),%xmm4
        movapd nb331_iyO(%esp),%xmm5
        movapd nb331_izO(%esp),%xmm6

        ## calc dr 
        subsd %xmm0,%xmm4
        subsd %xmm1,%xmm5
        subsd %xmm2,%xmm6

        ## store dr 
        movapd %xmm4,nb331_dxO(%esp)
        movapd %xmm5,nb331_dyO(%esp)
        movapd %xmm6,nb331_dzO(%esp)
        ## square it 
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        mulsd %xmm6,%xmm6
        addsd %xmm5,%xmm4
        addsd %xmm6,%xmm4
        movapd %xmm4,%xmm7
        ## rsqO in xmm7 

        ## move ixH1-izH1 to xmm4-xmm6 
        movapd nb331_ixH1(%esp),%xmm4
        movapd nb331_iyH1(%esp),%xmm5
        movapd nb331_izH1(%esp),%xmm6

        ## calc dr 
        subsd %xmm0,%xmm4
        subsd %xmm1,%xmm5
        subsd %xmm2,%xmm6

        ## store dr 
        movapd %xmm4,nb331_dxH1(%esp)
        movapd %xmm5,nb331_dyH1(%esp)
        movapd %xmm6,nb331_dzH1(%esp)
        ## square it 
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        mulsd %xmm6,%xmm6
        addsd %xmm5,%xmm6
        addsd %xmm4,%xmm6
        ## rsqH1 in xmm6 

        ## move ixH2-izH2 to xmm3-xmm5  
        movapd nb331_ixH2(%esp),%xmm3
        movapd nb331_iyH2(%esp),%xmm4
        movapd nb331_izH2(%esp),%xmm5

        ## calc dr 
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5

        ## store dr 
        movapd %xmm3,nb331_dxH2(%esp)
        movapd %xmm4,nb331_dyH2(%esp)
        movapd %xmm5,nb331_dzH2(%esp)
        ## square it 
        mulsd %xmm3,%xmm3
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        addsd %xmm4,%xmm5
        addsd %xmm3,%xmm5
        ## rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

        ## start with rsqO - put seed in xmm2 
        cvtsd2ss %xmm7,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb331_three(%esp),%xmm4
        mulsd   %xmm7,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb331_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm7,%xmm2
        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm2       ## rsq*lu*lu 
        movapd nb331_three(%esp),%xmm4
        subsd %xmm2,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb331_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,nb331_rinvO(%esp)         ## rinvO in xmm4 
        mulsd   %xmm4,%xmm7
        movapd  %xmm7,nb331_rO(%esp)    ## r in xmm7 

        ## rsqH1 - seed in xmm2 
        cvtsd2ss %xmm6,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb331_three(%esp),%xmm4
        mulsd   %xmm6,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb331_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm6,%xmm2
        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm2       ## rsq*lu*lu 
        movapd nb331_three(%esp),%xmm4
        subsd %xmm2,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb331_half(%esp),%xmm4   ## rinv 
        movapd %xmm4,nb331_rinvH1(%esp)         ## rinvH1 
        mulsd  %xmm4,%xmm6
        movapd %xmm6,nb331_rH1(%esp)    ## rH1 

        ## rsqH2 - seed in xmm2 
        cvtsd2ss %xmm5,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb331_three(%esp),%xmm4
        mulsd   %xmm5,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb331_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm5,%xmm2
        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm2       ## rsq*lu*lu 
        movapd nb331_three(%esp),%xmm4
        subsd %xmm2,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb331_half(%esp),%xmm4   ## rinv 
        movapd %xmm4,nb331_rinvH2(%esp)   ## rinv 
        mulsd %xmm4,%xmm5
        movapd %xmm5,nb331_rH2(%esp)   ## r 

        ## do O interactions 
        movd %eax,%mm0
        ## rO is still in xmm7 
        mulsd   nb331_tsc(%esp),%xmm7
        cvttsd2si %xmm7,%eax    ## lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm7,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax
        movl nb331_VFtab(%ebp),%esi
        leal  (%eax,%eax,2),%eax ## idx *= 3 (total *=12 now)   

        movsd 0(%esi,%eax,8),%xmm4      ## Y1   
        movsd 8(%esi,%eax,8),%xmm5      ## F1   
        movsd 16(%esi,%eax,8),%xmm6     ## G1   
        movsd 24(%esi,%eax,8),%xmm7     ## H1   
        ## coulomb table ready, in xmm4-xmm7  
    mulsd  %xmm1,%xmm6      ## xmm6=Geps 
    mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
    addsd  %xmm6,%xmm5  ## F+Geps 
    addsd  %xmm7,%xmm5      ## xmm5=Fp=F+Geps+Heps2        
    mulsd  nb331_two(%esp),%xmm7         ## two*Heps2 
    addsd  %xmm6,%xmm7  ## Geps + 2*Heps2 
    addsd  %xmm5,%xmm7 ## xmm7=FF 
    mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addsd  %xmm4,%xmm5 ## xmm5=VV 

    movapd  nb331_qqO(%esp),%xmm0
    mulsd  %xmm0,%xmm5 ## vcoul=qq*VV  
    mulsd  %xmm7,%xmm0 ## fijC=FF*qq 
        ## at this point xmm5 contains vcoul and xmm0 fijC 
        ## increment vcoul - then we can get rid of mm5 
        ## update vctot 
        addsd  nb331_vctot(%esp),%xmm5
        movlpd %xmm5,nb331_vctot(%esp)

        ## Dispersion 
        movsd 32(%esi,%eax,8),%xmm4     ## Y1   
        movsd 40(%esi,%eax,8),%xmm5     ## F1   
        movsd 48(%esi,%eax,8),%xmm6     ## G1   
        movsd 56(%esi,%eax,8),%xmm7     ## H1   
        ## Dispersion table ready, in xmm4-xmm7                 
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb331_two(%esp),%xmm7    ## two*Heps2 
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 

        movapd nb331_c6(%esp),%xmm4
        mulsd  %xmm4,%xmm7       ## fijD 
        mulsd  %xmm4,%xmm5       ## Vvdw6 
        addsd  %xmm7,%xmm0      ## fscal summation 

        addsd  nb331_Vvdwtot(%esp),%xmm5
        movsd %xmm5,nb331_Vvdwtot(%esp)

        ## Repulsion 
        movsd 64(%esi,%eax,8),%xmm4     ## Y1   
        movsd 72(%esi,%eax,8),%xmm5     ## F1   
        movsd 80(%esi,%eax,8),%xmm6     ## G1   
        movsd 88(%esi,%eax,8),%xmm7     ## H1   
        ## Repulsion table ready, in xmm4-xmm7                  
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb331_two(%esp),%xmm7    ## two*Heps2 
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 

        movapd nb331_c12(%esp),%xmm4
        mulsd  %xmm4,%xmm7 ## fijR 
        mulsd  %xmm4,%xmm5 ## Vvdw12 
        addsd  %xmm0,%xmm7      ## fscal summation 
        mulsd nb331_tsc(%esp),%xmm7
        mulsd nb331_rinvO(%esp),%xmm7

        addsd  nb331_Vvdwtot(%esp),%xmm5
        movsd %xmm5,nb331_Vvdwtot(%esp)

        xorpd %xmm4,%xmm4
        subsd %xmm7,%xmm4

        movapd nb331_dxO(%esp),%xmm0
        movapd nb331_dyO(%esp),%xmm1
        movapd nb331_dzO(%esp),%xmm2
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2      ## tx in xmm0-xmm2 

        ## update O forces 
        movapd nb331_fixO(%esp),%xmm3
        movapd nb331_fiyO(%esp),%xmm4
        movapd nb331_fizO(%esp),%xmm7
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm7
        movlpd %xmm3,nb331_fixO(%esp)
        movlpd %xmm4,nb331_fiyO(%esp)
        movlpd %xmm7,nb331_fizO(%esp)
        ## update j forces with water O 
        movlpd %xmm0,nb331_fjx(%esp)
        movlpd %xmm1,nb331_fjy(%esp)
        movlpd %xmm2,nb331_fjz(%esp)

        ## Done with O interactions - now H1! 
        movapd nb331_rH1(%esp),%xmm7
        mulpd nb331_tsc(%esp),%xmm7
        cvttsd2si %xmm7,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subpd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb331_VFtab(%ebp),%esi
        leal  (%eax,%eax,2),%eax ## idx *= 3 (total *=12 now)   

        movsd (%esi,%eax,8),%xmm4       ## Y1   
        movsd 8(%esi,%eax,8),%xmm5      ## F1   
        movsd 16(%esi,%eax,8),%xmm6     ## G1   
        movsd 24(%esi,%eax,8),%xmm7     ## H1   

        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb331_two(%esp),%xmm7    ## two*Heps2 
        movapd nb331_qqH(%esp),%xmm3
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulsd  %xmm7,%xmm3 ## fijC=FF*qq 
    ## at this point mm5 contains vcoul and xmm3 fijC 
    ## increment vcoul 
        xorpd  %xmm4,%xmm4
    addsd  nb331_vctot(%esp),%xmm5
        mulsd  nb331_rinvH1(%esp),%xmm3
    movlpd %xmm5,nb331_vctot(%esp)
        mulsd  nb331_tsc(%esp),%xmm3
        subsd %xmm3,%xmm4

        movapd nb331_dxH1(%esp),%xmm0
        movapd nb331_dyH1(%esp),%xmm1
        movapd nb331_dzH1(%esp),%xmm2
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2

        ## update H1 forces 
        movapd nb331_fixH1(%esp),%xmm3
        movapd nb331_fiyH1(%esp),%xmm4
        movapd nb331_fizH1(%esp),%xmm7
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm7
        movlpd %xmm3,nb331_fixH1(%esp)
        movlpd %xmm4,nb331_fiyH1(%esp)
        movlpd %xmm7,nb331_fizH1(%esp)
        ## update j forces with water H1 
        addsd  nb331_fjx(%esp),%xmm0
        addsd  nb331_fjy(%esp),%xmm1
        addsd  nb331_fjz(%esp),%xmm2
        movlpd %xmm0,nb331_fjx(%esp)
        movlpd %xmm1,nb331_fjy(%esp)
        movlpd %xmm2,nb331_fjz(%esp)

        ## Done with H1, finally we do H2 interactions 
        movapd nb331_rH2(%esp),%xmm7
        mulsd   nb331_tsc(%esp),%xmm7
        cvttsd2si %xmm7,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb331_VFtab(%ebp),%esi
        leal  (%eax,%eax,2),%eax ## idx *= 3 (total *=12 now)   

        movsd (%esi,%eax,8),%xmm4       ## Y1   
        movsd 8(%esi,%eax,8),%xmm5      ## F1   
        movsd 16(%esi,%eax,8),%xmm6     ## G1   
        movsd 24(%esi,%eax,8),%xmm7     ## H1   
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb331_two(%esp),%xmm7    ## two*Heps2 
        movapd nb331_qqH(%esp),%xmm3
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulsd  %xmm7,%xmm3 ## fijC=FF*qq 
    ## at this point mm5 contains vcoul and xmm3 fijC 
    ## increment vcoul 
        xorpd  %xmm4,%xmm4
    addsd  nb331_vctot(%esp),%xmm5
        mulsd  nb331_rinvH2(%esp),%xmm3
    movlpd %xmm5,nb331_vctot(%esp)
        mulsd  nb331_tsc(%esp),%xmm3
        subsd  %xmm3,%xmm4

        movapd nb331_dxH2(%esp),%xmm0
        movapd nb331_dyH2(%esp),%xmm1
        movapd nb331_dzH2(%esp),%xmm2
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2

    movd %mm0,%eax

        ## update H2 forces 
        movapd nb331_fixH2(%esp),%xmm3
        movapd nb331_fiyH2(%esp),%xmm4
        movapd nb331_fizH2(%esp),%xmm7
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm7
        movlpd %xmm3,nb331_fixH2(%esp)
        movlpd %xmm4,nb331_fiyH2(%esp)
        movlpd %xmm7,nb331_fizH2(%esp)

        movl nb331_faction(%ebp),%edi
        ## update j forces 
        addsd  nb331_fjx(%esp),%xmm0
        addsd  nb331_fjy(%esp),%xmm1
        addsd  nb331_fjz(%esp),%xmm2

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

_nb_kernel331_ia32_sse2.nb331_updateouterdata: 
        movl  nb331_ii3(%esp),%ecx
        movl  nb331_faction(%ebp),%edi
        movl  nb331_fshift(%ebp),%esi
        movl  nb331_is3(%esp),%edx

        ## accumulate  Oi forces in xmm0, xmm1, xmm2 
        movapd nb331_fixO(%esp),%xmm0
        movapd nb331_fiyO(%esp),%xmm1
        movapd nb331_fizO(%esp),%xmm2

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addsd  %xmm3,%xmm0
        addsd  %xmm4,%xmm1
        addsd  %xmm5,%xmm2 ## sum is in low xmm0-xmm2 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5

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

        ## accumulate force in xmm6/xmm7 for fshift 
        movapd %xmm0,%xmm6
        movsd %xmm2,%xmm7
        unpcklpd %xmm1,%xmm6

        ## accumulate H1i forces in xmm0, xmm1, xmm2 
        movapd nb331_fixH1(%esp),%xmm0
        movapd nb331_fiyH1(%esp),%xmm1
        movapd nb331_fizH1(%esp),%xmm2

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addsd  %xmm3,%xmm0
        addsd  %xmm4,%xmm1
        addsd  %xmm5,%xmm2 ## sum is in low xmm0-xmm2 

        ## increment i force 
        movsd  24(%edi,%ecx,8),%xmm3
        movsd  32(%edi,%ecx,8),%xmm4
        movsd  40(%edi,%ecx,8),%xmm5
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm5
        movsd  %xmm3,24(%edi,%ecx,8)
        movsd  %xmm4,32(%edi,%ecx,8)
        movsd  %xmm5,40(%edi,%ecx,8)

        ## accumulate force in xmm6/xmm7 for fshift 
        addsd %xmm2,%xmm7
        unpcklpd %xmm1,%xmm0
        addpd %xmm0,%xmm6

        ## accumulate H2i forces in xmm0, xmm1, xmm2 
        movapd nb331_fixH2(%esp),%xmm0
        movapd nb331_fiyH2(%esp),%xmm1
        movapd nb331_fizH2(%esp),%xmm2

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addsd  %xmm3,%xmm0
        addsd  %xmm4,%xmm1
        addsd  %xmm5,%xmm2 ## sum is in low xmm0-xmm2 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5

        ## increment i force 
        movsd  48(%edi,%ecx,8),%xmm3
        movsd  56(%edi,%ecx,8),%xmm4
        movsd  64(%edi,%ecx,8),%xmm5
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm5
        movsd  %xmm3,48(%edi,%ecx,8)
        movsd  %xmm4,56(%edi,%ecx,8)
        movsd  %xmm5,64(%edi,%ecx,8)

        ## accumulate force in xmm6/xmm7 for fshift 
        addsd %xmm2,%xmm7
        unpcklpd %xmm1,%xmm0
        addpd %xmm0,%xmm6

        ## increment fshift force 
        movlpd (%esi,%edx,8),%xmm3
        movhpd 8(%esi,%edx,8),%xmm3
        movsd  16(%esi,%edx,8),%xmm4
        addpd  %xmm6,%xmm3
        addsd  %xmm7,%xmm4
        movlpd %xmm3,(%esi,%edx,8)
        movhpd %xmm3,8(%esi,%edx,8)
        movsd  %xmm4,16(%esi,%edx,8)

        ## get n from stack
        movl nb331_n(%esp),%esi
        ## get group index for i particle 
        movl  nb331_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb331_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb331_Vc(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## accumulate total lj energy and update it 
        movapd nb331_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb331_Vvdw(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## finish if last 
        movl nb331_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel331_ia32_sse2.nb331_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb331_n(%esp)
        jmp _nb_kernel331_ia32_sse2.nb331_outer
_nb_kernel331_ia32_sse2.nb331_outerend: 
        ## check if more outer neighborlists remain
        movl  nb331_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel331_ia32_sse2.nb331_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel331_ia32_sse2.nb331_threadloop
_nb_kernel331_ia32_sse2.nb331_end: 
        emms

        movl nb331_nouter(%esp),%eax
        movl nb331_ninner(%esp),%ebx
        movl nb331_outeriter(%ebp),%ecx
        movl nb331_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb331_salign(%esp),%eax
        addl %eax,%esp
        addl $812,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret




.globl nb_kernel331nf_ia32_sse2
.globl _nb_kernel331nf_ia32_sse2
nb_kernel331nf_ia32_sse2:       
_nb_kernel331nf_ia32_sse2:      
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
.set nb331nf_argkrf, 56
.set nb331nf_argcrf, 60
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
        ## bottom of stack is cache-aligned for sse use 
.set nb331nf_ixO, 0
.set nb331nf_iyO, 16
.set nb331nf_izO, 32
.set nb331nf_ixH1, 48
.set nb331nf_iyH1, 64
.set nb331nf_izH1, 80
.set nb331nf_ixH2, 96
.set nb331nf_iyH2, 112
.set nb331nf_izH2, 128
.set nb331nf_iqO, 144
.set nb331nf_iqH, 160
.set nb331nf_qqO, 176
.set nb331nf_qqH, 192
.set nb331nf_rinvO, 208
.set nb331nf_rinvH1, 224
.set nb331nf_rinvH2, 240
.set nb331nf_rO, 256
.set nb331nf_rH1, 272
.set nb331nf_rH2, 288
.set nb331nf_tsc, 304
.set nb331nf_c6, 320
.set nb331nf_c12, 336
.set nb331nf_vctot, 352
.set nb331nf_Vvdwtot, 368
.set nb331nf_half, 384
.set nb331nf_three, 400
.set nb331nf_is3, 416
.set nb331nf_ii3, 420
.set nb331nf_ntia, 424
.set nb331nf_innerjjnr, 428
.set nb331nf_innerk, 432
.set nb331nf_n, 436
.set nb331nf_nn1, 440
.set nb331nf_nri, 444
.set nb331nf_nouter, 448
.set nb331nf_ninner, 452
.set nb331nf_salign, 456
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $460,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb331nf_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb331nf_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb331nf_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb331nf_nouter(%esp)
        movl %eax,nb331nf_ninner(%esp)


        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double 0.5 IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb331nf_half(%esp)
        movl %ebx,nb331nf_half+4(%esp)
        movsd nb331nf_half(%esp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## 1.0
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## 2.0
        addpd  %xmm2,%xmm3      ## 3.0
        movapd %xmm1,nb331nf_half(%esp)
        movapd %xmm3,nb331nf_three(%esp)
        movl nb331nf_p_tabscale(%ebp),%eax
        movsd (%eax),%xmm5
        shufpd $0,%xmm5,%xmm5
        movapd %xmm5,nb331nf_tsc(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb331nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx),%ebx           ## ebx =ii 

        movl  nb331nf_charge(%ebp),%edx
        movsd (%edx,%ebx,8),%xmm3
        movsd 8(%edx,%ebx,8),%xmm4
        movl nb331nf_p_facel(%ebp),%esi
        movsd (%esi),%xmm5
        mulsd  %xmm5,%xmm3
        mulsd  %xmm5,%xmm4

        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        movapd %xmm3,nb331nf_iqO(%esp)
        movapd %xmm4,nb331nf_iqH(%esp)

        movl  nb331nf_type(%ebp),%edx
        movl  (%edx,%ebx,4),%ecx
        shll  %ecx
        movl nb331nf_p_ntype(%ebp),%edi
        imull (%edi),%ecx     ## ecx = ntia = 2*ntype*type[ii0] 
        movl  %ecx,nb331nf_ntia(%esp)
_nb_kernel331nf_ia32_sse2.nb331nf_threadloop: 
        movl  nb331nf_count(%ebp),%esi          ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel331nf_ia32_sse2.nb331nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                           ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel331nf_ia32_sse2.nb331nf_spinlock

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
        jg  _nb_kernel331nf_ia32_sse2.nb331nf_outerstart
        jmp _nb_kernel331nf_ia32_sse2.nb331nf_end

_nb_kernel331nf_ia32_sse2.nb331nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb331nf_nouter(%esp),%ebx
        movl %ebx,nb331nf_nouter(%esp)

_nb_kernel331nf_ia32_sse2.nb331nf_outer: 
        movl  nb331nf_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 

        movl  nb331nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movsd (%eax,%ebx,8),%xmm0
        movsd 8(%eax,%ebx,8),%xmm1
        movsd 16(%eax,%ebx,8),%xmm2

        movl  nb331nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx,%esi,4),%ebx            ## ebx =ii 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb331nf_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb331nf_ii3(%esp)

        addsd (%eax,%ebx,8),%xmm3
        addsd 8(%eax,%ebx,8),%xmm4
        addsd 16(%eax,%ebx,8),%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb331nf_ixO(%esp)
        movapd %xmm4,nb331nf_iyO(%esp)
        movapd %xmm5,nb331nf_izO(%esp)

        movsd %xmm0,%xmm3
        movsd %xmm1,%xmm4
        movsd %xmm2,%xmm5
        addsd 24(%eax,%ebx,8),%xmm0
        addsd 32(%eax,%ebx,8),%xmm1
        addsd 40(%eax,%ebx,8),%xmm2
        addsd 48(%eax,%ebx,8),%xmm3
        addsd 56(%eax,%ebx,8),%xmm4
        addsd 64(%eax,%ebx,8),%xmm5

        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm0,nb331nf_ixH1(%esp)
        movapd %xmm1,nb331nf_iyH1(%esp)
        movapd %xmm2,nb331nf_izH1(%esp)
        movapd %xmm3,nb331nf_ixH2(%esp)
        movapd %xmm4,nb331nf_iyH2(%esp)
        movapd %xmm5,nb331nf_izH2(%esp)

        ## clear vctot 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb331nf_vctot(%esp)
        movapd %xmm4,nb331nf_Vvdwtot(%esp)

        movl  nb331nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb331nf_pos(%ebp),%esi
        movl  nb331nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb331nf_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb331nf_ninner(%esp),%ecx
        movl  %ecx,nb331nf_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb331nf_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel331nf_ia32_sse2.nb331nf_unroll_loop
        jmp   _nb_kernel331nf_ia32_sse2.nb331nf_checksingle
_nb_kernel331nf_ia32_sse2.nb331nf_unroll_loop: 
        ## twice unrolled innerloop here 
        movl  nb331nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx

        addl $8,nb331nf_innerjjnr(%esp)             ## advance pointer (unrolled 2) 

        movl nb331nf_charge(%ebp),%esi     ## base of charge[] 

        movlpd (%esi,%eax,8),%xmm3
        movhpd (%esi,%ebx,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb331nf_iqO(%esp),%xmm3
        mulpd  nb331nf_iqH(%esp),%xmm4
        movapd  %xmm3,nb331nf_qqO(%esp)
        movapd  %xmm4,nb331nf_qqH(%esp)

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1
        movl nb331nf_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl (%esi,%ebx,4),%ebx
        movl nb331nf_vdwparam(%ebp),%esi
        shll %eax
        shll %ebx
        movl nb331nf_ntia(%esp),%edi
        addl %edi,%eax
        addl %edi,%ebx

        movlpd (%esi,%eax,8),%xmm6 # # c6a
        movlpd (%esi,%ebx,8),%xmm7 # # c6b
        movhpd 8(%esi,%eax,8),%xmm6    # # c6a c12a
        movhpd 8(%esi,%ebx,8),%xmm7    # # c6b c12b


        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movd  %mm0,%eax
        movd  %mm1,%ebx
        movapd %xmm4,nb331nf_c6(%esp)
        movapd %xmm6,nb331nf_c12(%esp)

        movl nb331nf_pos(%ebp),%esi        ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx

        ## move two coordinates to xmm0-xmm2    
        movlpd (%esi,%eax,8),%xmm0
        movlpd 8(%esi,%eax,8),%xmm1
        movlpd 16(%esi,%eax,8),%xmm2
        movhpd (%esi,%ebx,8),%xmm0
        movhpd 8(%esi,%ebx,8),%xmm1
        movhpd 16(%esi,%ebx,8),%xmm2

        ## move ixO-izO to xmm4-xmm6 
        movapd nb331nf_ixO(%esp),%xmm4
        movapd nb331nf_iyO(%esp),%xmm5
        movapd nb331nf_izO(%esp),%xmm6

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
        movapd %xmm4,%xmm7
        ## rsqO in xmm7 

        ## move ixH1-izH1 to xmm4-xmm6 
        movapd nb331nf_ixH1(%esp),%xmm4
        movapd nb331nf_iyH1(%esp),%xmm5
        movapd nb331nf_izH1(%esp),%xmm6

        ## calc dr 
        subpd %xmm0,%xmm4
        subpd %xmm1,%xmm5
        subpd %xmm2,%xmm6

        ## square it 
        mulpd %xmm4,%xmm4
        mulpd %xmm5,%xmm5
        mulpd %xmm6,%xmm6
        addpd %xmm5,%xmm6
        addpd %xmm4,%xmm6
        ## rsqH1 in xmm6 

        ## move ixH2-izH2 to xmm3-xmm5  
        movapd nb331nf_ixH2(%esp),%xmm3
        movapd nb331nf_iyH2(%esp),%xmm4
        movapd nb331nf_izH2(%esp),%xmm5

        ## calc dr 
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5

        ## square it 
        mulpd %xmm3,%xmm3
        mulpd %xmm4,%xmm4
        mulpd %xmm5,%xmm5
        addpd %xmm4,%xmm5
        addpd %xmm3,%xmm5
        ## rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

        ## start with rsqO - put seed in xmm2 
        cvtpd2ps %xmm7,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb331nf_three(%esp),%xmm4
        mulpd   %xmm7,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulpd   nb331nf_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm7,%xmm2
        movapd %xmm4,%xmm3
        mulpd %xmm4,%xmm4       ## lu*lu 
        mulpd %xmm4,%xmm2       ## rsq*lu*lu 
        movapd nb331nf_three(%esp),%xmm4
        subpd %xmm2,%xmm4       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulpd nb331nf_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,nb331nf_rinvO(%esp)       ## rinvO in xmm4 
        mulpd   %xmm4,%xmm7
        movapd  %xmm7,nb331nf_rO(%esp)          ## r in xmm7 

        ## rsqH1 - seed in xmm2 
        cvtpd2ps %xmm6,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb331nf_three(%esp),%xmm4
        mulpd   %xmm6,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulpd   nb331nf_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm6,%xmm2
        movapd %xmm4,%xmm3
        mulpd %xmm4,%xmm4       ## lu*lu 
        mulpd %xmm4,%xmm2       ## rsq*lu*lu 
        movapd nb331nf_three(%esp),%xmm4
        subpd %xmm2,%xmm4       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulpd nb331nf_half(%esp),%xmm4   ## rinv 
        movapd %xmm4,nb331nf_rinvH1(%esp)       ## rinvH1 
        mulpd  %xmm4,%xmm6
        movapd %xmm6,nb331nf_rH1(%esp)          ## rH1 

        ## rsqH2 - seed in xmm2 
        cvtpd2ps %xmm5,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb331nf_three(%esp),%xmm4
        mulpd   %xmm5,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulpd   nb331nf_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm5,%xmm2
        movapd %xmm4,%xmm3
        mulpd %xmm4,%xmm4       ## lu*lu 
        mulpd %xmm4,%xmm2       ## rsq*lu*lu 
        movapd nb331nf_three(%esp),%xmm4
        subpd %xmm2,%xmm4       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulpd nb331nf_half(%esp),%xmm4   ## rinv 
        movapd %xmm4,nb331nf_rinvH2(%esp)   ## rinv 
        mulpd %xmm4,%xmm5
        movapd %xmm5,nb331nf_rH2(%esp)   ## r 

        ## do O interactions 
        ## rO is still in xmm7 
        movapd nb331nf_rinvO(%esp),%xmm0
        mulpd   nb331nf_tsc(%esp),%xmm7
        cvttpd2pi %xmm7,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movd %eax,%mm0
        movd %ebx,%mm1
        movl nb331nf_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 
        leal  (%eax,%eax,2),%eax ## idx *= 3 (total *=12 now) 
        leal  (%ebx,%ebx,2),%ebx

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
    movapd nb331nf_qqO(%esp),%xmm0
    mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addpd  %xmm4,%xmm5 ## xmm5=VV 
    mulpd  %xmm0,%xmm5 ## vcoul=qq*VV  
        ## at this point xmm5 contains vcoul 
        ## increment vcoul - then we can get rid of mm5 
        ## update vctot 
        addpd  nb331nf_vctot(%esp),%xmm5
        movapd %xmm5,nb331nf_vctot(%esp)

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

        mulpd  nb331nf_c6(%esp),%xmm5    ## Vvdw6 

        addpd  nb331nf_Vvdwtot(%esp),%xmm5
        movapd %xmm5,nb331nf_Vvdwtot(%esp)

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

        mulpd  nb331nf_c12(%esp),%xmm5   ## Vvdw12 
        addpd  nb331nf_Vvdwtot(%esp),%xmm5
        movapd %xmm5,nb331nf_Vvdwtot(%esp)

        ## Done with O interactions - now H1! 
        movapd nb331nf_rH1(%esp),%xmm7
        mulpd nb331nf_tsc(%esp),%xmm7
        cvttpd2pi %xmm7,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb331nf_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 
        leal  (%eax,%eax,2),%eax ## idx *= 3 (total *=12 now)   
        leal  (%ebx,%ebx,2),%ebx

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
        movapd nb331nf_qqH(%esp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul
    ## increment vcoul 
    addpd  nb331nf_vctot(%esp),%xmm5
    movapd %xmm5,nb331nf_vctot(%esp)

        ## Done with H1, finally we do H2 interactions 
        movapd nb331nf_rH2(%esp),%xmm7
        mulpd   nb331nf_tsc(%esp),%xmm7
        cvttpd2pi %xmm7,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb331nf_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 
        leal  (%eax,%eax,2),%eax ## idx *= 3 (total *=12 now)
        leal  (%ebx,%ebx,2),%ebx

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
        movapd nb331nf_qqH(%esp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 
    ## increment vcoul 
    addpd  nb331nf_vctot(%esp),%xmm5
    movapd %xmm5,nb331nf_vctot(%esp)

        ## should we do one more iteration? 
        subl $2,nb331nf_innerk(%esp)
        jl    _nb_kernel331nf_ia32_sse2.nb331nf_checksingle
        jmp   _nb_kernel331nf_ia32_sse2.nb331nf_unroll_loop
_nb_kernel331nf_ia32_sse2.nb331nf_checksingle: 
        movl  nb331nf_innerk(%esp),%edx
        andl  $1,%edx
        jnz   _nb_kernel331nf_ia32_sse2.nb331nf_dosingle
        jmp   _nb_kernel331nf_ia32_sse2.nb331nf_updateouterdata
_nb_kernel331nf_ia32_sse2.nb331nf_dosingle: 
        movl  nb331nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax

        movl nb331nf_charge(%ebp),%esi     ## base of charge[] 
        xorpd %xmm3,%xmm3
        movlpd (%esi,%eax,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb331nf_iqO(%esp),%xmm3
        mulpd  nb331nf_iqH(%esp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movapd  %xmm3,nb331nf_qqO(%esp)
        movapd  %xmm4,nb331nf_qqH(%esp)

        movl nb331nf_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl nb331nf_vdwparam(%ebp),%esi
        shll %eax
        movl nb331nf_ntia(%esp),%edi
        addl %edi,%eax

        movlpd (%esi,%eax,8),%xmm6      ## c6a
        movhpd 8(%esi,%eax,8),%xmm6     ## c6a c12a 

        xorpd %xmm7,%xmm7
        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movd  %mm0,%eax
        movapd %xmm4,nb331nf_c6(%esp)
        movapd %xmm6,nb331nf_c12(%esp)

        movl nb331nf_pos(%ebp),%esi        ## base of pos[] 
        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 

        ## move coords to xmm0-xmm2 
        movlpd (%esi,%eax,8),%xmm0
        movlpd 8(%esi,%eax,8),%xmm1
        movlpd 16(%esi,%eax,8),%xmm2

        ## move ixO-izO to xmm4-xmm6 
        movapd nb331nf_ixO(%esp),%xmm4
        movapd nb331nf_iyO(%esp),%xmm5
        movapd nb331nf_izO(%esp),%xmm6

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
        movapd %xmm4,%xmm7
        ## rsqO in xmm7 

        ## move ixH1-izH1 to xmm4-xmm6 
        movapd nb331nf_ixH1(%esp),%xmm4
        movapd nb331nf_iyH1(%esp),%xmm5
        movapd nb331nf_izH1(%esp),%xmm6

        ## calc dr 
        subsd %xmm0,%xmm4
        subsd %xmm1,%xmm5
        subsd %xmm2,%xmm6

        ## square it 
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        mulsd %xmm6,%xmm6
        addsd %xmm5,%xmm6
        addsd %xmm4,%xmm6
        ## rsqH1 in xmm6 

        ## move ixH2-izH2 to xmm3-xmm5  
        movapd nb331nf_ixH2(%esp),%xmm3
        movapd nb331nf_iyH2(%esp),%xmm4
        movapd nb331nf_izH2(%esp),%xmm5

        ## calc dr 
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5

        ## square it 
        mulsd %xmm3,%xmm3
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        addsd %xmm4,%xmm5
        addsd %xmm3,%xmm5
        ## rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

        ## start with rsqO - put seed in xmm2 
        cvtsd2ss %xmm7,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb331nf_three(%esp),%xmm4
        mulsd   %xmm7,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb331nf_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm7,%xmm2
        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm2       ## rsq*lu*lu 
        movapd nb331nf_three(%esp),%xmm4
        subsd %xmm2,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb331nf_half(%esp),%xmm4   ## rinv 
        movapd  %xmm4,nb331nf_rinvO(%esp)       ## rinvO in xmm4 
        mulsd   %xmm4,%xmm7
        movapd  %xmm7,nb331nf_rO(%esp)          ## r in xmm7 

        ## rsqH1 - seed in xmm2 
        cvtsd2ss %xmm6,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb331nf_three(%esp),%xmm4
        mulsd   %xmm6,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb331nf_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm6,%xmm2
        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm2       ## rsq*lu*lu 
        movapd nb331nf_three(%esp),%xmm4
        subsd %xmm2,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb331nf_half(%esp),%xmm4   ## rinv 
        movapd %xmm4,nb331nf_rinvH1(%esp)       ## rinvH1 
        mulsd  %xmm4,%xmm6
        movapd %xmm6,nb331nf_rH1(%esp)          ## rH1 

        ## rsqH2 - seed in xmm2 
        cvtsd2ss %xmm5,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb331nf_three(%esp),%xmm4
        mulsd   %xmm5,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb331nf_half(%esp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm5,%xmm2
        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm2       ## rsq*lu*lu 
        movapd nb331nf_three(%esp),%xmm4
        subsd %xmm2,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb331nf_half(%esp),%xmm4   ## rinv 
        movapd %xmm4,nb331nf_rinvH2(%esp)   ## rinv 
        mulsd %xmm4,%xmm5
        movapd %xmm5,nb331nf_rH2(%esp)   ## r 

        ## do O interactions 
        movd %eax,%mm0
        ## rO is still in xmm7 
        mulsd   nb331nf_tsc(%esp),%xmm7
        cvttsd2si %xmm7,%eax    ## lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm7,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax
        movl nb331nf_VFtab(%ebp),%esi
        leal  (%eax,%eax,2),%eax ## idx *= 3 (total *=12 now)   

        movsd (%esi,%eax,8),%xmm4       ## Y1   
        movsd 8(%esi,%eax,8),%xmm5      ## F1   
        movsd 16(%esi,%eax,8),%xmm6     ## G1   
        movsd 24(%esi,%eax,8),%xmm7     ## H1   

        ## coulomb table ready, in xmm4-xmm7  
    mulsd  %xmm1,%xmm6      ## xmm6=Geps 
    mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
    addsd  %xmm6,%xmm5  ## F+Geps 
    addsd  %xmm7,%xmm5      ## xmm5=Fp=F+Geps+Heps2        
    mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addsd  %xmm4,%xmm5 ## xmm5=VV 

    movapd  nb331nf_qqO(%esp),%xmm0
    mulsd  %xmm0,%xmm5 ## vcoul=qq*VV  
        ## at this point xmm5 contains vcoul 
        ## increment vcoul - then we can get rid of mm5 
        ## update vctot 
        addsd  nb331nf_vctot(%esp),%xmm5
        movlpd %xmm5,nb331nf_vctot(%esp)

        ## Dispersion 
        movsd 32(%esi,%eax,8),%xmm4     ## Y1   
        movsd 40(%esi,%eax,8),%xmm5     ## F1   
        movsd 48(%esi,%eax,8),%xmm6     ## G1   
        movsd 56(%esi,%eax,8),%xmm7     ## H1   
        ## Dispersion table ready, in xmm4-xmm7                 
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 

        mulsd  nb331nf_c6(%esp),%xmm5    ## Vvdw6 

        addsd  nb331nf_Vvdwtot(%esp),%xmm5
        movsd %xmm5,nb331nf_Vvdwtot(%esp)

        ## Repulsion 
        movsd 64(%esi,%eax,8),%xmm4     ## Y1   
        movsd 72(%esi,%eax,8),%xmm5     ## F1   
        movsd 80(%esi,%eax,8),%xmm6     ## G1
        movsd 88(%esi,%eax,8),%xmm7     ## H1   
        ## Dispersion table ready, in xmm4-xmm7                 
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 

        mulsd  nb331nf_c12(%esp),%xmm5   ## Vvdw12 

        addsd  nb331nf_Vvdwtot(%esp),%xmm5
        movsd %xmm5,nb331nf_Vvdwtot(%esp)

        ## Done with O interactions - now H1! 
        movapd nb331nf_rH1(%esp),%xmm7
        mulpd nb331nf_tsc(%esp),%xmm7
        cvttsd2si %xmm7,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subpd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb331nf_VFtab(%ebp),%esi
        leal  (%eax,%eax,2),%eax ## idx *= 3 (total *=12 now)   

        movsd (%esi,%eax,8),%xmm4       ## Y1   
        movsd 8(%esi,%eax,8),%xmm5      ## F1   
        movsd 16(%esi,%eax,8),%xmm6     ## G1   
        movsd 24(%esi,%eax,8),%xmm7     ## H1   
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb331nf_qqH(%esp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 
    ## increment vcoul 
    addsd  nb331nf_vctot(%esp),%xmm5
    movlpd %xmm5,nb331nf_vctot(%esp)

        ## Done with H1, finally we do H2 interactions 
        movapd nb331nf_rH2(%esp),%xmm7
        mulsd   nb331nf_tsc(%esp),%xmm7
        cvttsd2si %xmm7,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb331nf_VFtab(%ebp),%esi
        leal  (%eax,%eax,2),%eax ## idx *= 3 (total *=12 now)   

        movsd (%esi,%eax,8),%xmm4       ## Y1   
        movsd 8(%esi,%eax,8),%xmm5      ## F1   
        movsd 16(%esi,%eax,8),%xmm6     ## G1   
        movsd 24(%esi,%eax,8),%xmm7     ## H1   
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb331nf_qqH(%esp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 
        ## increment vcoul 
        addsd  nb331nf_vctot(%esp),%xmm5
        movlpd %xmm5,nb331nf_vctot(%esp)

_nb_kernel331nf_ia32_sse2.nb331nf_updateouterdata: 
        ## get n from stack
        movl nb331nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb331nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb331nf_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb331nf_Vc(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## accumulate total lj energy and update it 
        movapd nb331nf_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb331nf_Vvdw(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## finish if last 
        movl nb331nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel331nf_ia32_sse2.nb331nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb331nf_n(%esp)
        jmp _nb_kernel331nf_ia32_sse2.nb331nf_outer
_nb_kernel331nf_ia32_sse2.nb331nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb331nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel331nf_ia32_sse2.nb331nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel331nf_ia32_sse2.nb331nf_threadloop
_nb_kernel331nf_ia32_sse2.nb331nf_end: 
        emms

        movl nb331nf_nouter(%esp),%eax
        movl nb331nf_ninner(%esp),%ebx
        movl nb331nf_outeriter(%ebp),%ecx
        movl nb331nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb331nf_salign(%esp),%eax
        addl %eax,%esp
        addl $460,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret



